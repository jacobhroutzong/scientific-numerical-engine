#include "sce/lu.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace sce {
namespace {
static_assert(std::numeric_limits<double>::is_iec559 &&
              std::numeric_limits<double>::digits == 53 && sizeof(double) == 8);
using I = index_type;
std::size_t offset(I row, I col, I ld) {
    return static_cast<std::size_t>(col) * static_cast<std::size_t>(ld) +
           static_cast<std::size_t>(row);
}
void swaps(double* a, I columns, I ld, std::span<I> pivots, I first, I end) {
    for (I k = first; k < end; ++k)
        if (pivots[k] != k)
            for (I col = 0; col < columns; ++col)
                std::swap(a[offset(k, col, ld)], a[offset(pivots[k], col, ld)]);
}

// Evidence: DGETRF2 at 0x2e066ace0; complete dossier and constants in docs.
I factor(double* a, I m, I n, I ld, std::span<I> pivots) {
    if (m == 0 || n == 0) return 0;
    if (m == 1) {
        pivots[0] = 0;
        return a[0] == 0.0 ? 1 : 0;
    }
    if (n == 1) {
        I best = 0;
        for (I row = 1; row < m; ++row)
            if (std::abs(a[row]) > std::abs(a[best])) best = row;
        pivots[0] = best;
        if (a[best] == 0.0) return 1;
        std::swap(a[0], a[best]);
        const double diagonal = a[0];
        // DLAMCH('S') recovered bits: 0x0010000000000000.
        if (std::abs(diagonal) >= std::numeric_limits<double>::min()) {
            const double reciprocal = 1.0 / diagonal;
            for (I row = 1; row < m; ++row) a[row] *= reciprocal;
        } else {
            for (I row = 1; row < m; ++row) a[row] /= diagonal;
        }
        return 0;
    }

    const I split = std::min(m, n) / 2;
    const I width = n - split;
    const I first_info = factor(a, m, split, ld, pivots.first(split));
    // DLASWP on trailing columns using the first panel's row swaps.
    swaps(a + offset(0, split, ld), width, ld, pivots, 0, split);

    // Scalar DTRSM(Left,Lower,NoTranspose,Unit): L11 * U12 = A12.
    for (I col = split; col < n; ++col)
        for (I k = 0; k < split; ++k) {
            const double value = a[offset(k, col, ld)];
            if (value != 0.0)
                for (I row = k + 1; row < split; ++row)
                    a[offset(row, col, ld)] -= value * a[offset(row, k, ld)];
        }

    // Scalar DGEMM update: A22 -= L21 * U12. SIMD reduction order is not copied.
    for (I col = split; col < n; ++col)
        for (I k = 0; k < split; ++k) {
            const double negative_upper = -a[offset(k, col, ld)];
            for (I row = split; row < m; ++row)
                a[offset(row, col, ld)] += a[offset(row, k, ld)] * negative_upper;
        }

    const I second_info = factor(a + offset(split, split, ld), m - split,
                                 width, ld, pivots.subspan(split));
    const I count = std::min(m, n);
    for (I k = split; k < count; ++k) pivots[k] += split;
    swaps(a, split, ld, pivots, split, count);
    return first_info != 0 ? first_info : (second_info == 0 ? 0 : split + second_info);
}
}

LuResult factor_lu(std::span<double> storage, I rows, I columns, I ld) {
    if (rows < 0 || columns < 0 || ld < std::max<I>(1, rows))
        throw std::invalid_argument("Invalid LU dimensions or leading dimension");
    if (rows != 0 && columns != 0) {
        const auto last_column = static_cast<std::size_t>(columns - 1);
        if (last_column > (std::numeric_limits<std::size_t>::max() - rows) /
                              static_cast<std::size_t>(ld))
            throw std::invalid_argument("LU storage size overflow");
        const auto required = last_column * static_cast<std::size_t>(ld) + rows;
        if (storage.size() < required) throw std::invalid_argument("LU storage too small");
        for (I col = 0; col < columns; ++col)
            for (I row = 0; row < rows; ++row)
                if (!std::isfinite(storage[offset(row, col, ld)]))
                    throw std::invalid_argument("LU supports finite input values only");
    }
    LuResult result;
    result.pivots.resize(std::min(rows, columns));
    result.info = factor(storage.data(), rows, columns, ld, result.pivots);
    return result;
}
}
