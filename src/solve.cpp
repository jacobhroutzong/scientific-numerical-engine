#include "sce/solve.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace sce {
namespace {
using I = index_type;
std::size_t at(I row, I col, I ld) {
  return static_cast<std::size_t>(col) * static_cast<std::size_t>(ld) + row;
}
void check(std::span<const double> a, I rows, I cols, I ld) {
  if (rows < 0 || cols < 0 || ld < std::max<I>(1, rows))
    throw std::invalid_argument("Invalid solve dimensions");
  if (rows == 0 || cols == 0)
    return;
  const auto c = static_cast<std::size_t>(cols - 1);
  if (c > (std::numeric_limits<std::size_t>::max() - rows) /
              static_cast<std::size_t>(ld) ||
      a.size() < c * static_cast<std::size_t>(ld) + rows)
    throw std::invalid_argument("Solve storage too small or overflowing");
  for (I j = 0; j < cols; ++j)
    for (I i = 0; i < rows; ++i)
      if (!std::isfinite(a[at(i, j, ld)]))
        throw std::invalid_argument("Solve requires finite logical cells");
}
double finite(double x) {
  if (!std::isfinite(x))
    throw std::overflow_error("Nonfinite solve arithmetic");
  return x;
}
} // namespace
std::vector<double> solve_lu(std::span<const double> a, I n, I ld,
                             const LuResult &f, std::span<const double> b,
                             I nrhs, I ldb) {
  check(a, n, n, ld);
  check(b, n, nrhs, ldb);
  if (f.pivots.size() != static_cast<std::size_t>(n) || f.info < 0 ||
      f.info > n)
    throw std::invalid_argument("Inconsistent LU metadata");
  for (I k = 0; k < n; ++k)
    if (f.pivots[k] < k || f.pivots[k] >= n)
      throw std::invalid_argument("Invalid sequential LU pivot");
  if (f.info > 0)
    throw std::domain_error("Exactly singular LU factors");
  for (I k = 0; k < n; ++k)
    if (a[at(k, k, ld)] == 0.0)
      throw std::domain_error("Exactly zero U diagonal");
  if (n != 0 && static_cast<std::size_t>(nrhs) >
                    std::numeric_limits<std::size_t>::max() / n)
    throw std::invalid_argument("Solve output size overflow");
  std::vector<double> x(static_cast<std::size_t>(n) * nrhs);
  for (I j = 0; j < nrhs; ++j) {
    for (I i = 0; i < n; ++i)
      x[at(i, j, n)] = b[at(i, j, ldb)];
    // P*A=L*U, so apply the same forward sequence of swaps to B.
    for (I k = 0; k < n; ++k)
      std::swap(x[at(k, j, n)], x[at(f.pivots[k], j, n)]);
    for (I i = 0; i < n; ++i)
      for (I k = 0; k < i; ++k)
        x[at(i, j, n)] =
            finite(x[at(i, j, n)] - a[at(i, k, ld)] * x[at(k, j, n)]);
    for (I i = n; i-- > 0;) {
      for (I k = i + 1; k < n; ++k)
        x[at(i, j, n)] =
            finite(x[at(i, j, n)] - a[at(i, k, ld)] * x[at(k, j, n)]);
      x[at(i, j, n)] = finite(x[at(i, j, n)] / a[at(i, i, ld)]);
    }
  }
  return x;
}
} // namespace sce
