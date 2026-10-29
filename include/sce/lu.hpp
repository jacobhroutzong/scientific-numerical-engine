#pragma once
#include <cstdint>
#include <span>
#include <vector>

namespace sce {
using index_type = std::int32_t;

struct LuResult {
    // Sequential zero-based row swaps, not a final permutation vector.
    std::vector<index_type> pivots;
    // 0: no exact zero pivot; positive: first zero U diagonal, one-based.
    index_type info = 0;
};

// In-place column-major recursive real LU, reconstructing DGETRF2's control flow.
// Packed output: strict lower part is L, its diagonal is implicit 1; upper is U.
// Applying pivots in order to input A gives P*A = L*U.
// Invalid dimensions/storage and non-finite inputs throw std::invalid_argument.
// That checked C++ boundary replaces the original XERBLA/error-argument protocol.
LuResult factor_lu(std::span<double> storage, index_type rows,
                   index_type columns, index_type leading_dimension);
}
