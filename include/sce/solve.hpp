#pragma once
#include "sce/lu.hpp"

namespace sce {
// Solve A*X=B using square packed LU and sequential pivots from factor_lu.
// Inputs are immutable column-major views; returned X is dense, ld=order.
// Invalid shape/pivots/nonfinite input: invalid_argument. Exact singularity:
// domain_error. Nonfinite arithmetic/output: overflow_error. No condition
// estimate. This is new mathematical implementation, not recovered DGETRS
// machine code.
std::vector<double> solve_lu(std::span<const double> packed, index_type order,
                             index_type leading_dimension,
                             const LuResult &factorization,
                             std::span<const double> rhs, index_type rhs_count,
                             index_type rhs_leading_dimension);
} // namespace sce
