#pragma once
#include "sce/lu.hpp"
namespace sce {
struct QrResult {
  index_type rows{}, columns{};
  // Thin Q (rows x columns) and R (columns x columns), column-major.
  std::vector<double> q, r;
};
// Independent Householder QR; 1 <= columns <= rows <= 256.
// Input is unchanged. Finite values and valid padded storage required.
// Arithmetic outside binary64 range throws overflow_error; no rank estimate.
QrResult factor_qr(std::span<const double> input, index_type rows,
                   index_type columns, index_type leading_dimension);
} // namespace sce
