#pragma once
#include "sce/qr.hpp"
#include <optional>
#include <stdexcept>

namespace sce {
struct LeastSquaresOptions {
  std::optional<double> relative_rank_tolerance;
  bool estimate_r_condition = false;
};
class RankDeficient : public std::domain_error {
public:
  using std::domain_error::domain_error;
};
struct FitDiagnostics {
  double residual_norm2{};
  double scaled_residual{};
  double normalized_stationarity{};
};
struct LeastSquaresResult {
  std::vector<double> x;                   // n x nrhs, dense column-major
  std::vector<FitDiagnostics> diagnostics; // one per RHS
  double minimum_relative_r_diagonal{};
  double relative_rank_tolerance{};
  std::optional<double> r_condition_inf;
};
// Full-column-rank policy; 1<=n<=m<=256, 1<=nrhs<=16, immutable inputs.
// Throws invalid_argument (boundary), RankDeficient (diagonal threshold),
// overflow_error (nonfinite arithmetic). See docs/phase5-least-squares.md.
LeastSquaresResult least_squares(std::span<const double> a, index_type m,
                                 index_type n, index_type lda,
                                 std::span<const double> b, index_type nrhs,
                                 index_type ldb,
                                 LeastSquaresOptions options = {});
} // namespace sce
