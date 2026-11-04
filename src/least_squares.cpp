#include "sce/least_squares.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace sce {
namespace {
double checked(double value) {
  if (!std::isfinite(value))
    throw std::overflow_error(
        "Least-squares arithmetic or diagnostic overflow");
  return value;
}
std::size_t at(int row, int col, int ld) {
  return static_cast<std::size_t>(col) * ld + row;
}
void validate(std::span<const double> v, int rows, int columns, int ld) {
  if (ld < rows || static_cast<std::size_t>(columns - 1) >
                       (std::numeric_limits<std::size_t>::max() - rows) /
                           static_cast<std::size_t>(ld))
    throw std::invalid_argument("Invalid least-squares leading dimension");
  if (v.size() < at(rows - 1, columns - 1, ld) + 1)
    throw std::invalid_argument("Insufficient least-squares storage");
  for (int j = 0; j < columns; ++j)
    for (int i = 0; i < rows; ++i)
      if (!std::isfinite(v[at(i, j, ld)]))
        throw std::invalid_argument(
            "Least squares requires finite active entries");
}
void back_solve(const std::vector<double> &r, int n, std::span<double> x) {
  for (int i = n - 1; i >= 0; --i) {
    double value = x[i];
    for (int j = i + 1; j < n; ++j)
      value = checked(value - checked(r[at(i, j, n)] * x[j]));
    x[i] = checked(value / r[at(i, i, n)]);
  }
}
double condition(const std::vector<double> &r, int n) {
  double norm = 0;
  std::vector<double> inverse_rows(n);
  for (int i = 0; i < n; ++i) {
    double row = 0;
    for (int j = i; j < n; ++j)
      row = checked(row + std::abs(r[at(i, j, n)]));
    norm = std::max(norm, row);
  }
  for (int j = 0; j < n; ++j) {
    std::vector<double> column(n);
    column[j] = 1;
    back_solve(r, n, column);
    for (int i = 0; i < n; ++i)
      inverse_rows[i] = checked(inverse_rows[i] + std::abs(column[i]));
  }
  return checked(norm *
                 *std::max_element(inverse_rows.begin(), inverse_rows.end()));
}
} // namespace
LeastSquaresResult least_squares(std::span<const double> a, index_type m,
                                 index_type n, index_type lda,
                                 std::span<const double> b, index_type nrhs,
                                 index_type ldb, LeastSquaresOptions options) {
  if (n < 1 || m < n || m > 256 || nrhs < 1 || nrhs > 16)
    throw std::invalid_argument(
        "Least squares requires 1<=n<=m<=256, 1<=RHS<=16");
  const double tolerance = options.relative_rank_tolerance.value_or(
      std::numeric_limits<double>::epsilon() * m);
  if (!std::isfinite(tolerance) || tolerance < 0 || tolerance >= 1)
    throw std::invalid_argument(
        "Relative rank tolerance must be finite in [0,1)");
  validate(a, m, n, lda);
  validate(b, m, nrhs, ldb);
  const auto qr = factor_qr(a, m, n, lda);
  double scale = 0, minimum = 1;
  for (double value : qr.r)
    scale = std::max(scale, std::abs(checked(value)));
  if (scale == 0)
    throw RankDeficient("Zero QR triangular factor");
  for (int i = 0; i < n; ++i) {
    const double relative = std::abs(qr.r[at(i, i, n)]) / scale;
    minimum = std::min(minimum, relative);
    if (relative <= tolerance)
      throw RankDeficient(
          "R diagonal rejected by unpivoted relative-rank policy");
  }
  LeastSquaresResult result;
  result.relative_rank_tolerance = tolerance;
  result.minimum_relative_r_diagonal = minimum;
  result.x.resize(static_cast<std::size_t>(n) * nrhs);
  for (int rhs = 0; rhs < nrhs; ++rhs) {
    auto x = std::span<double>(result.x).subspan(
        static_cast<std::size_t>(rhs) * n, n);
    for (int i = 0; i < n; ++i) {
      double dot = 0;
      for (int k = 0; k < m; ++k)
        dot = checked(dot + checked(qr.q[at(k, i, m)] * b[at(k, rhs, ldb)]));
      x[i] = dot;
    }
    back_solve(qr.r, n, x);
  }
  double anorm = 0;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i)
      anorm = checked(std::hypot(anorm, a[at(i, j, lda)]));
  for (int rhs = 0; rhs < nrhs; ++rhs) {
    double bnorm = 0, xnorm = 0, rnorm = 0, gradient = 0;
    std::vector<double> residual(m);
    for (int j = 0; j < n; ++j)
      xnorm = checked(std::hypot(xnorm, result.x[at(j, rhs, n)]));
    for (int i = 0; i < m; ++i) {
      double value = -b[at(i, rhs, ldb)];
      bnorm = checked(std::hypot(bnorm, b[at(i, rhs, ldb)]));
      for (int j = 0; j < n; ++j)
        value = checked(value +
                        checked(a[at(i, j, lda)] * result.x[at(j, rhs, n)]));
      residual[i] = value;
      rnorm = checked(std::hypot(rnorm, value));
    }
    for (int j = 0; j < n; ++j) {
      double dot = 0;
      // Normalize before the dot product to avoid squaring the input scale.
      if (rnorm != 0)
        for (int i = 0; i < m; ++i)
          dot =
              checked(dot + (a[at(i, j, lda)] / anorm) * (residual[i] / rnorm));
      gradient = checked(std::hypot(gradient, dot));
    }
    // Compute r/(anorm*xnorm+bnorm) without forming an overflowing product.
    double scaled;
    if (xnorm > 1)
      scaled = checked((rnorm / xnorm) / checked(anorm + bnorm / xnorm));
    else {
      const double denominator = checked(anorm * xnorm + bnorm);
      scaled = denominator == 0 ? 0 : checked(rnorm / denominator);
    }
    result.diagnostics.push_back({rnorm, scaled, gradient});
  }
  if (options.estimate_r_condition)
    result.r_condition_inf = condition(qr.r, n);
  return result;
}
} // namespace sce
