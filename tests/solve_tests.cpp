#include "sce/solve.hpp"
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <random>
#include <stdexcept>
using V = std::vector<double>;
namespace {
V multiply(const V &a, const V &x, int n, int r) {
  V b(n * r, 0);
  for (int j = 0; j < r; ++j)
    for (int i = 0; i < n; ++i)
      for (int k = 0; k < n; ++k)
        b[j * n + i] += a[k * n + i] * x[j * n + k];
  return b;
}
} // namespace
TEST(Solve, PivotedMultipleRightHandSides) {
  V a{0, 2, 1, 1, 0, 3, 4, 1, 2};
  const V original = a;
  V expected{1, 2, 3, -2, 0, 4, 0, 0, 0};
  auto b = multiply(a, expected, 3, 3);
  auto f = sce::factor_lu(a, 3, 3, 3);
  const V factors = a, before = b;
  auto x = sce::solve_lu(a, 3, 3, f, b, 3, 3);
  for (std::size_t i = 0; i < x.size(); ++i)
    EXPECT_NEAR(x[i], expected[i], 1e-13);
  EXPECT_EQ(a, factors);
  EXPECT_EQ(b, before);
  for (int j = 0; j < 3; ++j) {
    auto single = sce::solve_lu(
        a, 3, 3, f, std::span<const double>(b).subspan(j * 3, 3), 1, 3);
    for (int i = 0; i < 3; ++i)
      EXPECT_DOUBLE_EQ(single[i], x[j * 3 + i]);
  }
}
TEST(Solve, SingularAndInvalidMetadata) {
  V a{1, 2, 2, 4}, b{1, 2};
  auto f = sce::factor_lu(a, 2, 2, 2);
  EXPECT_THROW(sce::solve_lu(a, 2, 2, f, b, 1, 2), std::domain_error);
  f.info = 0;
  EXPECT_THROW(sce::solve_lu(a, 2, 2, f, b, 1, 2), std::domain_error);
  f.pivots[0] = 2;
  EXPECT_THROW(sce::solve_lu(a, 2, 2, f, b, 1, 2), std::invalid_argument);
  f.pivots.clear();
  EXPECT_THROW(sce::solve_lu(a, 2, 2, f, b, 1, 2), std::invalid_argument);
}
TEST(Solve, CheckedInputsEmptyPaddingAndOverflow) {
  V empty;
  sce::LuResult f;
  EXPECT_TRUE(sce::solve_lu(empty, 0, 1, f, empty, 3, 1).empty());
  EXPECT_THROW(sce::solve_lu(empty, -1, 1, f, empty, 1, 1),
               std::invalid_argument);
  V a{2, 0, 99, 0, 4};
  auto lu = sce::factor_lu(a, 2, 2, 3);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  V b{2, 8, nan, 4, 16};
  auto x = sce::solve_lu(a, 2, 3, lu, b, 2, 3);
  EXPECT_EQ(x, (V{1, 2, 2, 4}));
  EXPECT_TRUE(std::isnan(b[2]));
  EXPECT_THROW(sce::solve_lu(a, 2, 3, lu, b, 2, 1), std::invalid_argument);
  EXPECT_THROW(
      sce::solve_lu(a, 2, 3, lu, std::span<const double>(b).first(1), 1, 2),
      std::invalid_argument);
  b[0] = nan;
  EXPECT_THROW(sce::solve_lu(a, 2, 3, lu, b, 2, 3), std::invalid_argument);
  V tiny{std::numeric_limits<double>::min()},
      huge{std::numeric_limits<double>::max()};
  auto tf = sce::factor_lu(tiny, 1, 1, 1);
  EXPECT_THROW(sce::solve_lu(tiny, 1, 1, tf, huge, 1, 1), std::overflow_error);
}
TEST(Solve, SeededResidualAndForwardError) {
  std::mt19937 gen(20261026);
  std::uniform_real_distribution<double> d(-1, 1);
  for (int n = 1; n <= 24; ++n) {
    V a(n * n), expected(n * 4);
    for (double &v : a)
      v = d(gen);
    // Controlled well-conditioned family for a meaningful forward error check.
    for (int i = 0; i < n; ++i)
      a[i * n + i] += 2 * n;
    for (double &v : expected)
      v = d(gen);
    auto b = multiply(a, expected, n, 4);
    V packed = a;
    auto f = sce::factor_lu(packed, n, n, n);
    auto x = sce::solve_lu(packed, n, n, f, b, 4, n);
    auto ax = multiply(a, x, n, 4);
    double anorm = 0;
    for (int i = 0; i < n; ++i) {
      double row = 0;
      for (int k = 0; k < n; ++k)
        row += std::abs(a[k * n + i]);
      anorm = std::max(anorm, row);
    }
    for (int j = 0; j < 4; ++j) {
      double error = 0, bnorm = 0, xnorm = 0;
      for (int i = 0; i < n; ++i) {
        error = std::max(error, std::abs(ax[j * n + i] - b[j * n + i]));
        bnorm = std::max(bnorm, std::abs(b[j * n + i]));
        xnorm = std::max(xnorm, std::abs(x[j * n + i]));
        EXPECT_NEAR(x[j * n + i], expected[j * n + i], 2e-14);
      }
      EXPECT_LE(error / (anorm * xnorm + bnorm),
                32 * n * std::numeric_limits<double>::epsilon());
    }
  }
}
TEST(Solve, IllConditionedResidualIsNotForwardAccuracy) {
  constexpr int n = 10;
  V a(n * n), ones(n, 1);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      a[j * n + i] = 1.0 / (i + j + 1);
  auto b = multiply(a, ones, n, 1);
  V packed = a;
  auto f = sce::factor_lu(packed, n, n, n);
  auto x = sce::solve_lu(packed, n, n, f, b, 1, n);
  auto ax = multiply(a, x, n, 1);
  for (int i = 0; i < n; ++i)
    EXPECT_NEAR(ax[i], b[i], 2e-14);
  // Deliberately no small forward-error promise for this ill-conditioned
  // matrix.
}
