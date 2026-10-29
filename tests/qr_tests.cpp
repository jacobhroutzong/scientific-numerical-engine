#include "sce/qr.hpp"
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <random>
#include <stdexcept>

namespace {
void check(const std::vector<double> &a, int m, int n,
           double tolerance = 2e-12) {
  const auto original = a;
  const auto qr = sce::factor_qr(a, m, n, m);
  EXPECT_EQ(a, original);
  double scale = 0;
  for (auto value : a)
    scale = std::max(scale, std::abs(value));
  if (scale == 0)
    scale = 1;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i) {
      double value = 0;
      for (int k = 0; k < n; ++k)
        value += qr.q[k * m + i] * qr.r[j * n + k];
      EXPECT_NEAR(value / scale, a[j * m + i] / scale, tolerance);
    }
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      double dot = 0;
      for (int k = 0; k < m; ++k)
        dot += qr.q[i * m + k] * qr.q[j * m + k];
      EXPECT_NEAR(dot, i == j ? 1.0 : 0.0, tolerance);
    }
}
} // namespace
TEST(Qr, ReconstructionAndOrthogonality) {
  std::mt19937 gen(2468);
  std::uniform_real_distribution<double> random(-1, 1);
  for (int m : {1, 2, 5, 16, 32})
    for (int n : {1, m}) {
      std::vector<double> a(m * n);
      for (auto &value : a)
        value = random(gen);
      check(a, m, n);
    }
  check({12, 6, -4, -51, 167, 24, 4, -68, -41}, 3, 3);
  check({1, 2, 3, 4, 2, 4, 6, 8}, 4, 2); // dependent columns
  check(std::vector<double>(12), 4, 3);
}
TEST(Qr, ScalingStorageAndValidation) {
  check({1e-300, 2e-300, 3e-300, 4e-300}, 2, 2);
  check({1e200, 2e200, 3e200, 4e200}, 2, 2);
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  auto padded = sce::factor_qr(std::vector<double>{1, 2, nan, 3, 4}, 2, 2, 3);
  EXPECT_EQ(padded.q.size(), 4u);
  EXPECT_THROW(sce::factor_qr(std::vector<double>{1}, 2, 2, 2),
               std::invalid_argument);
  EXPECT_THROW(sce::factor_qr(std::vector<double>{nan}, 1, 1, 1),
               std::invalid_argument);
  EXPECT_THROW(sce::factor_qr({}, 1, 2, 1), std::invalid_argument);
  EXPECT_THROW(sce::factor_qr(std::vector<double>{1.7e308, 1.7e308}, 2, 1, 2),
               std::overflow_error);
}
