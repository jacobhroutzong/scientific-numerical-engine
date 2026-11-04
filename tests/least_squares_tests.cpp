#include "sce/least_squares.hpp"
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
TEST(LeastSquares, InconsistentFitAndMultipleRhs) {
  // Design [1,t], t=0,1,2. Fits b=[1,2,2] with x=[7/6,1/2].
  std::vector<double> a{1, 1, 1, 0, 1, 2}, b{1, 2, 2, 2, 4, 4};
  const auto before = a;
  auto r = sce::least_squares(a, 3, 2, 3, b, 2, 3, {std::nullopt, true});
  EXPECT_EQ(a, before);
  EXPECT_EQ(b, (std::vector<double>{1, 2, 2, 2, 4, 4}));
  ASSERT_EQ(r.x.size(), 4u);
  EXPECT_NEAR(r.x[0], 7. / 6, 1e-14);
  EXPECT_NEAR(r.x[1], .5, 1e-14);
  EXPECT_NEAR(r.x[2], 7. / 3, 1e-14);
  EXPECT_NEAR(r.x[3], 1, 1e-14);
  EXPECT_NEAR(r.diagnostics[0].residual_norm2, std::sqrt(1. / 6), 1e-14);
  EXPECT_GT(r.diagnostics[0].scaled_residual, 0);
  EXPECT_LT(r.diagnostics[0].normalized_stationarity, 1e-13);
  ASSERT_TRUE(r.r_condition_inf);
  EXPECT_GT(*r.r_condition_inf, 1);
}
TEST(LeastSquares, RankPolicyAndConditionAreDistinct) {
  EXPECT_THROW(sce::least_squares(std::vector<double>{1, 2, 2, 4}, 2, 2, 2,
                                  std::vector<double>{1, 2}, 1, 2),
               sce::RankDeficient);
  const std::vector<double> a{1, 0, 0, 1e-17}, b{1, 1e-17};
  EXPECT_THROW(sce::least_squares(a, 2, 2, 2, b, 1, 2), sce::RankDeficient);
  auto r = sce::least_squares(a, 2, 2, 2, b, 1, 2, {0., true});
  EXPECT_DOUBLE_EQ(r.x[0], 1);
  EXPECT_DOUBLE_EQ(r.x[1], 1);
  EXPECT_NEAR(*r.r_condition_inf / 1e17, 1, 1e-14);
}
TEST(LeastSquares, PaddingZeroRhsAndScale) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  auto r = sce::least_squares(std::vector<double>{1, 0, nan, 0, 1}, 2, 2, 3,
                              std::vector<double>{0, 0}, 1, 2);
  EXPECT_EQ(r.x, (std::vector<double>{0, 0}));
  EXPECT_DOUBLE_EQ(r.diagnostics[0].scaled_residual, 0);
  for (double scale : {1e-150, 1e150}) {
    auto q =
        sce::least_squares(std::vector<double>{scale, 0, 0, 2 * scale}, 2, 2, 2,
                           std::vector<double>{3 * scale, 4 * scale}, 1, 2);
    EXPECT_NEAR(q.x[0], 3, 1e-14);
    EXPECT_NEAR(q.x[1], 2, 1e-14);
  }
}
TEST(LeastSquares, InvalidInputsAndOverflow) {
  auto a = std::vector<double>{1}, b = a;
  EXPECT_THROW(sce::least_squares(a, 1, 2, 1, b, 1, 1), std::invalid_argument);
  EXPECT_THROW(sce::least_squares({}, 1, 1, 1, b, 1, 1), std::invalid_argument);
  EXPECT_THROW(sce::least_squares(a, 1, 1, 0, b, 1, 1), std::invalid_argument);
  EXPECT_THROW(sce::least_squares(a, 1, 1, 1, b, 0, 1), std::invalid_argument);
  EXPECT_THROW(sce::least_squares(a, 1, 1, 1, b, 1, 1, {-1., false}),
               std::invalid_argument);
  b[0] = std::numeric_limits<double>::infinity();
  EXPECT_THROW(sce::least_squares(a, 1, 1, 1, b, 1, 1), std::invalid_argument);
  EXPECT_THROW(sce::least_squares(std::vector<double>{1e-300}, 1, 1, 1,
                                  std::vector<double>{1e300}, 1, 1),
               std::overflow_error);
}
