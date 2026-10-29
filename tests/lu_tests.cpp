#include "sce/lu.hpp"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace {
using I = sce::index_type;
using V = std::vector<double>;
// Independent mathematical oracle: form explicit P*A and L*U in long double.
void residual(const V& original, const V& packed, I m, I n, I ld,
              const sce::LuResult& result) {
    const I k = std::min(m, n);
    ASSERT_EQ(result.pivots.size(), static_cast<std::size_t>(k));
    V permuted = original;
    for (I j = 0; j < k; ++j) {
        ASSERT_GE(result.pivots[j], j);
        ASSERT_LT(result.pivots[j], m);
        for (I col = 0; col < n; ++col)
            std::swap(permuted[static_cast<std::size_t>(col)*ld+j],
                      permuted[static_cast<std::size_t>(col)*ld+result.pivots[j]]);
    }
    long double error = 0, norm = 0;
    for (I row = 0; row < m; ++row) {
        long double row_error = 0, row_norm = 0;
        for (I col = 0; col < n; ++col) {
            long double product = 0;
            for (I inner = 0; inner < k; ++inner) {
                const long double lower = row == inner ? 1.0L :
                    (row > inner ? packed[static_cast<std::size_t>(inner)*ld+row] : 0.0L);
                const long double upper = inner <= col ?
                    packed[static_cast<std::size_t>(col)*ld+inner] : 0.0L;
                product += lower * upper;
            }
            row_error += std::abs(permuted[static_cast<std::size_t>(col)*ld+row] - product);
            row_norm += std::abs(original[static_cast<std::size_t>(col)*ld+row]);
        }
        error = std::max(error, row_error); norm = std::max(norm, row_norm);
    }
    const auto tolerance = 128.0L * std::numeric_limits<double>::epsilon() *
        std::max<I>(1, k) * std::max(norm, static_cast<long double>(std::numeric_limits<double>::min()));
    EXPECT_LE(error, tolerance) << "m=" << m << " n=" << n;
}

TEST(Lu, KnownPivotedExample) {
    V a{1,3,2,4}; const V original=a;
    auto r=sce::factor_lu(a,2,2,2);
    EXPECT_EQ(r.info,0); EXPECT_EQ(r.pivots,(std::vector<I>{1,1}));
    EXPECT_DOUBLE_EQ(a[0],3); EXPECT_DOUBLE_EQ(a[1],1.0/3);
    EXPECT_DOUBLE_EQ(a[2],4); EXPECT_NEAR(a[3],2.0/3,1e-15);
    residual(original,a,2,2,2,r);
}
TEST(Lu, EmptyShapes) {
    V a; EXPECT_TRUE(sce::factor_lu(a,0,5,1).pivots.empty());
    EXPECT_TRUE(sce::factor_lu(a,4,0,4).pivots.empty());
    EXPECT_EQ(sce::factor_lu(a,0,0,1).info,0);
}
TEST(Lu, OneRowSingularity) {
    V a{0,4,7}; auto r=sce::factor_lu(a,1,3,1);
    EXPECT_EQ(r.info,1); EXPECT_EQ(r.pivots[0],0); EXPECT_EQ(a,(V{0,4,7}));
}
TEST(Lu, ColumnPivotAndTies) {
    V a{-1,-4,4,2}; auto r=sce::factor_lu(a,4,1,4);
    EXPECT_EQ(r.info,0); EXPECT_EQ(r.pivots[0],1); // Chosen first-maximum policy.
    EXPECT_EQ(a,(V{-4,.25,-1,-.5}));
}
TEST(Lu, FirstSingularityAndContinuation) {
    V a{0,0,1,2}; const V original=a; auto r=sce::factor_lu(a,2,2,2);
    EXPECT_EQ(r.info,1); residual(original,a,2,2,2,r);
    V b{1,2,2,4}; auto s=sce::factor_lu(b,2,2,2); EXPECT_EQ(s.info,2);
}
TEST(Lu, TinyPivotAvoidsReciprocalOverflow) {
    const double t=std::numeric_limits<double>::denorm_min()*16;
    ASSERT_GT(t,0); V a{t,2*t,t}; auto r=sce::factor_lu(a,3,1,3);
    EXPECT_EQ(r.info,0); EXPECT_EQ(r.pivots[0],1);
    EXPECT_DOUBLE_EQ(a[1],.5); EXPECT_DOUBLE_EQ(a[2],.5);
    for(double x:a) EXPECT_TRUE(std::isfinite(x));
}
TEST(Lu, SafeMinimumBoundary) {
    const double t=std::numeric_limits<double>::min();
    V a{t,t/2}; auto r=sce::factor_lu(a,2,1,2);
    EXPECT_EQ(r.info,0); EXPECT_DOUBLE_EQ(a[1],.5);
}
TEST(Lu, LeadingDimensionPaddingPreserved) {
    V a{2,4,99,99,1,3,99,99,5,7,99,99}; const V original=a;
    auto r=sce::factor_lu(a,2,3,4); residual(original,a,2,3,4,r);
    for(auto p:{2,3,6,7,10,11}) EXPECT_EQ(a[p],99);
}
TEST(Lu, CheckedBoundaryRejectsInvalidStorageAndNonfiniteInput) {
    V a(4,1); EXPECT_THROW(sce::factor_lu(a,-1,2,2),std::invalid_argument);
    EXPECT_THROW(sce::factor_lu(a,2,-1,2),std::invalid_argument);
    EXPECT_THROW(sce::factor_lu(a,2,2,1),std::invalid_argument);
    EXPECT_THROW(sce::factor_lu(a,3,2,3),std::invalid_argument);
    a[0]=std::numeric_limits<double>::infinity();
    EXPECT_THROW(sce::factor_lu(a,2,2,2),std::invalid_argument);
    a[0]=std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(sce::factor_lu(a,2,2,2),std::invalid_argument);
}
TEST(Lu, SeededRectangularResiduals) {
    std::mt19937_64 generator(0x4c55474554524632ULL);
    std::uniform_real_distribution<double> distribution(-1,1);
    for(I m=1;m<=17;++m) for(I n=1;n<=19;++n) {
        const I ld=m+2; V a(static_cast<std::size_t>(ld)*n,123);
        for(I col=0;col<n;++col) for(I row=0;row<m;++row)
            a[static_cast<std::size_t>(col)*ld+row]=distribution(generator);
        const V original=a; auto r=sce::factor_lu(a,m,n,ld);
        EXPECT_EQ(r.info,0); residual(original,a,m,n,ld,r);
    }
}
TEST(Lu, HilbertResidual) {
    constexpr I n=8; V a(n*n);
    for(I col=0;col<n;++col) for(I row=0;row<n;++row) a[col*n+row]=1.0/(1+row+col);
    const V original=a; auto r=sce::factor_lu(a,n,n,n);
    EXPECT_EQ(r.info,0); residual(original,a,n,n,n,r);
}
}

// Audit coverage: singular status must survive recursive panel offsets.
TEST(LuAudit, RecursiveSingularityAndContinuation) {
    V a{2,0,0,0, 0,3,0,0, 0,0,0,0, 0,0,0,5};
    const V original=a; auto r=sce::factor_lu(a,4,4,4);
    EXPECT_EQ(r.info,3); EXPECT_EQ(r.pivots,(std::vector<I>{0,1,2,3}));
    EXPECT_EQ(a[15],5); residual(original,a,4,4,4,r);
}
TEST(LuAudit, MinimalSpanAndNonfinitePadding) {
    const double nan=std::numeric_limits<double>::quiet_NaN();
    V a{1,3,nan,nan,2,4}; const V original=a;
    auto r=sce::factor_lu(a,2,2,4);
    EXPECT_EQ(r.info,0); residual(original,a,2,2,4,r);
    EXPECT_TRUE(std::isnan(a[2])); EXPECT_TRUE(std::isnan(a[3]));
    V invalid{1,3,2,nan}; const V before=invalid;
    EXPECT_THROW(sce::factor_lu(invalid,2,2,2),std::invalid_argument);
    for(int i=0;i<3;++i) EXPECT_EQ(invalid[i],before[i]);
    EXPECT_TRUE(std::isnan(invalid[3]));
}
TEST(LuAudit, WilkinsonGrowthIsNotAStabilityGuarantee) {
    constexpr I n=20; V a(n*n,0);
    for(I col=0;col<n;++col) for(I row=0;row<n;++row)
        a[col*n+row] = col==n-1 ? 1.0 : (row==col ? 1.0 : (row>col ? -1.0 : 0.0));
    const V original=a; auto r=sce::factor_lu(a,n,n,n);
    EXPECT_EQ(r.info,0);
    EXPECT_DOUBLE_EQ(a[n*n-1],std::ldexp(1.0,n-1));
    residual(original,a,n,n,n,r);
}
