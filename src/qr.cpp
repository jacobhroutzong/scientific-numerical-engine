#include "sce/qr.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace sce {
namespace {
double finite(double value) {
  if (!std::isfinite(value))
    throw std::overflow_error("QR arithmetic overflow");
  return value;
}
} // namespace
QrResult factor_qr(std::span<const double> input, index_type m, index_type n,
                   index_type ld) {
  if (n < 1 || m < n || m > 256 || ld < m)
    throw std::invalid_argument(
        "QR requires 1 <= columns <= rows <= 256 and ld >= rows");
  const auto required = static_cast<std::size_t>(n - 1) * ld + m;
  if (input.size() < required)
    throw std::invalid_argument("QR storage is too small");
  std::vector<double> a(static_cast<std::size_t>(m) * n), tau(n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i) {
      const double value = input[static_cast<std::size_t>(j) * ld + i];
      if (!std::isfinite(value))
        throw std::invalid_argument("QR requires finite active entries");
      a[j * m + i] = value;
    }
  for (int k = 0; k < n; ++k) {
    double tail = 0;
    for (int i = k + 1; i < m; ++i)
      tail = std::hypot(tail, a[k * m + i]);
    const double alpha = a[k * m + k];
    if (tail == 0)
      continue;
    const double beta = -std::copysign(finite(std::hypot(alpha, tail)), alpha);
    tau[k] = 1 - alpha / beta;
    // Scaling avoids alpha-beta overflow and a tiny reciprocal overflow.
    const double scale = std::max(std::abs(alpha), tail);
    const double denominator = alpha / scale - beta / scale;
    for (int i = k + 1; i < m; ++i)
      a[k * m + i] = (a[k * m + i] / scale) / denominator;
    a[k * m + k] = beta;
    for (int j = k + 1; j < n; ++j) {
      double dot = a[j * m + k];
      for (int i = k + 1; i < m; ++i)
        dot = finite(dot + a[k * m + i] * a[j * m + i]);
      dot = finite(tau[k] * dot);
      a[j * m + k] = finite(a[j * m + k] - dot);
      for (int i = k + 1; i < m; ++i)
        a[j * m + i] = finite(a[j * m + i] - a[k * m + i] * dot);
    }
  }
  QrResult result{m, n, std::vector<double>(static_cast<std::size_t>(m) * n),
                  std::vector<double>(static_cast<std::size_t>(n) * n)};
  for (int j = 0; j < n; ++j) {
    result.q[j * m + j] = 1;
    for (int i = 0; i <= j; ++i)
      result.r[j * n + i] = a[j * m + i];
  }
  for (int k = n - 1; k >= 0; --k) {
    if (tau[k] == 0)
      continue;
    for (int j = 0; j < n; ++j) {
      double dot = result.q[j * m + k];
      for (int i = k + 1; i < m; ++i)
        dot += a[k * m + i] * result.q[j * m + i];
      dot *= tau[k];
      result.q[j * m + k] -= dot;
      for (int i = k + 1; i < m; ++i)
        result.q[j * m + i] -= a[k * m + i] * dot;
    }
  }
  return result;
}
} // namespace sce
