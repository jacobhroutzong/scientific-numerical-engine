#include "sce/lu.hpp"
#include "sce/solve.hpp"
#include <algorithm>
#include <cmath>
#include <emscripten/bind.h>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace {
std::string quote(const std::string &s) {
  std::string out = "\"";
  for (unsigned char c : s) {
    if (c == '"' || c == '\\')
      out += '\\';
    if (c >= 32)
      out += static_cast<char>(c);
  }
  return out + '"';
}
std::string solve(const std::string &wire) {
  try {
    if (wire.size() > 32768)
      throw std::invalid_argument("Input exceeds the browser request limit");
    std::istringstream in(wire);
    int n = 0, r = 0;
    if (!(in >> n >> r) || n < 1 || n > 16 || r < 1 || r > 4)
      throw std::invalid_argument(
          "Use an order of 1–16 and 1–4 right-hand sides");
    std::vector<double> a(n * n), b(n * r);
    auto read = [&](std::vector<double> &v, int columns) {
      for (int row = 0; row < n; ++row)
        for (int col = 0; col < columns; ++col) {
          double x;
          if (!(in >> x) || !std::isfinite(x) || std::abs(x) > 1e100)
            throw std::invalid_argument(
                "Use finite numbers with magnitude at most 1e100");
          v[col * n + row] = x;
        }
    };
    read(a, n);
    read(b, r);
    std::string extra;
    if (in >> extra)
      throw std::invalid_argument("Unexpected extra matrix values");
    auto packed = a;
    auto f = sce::factor_lu(packed, n, n, n);
    if (f.info)
      throw std::domain_error(
          "Exactly singular matrix: zero U diagonal at position " +
          std::to_string(f.info));
    auto x = sce::solve_lu(packed, n, n, f, b, r, n);
    double anorm = 0, xnorm = 0, bnorm = 0, error = 0;
    for (int i = 0; i < n; ++i) {
      double ar = 0, xr = 0, br = 0, er = 0;
      for (int k = 0; k < n; ++k)
        ar += std::abs(a[k * n + i]);
      for (int j = 0; j < r; ++j) {
        xr += std::abs(x[j * n + i]);
        br += std::abs(b[j * n + i]);
        double ax = 0;
        for (int k = 0; k < n; ++k)
          ax += a[k * n + i] * x[j * n + k];
        er += std::abs(ax - b[j * n + i]);
      }
      anorm = std::max(anorm, ar);
      xnorm = std::max(xnorm, xr);
      bnorm = std::max(bnorm, br);
      error = std::max(error, er);
    }
    const double scale = anorm * xnorm + bnorm;
    const double backward = scale == 0 ? 0 : error / scale;
    if (!std::isfinite(error) || !std::isfinite(scale) ||
        !std::isfinite(backward))
      throw std::overflow_error("Residual overflow; rescale the matrix values");
    std::ostringstream out;
    out << std::setprecision(17) << "{\"ok\":true,\"x\":[";
    for (int i = 0; i < n; ++i) {
      if (i)
        out << ',';
      out << '[';
      for (int j = 0; j < r; ++j) {
        if (j)
          out << ',';
        out << x[j * n + i];
      }
      out << ']';
    }
    out << "],\"pivots\":[";
    for (int i = 0; i < n; ++i) {
      if (i)
        out << ',';
      out << f.pivots[i];
    }
    out << "],\"backwardError\":" << backward << ",\"info\":0}";
    return out.str();
  } catch (const std::exception &e) {
    return "{\"ok\":false,\"error\":" + quote(e.what()) + "}";
  }
}
} // namespace
EMSCRIPTEN_BINDINGS(numerical_engine) { emscripten::function("solve", &solve); }
