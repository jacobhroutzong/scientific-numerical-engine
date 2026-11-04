#include "sce/least_squares.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
int main() {
  std::cout << std::setprecision(17);
  int m, n, r, repetitions, estimate;
  double tolerance;
  while (std::cin >> m >> n >> r >> tolerance >> estimate >> repetitions) {
    if (m < 1 || m > 256 || n < 1 || n > 256 || r < 1 || r > 16 ||
        repetitions < 1 || repetitions > 10000)
      return 2;
    std::vector<double> a(static_cast<std::size_t>(m) * n),
        b(static_cast<std::size_t>(m) * r);
    for (auto &x : a)
      if (!(std::cin >> x))
        return 3;
    for (auto &x : b)
      if (!(std::cin >> x))
        return 3;
    try {
      sce::LeastSquaresOptions options;
      if (tolerance != -1)
        options.relative_rank_tolerance = tolerance;
      options.estimate_r_condition = estimate != 0;
      const auto start = std::chrono::steady_clock::now();
      sce::LeastSquaresResult result;
      for (int i = 0; i < repetitions; ++i)
        result = sce::least_squares(a, m, n, m, b, r, m, options);
      const auto elapsed = std::chrono::duration<double, std::nano>(
                               std::chrono::steady_clock::now() - start)
                               .count() /
                           repetitions;
      std::cout << "{\"status\":\"ok\",\"x\":[";
      for (std::size_t i = 0; i < result.x.size(); ++i) {
        if (i)
          std::cout << ',';
        std::cout << result.x[i];
      }
      std::cout << "],\"diagnostics\":[";
      for (std::size_t i = 0; i < result.diagnostics.size(); ++i) {
        if (i)
          std::cout << ',';
        const auto &d = result.diagnostics[i];
        std::cout << "[" << d.residual_norm2 << ',' << d.scaled_residual << ','
                  << d.normalized_stationarity << ']';
      }
      std::cout << "],\"min_relative_diagonal\":"
                << result.minimum_relative_r_diagonal
                << ",\"rank_tolerance\":" << result.relative_rank_tolerance
                << ",\"condition_inf_r\":";
      if (result.r_condition_inf)
        std::cout << *result.r_condition_inf;
      else
        std::cout << "null";
      std::cout << ",\"nanoseconds_per_call\":" << elapsed << "}\n";
    } catch (const sce::RankDeficient &) {
      std::cout << "{\"status\":\"rank_rejected\"}\n";
    } catch (const std::overflow_error &) {
      std::cout << "{\"status\":\"overflow\"}\n";
    } catch (const std::invalid_argument &) {
      std::cout << "{\"status\":\"invalid\"}\n";
    }
  }
}
