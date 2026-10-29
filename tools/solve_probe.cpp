#include "sce/solve.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
int main() {
  std::cout << std::setprecision(17);
  int n, r;
  while (std::cin >> n >> r) {
    if (n < 0 || n > 256 || r < 0 || r > 16)
      return 2;
    std::vector<double> a(static_cast<std::size_t>(n) * n),
        b(static_cast<std::size_t>(n) * r);
    for (double &v : a)
      if (!(std::cin >> v))
        return 3;
    for (double &v : b)
      if (!(std::cin >> v))
        return 3;
    try {
      auto f = sce::factor_lu(a, n, n, std::max(1, n));
      auto x = sce::solve_lu(a, n, std::max(1, n), f, b, r, std::max(1, n));
      std::cout << "{\"x\":[";
      for (std::size_t i = 0; i < x.size(); ++i) {
        if (i)
          std::cout << ',';
        std::cout << x[i];
      }
      std::cout << "]}\n";
    } catch (const std::domain_error &) {
      std::cout << "{\"singular\":true}\n";
    }
  }
}
