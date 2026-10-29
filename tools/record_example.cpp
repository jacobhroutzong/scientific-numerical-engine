#include "sce/solve.hpp"
#include <iomanip>
#include <iostream>
int main() {
  std::vector<double> a{1, 3, 2, 4}, b{5, 11, 5, 9};
  auto f = sce::factor_lu(a, 2, 2, 2);
  auto x = sce::solve_lu(a, 2, 2, f, b, 2, 2);
  std::cout << std::setprecision(17)
            << "{\"kind\":\"recorded-cpp-fixture\",\"date\":\"2026-10-26\","
               "\"order\":2,\"storage\":\"column-major\",\"a\":[1,3,2,4],\"b\":"
               "[5,11,5,9],\"packed\":[";
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (i)
      std::cout << ',';
    std::cout << a[i];
  }
  std::cout << "],\"pivots\":[" << f.pivots[0] << ',' << f.pivots[1]
            << "],\"info\":" << f.info << ",\"x\":[";
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (i)
      std::cout << ',';
    std::cout << x[i];
  }
  std::cout << "]}\n";
}
