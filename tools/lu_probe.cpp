#include "sce/lu.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// Local validation/benchmark adapter. This is new application code, not recovered.
int main() {
    std::cout << std::setprecision(17);
    int n, repetitions, packed_output;
    while (std::cin >> n >> repetitions >> packed_output) {
        if (n < 0 || n > 1024 || repetitions < 1 || repetitions > 100) return 2;
        std::vector<double> original(static_cast<std::size_t>(n)*n);
        for (double& x : original) if (!(std::cin >> x)) return 3;
        std::vector<double> packed;
        sce::LuResult result;
        std::vector<long long> times;
        int sign = 1;
        double logarithm = 0;
        for (int trial = 0; trial <= repetitions; ++trial) {
            const auto start = std::chrono::steady_clock::now();
            packed = original;
            result = sce::factor_lu(packed,n,n,std::max(1,n));
            sign = result.info == 0 ? 1 : 0;
            logarithm = 0;
            if (sign != 0) for (int j = 0; j < n; ++j) {
                const double diagonal = packed[static_cast<std::size_t>(j)*n+j];
                if (result.pivots[j] != j) sign = -sign;
                if (std::signbit(diagonal)) sign = -sign;
                logarithm += std::log(std::abs(diagonal));
            }
            const auto end = std::chrono::steady_clock::now();
            if (trial != 0) times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count());
        }
        const auto samples = times; // Preserve acquisition order for reproducibility.
        std::sort(times.begin(),times.end());
        std::cout << "{\"n\":" << n << ",\"info\":" << result.info << ",\"sign\":" << sign
                  << ",\"logabsdet\":";
        if (sign == 0) std::cout << "null"; else std::cout << logarithm;
        std::cout << ",\"median_ns\":" << times[times.size()/2];
        std::cout << ",\"samples_ns\":[";
        for (std::size_t j=0;j<samples.size();++j) { if (j) std::cout << ','; std::cout << samples[j]; }
        std::cout << ']';
        if (packed_output) {
            std::cout << ",\"pivots\":[";
            for (std::size_t j=0;j<result.pivots.size();++j) {
                if (j) std::cout << ','; std::cout << result.pivots[j];
            }
            std::cout << "],\"packed\":[";
            for (std::size_t j=0;j<packed.size();++j) { if (j) std::cout << ','; std::cout << packed[j]; }
            std::cout << ']';
        }
        std::cout << "}\n";
    }
}
