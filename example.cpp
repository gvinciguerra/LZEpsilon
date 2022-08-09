#include <iostream>
#include <random>
#include "LZEpsilon.hpp"

int main() {
    std::mt19937 engine;
    std::geometric_distribution<uint32_t> distribution(0.85);
    sdsl::int_vector<32> data(1000000);
    for (auto i = 1; i < data.size(); ++i)
        data[i] = data[i - 1] + distribution(engine) + 1;

    LZEpsilon lze(data, 3);

    std::cout << "Bits per integer:      " << 8. * lze.size_in_bytes() / data.size() << std::endl
              << "Phrases count:         " << lze.phrases_count() << std::endl
              << "# of elements <= 500:  " << lze.rank(500) << std::endl
              << "10th smallest element: " << lze.select(10) << std::endl;

    return 0;
}