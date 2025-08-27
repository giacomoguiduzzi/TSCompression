//
// Created by giacomo on 20/06/25.
//

#include "ts_compression.h"

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <size_ref_MB> <size_block>" << std::endl;
        return 1;
    }
    std::string input = argv[1];
    auto size_ref = static_cast<uint64_t>(strtol(
        argv[2], nullptr, 10
    ));
    auto size_block = static_cast<uint64_t>(strtol(argv[3], nullptr, 10));
    auto size = util::file::file_size(input) / sizeof(uint64_t);
    std::vector<int64_t> time_series(size);
    std::ifstream in(input);
    in.read(reinterpret_cast<char *>(time_series.data()), static_cast<long>(time_series.size() * sizeof(int64_t)));
    in.close();

    // Add 1 to the time series if it contains a 0, to avoid issues with RLZ compression
    if (std::ranges::find(time_series.begin(), time_series.end(), 0) != time_series.end()) {
        std::ranges::transform(time_series.begin(), time_series.end(), time_series.begin(),
                               [](const int64_t value) { return value + 1; });
    }

    NaiveRLZ rlz(time_series, size_ref, size_block, 0, true);
    rlz.init_factorization(&time_series);

    std::vector<NaiveRLZ::factor_type> factors;
    uint64_t max_factors_length = 0;

    while (rlz.has_next()) {
        auto factor = rlz.next();

        if (factor.length > max_factors_length) {
            max_factors_length = factor.length;
        }

        factors.push_back(factor);
    }

    std::vector<int64_t> result;
    rlz.decompress(factors, result);
    if (result.size() != time_series.size()) {
        std::cout << "Error size" << std::endl;
        exit(1);
    }
    for (uint64_t i = 0; i < result.size(); ++i) {
        if (result[i] != time_series[i]) {
            std::cout << "Error in time series at position: " << i << std::endl;
            exit(1);
        }
    }

    const std::string reference_filename = std::to_string(getpid()) + ".rev_ref";
    std::remove(reference_filename.c_str());

    std::cout << max_factors_length << std::flush;

    return 0;
}
