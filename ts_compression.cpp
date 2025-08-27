#include "ts_compression.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <filesystem>
#include <stdexcept>

// Helper function to read a single byte
uint8_t read_u8(std::ifstream &reader) {
    uint8_t value;
    if (!reader.read(reinterpret_cast<char *>(&value), sizeof(value))) {
        throw std::runtime_error("Failed to read uint8_t");
    }
    return value;
}

// Helper function to read an integer of a given size
int64_t read_int(std::ifstream &reader, const int64_t size) {
    std::vector<uint8_t> buffer(size);
    if (!reader.read(reinterpret_cast<char *>(buffer.data()), size)) {
        throw std::runtime_error("Failed to read integer");
    }
    int64_t result = 0;
    for (size_t i = 0; i < size; ++i) {
        result |= static_cast<int64_t>(buffer[i]) << (8 * i);
    }
    return result;
}

// Helper function to read a signed integer of a given size
int64_t read_signed_int(std::ifstream &reader, const int64_t size) {
    std::vector<uint8_t> buffer(size);
    if (!reader.read(reinterpret_cast<char *>(buffer.data()), size)) {
        throw std::runtime_error("Failed to read signed integer");
    }
    int64_t result = 0;
    for (size_t i = 0; i < size; ++i) {
        result |= static_cast<int64_t>(buffer[i]) << (8 * i);
    }
    // Sign extension for negative numbers
    if (buffer[size - 1] & 0x80) {
        for (size_t i = size; i < sizeof(int64_t); ++i) {
            result |= static_cast<int64_t>(0xFF) << (8 * i);
        }
    }
    return result;
}

/**
 * @brief Reconstructs reference data and timestamps from a binary file.
 *
 * This function reads a binary file containing compressed time series data and reconstructs
 * the reference, timestamps, and other metadata. It allocates memory for the output arrays
 * and populates them with the data read from the file. The time series and reference values are **normalized**, thus
 * to retrieve the real values it is necessary to perform the classical min-max inverse operation,
 * once ts_max and ts_min have been multiplied by @code 10^value_multiplier@endcode and casted to double:
 *
 * @code
 *      double real_ts = new double[ts_size];
 *      for (auto counter = 0; counter < ts_size; counter++) {
 *           auto real_value = static_cast<double>(ts[counter] * math.pow(10, value_multiplier);
 *           real_ts[counter] = real_value * (ts_max - ts_min) + ts_min;
 *      }
 * @endcode
 *
 * @param file_path The path to the binary file to be read.
 * @param timestamps A reference to a pointer that will be allocated and filled with the timestamps.
 * @param reference A reference to a pointer that will be allocated and filled with the reference data.
 * @param ts A reference to a pointer that will be allocated and filled with the time series data.
 * @param ts_size A reference to a variable that will hold the size of the time series data.
 * @param reference_size A reference to a variable that will hold the size of the reference data.
 * @param ts_min A reference to a variable that will hold the minimum value in the time series.
 * @param ts_max A reference to a variable that will hold the maximum value in the time series.
 * @param min_decimals A reference to a variable that will hold the number of decimals for the minimum value.
 * @param max_decimals A reference to a variable that will hold the number of decimals for the maximum value.
 * @param value_multiplier A reference to a variable that will hold the multiplier used for scaling values.
 *
 * @throws std::runtime_error If the file cannot be opened or if there is an error reading the file.
 */
void reconstruct_reference_from_binary(const std::string &file_path, std::vector<uint64_t> &timestamps,
                                       std::vector<int64_t> &reference, std::vector<int64_t> &ts,
                                       uint64_t &ts_size, uint64_t &reference_size,
                                       int64_t &ts_min, int64_t &ts_max, uint64_t &min_decimals,
                                       uint64_t &max_decimals, uint64_t &value_multiplier) {
    std::ifstream reader(file_path, std::ios::binary);
    if (!reader) {
        throw std::runtime_error("Failed to open input file.");
    }

    auto integer_ts = static_cast<bool>(read_u8(reader));

    // Read the number of decimals for min and max
    min_decimals = read_u8(reader);
    max_decimals = read_u8(reader);

    // Read subsequence size information
    int64_t subsequence_size_bytes_size = read_u8(reader);
    int64_t subsequence_size = read_int(reader, subsequence_size_bytes_size);

    // Read time series length information
    int64_t ts_length_size_size = read_u8(reader);
    int64_t ts_length = read_int(reader, ts_length_size_size);

    // Read reference length information
    int64_t reference_length_size_size = read_u8(reader);
    int64_t reference_length = read_int(reader, reference_length_size_size);

    // Read min and max values
    // int64_t ts_min_scaled = read_signed_int(reader, subsequence_size);
    // int64_t ts_max_scaled = read_signed_int(reader, subsequence_size);
    ts_min = read_signed_int(reader, subsequence_size);
    ts_max = read_signed_int(reader, subsequence_size);

    // ts_min = static_cast<double>(ts_min_scaled) / std::pow(10.0, min_decimals);
    // ts_max = static_cast<double>(ts_max_scaled) / std::pow(10.0, max_decimals);

    // Read timestamp size information
    int64_t timestamp_size_bytes_size = read_u8(reader);
    int64_t timestamp_size = read_int(reader, timestamp_size_bytes_size);

    // Allocate memory for timestamps array
    // timestamps = new int64_t[ts_length];

    // Read timestamps
    timestamps.resize(ts_length);
    for (size_t i = 0; i < ts_length; ++i) {
        timestamps[i] = read_signed_int(reader, timestamp_size);
    }

    // Read the subsequence multiplier for references and ts
    value_multiplier = read_u8(reader);

    // Allocate memory for reference array
    // reference = new int64_t[reference_length];

    // Read reference
    reference.resize(reference_length);
    for (size_t i = 0; i < reference_length; ++i) {
        int64_t value_scaled = read_signed_int(reader, subsequence_size);
        // reference[i] = static_cast<double>(value_scaled) / std::pow(10.0, value_multiplier);
        // if(integer_ts){
        //     value_scaled = static_cast<int64_t>(static_cast<double>(value_scaled) / std::pow(10.0, value_multiplier));
        // }
        reference[i] = value_scaled;
    }

    // Allocate memory for time series array
    // ts = new int64_t[ts_length];

    // Read time series data
    ts.resize(ts_length);
    for (size_t i = 0; i < ts_length; ++i) {
        int64_t ts_value = read_signed_int(reader, subsequence_size);
        // ts[i] = static_cast<double>(ts_value) / std::pow(10.0, integer_ts ? 0 : value_multiplier);
        ts[i] = ts_value;
    }

    ts_size = ts_length;
    reference_size = reference_length;
}

void search_for_subsequences(const std::vector<int64_t> &reference, const std::vector<int64_t> &ts,
                             const int64_t ts_size, const int64_t reference_size) {
    typedef struct {
        std::vector<int64_t> subseq;
        uint64_t ts_index = 0;
    } subsequence_t;

    std::vector<subsequence_t> subsequences;

    // for each value in the time series
    auto ts_index = -1;
    while (ts_index < ts_size) {
        ts_index++;
        if (ts_index % 10000 == 0) {
            std::cout << "Progress: index " << ts_index << " / " << ts_size << std::endl;
        }
        auto go_to_next_ts_value = false;
        // for each value in the reference
        auto ref_index = -1;

        while (ref_index < reference_size) {
            if (go_to_next_ts_value) {
                break;
            }

            ref_index++;

            auto ts_value = ts[ts_index];

            // if the value is the same as the reference, start saving the subsequence
            if (const auto ref_value = reference[ref_index]; ts_value == ref_value) {
                // std::cout << "Found subsequence starting at index " << ts_index << std::endl;
                subsequence_t subsequence;
                subsequence.subseq.push_back(ts_value);
                subsequence.ts_index = ts_index;

                auto ts_value_offset = 1;
                if (ts_index + ts_value_offset < ts_size) {
                    auto next_ts_value = ts[ts_index + ts_value_offset];
                    // while the time series does not change, it keeps being equal to the reference.
                    // add each point to the subsequence until there is a difference
                    while (next_ts_value == ref_value) {
                        if (ts_index + ts_value_offset < ts_size) {
                            subsequence.subseq.push_back(next_ts_value);
                            ts_value_offset++;
                            next_ts_value = ts[ts_index + ts_value_offset];
                        }
                    }
                } else // if the time series is at the end, then the next value is not available
                    break;

                // a difference was found, so start looking from the next value in the reference
                auto subseq_index = ref_index;
                while (subseq_index < reference_size) {
                    subseq_index++;

                    if (subseq_index < reference_size) {
                        auto next_ref_value = reference[subseq_index];
                        // if the next value in the reference is equal to the time series, add it to the subsequence
                        if (auto next_ts_value = ts[ts_index + ts_value_offset];
                            ts_index + ts_value_offset < ts_size && next_ts_value == next_ref_value
                        ) {
                            subsequence.subseq.push_back(next_ref_value);

                            ts_value_offset++;
                            next_ts_value = ts[ts_index + ts_value_offset];

                            // while the time series does not change, it keeps being equal to the reference.
                            while (next_ts_value == next_ref_value) {
                                if (ts_index + ts_value_offset < ts_size) {
                                    subsequence.subseq.push_back(next_ts_value);
                                    ts_value_offset++;
                                    next_ts_value = ts[ts_index + ts_value_offset];
                                }
                            }
                        } else {
                            // if the next value in the reference is not equal to the time series, break the loop.
                            // the subsequence is complete, so the most external for loop can go on with the next
                            // value in the time series
                            go_to_next_ts_value = true;
                            ts_index += ts_value_offset - 1;
                            break;
                        }
                    }
                }

                subsequences.push_back(subsequence);
                // std::cout << "Total subsequence length: " << subsequence.subseq.size() << std::endl;
            }
        }
    }

    std::cout << "Found " << subsequences.size() << " subsequences." << std::endl;
    auto length_one_subsequences = 0, bigger_subsequences = 0;
    for (const auto &[subseq, ts_index]: subsequences) {
        // std::cout << "Subsequence of length " << subseq.size() << " at index: " << ts_index << std::endl;
        if (subseq.size() == 1) {
            length_one_subsequences++;
        } else {
            bigger_subsequences++;
        }
    }
    std::cout << "Found " << length_one_subsequences << " subsequences of length 1." << std::endl;
    std::cout << "Found " << bigger_subsequences << " subsequences of length > 1." << std::endl;
}

void search_for_duplicates(const std::vector<int64_t> &reference) {
    std::set<int64_t> unique_values;
    std::set<int64_t> duplicates;

    for (const auto &value: reference) {
        if (unique_values.contains(value)) {
            duplicates.insert(value);
        } else {
            unique_values.insert(value);
        }
    }

    if (!duplicates.empty()) {
        std::cout << "Duplicates found in the reference:" << std::endl;
        for (const auto &duplicate: duplicates) {
            std::cout << duplicate << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "No duplicates found in the reference." << std::endl;
    }
}

int main(int argc, char *argv[]) {
    try {
        const std::string binary_filename = argc < 2 ? "./data/1week-10190.bin" : argv[1];

        auto timestamps = std::vector<uint64_t>();
        auto reference = std::vector<int64_t>();
        auto ts = std::vector<int64_t>();
        int64_t ts_min, ts_max;
        uint64_t ts_size, reference_size, min_decimals, max_decimals, value_multiplier;


        reconstruct_reference_from_binary(binary_filename, timestamps, reference, ts, ts_size,
                                          reference_size, ts_min, ts_max,
                                          min_decimals, max_decimals, value_multiplier);

        std::cout << "Read " << ts_size << " timestamps and values for a series "
                << "with minimum value " << ts_min << " and maximum value " << ts_max
                << ". The reference length is " << reference_size << std::endl;

        // check for duplicates in reference
        // search_for_duplicates(reference);
        // check for the reference validity
        // search_for_subsequences(reference, ts, static_cast<int64_t>(ts_size), static_cast<int64_t>(reference_size));

        // Invert normalization, if needed
        // const auto *denormalized_reference = new double[reference_size];
        // std::transform(reference[0], reference[reference_size - 1], denormalized_reference[0],
        //                [ts_min, ts_max](const double x) {
        //                    return x * (ts_max - ts_min) + ts_min;
        //                });

        // before saving the values as bytes, it is necessary to multiply each value by the value_multiplier,
        // so that the values are stored as integers
        // auto integer_reference = std::vector<int64_t>(reference_size);
        //
        // // Convert references to bytes
        // const uint64_t bytes_reference_size = reference_size * sizeof(int64_t);
        // auto reference_bytes = std::vector<uint8_t>(bytes_reference_size);
        // constexpr size_t num_bytes_per_value = sizeof(int64_t) / sizeof(uint8_t);
        //
        // for (size_t reference_counter = 0; reference_counter < reference_size; reference_counter++) {
        //     integer_reference[reference_counter] =
        //             reference[reference_counter] * static_cast<int64_t>(std::pow(10.0, value_multiplier));
        //
        //     // handle zero values as they arean issue for the creation of the compressed file through RLZ
        //     if (integer_reference[reference_counter] == 0) {
        //         integer_reference[reference_counter] = -1;
        //     }
        //     // const auto ref_bytes = reinterpret_cast<const uint8_t *>(&denormalized_reference[reference_counter]);
        //     const auto ref_bytes = reinterpret_cast<uint8_t *>(&integer_reference[reference_counter]);
        //
        //     for (size_t byte_idx = 0; byte_idx < num_bytes_per_value; byte_idx++) {
        //         reference_bytes[reference_counter * sizeof(int64_t) + byte_idx] = ref_bytes[byte_idx];
        //     }
        // }

        // // the same applies for the time series
        // auto integer_ts = std::vector<int64_t>(ts_size);
        //
        // // Convert time series to bytes
        // const uint64_t bytes_ts_size = ts_size * sizeof(int64_t);
        // auto ts_bytes = std::vector<uint8_t>(bytes_ts_size);
        // constexpr size_t num_bytes_per_value_ts = sizeof(int64_t) / sizeof(uint8_t);
        // for (size_t ts_counter = 0; ts_counter < ts_size; ts_counter++) {
        //     integer_ts[ts_counter] = ts[ts_counter] * static_cast<int64_t>(std::pow(10.0, value_multiplier));
        //
        //     // handle zero values as they arean issue for the creation of the compressed file through RLZ
        //     if (integer_ts[ts_counter] == 0) {
        //         integer_ts[ts_counter] = -1;
        //     }
        //
        //     // const auto ref_bytes = reinterpret_cast<const uint8_t *>(&ts[ts_counter]);
        //     const auto ref_bytes = reinterpret_cast<uint8_t *>(&integer_ts[ts_counter]);
        //
        //     for (size_t byte_idx = 0; byte_idx < num_bytes_per_value_ts; byte_idx++) {
        //         ts_bytes[ts_counter * sizeof(int64_t) + byte_idx] = ref_bytes[byte_idx];
        //     }
        // }

        // before creating the RLZ object, we need to check if the 0 value is contained in the reference and
        // time series object.
        // RLZ needs to use the 0 value to create an efficient compressed version of the original data, thus it needs
        // to be replaced with another value that is easily identifiable.
        // int64_t replacement_value = std::numeric_limits<int64_t>::max();
        auto had_zeroes = false;
        if (std::ranges::find(reference.begin(), reference.end(), 0) != reference.end()) {
            // std::cout << "Reference contains 0 value, replacing it with the maximum representable value" << std::endl;
            std::cout << "Reference contains 0 value, adding 1 to reference and time series." << std::endl;
            had_zeroes = true;
            // std::ranges::replace(reference.begin(), reference.end(), 0, replacement_value);
            // std::ranges::replace(ts.begin(), ts.end(), 0, replacement_value);
            std::ranges::transform(reference.begin(), reference.end(), reference.begin(),
                                   [](const int64_t value) { return value + 1; });
            std::ranges::transform(ts.begin(), ts.end(), ts.begin(),
                                   [](const int64_t value) { return value + 1; });
        }

        // BinaryRLZ rlz(reference_bytes, bytes_reference_size, sizeof(uint8_t));
        // NaiveRLZ rlz(integer_reference, integer_reference.size());
        CustomRLZ rlz(reference, reference_size, 0, 0.0, false);
        // este usa una block size de tamaño ts_size al lugar de 1MB.
        // NaiveRLZ rlz_naive(ts, ts_size, ts_size, 0.0, true);
        // este nos da una compresion naive con una referencia hecha con bloques random de 1MB

        uintmax_t compressed_custom_file_size;
        try {
            compressed_custom_file_size = std::filesystem::file_size(std::to_string(getpid()) + ".rev_ref");
        } catch (const std::filesystem::filesystem_error &e) {
            std::cerr << "Error retrieving the file system size of the original and compressed files: " << e.what() <<
                    std::endl;

            return 4;
        }

        // TODO: reference_size probably needs to be in bytes, so change reference_size with
        //  reference_size * sizeof(uint64_t) and ts_size with ts_size * sizeof(uint64_t)
        NaiveRLZ rlz_naive(
            ts,
            std::min(reference_size * sizeof(int64_t), ts_size * sizeof(int64_t)),
            1024,
            0.0,
            true
        );
        // implicit or lazy encoding
        // rlz.init_factorization(&ts_bytes);
        // rlz.init_factorization(&integer_ts);
        rlz.init_factorization(&ts);
        rlz_naive.init_factorization(&ts);

        // search_for_duplicates(rlz.reference);
        // std::vector<int64_t> rlz_naive_reference(rlz.reference.size());
        // std::ranges::copy(rlz.reference, rlz_naive_reference.begin());
        // search_for_duplicates(rlz_naive_reference);

        // decode
        std::vector<CustomRLZ::factor_type> custom_factors;
        // factor has offset and length. we need to save the max length to compute the size of the reference
        uint64_t max_factors_length = 0, max_naive_factors_length = 0;

        while (rlz.has_next()) {
            auto factor = rlz.next();

            if (factor.length > max_factors_length) {
                max_factors_length = factor.length;
            }

            custom_factors.push_back(factor);
        }

        std::vector<NaiveRLZ::factor_type> naive_factors;
        while (rlz_naive.has_next()) {
            auto factor = rlz_naive.next();

            if (factor.length > max_naive_factors_length) {
                max_naive_factors_length = factor.length;
            }

            naive_factors.push_back(factor);
        }

        std::vector<int64_t> recovered_values;
        std::cout << "Decompressing " << std::endl << std::flush;
        rlz.decompress(custom_factors, recovered_values);

        std::cout << "Decompressing naive" << std::endl << std::flush;
        std::vector<int64_t> naive_recovered_values;
        rlz_naive.decompress(naive_factors, naive_recovered_values);

        if (recovered_values.size() != ts.size()) {
            std::cout << "Error retrieving custom data, sizes do not match." << std::endl;
            exit(2);
        }

        if (naive_recovered_values.size() != ts.size()) {
            std::cout << "Error retrieving naive data, sizes do not match." << std::endl;
            exit(2);
        }

        // replace placeholder value
        if (had_zeroes) {
            // std::ranges::replace(recovered_values.begin(), recovered_values.end(), replacement_value, 0);
            // std::ranges::replace(ts.begin(), ts.end(), replacement_value, 0);
            std::ranges::transform(recovered_values.begin(), recovered_values.end(), recovered_values.begin(),
                                   [](const int64_t value) { return value - 1; });
            std::ranges::transform(naive_recovered_values.begin(), naive_recovered_values.end(),
                                   naive_recovered_values.begin(),
                                   [](const int64_t value) { return value - 1; }
            );
            std::ranges::transform(ts.begin(), ts.end(), ts.begin(),
                                   [](const int64_t value) { return value - 1; });
        }
        std::cout << "Checking for data mismatches" << std::endl;
        // convert back to int, then double
        // auto recovered_ts = std::vector<int64_t>(ts_size);
        // for (size_t ts_counter = 0; ts_counter < ts_size; ts_counter++) {
        //     const auto ref_bytes = reinterpret_cast<const int64_t *>(ts.data());
        //     recovered_ts[ts_counter] = ref_bytes[ts_counter] / static_cast<int64_t>(std::pow(10.0, value_multiplier));
        // }
        auto recovered_ts = recovered_values;

        // define tolerance as 5e-6 because of the precision limitations of
        // floating point representation.
        // for example, if there are 6 decimals, there might be small mismatches at the 5th and 6th decimals.
        // if the mismatch is bigger than 5e-6, it means that when rounded, it could actually cause a mismatch at the
        // previous decimal; else, it could just be an issue with floating point representation.
        for (size_t ts_counter = 0; ts_counter < ts_size; ts_counter++) {
            if (constexpr double tolerance = 5e-6;
                std::abs(static_cast<double>(recovered_ts[ts_counter] - ts[ts_counter])) > tolerance
            ) {
                std::cout << "Values do not match between original and recovered time series at index "
                        << ts_counter << std::endl;
                std::cout << "Original value: " << ts[ts_counter] << std::endl;
                std::cout << "Recovered value: " << recovered_ts[ts_counter] << std::endl;
                exit(3);
            }
        }

        std::uintmax_t original_file_size, compressed_naive_file_size, csv_file_size;

        try {
            original_file_size = std::filesystem::file_size(binary_filename);
            compressed_naive_file_size = std::filesystem::file_size(std::to_string(getpid()) + ".rev_ref");
            csv_file_size = std::filesystem::file_size(
                "/home/lucasala/LBD/BuildTSReference/data/1week-10190.csv"
            );
        } catch (const std::filesystem::filesystem_error &e) {
            std::cerr << "Error retrieving the file system size of the original and compressed files: " << e.what() <<
                    std::endl;

            return 5;
        }

        // Timestamps compression
        sdsl::int_vector<sizeof(int64_t)> timestamps_vector(ts_size);
        for (size_t timestamp_counter = 0; timestamp_counter < ts_size; timestamp_counter++) {
            timestamps_vector[timestamp_counter] = timestamps[timestamp_counter];
        }
        sdsl::bit_vector timestamps_bit_vector(timestamps_vector.size() * sizeof(int64_t));

        for (size_t timestamp_counter = 0; timestamp_counter < ts_size; timestamp_counter++) {
            const auto ref_bytes = reinterpret_cast<const uint8_t *>(&timestamps_vector[timestamp_counter]);
            for (size_t byte_idx = 0; byte_idx < sizeof(int64_t); byte_idx++) {
                timestamps_bit_vector[timestamp_counter * sizeof(int64_t) + byte_idx] = ref_bytes[byte_idx];
            }
        }

        sdsl::sd_vector<> compressed_timestamps(timestamps_bit_vector);


        const std::uintmax_t compressed_timestamps_size = sdsl::size_in_bytes(compressed_timestamps);
        // compressed_timestamps.high.size() / 8 // high bitvector size in bytes
        // + compressed_timestamps.low.size() * sizeof(compressed_timestamps.low[0]) // low array size
        // + sizeof(compressed_timestamps); // object overhead

        // sizeof(RLZ::factor_type)
        // sizeof(uint8_t)
        // const std::uintmax_t total_compressed_size = compressed_file_size + factors.size() * 8;

        // const std::uintmax_t total_compressed_size = compressed_file_size
        //                                              + factors.size() * sizeof(RLZ::factor_type)
        //                                              + compressed_timestamps_size;
        //
        // const std::uintmax_t total_uncompressed_size = ts_memory_size
        //                                                + reference_size * sizeof(int64_t)
        //                                                + sizeof(int64_t) * ts_size;

        const uint64_t uncompressed_reference_size = reference_size * sizeof(int64_t);
        const uint64_t effective_custom_reference_size =
                reference.size()
                * static_cast<uint64_t>(ceil(log2(pow(10, value_multiplier)))) / 8; // convert bits to bytes
        const uint64_t effective_naive_reference_size =
                rlz_naive.reference.size()
                * static_cast<uint64_t>(ceil(log2(pow(10, value_multiplier)))) / 8; // convert bits to bytes
        const uint64_t custom_compressed_ts_size_bits = (
                                                            static_cast<uint64_t>(ceil(log2(max_factors_length)))
                                                            + static_cast<uint64_t>(ceil(log2(reference.size())))
                                                        ) * custom_factors.size();
        const uint64_t naive_compressed_ts_size_bits = (
                                                           static_cast<uint64_t>(ceil(log2(max_naive_factors_length)))
                                                           + static_cast<uint64_t>(ceil(
                                                               log2(rlz_naive.reference.size())))
                                                       ) * naive_factors.size();
        const uint64_t compressed_custom_ts_size = custom_compressed_ts_size_bits / 8; // convert bits to bytes
        const uint64_t compressed_naive_ts_size = naive_compressed_ts_size_bits / 8; // convert bits to bytes
        const uint64_t uncompressed_ts_size = ts_size * sizeof(int64_t);
        // uint64_t effective_ts_size = integer_ts.size() * static_cast<uint64_t>(round(log2(10^value_multiplier)));

        const uint64_t total_custom_compressed_size = compressed_custom_ts_size + effective_custom_reference_size;
        const uint64_t total_naive_compressed_size = compressed_naive_ts_size + effective_naive_reference_size;
        const uint64_t total_uncompressed_size = uncompressed_ts_size + uncompressed_reference_size;

        const double compression_ratio_with_csv =
                static_cast<double>(total_custom_compressed_size) / static_cast<double>(csv_file_size) * 100.0;
        const double compression_ratio_with_python_compressed =
                static_cast<double>(total_custom_compressed_size) / static_cast<double>(original_file_size) * 100.0;
        const double compression_ratio_custom_reference_vs_uncompressed =
                static_cast<double>(total_custom_compressed_size) / static_cast<double>(total_uncompressed_size) *
                100.0;
        const double compression_ratio_naive_reference_vs_uncompressed =
                static_cast<double>(total_naive_compressed_size) / static_cast<double>(total_uncompressed_size) * 100.0;

        std::cout << "Data matches once rebuilt." << std::endl;
        std::cout << "CSV file size: " << csv_file_size << std::endl;
        std::cout << "Time Series memory size: " << uncompressed_ts_size << std::endl;
        std::cout << "Custom reference memory size: " << uncompressed_reference_size << std::endl;
        std::cout << "Naive reference memory size: " << rlz_naive.reference.size() * sizeof(int64_t) <<
                std::endl;
        std::cout << "Python compressed file size: " << original_file_size << std::endl;
        std::cout << "Compressed TS file size with Custom RLZ: " << compressed_custom_file_size << std::endl;
        std::cout << "Compressed TS file size with Naive RLZ: " << compressed_naive_file_size << std::endl;
        std::cout << "Effective (bit-wise) custom compressed TS size: " << compressed_custom_ts_size <<
                std::endl;
        std::cout << "Effective (bit-wise) naive compressed TS size: " << compressed_naive_ts_size <<
                std::endl;
        std::cout << "Effective (bit-wise) custom reference memory size: " << effective_custom_reference_size <<
                std::endl;
        std::cout << "Effective (bit-wise) naive reference memory size: " << effective_naive_reference_size <<
                std::endl;
        std::cout << "Largest custom factor: " << max_factors_length << std::endl;
        std::cout << "Largest naive factor: " << max_naive_factors_length << std::endl;
        std::cout << "Custom factor count: " << custom_factors.size() << std::endl;
        std::cout << "Naive factor count: " << naive_factors.size() << std::endl;
        std::cout << "Compressed timestamps size: " << compressed_timestamps_size << std::endl;
        std::cout << "Total compressed size (TS + Reference): " <<
                total_custom_compressed_size << std::endl;
        std::cout << "Compression ratio with CSV: " <<
                compression_ratio_with_csv <<
                "% (-" << 100.0 - compression_ratio_with_csv << "%)" << std::endl << std::flush;
        std::cout << "Compression ratio with Python compressed: " <<
                compression_ratio_with_python_compressed <<
                "% (-" << 100.0 - compression_ratio_with_python_compressed <<
                "%)" << std::endl << std::flush;
        std::cout << "Compression ratio of custom reference vs uncompressed memory usage: " <<
                compression_ratio_custom_reference_vs_uncompressed
                <<
                "% (-" << 100.0 - compression_ratio_custom_reference_vs_uncompressed << "%)" << std::endl << std::flush;
        std::cout << "Compression ratio with naive reference vs uncompressed memory usage: " <<
                compression_ratio_naive_reference_vs_uncompressed <<
                "% (-" << 100.0 - compression_ratio_naive_reference_vs_uncompressed << "%)" << std::endl << std::flush;

        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
