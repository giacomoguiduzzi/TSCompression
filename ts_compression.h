//
// Created by giacomo on 23/04/25.
//

#ifndef TS_COMPRESS_H
#define TS_COMPRESS_H
#include <rlz_naive.hpp>
typedef rct::rlz_naive<
    uint8_t,
    rct::reference_uniform_sample<uint8_t>,
    sdsl::csa_bitcompressed<sdsl::int_alphabet<> >
>
BinaryRLZ;
typedef rct::rlz_naive<
    int64_t,
    rct::reference_uniform_sample<int64_t>,
    sdsl::csa_bitcompressed<sdsl::int_alphabet<> >
>
NaiveRLZ;

typedef rct::rlz_naive<
    int64_t,
    std::vector<int64_t>,
    sdsl::csa_bitcompressed<sdsl::int_alphabet<> >
>
CustomRLZ;

#endif //TS_COMPRESS_H
