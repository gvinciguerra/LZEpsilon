// This file is part of LZEpsilon <https://github.com/gvinciguerra/LZEpsilon>.
// Copyright (c) 2022 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include "piecewise_linear_model.hpp"
#include <sdsl/int_vector.hpp>
#include <cstdint>
#include <type_traits>

/** Computes (bits_per_correction > 0 ? 2^(bits_per_correction-1) - 1 : 0) without the conditional operator. */
#define BPC_TO_EPSILON(bits_per_correction) (((1ul << (bits_per_correction)) + 1) / 2 - 1)

/** Computes the number of bits needed to store x, that is, 0 if x is 0, 1 + floor(log2(x)) otherwise. */
#define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))

/** Computes the smallest integral value not less than x / y, where x and y must be positive integers. */
#define CEIL_UINT_DIV(x, y) ((x) / (y) + ((x) % (y) != 0))

#pragma pack(push, 1)

template<typename K>
struct TailSegment {
    K slope_numerator;
    K slope_denominator;
    std::make_signed_t<K> intercept;
    using larger_key_type = typename std::conditional_t<sizeof(K) <= 4, int64_t, __int128>;
    using size_type = size_t;

    TailSegment() = default;

    template<typename RandomIt>
    TailSegment(RandomIt begin, RandomIt end, uint8_t bpc, size_t &length) {
        const auto n = (size_t) std::distance(begin, end);
        const auto epsilon = BPC_TO_EPSILON(bpc);

        OptimalPiecewiseLinearModel<K, K> opt(epsilon);
        opt.add_point(0, begin[0]);

        for (length = 1; length < n; ++length)
            if (!opt.add_point(length, begin[length]))
                break;

        auto cs = opt.get_segment();
        auto max_slope = cs.rectangle[3] - cs.rectangle[1];
        auto intercept_numerator = cs.rectangle[3].x * cs.rectangle[1].y - cs.rectangle[1].x * cs.rectangle[3].y;
        slope_numerator = max_slope.dy;
        slope_denominator = std::max<K>(1, max_slope.dx);
        intercept = max_slope.dx == 0 ? begin[0] : intercept_numerator / max_slope.dx;
    }

    larger_key_type approximate(size_t i) const {
        return larger_key_type(slope_numerator * i) / slope_denominator + intercept;
    }

    std::pair<size_t, size_t> approximate_position(const K &value, uint64_t epsilon, size_t n) const {
        auto numerator = std::max<larger_key_type>(1, slope_numerator);
        auto position = ((larger_key_type(value) - intercept) * slope_denominator) / numerator;
        auto bound = 1 + (epsilon * slope_denominator) / numerator;
        return {std::clamp<larger_key_type>(position, 0, n), bound};
    }

    uint64_t get_int(size_t i, const sdsl::int_vector<> &corrections, size_t offset, uint64_t epsilon) const {
        if (epsilon == 0)
            return approximate(i);
        return approximate(i) + corrections[offset + i] - epsilon;
    }

    size_t rank(const K &value, const sdsl::int_vector<> &corrections, size_t offset, uint64_t epsilon, size_t n) const {
        auto[pos, bound] = approximate_position(value, epsilon, n);
        auto lo = pos <= bound ? 0 : pos - bound;
        auto hi = std::min<size_t>(pos + bound + 1, n);

        if (epsilon == 0) {
            while (lo < hi) {
                auto mid = lo + (hi - lo) / 2;
                if (approximate(mid) < value)
                    lo = mid + 1;
                else
                    hi = mid;
            }
        } else {
            while (lo < hi) {
                auto mid = lo + (hi - lo) / 2;
                if (approximate(mid) + corrections[offset + mid] - epsilon < value)
                    lo = mid + 1;
                else
                    hi = mid;
            }
        }

        return lo;
    }
};

#pragma pack(pop)