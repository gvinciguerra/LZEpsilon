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

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp_bitcompressed.hpp>
#include <sdsl/rmq_support.hpp>
#include "TailSegment.hpp"

#include <string_view>
#include <map>

using TailType = TailSegment<uint32_t>;

class LZEpsilon {
    /** Elias-Fano representation without efficient rank support. */
    using elias_fano = sdsl::sd_vector<sdsl::bit_vector, sdsl::bit_vector::select_1_type, sdsl::select_support_scan<0>>;

    size_t text_length;
    uint64_t epsilon;
    sdsl::sd_vector<> ending_positions; ///< Stores the left and right boundaries of the segments.
    sdsl::sd_vector<>::select_1_type ending_positions_select;
    sdsl::sd_vector<>::rank_1_type ending_positions_rank;
    sdsl::int_vector<> head_sources;
    std::vector<TailType> tails;
    sdsl::int_vector<> corrections;
    elias_fano offsets; ///< Marks the beginning of each segment's corrections.
    elias_fano::select_1_type offsets_select;
    std::map<std::string, std::string> metadata;

    struct LZEpsilonPhrase {
        uint32_t source;      ///< The phrase identifier where the source of this phrase ends.
        uint32_t head_length; ///< The length of the head of the phrase.
        uint32_t tail_length; ///< The length of the tail of the phrase.
        TailType tail;        ///< The tail segment, which covers the last head position if the head is nonempty.
        LZEpsilonPhrase(uint32_t source, uint32_t head_length, uint32_t tail_length, TailType tail) :
            source(source), head_length(head_length), tail_length(tail_length), tail(tail) {}
    };

public:

    using size_type = size_t;

    LZEpsilon() = default;

    template<class IntVector>
    explicit LZEpsilon(const IntVector &v, uint8_t bpc) : text_length(v.size()), epsilon(BPC_TO_EPSILON(bpc)) {
        auto t0 = std::chrono::high_resolution_clock::now();
        std::vector<LZEpsilonPhrase> phrases;
        phrases.reserve(text_length / (1 + sdsl::bits::hi(text_length)));

        using rk_t = sdsl::rank_support_v<1>;
        using sl1_t = sdsl::select_support_scan<1>;
        using sl0_t = sdsl::select_support_scan<0>;
        sdsl::csa_wt_int<sdsl::wt_huff_int<sdsl::bit_vector, rk_t, sl1_t, sl0_t>, 1, 1> sa;
        sdsl::rmq_support_sparse_table<decltype(sa), false> rmq;
        construct_im_reversed(v, sa, rmq);

        std::map<size_t, size_t> position_to_phrase; // maps ending positions in the text to their phrase id
        position_to_phrase.emplace(text_length + 1, 0);
        size_t phrase_start = 0;
        size_t total_corrections = 0;
        size_t prev_segment_start = 0;
        size_t short_tails_count = 0;
        size_t total_head_length = 0;
        size_t total_tail_length = 0;

        while (phrase_start < text_length) {
            // Compute the phrase head
            typename decltype(sa)::size_type sp = 0, ep = text_length;
            size_t ip = phrase_start;
            size_t head_past_end = phrase_start;
            size_t source_end_pos;
            size_t q = 0;

            while (ip < text_length) {
                auto matches = sdsl::backward_search(sa, sp, ep, v[ip] - (ip > 0 ? v[ip - 1] : 0), sp, ep);
                if (!matches || sa[rmq(sp, ep)] <= text_length - phrase_start - 1)
                    break;
                ++ip;
                auto [fpos, qp] = *position_to_phrase.lower_bound(sp);
                if (fpos <= ep) {
                    head_past_end = ip;
                    q = qp;
                    source_end_pos = fpos;
                }
            }

            // Possibly shorten the previous tail by extending this phrase head leftwards
            auto head_length = head_past_end - phrase_start;
            auto is_head_nonempty = head_length != 0;

            if (bpc > 0 && is_head_nonempty && phrases.back().tail_length > 2
                && (source_end_pos = text_length - sa[source_end_pos] - 1) >= head_length) {
                auto a = phrase_start - 1;             // index before the start of the phrase
                auto b = source_end_pos - head_length; // index before the start of the source phrase
                //for (auto k = 1; k <= head_length; ++k) {
                //    auto left = v[a + k] - v[a + k - 1];
                //    auto right = v[b + k] - v[b + k - 1];
                //    assert(left == right);
                //}
                size_t shortened_length = 0;
                for (auto k = 0; k + 2 < phrases.back().tail_length && b >= k && a - k > source_end_pos; ++k) {
                    auto left = v[a - k] - v[a - k - 1];
                    auto right = b - k == 0 ? v[0] : v[b - k] - v[b - k - 1];
                    if (left != right)
                        break;
                    else
                        ++shortened_length;
                }

                if (shortened_length > 0) {
                    auto f = position_to_phrase.find(sa.isa[text_length - phrase_start]);
                    position_to_phrase.erase(sa.isa[text_length - phrase_start]);
                    head_length += shortened_length;
                    phrase_start -= shortened_length;
                    total_corrections -= shortened_length;
                    phrases.back().tail_length -= shortened_length;
                    if (phrases.back().tail_length + (phrases.back().head_length > 0) == 2) {
                        size_t tmp;
                        phrases.back().tail = TailType(v.begin() + prev_segment_start, v.end(), 0, tmp);
                        total_corrections -= 2;
                    }
                    position_to_phrase.emplace(sa.isa[text_length - phrase_start], phrases.size() - 1);
                }
            }

            // Compute the phrase tail
            prev_segment_start = head_past_end - is_head_nonempty;
            size_t segment_length;
            TailType tail(v.begin() + prev_segment_start, v.end(), bpc, segment_length);
            if (segment_length <= 2) {
                tail = TailType(v.begin() + prev_segment_start, v.end(), 0, segment_length);
                ++short_tails_count;
            }

            auto tail_length = segment_length - is_head_nonempty;
            position_to_phrase.emplace(sa.isa[text_length - head_past_end - tail_length], phrases.size());
            phrases.emplace_back(q, head_length, tail_length, tail);
            total_corrections += segment_length <= 2 ? 0 : segment_length;
            phrase_start = head_past_end + tail_length;

            total_head_length += head_length;
            total_tail_length += tail_length;
        }

        // Store the phrases and the corrections
        phrase_start = 0;
        corrections.width(bpc);
        corrections.reserve(total_corrections);
        sdsl::int_vector<> tmp_offsets;
        tmp_offsets.reserve(phrases.size());
        sdsl::int_vector<> tmp_ending_positions;
        tmp_ending_positions.reserve(phrases.size() * 2);
        for (auto &p: phrases) {
            tmp_offsets.push_back(corrections.size());
            head_sources.push_back(p.source);
            tails.push_back(p.tail);

            auto is_head_nonempty = p.head_length > 0;
            auto segment_start = phrase_start + p.head_length - is_head_nonempty;
            auto segment_length = p.tail_length + is_head_nonempty;

            if (bpc > 0 && segment_length > 2) {
                for (size_t i = 0; i < segment_length; ++i) {
                    auto error = v[segment_start + i] - p.tail.approximate(i);
                    auto correction = uint64_t(error + epsilon);
                    if (BIT_WIDTH(correction) > bpc)
                        throw std::overflow_error("Segment correction too large");
                    corrections.push_back(correction);
                }
            }

            tmp_ending_positions.push_back(phrase_start + p.head_length - is_head_nonempty);
            if (&p != &phrases.back())
                tmp_ending_positions.push_back(phrase_start + p.head_length + p.tail_length - 1);

            phrase_start += p.head_length + p.tail_length;
        }

        tmp_ending_positions.push_back(text_length - 1);

        sdsl::util::bit_compress(head_sources);

        if (bpc > 0)
            offsets = decltype(offsets)(tmp_offsets.begin(), tmp_offsets.end());
        sdsl::util::init_support(offsets_select, &offsets);
        ending_positions = decltype(ending_positions)(tmp_ending_positions.begin(), tmp_ending_positions.end());
        sdsl::util::init_support(ending_positions_select, &ending_positions);
        sdsl::util::init_support(ending_positions_rank, &ending_positions);

        // Store metadata
        auto t1 = std::chrono::high_resolution_clock::now();
        auto construction_s = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() / 1000.;
        auto to_bpi = [&](auto bytes) { return std::to_string(bytes * 8. / text_length); };
        metadata["tails_bpi"] = to_bpi(tails.size() * sizeof(TailType));
        metadata["ending_positions_bpi"] = to_bpi(sdsl::size_in_bytes(ending_positions));
        metadata["head_sources_bpi"] = to_bpi(sdsl::size_in_bytes(head_sources));
        metadata["corrections_bpi"] = to_bpi(sdsl::size_in_bytes(corrections));
        metadata["corrections_offsets_bpi"] = to_bpi(sdsl::size_in_bytes(offsets));
        metadata["phrases_count"] = std::to_string(tails.size());
        metadata["avg_head_length"] = std::to_string(total_head_length / double(tails.size()));
        metadata["avg_tail_length"] = std::to_string(total_tail_length / double(tails.size()));
        metadata["construction_s"] = std::to_string(construction_s);
        metadata["short_tails_count"] = std::to_string(short_tails_count);
    }

    std::map<std::string, std::string> get_metadata() { return metadata; }

    uint64_t select(size_t i) const {
        assert(i > 0);
        --i;
        return select_aux(i, ending_positions_rank(i) / 2);
    }

    size_t rank(uint64_t x) const {
        size_t phrase_id = 0;
        auto count = phrases_count();
        while (count > 0) {
            auto step = count / 2;
            auto value = get_end_value(phrase_id + step);
            if (value < x) {
                phrase_id = phrase_id + step + 1;
                count -= step + 1;
            } else
                count = step;
        }
        return rank_aux(x, phrase_id);
    }

    size_t size() const { return text_length; }
    size_t phrases_count() const { return tails.size(); }
    size_t size_in_bytes() const {
        return sdsl::size_in_bytes(ending_positions) + sdsl::size_in_bytes(ending_positions_select)
            + sdsl::size_in_bytes(ending_positions_rank) + sdsl::size_in_bytes(head_sources)
            + sdsl::size_in_bytes(corrections) + sdsl::size_in_bytes(offsets)
            + sizeof(tails[0]) * tails.size();
    }

    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, const std::string &name = "") const {
        auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        written_bytes += sdsl::write_member(text_length, out, child, "text_length");
        written_bytes += sdsl::write_member(epsilon, out, child, "epsilon");
        written_bytes += sdsl::serialize(ending_positions, out, child, "ending_positions");
        written_bytes += sdsl::serialize(ending_positions_select, out, child, "ending_positions_select");
        written_bytes += sdsl::serialize(ending_positions_rank, out, child, "ending_positions_rank");
        written_bytes += sdsl::serialize(head_sources, out, child, "head_sources");
        written_bytes += sdsl::serialize(corrections, out, child, "corrections");
        written_bytes += sdsl::serialize(offsets, out, child, "offsets");
        written_bytes += sdsl::serialize(offsets_select, out, child, "offsets_select");
        written_bytes += sdsl::serialize(tails, out, child, "tails");

        sdsl::write_member(metadata.at("tails_bpi"), out, child, "tails_bpi");
        sdsl::write_member(metadata.at("ending_positions_bpi"), out, child, "ending_positions_bpi");
        sdsl::write_member(metadata.at("head_sources_bpi"), out, child, "head_sources_bpi");
        sdsl::write_member(metadata.at("corrections_bpi"), out, child, "corrections_bpi");
        sdsl::write_member(metadata.at("corrections_offsets_bpi"), out, child, "corrections_offsets_bpi");
        sdsl::write_member(metadata.at("phrases_count"), out, child, "phrases_count");
        sdsl::write_member(metadata.at("avg_head_length"), out, child, "avg_head_length");
        sdsl::write_member(metadata.at("avg_tail_length"), out, child, "avg_tail_length");
        sdsl::write_member(metadata.at("construction_s"), out, child, "construction_s");
        sdsl::write_member(metadata.at("short_tails_count"), out, child, "short_tails_count");

        return written_bytes;
    }

    void load(std::istream &in) {
        sdsl::read_member(text_length, in);
        sdsl::read_member(epsilon, in);
        sdsl::load(ending_positions, in);
        ending_positions_select.load(in, &ending_positions);
        ending_positions_rank.load(in, &ending_positions);
        sdsl::load(head_sources, in);
        sdsl::load(corrections, in);
        sdsl::load(offsets, in);
        offsets_select.load(in, &offsets);
        sdsl::load(tails, in);

        sdsl::read_member(metadata["tails_bpi"], in);
        sdsl::read_member(metadata["ending_positions_bpi"], in);
        sdsl::read_member(metadata["head_sources_bpi"], in);
        sdsl::read_member(metadata["corrections_bpi"], in);
        sdsl::read_member(metadata["corrections_offsets_bpi"], in);
        sdsl::read_member(metadata["phrases_count"], in);
        sdsl::read_member(metadata["avg_head_length"], in);
        sdsl::read_member(metadata["avg_tail_length"], in);
        sdsl::read_member(metadata["construction_s"], in);
        sdsl::read_member(metadata["short_tails_count"], in);
    }

private:

    size_t get_offset(size_t phrase_id, uint64_t eps) const { return eps ? offsets_select(phrase_id + 1) : 0; }

    uint64_t select_aux(size_t i, size_t q) const {
        assert(q == 0 || i > get_phrase_end(q - 1));
        assert(i <= get_phrase_end(q));

        auto get_tail_value = [&](size_t phrase, size_t pos, size_t eps) {
            return tails[phrase].get_int(pos, corrections, get_offset(phrase, eps), eps);
        };

        auto [phrase_start, segment_start, phrase_end, q_epsilon] = get_phrase_endpoints(q);
        if (i >= segment_start)
            return get_tail_value(q, i - segment_start, q_epsilon);

        auto head_end = segment_start;
        auto value_at_head_end = get_tail_value(q, 0, q_epsilon);

        auto r = head_sources[q];
        auto [source_start, source_segment_start, source_end, r_epsilon] = get_phrase_endpoints(r);
        auto value_at_source_end = get_tail_value(r, source_end - source_segment_start, r_epsilon);
        auto delta = value_at_head_end - value_at_source_end;

        auto ip = i - (head_end - source_end);
        auto qp = r;
        while (qp > 0 && ip <= get_phrase_end(qp - 1))
            --qp;
        return delta + select_aux(ip, qp);
    }

    size_t rank_aux(uint64_t x, size_t q) const {
        assert(q == 0 || get_end_value(q - 1) < x);
        assert(x <= get_end_value(q));

        auto get_tail_value = [&](size_t phrase, size_t pos, size_t eps) {
            return tails[phrase].get_int(pos, corrections, get_offset(phrase, eps), eps);
        };

        auto [phrase_start, segment_start, phrase_end, q_epsilon] = get_phrase_endpoints(q);
        auto head_end = segment_start;
        auto value_at_head_end = get_tail_value(q, 0, q_epsilon);

        if (x >= value_at_head_end || phrase_start == segment_start)
            return segment_start
                + tails[q].rank(x, corrections, get_offset(q, q_epsilon), q_epsilon, phrase_end - segment_start + 1);

        auto r = head_sources[q];
        auto [source_start, source_segment_start, source_end, r_epsilon] = get_phrase_endpoints(r);
        auto value_at_source_end = get_tail_value(r, source_end - source_segment_start, r_epsilon);
        auto delta = value_at_head_end - value_at_source_end;

        auto xp = x - delta;
        auto qp = r;
        while (qp > 0 && xp <= get_end_value(qp - 1))
            --qp;
        return head_end - source_end + rank_aux(xp, qp);
    }

    /** Returns the value at the last position of a phrase. */
    uint64_t get_end_value(size_t phrase_id) const {
        auto [phrase_start, segment_start, phrase_end, eps] = get_phrase_endpoints(phrase_id);
        return tails[phrase_id].get_int(phrase_end - segment_start, corrections, get_offset(phrase_id, eps), eps);
    }

    /** Returns a tuple with the first position of the phrase, the first position of the segment, the last position
     * of the phrase, and the epsilon value of the segment. */
    std::tuple<size_t, size_t, size_t, uint64_t> get_phrase_endpoints(size_t phrase_id) const {
//        auto phrase_start = phrase_id > 0 ? ending_positions_select(phrase_id * 2) + 1 : 0;
//        auto segment_start = ending_positions_select(phrase_id * 2 + 1);
//        auto phrase_end = ending_positions_select(phrase_id * 2 + 2);

        auto lo_ptr = phrase_id * 2 + 2;
        auto hi_ptr = ending_positions.high_1_select(phrase_id * 2 + 2);
        auto phrase_end = ending_positions.low[lo_ptr - 1] + ((hi_ptr + 1 - lo_ptr) << ending_positions.wl);
        assert(phrase_end == ending_positions_select(phrase_id * 2 + 2));

        --lo_ptr;
        hi_ptr = sdsl::util::prev_bit(ending_positions.high, hi_ptr - 1);
        auto segment_start = ending_positions.low[lo_ptr - 1] + ((hi_ptr + 1 - lo_ptr) << ending_positions.wl);
        assert(segment_start == ending_positions_select(phrase_id * 2 + 1));

        size_t phrase_start = 0;
        if (phrase_id > 0) {
            --lo_ptr;
            hi_ptr = sdsl::util::prev_bit(ending_positions.high, hi_ptr - 1);
            phrase_start = 1 + ending_positions.low[lo_ptr - 1] + ((hi_ptr + 1 - lo_ptr) << ending_positions.wl);
            assert(phrase_start == 1 + ending_positions_select(phrase_id * 2));
        }

        auto eps = phrase_end - segment_start <= 1 ? 0 : epsilon;
        return {phrase_start, segment_start, phrase_end, eps};
    }

    /** Returns the last position of the given phrase, i.e. the index where its tail ends. */
    size_t get_phrase_end(size_t phrase_id) const { return ending_positions_select(phrase_id * 2 + 2); }

    template<class IntVector, class SA, class RMQ>
    static void construct_im_reversed(IntVector &s, SA &sa, RMQ &rmq) {
        using namespace sdsl;
        auto tmp_file = ram_file_name(util::to_string(util::pid()) + "_" + util::to_string(util::id()));
        store_reversed_gap_string_to_file(s, tmp_file);
        cache_config config(false, "@");
        construct(sa, tmp_file, config, 0);
        config.delete_files = true;
        config.delete_data = true;
        rmq = {&sa};
        ram_fs::remove(tmp_file);
    }

    template<uint8_t t_width>
    static bool store_reversed_gap_string_to_file(const sdsl::int_vector<t_width> &v, const std::string &file) {
        using namespace sdsl;
        osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
        if (!out) {
            std::cerr << "ERROR: store_to_file_reversed: Could not open file `" << file << "`" << std::endl;
            return false;
        }

        sdsl::int_vector<t_width> v_reversed;
        v_reversed.reserve(v.size());
        for (auto i = v.size() - 1; i >= 1; --i)
            v_reversed.push_back(uint64_t(v[i] - v[i - 1]));
        v_reversed.push_back(v[0] > 0 ? v[0] : 1);
        v_reversed.serialize(out, nullptr, "");
        out.close();
        return true;
    }
};
