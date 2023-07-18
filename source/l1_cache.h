#ifndef L1_CACHE_H
#define L1_CACHE_H

#include "types.h"
#include "address.h"
#include "utils.h"
#include <ap_int.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wpedantic"
#pragma GCC diagnostic error "-Wall"
#pragma GCC diagnostic error "-Wextra"
#pragma GCC diagnostic ignored "-Wunused-label"

using namespace types;

template <typename WORD_TYPE, size_t MAIN_SIZE, size_t N_SETS, size_t N_WAYS,
	 size_t N_WORDS_PER_LINE, bool SWAP_TAG_SET, storage_impl_type STORAGE_IMPL>
class l1_cache {
	private:
		static const size_t ADDR_SIZE = utils::log2_ceil(MAIN_SIZE);
		static const size_t SET_SIZE = utils::log2_ceil(N_SETS);
		static const size_t OFF_SIZE = utils::log2_ceil(N_WORDS_PER_LINE);
		static const size_t TAG_SIZE = (ADDR_SIZE - (SET_SIZE + OFF_SIZE));
		static const size_t WAY_SIZE = utils::log2_ceil(N_WAYS);
		static const size_t N_LINES = (((N_SETS * N_WAYS) > 0) ?
				(N_SETS * N_WAYS) : 1);

		static_assert(((MAIN_SIZE > 0) && ((1 << ADDR_SIZE) == MAIN_SIZE)),
				"MAIN_SIZE must be a power of 2 greater than 0");
		static_assert(((N_SETS == 0) || (1 << SET_SIZE) == N_SETS),
				"N_SETS must be a power of 2");
		static_assert(((N_WAYS == 0) || ((1 << WAY_SIZE) == N_WAYS)),
				"N_WAYS must be a power of 2");
		static_assert(((N_WORDS_PER_LINE > 0) &&
					((1 << OFF_SIZE) == N_WORDS_PER_LINE)),
				"N_WORDS_PER_LINE must be a power of 2 greater than 0");
		static_assert((MAIN_SIZE >= (N_SETS * N_WAYS * N_WORDS_PER_LINE)),
				"N_SETS and/or N_WAYS and/or N_WORDS_PER_LINE are too big for the specified MAIN_SIZE");

		typedef WORD_TYPE line_type[N_WORDS_PER_LINE];
		typedef address<ADDR_SIZE, TAG_SIZE, SET_SIZE, WAY_SIZE, SWAP_TAG_SET>
			addr_type;
		typedef replacer<false, addr_type, ((N_SETS > 0) ? N_SETS : 1),
			((N_WAYS > 0) ? N_WAYS : 1), N_WORDS_PER_LINE> replacer_type;

		ap_uint<(TAG_SIZE > 0) ? TAG_SIZE + 1 : 1> m_tag[N_LINES] = {0};	// 1
		// ap_uint<N_LINES> m_valid;				// 2
		WORD_TYPE m_cache_mem[N_LINES][N_WORDS_PER_LINE];	// 3
		replacer_type m_replacer;				// 4

	public:
		l1_cache() {
#pragma HLS array_partition variable=m_cache_mem type=complete dim=2
#pragma HLS array_partition variable=m_tag type=complete dim=0

			switch (STORAGE_IMPL) {
				case URAM:
#pragma HLS bind_storage variable=m_cache_mem type=RAM_2P impl=URAM
					break;
				case BRAM:
#pragma HLS bind_storage variable=m_cache_mem type=RAM_2P impl=BRAM latency=3
					break;
				case LUTRAM:
#pragma HLS bind_storage variable=m_cache_mem type=RAM_2P impl=LUTRAM
					break;
				default:
					break;
			}
		}

		void init() {
#pragma HLS inline
			for(int g = 0; g < (int)N_LINES; g++){
				#pragma HLS unroll
				m_tag[g] = 0;
			}
			m_replacer.init();
		}

		void set_tag(const ap_uint<ADDR_SIZE> addr_main){
			#pragma HLS inline
			addr_type addr(addr_main);
			// addr.set_way(0);
			// m_tag[addr.m_set] = ((ap_uint<1>)1, addr.m_tag);
			m_tag[addr.m_set] = addr.m_tag;
			m_tag[addr.m_set][TAG_SIZE] = 1;
		}

		bool get_line(const ap_uint<ADDR_SIZE> addr_main, line_type line) const {
#pragma HLS inline
			addr_type addr(addr_main);
			const auto way = hit(addr);

			if (way == -1)
				return false;

			addr.set_way(way);
			for (size_t off = 0; off < N_WORDS_PER_LINE; off++) {
#pragma HLS unroll
				line[off] = m_cache_mem[addr.m_addr_line][off];
			}

			return true;
		}

		void set_line(const ap_uint<ADDR_SIZE> addr_main,
				const line_type line) {
#pragma HLS inline
			addr_type addr(addr_main);

			addr.set_way(m_replacer.get_way(addr));

			for (size_t off = 0; off < N_WORDS_PER_LINE; off++) {
#pragma HLS unroll
				m_cache_mem[addr.m_addr_line][off] = line[off];
			}
			// m_valid[addr.m_addr_line] = true;
			// m_tag[addr.m_addr_line] = addr.m_tag;

			m_replacer.notify_insertion(addr);
		}

		void notify_write(const ap_uint<ADDR_SIZE> addr_main) {
#pragma HLS inline
			const addr_type addr(addr_main);

			if (hit(addr) != -1)
				m_tag[addr.m_addr_line] = 0;
				// m_valid[addr.m_addr_line] = false;
		}

	private:
		inline int hit(const addr_type &addr) const {
#pragma HLS inline
			addr_type addr_tmp = addr;
			auto hit_way = -1;
			ap_uint<TAG_SIZE + 1> tag_search = addr_tmp.m_tag;
			tag_search[TAG_SIZE] = 1;
			for (size_t way = 0; way < N_WAYS; way++) {
				addr_tmp.set_way(way);
				if (tag_search == m_tag[addr_tmp.m_addr_line]) {
					hit_way = way;
				}
			}
			// std::cout << ((ap_uint<1>)1, addr_tmp.m_tag) << " " << m_tag[addr_tmp.m_addr_line] << std::endl;

			return hit_way;
		}
};

#pragma GCC diagnostic pop

#endif /* L1_CACHE_H */

