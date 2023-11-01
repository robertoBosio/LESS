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
	 size_t N_WORDS_PER_LINE, size_t PORTS, bool SWAP_TAG_SET, storage_impl_type STORAGE_IMPL>
class l1_cache {
	private:
		static const size_t ADDR_SIZE = utils::log2_ceil(MAIN_SIZE);
		static const size_t SET_SIZE = utils::log2_ceil(N_SETS);
		static const size_t OFF_SIZE = utils::log2_ceil(N_WORDS_PER_LINE);
		static const size_t TAG_SIZE = (ADDR_SIZE - (SET_SIZE + OFF_SIZE));
		static const size_t WAY_SIZE = utils::log2_ceil(N_WAYS);
		static const size_t N_LINES = (((N_SETS * N_WAYS) > 0) ?
				(N_SETS * N_WAYS) : 1);
		static const size_t BRAM_LAT = 2;

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

		ap_uint<(TAG_SIZE > 0) ? TAG_SIZE + 1 : 1> m_tag[N_LINES];	// 1
		// ap_uint<N_LINES> m_valid;				// 2
		WORD_TYPE m_cache_mem[N_LINES][N_WORDS_PER_LINE];	// 3
		replacer_type m_replacer;				// 4


		typedef struct{
            ap_uint<((SET_SIZE + WAY_SIZE) > 0) ? (SET_SIZE + WAY_SIZE) : 1> addr_line;
            ap_uint<(TAG_SIZE > 0) ? TAG_SIZE : 1> tag;
		} block_t;

        block_t m_local_raw_cache[BRAM_LAT];
        ap_uint<BRAM_LAT> m_local_valid;
	
	
	public:
		l1_cache() {
			if (PORTS != 1){
#pragma HLS array_partition variable=m_cache_mem type=complete dim=2
			} else {
#pragma HLS array_partition variable=m_cache_mem type=complete dim=3
			}
#pragma HLS array_partition variable=m_local_raw_cache type=complete dim=1
// #pragma HLS array_partition variable=m_tag type=complete dim=0
#pragma HLS bind_storage variable=m_tag type=RAM_2P impl=BRAM latency=1

			switch (STORAGE_IMPL) {
				case URAM:
#pragma HLS bind_storage variable=m_cache_mem type=RAM_2P impl=URAM
					break;
				case BRAM:
#pragma HLS bind_storage variable=m_cache_mem type=RAM_2P impl=BRAM latency=1
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
            // m_valid = 0;
            for (int g = 0; g < (int)N_LINES; g++){
                m_tag[g] = 0;
            }
			m_local_valid = 0;
			m_replacer.init();
		}
        
		void foo(){
#pragma HLS dependence variable=m_tag type=inter direction=RAW false
#pragma HLS dependence variable=m_tag type=inter direction=RAW distance=3 true
        }

		// void set_tag(const ap_uint<ADDR_SIZE> addr_main){
		// 	#pragma HLS inline
		// 	addr_type addr(addr_main);
		// 	addr.set_way(0);
		// 	m_tag[addr.m_set] = ((ap_uint<1>)1, addr.m_tag);
		// 	m_tag[addr.m_set] = addr.m_tag;
		// 	m_tag[addr.m_set][TAG_SIZE] = 1;
		// }

		bool get_line(const ap_uint<ADDR_SIZE> addr_main, line_type line) {
#pragma HLS inline
			addr_type addr(addr_main);
            addr.set_way(0);
            bool hit_last_used = false;
            bool hit_local_storage = false;

            for (size_t s = 0; s < BRAM_LAT; s++){
#pragma HLS unroll
                auto g = BRAM_LAT - s - 1;
                if (m_local_raw_cache[g].addr_line == addr.m_addr_line && m_local_valid[g]){
                    hit_local_storage = true;
                    hit_last_used = (m_local_raw_cache[g].tag == addr.m_tag);
                }
            }
            for (size_t s = 0; s < BRAM_LAT - 1; s++){
#pragma HLS unroll
                auto g = BRAM_LAT - s - 1;
                m_local_raw_cache[g].addr_line = m_local_raw_cache[g - 1].addr_line;
                m_local_raw_cache[g].tag = m_local_raw_cache[g - 1].tag;
				m_local_valid[g] = m_local_valid[g - 1];
            }

            m_local_raw_cache[0].addr_line = addr.m_addr_line;
            m_local_raw_cache[0].tag = addr.m_tag;
            m_local_valid[0] = true;

            int way = -1;
            if (hit_local_storage){
                if (hit_last_used) {
                    way = 0;
                }
            } else {
                way = hit(addr);
            }

			if (way == -1){
				return false;
	        }

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
			m_tag[addr.m_addr_line] = addr.m_tag;
			m_tag[addr.m_addr_line].set(TAG_SIZE);

			m_replacer.notify_insertion(addr);
		}

		void notify_write(const ap_uint<ADDR_SIZE> addr_main) {
#pragma HLS inline
			const addr_type addr(addr_main);

			if (hit(addr) != -1){
				// m_valid[addr.m_addr_line] = false;
                m_tag[addr.m_addr_line].clear(TAG_SIZE);
			}
		}

	private:
		inline int hit(const addr_type &addr) const {
#pragma HLS inline
			addr_type addr_tmp = addr;
			auto hit_way = -1;
			for (size_t way = 0; way < N_WAYS; way++) {
				addr_tmp.set_way(way);
				// if (m_valid[addr_tmp.m_addr_line] && (addr_tmp.m_tag == m_tag[addr_tmp.m_addr_line])) {
				if ((addr_tmp.m_tag == m_tag[addr_tmp.m_addr_line].range(TAG_SIZE - 1, 0)) && m_tag[addr_tmp.m_addr_line].test(TAG_SIZE)) {
					hit_way = way;
				}
			}

			return hit_way;
		}
};

#pragma GCC diagnostic pop

#endif /* L1_CACHE_H */

