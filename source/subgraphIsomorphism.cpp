#include "subgraphIsomorphism.hpp"

#ifndef __SYNTHESIS__
#include <cassert>
#include <fstream>
#include <unordered_map>
#include <thread>
#endif

#include "cache.h"
#include <hls_stream.h>
#include <ap_int.h>
#include <ap_axi_sdata.h>
#include "hls_task.h"
#include "hls_print.h"
#include "Parameters.hpp"
#include "QueryVertex.hpp"
#include "Trie.hpp"
#include "hash_lookup3.hpp"
#include "dynfifo_utils.hpp"
#include "stack_utils.hpp"

#if !SOFTWARE_PREPROC
#include "preprocess.hpp"
#endif /* SOFTWARE_PREPROC */

#if DEBUG_STATS
#include "debug.hpp"
#endif /* DEBUG_STATS */

#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wpedantic"
#pragma GCC diagnostic error "-Wall"
#pragma GCC diagnostic error "-Wextra"
#pragma GCC diagnostic ignored "-Wunused-label"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

#define STOP_S      9    
#define V_ID_W      VERTEX_WIDTH_BIT
#define V_L_W       LABEL_WIDTH
#define MAX_QV      MAX_QUERY_VERTICES
#define MAX_CL      MAX_COLLISIONS
#define MAX_TB      MAX_TABLES
#define C_W         COUNTER_WIDTH
#define E_W         EDGE_WIDTH
#define S_D         DEFAULT_STREAM_DEPTH    
#define DDR_WORDS   RES_WIDTH
#define DDR_W       DDR_WORD

#if CACHE_ENABLE
typedef cache< ap_uint<DDR_W>, true, false, 2,
        HASHTABLES_SPACE, 0, 0, (1UL << CACHE_WORDS_PER_LINE), false, 512, 1,
        false, 1, AUTO, BRAM> htb_cache_t;

// typedef cache< bloom_t, true, false, 1,
//         BLOOM_SPACE, 1, 1, (1UL << K_FUNCTIONS), false, 128, 1,
//         false, 1, AUTO, BRAM> bloom_cache_t;
#endif /* CACHE_ENABLE */

// std::ofstream f("../../../../restot_stack.txt");
/******** Tuple definition ********/
enum edge_flag
{
    MIN_SET = 1,
    NO_EDGE = 2,
    CHECK = 0
};
typedef ap_uint<2> edge_flag_type;

typedef struct
{
    ap_uint<V_ID_W> indexing_v;
    unsigned char tb_index;
    unsigned char iv_pos; // Query indexing vertex position.
    unsigned int address;
    bool reset;
    bool last;
} findmin_tuple_t;

typedef struct
{
    ap_uint<V_ID_W> indexing_v;
    unsigned char tb_index;
    unsigned char iv_pos; // Query indexing vertex position.
} readmin_counter_tuple_t;

typedef struct
{
    ap_uint<V_ID_W> indexing_v;
    unsigned char tb_index;
    unsigned char iv_pos; // Query indexing vertex position.
    unsigned int rowstart;
    unsigned int rowend;
} readmin_edge_tuple_t;

typedef struct
{
    unsigned char tb_index;
    unsigned char iv_pos; // Query indexing vertex position.
} tuplebuild_tuple_t;

typedef struct
{
    ap_uint<V_ID_W> indexing_v;
    ap_uint<V_ID_W> indexed_v;
    ap_uint<64> addr_counter;
    unsigned char tb_index;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    bool bit_last_edge;
    bool skip_counter;
    edge_flag_type flag;
} intersect_tuple_t;

typedef struct
{
    ap_uint<V_ID_W> indexing_v;
    ap_uint<V_ID_W> indexed_v;
    unsigned char tb_index;
    ap_uint<(1UL << C_W)> offset;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    bool bit_last_edge;
    edge_flag_type flag;
} offset_tuple_t;

typedef struct
{
    ap_uint<V_ID_W> indexing_v;
    ap_uint<V_ID_W> indexed_v;
    unsigned char tb_index;
    unsigned int first_block;
    unsigned int end_block;
    ap_uint<(1UL << C_W)> start_off;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    bool bit_last_edge;
    edge_flag_type flag;
} split_tuple_t;

typedef struct
{
    ap_uint<V_ID_W> indexing_v;
    ap_uint<V_ID_W> indexed_v;
    unsigned char tb_index;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    ap_uint<(1UL << C_W)> address;
    bool bit_last_address;
    bool bit_last_edge;
    edge_flag_type flag;
} verify_tuple_t;

typedef struct {
    ap_uint<V_ID_W>         indexed_v;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    bool                    bit_last_address;
    bool                    bit_last_edge;
    bool                    bit_equal;
} compact_tuple_t;

typedef struct {
    ap_uint<V_ID_W>         indexed_v;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    bool                    bit_last_edge;
    bool                    bit_checked;
} assembly_tuple_t;

/******** End tuple definition ********/
void mwj_propose(
    hls::stream<ap_uint<V_ID_W> > &stream_fifo_in,

    hls::stream<ap_uint<V_ID_W> > &stream_sol_out,
    hls::stream<bool> &stream_sol_end_out)
{
    // const ap_uint<V_ID_W> MASK_NEW_SOLUTION = ~(1UL << (V_ID_W - 1));
    const ap_uint<V_ID_W> MASK_END_EXTENSION = ~(1UL << (V_ID_W - 2));
    static ap_uint<V_ID_W> curSol_fifo[MAX_QV];
    static ap_uint<V_ID_W> curQV = 0;
    ap_uint<V_ID_W> readv;

    // Read BFS solutions
    if (stream_fifo_in.read_nb(readv)){
        if (readv.test(V_ID_W - 1))
        {
            // Delimiter, new solution
            curQV = readv.range(V_ID_W - 2, 0);
PROPOSE_FIFO_NEW_SOLUTION_LOOP:
            for (int g = 0; g < curQV; g++)
            {
#pragma HLS pipeline II = 1
                curSol_fifo[g] = stream_fifo_in.read();
            }
            curSol_fifo[curQV] = stream_fifo_in.read() & MASK_END_EXTENSION;
        }
        else
        {
            // Updating only the last vertex
            curSol_fifo[curQV] = readv & MASK_END_EXTENSION;
        }

        for (int g = 0; g < curQV + 1; g++)
        {
#pragma HLS pipeline II = 1
            stream_sol_out.write(curSol_fifo[g]);
            stream_sol_end_out.write(false);
        }
        stream_sol_end_out.write(true);
    }
}

template<size_t LKP3_HASH_W, size_t MAX_HASH_W, size_t FULL_HASH_W>
void
mwj_edgebuild(const unsigned char hash1_w,
              QueryVertex* qVertices,
              hls::stream<ap_uint<V_ID_W> >& stream_sol_in,
              hls::stream<bool>& stream_sol_end_in,
              hls::stream<bool>& stream_stop,

              hls::stream<findmin_tuple_t>& stream_tuple_out,
              hls::stream<ap_uint<V_ID_W> >& stream_sol_out,
              hls::stream<bool>& stream_sol_end_out)
{
    ap_uint<V_ID_W> curEmb[MAX_QV];
    unsigned char curQV = 0;
    bool stop, last;
    findmin_tuple_t tuple_out;

    //Initializing filter in findmin
    tuple_out.reset = true;
    tuple_out.last = false;
    tuple_out.address = 0;
    stream_tuple_out.write(tuple_out);
   
    while(true) {
        if (stream_sol_end_in.read_nb(last)){ 
            curQV = 0;

EDGEBUILD_COPYING_EMBEDDING_LOOP:
            while(!last){
                curEmb[curQV] = stream_sol_in.read();
                stream_sol_out.write(curEmb[curQV]);
                stream_sol_end_out.write(false);
                curQV++;
                last = stream_sol_end_in.read();
            }
            stream_sol_end_out.write(true);

EDGEBUILD_MAIN_LOOP:
            for(int g = 0; g < qVertices[curQV].numTablesIndexed; g++){
                unsigned char tb_index = qVertices[curQV].tables_indexed[g];
                unsigned char iv_pos = qVertices[curQV].vertex_indexing[g];

                //Computing addresses of indexed sets
                ap_uint<LKP3_HASH_W> hash_out;
                ap_uint<MAX_HASH_W> hash_trimmed;
                xf::database::details::hashlookup3_core<V_ID_W>(curEmb[iv_pos], hash_out);
                hash_trimmed = hash_out;
                hash_trimmed = hash_trimmed.range(hash1_w - 1, 0);
                unsigned int address = (tb_index * (1UL << hash1_w)) + hash_trimmed;
                tuple_out.indexing_v = curEmb[iv_pos];
                tuple_out.iv_pos = iv_pos;
                tuple_out.tb_index = tb_index;
                tuple_out.address = address;
                tuple_out.reset = false;
                tuple_out.last = (g == (qVertices[curQV].numTablesIndexed - 1));
                stream_tuple_out.write(tuple_out);
            }
            tuple_out.reset = true;
            tuple_out.last = false;
            stream_tuple_out.write(tuple_out);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

template <typename T_BLOOM,
         size_t BLOOM_LOG,
        size_t K_FUN_LOG>
unsigned int bloom_bitset(T_BLOOM filter)                         
{
    unsigned int count {0};
    for (int c = 0; c < (1UL << (BLOOM_LOG - 5)); c++){
#pragma HLS unroll
        unsigned int u = filter.range(((c + 1) * 32) - 1, c * 32);
        u = u - ((u >> 1) & 0x55555555);
        u = (u & 0x33333333) + ((u >> 2) & 0x33333333);
        count += (((u + (u >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
    }
    return (count >> K_FUN_LOG);
}


template<typename T_BLOOM, size_t BLOOM_LOG>
unsigned short
bloom_intersect(T_BLOOM &filter, T_BLOOM set_bloom)
{
#pragma HLS inline off
#pragma HLS pipeline II=1
    unsigned short bloom_s = bloom_bitset<T_BLOOM, BLOOM_LOG, 0>(set_bloom);
    filter = filter & set_bloom;
    return bloom_s;
}

template<typename T_BLOOM, size_t BLOOM_LOG, size_t K_FUN_LOG>
void
mwj_findmin(bloom_t* bloom_p,
            hls::stream<findmin_tuple_t>& stream_tuple_in,
            hls::stream<bool>& stream_stop,

            hls::stream<readmin_counter_tuple_t>& stream_tuple_out,
            hls::stream<T_BLOOM>& stream_filter_out)
{
    constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    T_BLOOM filter[K_FUN];
    readmin_counter_tuple_t tuple_out;
    findmin_tuple_t tuple_in;
    unsigned int min_size;
    bool stop;

#pragma HLS array_partition variable = filter type = complete dim = 1
#pragma HLS allocation function instances=bloom_intersect<T_BLOOM, BLOOM_LOG> limit=1

FINDMIN_TASK_LOOP:
    while (1) {
#pragma HLS pipeline II = 8

        if (stream_tuple_in.read_nb(tuple_in)) {
            unsigned int address = tuple_in.address;
            address <<= K_FUN_LOG;
            unsigned short bloom_s = 0;

            if (tuple_in.reset) {
                for (int s = 0; s < K_FUN; s++) {
#pragma HLS unroll
                    filter[s] = ~0;
                }
            } else {
                for (int s = 0; s < K_FUN; s++) {
#pragma HLS unroll
                    T_BLOOM set_bloom = bloom_p[address + s];
                    bloom_s +=
                      bloom_intersect<T_BLOOM, BLOOM_LOG>(filter[s], set_bloom);
                }
#if DEBUG_STATS
                debug::findmin_reads++;
#endif
            }
            bloom_s >>= K_FUN_LOG;

            if (tuple_in.reset) {
                min_size = ~0;
            } else {
                if (bloom_s < min_size) {
                  min_size = bloom_s;
                  tuple_out.indexing_v = tuple_in.indexing_v;
                  tuple_out.tb_index = tuple_in.tb_index;
                  tuple_out.iv_pos = tuple_in.iv_pos;
                }
            }

            if (tuple_in.last) {
                stream_tuple_out.write(tuple_out);

                for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
                  stream_filter_out.write(filter[g]);
                }
            }
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

template<typename T_BLOOM,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t FULL_HASH_W>
void
mwj_readmin_counter(const unsigned char hash1_w,
                    const unsigned char hash2_w,
                    AdjHT* hTables,
                    row_t* m_axi,
                    hls::stream<readmin_counter_tuple_t>& stream_tuple_in,
                    hls::stream<T_BLOOM>& stream_filter_in,
                    hls::stream<bool>& stream_stop,

                    hls::stream<readmin_edge_tuple_t>& stream_tuple_out,
                    hls::stream<T_BLOOM>& stream_filter_out)
{
    constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    readmin_counter_tuple_t tuple_in;
    readmin_edge_tuple_t tuple_out;
    bool stop;

READMIN_COUNTER_TASK_LOOP:
    while (true) {
        if (stream_tuple_in.read_nb(tuple_in)) {

            for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
                stream_filter_out.write(stream_filter_in.read());
            }

            ap_uint<LKP3_HASH_W> hash_out;
            ap_uint<MAX_HASH_W> hash_trimmed;
            xf::database::details::hashlookup3_core<V_ID_W>(tuple_in.indexing_v,
                                                            hash_out);
            volatile unsigned int start_off = 0;
            volatile unsigned int end_off;
            ap_uint<DDR_BIT - C_W> addr_inrow;
            ap_uint<DDR_W> ram_row;
            unsigned long addr_row;
            ap_uint<64> addr_counter;
            hash_trimmed = hash_out;
            hash_trimmed = hash_trimmed.range(hash1_w - 1, 0);

            if (hash_trimmed != 0) {
                addr_counter = hash_trimmed - 1;
                addr_counter <<= hash2_w;
                addr_counter += (1UL << hash2_w) - 1;

                /* Compute address of row storing the counter */
                addr_row = hTables[tuple_in.tb_index].start_offset +
                           (addr_counter >> (DDR_BIT - C_W));

                /* Compute address of data inside the row */
                addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);

                /* Read the data */
                ram_row = m_axi[addr_row];
                if (addr_inrow == 0) {
                  start_off = ram_row.range((1UL << C_W) - 1, 0);
                } else if (addr_inrow == 1) {
                  start_off = ram_row.range((2UL << C_W) - 1, 1UL << C_W);
                } else if (addr_inrow == 2) {
                  start_off = ram_row.range((3UL << C_W) - 1, 2UL << C_W);
                } else {
                  start_off = ram_row.range((4UL << C_W) - 1, 3UL << C_W);
                }
            }

            addr_counter = hash_trimmed;
            addr_counter <<= hash2_w;
            addr_counter += (1UL << hash2_w) - 1;

            /* Compute address of row storing the counter */
            addr_row = hTables[tuple_in.tb_index].start_offset +
                       (addr_counter >> (DDR_BIT - C_W));

            /* Compute address of data inside the row */
            addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);

            /* Read the data */
            ram_row = m_axi[addr_row];
            if (addr_inrow == 0) {
                end_off = ram_row.range((1UL << C_W) - 1, 0);
            } else if (addr_inrow == 1) {
                end_off = ram_row.range((2UL << C_W) - 1, 1UL << C_W);
            } else if (addr_inrow == 2) {
                end_off = ram_row.range((3UL << C_W) - 1, 2UL << C_W);
            } else {
                end_off = ram_row.range((4UL << C_W) - 1, 3UL << C_W);
            }

            unsigned int rowstart = hTables[tuple_in.tb_index].start_edges +
                                    (start_off >> (DDR_BIT - E_W));
            unsigned int rowend = hTables[tuple_in.tb_index].start_edges +
                                  (end_off >> (DDR_BIT - E_W));

            tuple_out.indexing_v = tuple_in.indexing_v;
            tuple_out.tb_index = tuple_in.tb_index;
            tuple_out.iv_pos = tuple_in.iv_pos;
            tuple_out.rowstart = rowstart;
            tuple_out.rowend = rowend;

            stream_tuple_out.write(tuple_out);

#if DEBUG_STATS
            debug::readmin_counter_reads += 2;
#endif
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

template <size_t W>
void hash_wrapper(ap_uint<W> key_val, ap_uint<64> &hash_val) {
#pragma HLS inline off
    xf::database::details::hashlookup3_core<W>(key_val, hash_val);
}

template<typename T_BLOOM,
         size_t BLOOM_LOG,
         size_t K_FUN,
         size_t FULL_HASH_W>
void bloom_test(T_BLOOM filter[K_FUN], ap_uint<64> hash_val, bool &test){
#pragma HLS inline off
    for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
        ap_uint<BLOOM_LOG> idx =
          hash_val.range((FULL_HASH_W / K_FUN) * (g + 1) - 1,
                         (FULL_HASH_W / K_FUN) * (g + 1) - BLOOM_LOG);
        test = test && (filter[g][idx] == 1);
    }
}

template<typename T_BLOOM,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t FULL_HASH_W>
void
mwj_readmin_edge(row_t* m_axi,
                 hls::stream<readmin_edge_tuple_t>& stream_tuple_in,
                 hls::stream<T_BLOOM>& stream_filter_in,
                 hls::stream<bool>& stream_stop,

                 hls::stream<ap_uint<V_ID_W> >& stream_set_out,
                 hls::stream<bool>& stream_set_end_out,
                 hls::stream<tuplebuild_tuple_t>& stream_tuple_out)
{
    constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    T_BLOOM filter[K_FUN];
#pragma HLS array_partition variable = filter type = complete dim = 1
    readmin_edge_tuple_t tuple_in;
    hls::stream<ap_uint<V_ID_W>, 2> hash_in_s;
    hls::stream<ap_uint<FULL_HASH_W>, 2> hash_out_s;
    tuplebuild_tuple_t tuple_out;
    ap_uint<V_ID_W> vertexCheck;
    ap_uint<V_ID_W> vertex;
    ap_uint<V_ID_W * 2> edge;
    ap_uint<FULL_HASH_W> hash_out;
    bool stop;

READMIN_EDGE_TASK_LOOP:
    while (true) {
        if (stream_tuple_in.read_nb(tuple_in)) {

            tuple_out.tb_index = tuple_in.tb_index;
            tuple_out.iv_pos = tuple_in.iv_pos;
            stream_tuple_out.write(tuple_out);

            for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
                filter[g] = stream_filter_in.read();
            }

        READMIN_EDGES_MAIN_LOOP:
            for (int g = 0; g <= tuple_in.rowend - tuple_in.rowstart; g++) {
#pragma HLS pipeline II = 2
#pragma HLS allocation function instances=hash_wrapper<V_ID_W> limit=1 
#pragma HLS allocation function instances=bloom_test<T_BLOOM, BLOOM_LOG, K_FUN, FULL_HASH_W> limit=1 
                row_t row = m_axi[tuple_in.rowstart + g];
                for (int i = 0; i < EDGE_ROW; i++) {
#pragma HLS unroll
                    edge = row.range(((i + 1) << E_W) - 1, i << E_W);
                    vertexCheck = edge.range(V_ID_W * 2 - 1, V_ID_W);
                    vertex = edge.range(V_ID_W - 1, 0);

                    if (tuple_in.indexing_v == vertexCheck) {

                      hash_wrapper<V_ID_W>(vertex, hash_out);
                      //   xf::database::details::hashlookup3_core<V_ID_W>(vertex,
                      //                                                   hash_out);
                      bool test = true;
                      bloom_test<T_BLOOM, BLOOM_LOG, K_FUN, FULL_HASH_W>(
                        filter, hash_out, test);

                      if (test) {
                        stream_set_out.write(vertex);
                        stream_set_end_out.write(false);
#if DEBUG_STATS
                        debug::readmin_vstream++;
                      } else {
                        debug::bloom_filter++;
#endif
                      }
                    }
                }
            }
#if DEBUG_STATS
            debug::readmin_edge_reads +=
              ceil(((tuple_in.rowend - tuple_in.rowstart) / 16.0));
#endif
            stream_set_end_out.write(true);
#if DEBUG_STATS
            debug::readmin_n_sets++;
#endif
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

void mwj_homomorphism(
        hls::stream< ap_uint<V_ID_W> >  &stream_set_in,
        hls::stream<bool>               &stream_set_end_in,
        hls::stream<tuplebuild_tuple_t> &stream_tuple_in,
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_in,
        hls::stream<bool>               &stream_sol_end_in,
        
        hls::stream< ap_uint<V_ID_W> >  &stream_set_out,
        hls::stream<bool>               &stream_set_end_out,
        hls::stream<tuplebuild_tuple_t> &stream_tuple_out,
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_out,
        hls::stream<bool>               &stream_sol_end_out)
{
    ap_uint<8> curQV {0};
    ap_uint<2> last_set;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool last_sol;

    last_sol = stream_sol_end_in.read();
HOMOMORPHISM_COPYING_EMBEDDING_LOOP:
    while(!last_sol){
        curEmb[curQV] = stream_sol_in.read();
        stream_sol_out.write(curEmb[curQV]);
        stream_sol_end_out.write(false);
        curQV++;
        last_sol = stream_sol_end_in.read();
    }
    stream_sol_end_out.write(true);

    stream_tuple_out.write(stream_tuple_in.read());
    last_set = stream_set_end_in.read();

HOMOMORPHISM_CHECK_LOOP:
    while(!last_set.test(0)){
        ap_uint<V_ID_W> vToVerify = stream_set_in.read();
        bool homomorphism = false;

        for (int g = 0; g < curQV; g++){
            if (vToVerify == curEmb[g])
                homomorphism = true;
        }

        if (!homomorphism){
            stream_set_out.write(vToVerify);
            stream_set_end_out.write(last_set);
        }
#if DEBUG_STATS
        else {
            debug::homomo_trashed++;
        }
#endif
        last_set = stream_set_end_in.read();
    }
    stream_set_end_out.write(last_set);
}

template <size_t MAX_BATCH_SIZE>
void mwj_batchbuild(
        hls::stream< ap_uint<V_ID_W> >  &stream_set_in,
        hls::stream<bool>               &stream_set_end_in,
        hls::stream<tuplebuild_tuple_t> &stream_tuple_in,
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_in,
        hls::stream<bool>               &stream_sol_end_in,
       
        hls::stream<bool>               &stream_req,
        hls::stream< ap_uint<V_ID_W> >  &stream_batch_out,
        hls::stream<bool>               &stream_batch_end_out,
        hls::stream<tuplebuild_tuple_t> &stream_tuple_out,
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_out,
        hls::stream<bool>               &stream_sol_end_out)
{
    ap_uint<8> curQV {0};
    ap_uint<V_ID_W> curEmb[MAX_QV];
    unsigned int batch_counter {0};
    bool last;
    
    last = stream_sol_end_in.read();
BATCHBUILD_COPYING_EMBEDDING_LOOP:
    while(!last){
        curEmb[curQV] = stream_sol_in.read();
        stream_sol_out.write(curEmb[curQV]);
        stream_sol_end_out.write(false);
        curQV++;
        last = stream_sol_end_in.read();
    }
    stream_sol_end_out.write(true);

    tuplebuild_tuple_t tuple_in = stream_tuple_in.read();
    stream_tuple_out.write(tuple_in);
    last = stream_set_end_in.read();

BATCHBUILD_MAIN_LOOP:
    while(!last){

BATCHBUILD_MOVING_SET_LOOP:
        while(!last && batch_counter != (MAX_BATCH_SIZE - 1)){
#pragma HLS pipeline II=1
            ap_uint<V_ID_W> node = stream_set_in.read();
            stream_batch_out.write(node);
            stream_batch_end_out.write(false);
            batch_counter++;
            last = stream_set_end_in.read();
        }

        // Stream again partial solution if max batch size is reached
        if (!last && batch_counter == (MAX_BATCH_SIZE - 1)){
            stream_req.write(true);
            stream_batch_end_out.write(true);
            batch_counter = 0;

BATCHBUILD_ADD_SOLUTION_LOOP:
            for (int g = 0; g < curQV; g++){
                stream_sol_out.write(curEmb[g]);
                stream_sol_end_out.write(false);
            }
            stream_sol_end_out.write(true);
            stream_tuple_out.write(tuple_in);
        }
    }
    stream_batch_end_out.write(true);
}

template<size_t BATCH_SIZE_LOG,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t FULL_HASH_W>
void
mwj_tuplebuild(const unsigned char hash1_w,
               const unsigned char hash2_w,
               QueryVertex* qVertices,
               hls::stream<ap_uint<V_ID_W> >& stream_set_in,
               hls::stream<bool>& stream_set_end_in,
               hls::stream<tuplebuild_tuple_t>& stream_tuple_in,
               hls::stream<ap_uint<V_ID_W> >& stream_sol_in,
               hls::stream<bool>& stream_sol_end_in,
               hls::stream<bool>& stream_stop,

               hls::stream<intersect_tuple_t>& stream_tuple_out,
               hls::stream<bool>& stream_tuple_end_out,
               hls::stream<ap_uint<V_ID_W> >& stream_sol_out,
               hls::stream<bool>& stream_sol_end_out)
{
    ap_uint<V_ID_W> buffer[(1UL << BATCH_SIZE_LOG)];
    ap_uint<BATCH_SIZE_LOG> buffer_size {0};
    ap_uint<BATCH_SIZE_LOG> buffer_p;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<V_ID_W> vToVerify;
    hls::stream<ap_uint<V_ID_W>, 4> hash_in0, hash_in1;
    hls::stream<ap_uint<LKP3_HASH_W>, 4> hash_out0, hash_out1;
    unsigned long addr_counter;
    intersect_tuple_t tuple_out;
    unsigned char curQV {0};
    bool stop, last;
   
    while(true) {
        if (stream_sol_end_in.read_nb(last)){ 
            curQV = 0;
            buffer_size = 0;

TUPLEBUILD_COPYING_EMBEDDING_LOOP:
            while(!last){
                curEmb[curQV] = stream_sol_in.read();
                stream_sol_out.write(curEmb[curQV]);
                stream_sol_end_out.write(false);
                curQV++;
                last = stream_sol_end_in.read();
            }
            stream_sol_end_out.write(true);
            
            tuplebuild_tuple_t tuple_in = stream_tuple_in.read();
            unsigned char cycles = qVertices[curQV].numTablesIndexed;

            last = stream_set_end_in.read();

TUPLEBUILD_SAVE_NODE_LOOP:
            while (!last) {
#pragma HLS pipeline II = 1
                vToVerify = stream_set_in.read();
                buffer[buffer_size++] = vToVerify;
                last = stream_set_end_in.read();
            }

TUPLEBUILD_EDGE_LOOP:
            for (int g = 0; g < cycles; g++){
#pragma HLS pipeline II=2
                uint8_t tableIndex = qVertices[curQV].tables_indexed[g];
                uint8_t ivPos = qVertices[curQV].vertex_indexing[g];
                bool bit_last = (g == cycles - 1);
                bool bit_min = (tuple_in.tb_index == tableIndex && 
                        tuple_in.iv_pos == ivPos);

                for (int buffer_p = 0; buffer_p < buffer_size; buffer_p++){
#pragma HLS loop_flatten
                    vToVerify = buffer[buffer_p];

                    hash_in0.write(vToVerify);
                    hash_in1.write(curEmb[ivPos]);
                    xf::database::hashLookup3<V_ID_W>(hash_in0, hash_out0);
                    xf::database::hashLookup3<V_ID_W>(hash_in1, hash_out1);
                    ap_uint<MAX_HASH_W> indexed_h = hash_out0.read();
                    ap_uint<MAX_HASH_W> indexing_h = hash_out1.read();
                    
                    addr_counter = indexing_h.range(hash1_w - 1, 0);
                    addr_counter <<= hash2_w;
                    addr_counter += indexed_h.range(hash2_w - 1, 0);
                    
                    tuple_out.indexed_v = vToVerify;
                    tuple_out.indexing_v = curEmb[ivPos];
                    tuple_out.addr_counter = addr_counter;
                    tuple_out.tb_index = tableIndex;
                    tuple_out.pos = buffer_p;
                    tuple_out.bit_last_edge = bit_last;
                    tuple_out.flag = (bit_min)? MIN_SET: CHECK;
                    tuple_out.skip_counter = false;
                    
                    stream_tuple_out.write(tuple_out);
                    stream_tuple_end_out.write(false);
                    
                    if (addr_counter == 0)
                        tuple_out.skip_counter = true; 
                    tuple_out.addr_counter = addr_counter - 1;
                    stream_tuple_out.write(tuple_out);
                    stream_tuple_end_out.write(false);

                    // std::cout << "( "
                    //           << tuple_out.indexing_v
                    //           << ", " << tuple_out.indexed_v
                    //           << ", " << tuple_out.addr_counter
                    //           << ", " << (int)tuple_out.tb_index
                    //           << ", " << tuple_out.pos
                    //           << ", " << tuple_out.bit_last_edge
                    //           << ", " << tuple_out.bit_min_set << ")\n" << std::endl;

                }
            }
            stream_tuple_end_out.write(true);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

void mwj_bypass_sol(
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_in,
        hls::stream<bool>               &stream_sol_end_in,

        hls::stream< ap_uint<V_ID_W> >  &stream_sol_out,
        hls::stream<bool>               &stream_sol_end_out)
{
#pragma HLS pipeline II=1
    bool last = stream_sol_end_in.read();
    if (!last){
        stream_sol_out.write(stream_sol_in.read());
        stream_sol_end_out.write(false);
    } else {
        stream_sol_end_out.write(true);
    }
}

template <size_t BATCH_SIZE_LOG>
void mwj_intersect(
    AdjHT *hTables,
    htb_cache_t &htb_buf,
    hls::stream<intersect_tuple_t> &stream_tuple_in,
    hls::stream<bool> &stream_tuple_end_in,
    hls::stream<bool> &stream_stop,

    hls::stream<offset_tuple_t> &stream_tuple_out,
    hls::stream<bool> &stream_tuple_end_out)
{
    ap_uint<V_ID_W> indexing_v;
    intersect_tuple_t tuple_in;
    offset_tuple_t tuple_out;
    unsigned char tableIndex;
    ap_uint<DDR_W> ram_row;
    unsigned long addr_row;
    ap_uint<DDR_BIT - C_W> addr_inrow;
    ap_uint<64> addr_counter;
    bool stop, last;

INTERSECT_TASK_LOOP:
    while (1) {
#pragma HLS pipeline II=1
        if (stream_tuple_end_in.read_nb(last)){

            if (!last){
                tuple_in = stream_tuple_in.read();
                ap_uint<(1UL << C_W)> offset = 0;
                tableIndex = tuple_in.tb_index;
                addr_counter = tuple_in.addr_counter;

                if (tuple_in.flag == CHECK && !tuple_in.skip_counter){
                    
                    /* Compute address of row storing the counter */
                    addr_row = hTables[tableIndex].start_offset + 
                        (addr_counter >> (DDR_BIT - C_W));

                    /* Compute address of data inside the row */
                    addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);
                    ram_row = htb_buf.get(addr_row, 1);
                    if (addr_inrow == 0) {
                        offset = ram_row.range((1UL << C_W) - 1, 0);
                    } else if (addr_inrow == 1) {
                        offset =
                          ram_row.range((2UL << C_W) - 1, 1UL << C_W);
                    } else if (addr_inrow == 2) {
                        offset =
                          ram_row.range((3UL << C_W) - 1, 2UL << C_W);
                    } else {
                        offset =
                          ram_row.range((4UL << C_W) - 1, 3UL << C_W);
                    }

#if DEBUG_STATS
                    debug::intersect_reads += 1;
#endif /* DEBUG_STATS */
                }

                tuple_out.indexed_v     = tuple_in.indexed_v;
                tuple_out.indexing_v    = tuple_in.indexing_v;
                tuple_out.tb_index      = tuple_in.tb_index;
                tuple_out.pos           = tuple_in.pos;
                tuple_out.offset        = offset;
                tuple_out.bit_last_edge = tuple_in.bit_last_edge;
                tuple_out.flag          = tuple_in.flag;

                // std::cout << "( "
                //           << tuple_out.indexing_v
                //           << ", " << tuple_out.indexed_v
                //           << ", " << (int)tuple_out.tb_index
                //           << ", " << tuple_out.pos
                //           << ", " << tuple_out.offset
                //           << ", " << tuple_out.bit_last_edge
                //           << ", " << tuple_out.bit_min_set << ")" << std::endl;
                /* bits[pos] = bits[pos] & (start_off < end_off); */
                stream_tuple_out.write(tuple_out);
            }
            stream_tuple_end_out.write(last);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

void mwj_offset(
    hls::stream<offset_tuple_t> &stream_tuple_in,
    hls::stream<bool> &stream_tuple_end_in,

    hls::stream<split_tuple_t> &stream_tuple_out,
    hls::stream<bool> &stream_tuple_end_out)
{
#pragma HLS pipeline II=2 style=flp
    constexpr size_t EDGE_BLOCK = (CACHE_WORDS_PER_LINE + DDR_BIT - E_W);
    offset_tuple_t tuple_in;
    split_tuple_t tuple_out;

    bool last = stream_tuple_end_in.read();
    // if (stream_tuple_end_in.read_nb(last)){
    if (!last)
    {
        tuple_in = stream_tuple_in.read();
        ap_uint<(1UL << C_W)> end_off = tuple_in.offset;
        tuple_out.indexed_v = tuple_in.indexed_v;
        tuple_out.indexing_v = tuple_in.indexing_v;
        tuple_out.pos = tuple_in.pos;
        tuple_out.tb_index = tuple_in.tb_index;
        tuple_out.bit_last_edge = tuple_in.bit_last_edge;

        last = stream_tuple_end_in.read();
        tuple_in = stream_tuple_in.read();

        tuple_out.start_off = tuple_in.offset;
        tuple_out.first_block = tuple_in.offset >> EDGE_BLOCK;
        tuple_out.end_block = (tuple_in.offset == end_off) ? (unsigned int)(end_off >> EDGE_BLOCK) : (unsigned int)((end_off - 1) >> EDGE_BLOCK);
        tuple_out.flag = (tuple_in.flag == MIN_SET) ? MIN_SET : ((tuple_in.offset == end_off) ? NO_EDGE : CHECK);

#if DEBUG_STATS
        debug::intersect_filter += (tuple_out.flag == NO_EDGE)? 1 : 0;
#endif /* DEBUG_STATS */

        stream_tuple_out.write(tuple_out);
    }
    stream_tuple_end_out.write(last);
    // }
}

template <size_t MAX_BATCH_SIZE>
void mwj_split(
    hls::stream<split_tuple_t> &stream_tuple_in,
    hls::stream<bool> &stream_tuple_end_in,

    hls::stream<verify_tuple_t> &stream_tuple_out,
    hls::stream<bool> &stream_tuple_end_out)
{
    constexpr size_t EDGE_BLOCK = (CACHE_WORDS_PER_LINE + DDR_BIT - E_W);
    split_tuple_t tuple_in;
    verify_tuple_t tuple_out;
    bool last;

    last = stream_tuple_end_in.read();
    if (!last){
        tuple_in = stream_tuple_in.read();
        
        tuple_out.indexing_v    = tuple_in.indexing_v;
        tuple_out.indexed_v     = tuple_in.indexed_v;
        tuple_out.tb_index      = tuple_in.tb_index;
        tuple_out.pos           = tuple_in.pos;
        tuple_out.bit_last_edge = tuple_in.bit_last_edge;
        tuple_out.flag          = tuple_in.flag;

        // First loop iteration 
        tuple_out.address = tuple_in.start_off;
        tuple_out.bit_last_address = (tuple_in.first_block == tuple_in.end_block);
        stream_tuple_out.write(tuple_out);
        stream_tuple_end_out.write(false);
        tuple_in.first_block++;
        
SPLIT_MAIN_LOOP:
        for (unsigned short s = (1UL << EDGE_BLOCK);
             tuple_in.first_block <= tuple_in.end_block;
             tuple_in.first_block++, s += (1UL << EDGE_BLOCK))
        {
#pragma HLS pipeline II=1
            tuple_out.address = tuple_in.start_off + s;
            tuple_out.bit_last_address = (tuple_in.first_block == tuple_in.end_block);
            stream_tuple_out.write(tuple_out);
            stream_tuple_end_out.write(false);

            // std::cout << "( "
            //           << tuple_out.indexing_v
            //           << ", " << tuple_out.indexed_v
            //           << ", " << (int)tuple_out.tb_index
            //           << ", " << tuple_out.pos
            //           << ", " << tuple_out.address
            //           << ", " << tuple_out.bit_last_address
            //           << ", " << tuple_out.bit_last_edge
            //           << ", " << tuple_out.bit_no_edge
            //           << ", " << tuple_out.bit_min_set << ")" << std::endl;
        }
    } else {
        stream_tuple_end_out.write(true);
    }
}

template <size_t BATCH_SIZE_LOG>
void mwj_verify(
        AdjHT *hTables,
        htb_cache_t                  &htb_buf,
        hls::stream<verify_tuple_t>  &stream_tuple_in,
        hls::stream<bool>            &stream_tuple_end_in,
        hls::stream<bool>            &stream_stop,

        hls::stream<compact_tuple_t> &stream_tuple_out,
        hls::stream<bool>            &stream_tuple_end_out)
{
    constexpr size_t EDGE_PER_WORD = (DDR_BIT - E_W);
    
    verify_tuple_t tuple_in;
    compact_tuple_t tuple_out;
    ap_uint<V_ID_W> candidate_v;
    ap_uint<V_ID_W> indexing_v;
    unsigned char tableIndex;
    bool stop, last;
    ap_uint<128> edge_block[(1UL << CACHE_WORDS_PER_LINE)];
#pragma HLS array_partition variable=edge_block type=complete

VERIFY_TASK_LOOP:
    while(1){
#pragma HLS pipeline II=1
        if (stream_tuple_end_in.read_nb(last)){

            if(!last){
                tuple_in = stream_tuple_in.read();
                candidate_v = tuple_in.indexed_v;
                tuple_out.bit_equal = (tuple_in.flag == MIN_SET);

                if (tuple_in.flag == CHECK){
                    indexing_v = tuple_in.indexing_v;
                    tableIndex = tuple_in.tb_index;

                    // 128 bit word address
                    unsigned long addr_row = hTables[tableIndex].start_edges +
                        (tuple_in.address >> EDGE_PER_WORD);

                    // Read the data
                    htb_buf.get_line(addr_row, 0, edge_block);
                    ap_uint<(1UL << E_W)> edge;
                    edge.range(V_ID_W - 1, 0) = candidate_v;
                    edge.range(2 * V_ID_W - 1, V_ID_W) = indexing_v;

                    for (int g = 0; g < (1UL << CACHE_WORDS_PER_LINE); g++){
#pragma HLS unroll
                        for (int s = 0; s < (1UL << (EDGE_PER_WORD)); s++){
#pragma HLS unroll
                            if (edge == edge_block[g].range(((s + 1) << E_W) - 1, s << E_W))
                                tuple_out.bit_equal = true;
                        }
                    }

#if DEBUG_STATS
                    debug::verify_reads++;
#endif
                }

                tuple_out.bit_last_edge     = tuple_in.bit_last_edge;
                tuple_out.bit_last_address  = tuple_in.bit_last_address;
                tuple_out.indexed_v         = candidate_v;
                tuple_out.pos               = tuple_in.pos;
                stream_tuple_out.write(tuple_out);
            }
            stream_tuple_end_out.write(last);
        }

        if (stream_stop.read_nb(stop)){
            break;
        }
    }
}

void mwj_compact(
        hls::stream<compact_tuple_t>  &stream_tuple_in,
        hls::stream<bool>             &stream_tuple_end_in,
        
        hls::stream<assembly_tuple_t> &stream_tuple_out,
        hls::stream<bool>             &stream_tuple_end_out)
{
#pragma HLS pipeline II=1
    static bool checked = false;
    
    compact_tuple_t tuple_in;
    assembly_tuple_t tuple_out;
    bool last;

    last = stream_tuple_end_in.read();
    if (!last) {
        tuple_in = stream_tuple_in.read();
        checked |= tuple_in.bit_equal;
        if (tuple_in.bit_last_address){
            tuple_out.bit_checked = checked;
            tuple_out.indexed_v = tuple_in.indexed_v;
            tuple_out.pos = tuple_in.pos;
            tuple_out.bit_last_edge = tuple_in.bit_last_edge;
            stream_tuple_out.write(tuple_out); 
            stream_tuple_end_out.write(false);
#if DEBUG_STATS
            if (!checked)
                debug::verify_filter++;
#endif /* DEBUG_STATS */
            checked = false;
        }
    } else {
        stream_tuple_end_out.write(true);
        checked = false; 
    }
}

template<size_t MAX_BATCH_SIZE>
void mwj_filter(
        hls::stream<assembly_tuple_t>   &stream_tuple_in,
        hls::stream<bool>               &stream_tuple_end_in,

        hls::stream< ap_uint<V_ID_W> >  &stream_set_out,
        hls::stream<bool>               &stream_set_end_out)
{
#pragma HLS pipeline II=1
    static ap_uint<MAX_BATCH_SIZE> bits = ~0;
    assembly_tuple_t tuple_in;
    bool last;
    
    last = stream_tuple_end_in.read();
    if (!last){
        tuple_in = stream_tuple_in.read();
        unsigned short p = tuple_in.pos;
        bits[p] = bits[p] && tuple_in.bit_checked;
        if (tuple_in.bit_last_edge && bits.test(p)){
            stream_set_out.write(tuple_in.indexed_v);
            stream_set_end_out.write(false);
        }
    } else {
        bits = ~0;
        stream_set_end_out.write(true);
    }
}

void mwj_assembly(
    unsigned short nQueryVer,
    hls::stream<ap_uint<V_ID_W> > &stream_inter_in,
    hls::stream<bool> &stream_end_inter_in,
    hls::stream<ap_uint<V_ID_W> > &stream_embed_in,
    hls::stream<bool> &stream_end_embed_in,
    hls::stream<ap_uint<V_ID_W> > &stream_batch,
    hls::stream<bool> &stream_batch_end,

    hls::stream<bool> streams_stop[STOP_S],
    hls::stream<bool> &stream_sol0,
    hls::stream<bool> &stream_req,
    hls::stream<ap_uint<V_ID_W> > &stream_partial_out,
#if COUNT_ONLY
    long unsigned int &result
#else
    hls::stream<T_NODE> &result
#endif
)
{
    const ap_uint<V_ID_W> MASK_NEW_SOLUTION = (1UL << (V_ID_W - 1));
    const ap_uint<V_ID_W> MASK_END_EXTENSION = (1UL << (V_ID_W - 2));
    ap_uint<V_ID_W> curQV;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool last_sol, last_set, stop, last_start;
    bool token_new_start;
    T_NODE node;

#if COUNT_ONLY
    unsigned long int counter {0};
#endif

    last_start = stream_batch_end.read();
    last_start = stream_batch_end.read();
    token_new_start = false;
    stream_partial_out.write(0 | MASK_NEW_SOLUTION);
    stream_partial_out.write(stream_batch.read() | MASK_END_EXTENSION);
    stream_req.write(true);
ASSEMBLY_TASK_LOOP:
    while(1) {
        if (stream_end_embed_in.read_nb(last_sol)){

            curQV = 0;
ASSEMBLY_COPYING_EMBEDDING_LOOP:
            while(!last_sol){
                curEmb[curQV] = stream_embed_in.read();
                curQV++;
                last_sol = stream_end_embed_in.read();
            }

            last_set = stream_end_inter_in.read();

            //Last bit used to signal new solution
            if (!last_set && (curQV != nQueryVer - 1)){
                stream_partial_out.write(curQV | MASK_NEW_SOLUTION);
ASSEMBLY_RADIX_LOOP:
                for (int g = 0; g < curQV; g++){
                    stream_partial_out.write(curEmb[g]);
                }
            }

ASSEMBLY_SET_LOOP:
            while(!last_set){
                last_set = stream_end_inter_in.read();
                ap_uint<V_ID_W> vToVerify = stream_inter_in.read();
                
                /* Write in the correct stream */
                if (curQV == nQueryVer - 1){
#if COUNT_ONLY
                    if (!last_start && token_new_start){
                        last_start = stream_batch_end.read();
                        token_new_start = false;
                        stream_partial_out.write(0 | MASK_NEW_SOLUTION);
                        stream_partial_out.write(stream_batch.read() | MASK_END_EXTENSION);
                        stream_req.write(true);
                    }
                    // for (int g = 0; g < nQueryVer - 1; g++){
                    //     f << curEmb[g] << " ";
                    // }
                    // f << vToVerify << std::endl;
                    counter++;
#else
ASSEMBLY_WRITE_FINAL_LOOP:
                    for (int g = 0; g < curQV; g++){
                        node.data = curEmb[g];
                        node.last = false;
                        node.keep = ~0;
                        result.write(node);
                    }
                    node.data = vToVerify;
                    node.last = false;
                    node.keep = ~0;
                    result.write(node);
#endif
                } else {
                    token_new_start = true;
                    stream_partial_out.write((last_set)? (vToVerify | MASK_END_EXTENSION) : vToVerify);
                    stream_req.write(true);
                }
            }

            // Last batch of a set 
            stream_req.write(false);

        }

        if (stream_sol0.read_nb(stop)){

#if DEBUG_STATS
                debug::miss_indexing++;
#endif

            // Test if there are some node from start batch 
            if (!last_start){
                last_start = stream_batch_end.read();
                token_new_start = true;
                stream_partial_out.write(0 | MASK_NEW_SOLUTION);
                stream_partial_out.write(stream_batch.read() | MASK_END_EXTENSION);
                stream_req.write(true);
            } else {
                break;
            }
        }
    }

#if COUNT_ONLY
    /* Write in output number of results */
    result += counter;
#else
    /* Write last node */
    node.data = 0;
    node.last = true;
    node.keep = ~0;
    result.write(node);
#endif
        
    for (int g = 0; g < STOP_S; g++){
#pragma HLS unroll
        streams_stop[g].write(true);
    }
}

template<typename T>
void mwj_merge(
        hls::stream<T> in[MERGE_IN_STREAMS],
        hls::stream<T> &out)
{
    T data;
    for(int g = 0; g < MERGE_IN_STREAMS; g++){
#pragma HLS unroll
        if (in[g].read_nb(data))
            out.write(data);
    }
}

void mwj_stop(
        hls::stream<bool> &stream_req,
        hls::stream<bool> &dynfifo_overflow,
        hls::stream<bool> &stream_sol0)
{
    static unsigned long sol {0};
    bool ovf {false};
    bool req = stream_req.read();
    
    if (req) {
        sol++;
    } else {
        sol--;
    }

    dynfifo_overflow.read_nb(ovf);

    if (sol == 0 || ovf == true){
        sol = 0;
        stream_sol0.write(true);
    }
}

template<size_t LKP3_HASH_W, size_t MAX_HASH_W, size_t FULL_HASH_W>
void
mwj_batch(const unsigned char hash1_w,
          AdjHT* hTables,
          QueryVertex* qVertices,
          ap_uint<DDR_W>* htb_buf,

          hls::stream<bool>& stream_batch_end,
          hls::stream<ap_uint<V_ID_W> >& stream_batch)
{
    ap_uint<8> tableIndex {0};
    ap_uint<32> minSize = (1UL << 32) - 1;
    ap_uint<32> minStart, minOff;
    ap_uint<8 + 1 + V_ID_W> minData;
    ap_uint<V_ID_W*2> edge;
    ap_uint<V_ID_W> vertex, vertexCheck;
    ap_uint<V_ID_W> set[MAX_CL];
    ap_uint<MAX_HASH_W> hash_buff, hash_new;
    unsigned char set_counter = 0;
    bool flag_buff = false;
    bool flag_new = true;
    hash_buff = hash_new = 0;
    unsigned int rm_start = 0;

PROPOSE_TBINDEXING_LOOP:
    for(int g = 0; g < qVertices[0].numTablesIndexing; g++){
        tableIndex = qVertices[0].tables_indexing[g];

        if (hTables[tableIndex].n_edges < minSize){
            minSize = hTables[tableIndex].n_edges;
            minOff = hTables[tableIndex].start_edges;
            minStart = 0;
            minData.range(7, 0) = tableIndex;
            minData.clear(8);
            minData.range(V_ID_W + 8, 9) = 0;
        }
    }
    
    unsigned int rowstart = minOff + (minStart >> (DDR_BIT - E_W));
    unsigned int rowend = minOff + (minSize >> (DDR_BIT - E_W));
    unsigned int window_right = minSize.range((DDR_BIT - E_W) - 1, 0) + 
        (rowend - rowstart) * EDGE_ROW;
    unsigned int cnt = 0;

PROPOSE_READ_MIN_INDEXING_LOOP:
    for (unsigned int g = 0; g <= rowend - rowstart; g++){
        row_t row = htb_buf[rowstart + g];
        for (unsigned int i = 0; i < EDGE_ROW; i++, cnt++){
#pragma HLS unroll
            if (cnt < window_right){
                edge = row.range((1UL << E_W) - 1, 0);
                vertex = edge.range(V_ID_W * 2 - 1, V_ID_W);
                ap_uint<LKP3_HASH_W> hash_out;
                xf::database::details::hashlookup3_core<V_ID_W>(
                        vertex,
                        hash_out);
                hash_new = hash_out.range(MAX_HASH_W - 1, 0);
                hash_new = hash_new.range(hash1_w - 1, 0);

                if (flag_buff && hash_buff == hash_new){
                    flag_new = true;
EXTRACT_BAGTOSET_SETCHECKER_LOOP:
                    for(int nSet = 0; nSet < set_counter; nSet++){
                        if (vertex == set[nSet]){
                            flag_new = false;
                            break;
                        }
                    }
                } else {
                    flag_new = true;
                    set_counter = 0;
                }

#ifndef __SYNTHESIS__
                assert(set_counter < MAX_CL);
#endif
                if (flag_new) {
                    set[set_counter++] = vertex;
#if DEBUG_STATS
                    debug::start_set++;
#endif
                    // unsigned int start_off = 0;
                    // unsigned int end_off;
                    // ap_uint<DDR_BIT - C_W> addr_inrow;
                    // ap_uint<DDR_W> ram_row;
                    // unsigned long addr_row;
                    // ap_uint<64> addr_counter;
                    // ap_uint<V_ID_W * 2> edge2;
                    // ap_uint<V_ID_W> vertex2;
                    // bool verified = true;

                    // for (int t = 0; t < qVertices[0].numTablesIndexing; t++) {
                    //     unsigned char tb_index = qVertices[0].tables_indexing[t];
                    //     bool tb_verified = false;

                    //     if (hash_buff != 0) {
                    //         addr_counter = hash_new - 1;
                    //         addr_counter <<= 5;
                    //         addr_counter += (1UL << 5) - 1;

                    //         /* Compute address of row storing the counter */
                    //         addr_row = hTables[tb_index].start_offset +
                    //                    (addr_counter >> (DDR_BIT - C_W));

                    //         /* Compute address of data inside the row */
                    //         addr_inrow =
                    //           addr_counter.range((DDR_BIT - C_W) - 1, 0);

                    //         /* Read the data */
                    //         ram_row = htb_buf[addr_row];
                    //         if (addr_inrow == 0) {
                    //             start_off = ram_row.range((1UL << C_W) - 1, 0);
                    //         } else if (addr_inrow == 1) {
                    //             start_off =
                    //               ram_row.range((2UL << C_W) - 1, 1UL << C_W);
                    //         } else if (addr_inrow == 2) {
                    //             start_off =
                    //               ram_row.range((3UL << C_W) - 1, 2UL << C_W);
                    //         } else {
                    //             start_off =
                    //               ram_row.range((4UL << C_W) - 1, 3UL << C_W);
                    //         }
                    //     }

                    //     addr_counter = hash_new;
                    //     addr_counter <<= 5;
                    //     addr_counter += (1UL << 5) - 1;

                    //     /* Compute address of row storing the counter */
                    //     addr_row = hTables[tb_index].start_offset +
                    //                (addr_counter >> (DDR_BIT - C_W));

                    //     /* Compute address of data inside the row */
                    //     addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);

                    //     /* Read the data */
                    //     ram_row = htb_buf[addr_row];
                    //     if (addr_inrow == 0) {
                    //         end_off = ram_row.range((1UL << C_W) - 1, 0);
                    //     } else if (addr_inrow == 1) {
                    //         end_off =
                    //           ram_row.range((2UL << C_W) - 1, 1UL << C_W);
                    //     } else if (addr_inrow == 2) {
                    //         end_off =
                    //           ram_row.range((3UL << C_W) - 1, 2UL << C_W);
                    //     } else {
                    //         end_off =
                    //           ram_row.range((4UL << C_W) - 1, 3UL << C_W);
                    //     }

                    //     unsigned int rowstart2 = hTables[tb_index].start_edges +
                    //                              (start_off >> (DDR_BIT - E_W));
                    //     unsigned int rowend2 = hTables[tb_index].start_edges +
                    //                            (end_off >> (DDR_BIT - E_W));

                    //     for (unsigned int s = 0; s <= rowend2 - rowstart2;
                    //          s++) {
                    //         row_t rowss = htb_buf[rowstart2 + s];
                    //         for (unsigned int i2 = 0; i2 < E_W; i2++) {
                    //             edge2 = rowss.range((1UL << E_W) - 1, 0);
                    //             vertex2 = edge2.range(V_ID_W * 2 - 1, V_ID_W);
                    //             if (vertex == vertex2) {
                    //               tb_verified = true;
                    //             }
                    //             rowss >>= (1UL << E_W);
                    //         }
                    //     }
                    //     verified = verified && tb_verified;

                    // }

                    // if (verified) {
                    //     stream_batch_end.write(false);
                    //     stream_batch.write(vertex);
                    // } else {
                    //     rm_start++;
                    // }
                    stream_batch_end.write(false);
                    stream_batch.write(vertex);
                }

                hash_buff = hash_new;
                flag_buff = true;
            }
            row >>= (1UL << E_W);
        }
    }

#if DEBUG_STATS
    debug::batch_reads += ceil((rowend - rowstart) / 16.0);
    std::cout << rm_start << " removed\n";
#endif
    stream_batch_end.write(true);
}

// template <typename T_BLOOM,
//           size_t BLOOM_LOG,
//           size_t K_FUN_LOG>
// void readmincache_wrapper(
//     AdjHT *hTables,
//     htb_cache_t &htb_buf,
//     hls::stream<readmin_tuple_t> &stream_tuple_in,
//     hls::stream<T_BLOOM> &stream_filter_in,
//     hls::stream<bool> &stream_stop,

//     hls::stream< ap_uint<V_ID_W> > &stream_set_out,
//     hls::stream<bool> &stream_set_end_out,
//     hls::stream<tuplebuild_tuple_t> &stream_tuple_out)
// {
//     htb_buf.init(2);
//     mwj_readmin<T_BLOOM,
//                 BLOOM_LOG,
//                 K_FUN_LOG>
//                 (
//         hTables,
//         htb_buf,
//         stream_tuple_in,
//         stream_filter_in,
//         stream_stop,
//         stream_set_out,
//         stream_set_end_out,
//         stream_tuple_out);
// }

template <size_t BATCH_SIZE_LOG>
void intersectcache_wrapper(
    AdjHT *hTables,
    htb_cache_t &htb_buf,
    hls::stream<intersect_tuple_t> &stream_tuple_in,
    hls::stream<bool> &stream_tuple_end_in,
    hls::stream<bool> &stream_stop,

    hls::stream<offset_tuple_t> &stream_tuple_out,
    hls::stream<bool> &stream_tuple_end_out)
{
    htb_buf.init(1);
    mwj_intersect<BATCH_SIZE_LOG>(
        hTables,
        htb_buf,
        stream_tuple_in,
        stream_tuple_end_in,
        stream_stop,
        stream_tuple_out,
        stream_tuple_end_out);
}

template <size_t BATCH_SIZE_LOG>
void verifycache_wrapper(
        AdjHT *hTables,
        htb_cache_t                  &htb_buf,
        hls::stream<verify_tuple_t>  &stream_tuple_in,
        hls::stream<bool>            &stream_tuple_end_in,
        hls::stream<bool>            &stream_stop,

        hls::stream<compact_tuple_t> &stream_tuple_out,
        hls::stream<bool>            &stream_tuple_end_out)
{
    htb_buf.init(0);
    mwj_verify<BATCH_SIZE_LOG>(
            hTables,
            htb_buf,
            stream_tuple_in,
            stream_tuple_end_in,          
            stream_stop,                  
            stream_tuple_out,             
            stream_tuple_end_out);
    htb_buf.stop();
}

template <typename T_BLOOM,
          size_t BLOOM_LOG,
          size_t K_FUN_LOG,
          size_t LKP3_HASH_W,
          size_t MAX_HASH_W,
          size_t FULL_HASH_W>
void multiwayJoin(
    ap_uint<DDR_W> *htb_buf0,
    ap_uint<DDR_W> *htb_buf1,
    ap_uint<DDR_W> *htb_buf2,
    T_BLOOM *bloom_p,
    row_t *res_buf,
    AdjHT *hTables0,
    AdjHT *hTables1,
    QueryVertex *qVertices0,
    const unsigned short nQueryVer,
    const unsigned char hash1_w,
    const unsigned char hash2_w,
    const unsigned long dynfifo_space,
    hls::stream< ap_uint<V_ID_W> > &stream_batch,
    hls::stream<bool> &stream_batch_end,
    unsigned long &dynfifo_diagnostic,

#if COUNT_ONLY
    long unsigned int &result
#else
    hls::stream<T_NODE> &result
#endif
)
{
#pragma HLS STABLE variable=htb_buf0
#pragma HLS STABLE variable=htb_buf1
#pragma HLS STABLE variable=htb_buf2
#pragma HLS STABLE variable=bloom_p
#pragma HLS STABLE variable=hTables0
#pragma HLS STABLE variable=hTables1
#pragma HLS STABLE variable=qVertices0
#pragma HLS STABLE variable=nQueryVer
#pragma HLS STABLE variable=hash1_w
#pragma HLS STABLE variable=hash2_w
#pragma HLS STABLE variable=dynfifo_space

#pragma HLS DATAFLOW

    /* Propose data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> p0_stream_sol
        ("Propose - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> p0_stream_sol_end
        ("Propose - partial solution end flag");
    
    /* Edgebuild data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> e_stream_sol
        ("Edgebuild - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> e_stream_sol_end
        ("Edgebuild - partial solution end flag");
    hls_thread_local hls::stream<findmin_tuple_t, S_D> e_stream_tuple
        ("Edgebuild - tuples");
    
    /* Findmin data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> p_stream_sol
        ("Propose - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> p_stream_sol_end
        ("Propose - partial solution end flag");
    hls_thread_local hls::stream<T_BLOOM, (1UL << K_FUN_LOG) * 4> p_stream_filter
        ("Propose - filter");
    hls_thread_local hls::stream<readmin_counter_tuple_t, S_D> p_stream_tuple
        ("Propose - tuples");

    /* Readmin counter data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> rc_stream_sol
        ("Readmin counter - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> rc_stream_sol_end
        ("Readmin counter - partial solution end flag");
    hls_thread_local hls::stream<T_BLOOM, (1UL << K_FUN_LOG) * 4> rc_stream_filter
        ("Readmin counter - filter");
    hls_thread_local hls::stream<readmin_edge_tuple_t, S_D> rc_stream_tuple
        ("Readmin counter - tuples");
    
    /* Readmin edge data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> re_stream_sol
        ("Readmin edge - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> re_stream_sol_end
        ("Readmin edge - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D + 32> re_stream_set
        ("Readmin edge - set nodes");
    hls_thread_local hls::stream<bool, S_D + 32> re_stream_set_end
        ("Readmin edge - set nodes end flag");
    hls_thread_local hls::stream<tuplebuild_tuple_t, S_D> re_stream_tuple
        ("Readmin edge - tuples");
    
    /* Homomorphism data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> h_stream_sol
        ("Homomorphism - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> h_stream_sol_end
        ("Homomorphism - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> h_stream_set
        ("Homomorphism - set nodes");
    hls_thread_local hls::stream<bool, S_D> h_stream_set_end
        ("Homomorphism - set nodes end flag");
    hls_thread_local hls::stream<tuplebuild_tuple_t, S_D> h_stream_tuple
        ("Homomorphism - tuples");
    
    /* Batchbuild data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> b_stream_sol
        ("Batchbuild - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> b_stream_sol_end
        ("Batchbuild - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> b_stream_set
        ("Batchbuild - set nodes");
    hls_thread_local hls::stream<bool, S_D> b_stream_set_end
        ("Batchbuild - set nodes end flag");
    hls_thread_local hls::stream<tuplebuild_tuple_t, S_D> b_stream_tuple
        ("Batchbuild - tuples");
    
    /* Tuplebuild data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> t_stream_sol
        ("Tuplebuild - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> t_stream_sol_end
        ("Tuplebuild - partial solution end flag");
    hls_thread_local hls::stream<intersect_tuple_t, S_D + 32> t_stream_tuple
        ("Tuplebuild - tuples");
    hls_thread_local hls::stream<bool, S_D + 32> t_stream_tuple_end
        ("Tuplebuild - tuples end flag");
    
    /* Intersect data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> i_stream_sol
        ("Intersect - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> i_stream_sol_end
        ("Intersect - partial solution end flag");
    hls_thread_local hls::stream<offset_tuple_t, S_D> i_stream_tuple
        ("Intersect - tuples");
    hls_thread_local hls::stream<bool, S_D> i_stream_tuple_end
        ("Intersect - tuples end flag");
    
    /* Offset data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> o_stream_sol
        ("Offset - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> o_stream_sol_end
        ("Offset - partial solution end flag");
    hls_thread_local hls::stream<split_tuple_t, S_D> o_stream_tuple
        ("Offset - tuples");
    hls_thread_local hls::stream<bool, S_D> o_stream_tuple_end
        ("Offset - tuples end flag");
    
    /* Split data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> s_stream_sol
        ("Split - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> s_stream_sol_end
        ("Split - partial solution end flag");
    hls_thread_local hls::stream<verify_tuple_t, S_D + 32> s_stream_tuple
        ("Split - tuples");
    hls_thread_local hls::stream<bool, S_D + 32> s_stream_tuple_end
        ("Split - tuples end flag");
    
    /* Verify data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> v_stream_sol
        ("Verify - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> v_stream_sol_end
        ("Verify - partial solution end flag");
    hls_thread_local hls::stream<compact_tuple_t, S_D> v_stream_tuple
        ("Verify - tuples");
    hls_thread_local hls::stream<bool, S_D> v_stream_tuple_end
        ("Verify - tuples end flag");
    
    /* Compact data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> c_stream_sol
        ("Compact - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> c_stream_sol_end
        ("Compact - partial solution end flag");
    hls_thread_local hls::stream<assembly_tuple_t, S_D> c_stream_tuple
        ("Compact - tuples");
    hls_thread_local hls::stream<bool, S_D> c_stream_tuple_end
        ("Compact - tuples end flag");

    /* Filter data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> f_stream_sol
        ("Filter - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> f_stream_sol_end
        ("Filter - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> f_stream_set
        ("Filter - set nodes");
    hls_thread_local hls::stream<bool, S_D> f_stream_set_end
        ("Filter - set nodes end flag");
    
    /* Assembly data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, DYN_FIFO_DEPTH> a_stream_sol
        ("Assembly - partial solution");

    /* Dynamic fifo data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, DYN_FIFO_DEPTH> dyn_stream_sol
        ("Dynamic fifo - partial solution");
    hls_thread_local hls::stream<bool, 4> dyn_stream_ovf
        ("Dynamic fifo - overflow");
 
    /* Stop signals */
    hls_thread_local hls::stream<bool, 4> streams_stop[STOP_S];
    hls_thread_local hls::stream<bool, S_D> merge_out;
    hls_thread_local hls::stream<bool, S_D> merge_in[MERGE_IN_STREAMS];
    hls_thread_local hls::stream<bool, 4> stream_sol0;
#pragma HLS array_partition variable=merge_in type=complete    

#ifndef __SYNTHESIS__
    for (int g = 0; g < STOP_S; g++) {
        char stream_name[10];
        sprintf(stream_name, "stop_%d", g);
        streams_stop[g].set_name(stream_name);
    } 
#endif /* __SYNTHESIS__ */

    htb_cache_t htb_cache(htb_buf0);
    // bloom_cache_t bloom_cache(bloom_p);

    dynfifo_init<
        ap_uint<V_ID_W>,    /* fifo data type */
        row_t,              /* fifo data type */
        DYN_FIFO_DEPTH,     /* in/out stream size */
        DYN_FIFO_BURST * 2, /* load/store stream size */
        DDR_WORD,           /* bitwidth ddr word */
        DYN_FIFO_BURST,     /* burst transaction size */
        RESULTS_SPACE>      /* memory words available */
        (res_buf,
         dynfifo_diagnostic,
         dynfifo_space,
         a_stream_sol,
         dyn_stream_sol,
         streams_stop[STOP_S - 2],
         streams_stop[STOP_S - 1],
         dyn_stream_ovf);

    hls_thread_local hls::task mwj_propose_t(
        mwj_propose,
        dyn_stream_sol,
        p0_stream_sol,
        p0_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_findmin(
        mwj_bypass_sol,
        e_stream_sol,
        e_stream_sol_end,
        p_stream_sol,
        p_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_readmin_counter(
        mwj_bypass_sol,
        p_stream_sol,
        p_stream_sol_end,
        rc_stream_sol,
        rc_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_readmin_edge(
        mwj_bypass_sol,
        rc_stream_sol,
        rc_stream_sol_end,
        re_stream_sol,
        re_stream_sol_end);

    hls_thread_local hls::task mwj_homomorphism_t(
        mwj_homomorphism,
        re_stream_set,
        re_stream_set_end,
        re_stream_tuple,
        re_stream_sol,
        re_stream_sol_end,
        h_stream_set,
        h_stream_set_end,
        h_stream_tuple,
        h_stream_sol,
        h_stream_sol_end);
    
    hls_thread_local hls::task mwj_batchbuild_t(
        mwj_batchbuild<(1UL << PROPOSE_BATCH_LOG)>,
        h_stream_set,
        h_stream_set_end,
        h_stream_tuple,
        h_stream_sol,
        h_stream_sol_end,
        merge_in[0],
        b_stream_set,
        b_stream_set_end,
        b_stream_tuple,
        b_stream_sol,
        b_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_intersect(
        mwj_bypass_sol,
        t_stream_sol,
        t_stream_sol_end,
        i_stream_sol,
        i_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_offset(
        mwj_bypass_sol,
        i_stream_sol,
        i_stream_sol_end,
        o_stream_sol,
        o_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_split(
        mwj_bypass_sol,
        o_stream_sol,
        o_stream_sol_end,
        s_stream_sol,
        s_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_verify(
        mwj_bypass_sol,
        s_stream_sol,
        s_stream_sol_end,
        v_stream_sol,
        v_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_compact(
        mwj_bypass_sol,
        v_stream_sol,
        v_stream_sol_end,
        c_stream_sol,
        c_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_filter(
        mwj_bypass_sol,
        c_stream_sol,
        c_stream_sol_end,
        f_stream_sol,
        f_stream_sol_end);

    hls_thread_local hls::task mwj_offset_t(
        mwj_offset,
        i_stream_tuple,
        i_stream_tuple_end,
        o_stream_tuple,
        o_stream_tuple_end);

    hls_thread_local hls::task mwj_split_t(
        mwj_split<(1UL << PROPOSE_BATCH_LOG)>,
        o_stream_tuple,
        o_stream_tuple_end,
        s_stream_tuple,
        s_stream_tuple_end);

    hls_thread_local hls::task mwj_compact_t(
        mwj_compact,
        v_stream_tuple,
        v_stream_tuple_end,
        c_stream_tuple,
        c_stream_tuple_end);

    hls_thread_local hls::task mwj_filter_t(
        mwj_filter<(1UL << PROPOSE_BATCH_LOG)>,
        c_stream_tuple,
        c_stream_tuple_end,
        f_stream_set,
        f_stream_set_end);

#ifdef __SYNTHESIS__

    // cache_wrapper(mwj_findmin<
    //                   T_BLOOM,
    //                   BLOOM_LOG,
    //                   K_FUN_LOG>,
    //               hash1_w,
    //               qVertices0,
    //               bloom_cache,
    //               dyn_stream_sol,
    //               streams_stop[0],
    //               p_stream_tuple,
    //               p_stream_filter,
    //               p_stream_sol,
    //               p_stream_sol_end);
    
    // stack<
    //     ap_uint<V_ID_W>,
    //     V_ID_W,
    //     MAX_QV,
    //     DYN_FIFO_DEPTH,
    //     8192>(
    //     a_stream_sol,
    //     dyn_stream_sol,
    //     streams_stop[STOP_S - 1],
    //     dyn_stream_ovf);

    mwj_edgebuild<LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>(hash1_w,
                                       qVertices0,
                                       p0_stream_sol,
                                       p0_stream_sol_end,
                                       streams_stop[0],
                                       e_stream_tuple,
                                       e_stream_sol,
                                       e_stream_sol_end);

    mwj_findmin<T_BLOOM, BLOOM_LOG, K_FUN_LOG>(bloom_p,
                                               e_stream_tuple,
                                               streams_stop[1],
                                               p_stream_tuple,
                                               p_stream_filter);

    mwj_readmin_counter<T_BLOOM, BLOOM_LOG, K_FUN_LOG, LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>(
      hash1_w,
      hash2_w,
      hTables0,
      htb_buf1,
      p_stream_tuple,
      p_stream_filter,
      streams_stop[2],
      rc_stream_tuple,
      rc_stream_filter);

    mwj_readmin_edge<T_BLOOM, BLOOM_LOG, K_FUN_LOG, FULL_HASH_W>(
      htb_buf2,
      rc_stream_tuple,
      rc_stream_filter,
      streams_stop[3],
      re_stream_set,
      re_stream_set_end,
      re_stream_tuple);

    mwj_tuplebuild<PROPOSE_BATCH_LOG, LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>(
        hash1_w,
        hash2_w,
        qVertices0,
        b_stream_set,
        b_stream_set_end,
        b_stream_tuple,
        b_stream_sol,
        b_stream_sol_end,
        streams_stop[4],
        t_stream_tuple,
        t_stream_tuple_end,
        t_stream_sol,
        t_stream_sol_end);

    intersectcache_wrapper<PROPOSE_BATCH_LOG>(
        hTables1,
        htb_cache,
        t_stream_tuple,
        t_stream_tuple_end,
        streams_stop[5],
        i_stream_tuple,
        i_stream_tuple_end);

    verifycache_wrapper<PROPOSE_BATCH_LOG>(
        hTables1,
        htb_cache,
        s_stream_tuple,
        s_stream_tuple_end,
        streams_stop[6],
        v_stream_tuple,
        v_stream_tuple_end);

    mwj_assembly(
        nQueryVer,
        f_stream_set,
        f_stream_set_end,
        f_stream_sol,
        f_stream_sol_end,
        stream_batch,
        stream_batch_end,
        streams_stop,
        stream_sol0,
        merge_in[1],
        a_stream_sol,
        result);

#else
  
    for (int g = 0; g < STOP_S - 2; g++) 
        hls::stream_globals::incr_task_counter();
    
    htb_cache.init();
    // bloom_cache.init();

    // std::thread stack_t(
    //     stack<ap_uint<V_ID_W>,
    //           V_ID_W,
    //           MAX_QV,
    //           DYN_FIFO_DEPTH,
    //           8192>,
    //     std::ref(a_stream_sol),
    //     std::ref(dyn_stream_sol),
    //     std::ref(streams_stop[STOP_S - 1]),
    //     std::ref(dyn_stream_ovf));

    std::thread mwj_edgebuild_t(mwj_edgebuild<LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>,
                                hash1_w,
                                qVertices0,
                                std::ref(p0_stream_sol),
                                std::ref(p0_stream_sol_end),
                                std::ref(streams_stop[0]),
                                std::ref(e_stream_tuple),
                                std::ref(e_stream_sol),
                                std::ref(e_stream_sol_end));

    std::thread mwj_findmin_t(mwj_findmin<T_BLOOM, BLOOM_LOG, K_FUN_LOG>,
                              bloom_p,
                              // std::ref(bloom_cache),
                              std::ref(e_stream_tuple),
                              std::ref(streams_stop[1]),
                              std::ref(p_stream_tuple),
                              std::ref(p_stream_filter));

    std::thread mwj_readmin_counter_t(mwj_readmin_counter<T_BLOOM,
                                                          BLOOM_LOG,
                                                          K_FUN_LOG,
                                                          LKP3_HASH_W,
                                                          MAX_HASH_W,
                                                          FULL_HASH_W>,
                                      hash1_w,
                                      hash2_w,
                                      hTables0,
                                      htb_buf1,
                                      std::ref(p_stream_tuple),
                                      std::ref(p_stream_filter),
                                      std::ref(streams_stop[2]),
                                      std::ref(rc_stream_tuple),
                                      std::ref(rc_stream_filter));

    std::thread mwj_readmin_edge_t(
      mwj_readmin_edge<T_BLOOM, BLOOM_LOG, K_FUN_LOG, FULL_HASH_W>,
      htb_buf2,
      std::ref(rc_stream_tuple),
      std::ref(rc_stream_filter),
      std::ref(streams_stop[3]),
      std::ref(re_stream_set),
      std::ref(re_stream_set_end),
      std::ref(re_stream_tuple));

    std::thread mwj_tuplebuild_t(
      mwj_tuplebuild<PROPOSE_BATCH_LOG, LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>,
      hash1_w,
      hash2_w,
      qVertices0,
      std::ref(b_stream_set),
      std::ref(b_stream_set_end),
      std::ref(b_stream_tuple),
      std::ref(b_stream_sol),
      std::ref(b_stream_sol_end),
      std::ref(streams_stop[4]),
      std::ref(t_stream_tuple),
      std::ref(t_stream_tuple_end),
      std::ref(t_stream_sol),
      std::ref(t_stream_sol_end));

    std::thread mwj_intersect_t(
        mwj_intersect<PROPOSE_BATCH_LOG>,
        hTables1,
        std::ref(htb_cache),
        std::ref(t_stream_tuple),
        std::ref(t_stream_tuple_end),
        std::ref(streams_stop[5]),
        std::ref(i_stream_tuple),
        std::ref(i_stream_tuple_end));

    std::thread mwj_verify_t(
        mwj_verify<PROPOSE_BATCH_LOG>,
        hTables1,
        std::ref(htb_cache),
        std::ref(s_stream_tuple),
        std::ref(s_stream_tuple_end),
        std::ref(streams_stop[6]),
        std::ref(v_stream_tuple),
        std::ref(v_stream_tuple_end));

    std::thread mwj_assembly_t(
        mwj_assembly,
        nQueryVer,
        std::ref(f_stream_set),
        std::ref(f_stream_set_end),
        std::ref(f_stream_sol),
        std::ref(f_stream_sol_end),
        std::ref(stream_batch),
        std::ref(stream_batch_end),
        std::ref(streams_stop),
        std::ref(stream_sol0),
        std::ref(merge_in[0]),
        std::ref(a_stream_sol),
        std::ref(result));

#endif /* __SYNTHESIS__ */

    hls_thread_local hls::task mwj_merge_t(
        mwj_merge<bool>,
        merge_in,
        merge_out);

    hls_thread_local hls::task mwj_stop_t(
        mwj_stop,
        merge_out,
        dyn_stream_ovf,
        stream_sol0);

#ifndef __SYNTHESIS__

    mwj_edgebuild_t.join();
    mwj_findmin_t.join();
    mwj_readmin_counter_t.join();
    mwj_readmin_edge_t.join();
    mwj_tuplebuild_t.join();
    mwj_intersect_t.join(); 
    mwj_verify_t.join();
    // stack_t.join();
    mwj_assembly_t.join();


#if DEBUG_STATS
    // debug::cache_hit_prop   = bloom_cache.get_n_hits(0);
    debug::cache_hit_inter  = htb_cache.get_n_l1_hits(1);
    debug::cache_hit_verify = htb_cache.get_n_l1_hits(0);
    // debug::cache_req_prop   = bloom_cache.get_n_reqs(0);
    debug::cache_req_inter  = htb_cache.get_n_l1_reqs(1);
    debug::cache_req_verify = htb_cache.get_n_l1_reqs(0);
#endif /* DEBUG_STATS */

    htb_cache.stop();
    // bloom_cache.stop();

#endif /* __SYNTHESIS__ */
}

template <typename T_BLOOM,
        size_t BLOOM_LOG, 
        size_t K_FUN_LOG,
        size_t LKP3_HASH_W,
        size_t MAX_HASH_W,
        size_t FULL_HASH_W>
void multiwayJoinWrap(
        ap_uint<DDR_W> *htb_buf0,
        ap_uint<DDR_W> *htb_buf1,
        ap_uint<DDR_W> *htb_buf2,
        ap_uint<DDR_W> *htb_buf3,
        T_BLOOM *bloom_p,
        row_t *res_buf,
        AdjHT *hTables0,
        AdjHT *hTables1,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        const unsigned short nQueryVer,
        const unsigned char hash1_w,
        const unsigned char hash2_w,
        const unsigned long dynfifo_space,
        unsigned long &dynfifo_diagnostic,

#if COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif
        )
{
#pragma HLS STABLE variable=htb_buf0
#pragma HLS STABLE variable=htb_buf1
#pragma HLS STABLE variable=htb_buf2
#pragma HLS STABLE variable=htb_buf3
#pragma HLS STABLE variable=hTables0
#pragma HLS STABLE variable=hTables1
#pragma HLS STABLE variable=qVertices0
#pragma HLS STABLE variable=qVertices1
#pragma HLS STABLE variable=nQueryVer
#pragma HLS STABLE variable=hash1_w
#pragma HLS STABLE variable=hash2_w
#pragma HLS STABLE variable=dynfifo_space
#pragma HLS dataflow

    hls::stream<bool, S_D> stream_batch_end("Stream batch end");
    hls::stream<ap_uint<V_ID_W>, S_D> stream_batch("Stream batch");

    mwj_batch<LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>(
      hash1_w, hTables0, qVertices0, htb_buf0, stream_batch_end, stream_batch);

    multiwayJoin<T_BLOOM, BLOOM_LOG, K_FUN_LOG, LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>(
      htb_buf1,
      htb_buf2,
      htb_buf3,
      bloom_p,
      res_buf,
      hTables0,
      hTables1,
      qVertices1,
      nQueryVer,
      hash1_w,
      hash2_w,
      dynfifo_space,
      stream_batch,
      stream_batch_end,
      dynfifo_diagnostic,
      result);
}

#if SOFTWARE_PREPROC
void subgraphIsomorphism(
        row_t htb_buf0[HASHTABLES_SPACE],
        row_t htb_buf1[HASHTABLES_SPACE],
        row_t htb_buf2[HASHTABLES_SPACE],
        row_t htb_buf3[HASHTABLES_SPACE],
        row_t bloom_p[BLOOM_SPACE],
        row_t res_buf[RESULTS_SPACE],
        const unsigned short numQueryVert,
        const unsigned char hash1_w,
        const unsigned char hash2_w,
        const unsigned long dynfifo_space,
        unsigned long &dynfifo_diagnostic,
        QueryVertex qVertices0[MAX_QV], 
        QueryVertex qVertices1[MAX_QV],
        AdjHT hTables0[MAX_TB],
        AdjHT hTables1[MAX_TB],

#if DEBUG_INTERFACE
        volatile unsigned int &debif_endpreprocess,
#endif /* DEBUG_INTERFACE */

#if COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif /* COUNT_ONLY */

        )
{

#pragma HLS INTERFACE mode=m_axi port=htb_buf0 bundle=prop_batch \
    max_widen_bitwidth=128 latency=1
#pragma HLS INTERFACE mode=m_axi port=htb_buf1 bundle=cache \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
    latency=1
#pragma HLS INTERFACE mode=m_axi port=htb_buf2 bundle=readmin_c \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
    latency=1
#pragma HLS INTERFACE mode=m_axi port=htb_buf3 bundle=readmin_e \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
    latency=1
#pragma HLS INTERFACE mode=m_axi port=res_buf bundle=fifo \
    max_widen_bitwidth=128 max_read_burst_length=32 max_write_burst_length=32 \
    latency=20
#pragma HLS INTERFACE mode=m_axi port=bloom_p bundle=bloom \
    max_widen_bitwidth=128 latency=20

#pragma HLS alias ports=htb_buf0,htb_buf1,htb_buf2,htb_buf3 distance=0

#pragma HLS INTERFACE mode=s_axilite port=numQueryVert
#pragma HLS INTERFACE mode=s_axilite port=hash1_w
#pragma HLS INTERFACE mode=s_axilite port=hash2_w
#pragma HLS INTERFACE mode=s_axilite port=dynfifo_space
#pragma HLS INTERFACE mode=s_axilite port=dynfifo_diagnostic
#pragma HLS INTERFACE mode=s_axilite port=return

#if DEBUG_INTERFACE
#pragma HLS INTERFACE mode=s_axilite port=debif_endpreprocess
#endif /* DEBUG_INTERFACE */

#if COUNT_ONLY
#pragma HLS INTERFACE mode=s_axilite port=result
#else
#pragma HLS INTERFACE mode=axis port=result
#endif /* COUNT_ONLY */

    unsigned long localResult = 0;
#if DEBUG_STATS
    debug::init();
#endif /* DEBUG_STATS */

#if DEBUG_INTERFACE
    ap_wait();
    debif_endpreprocess = 1;
    ap_wait();
#endif /* DEBUG_INTERFACE */

    multiwayJoinWrap<bloom_t,
                     BLOOM_FILTER_WIDTH,
                     K_FUNCTIONS,
                     HASH_LOOKUP3_BIT,
                     MAX_HASH_TABLE_BIT,
                     64>(htb_buf0,
                         htb_buf1,
                         htb_buf2,
                         htb_buf3,
                         bloom_p,
                         res_buf,
                         hTables0,
                         hTables1,
                         qVertices0,
                         qVertices1,
                         numQueryVert,
                         hash1_w,
                         hash2_w,
                         dynfifo_space,
                         dynfifo_diagnostic,
                         localResult);

    result = localResult;

#if DEBUG_STATS
    debug::print(hash1_w, hash2_w);
#endif /* DEBUG_STATS */
    // f.close();
}

#else

void subgraphIsomorphism(
        row_t htb_buf0[HASHTABLES_SPACE],
        row_t htb_buf1[HASHTABLES_SPACE],
        row_t htb_buf2[HASHTABLES_SPACE],
        row_t htb_buf3[HASHTABLES_SPACE],
        row_t bloom_p[BLOOM_SPACE],
        row_t res_buf[RESULTS_SPACE],
        const unsigned short numQueryVert,
        const unsigned short numQueryEdges,
        const unsigned long numDataEdges,
        const unsigned char hash1_w,
        const unsigned char hash2_w,
        const unsigned long dynfifo_space,
        unsigned long &dynfifo_diagnostic,

#if DEBUG_INTERFACE
        volatile unsigned int &debif_endpreprocess,
#endif /* DEBUG_INTERFACE */

#if COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif /* COUNT_ONLY */

        )
{

#pragma HLS INTERFACE mode=m_axi port=htb_buf0 bundle=prop_batch \
    max_widen_bitwidth=128 latency=1
#pragma HLS INTERFACE mode=m_axi port=htb_buf1 bundle=cache \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
    latency=1
#pragma HLS INTERFACE mode=m_axi port=htb_buf2 bundle=readmin_c \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
    latency=1
#pragma HLS INTERFACE mode=m_axi port=htb_buf3 bundle=readmin_e \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
    latency=20
#pragma HLS INTERFACE mode=m_axi port=res_buf bundle=fifo \
    max_widen_bitwidth=128 max_read_burst_length=32 max_write_burst_length=32 \
    latency=1
#pragma HLS INTERFACE mode=m_axi port=bloom_p bundle=bloom \
    max_widen_bitwidth=128 latency=20

#pragma HLS alias ports=htb_buf0,htb_buf1,htb_buf2,htb_buf3 distance=0

#pragma HLS INTERFACE mode=s_axilite port=numQueryVert
#pragma HLS INTERFACE mode=s_axilite port=numQueryEdges
#pragma HLS INTERFACE mode=s_axilite port=numDataEdges
#pragma HLS INTERFACE mode=s_axilite port=hash1_w
#pragma HLS INTERFACE mode=s_axilite port=hash2_w
#pragma HLS INTERFACE mode=s_axilite port=dynfifo_space
#pragma HLS INTERFACE mode=s_axilite port=dynfifo_diagnostic
#pragma HLS INTERFACE mode=s_axilite port=return

#if DEBUG_INTERFACE
#pragma HLS INTERFACE mode=s_axilite port=debif_endpreprocess
#endif /* DEBUG_INTERFACE */

#if COUNT_ONLY
#pragma HLS INTERFACE mode=s_axilite port=result
#else
#pragma HLS INTERFACE mode=axis port=result
#endif /* COUNT_ONLY */

#if DEBUG_STATS
    debug::init();
#endif /* DEBUG_STATS */

    QueryVertex qVertices0[MAX_QV], qVertices1[MAX_QV];
    AdjHT hTables0[MAX_TB], hTables1[MAX_TB];
    unsigned long localResult = 0;

    preprocess<row_t,
               bloom_t,
               EDGE_WIDTH,
               COUNTER_WIDTH,
               BLOOM_FILTER_WIDTH,
               K_FUNCTIONS,
               DDR_BIT,
               VERTEX_WIDTH_BIT,
               HASH_LOOKUP3_BIT,
               MAX_HASH_TABLE_BIT,
               64,
               LABEL_WIDTH,
               DEFAULT_STREAM_DEPTH,
               HASHTABLES_SPACE,
               MAX_QUERY_VERTICES,
               MAX_TABLES>(&res_buf[dynfifo_space],
                           htb_buf0,
                           bloom_p,
                           qVertices0,
                           qVertices1,
                           hTables0,
                           hTables1,
                           numQueryVert,
                           numQueryEdges,
                           numDataEdges,
                           hash1_w,
                           hash2_w);

#if DEBUG_INTERFACE
    ap_wait();
    debif_endpreprocess = 1;
    ap_wait();
#endif /* DEBUG_INTERFACE */

    multiwayJoinWrap<bloom_t,
                     BLOOM_FILTER_WIDTH,
                     K_FUNCTIONS,
                     HASH_LOOKUP3_BIT,
                     MAX_HASH_TABLE_BIT,
                     64>(htb_buf0,
                         htb_buf1,
                         htb_buf2,
                         htb_buf3,
                         bloom_p,
                         res_buf,
                         hTables0,
                         hTables1,
                         qVertices0,
                         qVertices1,
                         numQueryVert,
                         hash1_w,
                         hash2_w,
                         dynfifo_space,
                         dynfifo_diagnostic,
                         localResult);

    result = localResult;

#if DEBUG_STATS
    debug::print(hash1_w, hash2_w);
#endif /* DEBUG_STATS */
}
#endif /* SOFTWARE_PREPROC */

#pragma GCC diagnostic pop
