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
#include "hls_np_channel.h"

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
#pragma GCC diagnostic ignored "-Wc++11-compat"
// #pragma GCC diagnostic ignored "-Wunused-variable"
// #pragma GCC diagnostic ignored "-Wunused-but-set-variable"

#define STOP_S      13    
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
#define BLOCKBUILD_NUM 2

#if DEBUG_INTERFACE
    unsigned long propose_empty = 0;
    unsigned long edgebuild_empty = 0;
    unsigned long findmin_empty = 0;
    unsigned long readmin_counter_empty = 0;
    unsigned long readmin_edge_empty = 0;
    unsigned long homomorphism_empty = 0;
    unsigned long batchbuild_empty = 0;
    unsigned long tuplebuild_empty = 0;
    unsigned long intersect_empty = 0;
    unsigned long offset_empty = 0;
    unsigned long split_empty = 0;
    unsigned long verify_empty = 0;
    unsigned long compact_empty = 0;
    unsigned long filter_empty = 0;
    unsigned long assembly_empty = 0;
#endif

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
  /* Set information (edge) */
  ap_uint<V_ID_W> indexing_v;   // Indexing vertex in the edge
  unsigned char tb_index;       // Table index
  unsigned char iv_pos;         // Query indexing vertex position.
  unsigned char num_tb_indexed; // Used by sequencebuild

  unsigned int address; // Used by findmin to read the bloom filter
  bool reset;           // Used by findmin to reset the result bloom filter

  bool last; // Delimeter flag
} findmin_tuple_t;

typedef struct
{
  /* Minimum size set information (edge), no need for the delimeter
  this time since the min_set can be only one */
  ap_uint<V_ID_W> indexing_v;   // Indexing vertex in the edge
  unsigned char tb_index;       // Table index
  unsigned char iv_pos;         // Query indexing vertex position.
  unsigned char num_tb_indexed; // Used by sequencebuild
} readmin_counter_tuple_t;

typedef struct
{
  /* Minimum size set information (edge), no need for the delimeter
  this time since the min_set can be only one */
  ap_uint<V_ID_W> indexing_v;   // Indexing vertex in the edge
  unsigned char tb_index;       // Table index
  unsigned char iv_pos;         // Query indexing vertex position.
  unsigned char num_tb_indexed; // Used by sequencebuild
  
  /* Minimum set starting and ending row addresses */
  unsigned int rowstart;
  unsigned int rowend;
} readmin_edge_tuple_t;

template<typename NODE_T>
struct homomorphism_set_t
{
    NODE_T node;
    bool last;
};

template<typename NODE_T>
struct sequencebuild_set_t
{
    /* Adding the flag for min set since homomorphism aggregate tuple and set in
     * a unique stream to let sequencebuild be a simple pipeline */
    NODE_T node;  // store the node id or the data of the min set
    bool last;    // last of the set flag
    bool sol;     // solution node flag
    bool min_set; // min set information flag
};

template<typename NODE_T>
struct sequencebuild_tuple_t
{
    NODE_T node;
    unsigned char tb_index;
    unsigned char iv_pos;
    unsigned char query_edge;
    unsigned char pos;
    bool bit_last;
    bool sol;

    /* We used two delimeter to respectively signal the end of a batch or the end
    of the entire set of extension */
    bool last_set;
    bool last_batch;
};

typedef struct
{
    unsigned char tb_index;
    unsigned char iv_pos; // Query indexing vertex position.
    unsigned char num_tb_indexed;
} minset_tuple_t;

typedef struct
{
    ap_uint<V_ID_W> indexing_v;
    ap_uint<V_ID_W> indexed_v;
    ap_uint<64> addr_counter;
    unsigned char tb_index;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    bool bit_last_edge;
    bool skip_counter;
    bool last_set;
    bool last_batch;
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
    bool last_set;
    bool last_batch;
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
    bool last_set;
    bool last_batch;
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
    bool last_set;
    bool last_batch;
    edge_flag_type flag;
} verify_tuple_t;

typedef struct {
    ap_uint<V_ID_W> indexed_v;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    bool bit_last_address;
    bool bit_last_edge;
    bool bit_equal;
    bool last_set;
    bool last_batch;
} compact_tuple_t;

typedef struct
{
    ap_uint<V_ID_W> indexed_v;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    bool bit_last_edge;
    bool bit_checked;
    bool last_set;
    bool last_batch;
} assembly_tuple_t;

typedef struct
{
    ap_uint<V_ID_W> node;
    bool last_set;
} assembly_set_t;

/******** End tuple definition ********/

void
mwj_propose(const unsigned short nQueryVer,
            hls::stream<ap_uint<V_ID_W> >& stream_fifo_in,
            hls::stream<bool>& stream_stop,

            hls::stream<ap_uint<V_ID_W> >& stream_sol_out,
            hls::stream<bool>& stream_sol_end_out)
{
    const ap_uint<V_ID_W> BIT_FINAL_SOLUTION = (1UL << (V_ID_W - 1));
    const ap_uint<V_ID_W> MASK_NEW_SOLUTION = ~(1UL << (V_ID_W - 1));
    // const ap_uint<V_ID_W> MASK_END_EXTENSION = ~(1UL << (V_ID_W - 2));
    static ap_uint<V_ID_W> curSol_fifo[MAX_QV];
    static ap_uint<V_ID_W> curQV = 0;
    ap_uint<V_ID_W> readv;
    bool stop = false;

    // Read BFS solutions
PROPOSE_TASK_LOOP:
    while (true) {
      if (stream_fifo_in.read_nb(readv)) {
        if (readv.test(V_ID_W - 1)) {
          /* A node with the 31st bit asserted is a vertex of
          a radix of a solution, it will probably not change in
          future iteration */
          curQV = 0;
        PROPOSE_READ_NEW_SOLUTION_LOOP:
          while (readv.test(V_ID_W - 1)) {
#pragma HLS pipeline II = 1
            curSol_fifo[curQV++] = readv & MASK_NEW_SOLUTION;
            readv = stream_fifo_in.read();
          }

          /* In case of starting solution, i.e. solution with only
          one node, a second node with 30th bit asserted is used as
          false extension. This is done to handle special case */
          if (!readv.test(V_ID_W - 2)) {
            curSol_fifo[curQV] = readv;
          } else {
            curQV--;
          }
        } else {
          // Updating only the last vertex
          curSol_fifo[curQV] = readv;
        }

      PROPOSE_NEW_SOLUTION_OUT_LOOP:
        for (int g = 0; g <= curQV; g++) {
#pragma HLS pipeline II = 1
          /* Mark as possibile final solution the ones with only one
          node missing, in a way that the mwj_assembly can recognize.
          We are reusing the same 31st bit. */
          stream_sol_out.write((curQV == nQueryVer - 2)
                                 ? curSol_fifo[g] | BIT_FINAL_SOLUTION
                                 : curSol_fifo[g]);
          stream_sol_end_out.write(false);
        //   std::cout << curSol_fifo[g] << " ";
        }
        // std::cout << std::endl;
        stream_sol_end_out.write(true);
      }
#if DEBUG_INTERFACE
      else {
        propose_empty++;
      }
#endif

      if (stream_stop.read_nb(stop))
        break;
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
    const ap_uint<V_ID_W> MASK_FINAL_SOLUTION = ~(1UL << (V_ID_W - 1));
    ap_uint<V_ID_W> curEmb[MAX_QV];
    unsigned char curQV = 0;
    bool stop, last;
    findmin_tuple_t tuple_out;

    //Initializing filter in findmin
    tuple_out.reset = true;
    tuple_out.last = false;
    tuple_out.address = 0;
    stream_tuple_out.write(tuple_out);

EDGEBUILD_TASK_LOOP:
    while (true) {
        if (stream_sol_end_in.read_nb(last)) {
            curQV = 0;

        EDGEBUILD_COPYING_EMBEDDING_LOOP:
            while (!last) {
                curEmb[curQV] = stream_sol_in.read();
                stream_sol_out.write(curEmb[curQV]);
                stream_sol_end_out.write(false);
                curQV++;
                last = stream_sol_end_in.read();
            }
            stream_sol_end_out.write(true);

        EDGEBUILD_MAIN_LOOP:
            for (int g = 0; g < qVertices[curQV].numTablesIndexed; g++) {
                unsigned char tb_index = qVertices[curQV].tables_indexed[g];
                unsigned char iv_pos = qVertices[curQV].vertex_indexing[g];

                // Computing addresses of indexed sets
                ap_uint<LKP3_HASH_W> hash_out;
                ap_uint<MAX_HASH_W> hash_trimmed;
                xf::database::details::hashlookup3_core<V_ID_W>(
                  curEmb[iv_pos] & MASK_FINAL_SOLUTION, hash_out);
                hash_trimmed = hash_out;
                hash_trimmed = hash_trimmed.range(hash1_w - 1, 0);
                unsigned int address =
                  (tb_index * (1UL << hash1_w)) + hash_trimmed;
                tuple_out.indexing_v = curEmb[iv_pos] & MASK_FINAL_SOLUTION;
                tuple_out.iv_pos = iv_pos;
                tuple_out.tb_index = tb_index;
                tuple_out.address = address;
                tuple_out.reset = false;
                tuple_out.num_tb_indexed = qVertices[curQV].numTablesIndexed;
                tuple_out.last = (g == (qVertices[curQV].numTablesIndexed - 1));
                stream_tuple_out.write(tuple_out);
            }
            tuple_out.reset = true;
            tuple_out.last = false;
            stream_tuple_out.write(tuple_out);
        }
#if DEBUG_INTERFACE
        else {
            edgebuild_empty++;
        }
#endif

        if (stream_stop.read_nb(stop))
            break;
    }
}

template<typename T_BLOOM, size_t BLOOM_LOG, size_t K_FUN_LOG>
unsigned int
bloom_bitset(T_BLOOM filter)
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
                  tuple_out.num_tb_indexed = tuple_in.num_tb_indexed;
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
#if DEBUG_INTERFACE
        else {
            findmin_empty++;
        }
#endif

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
#pragma HLS pipeline II = 8
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
            tuple_out.num_tb_indexed = tuple_in.num_tb_indexed;
            tuple_out.rowstart = rowstart;
            tuple_out.rowend = rowend;

            stream_tuple_out.write(tuple_out);

#if DEBUG_STATS
            debug::readmin_counter_reads += 2;
#endif
        }
#if DEBUG_INTERFACE
        else {
            readmin_counter_empty++;
        }
#endif

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

                 hls::stream<homomorphism_set_t<ap_uint<V_ID_W> > >& stream_set_out,
                 hls::stream<minset_tuple_t>& stream_tuple_out)
{
    constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    T_BLOOM filter[K_FUN];
#pragma HLS array_partition variable = filter type = complete dim = 1
    readmin_edge_tuple_t tuple_in;
    hls::stream<ap_uint<V_ID_W>, 2> hash_in_s;
    hls::stream<ap_uint<FULL_HASH_W>, 2> hash_out_s;
    minset_tuple_t tuple_out;
    homomorphism_set_t< ap_uint<V_ID_W> > set_out;
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
            tuple_out.num_tb_indexed = tuple_in.num_tb_indexed;
            stream_tuple_out.write(tuple_out);

            for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
                filter[g] = stream_filter_in.read();
            }
        
        unsigned int cycles = tuple_in.rowend - tuple_in.rowstart;
        READMIN_EDGES_MAIN_LOOP:
            for (int g = 0; g <= cycles; g++) {
#pragma HLS pipeline II = 2
#pragma HLS allocation function instances=hash_wrapper<V_ID_W> limit=1 
#pragma HLS allocation function instances=bloom_test<T_BLOOM, BLOOM_LOG, K_FUN, FULL_HASH_W> limit=1 
                row_t row = m_axi[tuple_in.rowstart + g];
                for (int i = 0; i < EDGE_ROW; i++) {
#pragma HLS unroll
                    edge = row.range(((i + 1) << E_W) - 1, i << E_W);
                    vertexCheck = edge.range(V_ID_W * 2 - 1, V_ID_W);
                    vertex = edge.range(V_ID_W - 1, 0);

                    hash_wrapper<V_ID_W>(vertex, hash_out);
                    bool test = true;
                    bloom_test<T_BLOOM, BLOOM_LOG, K_FUN, FULL_HASH_W>(
                      filter, hash_out, test);

                    set_out.node = vertex;
                    set_out.last = false;
                    if (test && (tuple_in.indexing_v == vertexCheck)) {
                      stream_set_out.write(set_out);
                    }
#if DEBUG_STATS
                    if (tuple_in.indexing_v == vertexCheck) {
                      if (test) {
                        debug::readmin_vstream++;
                      } else {
                        debug::bloom_filter++;
                      }
                    }
#endif
                }
            }
            set_out.last = true;
            stream_set_out.write(set_out);

#if DEBUG_STATS
            debug::readmin_edge_reads +=
              ceil(((tuple_in.rowend - tuple_in.rowstart) / 16.0));
            debug::readmin_n_sets++;
#endif
        }
#if DEBUG_INTERFACE
        else {
            readmin_edge_empty++;
        }
#endif

        if (stream_stop.read_nb(stop))
            break;
    }
}

void
mwj_homomorphism(
  hls::stream<homomorphism_set_t<ap_uint<V_ID_W> > >& stream_set_in,
  hls::stream<minset_tuple_t>& stream_tuple_in, 
  hls::stream<ap_uint<V_ID_W>>& stream_sol_in,
  hls::stream<bool>& stream_sol_end_in,
  hls::stream<bool>& stream_stop,

  hls::stream<sequencebuild_set_t<ap_uint<V_ID_W> > >& stream_set_out)
{
    const ap_uint<V_ID_W> MASK_FINAL_SOLUTION = ~(1UL << (V_ID_W - 1));
    ap_uint<8> curQV;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<V_ID_W> fake_node;
    ap_uint<MAX_QV> equal_bits;
    ap_uint<MAX_QV> valid_bits;
    homomorphism_set_t<ap_uint<V_ID_W> > set_in;
    sequencebuild_set_t<ap_uint<V_ID_W> > set_out;
    bool last_sol, stop;

    while (true) {
        if (stream_sol_end_in.read_nb(last_sol)) {
            curQV = 0;

        HOMOMORPHISM_COPYING_EMBEDDING_LOOP:
            do {
                ap_uint<V_ID_W> node = stream_sol_in.read();
                set_out.node = node;
                set_out.last = false;
                set_out.min_set = false;
                set_out.sol = true;
                stream_set_out.write(set_out);
                curEmb[curQV] = node & MASK_FINAL_SOLUTION;
                curQV++;
                last_sol = stream_sol_end_in.read();
                // std::cout << (int)curEmb[curQV - 1] << " sol in " << select << std::endl;
            }while (!last_sol);
            
            /* Fake node is the tuple about min_set data */
            minset_tuple_t tuple = stream_tuple_in.read();
            fake_node.range(7, 0) = tuple.tb_index;
            fake_node.range(15, 8) = tuple.iv_pos;
            fake_node.range(31, 16) = tuple.num_tb_indexed;
            set_out.last = false;
            set_out.min_set = true;
            set_out.sol = false;
            set_out.node = fake_node;
            stream_set_out.write(set_out);
            valid_bits = (1UL << curQV) - 1;
        
        HOMOMORPHISM_CHECK_LOOP:
            do {
#pragma HLS pipeline II = 1
                set_in = stream_set_in.read();
                ap_uint<V_ID_W> vToVerify = set_in.node;
                equal_bits = 0;

                for (int g = 0; g < MAX_QV; g++) {
#pragma HLS unroll
                    if (vToVerify == curEmb[g]){
                        equal_bits[g] = 1;
                    }
                }

                set_out.node = vToVerify;
                set_out.last = set_in.last;
                set_out.min_set = false;
                set_out.sol = false;

                /* Write out in case of not duplicate in current solution or
                delimeter node */
                if ((equal_bits & valid_bits) == 0 || set_in.last) {
                    stream_set_out.write(set_out);
                }
#if DEBUG_STATS
                else {
                    debug::homomo_trashed++;
                }
#endif
            } while (!set_in.last); 
        }
#if DEBUG_INTERFACE
        else {
            homomorphism_empty++;
        }
#endif

        if (stream_stop.read_nb(stop))
            break;
    }
}

template<size_t BATCH_SIZE_LOG>
void
mwj_sequencebuild(
  hls::stream<sequencebuild_set_t<ap_uint<V_ID_W>>>& stream_tuple_in,
  hls::stream<sequencebuild_tuple_t<ap_uint<V_ID_W>>>& stream_tuple_out)
{
#pragma HLS pipeline II = 1
    static ap_uint<V_ID_W> buffer[(1UL << BATCH_SIZE_LOG)];
    static ap_uint<BATCH_SIZE_LOG + 1> buffer_size = 0;
    static ap_uint<BATCH_SIZE_LOG + 1> buffer_p = 0;
    static unsigned short cycles = 0;
    static unsigned char query_edge = 0;
    static unsigned char tb_index = 0;
    static unsigned char iv_pos = 0;
    ap_uint<V_ID_W> vToVerify;
    sequencebuild_set_t< ap_uint<V_ID_W> > tuple_in;
    sequencebuild_tuple_t< ap_uint<V_ID_W> > tuple_out;
    bool last_inner;
    
    if (query_edge == 0) {
        tuple_in = stream_tuple_in.read();

        /* Two different cases of termination, one is given by the last
        delimeter and the other is given by the overflow of the batch size*/
        last_inner =
          tuple_in.last || (buffer_p == ((1UL << BATCH_SIZE_LOG) - 1));

        vToVerify = tuple_in.node;
        buffer[buffer_p] = vToVerify;

        if (tuple_in.min_set) {

            /* Passing information of min_set through the node variable to
            save some bits. */
            tb_index = vToVerify.range(7, 0);
            iv_pos = vToVerify.range(15, 8);
            cycles = vToVerify.range(31, 16);
            tuple_out.sol = true;
            tuple_out.last_set = true;
            stream_tuple_out.write(tuple_out);

        } else if (tuple_in.sol) {
            tuple_out.last_set = false;
            tuple_out.sol = true;
            tuple_out.node = tuple_in.node;
            stream_tuple_out.write(tuple_out);

        } else {
            tuple_out.node = vToVerify;
            tuple_out.tb_index = tb_index;
            tuple_out.iv_pos = iv_pos;
            tuple_out.pos = buffer_p;
            tuple_out.query_edge = query_edge;
            tuple_out.bit_last = (cycles == 1);
            tuple_out.last_set = tuple_in.last;
            tuple_out.sol = false;
            tuple_out.last_batch =
              (buffer_p == ((1UL << BATCH_SIZE_LOG) - 1)) && (cycles == 1);

            /* Last is not a real node, so do not increment */
            buffer_p += (tuple_in.last) ? 0 : 1;
            if (!tuple_in.last || cycles == 1) {
                stream_tuple_out.write(tuple_out);
            }

            if (last_inner) {
                query_edge++;
                buffer_size = buffer_p;
                // std::cout << "last inner cycle 0 b_s: " << (int)buffer_size
                // << " b_p: " << (int)buffer_p << std::endl;
                buffer_p = 0;
                if (cycles == 1) {
                    /* Checking outer loop end condition */
                    buffer_size = 0;
                    query_edge = 0;
                    // std::cout << "last outer" << std::endl;
                }
            }
        }

    } else {
        last_inner = (buffer_p == buffer_size);
        vToVerify = buffer[buffer_p];

        if (!last_inner) {
            tuple_out.node = vToVerify;
            tuple_out.tb_index = tb_index;
            tuple_out.iv_pos = iv_pos;
            tuple_out.pos = buffer_p;
            tuple_out.query_edge = query_edge;
            tuple_out.bit_last = (query_edge == (cycles - 1));
            tuple_out.sol = false;
            tuple_out.last_set = false;
            tuple_out.last_batch =
              (buffer_p == ((1UL << BATCH_SIZE_LOG) - 1)) &&
              query_edge == (cycles - 1);
            stream_tuple_out.write(tuple_out);
            buffer_p++;
            // std::cout << "read from mem " << (int)vToVerify << " b_p: "
            // << (int)buffer_p  << " last: "
            //               << tuple_in.last << " batch: " <<
            //               tuple_out.last_batch << std::endl;
        } else {

            /* Checking inner loop end condition.
            Inner loop cycles on nodes over an edge; Doing single node over
            the edges would have been much more easier but this rearrangment
            achieve better hit ratios since each edge is a single table */
            query_edge++;
            // std::cout << "last inner cycle >0 b_s: " << (int)buffer_size
            // << " b_p: " << (int)buffer_p << std::endl;
            buffer_p = 0;
            if (query_edge == cycles) {
                /* Checking outer loop end condition */
                query_edge = 0;
                if (buffer_size != (1UL << BATCH_SIZE_LOG)) {
                    // std::cout << "write last " << (int)vToVerify <<
                    // std::endl;
                    tuple_out.last_batch = false;
                    tuple_out.last_set = true;
                    tuple_out.sol = false;
                    stream_tuple_out.write(tuple_out);
                }
                buffer_size = 0;
            }
        }
    }
}

template<size_t BATCH_SIZE_LOG,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t FULL_HASH_W>
void
mwj_tuplebuild(
  const unsigned char hash1_w,
  const unsigned char hash2_w,
  QueryVertex* qVertices,
  hls::stream<sequencebuild_tuple_t<ap_uint<V_ID_W>>>& stream_tuple_in,
  hls::stream<bool>& stream_stop,

  hls::stream<intersect_tuple_t> stream_tuple_out[2],
  hls::stream<ap_uint<V_ID_W>>& stream_sol_out,
  hls::stream<bool>& stream_sol_end_out)
{
    const ap_uint<V_ID_W> MASK_FINAL_SOLUTION = ~(1UL << (V_ID_W - 1));
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<V_ID_W> vToVerify;
    hls::stream<ap_uint<V_ID_W>, 4> hash_in0, hash_in1;
    hls::stream<ap_uint<LKP3_HASH_W>, 4> hash_out0, hash_out1;
    unsigned long addr_counter;
    intersect_tuple_t tuple_out;
    sequencebuild_tuple_t< ap_uint<V_ID_W> > tuple_in;
    unsigned char curQV = 0;
    bool stop;

TUPLEBUILD_TASK_LOOP:
    while (true) {
#pragma HLS pipeline II = 1
        if (stream_tuple_in.read_nb(tuple_in)) {

            if (tuple_in.sol) {
                curEmb[curQV] = tuple_in.node;
                stream_sol_end_out.write(tuple_in.last_set);
                if (!tuple_in.last_set) {
                    stream_sol_out.write(tuple_in.node);
                    curQV++;
                }
            } else {
                uint8_t tableIndex =
                  qVertices[curQV].tables_indexed[tuple_in.query_edge];
                uint8_t ivPos =
                  qVertices[curQV].vertex_indexing[tuple_in.query_edge];
                
                if (tuple_in.last_set) {
                    curQV = 0;
                }

                bool bit_last = tuple_in.bit_last;
                bool bit_min =
                  (tuple_in.tb_index == tableIndex && tuple_in.iv_pos == ivPos);

                vToVerify = tuple_in.node;
                hash_in0.write(vToVerify);
                hash_in1.write(curEmb[ivPos] & MASK_FINAL_SOLUTION);
                xf::database::hashLookup3<V_ID_W>(hash_in0, hash_out0);
                xf::database::hashLookup3<V_ID_W>(hash_in1, hash_out1);
                ap_uint<MAX_HASH_W> indexed_h = hash_out0.read();
                ap_uint<MAX_HASH_W> indexing_h = hash_out1.read();

                addr_counter = indexing_h.range(hash1_w - 1, 0);
                addr_counter <<= hash2_w;
                addr_counter += indexed_h.range(hash2_w - 1, 0);

                tuple_out.indexed_v = vToVerify;
                tuple_out.indexing_v = curEmb[ivPos] & MASK_FINAL_SOLUTION;
                tuple_out.addr_counter = addr_counter;
                tuple_out.tb_index = tableIndex;
                tuple_out.pos = tuple_in.pos;
                tuple_out.bit_last_edge = bit_last;
                tuple_out.flag = (bit_min) ? MIN_SET : CHECK;
                tuple_out.skip_counter = false;
                tuple_out.last_set = tuple_in.last_set;
                tuple_out.last_batch = tuple_in.last_batch;

                stream_tuple_out[0].write(tuple_out);

                if (addr_counter == 0)
                    tuple_out.skip_counter = true;
                tuple_out.addr_counter = addr_counter - 1;

                stream_tuple_out[1].write(tuple_out);
            }
        }
#if DEBUG_INTERFACE
        else {
            tuplebuild_empty++;
        }
#endif

        if (stream_stop.read_nb(stop)){
            break;
        }
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
    hls::stream<intersect_tuple_t> stream_tuple_in[2],
    hls::stream<bool> &stream_stop,

    hls::stream<offset_tuple_t> &stream_tuple_out)
{
    ap_uint<V_ID_W> indexing_v;
    intersect_tuple_t tuple_in;
    offset_tuple_t tuple_out;
    unsigned char tableIndex;
    ap_uint<DDR_W> ram_row;
    unsigned long addr_row;
    ap_uint<DDR_BIT - C_W> addr_inrow;
    ap_uint<64> addr_counter;
    unsigned char stream_p = 0;
    bool stop;

INTERSECT_TASK_LOOP:
    while (1) {
#pragma HLS pipeline II=1
        if (stream_tuple_in[stream_p].read_nb(tuple_in)){
            stream_p = (stream_p + 1) % 2;

            if (!tuple_in.last_set){
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

                tuple_out.indexed_v = tuple_in.indexed_v;
                tuple_out.indexing_v = tuple_in.indexing_v;
                tuple_out.tb_index = tuple_in.tb_index;
                tuple_out.pos = tuple_in.pos;
                tuple_out.offset = offset;
                tuple_out.bit_last_edge = tuple_in.bit_last_edge;
                tuple_out.flag = tuple_in.flag;
            }
            tuple_out.last_batch = tuple_in.last_batch;
            tuple_out.last_set = tuple_in.last_set;
            stream_tuple_out.write(tuple_out);
        }
#if DEBUG_INTERFACE
        else {
            intersect_empty++;
        }
#endif

        if (stream_stop.read_nb(stop))
            break;
    }
}

template<size_t NUM_SPLIT>
void mwj_offset(
    hls::stream<offset_tuple_t> &stream_tuple_in,

    hls::stream<split_tuple_t> stream_tuple_out[NUM_SPLIT])
{
#pragma HLS pipeline II = 2 style = flp
    constexpr size_t EDGE_BLOCK = (CACHE_WORDS_PER_LINE + DDR_BIT - E_W);
    static unsigned char pointer = 0;
    offset_tuple_t tuple_in;
    split_tuple_t tuple_out;

    tuple_in = stream_tuple_in.read();
    ap_uint<(1UL << C_W)> end_off = tuple_in.offset;
    tuple_out.indexed_v = tuple_in.indexed_v;
    tuple_out.indexing_v = tuple_in.indexing_v;
    tuple_out.pos = tuple_in.pos;
    tuple_out.tb_index = tuple_in.tb_index;
    tuple_out.bit_last_edge = tuple_in.bit_last_edge;

    tuple_in = stream_tuple_in.read();

    tuple_out.last_set = tuple_in.last_set;
    tuple_out.last_batch = tuple_in.last_batch;
    tuple_out.start_off = tuple_in.offset;
    tuple_out.first_block = tuple_in.offset >> EDGE_BLOCK;
    tuple_out.end_block = (tuple_in.offset == end_off)
                            ? (unsigned int)(end_off >> EDGE_BLOCK)
                            : (unsigned int)((end_off - 1) >> EDGE_BLOCK);
    tuple_out.flag = (tuple_in.flag == MIN_SET)
                       ? MIN_SET
                       : ((tuple_in.offset == end_off) ? NO_EDGE : CHECK);

#if DEBUG_STATS
    debug::intersect_filter += (tuple_out.flag == NO_EDGE) ? 1 : 0;
#endif /* DEBUG_STATS */
    stream_tuple_out[pointer].write(tuple_out);
    pointer = (pointer + 1) % NUM_SPLIT;
}

template<size_t MAX_BATCH_SIZE, size_t ID>
void
mwj_blockbuild(hls::stream<split_tuple_t>& stream_tuple_in,

          hls::stream<verify_tuple_t>& stream_tuple_out)
{
    constexpr size_t EDGE_BLOCK = (CACHE_WORDS_PER_LINE + DDR_BIT - E_W);
    split_tuple_t tuple_in;
    verify_tuple_t tuple_out;

    tuple_in = stream_tuple_in.read();
    if (!tuple_in.last_set) {

        tuple_out.indexing_v = tuple_in.indexing_v;
        tuple_out.indexed_v = tuple_in.indexed_v;
        tuple_out.tb_index = tuple_in.tb_index;
        tuple_out.pos = tuple_in.pos;
        tuple_out.bit_last_edge = tuple_in.bit_last_edge;
        tuple_out.flag = tuple_in.flag;
        tuple_out.last_batch = tuple_in.last_batch;
        tuple_out.last_set = false;
        unsigned int s = 0;

    SPLIT_MAIN_LOOP:
        do {
#pragma HLS pipeline II = 1
            tuple_out.address = tuple_in.start_off + s;
            tuple_out.bit_last_address =
              (tuple_in.first_block == tuple_in.end_block);
            stream_tuple_out.write(tuple_out);
            
            tuple_in.first_block++;
            s += (1UL << EDGE_BLOCK);
        } while (tuple_in.first_block <= tuple_in.end_block);
    } else {
        /* Setting bit_last_address even for delimeter nodes for
        the merge_split to change stream on which reads. In non delimeter
        nodes, bit_last_address is used to identify the last tuple generated */
        tuple_out.last_set = tuple_in.last_set;
        tuple_out.last_batch = tuple_in.last_batch;
        tuple_out.bit_last_address = true;
        stream_tuple_out.write(tuple_out);
    }
}

template<size_t NUM_SPLIT>
void
mwj_split_merge(hls::stream<verify_tuple_t> stream_tuple_in[NUM_SPLIT],

                hls::stream<verify_tuple_t>& stream_tuple_out)
{
#pragma HLS pipeline II = 1
    static unsigned char pointer = 0;
    verify_tuple_t tuple_in;

    tuple_in = stream_tuple_in[pointer].read(); 
    stream_tuple_out.write(tuple_in);

    /* We need to guarantee to read from the same stream until 
    the tuples produced from an edge verification job are concluded,
    thus not mixing tuple coming from two different edge. */
    if (tuple_in.bit_last_address){
        pointer = (pointer + 1) % NUM_SPLIT;
    }
}

template <size_t BATCH_SIZE_LOG>
void mwj_verify(
        AdjHT *hTables,
        htb_cache_t                  &htb_buf,
        hls::stream<verify_tuple_t>  &stream_tuple_in,
        hls::stream<bool>            &stream_stop,

        hls::stream<compact_tuple_t> &stream_tuple_out)
{
    constexpr size_t EDGE_PER_WORD = (DDR_BIT - E_W);
    
    verify_tuple_t tuple_in;
    compact_tuple_t tuple_out;
    ap_uint<V_ID_W> candidate_v;
    ap_uint<V_ID_W> indexing_v;
    unsigned char tableIndex;
    bool stop;
    ap_uint<128> edge_block[(1UL << CACHE_WORDS_PER_LINE)];
#pragma HLS array_partition variable=edge_block type=complete

VERIFY_TASK_LOOP:
    while(1){
#pragma HLS pipeline II=1
        if (stream_tuple_in.read_nb(tuple_in)){

            // std::cout <<
            // (int)tuple_in.indexed_v << " " <<
            // (int)tuple_in.indexing_v << " " <<
            // tuple_in.bit_last_address << " " <<
            // tuple_in.bit_last_edge << " " <<
            // tuple_in.last <<
            // std::endl;

            if(!tuple_in.last_set){
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

                tuple_out.bit_last_edge = tuple_in.bit_last_edge;
                tuple_out.bit_last_address = tuple_in.bit_last_address;
                tuple_out.indexed_v = candidate_v;
                tuple_out.pos = tuple_in.pos;
            }
            tuple_out.last_batch = tuple_in.last_batch;
            tuple_out.last_set = tuple_in.last_set;
            stream_tuple_out.write(tuple_out);
        }
#if DEBUG_INTERFACE
        else {
            verify_empty++;
        }
#endif

        if (stream_stop.read_nb(stop)){
            break;
        }
    }
}

void mwj_compact(
        hls::stream<compact_tuple_t>  &stream_tuple_in,
        hls::stream<bool>             &stream_stop,
        
        hls::stream<assembly_tuple_t> &stream_tuple_out)
{
    bool checked = false;

    compact_tuple_t tuple_in;
    assembly_tuple_t tuple_out;
    bool stop;

COMPACT_TASK_LOOP:
    while (true) {
#pragma HLS pipeline II = 1

        if (stream_tuple_in.read_nb(tuple_in)) {
            if (!tuple_in.last_set) {
                checked |= tuple_in.bit_equal;
                if (tuple_in.bit_last_address) {
                    tuple_out.bit_checked = checked;
                    tuple_out.indexed_v = tuple_in.indexed_v;
                    tuple_out.pos = tuple_in.pos;
                    tuple_out.bit_last_edge = tuple_in.bit_last_edge;
                    tuple_out.last_batch = tuple_in.last_batch; 
                    tuple_out.last_set = false; 
                    stream_tuple_out.write(tuple_out);
#if DEBUG_STATS
                    if (!checked)
                        debug::verify_filter++;
#endif /* DEBUG_STATS */
                    checked = false;
                }
            } else {
                tuple_out.last_batch = tuple_in.last_batch;
                tuple_out.last_set = tuple_in.last_set;
                stream_tuple_out.write(tuple_out);
                checked = false;
            }
        }
#if DEBUG_INTERFACE
        else {
            compact_empty++;
        }
#endif

        if (stream_stop.read_nb(stop)) {
            break;
        }
    }
}

template<size_t MAX_BATCH_SIZE>
void mwj_filter(
        hls::stream<assembly_tuple_t>   &stream_tuple_in,
        hls::stream<bool>               &stream_stop,

        hls::stream<assembly_set_t > &stream_set_out)
{
    static ap_uint<MAX_BATCH_SIZE> bits = ~0;
    assembly_tuple_t tuple_in;
    assembly_set_t tuple_out;
    bool stop;

    while (true) {
#pragma HLS pipeline II = 1

        if (stream_tuple_in.read_nb(tuple_in)) {
            if (!tuple_in.last_set) {
                unsigned short p = tuple_in.pos;
                bits[p] = bits[p] && tuple_in.bit_checked;
                tuple_out.node = tuple_in.indexed_v;
                tuple_out.last_set = false;
                if (tuple_in.bit_last_edge && bits.test(p)) {
                    stream_set_out.write(tuple_out);
                }
                if (tuple_in.last_batch && tuple_in.bit_last_edge){
                    bits = ~0;
                }
            } else {
                bits = ~0;
                tuple_out.last_set = tuple_in.last_set;
                stream_set_out.write(tuple_out);
            }
        }
#if DEBUG_INTERFACE
        else {
            filter_empty++;
        }
#endif

        if (stream_stop.read_nb(stop)) {
            break;
        }
    }
}

void mwj_assembly(
    const unsigned short nQueryVer,
    hls::stream<assembly_set_t > &stream_inter_in,
    hls::stream<ap_uint<V_ID_W> > &stream_embed_in,
    hls::stream<bool> &stream_end_embed_in,
    hls::stream<ap_uint<V_ID_W> > &stream_batch,
    hls::stream<bool> &stream_batch_end,

    hls::stream<bool> streams_stop[STOP_S],
    // hls::stream<bool> &stream_sol0,
    // hls::stream<bool> &stream_req,
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
    assembly_set_t tuple_in;
    bool last_sol, last_start;
    unsigned long partial_sol = 0;
    // bool token_new_start;
    T_NODE node;

#if COUNT_ONLY
    unsigned long int counter {0};
#endif

    last_start = stream_batch_end.read();
    last_start = stream_batch_end.read();
    // token_new_start = false;
    stream_partial_out.write(stream_batch.read() | MASK_NEW_SOLUTION);
    /* False extension for single node solutions */
    stream_partial_out.write(0 | MASK_END_EXTENSION);

    partial_sol++;
    // stream_req.write(true);
ASSEMBLY_TASK_LOOP:
    while(1) {
        if (stream_inter_in.read_nb(tuple_in)){
            curQV = 0;

            last_sol = stream_end_embed_in.read();
ASSEMBLY_COPYING_EMBEDDING_LOOP:
            while(!last_sol){
#pragma HLS pipeline II = 1
                ap_uint<V_ID_W> node = stream_embed_in.read();

                /* If there is at least an extensions for this solution, 
                and the solution is not a one with only one node missing */
                if (!tuple_in.last_set && !node.test(V_ID_W - 1)){
                    stream_partial_out.write(node | MASK_NEW_SOLUTION);
                }

                curQV++;
                last_sol = stream_end_embed_in.read();
            }

ASSEMBLY_SET_LOOP:
            while(!tuple_in.last_set){
#pragma HLS pipeline II = 1
                
                /* Write in the correct stream */
                if (curQV == nQueryVer - 1) {
#if COUNT_ONLY
                    // if (!last_start && token_new_start){
                    //     last_start = stream_batch_end.read();
                    //     token_new_start = false;
                    //     stream_partial_out.write(0 | MASK_NEW_SOLUTION);
                    //     stream_partial_out.write(stream_batch.read() |
                    //     MASK_END_EXTENSION); stream_req.write(true);
                    // }
                    // for (int g = 0; g < nQueryVer - 1; g++){
                    //     f << curEmb[g] << " ";
                    // }
                    // f << vToVerify << std::endl;
                    counter++;
#else
                ASSEMBLY_WRITE_FINAL_LOOP:
                    for (int g = 0; g < curQV; g++) {
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
                    // token_new_start = true;
                    stream_partial_out.write(tuple_in.node);
                    partial_sol++;
                    // stream_req.write(true);
                }

                /* Useless but introduced to solve a bug in Vitis HLS 2022.2 */
                ap_wait();
                ap_wait();
                tuple_in = stream_inter_in.read();
            }

            // Last batch of a set 
            partial_sol--;
            // stream_req.write(false);

        }
#if DEBUG_INTERFACE
        else {
            assembly_empty++;
        }
#endif

        if (partial_sol == 0){

#if DEBUG_STATS
                debug::miss_indexing++;
#endif

            // Test if there are some node from start batch 
            if (!last_start){
                last_start = stream_batch_end.read();
                // token_new_start = true;
                stream_partial_out.write(stream_batch.read() |
                                         MASK_NEW_SOLUTION);

                /* False extension for single node solutions */
                stream_partial_out.write(0 | MASK_END_EXTENSION);
                partial_sol++;
                // stream_req.write(true);
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
    ap_uint<V_ID_W * 2> edge;
    ap_uint<V_ID_W> vertex, vertexCheck;
    ap_uint<V_ID_W> set[MAX_CL];
    ap_uint<MAX_HASH_W> hash_buff, hash_new;
    unsigned char set_counter = 0;
    bool flag_buff = false;
    bool flag_new = true;
    hash_buff = hash_new = 0;
    // unsigned int rm_start = 0;

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
    // std::cout << rm_start << " removed\n";
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
    hls::stream<intersect_tuple_t> stream_tuple_in[2],
    hls::stream<bool> &stream_stop,

    hls::stream<offset_tuple_t> &stream_tuple_out)
{
    htb_buf.init(1);
    mwj_intersect<BATCH_SIZE_LOG>(
        hTables,
        htb_buf,
        stream_tuple_in,
        stream_stop,
        stream_tuple_out);
}

template <size_t BATCH_SIZE_LOG>
void verifycache_wrapper(
        AdjHT *hTables,
        htb_cache_t                  &htb_buf,
        hls::stream<verify_tuple_t>  &stream_tuple_in,
        hls::stream<bool>            &stream_stop,

        hls::stream<compact_tuple_t> &stream_tuple_out)
{
    htb_buf.init(0);
    mwj_verify<BATCH_SIZE_LOG>(
            hTables,
            htb_buf,
            stream_tuple_in,
            stream_stop,                  
            stream_tuple_out);       
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
    hls_thread_local hls::stream<homomorphism_set_t< ap_uint<V_ID_W> >, S_D> re_stream_set
        ("Readmin edge - set nodes");
    hls_thread_local hls::stream<minset_tuple_t, S_D> re_stream_tuple
        ("Readmin edge - tuples");

    /* Homomorphism data out */
    hls_thread_local hls::stream<sequencebuild_set_t<ap_uint<V_ID_W>>, S_D>
      h_stream_set("Homomorphism - set nodes");

    /* Sequencebuild data out */
    hls_thread_local hls::stream<sequencebuild_tuple_t<ap_uint<V_ID_W>>, S_D>
      sb_stream_set("Sequencebuild - set nodes");

    /* Tuplebuild data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> t_stream_sol
        ("Tuplebuild - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> t_stream_sol_end
        ("Tuplebuild - partial solution end flag");
    hls_thread_local hls::stream<intersect_tuple_t, S_D> t_stream_tuple[2];
    
    /* Intersect data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> i_stream_sol
        ("Intersect - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> i_stream_sol_end
        ("Intersect - partial solution end flag");
    hls_thread_local hls::stream<offset_tuple_t, S_D> i_stream_tuple
        ("Intersect - tuples");
    
    /* Offset data out */    
    hls_thread_local hls::stream<split_tuple_t, S_D> o_stream_tuple[BLOCKBUILD_NUM];
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> o_stream_sol
        ("Offset - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> o_stream_sol_end
        ("Offset - partial solution end flag");
    
    /* Blockbuild data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> bb_stream_sol
        ("Blockbuild - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> bb_stream_sol_end
        ("Blockbuild - partial solution end flag");
    hls_thread_local hls::stream<verify_tuple_t, S_D> bb_stream_tuple[BLOCKBUILD_NUM];

    /* Blockbuild merge data out */
    hls_thread_local hls::stream<verify_tuple_t, S_D> bb_merge_stream_tuple;
    
    /* Verify data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> v_stream_sol
        ("Verify - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> v_stream_sol_end
        ("Verify - partial solution end flag");
    hls_thread_local hls::stream<compact_tuple_t, S_D> v_stream_tuple
        ("Verify - tuples");
    
    /* Compact data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> c_stream_sol
        ("Compact - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> c_stream_sol_end
        ("Compact - partial solution end flag");
    hls_thread_local hls::stream<assembly_tuple_t, S_D> c_stream_tuple
        ("Compact - tuples");

    /* Filter data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> f_stream_sol
        ("Filter - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> f_stream_sol_end
        ("Filter - partial solution end flag");
    hls_thread_local hls::stream<assembly_set_t, S_D> f_stream_set
        ("Filter - set nodes");
    
    /* Assembly data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, DYN_FIFO_DEPTH> a_stream_sol
        ("Assembly - partial solution");

    /* Dynamic fifo data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, DYN_FIFO_DEPTH> dyn_stream_sol
        ("Dynamic fifo - partial solution");
    // hls_thread_local hls::stream<bool, 4> dyn_stream_ovf
    //     ("Dynamic fifo - overflow");
 
    /* Stop signals */
    hls_thread_local hls::stream<bool, 4> streams_stop[STOP_S];
    // hls_thread_local hls::stream<bool, S_D> merge_out;
    // hls_thread_local hls::stream<bool, S_D> merge_in[MERGE_IN_STREAMS];
    // hls_thread_local hls::stream<bool, 4> stream_sol0;
// #pragma HLS array_partition variable=merge_in type=complete    

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
         streams_stop[STOP_S - 1]);
        //  dyn_stream_ovf);

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

    hls_thread_local hls::task mwj_sequencebuild_t(
        mwj_sequencebuild<PROPOSE_BATCH_LOG>,
        h_stream_set,
        sb_stream_set);

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
        bb_stream_sol,
        bb_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_verify(
        mwj_bypass_sol,
        bb_stream_sol,
        bb_stream_sol_end,
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
      mwj_offset<BLOCKBUILD_NUM>, i_stream_tuple, o_stream_tuple);

    hls_thread_local hls::task mwj_blockbuild_t[BLOCKBUILD_NUM];

    for (int g = 0; g < BLOCKBUILD_NUM; g++) {
#pragma HLS unroll
        mwj_blockbuild_t[g](
          mwj_blockbuild<(1UL << PROPOSE_BATCH_LOG), __COUNTER__>,
          o_stream_tuple[g],
          bb_stream_tuple[g]);
    }

    hls_thread_local hls::task mwj_blockbuild_merge_t(
      mwj_split_merge<BLOCKBUILD_NUM>, bb_stream_tuple, bb_merge_stream_tuple);

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

    mwj_propose(nQueryVer,
                dyn_stream_sol,
                streams_stop[0],
                p0_stream_sol,
                p0_stream_sol_end);

    mwj_edgebuild<LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>(hash1_w,
                                       qVertices0,
                                       p0_stream_sol,
                                       p0_stream_sol_end,
                                       streams_stop[1],
                                       e_stream_tuple,
                                       e_stream_sol,
                                       e_stream_sol_end);

    mwj_findmin<T_BLOOM, BLOOM_LOG, K_FUN_LOG>(bloom_p,
                                               e_stream_tuple,
                                               streams_stop[2],
                                               p_stream_tuple,
                                               p_stream_filter);

    mwj_readmin_counter<T_BLOOM, BLOOM_LOG, K_FUN_LOG, LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>(
      hash1_w,
      hash2_w,
      hTables0,
      htb_buf1,
      p_stream_tuple,
      p_stream_filter,
      streams_stop[3],
      rc_stream_tuple,
      rc_stream_filter);

    mwj_readmin_edge<T_BLOOM, BLOOM_LOG, K_FUN_LOG, FULL_HASH_W>(
      htb_buf2,
      rc_stream_tuple,
      rc_stream_filter,
      streams_stop[4],
      re_stream_set,
      re_stream_tuple);

    mwj_homomorphism(re_stream_set,
                     re_stream_tuple,
                     re_stream_sol,
                     re_stream_sol_end,
                     streams_stop[5],
                     h_stream_set);

    mwj_tuplebuild<PROPOSE_BATCH_LOG, LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>(
      hash1_w,
      hash2_w,
      qVertices0,
      sb_stream_set,
      streams_stop[6],
      t_stream_tuple,
      t_stream_sol,
      t_stream_sol_end);

    intersectcache_wrapper<PROPOSE_BATCH_LOG>(
        hTables1,
        htb_cache,
        t_stream_tuple,
        streams_stop[7],
        i_stream_tuple);

    mwj_compact(v_stream_tuple,
                streams_stop[8],
                c_stream_tuple);

    mwj_filter<(1UL << PROPOSE_BATCH_LOG)>(c_stream_tuple,
                                           streams_stop[9],
                                           f_stream_set);

    verifycache_wrapper<PROPOSE_BATCH_LOG>(
        hTables1,
        htb_cache,
        bb_merge_stream_tuple,
        streams_stop[10],
        v_stream_tuple);


    mwj_assembly(
        nQueryVer,
        f_stream_set,
        f_stream_sol,
        f_stream_sol_end,
        stream_batch,
        stream_batch_end,
        streams_stop,
        // stream_sol0,
        // merge_in[1],
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
    std::thread mwj_propose_t(mwj_propose,
                              nQueryVer,
                              std::ref(dyn_stream_sol),
                              std::ref(streams_stop[0]),
                              std::ref(p0_stream_sol),
                              std::ref(p0_stream_sol_end));

    std::thread mwj_edgebuild_t(mwj_edgebuild<LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>,
                                hash1_w,
                                qVertices0,
                                std::ref(p0_stream_sol),
                                std::ref(p0_stream_sol_end),
                                std::ref(streams_stop[1]),
                                std::ref(e_stream_tuple),
                                std::ref(e_stream_sol),
                                std::ref(e_stream_sol_end));

    std::thread mwj_findmin_t(mwj_findmin<T_BLOOM, BLOOM_LOG, K_FUN_LOG>,
                              bloom_p,
                              // std::ref(bloom_cache),
                              std::ref(e_stream_tuple),
                              std::ref(streams_stop[2]),
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
                                      std::ref(streams_stop[3]),
                                      std::ref(rc_stream_tuple),
                                      std::ref(rc_stream_filter));

    std::thread mwj_readmin_edge_t(
      mwj_readmin_edge<T_BLOOM, BLOOM_LOG, K_FUN_LOG, FULL_HASH_W>,
      htb_buf2,
      std::ref(rc_stream_tuple),
      std::ref(rc_stream_filter),
      std::ref(streams_stop[4]),
      std::ref(re_stream_set),
      std::ref(re_stream_tuple));
    
    std::thread mwj_homomorphism_t(
        mwj_homomorphism,
        std::ref(re_stream_set),
        std::ref(re_stream_tuple),
        std::ref(re_stream_sol),
        std::ref(re_stream_sol_end),
        std::ref(streams_stop[5]),
        std::ref(h_stream_set));

    std::thread mwj_tuplebuild_t(
      mwj_tuplebuild<PROPOSE_BATCH_LOG, LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>,
      hash1_w,
      hash2_w,
      qVertices0,
      std::ref(sb_stream_set),
      std::ref(streams_stop[6]),
      std::ref(t_stream_tuple),
      std::ref(t_stream_sol),
      std::ref(t_stream_sol_end));

    std::thread mwj_intersect_t(
        mwj_intersect<PROPOSE_BATCH_LOG>,
        hTables1,
        std::ref(htb_cache),
        std::ref(t_stream_tuple),
        std::ref(streams_stop[7]),
        std::ref(i_stream_tuple));
    
    std::thread mwj_compact_t(
        mwj_compact,
        std::ref(v_stream_tuple),
        std::ref(streams_stop[8]),
        std::ref(c_stream_tuple));

    std::thread mwj_filter_t(
        mwj_filter<(1UL << PROPOSE_BATCH_LOG)>,
        std::ref(c_stream_tuple),
        std::ref(streams_stop[9]),
        std::ref(f_stream_set));

    std::thread mwj_verify_t(
        mwj_verify<PROPOSE_BATCH_LOG>,
        hTables1,
        std::ref(htb_cache),
        std::ref(bb_merge_stream_tuple),
        std::ref(streams_stop[10]),
        std::ref(v_stream_tuple));

    std::thread mwj_assembly_t(
        mwj_assembly,
        nQueryVer,
        std::ref(f_stream_set),
        std::ref(f_stream_sol),
        std::ref(f_stream_sol_end),
        std::ref(stream_batch),
        std::ref(stream_batch_end),
        std::ref(streams_stop),
        std::ref(a_stream_sol),
        std::ref(result));

    mwj_propose_t.join();
    mwj_edgebuild_t.join();
    mwj_findmin_t.join();
    mwj_readmin_counter_t.join();
    mwj_readmin_edge_t.join();
    mwj_homomorphism_t.join();
    mwj_tuplebuild_t.join();
    mwj_intersect_t.join(); 
    mwj_verify_t.join();
    mwj_compact_t.join();
    mwj_filter_t.join();
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

void
subgraphIsomorphism(row_t htb_buf0[HASHTABLES_SPACE],
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
                    unsigned long& dynfifo_diagnostic,

#if DEBUG_INTERFACE
                    volatile unsigned int& debif_endpreprocess,
                    unsigned long& p_propose_empty,
                    unsigned long& p_edgebuild_empty,
                    unsigned long& p_findmin_empty,
                    unsigned long& p_readmin_counter_empty,
                    unsigned long& p_readmin_edge_empty,
                    unsigned long& p_homomorphism_empty,
                    unsigned long& p_batchbuild_empty,
                    unsigned long& p_tuplebuild_empty,
                    unsigned long& p_intersect_empty,
                    unsigned long& p_offset_empty,
                    unsigned long& p_split_empty,
                    unsigned long& p_verify_empty,
                    unsigned long& p_compact_empty,
                    unsigned long& p_filter_empty,
                    unsigned long& p_assembly_empty,
#endif /* DEBUG_INTERFACE */

#if COUNT_ONLY
                    long unsigned int& result
#else
                     hls::stream<T_NODE>& result
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
    latency=1 max_read_burst_length=32
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
#pragma HLS INTERFACE mode=s_axilite port=p_propose_empty
#pragma HLS INTERFACE mode=s_axilite port=p_edgebuild_empty
#pragma HLS INTERFACE mode=s_axilite port=p_findmin_empty
#pragma HLS INTERFACE mode=s_axilite port=p_readmin_counter_empty
#pragma HLS INTERFACE mode=s_axilite port=p_readmin_edge_empty
#pragma HLS INTERFACE mode=s_axilite port=p_tuplebuild_empty
#pragma HLS INTERFACE mode=s_axilite port=p_intersect_empty
#pragma HLS INTERFACE mode=s_axilite port=p_verify_empty
#pragma HLS INTERFACE mode=s_axilite port=p_assembly_empty
#pragma HLS INTERFACE mode=s_axilite port=p_homomorphism_empty
#pragma HLS INTERFACE mode=s_axilite port=p_batchbuild_empty
#pragma HLS INTERFACE mode=s_axilite port=p_offset_empty 
#pragma HLS INTERFACE mode=s_axilite port=p_split_empty
#pragma HLS INTERFACE mode=s_axilite port=p_compact_empty
#pragma HLS INTERFACE mode=s_axilite port=p_filter_empty 
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
               MAX_TABLES>(res_buf,
                           htb_buf0,
                           bloom_p,
                           qVertices0,
                           qVertices1,
                           hTables0,
                           hTables1,
                           dynfifo_space,
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

    ap_wait();
    result = localResult;
    p_propose_empty = propose_empty;
    p_edgebuild_empty = edgebuild_empty;
    p_findmin_empty = findmin_empty;
    p_readmin_counter_empty = readmin_counter_empty;
    p_readmin_edge_empty = readmin_edge_empty;
    p_homomorphism_empty = homomorphism_empty;
    p_batchbuild_empty = batchbuild_empty;
    p_tuplebuild_empty = tuplebuild_empty;
    p_intersect_empty = intersect_empty;
    p_offset_empty = offset_empty;
    p_split_empty = split_empty;
    p_verify_empty = verify_empty;
    p_compact_empty = compact_empty;
    p_filter_empty = filter_empty;
    p_assembly_empty = assembly_empty;

#if DEBUG_STATS
    debug::print(hash1_w, hash2_w);
#endif /* DEBUG_STATS */
}
#endif /* SOFTWARE_PREPROC */

#pragma GCC diagnostic pop
