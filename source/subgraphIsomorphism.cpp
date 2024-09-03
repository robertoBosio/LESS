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
/* #pragma GCC diagnostic error "-Wall" */
#pragma GCC diagnostic error "-Wextra"
#pragma GCC diagnostic ignored "-Wunused-label"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wc++11-compat"
// #pragma GCC diagnostic ignored "-Wunused-variable"
// #pragma GCC diagnostic ignored "-Wunused-but-set-variable"

#define STOP_S      2   
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
unsigned long reqs_findmin = 0;
unsigned long reqs_readmin_counter = 0;
unsigned long reqs_readmin_edge = 0;
unsigned long reqs_verify = 0;
unsigned long reqs_intersect = 0;
unsigned long hits_findmin = 0;
unsigned long hits_readmin_counter = 0;
unsigned long hits_readmin_edge = 0;
unsigned long hits_intersect = 0;
unsigned long hits_verify = 0;
unsigned long bloom_filtered = 0;
#endif

typedef cache<ap_uint<DDR_W>,
              true,
              false,
              2,
              HASHTABLES_SPACE,
              0,
              0,
              (1UL << CACHE_WORDS_PER_LINE),
              false,
              4096,
              1,
              false,
              1,
              AUTO,
              URAM>
  htb_cache_t;

typedef cache<ap_uint<DDR_W>,
              true,
              false,
              1,
              HASHTABLES_SPACE,
              0,
              0,
              16,
              false,
              512,
              1,
              false,
              1,
              AUTO,
              BRAM>
  htb_cache2_t;

/******** Tuple definition ********/
enum edge_flag
{
  MIN_SET = 1,
  NO_EDGE = 2,
  CHECK = 0
};
typedef ap_uint<2> edge_flag_type;

template<typename NODE_T>
struct sol_node_t
{
  NODE_T node;
  unsigned char pos;
  bool last;
  bool stop;
};

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
  bool stop; // Stop flag
} findmin_tuple_t;

typedef struct
{
  /* Minimum size set information (edge), no need for the delimeter
  this time since the min_set can be only one */
  ap_uint<V_ID_W> indexing_v;   // Indexing vertex in the edge
  unsigned char tb_index;       // Table index
  unsigned char iv_pos;         // Query indexing vertex position.
  unsigned char num_tb_indexed; // Used by sequencebuild
  ap_uint<64> addr_counter;
  bool skip_counter;
  bool stop;
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
  unsigned int cycles;
  bool stop;
} readmin_edge_tuple_t;

template<typename NODE_T>
struct homomorphism_set_t
{
  NODE_T node;
  bool last;
  bool valid;
};

template<typename NODE_T>
struct sequencebuild_set_t
{
  /* Adding the flag for min set since homomorphism aggregate tuple and set in
   * a unique stream to let sequencebuild be a simple pipeline */
  NODE_T node;          // store the node id or the data of the min set
  unsigned char pos;
  bool last;            // last of the set flag
  bool sol;             // solution node flag
  bool min_set;         // min set information flag
  bool stop;            // stop flag
};

template<typename NODE_T>
struct sequencebuild_tuple_t
{
  NODE_T node;
  unsigned char tb_index;
  unsigned char iv_pos;
  unsigned char query_edge;
  unsigned char pos;
  bool last_edge;       // last edge of the node flag
  bool sol;             // solution node flag
  bool stop;            // stop flag

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
  bool stop;
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
  bool stop;
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
  bool stop;
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
  bool stop;
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
  bool stop;
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
  bool stop;
} compact_tuple_t;

typedef struct
{
  ap_uint<V_ID_W> indexed_v;
  ap_uint<PROPOSE_BATCH_LOG> pos;
  bool bit_last_edge;
  bool bit_checked;
  bool last_set;
  bool last_batch;
  bool stop;
} assembly_tuple_t;

typedef struct
{
    ap_uint<V_ID_W> node;
    bool last_set;
    bool stop;
} assembly_set_t;

template<typename NODE_T>
struct assembly_node_t
{
  NODE_T node;
  unsigned char pos;
  bool last;
  bool sol;
  bool stop;
};

/******** End tuple definition ********/

/* Count the number of 1s in a bloom filter */
template<typename T_BLOOM, size_t BLOOM_LOG, size_t K_FUN_LOG>
unsigned int
bloom_bitset(T_BLOOM filter)
{
  unsigned int count{ 0 };
  for (int c = 0; c < (1UL << (BLOOM_LOG - 5)); c++) {
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
bloom_intersect(T_BLOOM& filter, T_BLOOM set_bloom)
{
#pragma HLS inline
  unsigned short bloom_s = bloom_bitset<T_BLOOM, BLOOM_LOG, 0>(set_bloom);
  filter = filter & set_bloom;
  return bloom_s;
}

template<size_t W>
void
hash_wrapper(ap_uint<W> key_val, ap_uint<64>& hash_val)
{
#pragma HLS inline off
  xf::database::details::hashlookup3_core<W>(key_val, hash_val);
}

template<typename T_BLOOM, size_t BLOOM_LOG, size_t K_FUN, size_t FULL_HASH_W>
void
bloom_test(T_BLOOM filter[K_FUN], ap_uint<64> hash_val, bool& test)
{
#pragma HLS inline off
  for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
    ap_uint<BLOOM_LOG> idx =
      hash_val.range((FULL_HASH_W / K_FUN) * (g + 1) - 1,
                     (FULL_HASH_W / K_FUN) * (g + 1) - BLOOM_LOG);
    test = test && (filter[g][idx] == 1);
  }
}

/* Reads the partial solutions in the format of radix followed by a set of
 * extensions, and streams them out registering which position each node has in
 * the query vertex order */
void
mwj_propose(hls::stream<ap_uint<V_ID_W>>& stream_fifo_in,
            hls::stream<sol_node_t<vertex_t>>& stream_sol_out)
{
#pragma HLS pipeline II = 1 style = flp
  const ap_uint<V_ID_W> MASK_RADIX = ~(1UL << (V_ID_W - 1));
  const vertex_t STOP_NODE = ~0;
  const vertex_t FAKE_NODE = ~0 - 1;
  static ap_uint<V_ID_W> n_nodes = 0;
  ap_uint<V_ID_W> vertex_read;
  sol_node_t<vertex_t> vertex;

  /* Read BFS solutions from dynfifo. Solution are read as radix and extension.
   * The radix is marked by the 31th bit being asserted. The motivation behind
   * is to not always resend all the new solution but instead send only the
   * changed vertex.
   * STOP_NODE stops the entire pipeline.
   * FAKE_NODE used in case of radix without extension: single node solutions */
  vertex_read = stream_fifo_in.read();
  if (vertex_read == FAKE_NODE) {
    n_nodes = 0;
  } else {
    bool radix = vertex_read.test(V_ID_W - 1);
    vertex.pos = n_nodes;
    vertex.node = vertex_read & MASK_RADIX;
    vertex.last = !radix || (vertex_read == STOP_NODE);
    vertex.stop = (vertex_read == STOP_NODE);
    if (radix)
      n_nodes++;
    stream_sol_out.write(vertex);
  }
}

/* Edgebuild read the next vertex to be mapped and extract the information of
 * the sets involved in the intersection, computing also the addresses of the
 * bloom filters */
template<size_t LKP3_HASH_W, size_t MAX_HASH_W, size_t FULL_HASH_W>
void
mwj_edgebuild(const unsigned char hash1_w,
              QueryVertex* qVertices,
              hls::stream<sol_node_t<vertex_t>>& stream_sol_in,
              hls::stream<findmin_tuple_t>& stream_tuple_out,
              hls::stream<sol_node_t<vertex_t>>& stream_sol_out)
{
  ap_uint<V_ID_W> curEmb[MAX_QV];
  unsigned char curQV = 0;
  sol_node_t<vertex_t> vertex;
  findmin_tuple_t tuple_out;

  /* Initializing local filter in findmin */
  tuple_out.reset = true;
  tuple_out.last = false;
  tuple_out.stop = false;
  tuple_out.address = 0;
  stream_tuple_out.write(tuple_out);

EDGEBUILD_TASK_LOOP:
  while (true) {

  EDGEBUILD_COPYING_EMBEDDING_LOOP:
    do {
#pragma HLS pipeline II = 1
      vertex = stream_sol_in.read();
      curEmb[vertex.pos] = vertex.node;
      stream_sol_out.write(vertex);
    } while (!vertex.last);

    if (vertex.stop) {
      break;
    }

    curQV = vertex.pos + 1;
  EDGEBUILD_MAIN_LOOP:
    for (int g = 0; g < qVertices[curQV].numTablesIndexed; g++) {
#pragma HLS pipeline II = 1
      unsigned char tb_index = qVertices[curQV].tables_indexed[g];
      unsigned char iv_pos = qVertices[curQV].vertex_indexing[g];

      // Computing addresses of indexed sets
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
      tuple_out.num_tb_indexed = qVertices[curQV].numTablesIndexed;
      tuple_out.last = (g == (qVertices[curQV].numTablesIndexed - 1));
      stream_tuple_out.write(tuple_out);
    }
    tuple_out.reset = true;
    tuple_out.last = false;
    stream_tuple_out.write(tuple_out);
  }

  /* Propagate stop node on two different routes */
  tuple_out.stop = true;
  stream_tuple_out.write(tuple_out);
}

/* Finds the set with the smallest cardinality in the intersection. Intersects
 * the bloom filters */
template<typename T_BLOOM,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t FULL_HASH_W>
void
mwj_findmin(bloom_t* bloom_p,
            const unsigned char hash1_w,
            const unsigned char hash2_w,
            hls::stream<findmin_tuple_t>& stream_tuple_in,
            hls::stream<readmin_counter_tuple_t> stream_tuple_out[2],
            hls::stream<T_BLOOM> stream_filter_out[1 << K_FUN_LOG])
{
  constexpr size_t K_FUN = (1UL << K_FUN_LOG);
  T_BLOOM filter[K_FUN];
  readmin_counter_tuple_t tuple_out;
  findmin_tuple_t tuple_in;
  unsigned int min_size = ~0;
  unsigned long addr_counter;
  ap_uint<LKP3_HASH_W> hash_out;
  ap_uint<MAX_HASH_W> hash_trimmed;
#pragma HLS array_partition variable = filter type = complete dim = 1

  tuple_out.stop = false;

FINDMIN_TASK_LOOP:
  while (true) {
#pragma HLS pipeline II = 4

    /* Non blocking read since Vitis 2022.2 does not want to apply flushable
     * pipeline */
    if (stream_tuple_in.read_nb(tuple_in)) {
      unsigned int address = tuple_in.address;
      address <<= K_FUN_LOG;
      unsigned short bloom_s = 0;

      if (tuple_in.stop) {
        break;
      } else if (tuple_in.reset) {
        for (int s = 0; s < K_FUN; s++) {
#pragma HLS unroll
          filter[s] = ~0;
        }
      } else {
        for (int s = 0; s < K_FUN; s++) {
#pragma HLS unroll
          T_BLOOM set_bloom = bloom_p[address + s];
          bloom_s += bloom_intersect<T_BLOOM, BLOOM_LOG>(filter[s], set_bloom);
        }
        reqs_findmin++;
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
        xf::database::details::hashlookup3_core<V_ID_W>(tuple_out.indexing_v,
                                                        hash_out);
        hash_trimmed = hash_out;
        hash_trimmed = hash_trimmed.range(hash1_w - 1, 0);

        addr_counter = hash_trimmed - 1;
        addr_counter <<= hash2_w;
        addr_counter += (1UL << hash2_w) - 1;
        tuple_out.addr_counter = addr_counter;
        if (hash_trimmed == 0) {
          tuple_out.skip_counter = true;
        }
        stream_tuple_out[0].write(tuple_out);

        addr_counter = hash_trimmed;
        addr_counter <<= hash2_w;
        addr_counter += (1UL << hash2_w) - 1;
        tuple_out.addr_counter = addr_counter;
        tuple_out.skip_counter = false;

        stream_tuple_out[1].write(tuple_out);

        for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
          stream_filter_out[g].write(filter[g]);
        }
      }
    }
  }

  /* Propagate stop node */
  tuple_out.stop = true;
  stream_tuple_out[0].write(tuple_out);
}

void
mwj_bypassfilter(hls::stream<bloom_t> stream_filter_in[1UL << K_FUNCTIONS],
                 hls::stream<bloom_t> stream_filter_out[1UL << K_FUNCTIONS])
{
  for (int g = 0; g < (1UL << K_FUNCTIONS); g++) {
#pragma HLS unroll
    stream_filter_out[g].write(stream_filter_in[g].read());
  }
}

void
mwj_readmin_counter(AdjHT* hTables,
                    row_t* m_axi,
                    hls::stream<readmin_counter_tuple_t> stream_tuple_in[2],
                    hls::stream<readmin_edge_tuple_t>& stream_tuple_out)
{
  readmin_counter_tuple_t tuple_in;
  readmin_edge_tuple_t tuple_out;
  unsigned char stream_p = 0;
  tuple_out.stop = false;

READMIN_COUNTER_TASK_LOOP:
  while (true) {
#pragma HLS pipeline II = 1
    
    if (stream_tuple_in[stream_p].read_nb(tuple_in)){
      
      if (tuple_in.stop){
        break;
      }

      volatile unsigned int offset = 0;
      if (!tuple_in.skip_counter)
      {
        ap_uint<DDR_BIT - C_W> addr_inrow;
        ap_uint<DDR_W> ram_row;
        unsigned long addr_row;

        /* Compute address of row storing the counter */
        addr_row = hTables[tuple_in.tb_index].start_offset +
                   (tuple_in.addr_counter >> (DDR_BIT - C_W));

        /* Compute address of data inside the row */
        addr_inrow = tuple_in.addr_counter.range((DDR_BIT - C_W) - 1, 0);

        /* Read the data */
        ram_row = m_axi[addr_row];
        if (addr_inrow == 0)
        {
          offset = ram_row.range((1UL << C_W) - 1, 0);
        }
        else if (addr_inrow == 1)
        {
          offset = ram_row.range((2UL << C_W) - 1, 1UL << C_W);
        }
        else if (addr_inrow == 2)
        {
          offset = ram_row.range((3UL << C_W) - 1, 2UL << C_W);
        }
        else
        {
          offset = ram_row.range((4UL << C_W) - 1, 3UL << C_W);
        }
        reqs_readmin_counter++;
#if DEBUG_STATS
        debug::readmin_counter_reads++;
#endif
      }

      unsigned int row =
        hTables[tuple_in.tb_index].start_edges + (offset >> (DDR_BIT - E_W));

      tuple_out.indexing_v = tuple_in.indexing_v;
      tuple_out.tb_index = tuple_in.tb_index;
      tuple_out.iv_pos = tuple_in.iv_pos;
      tuple_out.num_tb_indexed = tuple_in.num_tb_indexed;

      if (stream_p == 0){
        tuple_out.rowstart = row;
      } else {
        tuple_out.cycles = row - tuple_out.rowstart;
        stream_tuple_out.write(tuple_out);
      }
      stream_p = (stream_p + 1) % 2;
    }
  }

  /* Propagate stop node */
  tuple_out.stop = true;
  stream_tuple_out.write(tuple_out);
}

template<typename T_BLOOM,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t FULL_HASH_W>
void
mwj_readmin_edge_pipelined(
  htb_cache2_t& m_axi,
  hls::stream<readmin_edge_tuple_t>& stream_tuple_in,
  hls::stream<T_BLOOM> stream_filter_in[1UL << K_FUN_LOG],
  hls::stream<homomorphism_set_t<ap_uint<V_ID_W>>> stream_set_out[2],
  hls::stream<minset_tuple_t>& stream_tuple_out)
{
  constexpr size_t K_FUN = (1UL << K_FUN_LOG);
  T_BLOOM filter[K_FUN];
  readmin_edge_tuple_t tuple_in;
  bool read_new = true;
  unsigned int cycles = 0;
  unsigned int counter = 0;
#pragma HLS array_partition variable = filter type = complete dim = 1

READMIN_EDGE_TASK_LOOP:
  while (true) {
#pragma HLS pipeline II = 1

    if (read_new) {

      // Non blocking read since Vitis 2022.2 does not want to apply flushable
      if (stream_tuple_in.read_nb(tuple_in)) {
        if (tuple_in.stop) {
          break;
        }
        
        minset_tuple_t tuple_out;
        tuple_out.stop = false;
        tuple_out.tb_index = tuple_in.tb_index;
        tuple_out.iv_pos = tuple_in.iv_pos;
        tuple_out.num_tb_indexed = tuple_in.num_tb_indexed;
        stream_tuple_out.write(tuple_out);
        
        for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
          filter[g] = stream_filter_in[g].read();
        }
        read_new = false;
        cycles = tuple_in.cycles;
        counter = 0;
      }

    } else {

      if (counter == cycles) {
        read_new = true;
      }

      row_t row = m_axi[tuple_in.rowstart + counter];
      for (int i = 0; i < EDGE_ROW; i++) {
#pragma HLS unroll
        ap_uint<V_ID_W> indexing_v, indexed_v;
        ap_uint<FULL_HASH_W> hash_out;
        ap_uint<V_ID_W * 2> edge = row.range(((i + 1) << E_W) - 1, i << E_W);
        indexing_v = edge.range(V_ID_W * 2 - 1, V_ID_W);
        indexed_v = edge.range(V_ID_W - 1, 0);

        hash_wrapper<V_ID_W>(indexed_v, hash_out);
        bool test = true;
        bloom_test<T_BLOOM, BLOOM_LOG, K_FUN, FULL_HASH_W>(
          filter, hash_out, test);

        homomorphism_set_t<ap_uint<V_ID_W>> set_out;
        set_out.node = indexed_v;
        set_out.last = (i == EDGE_ROW - 1) && (counter == cycles);
        set_out.valid = test && tuple_in.indexing_v == indexing_v;
        stream_set_out[i].write(set_out);
        
        if (tuple_in.indexing_v == indexing_v) {
          if (!test) {
            bloom_filtered++;
          }
        }
      }
      counter++;

#if DEBUG_STATS
      debug::readmin_edge_reads += cycles + 1;
      debug::readmin_n_sets++;
#endif /* DEBUG_STATS */
    }
  }
}

void mwj_homomorphism(
    hls::stream<homomorphism_set_t<ap_uint<V_ID_W>>> stream_set_in[2],
    hls::stream<minset_tuple_t> &stream_tuple_in,
    hls::stream<sol_node_t<vertex_t>> &stream_sol_in,
    hls::stream<sequencebuild_set_t<ap_uint<V_ID_W>>> &stream_set_out)
{
  ap_uint<8> curQV;
  ap_uint<V_ID_W> curEmb[MAX_QV];
  ap_uint<V_ID_W> fake_node;
  ap_uint<MAX_QV> equal_bits;
  ap_uint<MAX_QV> valid_bits;
  homomorphism_set_t<ap_uint<V_ID_W>> set_in;
  sequencebuild_set_t<ap_uint<V_ID_W>> set_out;
  sol_node_t<vertex_t> vertex;

  while (true)
  {
    ap_uint<1> select = 0;

  HOMOMORPHISM_COPYING_EMBEDDING_LOOP:
    do
    {
#pragma HLS pipeline II = 1
      vertex = stream_sol_in.read();
      set_out.node = vertex.node;
      set_out.last = vertex.last;
      set_out.pos = vertex.pos;
      set_out.min_set = false;
      set_out.sol = true;
      set_out.stop = vertex.stop;
      stream_set_out.write(set_out);
      curEmb[vertex.pos] = vertex.node;
    } while (!vertex.last);
    curQV = vertex.pos + 1;

    if (vertex.stop)
    {
      break;
    }

    /* Fake node is the tuple about min_set data */
    minset_tuple_t tuple = stream_tuple_in.read();
    fake_node.range(7, 0) = tuple.tb_index;
    fake_node.range(15, 8) = tuple.iv_pos;
    fake_node.range(31, 16) = tuple.num_tb_indexed;
    set_out.node = fake_node;
    set_out.last = false;
    set_out.min_set = true;
    set_out.sol = false;
    set_out.stop = false;
    stream_set_out.write(set_out);
    valid_bits = (1UL << curQV) - 1;

  HOMOMORPHISM_CHECK_LOOP:
    do
    {
#pragma HLS pipeline II = 1
      set_in = stream_set_in[select].read();
      select++;
      ap_uint<V_ID_W> vToVerify = set_in.node;
      equal_bits = 0;

      for (int g = 0; g < MAX_QV; g++)
      {
#pragma HLS unroll
        if (vToVerify == curEmb[g])
        {
          equal_bits[g] = 1;
        }
      }
      set_out.node = vToVerify;
      set_out.last = false;
      set_out.min_set = false;
      set_out.sol = false;

      /* Write out in case of not duplicate in current solution or
      delimeter node */
      if ((equal_bits & valid_bits) == 0 && set_in.valid)
      {
        stream_set_out.write(set_out);
        // std::cout << "Sending node " << vToVerify << std::endl;
      }
#if DEBUG_STATS
      else
      {
        debug::homomo_trashed++;
      }
#endif
    } while (!set_in.last);
    
    // std::cout << "Propagate last set" << std::endl;
    set_out.last = true;
    stream_set_out.write(set_out);
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
  sequencebuild_set_t<ap_uint<V_ID_W>> tuple_in;
  sequencebuild_tuple_t<ap_uint<V_ID_W>> tuple_out;
  bool last_inner;

  if (query_edge == 0) {
    tuple_in = stream_tuple_in.read();

    /* Two different cases of termination, one is given by the last
    delimeter and the other is given by the overflow of the batch size*/
    last_inner = tuple_in.last || (buffer_p == ((1UL << BATCH_SIZE_LOG) - 1));

    vToVerify = tuple_in.node;
    buffer[buffer_p] = vToVerify;

    if (tuple_in.min_set) {
      tb_index = vToVerify.range(7, 0);
      iv_pos = vToVerify.range(15, 8);
      cycles = vToVerify.range(31, 16);
    } else if (tuple_in.sol) {
      tuple_out.node = tuple_in.node;
      tuple_out.tb_index = 0;
      tuple_out.iv_pos = 0;
      tuple_out.query_edge = 0;
      tuple_out.pos = tuple_in.pos;
      tuple_out.last_edge = tuple_in.last;
      tuple_out.sol = true;
      tuple_out.stop = tuple_in.stop;
      tuple_out.last_set = false;
      tuple_out.last_batch = false;
      stream_tuple_out.write(tuple_out);
    } else {
      tuple_out.node = vToVerify;
      tuple_out.tb_index = tb_index;
      tuple_out.iv_pos = iv_pos;
      tuple_out.query_edge = query_edge;
      tuple_out.pos = buffer_p;
      tuple_out.last_edge = (cycles == 1);
      tuple_out.sol = false;
      tuple_out.stop = tuple_in.stop;
      tuple_out.last_set = tuple_in.last;
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
        buffer_p = 0;
        if (cycles == 1) {
          /* Checking outer loop end condition */
          buffer_size = 0;
          query_edge = 0;
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
      tuple_out.query_edge = query_edge;
      tuple_out.pos = buffer_p;
      tuple_out.last_edge = (query_edge == (cycles - 1));
      tuple_out.sol = false;
      tuple_out.stop = false;
      tuple_out.last_set = false;
      tuple_out.last_batch = (buffer_p == ((1UL << BATCH_SIZE_LOG) - 1)) &&
                             query_edge == (cycles - 1);
      stream_tuple_out.write(tuple_out);
      buffer_p++;
    } else {

      /* Checking inner loop end condition.
      Inner loop cycles on nodes over an edge; Doing single node over
      the edges would have been much more easier but this rearrangment
      achieve better hit ratios since each edge is a single table */
      query_edge++;
      buffer_p = 0;
      if (query_edge == cycles) {
        /* Checking outer loop end condition */
        query_edge = 0;
        if (buffer_size != (1UL << BATCH_SIZE_LOG)) {
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
  hls::stream<intersect_tuple_t> stream_tuple_out[2],
  hls::stream<sol_node_t<vertex_t>>& stream_sol_out)
{
  ap_uint<V_ID_W> curEmb[MAX_QV];
  ap_uint<V_ID_W> vToVerify;
  hls::stream<ap_uint<V_ID_W>, 4> hash_in0, hash_in1;
  hls::stream<ap_uint<LKP3_HASH_W>, 4> hash_out0, hash_out1;
  unsigned long addr_counter;
  intersect_tuple_t tuple_out;
  sequencebuild_tuple_t<ap_uint<V_ID_W>> tuple_in;
  unsigned char curQV = 0;

  tuple_out.stop = false;

TUPLEBUILD_TASK_LOOP:
  while (true) {
#pragma HLS pipeline II = 1
    if (stream_tuple_in.read_nb(tuple_in)) {

      if (tuple_in.sol) {
        curEmb[tuple_in.pos] = tuple_in.node;
        stream_sol_out.write(
          { tuple_in.node, tuple_in.pos, tuple_in.last_edge, tuple_in.stop });
        if (tuple_in.stop) {
          break;
        }
        curQV = tuple_in.pos + 1;
      } else {
        uint8_t tableIndex =
          qVertices[curQV].tables_indexed[tuple_in.query_edge];
        uint8_t ivPos =
          qVertices[curQV].vertex_indexing[tuple_in.query_edge];

        bool bit_last = tuple_in.last_edge;
        bool bit_min =
          (tuple_in.tb_index == tableIndex && tuple_in.iv_pos == ivPos);

        vToVerify = tuple_in.node;
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
  }

  /* Propagate stop node */
  tuple_out.stop = true;
  stream_tuple_out[0].write(tuple_out);
}

void
mwj_bypass_sol(hls::stream<sol_node_t<vertex_t>>& stream_sol_in,
               hls::stream<sol_node_t<vertex_t>>& stream_sol_out)
{
#pragma HLS pipeline II = 1
  stream_sol_out.write(stream_sol_in.read());
}

template<size_t BATCH_SIZE_LOG>
void
mwj_intersect(AdjHT* hTables,
              htb_cache_t& htb_buf,
              hls::stream<intersect_tuple_t> stream_tuple_in[2],
              hls::stream<split_tuple_t>& stream_tuple_out)
{
  constexpr size_t EDGE_BLOCK = (CACHE_WORDS_PER_LINE + DDR_BIT - E_W);
  intersect_tuple_t tuple_in;
  split_tuple_t tuple_out;
  ap_uint<(1UL << C_W)> end_off;
  unsigned char tableIndex;
  ap_uint<DDR_W> ram_row;
  unsigned long addr_row;
  ap_uint<DDR_BIT - C_W> addr_inrow;
  ap_uint<64> addr_counter;
  unsigned char stream_p_in = 0;
  tuple_out.stop = false;

INTERSECT_TASK_LOOP:
  while (true)
  {
#pragma HLS pipeline II = 1
    if (stream_tuple_in[stream_p_in].read_nb(tuple_in))
    {

      ap_uint<(1UL << C_W)> offset = 0;
      if (tuple_in.stop)
      {
        break;
      }
      else if (!tuple_in.last_set)
      {
        tableIndex = tuple_in.tb_index;
        addr_counter = tuple_in.addr_counter;

        if (tuple_in.flag == CHECK && !tuple_in.skip_counter)
        {

          /* Compute address of row storing the counter */
          addr_row = hTables[tableIndex].start_offset +
                     (addr_counter >> (DDR_BIT - C_W));

          /* Compute address of data inside the row */
          addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);
          ram_row = htb_buf.get(addr_row, 1);
          if (addr_inrow == 0)
          {
            offset = ram_row.range((1UL << C_W) - 1, 0);
          }
          else if (addr_inrow == 1)
          {
            offset = ram_row.range((2UL << C_W) - 1, 1UL << C_W);
          }
          else if (addr_inrow == 2)
          {
            offset = ram_row.range((3UL << C_W) - 1, 2UL << C_W);
          }
          else
          {
            offset = ram_row.range((4UL << C_W) - 1, 3UL << C_W);
          }

#if DEBUG_STATS
          debug::intersect_reads += 1;
#endif /* DEBUG_STATS */
        }

        tuple_out.indexed_v = tuple_in.indexed_v;
        tuple_out.indexing_v = tuple_in.indexing_v;
        tuple_out.tb_index = tuple_in.tb_index;
        tuple_out.pos = tuple_in.pos;
        tuple_out.bit_last_edge = tuple_in.bit_last_edge;
        tuple_out.flag = tuple_in.flag;
      }
      tuple_out.last_batch = tuple_in.last_batch;
      tuple_out.last_set = tuple_in.last_set;

      if (stream_p_in == 0)
      {
        end_off = offset;
      }
      else
      {
        tuple_out.start_off = offset;
        tuple_out.first_block = offset >> EDGE_BLOCK;
        tuple_out.end_block = (offset == end_off)
                                  ? (unsigned int)(end_off >> EDGE_BLOCK)
                                  : (unsigned int)((end_off - 1) >> EDGE_BLOCK);
        tuple_out.flag = (tuple_in.flag == MIN_SET)
                             ? MIN_SET
                             : ((offset == end_off) ? NO_EDGE : CHECK);
        
        // Writing and switching output stream
        stream_tuple_out.write(tuple_out);

#if DEBUG_STATS
        debug::intersect_filter += (tuple_out.flag == NO_EDGE) ? 1 : 0;
#endif /* DEBUG_STATS */
      }

      // Switching input stream
      stream_p_in = (stream_p_in + 1) % 2;
    }
  }

  /* Propagating stop node */
  tuple_out.stop = true;
  stream_tuple_out.write(tuple_out);
}

template<size_t MAX_BATCH_SIZE, size_t ID>
void
mwj_blockbuild_pipelined(hls::stream<split_tuple_t>& stream_tuple_in,
                         hls::stream<verify_tuple_t>& stream_tuple_out)
{
#pragma HLS pipeline II = 1
  constexpr size_t EDGE_BLOCK = (CACHE_WORDS_PER_LINE + DDR_BIT - E_W);
  static split_tuple_t tuple_in;
  static bool read_new = true;
  static unsigned int s = 0;

  if (read_new) {
    tuple_in = stream_tuple_in.read();
    s = 0;
    read_new = false;
  } else {
    if (!tuple_in.last_set) {
      verify_tuple_t tuple_out;
      tuple_out.indexing_v = tuple_in.indexing_v;
      tuple_out.indexed_v = tuple_in.indexed_v;
      tuple_out.tb_index = tuple_in.tb_index;
      tuple_out.pos = tuple_in.pos;
      tuple_out.bit_last_edge = tuple_in.bit_last_edge;
      tuple_out.flag = tuple_in.flag;
      tuple_out.last_batch = tuple_in.last_batch;
      tuple_out.last_set = false;
      tuple_out.stop = tuple_in.stop;

      tuple_out.address = tuple_in.start_off + s;
      tuple_out.bit_last_address = (tuple_in.first_block == tuple_in.end_block);
      stream_tuple_out.write(tuple_out);
      if (tuple_in.first_block == tuple_in.end_block) {
        read_new = true;
      }
      tuple_in.first_block++;
      s += (1UL << EDGE_BLOCK);
    } else {
      /* Setting bit_last_address even for delimeter nodes for
      the merge_split to change stream on which reads. In non delimeter
      nodes, bit_last_address is used to identify the last tuple generated */
      verify_tuple_t tuple_out;
      tuple_out.stop = tuple_in.stop;
      tuple_out.last_set = tuple_in.last_set;
      tuple_out.last_batch = tuple_in.last_batch;
      tuple_out.bit_last_address = true;
      stream_tuple_out.write(tuple_out);
      read_new = true;
    }
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
  if (tuple_in.bit_last_address) {
    pointer = (pointer + 1) % NUM_SPLIT;
  }
}

template<size_t BATCH_SIZE_LOG>
void
mwj_verify(AdjHT* hTables,
           htb_cache_t& htb_buf,
           hls::stream<verify_tuple_t>& stream_tuple_in,
           hls::stream<compact_tuple_t>& stream_tuple_out)
{
  constexpr size_t EDGE_PER_WORD = (DDR_BIT - E_W);

  verify_tuple_t tuple_in;
  compact_tuple_t tuple_out;
  ap_uint<V_ID_W> candidate_v;
  ap_uint<V_ID_W> indexing_v;
  unsigned char tableIndex;
  ap_uint<128> edge_block[(1UL << CACHE_WORDS_PER_LINE)];
#pragma HLS array_partition variable = edge_block type = complete
  tuple_out.stop = false;

VERIFY_TASK_LOOP:
  while (true) {
#pragma HLS pipeline II = 1

    if (stream_tuple_in.read_nb(tuple_in)) {
      if (tuple_in.stop) {
        break;
      } else if (!tuple_in.last_set) {
        candidate_v = tuple_in.indexed_v;
        tuple_out.bit_equal = (tuple_in.flag == MIN_SET);

        if (tuple_in.flag == CHECK) {
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
          for (int g = 0; g < (1UL << CACHE_WORDS_PER_LINE); g++) {
#pragma HLS unroll
            for (int s = 0; s < (1UL << (EDGE_PER_WORD)); s++) {
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
  }
}

/* Or reduce of bits coming from verify, compact edge blocks and output if a
 * single edge is verified or not */
void
mwj_compact(hls::stream<compact_tuple_t>& stream_tuple_in,
            hls::stream<assembly_tuple_t>& stream_tuple_out)
{
#pragma HLS pipeline II = 1
  static bool checked = false;
  compact_tuple_t tuple_in;
  assembly_tuple_t tuple_out;

  tuple_in = stream_tuple_in.read();
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
    tuple_out.stop = tuple_in.stop;
    stream_tuple_out.write(tuple_out);
    checked = false;
  }
}

template<size_t MAX_BATCH_SIZE>
void
mwj_filter(hls::stream<assembly_tuple_t>& stream_tuple_in,
           hls::stream<assembly_set_t>& stream_set_out)
{
#pragma HLS pipeline II = 1
  static ap_uint<MAX_BATCH_SIZE> bits = ~0;
  assembly_tuple_t tuple_in;
  assembly_set_t tuple_out;

  tuple_in = stream_tuple_in.read();
  if (!tuple_in.last_set) {
    unsigned short p = tuple_in.pos;
    bits[p] = bits[p] && tuple_in.bit_checked;
    tuple_out.node = tuple_in.indexed_v;
    tuple_out.last_set = false;
    if (tuple_in.bit_last_edge && bits.test(p)) {
      stream_set_out.write(tuple_out);
    }
    if (tuple_in.last_batch && tuple_in.bit_last_edge) {
      bits = ~0;
    }
  } else {
    bits = ~0;
    tuple_out.last_set = tuple_in.last_set;
    tuple_out.stop = tuple_in.stop;
    stream_set_out.write(tuple_out);
  }
}

void
mwj_enlarge_sol(hls::stream<sol_node_t<vertex_t>>& stream_sol_in,
               hls::stream<sol_node_t<vertex_t>>& stream_sol_out)
{
  static vertex_t curEmb[MAX_QV];
  sol_node_t<vertex_t> vertex;
  bool stop;
ENLARGE_SOL_STORING_EMBEDDINGS_LOOP:
  do {
#pragma HLS pipeline II = 1
    vertex = stream_sol_in.read();
    curEmb[vertex.pos] = vertex.node;
    stop = vertex.stop;
  } while (!vertex.last);
  unsigned char curQV = vertex.pos + 1;

ENLARGE_SOL_WRITING_EMBEDDINGS_LOOP:
  for (int g = 0; g < curQV; g++) {
#pragma HLS pipeline II = 1
    vertex.node = curEmb[g];
    vertex.last = (g == (curQV - 1));
    vertex.pos = curQV;
    vertex.stop = stop;
    stream_sol_out.write(vertex);
  }
}

void
mwj_merge_solandset(hls::stream<sol_node_t<vertex_t>>& stream_sol_in,
                    hls::stream<assembly_set_t>& stream_set_in,
                    hls::stream<assembly_node_t<vertex_t>>& stream_sol_out)
{
#pragma HLS PIPELINE II = 1
  static bool select = false;
  static unsigned char pos = 0;
  sol_node_t<vertex_t> sol_vertex;
  assembly_set_t set_vertex;
  static bool printing_sol = false;

  bool last, sol, stop;
  vertex_t node;
  if (select) {
    set_vertex = stream_set_in.read();
    select = !set_vertex.last_set;
    last = set_vertex.last_set;
    node = set_vertex.node;
    sol = false;
    stop = false;

    // if (printing_sol) {
    //   if (last) {
    //     std::cout << std::endl;
    //     printing_sol = false;
    //   } else {
    //     std::cout << node << ", ";
    //   }
    // }
  } else {
    sol_vertex = stream_sol_in.read();
    select = sol_vertex.last;
    pos = sol_vertex.pos;
    last = false;
    node = sol_vertex.node;
    sol = true;
    stop = sol_vertex.stop;

    // if (sol_vertex.pos == 3){
    //   std::cout << sol_vertex.node << " ";
    //   printing_sol = true;
    // }
  }

  stream_sol_out.write({ node, pos, last, sol, stop });
}

void
mwj_assembly(row_t* m_axi,
             const unsigned int n_candidate,
             const unsigned int start_candidate,
             const unsigned int n_queryv,
             hls::stream<assembly_node_t<vertex_t>>& stream_sol_in,
             hls::stream<bool> streams_stop[STOP_S],
             hls::stream<ap_uint<V_ID_W>>& stream_partial_out,
             long unsigned int& result)
{
  const ap_uint<V_ID_W> MASK_NEW_SOLUTION = (1UL << (V_ID_W - 1));
  const ap_uint<V_ID_W> STOP_NODE = ~0;
  const ap_uint<V_ID_W> FAKE_NODE = ~0 - 1;
  unsigned long partial_sol = 0;
  // bool token_new_start;
  unsigned long int counter = 0;
  ap_uint<32> nodes_read = 0;
  bool stop = false;

ASSEMBLY_TASK_LOOP:
  do {
    // Test if there are some node from start batch
    if (nodes_read < n_candidate) {
      ap_uint<2> word_select = nodes_read.range(1, 0);
      ap_uint<32> row_select = nodes_read >> 2;
      row_t row = m_axi[start_candidate + row_select];
      ap_uint<V_ID_W> node = row.range((V_ID_W * (word_select + 1)) - 1,
                                       V_ID_W * word_select);
      
      /* False extension for single node solutions */
      stream_partial_out.write(FAKE_NODE);
      stream_partial_out.write(node);

      partial_sol = 1;
      nodes_read++;
    } else {
      stream_partial_out.write(STOP_NODE);
    }

  ASSEMBLY_SET_LOOP:
    do {
#pragma HLS PIPELINE II = 1
      assembly_node_t<vertex_t> vertex = stream_sol_in.read();
      if (vertex.stop){
        stop = true;
        break;
      }

      vertex_t dynfifo_node;
      if (vertex.sol){
        dynfifo_node = vertex.node | MASK_NEW_SOLUTION;
      } else if (vertex.last) {
        dynfifo_node = FAKE_NODE;
      } else {
        dynfifo_node = vertex.node;
      }

      if (vertex.pos < (n_queryv - 1)){
        stream_partial_out.write(dynfifo_node);
      } else if (!vertex.sol && !vertex.last) {
        counter++;
      }

      if (!vertex.sol){
        if (vertex.last){
          if (partial_sol == 1){
            break;
          }
          partial_sol--;
        } else if (vertex.pos < (n_queryv - 1)){
          partial_sol++;
        }
      }
    } while (true);
  } while (!stop);

  /* Write in output number of results */
  result = counter;
  for (int g = 0; g < STOP_S; g++) {
#pragma HLS unroll
    streams_stop[g].write(true);
  }
}

template<size_t BATCH_SIZE_LOG>
void
intersectcache_wrapper(AdjHT* hTables,
                       htb_cache_t& htb_buf,
                       hls::stream<intersect_tuple_t> stream_tuple_in[2],
                       hls::stream<split_tuple_t>& stream_tuple_out)
{
  htb_buf.init(1);
  mwj_intersect<BATCH_SIZE_LOG>(
    hTables, htb_buf, stream_tuple_in, stream_tuple_out);
#if DEBUG_INTERFACE && __SYNTHESIS__
  htb_buf.get_l1_stats(1, hits_intersect, reqs_intersect);
#endif
}

template<size_t BATCH_SIZE_LOG>
void
verifycache_wrapper(AdjHT* hTables,
                    htb_cache_t& htb_buf,
                    hls::stream<verify_tuple_t>& stream_tuple_in,
                    hls::stream<compact_tuple_t>& stream_tuple_out)
{
  htb_buf.init(0);
  mwj_verify<BATCH_SIZE_LOG>(
    hTables, htb_buf, stream_tuple_in, stream_tuple_out);
#if DEBUG_INTERFACE && __SYNTHESIS__
  htb_buf.get_l1_stats(0, hits_verify, reqs_verify);
#endif
  htb_buf.stop();
}

template<typename T_BLOOM,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t FULL_HASH_W>
void
multiwayJoin(ap_uint<DDR_W>* htb_buf0,
             ap_uint<DDR_W>* htb_buf1,
             ap_uint<DDR_W>* htb_buf2,
             ap_uint<DDR_W>* htb_buf3,
             T_BLOOM* bloom_p,
             row_t* res_buf,
             AdjHT* hTables0,
             AdjHT* hTables1,
             QueryVertex* qVertices,
             const unsigned int n_candidate,
             const unsigned int start_candidate,
             const unsigned short nQueryVer,
             const unsigned char hash1_w,
             const unsigned char hash2_w,
             const unsigned long dynfifo_space,
             unsigned int& dynfifo_overflow,
             unsigned long& reqs_dynfifo,
             long unsigned int& result)
{
#pragma HLS STABLE variable=htb_buf0
#pragma HLS STABLE variable=htb_buf1
#pragma HLS STABLE variable=htb_buf2
#pragma HLS STABLE variable=htb_buf3
#pragma HLS STABLE variable=bloom_p
#pragma HLS STABLE variable=hTables0
#pragma HLS STABLE variable=hTables1
#pragma HLS STABLE variable=qVertices

#pragma HLS DATAFLOW

    /* Propose data out */    
    hls_thread_local hls::stream<sol_node_t<vertex_t>, MAX_QV> p0_stream_sol
        ("Propose - partial solution");
    
    /* Edgebuild data out */
    hls_thread_local hls::stream<sol_node_t<vertex_t>, MAX_QV * 3> e_stream_sol
        ("Edgebuild - partial solution");
    hls_thread_local hls::stream<findmin_tuple_t, S_D> e_stream_tuple
        ("Edgebuild - tuples");
    
    /* Findmin data out */
    hls_thread_local hls::stream<bloom_t, 4>
      p_stream_filter[1UL << K_FUNCTIONS];
    hls_thread_local hls::stream<readmin_counter_tuple_t, S_D> p_stream_tuple[2];

    /* Readmin counter data out */    
    hls_thread_local hls::stream<bloom_t, 4>
      rc_stream_filter[1UL << K_FUNCTIONS];
    hls_thread_local hls::stream<readmin_edge_tuple_t, S_D> rc_stream_tuple
        ("Readmin counter - tuples");
    
    /* Readmin edge data out */    
    hls_thread_local hls::stream<sol_node_t<vertex_t>, MAX_QV> re_stream_sol
        ("Readmin edge - partial solution");
    hls_thread_local hls::stream<homomorphism_set_t<ap_uint<V_ID_W>>, S_D>
      re_stream_set[2];
    hls_thread_local hls::stream<minset_tuple_t, S_D> re_stream_tuple
        ("Readmin edge - tuples");

    /* Homomorphism data out */
    hls_thread_local hls::stream<sequencebuild_set_t<ap_uint<V_ID_W>>, S_D>
      h_stream_set("Homomorphism - set nodes");

    /* Sequencebuild data out */
    hls_thread_local hls::stream<sequencebuild_tuple_t<ap_uint<V_ID_W>>, S_D>
      sb_stream_set("Sequencebuild - set nodes");

    /* Tuplebuild data out */    
    hls_thread_local hls::stream<sol_node_t<vertex_t>, MAX_QV * 7> t_stream_sol
        ("Tuplebuild - partial solution");
    hls_thread_local hls::stream<intersect_tuple_t, S_D> t_stream_tuple[2];
    
    /* Offset data out */    
    hls_thread_local hls::stream<split_tuple_t, S_D> o_stream_tuple;
    
    /* Blockbuild data out */    
    hls_thread_local hls::stream<verify_tuple_t, S_D> bb_stream_tuple;
    
    /* Verify data out */
    hls_thread_local hls::stream<compact_tuple_t, S_D> v_stream_tuple
        ("Verify - tuples");
    
    /* Compact data out */
    hls_thread_local hls::stream<assembly_tuple_t, S_D> c_stream_tuple
        ("Compact - tuples");

    /* Filter data out */    
    hls_thread_local hls::stream<sol_node_t<vertex_t>, MAX_QV> f_stream_sol
        ("Filter - partial solution");
    hls_thread_local hls::stream<assembly_set_t, S_D> f_stream_set
        ("Filter - set nodes");
    
    /* Merge solandset data out */    
    hls_thread_local hls::stream<assembly_node_t<vertex_t>, S_D> mss_stream_sol
        ("Merge solandset - nodes");
    
    /* Assembly data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, DYN_FIFO_DEPTH> a_stream_sol
        ("Assembly - partial solution");

    /* Dynamic fifo data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, DYN_FIFO_DEPTH> dyn_stream_sol
        ("Dynamic fifo - partial solution");
 
    /* Stop signals */
    hls_thread_local hls::stream<bool, 4> streams_stop[STOP_S];

    htb_cache_t htb_cache(htb_buf0);
    htb_cache2_t htb_cache2(htb_buf2);

    dynfifo_init<ap_uint<V_ID_W>,    /* fifo data type */
                 row_t,              /* fifo data type */
                 DYN_FIFO_DEPTH,     /* in/out stream size */
                 DYN_FIFO_BURST * 2, /* load/store stream size */
                 DDR_WORD,           /* bitwidth ddr word */
                 DYN_FIFO_BURST,     /* burst transaction size */
                 RESULTS_SPACE>      /* memory words available */
      (res_buf,
       dynfifo_space,
       dynfifo_overflow,
       reqs_dynfifo,
       a_stream_sol,
       dyn_stream_sol,
       streams_stop[STOP_S - 2],
       streams_stop[STOP_S - 1]);

    hls_thread_local hls::task mwj_propose_t(
      mwj_propose, dyn_stream_sol, p0_stream_sol);

    hls_thread_local hls::task mwj_bypassfilter_t(
      mwj_bypassfilter, (hls::stream<bloom_t, 4> *)p_stream_filter, (hls::stream<bloom_t, 4>*)rc_stream_filter);

    hls_thread_local hls::task mwj_bypass_sol_t(
      mwj_bypass_sol, e_stream_sol, re_stream_sol);

    hls_thread_local hls::task mwj_sequencebuild_t(
      mwj_sequencebuild<PROPOSE_BATCH_LOG>, h_stream_set, sb_stream_set);

    hls_thread_local hls::task mwj_enlarge_sol_t(
      mwj_enlarge_sol, t_stream_sol, f_stream_sol);
    
    hls_thread_local hls::task mwj_blockbuild_t(
      mwj_blockbuild_pipelined<(1UL << PROPOSE_BATCH_LOG), __COUNTER__>, o_stream_tuple, bb_stream_tuple);

    hls_thread_local hls::task mwj_compact_t(
      mwj_compact, v_stream_tuple, c_stream_tuple);

    hls_thread_local hls::task mwj_filter_t(
      mwj_filter<(1UL << PROPOSE_BATCH_LOG)>, c_stream_tuple, f_stream_set);
    
    hls_thread_local hls::task mwj_merge_solandset_t(
      mwj_merge_solandset, f_stream_sol, f_stream_set, mss_stream_sol);

#ifdef __SYNTHESIS__

    mwj_edgebuild<LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>(
      hash1_w, qVertices, p0_stream_sol, e_stream_tuple, e_stream_sol);

    mwj_findmin<T_BLOOM,
                BLOOM_LOG,
                K_FUN_LOG,
                LKP3_HASH_W,
                MAX_HASH_W,
                FULL_HASH_W>(bloom_p,
                             hash1_w,
                             hash2_w,
                             e_stream_tuple,
                             p_stream_tuple,
                             p_stream_filter);

    mwj_readmin_counter(hTables0, htb_buf1, p_stream_tuple, rc_stream_tuple);

    cache_wrapper(mwj_readmin_edge_pipelined<T_BLOOM, BLOOM_LOG, K_FUN_LOG, FULL_HASH_W>,
                  htb_cache2,
                  rc_stream_tuple,
                  rc_stream_filter,
                  re_stream_set,
                  re_stream_tuple);

    mwj_homomorphism(
      re_stream_set, re_stream_tuple, re_stream_sol, h_stream_set);

    mwj_tuplebuild<PROPOSE_BATCH_LOG, LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>(
      hash1_w,
      hash2_w,
      qVertices,
      sb_stream_set,
      t_stream_tuple,
      t_stream_sol);

    intersectcache_wrapper<PROPOSE_BATCH_LOG>(
      hTables1, htb_cache, t_stream_tuple, o_stream_tuple);

    verifycache_wrapper<PROPOSE_BATCH_LOG>(
      hTables1, htb_cache, bb_stream_tuple, v_stream_tuple);

    mwj_assembly(htb_buf3,
                 n_candidate,
                 start_candidate,
                 nQueryVer,
                 mss_stream_sol,
                 streams_stop,
                 a_stream_sol,
                 result);

    htb_cache2.get_l1_stats(0, hits_readmin_edge, reqs_readmin_edge);
#else

    for (int g = 0; g < STOP_S + 2; g++)
      // hls::stream_globals<false>::incr_task_counter();
      hls::stream_globals::incr_task_counter();

    htb_cache.init();
    htb_cache2.init();

    std::thread mwj_edgebuild_t(mwj_edgebuild<LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>,
                                hash1_w,
                                qVertices,
                                std::ref(p0_stream_sol),
                                std::ref(e_stream_tuple),
                                std::ref(e_stream_sol));

    std::thread mwj_findmin_t(mwj_findmin<T_BLOOM, BLOOM_LOG, K_FUN_LOG, LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>,
                              bloom_p,
                              hash1_w,
                              hash2_w,
                              std::ref(e_stream_tuple),
                              p_stream_tuple,
                              p_stream_filter);

    std::thread mwj_readmin_counter_t(
        mwj_readmin_counter,
        hTables0,
        htb_buf1,
        p_stream_tuple,
        std::ref(rc_stream_tuple));

    std::thread mwj_readmin_edge_t(
      mwj_readmin_edge_pipelined<T_BLOOM, BLOOM_LOG, K_FUN_LOG, FULL_HASH_W>,
      std::ref(htb_cache2),
      std::ref(rc_stream_tuple),
      rc_stream_filter,
      re_stream_set,
      std::ref(re_stream_tuple));

    std::thread mwj_homomorphism_t(mwj_homomorphism,
                                   re_stream_set,
                                   std::ref(re_stream_tuple),
                                   std::ref(re_stream_sol),
                                   std::ref(h_stream_set));

    std::thread mwj_tuplebuild_t(
      mwj_tuplebuild<PROPOSE_BATCH_LOG, LKP3_HASH_W, MAX_HASH_W, FULL_HASH_W>,
      hash1_w,
      hash2_w,
      qVertices,
      std::ref(sb_stream_set),
      std::ref(t_stream_tuple),
      std::ref(t_stream_sol));

    std::thread mwj_intersect_t(mwj_intersect<PROPOSE_BATCH_LOG>,
                                hTables1,
                                std::ref(htb_cache),
                                std::ref(t_stream_tuple),
                                std::ref(o_stream_tuple));

    std::thread mwj_verify_t(mwj_verify<PROPOSE_BATCH_LOG>,
                             hTables1,
                             std::ref(htb_cache),
                             std::ref(bb_stream_tuple),
                             std::ref(v_stream_tuple));

    std::thread mwj_assembly_t(mwj_assembly,
                               htb_buf3,
                               n_candidate,
                               start_candidate,
                               nQueryVer,
                               std::ref(mss_stream_sol),
                               std::ref(streams_stop),
                               std::ref(a_stream_sol),
                               std::ref(result));

    mwj_edgebuild_t.join();
    mwj_findmin_t.join();
    mwj_readmin_counter_t.join();
    mwj_readmin_edge_t.join();
    mwj_homomorphism_t.join();
    mwj_tuplebuild_t.join();
    mwj_intersect_t.join();
    mwj_verify_t.join();
    mwj_assembly_t.join();

#if DEBUG_STATS
    debug::cache_hit_findmin = 0;
    debug::cache_hit_read_c = 0;
    debug::cache_hit_read_e = htb_cache2.get_n_l1_hits(0);
    debug::cache_hit_inter = htb_cache.get_n_l1_hits(1);
    debug::cache_hit_verify = htb_cache.get_n_l1_hits(0);
    debug::cache_req_findmin = reqs_findmin;
    debug::cache_req_read_c = reqs_readmin_counter;
    debug::cache_req_read_e = htb_cache2.get_n_l1_reqs(0);
    debug::cache_req_inter = htb_cache.get_n_l1_reqs(1);
    debug::cache_req_verify = htb_cache.get_n_l1_reqs(0);
#endif /* DEBUG_STATS */

    htb_cache.stop();
    htb_cache2.stop();

#endif /* __SYNTHESIS__ */
}

#if SOFTWARE_PREPROC
void subgraphIsomorphism(row_t htb_buf0[HASHTABLES_SPACE],
                         row_t htb_buf1[HASHTABLES_SPACE],
                         row_t htb_buf2[HASHTABLES_SPACE],
                         row_t htb_buf3[HASHTABLES_SPACE],
                         row_t bloom_p[BLOOM_SPACE],
                         row_t res_buf[RESULTS_SPACE],
                         const unsigned short numQueryVert,
                         const unsigned char hash1_w,
                         const unsigned char hash2_w,
                         const unsigned long dynfifo_space,
                         unsigned int &dynfifo_overflow,
                         const unsigned int n_candidate,
                         const unsigned int start_candidate,
                         QueryVertex qVertices[MAX_QV],
                         AdjHT hTables0[MAX_TB],
                         AdjHT hTables1[MAX_TB],

#if DEBUG_INTERFACE
                         volatile unsigned int &debif_endpreprocess,
                         unsigned long &p_propose_empty,
                         unsigned long &p_edgebuild_empty,
                         unsigned long &p_findmin_empty,
                         unsigned long &p_readmin_counter_empty,
                         unsigned long &p_readmin_edge_empty,
                         unsigned long &p_homomorphism_empty,
                         unsigned long &p_batchbuild_empty,
                         unsigned long &p_tuplebuild_empty,
                         unsigned long &p_intersect_empty,
                         unsigned long &p_offset_empty,
                         unsigned long &p_split_empty,
                         unsigned long &p_verify_empty,
                         unsigned long &p_compact_empty,
                         unsigned long &p_filter_empty,
                         unsigned long &p_assembly_empty,
                         unsigned long &p_hits_findmin,
                         unsigned long &p_hits_readmin_counter,
                         unsigned long &p_hits_readmin_edge,
                         unsigned long &p_hits_intersect,
                         unsigned long &p_hits_verify,
                         unsigned long &p_reqs_findmin,
                         unsigned long &p_reqs_readmin_counter,
                         unsigned long &p_reqs_readmin_edge,
                         unsigned long &p_reqs_intersect,
                         unsigned long &p_reqs_verify,
#endif /* DEBUG_INTERFACE */
                         long unsigned int &result)
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
    latency=1
#pragma HLS INTERFACE mode=m_axi port=bloom_p bundle=bloom \
    max_widen_bitwidth=128 latency=20

#pragma HLS alias ports = htb_buf0, htb_buf1, htb_buf2, htb_buf3 distance = 0

#pragma HLS INTERFACE mode = s_axilite port = numQueryVert
#pragma HLS INTERFACE mode = s_axilite port = hash1_w
#pragma HLS INTERFACE mode = s_axilite port = hash2_w
#pragma HLS INTERFACE mode = s_axilite port = dynfifo_space
#pragma HLS INTERFACE mode = s_axilite port = dynfifo_overflow
#pragma HLS INTERFACE mode = s_axilite port = return

#if DEBUG_INTERFACE
#pragma HLS INTERFACE mode = s_axilite port = debif_endpreprocess
#pragma HLS INTERFACE mode = s_axilite port = p_propose_empty
#pragma HLS INTERFACE mode = s_axilite port = p_edgebuild_empty
#pragma HLS INTERFACE mode = s_axilite port = p_findmin_empty
#pragma HLS INTERFACE mode = s_axilite port = p_readmin_counter_empty
#pragma HLS INTERFACE mode = s_axilite port = p_readmin_edge_empty
#pragma HLS INTERFACE mode = s_axilite port = p_tuplebuild_empty
#pragma HLS INTERFACE mode = s_axilite port = p_intersect_empty
#pragma HLS INTERFACE mode = s_axilite port = p_verify_empty
#pragma HLS INTERFACE mode = s_axilite port = p_assembly_empty
#pragma HLS INTERFACE mode = s_axilite port = p_homomorphism_empty
#pragma HLS INTERFACE mode = s_axilite port = p_batchbuild_empty
#pragma HLS INTERFACE mode = s_axilite port = p_offset_empty
#pragma HLS INTERFACE mode = s_axilite port = p_split_empty
#pragma HLS INTERFACE mode = s_axilite port = p_compact_empty
#pragma HLS INTERFACE mode = s_axilite port = p_filter_empty
#pragma HLS INTERFACE mode = s_axilite port = p_reqs_findmin
#pragma HLS INTERFACE mode = s_axilite port = p_reqs_readmin_counter
#pragma HLS INTERFACE mode = s_axilite port = p_reqs_readmin_edge
#pragma HLS INTERFACE mode = s_axilite port = p_reqs_intersect
#pragma HLS INTERFACE mode = s_axilite port = p_reqs_verify
#pragma HLS INTERFACE mode = s_axilite port = p_hits_findmin
#pragma HLS INTERFACE mode = s_axilite port = p_hits_readmin_counter
#pragma HLS INTERFACE mode = s_axilite port = p_hits_readmin_edge
#pragma HLS INTERFACE mode = s_axilite port = p_hits_intersect
#pragma HLS INTERFACE mode = s_axilite port = p_hits_verify
#endif /* DEBUG_INTERFACE */

#pragma HLS INTERFACE mode = s_axilite port = result

  unsigned long local_result = 0;
  unsigned int local_dynfifo_overflow = 0;

#if DEBUG_STATS
  debug::init();
#endif /* DEBUG_STATS */

#if DEBUG_INTERFACE
  ap_wait();
  debif_endpreprocess = 1;
  ap_wait();
#endif /* DEBUG_INTERFACE */

  multiwayJoin<bloom_t,
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
                   qVertices,
                   n_candidate,
                   start_candidate,
                   numQueryVert,
                   hash1_w,
                   hash2_w,
                   dynfifo_space,
                   local_dynfifo_overflow,
                   local_result);

  ap_wait();
  result = local_result;
  dynfifo_overflow = local_dynfifo_overflow;
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
  p_hits_verify = hits_verify;
  p_hits_intersect = hits_intersect;
  p_reqs_findmin = reqs_findmin;
  p_reqs_readmin_counter = reqs_readmin_counter;
  p_reqs_readmin_edge = reqs_readmin_edge;
  p_reqs_verify = reqs_verify;
  p_reqs_intersect = reqs_intersect;

#if DEBUG_STATS
  debug::print(hash1_w, hash2_w);
#endif /* DEBUG_STATS */
       // f.close();
}

#else

void subgraphIsomorphism(row_t htb_buf0[HASHTABLES_SPACE],
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
                         unsigned int &dynfifo_overflow,

#if DEBUG_INTERFACE
                         volatile unsigned int &debif_endpreprocess,
                         unsigned long &p_hits_findmin,
                         unsigned long &p_hits_readmin_counter,
                         unsigned long &p_hits_readmin_edge,
                         unsigned long &p_hits_intersect,
                         unsigned long &p_hits_verify,
                         unsigned long &p_reqs_findmin,
                         unsigned long &p_reqs_readmin_counter,
                         unsigned long &p_reqs_readmin_edge,
                         unsigned long &p_reqs_intersect,
                         unsigned long &p_reqs_verify,
                         unsigned long &p_reqs_dynfifo,
                         unsigned long &p_bloom_filtered,
#endif /* DEBUG_INTERFACE */
                         long unsigned int &result)
{

#pragma HLS INTERFACE mode=m_axi port=htb_buf0 bundle=cache \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
    latency=0
// #pragma HLS INTERFACE mode=m_axi port=htb_buf1 bundle=readmin_c \
//     max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
//     latency=0
#pragma HLS INTERFACE mode=m_axi port=htb_buf2 bundle=readmin_e \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
    latency=0 max_read_burst_length=16
// #pragma HLS INTERFACE mode=m_axi port=htb_buf3 bundle=prop_batch \
//     max_widen_bitwidth=128 latency=0
// #pragma HLS INTERFACE mode=m_axi port=htb_buf2 bundle=readmin_e \
//     max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
//     latency=0 max_read_burst_length=16 channel=0
#pragma HLS INTERFACE mode=m_axi port=htb_buf1 bundle=readmin_c \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
    latency=0 max_read_burst_length=16 channel=1
#pragma HLS INTERFACE mode=m_axi port=htb_buf3 bundle=readmin_c \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
    latency=0 max_read_burst_length=16 channel=2
// #pragma HLS INTERFACE mode=m_axi port=htb_buf0 bundle=readmin_e \
//     max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2 \
//     latency=0 max_read_burst_length=16 channel=3
#pragma HLS INTERFACE mode=m_axi port=res_buf bundle=fifo \
    max_widen_bitwidth=128 max_read_burst_length=32 max_write_burst_length=32 \
    latency=0
#pragma HLS INTERFACE mode=m_axi port=bloom_p bundle=bloom \
    max_widen_bitwidth=128 latency=20

// #pragma HLS alias ports=htb_buf0,htb_buf1,htb_buf2,htb_buf3 distance=0

#pragma HLS INTERFACE mode=s_axilite port=numQueryVert
#pragma HLS INTERFACE mode=s_axilite port=numQueryEdges
#pragma HLS INTERFACE mode=s_axilite port=numDataEdges
#pragma HLS INTERFACE mode=s_axilite port=hash1_w
#pragma HLS INTERFACE mode=s_axilite port=hash2_w
#pragma HLS INTERFACE mode=s_axilite port=dynfifo_space
#pragma HLS INTERFACE mode=s_axilite port=dynfifo_overflow
#pragma HLS INTERFACE mode=s_axilite port=return

// #pragma HLS cache port=htb_buf2 lines=512 depth=32

#if DEBUG_INTERFACE
#pragma HLS INTERFACE mode=s_axilite port=debif_endpreprocess
#pragma HLS INTERFACE mode=s_axilite port=p_reqs_findmin
#pragma HLS INTERFACE mode=s_axilite port=p_reqs_readmin_counter
#pragma HLS INTERFACE mode=s_axilite port=p_reqs_readmin_edge
#pragma HLS INTERFACE mode=s_axilite port=p_reqs_intersect
#pragma HLS INTERFACE mode=s_axilite port=p_reqs_verify
#pragma HLS INTERFACE mode=s_axilite port=p_reqs_dynfifo
#pragma HLS INTERFACE mode=s_axilite port=p_hits_findmin
#pragma HLS INTERFACE mode=s_axilite port=p_hits_readmin_counter
#pragma HLS INTERFACE mode=s_axilite port=p_hits_readmin_edge
#pragma HLS INTERFACE mode=s_axilite port=p_hits_intersect
#pragma HLS INTERFACE mode=s_axilite port=p_hits_verify
#pragma HLS INTERFACE mode=s_axilite port=p_bloom_filtered
#endif /* DEBUG_INTERFACE */

#pragma HLS INTERFACE mode=s_axilite port=result

#if DEBUG_STATS
    debug::init();
#endif /* DEBUG_STATS */

    QueryVertex qVertices[MAX_QV];
    AdjHT hTables0[MAX_TB], hTables1[MAX_TB];
    unsigned long localResult = 0;
    unsigned int local_dynfifo_overflow = 0;
    unsigned int n_candidate = 0;
    unsigned int start_candidate = 0;
    unsigned long reqs_dynfifo = 0;

    preprocess<row_t,
               bloom_t,
               EDGE_WIDTH,
               COUNTER_WIDTH,
               BLOOM_FILTER_WIDTH,
               K_FUNCTIONS,
               DDR_BIT,
               VERTEX_WIDTH_BIT,
               VERTEX_WIDTH,
               HASH_LOOKUP3_BIT,
               MAX_HASH_TABLE_BIT,
               64,
               LABEL_WIDTH,
               DEFAULT_STREAM_DEPTH,
               HASHTABLES_SPACE,
               MAX_QUERY_VERTICES,
               MAX_TABLES,
               MAX_COLLISIONS>(res_buf,
                               htb_buf0,
                               htb_buf1,
                               bloom_p,
                               qVertices,
                               hTables0,
                               hTables1,
                               dynfifo_space,
                               n_candidate,
                               start_candidate,
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

    multiwayJoin<bloom_t,
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
                     qVertices,
                     n_candidate,
                     start_candidate,
                     numQueryVert,
                     hash1_w,
                     hash2_w,
                     dynfifo_space,
                     local_dynfifo_overflow,
                     reqs_dynfifo,
                     localResult);
    
    ap_wait();
    result = localResult;
    dynfifo_overflow = local_dynfifo_overflow;
    p_reqs_findmin = reqs_findmin;
    p_reqs_readmin_counter = reqs_readmin_counter;
    p_reqs_readmin_edge = reqs_readmin_edge;
    p_reqs_intersect = reqs_intersect;
    p_reqs_verify = reqs_verify;
    p_reqs_dynfifo = reqs_dynfifo;
    p_hits_findmin = hits_findmin;
    p_hits_readmin_counter = hits_readmin_counter;
    p_hits_readmin_edge = hits_readmin_edge;
    p_hits_verify = hits_verify;
    p_hits_intersect = hits_intersect;
    p_bloom_filtered = bloom_filtered;

#if DEBUG_STATS
    debug::print(hash1_w, hash2_w);
#endif /* DEBUG_STATS */
}
#endif /* SOFTWARE_PREPROC */

#pragma GCC diagnostic pop
