#pragma once

#define HLS_STREAM_THREAD_SAFE
#ifndef __SYNTHESIS__
#include <cassert>
#include <fstream>
#include <limits.h>
#endif

#include <hls_stream.h>
#include <ap_int.h>

#include "Parameters.hpp"
#include "QueryVertex.hpp"
#include "Trie.hpp"
#include "hash_lookup3.hpp"

#if DEBUG_STATS
#include "debug.hpp"
#endif /* DEBUG_STATS */

#pragma GCC diagnostic push
// #pragma GCC diagnostic error "-Wpedantic"
#pragma GCC diagnostic error "-Wall"
// #pragma GCC diagnostic error "-Wextra"
#pragma GCC diagnostic ignored "-Wunused-label"
#pragma GCC diagnostic ignored "-Wsign-compare"
// #pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"

struct bloom_write_tuple_t
{
  unsigned int address;
  bool last;
};

template<size_t HASH_W>
struct bloom_update_tuple_t
{
  unsigned int address;
  ap_uint<HASH_W> indexed_h;
  bool write;
  bool last;
};

template<size_t NODE_W>
struct bagtoset_tuple_t
{
  ap_uint<NODE_W> indexing_v;
  bool write;
  bool last;
  bool valid;
};

template<size_t NODE_W>
struct batch_tuple_t
{
  ap_uint<NODE_W> indexing_v;
  bool last;
};

struct counter_tuple_t
{
  unsigned int address;
  bool stop;
};

template<typename ROW_T>
struct store_tuple_t
{
  unsigned int address;
  ROW_T edge;
  bool stop;
};

/* Builds the table descriptors based on the information
 * from the query graph. */
template<size_t MAX_QV,
         size_t MAX_TB,
         size_t NODE_W,
         size_t LAB_W,
         size_t MAX_LABELS>
void
buildTableDescriptors(row_t* edge_buf,
                      QueryVertex* qVertices,
                      ap_uint<8> labelToTable[MAX_LABELS][MAX_LABELS],
                      unsigned short& numTables,
                      unsigned short numQueryVert,
                      unsigned short numQueryEdge)
{
    /* Translate from id of vertex to position in the order */
    unsigned short fromNumToPos[MAX_QV];
    constexpr size_t SRC_NODE = 0;
    constexpr size_t DST_NODE = 32;
    constexpr size_t LABELSRC_NODE = 64;
    constexpr size_t LABELDST_NODE = 96;

    /* Filling information about query vertices and coping
     * the vertex order needed by multiway join */
FILL_ORDER_LOOP:
    for (int g = 0; g < numQueryVert; g++) {
#pragma HLS pipeline II = 1
      ap_uint<NODE_W> nodesrc =
        edge_buf[g].range(SRC_NODE + NODE_W - 1, SRC_NODE);
#ifndef __SYNTHESIS__
      assert(numQueryVert < (MAX_QV));
#endif
      fromNumToPos[nodesrc] = g;
    }

    /* Creating table descriptors */
CREATE_TABDESC_LOOP:
    for (int s = 0; s < numQueryEdge; s++) {
#pragma HLS pipeline II = 4
      bool dirEdge = false;
      ap_uint<8> index = 0;
      row_t edge = edge_buf[s + numQueryVert];
      ap_uint<LAB_W> labeldst =
        edge.range(LABELDST_NODE + LAB_W - 1, LABELDST_NODE);
      ap_uint<LAB_W> labelsrc =
        edge.range(LABELSRC_NODE + LAB_W - 1, LABELSRC_NODE);
      ap_uint<NODE_W> nodedst = edge.range(DST_NODE + NODE_W - 1, DST_NODE);
      ap_uint<NODE_W> nodesrc = edge.range(SRC_NODE + NODE_W - 1, SRC_NODE);
      unsigned short nodeSrcPos = fromNumToPos[nodesrc];
      unsigned short nodeDstPos = fromNumToPos[nodedst];

      // Direction of the table is used to understand if the
      // source vertex is indexed or indexing the table
      if (nodeSrcPos < nodeDstPos)
        dirEdge = true;

#ifndef __SYNTHESIS__
      std::cout << (unsigned int)nodesrc << "(" << (int)labelsrc << ")"
                << " -> " << (unsigned int)nodedst << "(" << (int)labeldst
                << ")" << std::endl;
#endif

      // Saving the index of the table in the labels matrix
      // which is indexed by [indexing label][indexed label]
      if (dirEdge) {
        index = labelToTable[labelsrc][labeldst];
        if (index == 0) {
          index = ++numTables;
        }
        labelToTable[labelsrc][labeldst] = index;
      } else {
        index = labelToTable[labeldst][labelsrc];
        if (index == 0) {
          index = ++numTables;
        }
        labelToTable[labeldst][labelsrc] = index;
      }

#ifndef __SYNTHESIS__
      if (dirEdge) {
        std::cout << "Table " << (int)index - 1 << ": " << (int)labelsrc
                  << " -> " << (int)labeldst << std::endl;
      } else {
        std::cout << "Table " << (int)index - 1 << ": " << (int)labeldst
                  << " <- " << (int)labelsrc << std::endl;
      }
#endif

      /* Linking vertices to tables */
      if (dirEdge) {
        unsigned char idx = qVertices[nodeSrcPos].numTablesIndexing;
        qVertices[nodeSrcPos].tables_indexing[idx] = index - 1;
        qVertices[nodeSrcPos].numTablesIndexing++;

        idx = qVertices[nodeDstPos].numTablesIndexed;
        qVertices[nodeDstPos].tables_indexed[idx] = index - 1;
        qVertices[nodeDstPos].vertex_indexing[idx] = nodeSrcPos;

        qVertices[nodeDstPos].numTablesIndexed++;
      } else {
        unsigned char idx = qVertices[nodeDstPos].numTablesIndexing;
        qVertices[nodeDstPos].tables_indexing[idx] = index - 1;
        qVertices[nodeDstPos].numTablesIndexing++;

        idx = qVertices[nodeSrcPos].numTablesIndexed;
        qVertices[nodeSrcPos].tables_indexed[idx] = index - 1;
        qVertices[nodeSrcPos].vertex_indexing[idx] = nodeDstPos;

        qVertices[nodeSrcPos].numTablesIndexed++;
      }
    }
}

/* Reads edges from each table and divide the indexed vertices based
on indexing hash to create bloom filters */
template <typename T_DDR,
          size_t CNT_LOG,
          size_t ROW_LOG,
          size_t NODE_W,
          size_t EDGE_LOG,
          size_t LKP3_HASH_W,
          size_t MAX_HASH_W,
          size_t FULL_HASH_W>
void bloomRead(AdjHT *hTables,
               QueryVertex *qVertices,
               T_DDR *htb_buf,
               const unsigned short numTables,
               const unsigned char hash1_w,
               hls::stream<bagtoset_tuple_t<NODE_W> > &stream_tuple_bagtoset_out,
               hls::stream<bloom_update_tuple_t<FULL_HASH_W>> &stream_tuple_bloom_out)
{
  constexpr size_t EDGE_W = 1UL << EDGE_LOG;
  hls::stream<ap_uint<NODE_W>, 4> hash_in0, hash_in1;
  hls::stream<ap_uint<FULL_HASH_W>, 4> hash_out0;
  hls::stream<ap_uint<LKP3_HASH_W>, 4> hash_out1;
  ap_uint<EDGE_W> edge;
  ap_uint<NODE_W> indexing_v, indexed_v, prev_indexing_v;
  ap_uint<FULL_HASH_W> indexed_h, prev_indexed_h;
  ap_uint<MAX_HASH_W> indexing_h, prev_indexing_h;
  bloom_update_tuple_t<FULL_HASH_W> tuple_out;
  T_DDR row;
  unsigned int counter;
  prev_indexing_h = 0;
  prev_indexed_h = 0;

// Select the table with the minimum number of edges
// to start the partial solutions.
  unsigned int minSize = UINT32_MAX;
  unsigned short minTableIndex;
PROPOSE_TBINDEXING_LOOP:
  for (int g = 0; g < qVertices[0].numTablesIndexing; g++) {
    unsigned short tableIndex = qVertices[0].tables_indexing[g];

    if (hTables[tableIndex].n_edges < minSize) {
      minSize = hTables[tableIndex].n_edges;
      minTableIndex = tableIndex;
    }
  }

BLOOM_READ_TASK_LOOP:
  for (unsigned int ntb = 0; ntb < numTables; ntb++)
  {

    /* During first iteration do not consider the difference between
    prev_indexing_h and indexing_h to be useful to write the bloom */
    bool first_it = true;
    counter = 0;
    unsigned int cycles = (hTables[ntb].n_edges >> 1) + 1;
    unsigned int offset = hTables[ntb].start_edges;

    /* Read all the edges in a table and divide them by hash1 */
  BLOOM_READ_EDGES_BLOCK:
    for (unsigned int start = 0; start < cycles; start++)
    {
#pragma HLS pipeline II = 2

      row = htb_buf[offset + start];

      for (int i = 0; i < EDGE_ROW; i++)
      {
#pragma HLS unroll
        edge = row.range(((i + 1) << EDGE_LOG) - 1, i << EDGE_LOG);
        indexing_v = edge.range(NODE_W * 2 - 1, NODE_W);
        indexed_v = edge.range(NODE_W - 1, 0);

        hash_in0.write(indexed_v);
        hash_in1.write(indexing_v);
        xf::database::hashLookup3<NODE_W>(hash_in0, hash_out0);
        xf::database::hashLookup3<NODE_W>(hash_in1, hash_out1);
        indexed_h = hash_out0.read();
        indexing_h = hash_out1.read();
        indexing_h = indexing_h.range(hash1_w - 1, 0);

        bool valid = (counter < hTables[ntb].n_edges);
        bool write = (indexing_h != prev_indexing_h);
        /* Writing edge of previous iteration */
        if (valid && !first_it)
        {
          tuple_out.address = ntb * (1UL << hash1_w) + prev_indexing_h;
          tuple_out.last = false;
          tuple_out.write = write;
          tuple_out.indexed_h = prev_indexed_h;
          stream_tuple_bloom_out.write(tuple_out);

          if (ntb == minTableIndex)
          {
            bagtoset_tuple_t<NODE_W> tuple_bagtoset_out;
            tuple_bagtoset_out.indexing_v = prev_indexing_v;
            tuple_bagtoset_out.write = write;
            tuple_bagtoset_out.last = false;
            tuple_bagtoset_out.valid = true;
            stream_tuple_bagtoset_out.write(tuple_bagtoset_out);
          }
        }

        if (valid)
        {
          prev_indexing_h = indexing_h;
          prev_indexing_v = indexing_v;
          prev_indexed_h = indexed_h;
        }
        counter++;
        first_it = false;
      }
    }

    /* Write explicitly the last bloom filter since
    the difference between prev_indexing_h and indexing_h
    does not work at the end of the table */
    tuple_out.address = ntb * (1UL << hash1_w) + prev_indexing_h;
    tuple_out.indexed_h = prev_indexed_h;
    tuple_out.write = true;
    tuple_out.last = (ntb == (numTables - 1));
    stream_tuple_bloom_out.write(tuple_out);

    bagtoset_tuple_t<NODE_W> tuple_bagtoset_out;
    tuple_bagtoset_out.indexing_v = prev_indexing_v;
    tuple_bagtoset_out.write = true;
    tuple_bagtoset_out.valid = ntb == minTableIndex;
    tuple_bagtoset_out.last = (ntb == (numTables - 1));
    stream_tuple_bagtoset_out.write(tuple_bagtoset_out);
  }
}

template <size_t NODE_W,
          size_t MAX_CL>
void bagtoset(hls::stream<bagtoset_tuple_t<NODE_W> > &stream_tuple_in,
              hls::stream<batch_tuple_t<NODE_W> > &stream_tuple_out)
{
  bagtoset_tuple_t<NODE_W> tuple_in;
  ap_uint<NODE_W> set[MAX_CL];
  unsigned char pointer = 0;
  ap_uint<MAX_CL> valid_bits = 0; // 1 if the element is present
  ap_uint<MAX_CL> equal_bits = 0; // 1 if the element is equal to the one in the bag

BAGTOSET_TASK_LOOP:
  do
  {
#pragma HLS pipeline II = 2
    tuple_in = stream_tuple_in.read();

    if (tuple_in.valid)
    {
      // Check if the element is present in the set
      for (int i = 0; i < MAX_CL; i++)
      {
#pragma HLS unroll
        equal_bits[i] = (set[i] == tuple_in.indexing_v);
      }

      // If the element is not present in the set, add it
      // and write it to the output stream
      if ((equal_bits & valid_bits) == 0)
      {
        set[pointer] = tuple_in.indexing_v;
        valid_bits = valid_bits | (1 << pointer);
        pointer++;
        batch_tuple_t<NODE_W> tuple_out;
        tuple_out.indexing_v = tuple_in.indexing_v;
        tuple_out.last = false;
        stream_tuple_out.write(tuple_out);
      }

      // Last element of the hash set
      if (tuple_in.write)
      {
        pointer = 0;
        valid_bits = 0;
      }
    }

  } while (!tuple_in.last);

  batch_tuple_t<NODE_W> tuple_out;
  tuple_out.indexing_v = tuple_in.indexing_v;
  tuple_out.last = true;
  stream_tuple_out.write(tuple_out);
}

template <size_t NODE_LOG,
          size_t NODE_W,
          size_t ROW_LOG>
void batch(unsigned int &n_candidate,
               const unsigned int start_address,
               row_t *htb_buf,
               hls::stream<batch_tuple_t<NODE_W>> &stream_tuple_in)
{
  constexpr size_t NODE_PER_WORD_LOG = ROW_LOG - NODE_LOG;
  row_t word;
  ap_uint<32> pointer = 0;
  unsigned int offset = 0;

  batch_tuple_t<NODE_W> tuple_in = stream_tuple_in.read();
  while (!tuple_in.last)
  {
#pragma HLS pipeline II = 1
    ap_uint<NODE_PER_WORD_LOG> in_word_pointer = pointer.range(NODE_PER_WORD_LOG - 1, 0);
    word.range(NODE_W * (in_word_pointer + 1) - 1, NODE_W * in_word_pointer) = tuple_in.indexing_v;
    if (in_word_pointer == (1UL << NODE_PER_WORD_LOG) - 1)
    {
      htb_buf[start_address + offset] = word;
      offset++;
    }
    pointer++;
    tuple_in = stream_tuple_in.read();
  };
  htb_buf[start_address + offset] = word;
  n_candidate = pointer;
}

template<typename T_BLOOM,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t FULL_HASH_W>
void
bloomUpdate(hls::stream<bloom_update_tuple_t<FULL_HASH_W> >& stream_tuple_in,
           hls::stream<bloom_write_tuple_t>& stream_address,
           hls::stream<T_BLOOM> stream_filter[(1UL << K_FUN_LOG)])
{
    constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    T_BLOOM filter[K_FUN];
    bloom_update_tuple_t<FULL_HASH_W> tuple_in;
#pragma HLS array_partition variable = filter type = complete

    /*Reset filter*/
    for (auto g = 0; g < K_FUN; g++) {
#pragma HLS unroll
        filter[g] = 0;
    }

    do {
#pragma HLS pipeline II = 1

        tuple_in = stream_tuple_in.read();
        for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
            ap_uint<BLOOM_LOG> idx = tuple_in.indexed_h.range(
              (FULL_HASH_W / K_FUN) * (g + 1) - 1,
              (FULL_HASH_W / K_FUN) * (g + 1) - BLOOM_LOG);
            filter[g][idx] = 1;
            if (tuple_in.write) {
                stream_filter[g].write(filter[g]);
                filter[g] = 0;
            }
        }

        if (tuple_in.write) {
            stream_address.write({ tuple_in.address, tuple_in.last });
        }
    } while (!tuple_in.last);
}

template <typename T_BLOOM, size_t K_FUN_LOG>
void bloomWrite(T_BLOOM *bloom_p,
                hls::stream<bloom_write_tuple_t> &stream_address,
                hls::stream<T_BLOOM> stream_filter[(1UL << K_FUN_LOG)])
{
  constexpr size_t K_FUN = (1UL << K_FUN_LOG);
  bloom_write_tuple_t tuple_in;
BLOOM_WRITE_TASK_LOOP:
  do
  {
#pragma HLS pipeline II = (1UL << K_FUN_LOG)
    tuple_in = stream_address.read();
    for (int g = 0; g < K_FUN; g++)
    {
#pragma HLS unroll
      bloom_p[(tuple_in.address << K_FUN_LOG) + g] =
          stream_filter[g].read();

#if DEBUG_STATS
      /* Computing the number of ones in each filter*/
      T_BLOOM row = bloom_p[(tuple_in.address << K_FUN_LOG) + g];
      while (row > 0)
      {
        debug::bloom_fullness++;
        row = row & (row - 1);
      }
#endif /* DEBUG_STATS */
    }
  } while (!tuple_in.last);
}

template <typename T_DDR,
          typename T_BLOOM,
          size_t CNT_LOG,
          size_t ROW_LOG,
          size_t NODE_W,
          size_t NODE_LOG,
          size_t EDGE_LOG,
          size_t MAX_CL,
          size_t LKP3_HASH_W,
          size_t MAX_HASH_W,
          size_t FULL_HASH_W,
          size_t BLOOM_LOG,
          size_t K_FUN_LOG,
          size_t STREAM_D>
void writeBloom(
    T_BLOOM *bloom_p,
    T_DDR *htb_p0,
    T_DDR *htb_p1,
    AdjHT *hTables,
    QueryVertex *qVertices,
    unsigned int &n_candidate,
    const unsigned int start_address,
    const unsigned char numTables,
    const unsigned char hash1_w)
{
#pragma HLS DATAFLOW

    hls::stream<bloom_update_tuple_t<FULL_HASH_W>, 32>
      stream_tuple("Bloom tuple");
    hls::stream<bloom_write_tuple_t, 8> stream_address(
      "Bloom address");
    hls::stream<T_BLOOM, 8> stream_filter[(1UL << K_FUN_LOG)];
    hls::stream<bagtoset_tuple_t<NODE_W>, 8> stream_tuple_bagtoset(
      "Bagtoset tuple");
    hls::stream<batch_tuple_t<NODE_W>, 8> stream_tuple_batch(
      "Batch tuple");

    /* Read edges in each table */
    bloomRead<T_DDR,
              CNT_LOG,
              ROW_LOG,
              NODE_W,
              EDGE_LOG,
              LKP3_HASH_W,
              MAX_HASH_W,
              FULL_HASH_W>(hTables,
                           qVertices,
                           htb_p0,
                           numTables,
                           hash1_w,
                           stream_tuple_bagtoset,
                           stream_tuple);

    bloomUpdate<T_BLOOM, BLOOM_LOG, K_FUN_LOG, FULL_HASH_W>(
      stream_tuple, stream_address, stream_filter);

    bagtoset<NODE_W, MAX_CL>(stream_tuple_bagtoset, stream_tuple_batch);

    bloomWrite<T_BLOOM, K_FUN_LOG>(bloom_p, stream_address, stream_filter);

    batch<NODE_LOG, NODE_W, ROW_LOG>(n_candidate, start_address, htb_p1, stream_tuple_batch);
}

template<size_t NODE_W,
         size_t LAB_W,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t MAX_LABELS>
void
readEdgesPerBlock(row_t* edge_buf,
                  const unsigned char hash1_w,
                  const unsigned char hash2_w,
                  const ap_uint<8> labelToTable[MAX_LABELS][MAX_LABELS],
                  const unsigned long numDataEdges,
                  hls::stream<counter_tuple_t> stream_address[2])
{
    constexpr size_t COUNTERS_PER_BLOCK = 14;
    constexpr size_t SRC_NODE = 0;
    constexpr size_t DST_NODE = 32;
    constexpr size_t LABELSRC_NODE = 64;
    constexpr size_t LABELDST_NODE = 96;
    const unsigned int block_per_table = hash1_w + hash2_w - COUNTERS_PER_BLOCK;
    bool inverted = false;

READ_EDGES_PER_BLOCK_LOOP:
    for (auto s = 0; s < numDataEdges; s++) {
#pragma HLS pipeline II = 1

        row_t edge = edge_buf[s];

        ap_uint<LAB_W> labeldst =
          edge.range(LABELDST_NODE + LAB_W - 1, LABELDST_NODE);
        ap_uint<LAB_W> labelsrc =
          edge.range(LABELSRC_NODE + LAB_W - 1, LABELSRC_NODE);
        ap_uint<NODE_W> nodedst = edge.range(DST_NODE + NODE_W - 1, DST_NODE);
        ap_uint<NODE_W> nodesrc = edge.range(SRC_NODE + NODE_W - 1, SRC_NODE);
        
        // Retrieve index of table with source as indexing vertex
        ap_uint<8> index0 = labelToTable[labelsrc][labeldst];
        // Retrieve index of table with destination as indexing vertex
        ap_uint<8> index1 = labelToTable[labeldst][labelsrc];

        /* Compute indices for hash table */
        ap_uint<LKP3_HASH_W> hash_out0;
        ap_uint<LKP3_HASH_W> hash_out1;
        xf::database::details::hashlookup3_core<NODE_W>(nodesrc, hash_out0);

        ap_uint<MAX_HASH_W> hashsrc = hash_out0.range(MAX_HASH_W - 1, 0);
        hashsrc = hashsrc.range(hash1_w - 1, 0);

        xf::database::details::hashlookup3_core<NODE_W>(nodedst, hash_out1);

        ap_uint<MAX_HASH_W> hashdst = hash_out1.range(MAX_HASH_W - 1, 0);
        hashdst = hashdst.range(hash1_w - 1, 0);

        /* Compute inside which block the edge will finish, in the meanwhile
        also adjust the edge as indexing -> indexed, and also add the hash
        values */
        unsigned short address_intable0 =
          hashsrc.range(hash1_w - 1, hash1_w - block_per_table);
        unsigned int address0 = ((index0 - 1) << block_per_table) + address_intable0;
        unsigned short address_intable1 =
          hashdst.range(hash1_w - 1, hash1_w - block_per_table);
        unsigned int address1 = ((index1 - 1) << block_per_table) + address_intable1;

        /* This useless if is to explain to Vitis HLS 2022.2 that two write in
         * the same stream cannot happen in one cycle */
        if (index0 != 0 && index1 != 0) {
          if (inverted){
            stream_address[1].write({ address0, false });
            stream_address[0].write({ address1, false });
          } else {
            stream_address[0].write({ address0, false });
            stream_address[1].write({ address1, false });
          }
        } else if (index0 != 0){
          if (inverted){
            stream_address[1].write({ address0, false });
          } else {
            stream_address[0].write({ address0, false });
          }
        } else if (index1 != 0){
          if (inverted){
            stream_address[1].write({ address1, false });
          } else {
            stream_address[0].write({ address1, false });
          }
        }

        if ((index0 != 0) ^ (index1 != 0)){
          inverted = !inverted;
        }
    }
    if (inverted){
      stream_address[1].write({ 0, true });
    } else {
      stream_address[0].write({ 0, true });
    }
}

template<size_t NODE_W,
         size_t LAB_W,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t MAX_LABELS>
void
countEdgesPerBlock(hls::stream<counter_tuple_t> stream_address[2],
                  unsigned int block_n_edges[4096])
{
    const size_t BRAM_LAT = 3;
    unsigned int local_cache_address[BRAM_LAT];
    unsigned int local_cache_counter[BRAM_LAT];
    ap_uint<BRAM_LAT> local_cache_valid = 0;
    counter_tuple_t tuple_in;
    auto select = 0;
    
    tuple_in = stream_address[select].read();
    select = (select + 1) % 2;
COUNT_EDGES_PER_BLOCK_LOOP:
    while (!tuple_in.stop) {
#pragma HLS dependence variable = block_n_edges type = inter direction =       \
  RAW false
#pragma HLS pipeline II = 1

        unsigned int address = tuple_in.address;

        bool hit = false;
        unsigned int local_value_counter = 0;

        /* Check if the counter has been used recently, by cycling backword to
         * catch the updated value */
        for (auto s = 0; s < BRAM_LAT; s++) {
#pragma HLS unroll
            auto g = BRAM_LAT - s - 1;
            if (local_cache_address[g] == address && local_cache_valid[g]) {
                hit = true;
                local_value_counter = local_cache_counter[g];
            }
        }

        /* Read from memory only if is not present in local cache, in this way
         * it is possible to remove the RAW dependency */
        if (hit){
          local_value_counter++;
        } else {
          local_value_counter = block_n_edges[address];
          local_value_counter++;
        }

        /* Shift everything by one position and writes the last one in memory */
        for (auto s = 0; s < BRAM_LAT - 1; s++) {
#pragma HLS unroll
          auto g = BRAM_LAT - s - 1;
          local_cache_address[g] = local_cache_address[g - 1];
          local_cache_counter[g] = local_cache_counter[g - 1];
          local_cache_valid[g] = local_cache_valid[g - 1];
        }
        local_cache_address[0] = address;
        local_cache_counter[0] = local_value_counter;
        local_cache_valid[0] = true;
        block_n_edges[address] = local_value_counter;

        tuple_in = stream_address[select].read();
        select = (select + 1) % 2;

#ifndef __SYNTHESIS__
        assert(local_value_counter < UINT32_MAX);
#endif
    }
}


template<size_t NODE_W,
         size_t LAB_W,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t MAX_LABELS>
void
countEdgesPerBlockWrap(row_t* edge_buf,
                const unsigned char hash1_w,
                const unsigned char hash2_w,
                const ap_uint<8> labelToTable[MAX_LABELS][MAX_LABELS],
                const unsigned long numDataEdges,
                unsigned int block_n_edges[4096])
{
#pragma HLS dataflow
    hls::stream<counter_tuple_t, 32> stream_address[2];

    readEdgesPerBlock<NODE_W, LAB_W, LKP3_HASH_W, MAX_HASH_W, MAX_LABELS>(
      edge_buf, hash1_w, hash2_w, labelToTable, numDataEdges, stream_address);

    countEdgesPerBlock<NODE_W, LAB_W, LKP3_HASH_W, MAX_HASH_W, MAX_LABELS>(
      stream_address, block_n_edges);
}

template<size_t NODE_W,
         size_t LAB_W,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t MAX_LABELS>
void
readAndStreamEdgesPerBlock(row_t* edge_buf,
                  const unsigned char hash1_w,
                  const unsigned char hash2_w,
                  const ap_uint<8> labelToTable[MAX_LABELS][MAX_LABELS],
                  const unsigned int numDataEdges,
                  hls::stream<store_tuple_t<row_t> > stream_edge[2])
{
    constexpr size_t COUNTERS_PER_BLOCK = 14;
    constexpr size_t SRC_NODE = 0;
    constexpr size_t DST_NODE = 32;
    constexpr size_t LABELSRC_NODE = 64;
    constexpr size_t LABELDST_NODE = 96;
    const unsigned int block_per_table = hash1_w + hash2_w - COUNTERS_PER_BLOCK;
    bool inverted = false;

STORE_EDGE_PER_BLOCK_LOOP:
    for (auto s = 0; s < numDataEdges; s++) {
#pragma HLS pipeline II = 1
        row_t edge = edge_buf[s];

        ap_uint<LAB_W> labeldst =
          edge.range(LABELDST_NODE + LAB_W - 1, LABELDST_NODE);
        ap_uint<LAB_W> labelsrc =
          edge.range(LABELSRC_NODE + LAB_W - 1, LABELSRC_NODE);
        ap_uint<NODE_W> nodedst = edge.range(DST_NODE + NODE_W - 1, DST_NODE);
        ap_uint<NODE_W> nodesrc = edge.range(SRC_NODE + NODE_W - 1, SRC_NODE);

        // Retrieve index of table with source as indexing vertex
        ap_uint<8> index0 = labelToTable[labelsrc][labeldst];
        // Retrieve index of table with destination as indexing vertex
        ap_uint<8> index1 = labelToTable[labeldst][labelsrc];

        /* Compute indices for hash table */
        ap_uint<LKP3_HASH_W> hash_out0;
        ap_uint<LKP3_HASH_W> hash_out1;
        xf::database::details::hashlookup3_core<NODE_W>(nodesrc, hash_out0);

        ap_uint<MAX_HASH_W> hashsrc = hash_out0.range(MAX_HASH_W - 1, 0);
        hashsrc = hashsrc.range(hash1_w - 1, 0);

        xf::database::details::hashlookup3_core<NODE_W>(nodedst, hash_out1);

        ap_uint<MAX_HASH_W> hashdst = hash_out1.range(MAX_HASH_W - 1, 0);
        hashdst = hashdst.range(hash1_w - 1, 0);

        /* Compute inside which block the edge will finish, in the meanwhile
        also adjust the edge as indexing -> indexed, and also add the hash
        values */
        unsigned short address_intable0 =
          hashsrc.range(hash1_w - 1, hash1_w - block_per_table);
        unsigned int address0 =
          ((index0 - 1) << block_per_table) + address_intable0;
        unsigned short address_intable1 =
          hashdst.range(hash1_w - 1, hash1_w - block_per_table);
        unsigned int address1 =
          ((index1 - 1) << block_per_table) + address_intable1;

        row_t table_edge0;
        table_edge0.range(31, 0) = nodesrc;
        table_edge0.range(63, 32) = nodedst;
        table_edge0.range(95, 64) = hashsrc;
        table_edge0.range(127, 96) = hashdst.range(hash2_w - 1, 0);

        row_t table_edge1;
        table_edge1.range(31, 0) = nodedst;
        table_edge1.range(63, 32) = nodesrc;
        table_edge1.range(95, 64) = hashdst;
        table_edge1.range(127, 96) = hashsrc.range(hash2_w - 1, 0);

        /* This useless if is to explain to Vitis HLS 2022.2 that two write in
         * the same stream cannot happen in one cycle */
        if (index0 != 0 && index1 != 0) {
          if (inverted) {
                stream_edge[1].write({ address0, table_edge0, false });
                stream_edge[0].write({ address1, table_edge1, false });
          } else {
                stream_edge[0].write({ address0, table_edge0, false });
                stream_edge[1].write({ address1, table_edge1, false });
          }
        } else if (index0 != 0) {
          if (inverted) {
                stream_edge[1].write({ address0, table_edge0, false });
          } else {
                stream_edge[0].write({ address0, table_edge0, false });
          }
        } else if (index1 != 0) {
          if (inverted) {
                stream_edge[1].write({ address1, table_edge1, false });
          } else {
                stream_edge[0].write({ address1, table_edge1, false });
          }
        }

        if ((index0 != 0) ^ (index1 != 0)) {
          inverted = !inverted;
        }
    }
    if (inverted) {
        stream_edge[1].write({ 0, 0, true });
    } else {
        stream_edge[0].write({ 0, 0, true });
    }
}

template<size_t NODE_W,
         size_t LAB_W,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t MAX_LABELS>
void
storeEdgesPerBlock(hls::stream<store_tuple_t<row_t> > stream_edge[2],
                  row_t* m_axi,
                  unsigned int block_n_edges[4096])
{
    const size_t BRAM_LAT = 3;
    unsigned int local_cache_address[BRAM_LAT];
    unsigned int local_cache_counter[BRAM_LAT];
    ap_uint<BRAM_LAT> local_cache_valid = 0;
    store_tuple_t<row_t> tuple_in;
    auto select = 0;
    
    tuple_in = stream_edge[select].read();
    select = (select + 1) % 2;
STORE_EDGES_PER_BLOCK_LOOP:
    while (!tuple_in.stop) {
#pragma HLS dependence variable = block_n_edges type = inter direction =       \
  RAW false
#pragma HLS pipeline II = 1

        unsigned int address = tuple_in.address;

        bool hit = false;
        unsigned int local_value_counter = 0;

        /* Check if the counter has been used recently, by cycling backword to
         * catch the updated value */
        for (auto s = 0; s < BRAM_LAT; s++) {
#pragma HLS unroll
            auto g = BRAM_LAT - s - 1;
            if (local_cache_address[g] == address && local_cache_valid[g]) {
                hit = true;
                local_value_counter = local_cache_counter[g];
            }
        }

        /* Read from memory only if is not present in local cache, in this way
         * it is possible to remove the RAW dependency */
        if (!hit){
          local_value_counter = block_n_edges[address];
        }

        /* Shift everything by one position and writes the last one in memory */
        for (auto s = 0; s < BRAM_LAT - 1; s++) {
#pragma HLS unroll
          auto g = BRAM_LAT - s - 1;
          local_cache_address[g] = local_cache_address[g - 1];
          local_cache_counter[g] = local_cache_counter[g - 1];
          local_cache_valid[g] = local_cache_valid[g - 1];
        }
        local_cache_address[0] = address;
        local_cache_counter[0] = local_value_counter + 1;
        local_cache_valid[0] = true;
        block_n_edges[address] = local_value_counter + 1;
        m_axi[local_value_counter] = tuple_in.edge;
        tuple_in = stream_edge[select].read();
        select = (select + 1) % 2;

#ifndef __SYNTHESIS__
        assert(local_value_counter < UINT32_MAX);
#endif
    }
}

template<size_t NODE_W,
         size_t LAB_W,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t MAX_LABELS>
void
storeEdgePerBlockWrap(row_t* edge_buf,
                      row_t* block_buf,
                      const unsigned char hash1_w,
                      const unsigned char hash2_w,
                      const ap_uint<8> labelToTable[MAX_LABELS][MAX_LABELS],
                      const unsigned int numDataEdges,
                      unsigned int block_n_edges[4096])
{
#pragma HLS dataflow
    hls::stream<store_tuple_t<row_t>, 30> stream_edge[2];

    readAndStreamEdgesPerBlock<NODE_W,
                               LAB_W,
                               LKP3_HASH_W,
                               MAX_HASH_W,
                               MAX_LABELS>(edge_buf,
                                           hash1_w,
                                           hash2_w,
                                           labelToTable,
                                           numDataEdges,
                                           stream_edge);

    storeEdgesPerBlock<NODE_W, LAB_W, LKP3_HASH_W, MAX_HASH_W, MAX_LABELS>(
      stream_edge, block_buf, block_n_edges);
}

template<size_t NODE_W,
         size_t ROW_LOG,
         size_t EDGE_LOG,
         size_t LAB_W,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t MAX_LABELS>
void
blockToHTB(row_t* edge_buf,
           row_t* htb_buf,
           AdjHT* hTables,
           const unsigned char hash1_w,
           const unsigned char hash2_w,
           const unsigned int numTables,
           const unsigned int block_n_edges[4096])
{
    constexpr size_t COUNTERS_PER_BLOCK = 14;
    constexpr size_t IXG_NODE = 0;
    constexpr size_t IXD_NODE = 32;
    constexpr size_t IXG_HASH = 64;
    constexpr size_t IXD_HASH = 96;
    const unsigned int block_per_table =
      (1UL << (hash1_w + hash2_w - COUNTERS_PER_BLOCK));
    ap_uint<64> block_counter0[4096];
    ap_uint<64> block_counter1[4096];
#pragma HLS bind_storage variable=block_counter0 type=RAM_2P impl=URAM
#pragma HLS bind_storage variable=block_counter1 type=RAM_2P impl=URAM

    ap_uint<64> *htb_p0 = (ap_uint<64> *)htb_buf;

/* Loop 2^(COUNTERS_PER_BLOCK - 2) since one block is split in two memories and
 * every memroy keep two counters per line*/
INITIALIZE_URAM_LOOP:
    for (auto g = 0; g < (1UL << (COUNTERS_PER_BLOCK - 2)); g++) {
#pragma HLS pipeline II = 1
        block_counter0[g] = 0;
        block_counter1[g] = 0;
    }

    auto prev_offset = 0;
    auto prev_ntb = 0;
    unsigned int base_address = 0;
BLOCK_HTB_TOP_LOOP:
    for (auto s = 0; s < block_per_table * numTables; s++) {
        auto block_edges = block_n_edges[s] - prev_offset;
        auto ntb = s >> (hash1_w + hash2_w - COUNTERS_PER_BLOCK);
        if (prev_ntb != ntb){
            base_address = 0;
        }

COUNT_EDGES_INSIDE_BLOCK_LOOP:
        for (auto g = 0; g < block_edges; g++) {
#pragma HLS pipeline II = 2
            row_t edge = edge_buf[g + prev_offset];

            ap_uint<NODE_W> indexing_hash =
              edge.range(IXG_HASH + NODE_W - 1, IXG_HASH);
            ap_uint<NODE_W> indexed_hash =
              edge.range(IXD_HASH + NODE_W - 1, IXD_HASH);

            /* Computing the bucket in which the edge will be stored, 
            restricted to the block of counter in memory */
            ap_uint<COUNTERS_PER_BLOCK> address = indexing_hash;
            address <<= hash2_w;
            address += indexed_hash;
            ap_uint<64> row_counter0;
            ap_uint<64> row_counter1;
            ap_uint<64> row_counter;

            /* The first bit select which counter in the 64-bit word, the second
             * one select on which memory read. In this way counters are
             * consecutive inside a word */
            row_counter1 = block_counter1[(address >> 2)];
            row_counter0 = block_counter0[(address >> 2)];

            if (address.test(1)){
                row_counter = row_counter1;
            } else {
                row_counter = row_counter0;
            }

            if (address.test(0)){
                ap_uint<32> counter = row_counter.range(63, 32);
                row_counter.range(63, 32) = counter + 1;
            } else {
                ap_uint<32> counter = row_counter.range(31, 0);
                row_counter.range(31, 0) = counter + 1;
            }

            if (address.test(1)){
                row_counter1 = row_counter;
            } else {
                row_counter0 = row_counter;
            }

            block_counter1[(address >> 2)] = row_counter1;
            block_counter0[(address >> 2)] = row_counter0;
        }

COUNTERS_TO_OFFSETS_URAM_LOOP:
        for (auto g = 0; g < (1UL << (COUNTERS_PER_BLOCK - 2)); g++){
#pragma HLS pipeline II = 2
            ap_uint<64> counter0 = block_counter0[g];
            ap_uint<64> counter1 = block_counter1[g];
            ap_uint<64> offset0, offset1;
            offset0.range(31, 0) = base_address;
            offset0.range(63, 32) = base_address + counter0.range(31, 0);
            offset1.range(31, 0) =
              base_address + counter0.range(31, 0) + counter0.range(63, 32);
            offset1.range(63, 32) = base_address + counter0.range(31, 0) +
                                    counter0.range(63, 32) +
                                    counter1.range(31, 0);
            base_address += counter0.range(31, 0) + counter0.range(63, 32) +
                            counter1.range(31, 0) + counter1.range(63, 32);
            block_counter0[g] = offset0;
            block_counter1[g] = offset1;
        }

STORE_EDGES_INSIDE_BLOCK_LOOP:
        for (auto g = 0; g < block_edges; g++) {
#pragma HLS pipeline II = 2
            row_t edge = edge_buf[g + prev_offset];
            
            ap_uint<NODE_W> indexing_hash =
              edge.range(IXG_HASH + NODE_W - 1, IXG_HASH);
            ap_uint<NODE_W> indexed_hash =
              edge.range(IXD_HASH + NODE_W - 1, IXD_HASH);
            ap_uint<NODE_W> indexed_node =
              edge.range(IXD_NODE + NODE_W - 1, IXD_NODE);
            ap_uint<NODE_W> indexing_node =
              edge.range(IXG_NODE + NODE_W - 1, IXG_NODE);

            /* Computing the bucket in which the edge will be stored, 
            restricted to the block of offset in memory */
            ap_uint<COUNTERS_PER_BLOCK> address = indexing_hash;
            address <<= hash2_w;
            address += indexed_hash;
            ap_uint<64> row_offset0;
            ap_uint<64> row_offset1;
            ap_uint<64> row_offset;
            ap_uint<32>  offset;
            
            /* The first bit select which offset in the 64-bit word, the second
             * one select on which memory read. In this way OFFSETS are
             * consecutive inside a word */
            row_offset1 = block_counter1[(address >> 2)];
            row_offset0 = block_counter0[(address >> 2)];

            if (address.test(1)){
                row_offset = row_offset1;
            } else {
                row_offset = row_offset0;
            }

            if (address.test(0)){
                offset = row_offset.range(63, 32);
                row_offset.range(63, 32) = offset + 1;
            } else {
                offset = row_offset.range(31, 0);
                row_offset.range(31, 0) = offset + 1;
            }

            if (address.test(1)){
                row_offset1 = row_offset;
            } else {
                row_offset0 = row_offset;
            }

            block_counter1[(address >> 2)] = row_offset1;
            block_counter0[(address >> 2)] = row_offset0;
            
            ap_uint<32> addr_row_offset = (hTables[ntb].start_edges << 1) + offset;
            htb_p0[addr_row_offset] = indexing_node.concat(indexed_node);
        }

        /* Store the block counters, packing them in a row */
        row_t row;
STORE_OFFSETS_BLOCK_LOOP:
        for (auto g = 0; g < (1UL << (COUNTERS_PER_BLOCK - 2)); g++) {
#pragma HLS pipeline II = 1
            row.range(63, 0) = block_counter0[g];
            row.range(127, 64) = block_counter1[g];
            block_counter0[g] = 0;
            block_counter1[g] = 0;
            htb_buf[g + (s * (1UL << (COUNTERS_PER_BLOCK - 2)))] = row;
        }
        prev_offset = block_n_edges[s];
        prev_ntb = ntb;
    }
}

template<typename T_CNT,
         size_t ROW_LOG,
         size_t EDGE_LOG,
         size_t NODE_W,
         size_t LAB_W,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t MAX_LABELS>
void
storeEdgesHTB(row_t* edge_buf,
              row_t* htb_buf,
              AdjHT* hTables,
              const unsigned char hash1_w,
              const unsigned char hash2_w,
              const unsigned int numTables,
              const unsigned int block_n_edges[4096])
{
    constexpr size_t OFFSETS_PER_BLOCK = 14;
    constexpr size_t IXG_NODE = 0;
    constexpr size_t IXD_NODE = 32;
    constexpr size_t IXG_HASH = 64;
    constexpr size_t IXD_HASH = 96;
    const unsigned int block_per_table = (1UL << (hash1_w + hash2_w - OFFSETS_PER_BLOCK));
    ap_uint<64> block_offset0[4096];
    ap_uint<64> block_offset1[4096];
#pragma HLS bind_storage variable=block_offset0 type=RAM_2P impl=URAM
#pragma HLS bind_storage variable=block_offset1 type=RAM_2P impl=URAM

    auto prev_offset = 0;
STORE_EDGES_TOP_LOOP:
    for (auto s = 0; s < block_per_table * numTables; s++) {
        auto block_edges = block_n_edges[s] - prev_offset;
        auto ntb = s >> (hash1_w + hash2_w - OFFSETS_PER_BLOCK);
        
        row_t row;
LOAD_URAM_OFFSETS_LOOP:
        for (auto g = 0; g < (1UL << (OFFSETS_PER_BLOCK - 2)); g++) {
#pragma HLS pipeline II = 1
            row = htb_buf[g + (s * (1UL << (OFFSETS_PER_BLOCK - 2)))];
            block_offset0[g] = row.range(63, 0);
            block_offset1[g] = row.range(127, 64);
        }

STORE_EDGES_BLOCK_LOOP:
        for (auto g = 0; g < block_edges; g++) {
#pragma HLS pipeline II = 2
            row_t edge = edge_buf[g + prev_offset];
            
            ap_uint<NODE_W> indexing_hash =
              edge.range(IXG_HASH + NODE_W - 1, IXG_HASH);
            ap_uint<NODE_W> indexed_hash =
              edge.range(IXD_HASH + NODE_W - 1, IXD_HASH);
            ap_uint<NODE_W> indexed_node =
              edge.range(IXD_NODE + NODE_W - 1, IXD_NODE);
            ap_uint<NODE_W> indexing_node =
              edge.range(IXG_NODE + NODE_W - 1, IXG_NODE);

            /* Computing the bucket in which the edge will be stored, 
            restricted to the block of offset in memory */
            ap_uint<OFFSETS_PER_BLOCK> address = indexing_hash;
            address <<= hash2_w;
            address += indexed_hash;
            ap_uint<64> row_offset0;
            ap_uint<64> row_offset1;
            ap_uint<64> row_offset;
            ap_uint<32>  offset;
            
            /* The first bit select which offset in the 64-bit word, the second
             * one select on which memory read. In this way OFFSETS are
             * consecutive inside a word */
            row_offset1 = block_offset1[(address >> 2)];
            row_offset0 = block_offset0[(address >> 2)];

            if (address.test(1)){
                row_offset = row_offset1;
            } else {
                row_offset = row_offset0;
            }

            if (address.test(0)){
                offset = row_offset.range(63, 32);
                row_offset.range(63, 32) = offset + 1;
            } else {
                offset = row_offset.range(31, 0);
                row_offset.range(31, 0) = offset + 1;
            }

            if (address.test(1)){
                row_offset1 = row_offset;
            } else {
                row_offset0 = row_offset;
            }

            block_offset1[(address >> 2)] = row_offset1;
            block_offset0[(address >> 2)] = row_offset0;
            
            /* Compute address of row that will store the edge */
            T_CNT addr_row_edge =
              hTables[ntb].start_edges + (offset >> (ROW_LOG - EDGE_LOG));

            /* Compute address of the edge inside the row */
            T_CNT addr_inrow = offset.range((ROW_LOG - EDGE_LOG) - 1, 0);

            /* Read, modify and write the edge */
            row_t row_edge = htb_buf[addr_row_edge];
            row_edge.range(((addr_inrow + 1) << EDGE_LOG) - 1,
                           addr_inrow << EDGE_LOG) =
              indexing_node.concat(indexed_node);

            /* Store offset and edge modified */
            htb_buf[addr_row_edge] = row_edge;
        }

        /* Store the block OFFSETS packing them in a row */
STORE_OFFSETS_BLOCK_LOOP:
        for (auto g = 0; g < (1UL << (OFFSETS_PER_BLOCK - 2)); g++) {
#pragma HLS pipeline II = 1
            row.range(63, 0) = block_offset0[g];
            row.range(127, 64) = block_offset1[g];
            htb_buf[g + (s * (1UL << (OFFSETS_PER_BLOCK - 2)))] = row;
        }
        prev_offset = block_n_edges[s];
    }
}

/* Function in charge of reading the starting vertices of partial solutions. 
 * While reading an indexing set, it is critical to transform it from a bag 
 * in which nodes are repeated in a set. */
template<size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t FULL_HASH_W,
         size_t NODE_W,
         size_t EDGE_LOG,
         size_t ROW_LOG,
         size_t MAX_CL>
void
mwj_batch(const unsigned char hash1_w,
          unsigned int &n_candidate,
          const unsigned int start_address,
          AdjHT* hTables,
          QueryVertex* qVertices,
          row_t* htb_buf)
{
  ap_uint<8> tableIndex = 0;
  ap_uint<32> minSize = (1UL << 32) - 1;
  ap_uint<32> minOff;
  ap_uint<NODE_W * 2> edge;
  ap_uint<NODE_W> vertex;
  ap_uint<NODE_W> set[MAX_CL];
  ap_uint<MAX_HASH_W> hash_buff, hash_new;
  unsigned char set_counter = 0;
  bool flag_buff = false;
  bool flag_new = true;
  hash_buff = hash_new = 0;
  // unsigned int rm_start = 0;

PROPOSE_TBINDEXING_LOOP:
  for (int g = 0; g < qVertices[0].numTablesIndexing; g++) {
    tableIndex = qVertices[0].tables_indexing[g];

    if (hTables[tableIndex].n_edges < minSize) {
      minSize = hTables[tableIndex].n_edges;
      minOff = hTables[tableIndex].start_edges;
    }
  }

  unsigned int rowstart = minOff;
  unsigned int rowend = minOff + (minSize >> (ROW_LOG - EDGE_LOG));
  unsigned int window_right =
    minSize.range((ROW_LOG - EDGE_LOG) - 1, 0) + ((rowend - rowstart) << (ROW_LOG - EDGE_LOG));
  unsigned int cnt = 0;
  unsigned int address = start_address;
  n_candidate = 0;

PROPOSE_READ_MIN_INDEXING_LOOP:
  for (unsigned int g = 0; g <= rowend - rowstart; g++) {
    row_t row = htb_buf[rowstart + g];
    for (unsigned int i = 0; i < EDGE_ROW; i++, cnt++) {
#pragma HLS unroll
      if (cnt < window_right) {
        edge = row.range((1UL << EDGE_LOG) - 1, 0);
        vertex = edge.range(NODE_W * 2 - 1, NODE_W);
        ap_uint<LKP3_HASH_W> hash_out;
        xf::database::details::hashlookup3_core<NODE_W>(vertex, hash_out);
        hash_new = hash_out.range(MAX_HASH_W - 1, 0);
        hash_new = hash_new.range(hash1_w - 1, 0);

        if (flag_buff && hash_buff == hash_new) {
          flag_new = true;
        EXTRACT_BAGTOSET_SETCHECKER_LOOP:
          for (int nSet = 0; nSet < set_counter; nSet++) {
            if (vertex == set[nSet]) {
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
          htb_buf[address++] = vertex;
          n_candidate++;
        }

        hash_buff = hash_new;
        flag_buff = true;
      }
      row >>= (1UL << EDGE_LOG);
    }
  }

#if DEBUG_STATS
    debug::batch_reads += ceil((rowend - rowstart) / 16.0);
    // std::cout << rm_start << " removed\n";
#endif
}


/* Reads two times the data graph and fills the data stuctures */
template<typename T_DDR,
         typename T_BLOOM,
         size_t EDGE_LOG,
         size_t CNT_LOG,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t ROW_LOG,
         size_t NODE_W,
         size_t NODE_LOG,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t FULL_HASH_W,
         size_t LAB_W,
         size_t STREAM_D,
         size_t HTB_SPACE,
         size_t MAX_LABELS,
         size_t MAX_CL>
void
fillTablesURAM(row_t* edge_buf,
               T_DDR* htb_buf,
               T_DDR* htb_buf1,
               T_DDR* bloom_p,
               QueryVertex* qVertices,
               AdjHT* hTables0,
               AdjHT* hTables1,
               const ap_uint<8> labelToTable[MAX_LABELS][MAX_LABELS],
               const unsigned long dynfifo_space,
               unsigned int &n_candidate,
               unsigned int &start_candidate,
               const unsigned long numDataEdges,
               const unsigned short numTables,
               const unsigned char hash1_w,
               const unsigned char hash2_w)
{

    /* Resetting portion of memory dedicated to counters
     * 1 << HASH1_W * HASH2_W is the number of counters needed
     * for each table, then it should be divided by the number
     * of counters stored in each row which is 1 << (ROW_LOG - CNT_LOG)*/
    constexpr size_t COUNTERS_PER_BLOCK = 14;
    const unsigned long htb_size =
      (1UL << (hash1_w + hash2_w - (DDR_BIT - COUNTER_WIDTH)));
    const unsigned int block_per_table =
      (1UL << (hash1_w + hash2_w - COUNTERS_PER_BLOCK));
    unsigned long start_addr = 0;
    unsigned int block_n_edges[4096];
#pragma HLS bind_storage variable = block_n_edges type = RAM_T2P impl = BRAM

#ifndef __SYNTHESIS__
    unsigned long end_addr = numTables * htb_size;
#endif

STORE_HASHTABLES_POINTER_LOOP:
    for (unsigned int ntb = 0; ntb < numTables; ntb++){
        hTables0[ntb].start_offset = start_addr;
        start_addr += htb_size;
    }
    
INITIALIZE_BRAM_LOOP:
    for (auto g = 0; g < numTables * block_per_table; g++) {
#pragma HLS pipeline II = 1
        block_n_edges[g] = 0;
    }

    countEdgesPerBlockWrap<NODE_W, LAB_W, LKP3_HASH_W, MAX_HASH_W, MAX_LABELS>(
      &edge_buf[dynfifo_space],
      hash1_w,
      hash2_w,
      labelToTable,
      numDataEdges,
      block_n_edges);

    auto base_addr = 0;
COUNTER_TO_OFFSET_BLOCK_LOOP:
    for (auto g = 0; g < numTables * block_per_table; g++) {
#pragma HLS pipeline II = 1
        auto data = block_n_edges[g];
        block_n_edges[g] = base_addr;
        base_addr += data;
    }

    storeEdgePerBlockWrap<NODE_W, LAB_W, LKP3_HASH_W, MAX_HASH_W, MAX_LABELS>(
      &edge_buf[dynfifo_space],
      bloom_p,
      hash1_w,
      hash2_w,
      labelToTable,
      numDataEdges,
      block_n_edges);
    
    start_addr = (start_addr + (1UL << CACHE_WORDS_PER_LINE)) &
                 ~((1UL << CACHE_WORDS_PER_LINE) - 1);
    unsigned int prev_offset = 0;
STORE_EDGES_POINTER_LOOP:
    for (unsigned short ntb = 0; ntb < numTables; ntb++) {
        hTables0[ntb].start_edges = start_addr;
        unsigned int offset = block_n_edges[((ntb + 1) * block_per_table) - 1];
        hTables0[ntb].n_edges = offset - prev_offset;
        prev_offset = offset;
        start_addr += (hTables0[ntb].n_edges >> (ROW_LOG - EDGE_LOG)) + 1;
        start_addr = (start_addr + (1UL << CACHE_WORDS_PER_LINE)) &
                     ~((1UL << CACHE_WORDS_PER_LINE) - 1);
#ifndef __SYNTHESIS__
        assert(start_addr < HTB_SPACE);
#endif
        hTables1[ntb].start_offset = hTables0[ntb].start_offset;
        hTables1[ntb].start_edges = hTables0[ntb].start_edges;
        hTables1[ntb].n_edges = hTables0[ntb].n_edges;
    }

    blockToHTB<NODE_W,
               ROW_LOG,
               EDGE_LOG,
               LAB_W,
               LKP3_HASH_W,
               MAX_HASH_W,
               MAX_LABELS>(
      bloom_p, htb_buf, hTables0, hash1_w, hash2_w, numTables, block_n_edges);

    writeBloom<T_DDR,
               T_BLOOM,
               CNT_LOG,
               ROW_LOG,
               NODE_W,
               NODE_LOG,
               EDGE_LOG,
               MAX_CL,
               LKP3_HASH_W,
               MAX_HASH_W,
               FULL_HASH_W,
               BLOOM_LOG,
               K_FUN_LOG,
               STREAM_D>(bloom_p, htb_buf, htb_buf1, hTables0, qVertices, n_candidate, start_addr, numTables, hash1_w);

    // mwj_batch<LKP3_HASH_W,
    //           MAX_HASH_W,
    //           FULL_HASH_W,
    //           NODE_W,
    //           EDGE_LOG,
    //           ROW_LOG,
    //           MAX_CL>(
    //   hash1_w, n_candidate, start_addr, hTables0, qVertices, htb_buf);

    start_candidate = start_addr;
#ifndef __SYNTHESIS__
    end_addr = start_addr * (1UL << (ROW_LOG - 3)) + ((numTables * ((1 << hash1_w) + 1)) << (BLOOM_LOG - 3));
    std::cout << "Occupied " << end_addr << " bytes, " << 
        end_addr / (float)(1UL << 20) << " MB. " << std::endl;
#endif

#if DEBUG_STATS
    debug::bloom_fullness /= numTables * (1UL << hash1_w) * (1UL << BLOOM_LOG) * (1UL << K_FUN_LOG);
#endif /* DEBUG_STATS */
}

template<typename T_DDR,
         typename T_BLOOM,
         size_t EDGE_LOG,
         size_t CNT_LOG,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t ROW_LOG,
         size_t NODE_W,
         size_t NODE_LOG,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t FULL_HASH_W,
         size_t LAB_W,
         size_t STREAM_D,
         size_t HTB_SPACE,
         size_t MAX_QV,
         size_t MAX_TB,
         size_t MAX_CL>
void
preprocess(row_t* edge_buf,
           T_DDR* htb_buf0,
           T_DDR* htb_buf1,
           T_DDR* bloom_p,
           QueryVertex* qVertices,
           AdjHT* hTables0,
           AdjHT* hTables1,
           const unsigned long dynfifo_space,
           unsigned int &n_candidate,
           unsigned int &start_candidate,
           unsigned short numQueryVert,
           unsigned short numQueryEdges,
           unsigned long numDataEdges,
           const unsigned char hash1_w,
           const unsigned char hash2_w)
{
    constexpr size_t MAX_LABELS = (1UL << LAB_W);
    unsigned short numTables = 0;
    ap_uint<8> labelToTable[MAX_LABELS][MAX_LABELS];

INITIALIZE_LABELTOTABLE_LOOP:
    for (int g = 0; g < MAX_LABELS; g++) {
        for (int s = 0; s < MAX_LABELS; s++) {
#pragma HLS pipeline II = 2
            labelToTable[g][s] = 0;
        }
    }

    buildTableDescriptors<MAX_QV, MAX_TB, NODE_W, LAB_W, MAX_LABELS>(
      &edge_buf[dynfifo_space + numDataEdges],
      qVertices,
      labelToTable,
      numTables,
      numQueryVert,
      numQueryEdges);

    fillTablesURAM<T_DDR,
                   T_BLOOM,
                   EDGE_LOG,
                   CNT_LOG,
                   BLOOM_LOG,
                   K_FUN_LOG,
                   ROW_LOG,
                   NODE_W,
                   NODE_LOG,
                   LKP3_HASH_W,
                   MAX_HASH_W,
                   FULL_HASH_W,
                   LAB_W,
                   STREAM_D,
                   HTB_SPACE,
                   MAX_LABELS,
                   MAX_CL>(edge_buf,
                               htb_buf0,
                               htb_buf1,
                               bloom_p,
                               qVertices,
                               hTables0,
                               hTables1,
                               labelToTable,
                               dynfifo_space,
                               n_candidate,
                               start_candidate,
                               numDataEdges,
                               numTables,
                               hash1_w,
                               hash2_w);

}

#pragma GCC diagnostic pop
