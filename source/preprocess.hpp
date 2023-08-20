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
#include "cache.h"

#if DEBUG_STATS
#include "debug.hpp"
#endif /* DEBUG_STATS */

#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wpedantic"
#pragma GCC diagnostic error "-Wall"
#pragma GCC diagnostic error "-Wextra"
#pragma GCC diagnostic ignored "-Wunused-label"
#pragma GCC diagnostic ignored "-Wsign-compare"
// #pragma GCC diagnostic ignored "-Wunused-parameter"

template<size_t HASH_W>
struct bloom_tuple_t
{
  ap_uint<HASH_W> indexed_h;
  unsigned int address;
  bool write;
  bool valid;
};

struct counter_tuple_t
{
  unsigned int address;
  bool stop;
};

template<size_t EDGE_W>
struct store_tuple_t
{
  ap_uint<EDGE_W> edge;
  unsigned int address;
  unsigned char tb_index;
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
                      QueryVertex* qVertices0,
                      QueryVertex* qVertices1,
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
        unsigned char idx = qVertices0[nodeSrcPos].numTablesIndexing;
        qVertices0[nodeSrcPos].tables_indexing[idx] = index - 1;
        qVertices1[nodeSrcPos].tables_indexing[idx] = index - 1;
        qVertices0[nodeSrcPos].numTablesIndexing++;
        qVertices1[nodeSrcPos].numTablesIndexing++;

        idx = qVertices0[nodeDstPos].numTablesIndexed;
        qVertices0[nodeDstPos].tables_indexed[idx] = index - 1;
        qVertices0[nodeDstPos].vertex_indexing[idx] = nodeSrcPos;
        qVertices1[nodeDstPos].tables_indexed[idx] = index - 1;
        qVertices1[nodeDstPos].vertex_indexing[idx] = nodeSrcPos;

        qVertices0[nodeDstPos].numTablesIndexed++;
        qVertices1[nodeDstPos].numTablesIndexed++;
      } else {
        unsigned char idx = qVertices0[nodeDstPos].numTablesIndexing;
        qVertices0[nodeDstPos].tables_indexing[idx] = index - 1;
        qVertices1[nodeDstPos].tables_indexing[idx] = index - 1;
        qVertices0[nodeDstPos].numTablesIndexing++;
        qVertices1[nodeDstPos].numTablesIndexing++;

        idx = qVertices0[nodeSrcPos].numTablesIndexed;
        qVertices0[nodeSrcPos].tables_indexed[idx] = index - 1;
        qVertices0[nodeSrcPos].vertex_indexing[idx] = nodeDstPos;
        qVertices1[nodeSrcPos].tables_indexed[idx] = index - 1;
        qVertices1[nodeSrcPos].vertex_indexing[idx] = nodeDstPos;

        qVertices0[nodeSrcPos].numTablesIndexed++;
        qVertices1[nodeSrcPos].numTablesIndexed++;
      }
    }
}

/* Reads edges from each table and divide the indexed vertices based
on indexing hash to create bloom filters */
template<typename T_DDR,
         size_t CNT_LOG,
         size_t ROW_LOG,
         size_t NODE_W,
         size_t EDGE_LOG,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t FULL_HASH_W>
void
bloom_read(AdjHT* hTables,
           T_DDR* htb_buf,
           const unsigned short numTables,
           const unsigned char hash1_w,
           hls::stream<bloom_tuple_t<FULL_HASH_W>>& stream_tuple_out,
           hls::stream<bool>& stream_tuple_end_out)
{
    constexpr size_t EDGE_W = 1UL << EDGE_LOG;
    hls::stream<ap_uint<NODE_W>, 4> hash_in0, hash_in1;
    hls::stream<ap_uint<FULL_HASH_W>, 4> hash_out0;
    hls::stream<ap_uint<LKP3_HASH_W>, 4> hash_out1;
    ap_uint<EDGE_W> edge;
    ap_uint<NODE_W> indexing_v, indexed_v;
    ap_uint<FULL_HASH_W> indexed_h, prev_indexed_h;
    ap_uint<MAX_HASH_W> indexing_h, prev_indexing_h;
    bloom_tuple_t<FULL_HASH_W> tuple_out;
    T_DDR row;
    unsigned int counter;
    prev_indexing_h = 0;
    prev_indexed_h = 0;

BLOOM_READ_TABLES_LOOP:
    for (unsigned int ntb = 0; ntb < numTables; ntb++) {

        /* During first iteration do not consider the difference between
        prev_indexing_h and indexing_h to be useful to write the bloom */
        bool first_it = true;
        counter = 0;
        unsigned int cycles = (hTables[ntb].n_edges >> 1) + 1;
        unsigned int offset = hTables[ntb].start_edges;

        /* Read all the edges in a table and divide them by hash1 */
    BLOOM_READ_EDGES_BLOCK:
        for (unsigned int start = 0; start < cycles; start++) {
#pragma HLS pipeline II = 2

            row = htb_buf[offset + start];

            for (int i = 0; i < EDGE_ROW; i++) {
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

                /* Writing edge of previous iteration */
                tuple_out.address = ntb * (1UL << hash1_w) + prev_indexing_h;
                tuple_out.indexed_h = prev_indexed_h;
                tuple_out.write = (indexing_h != prev_indexing_h) && !first_it;
                tuple_out.valid = (counter < hTables[ntb].n_edges) && !first_it;
                stream_tuple_out.write(tuple_out);
                stream_tuple_end_out.write(false);

                if (counter < hTables[ntb].n_edges){
                  prev_indexing_h = indexing_h;
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
        tuple_out.valid = true;
        stream_tuple_out.write(tuple_out);
        stream_tuple_end_out.write(false);
    }

    stream_tuple_end_out.write(true);
}

template<typename T_BLOOM,
         size_t BLOOM_LOG,
         size_t FULL_HASH_W,
         size_t K_FUN_LOG>
void
bloom_write(T_BLOOM* bloom_p,
            hls::stream<bloom_tuple_t<FULL_HASH_W>>& stream_tuple_in,
            hls::stream<bool>& stream_tuple_end_in)
{
    constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    T_BLOOM filter[K_FUN];
    bloom_tuple_t<FULL_HASH_W> tuple_in;
#pragma HLS array_partition variable = filter type = complete

    for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
        filter[g] = 0;
    }

    bool last = stream_tuple_end_in.read();
BLOOM_WRITE_TASK_LOOP:
    while (!last) {
#pragma HLS pipeline II = 8
        tuple_in = stream_tuple_in.read();

        if (tuple_in.valid) {

            /* Forcing Vitis_HLS to insert burst write */
            if (tuple_in.write) {
                for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
                  ap_uint<BLOOM_LOG> idx = tuple_in.indexed_h.range(
                    (FULL_HASH_W / K_FUN) * (g + 1) - 1,
                    (FULL_HASH_W / K_FUN) * (g + 1) - BLOOM_LOG);
                  filter[g][idx] = 1;
                  bloom_p[(tuple_in.address << K_FUN_LOG) + g] = filter[g];
                  filter[g] = 0;
                }
            } else {
                for (int g = 0; g < K_FUN; g++) {
#pragma HLS unroll
                  ap_uint<BLOOM_LOG> idx = tuple_in.indexed_h.range(
                    (FULL_HASH_W / K_FUN) * (g + 1) - 1,
                    (FULL_HASH_W / K_FUN) * (g + 1) - BLOOM_LOG);
                  filter[g][idx] = 1;
                }
            }
        }

        last = stream_tuple_end_in.read();
    }
}

template <typename T_DDR,
          typename T_BLOOM,
          size_t CNT_LOG,
          size_t ROW_LOG,
          size_t NODE_W,
          size_t EDGE_LOG,
          size_t LKP3_HASH_W,
          size_t MAX_HASH_W,
          size_t FULL_HASH_W,
          size_t BLOOM_LOG,
          size_t K_FUN_LOG,
          size_t STREAM_D>
void writeBloom(
    T_BLOOM *bloom_p,
    T_DDR *htb_p,
    AdjHT *hTables,
    const unsigned char numTables,
    const unsigned char hash1_w)
{
#pragma HLS DATAFLOW

    hls::stream<bloom_tuple_t<FULL_HASH_W>, STREAM_D> stream_tuple("Bloom tuples");
    hls::stream<bool, STREAM_D> stream_tuple_end("Bloom tuple end");

    /* Read edges in each table */
    bloom_read<T_DDR,
               CNT_LOG,
               ROW_LOG,
               NODE_W,
               EDGE_LOG,
               LKP3_HASH_W,
               MAX_HASH_W,
               FULL_HASH_W>(
      hTables, htb_p, numTables, hash1_w, stream_tuple, stream_tuple_end);

    bloom_write<T_BLOOM, BLOOM_LOG, FULL_HASH_W, K_FUN_LOG>(
      bloom_p, stream_tuple, stream_tuple_end);
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
                  unsigned int block_n_edges[2048])
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
    };
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
                unsigned int block_n_edges[2048])
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
storeEdgePerBlock(row_t* edge_buf,
                const unsigned char hash1_w,
                const unsigned char hash2_w,
                const ap_uint<8> labelToTable[MAX_LABELS][MAX_LABELS],
                const unsigned long dynfifo_space,
                const unsigned int numDataEdges,
                unsigned int block_n_edges[2048])
{
    constexpr size_t COUNTERS_PER_BLOCK = 14;
    constexpr size_t SRC_NODE = 0;
    constexpr size_t DST_NODE = 32;
    constexpr size_t LABELSRC_NODE = 64;
    constexpr size_t LABELDST_NODE = 96;
    const unsigned int block_per_table = hash1_w + hash2_w - COUNTERS_PER_BLOCK;

STORE_EDGE_PER_BLOCK_LOOP:
    for (auto s = 0; s < numDataEdges; s++) {
        /* We can remove this dependency since writing and reading are on
        different portion of memory */
#pragma HLS dependence variable = edge_buf type = inter direction = RAW false
#pragma HLS dependence variable = block_n_edges type = intra false
#pragma HLS pipeline II = 4
        row_t edge = edge_buf[dynfifo_space + s];

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
        unsigned char step0 = 1; 
        unsigned char step1 = 1; 
        step0 = (index0 != 0) ? step0 : 0;
        step1 = (index1 != 0) ? step1 : 0;

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
        unsigned int address0 = (1UL << block_per_table) * (index0 - 1) + address_intable0;
        unsigned short address_intable1 =
          hashdst.range(hash1_w - 1, hash1_w - block_per_table);
        unsigned int address1 = (1UL << block_per_table) * (index1 - 1) + address_intable1;
#ifndef __SYNHTESIS__
        address0 = (index0 == 0) ? 0 : address0;
        address1 = (index1 == 0) ? 0 : address1;
#endif
        unsigned int offset0 = block_n_edges[address0];
        unsigned int offset1 = block_n_edges[address1];

        /* The case in which the nodes have the same labels and ends up in the same block 
        will write the same edge two times, but with different vertex indexing */
        if (address0 == address1 && index0 != 0 && index1 != 0){
            step0 = step0 << 1;
            offset1 += 1;
        }

        if (index0 != 0) {
            row_t table_edge;
            table_edge.range(31, 0) = nodesrc;
            table_edge.range(63, 32) = nodedst;
            table_edge.range(95, 64) = hashsrc;
            table_edge.range(127, 96) = hashdst.range(hash2_w - 1, 0);
            edge_buf[offset0] = table_edge;
            block_n_edges[address0] = offset0 + step0;
        }

        if (index1 != 0) {
            row_t table_edge;
            table_edge.range(31, 0) = nodedst;
            table_edge.range(63, 32) = nodesrc;
            table_edge.range(95, 64) = hashdst;
            table_edge.range(127, 96) = hashsrc.range(hash2_w - 1, 0);
            edge_buf[offset1] = table_edge;
            block_n_edges[address1] = offset1 + step1;
        }


#ifndef __SYNTHESIS__
        assert(block_n_edges[address0] < UINT32_MAX);
        assert(block_n_edges[address1] < UINT32_MAX);
#endif
    }
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
           const unsigned int block_n_edges[2048])
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
            
            /* Compute address of row that will store the edge */
            ap_uint<32> addr_row_edge =
              hTables[ntb].start_edges + (offset >> (ROW_LOG - EDGE_LOG));
            
            /* Compute address of the edge inside the row */
            ap_uint<32> addr_inrow = offset.range((ROW_LOG - EDGE_LOG) - 1, 0);

            /* Read, modify and write the edge */
            row_t row_edge = htb_buf[addr_row_edge];
            row_edge.range(((addr_inrow + 1) << EDGE_LOG) - 1,
                           addr_inrow << EDGE_LOG) =
              indexing_node.concat(indexed_node);

            /* Store offset and edge modified */
            htb_buf[addr_row_edge] = row_edge;
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
              const unsigned int block_n_edges[2048])
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

    /* Reads two times the data graph and fills the data stuctures */
template<typename T_DDR,
         typename T_BLOOM,
         size_t EDGE_LOG,
         size_t CNT_LOG,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t ROW_LOG,
         size_t NODE_W,
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t FULL_HASH_W,
         size_t LAB_W,
         size_t STREAM_D,
         size_t HTB_SPACE,
         size_t MAX_LABELS>
void
fillTablesURAM(row_t* edge_buf,
            T_DDR* htb_buf,
            T_DDR* bloom_p,
            AdjHT* hTables0,
            AdjHT* hTables1,
            const ap_uint<8> labelToTable[MAX_LABELS][MAX_LABELS],
            const unsigned long dynfifo_space,
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
    unsigned int block_n_edges[2048];
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

    storeEdgePerBlock<NODE_W, LAB_W, LKP3_HASH_W, MAX_HASH_W, MAX_LABELS>(
      edge_buf,
      hash1_w,
      hash2_w,
      labelToTable,
      dynfifo_space,
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
      edge_buf, htb_buf, hTables0, hash1_w, hash2_w, numTables, block_n_edges);

    /* From counts to offsets */
    // counterToOffset<T_DDR, CNT_LOG, ROW_LOG, MAX_HASH_W>(
    //   htb_size, numTables, hTables0, htb_buf);


    // storeEdgesHTB<ap_uint<(1UL << CNT_LOG)>,
    //               ROW_LOG,
    //               EDGE_LOG,
    //               NODE_W,
    //               LAB_W,
    //               LKP3_HASH_W,
    //               MAX_HASH_W,
    //               MAX_LABELS>(edge_buf,
    //                           htb_buf,
    //                           hTables0,
    //                           hash1_w,
    //                           hash2_w,
    //                           numTables,
    //                           block_n_edges);

    writeBloom<T_DDR,
               T_BLOOM,
               CNT_LOG,
               ROW_LOG,
               NODE_W,
               EDGE_LOG,
               LKP3_HASH_W,
               MAX_HASH_W,
               FULL_HASH_W,
               BLOOM_LOG,
               K_FUN_LOG,
               STREAM_D>(bloom_p, htb_buf, hTables0, numTables, hash1_w);

#ifndef __SYNTHESIS__
    end_addr = start_addr * (1UL << (ROW_LOG - 3)) + ((numTables * ((1 << hash1_w) + 1)) << (BLOOM_LOG - 3));
    std::cout << "Occupied " << end_addr << " bytes, " << 
        end_addr / (float)(1UL << 20) << " MB. " << std::endl;
#endif

#if DEBUG_STATS
    constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    for (unsigned int tab = 0; tab < numTables; tab++) {
        for (unsigned long addr = 0; addr < (1UL << hash1_w); addr++) {
            T_BLOOM row = bloom_p[addr + tab * (1UL << hash1_w)];
            for (int g = 0; g < K_FUN; g++) {
                ap_uint<(1UL << (BLOOM_LOG - K_FUN_LOG))> bloom_f =
                  row.range(((g + 1) << (BLOOM_LOG - K_FUN_LOG)) - 1,
                            g << (BLOOM_LOG - K_FUN_LOG));
                unsigned int count = 0;
                while (bloom_f > 0) {
                    count++;
                    bloom_f = bloom_f & (bloom_f - 1);
                }
                debug::bloom_fullness +=
                  (float)count / ((1UL << hash1_w) * numTables);
            }
        }
    }
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
         size_t LKP3_HASH_W,
         size_t MAX_HASH_W,
         size_t FULL_HASH_W,
         size_t LAB_W,
         size_t STREAM_D,
         size_t HTB_SPACE,
         size_t MAX_QV,
         size_t MAX_TB>
void
preprocess(row_t* edge_buf,
           T_DDR* htb_buf,
           T_DDR* bloom_p,
           QueryVertex* qVertices0,
           QueryVertex* qVertices1,
           AdjHT* hTables0,
           AdjHT* hTables1,
           const unsigned long dynfifo_space,
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
    for (int g = 0; g < MAX_TABLES; g++) {
        for (int s = 0; s < MAX_TABLES; s++) {
#pragma HLS pipeline II = 2
            labelToTable[g][s] = 0;
        }
    }

    buildTableDescriptors<MAX_QV, MAX_TB, NODE_W, LAB_W, MAX_LABELS>(
      &edge_buf[dynfifo_space + numDataEdges],
      qVertices0,
      qVertices1,
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
                   LKP3_HASH_W,
                   MAX_HASH_W,
                   FULL_HASH_W,
                   LAB_W,
                   STREAM_D,
                   HTB_SPACE,
                   MAX_LABELS>(edge_buf,
                               htb_buf,
                               bloom_p,
                               hTables0,
                               hTables1,
                               labelToTable,
                               dynfifo_space,
                               numDataEdges,
                               numTables,
                               hash1_w,
                               hash2_w);
    
    // fillTables<T_DDR,
    //                T_BLOOM,
    //                EDGE_LOG,
    //                CNT_LOG,
    //                BLOOM_LOG,
    //                K_FUN_LOG,
    //                ROW_LOG,
    //                NODE_W,
    //                LKP3_HASH_W,
    //                MAX_HASH_W,
    //                FULL_HASH_W,
    //                LAB_W,
    //                STREAM_D,
    //                HTB_SPACE,
    //                MAX_LABELS>(&edge_buf[dynfifo_space],
    //                            htb_buf,
    //                            bloom_p,
    //                            hTables0,
    //                            hTables1,
    //                            labelToTable,
    //                            numDataEdges,
    //                            numTables,
    //                            hash1_w,
    //                            hash2_w);

    // exit(-1);
}

#pragma GCC diagnostic pop
