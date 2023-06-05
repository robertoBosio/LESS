#pragma once

#define HLS_STREAM_THREAD_SAFE
#ifndef __SYNTHESIS__
#include <cassert>
#include <fstream>
#endif

#include <hls_stream.h>
#include <ap_int.h>

#include "Parameters.hpp"
#include "QueryVertex.hpp"
#include "Trie.hpp"
#include "hash_lookup3.hpp"
#include "hls_burst_maxi.h"

#if DEBUG_STATS
#include "debug.hpp"
#endif /* DEBUG_STATS */

/* Builds the table descriptors based on the information
 * from the query graph. */
template <size_t MAX_QV,
         size_t MAX_TB,
         size_t NODE_W,
         size_t LAB_W>
void buildTableDescriptors(
        edge_t *edge_buf,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        TableDescriptor *tDescriptors,
        unsigned short &numTables,
        unsigned short numQueryVert,
        unsigned short numQueryEdge)
{
    /* Translate from id of vertex to position in the order */
    unsigned short fromNumToPos[MAX_QV];

    /* Filling information about query vertices and coping
     * the vertex order needed by multiway join */
FILL_ORDER_LOOP:
    for (int g = 0; g < numQueryVert; g++){
        ap_uint<NODE_W> node = edge_buf[g].src;
#ifndef __SYNTHESIS__
        assert(numQueryVert < (MAX_QV));
#endif
        fromNumToPos[node] = g;
    }

    /* Creating table descriptors */
CREATE_TABDESC_LOOP:
    for (int s = 0; s < numQueryEdge; s++) {
        bool dirEdge = false;
        edge_t edge = edge_buf[s + numQueryVert];
        int temp = LAB_W;
        ap_uint<LAB_W> labeldst = edge.dst.range(LAB_W - 1, 0);
		ap_uint<LAB_W> labelsrc = edge.src.range(LAB_W - 1, 0);
		ap_uint<NODE_W> nodedst = edge.dst >> LAB_W;
		ap_uint<NODE_W> nodesrc = edge.src >> LAB_W;
        unsigned short nodeSrcPos = fromNumToPos[nodesrc];
        unsigned short nodeDstPos = fromNumToPos[nodedst];

        if (nodeSrcPos < nodeDstPos)
            dirEdge = true;

#ifndef __SYNTHESIS__
        std::cout << (unsigned int)nodesrc << "(" << (int)labelsrc << ")"
            << " -> " <<  (unsigned int)nodedst << "(" << (int)labeldst << ")" 
            << std::endl;
#endif

        /* Understanding if table already exists */
        unsigned int g = 0;
FIND_CORRECT_TABLE_LOOP:
        for (; g < numTables; g++){
            if (tDescriptors[g].src_label == labelsrc 
                    && tDescriptors[g].dst_label == labeldst 
                    && tDescriptors[g].dir == dirEdge)
                break;
        }
     
        /* Add new table descriptor */   
        if (g == numTables){
#ifndef __SYNTHESIS__
            assert(numTables < MAX_TB);
#endif
            tDescriptors[g].src_label = labelsrc;
            tDescriptors[g].dst_label = labeldst;
            tDescriptors[g].dir = dirEdge;
            numTables++;
        }

#ifndef __SYNTHESIS__
        if (dirEdge){
        std::cout << "Table " << (int)g << ": " << (int)labelsrc << 
             " -> " << (int)labeldst << std::endl;
        } else {
        std::cout << "Table " << (int)g << ": " << (int)labeldst << 
             " <- " << (int)labelsrc << std::endl;
        }
#endif

        /* Linking vertices to tables */
        if (dirEdge){
            unsigned char idx = qVertices0[nodeSrcPos].numTablesIndexing;
            qVertices0[nodeSrcPos].tables_indexing[idx] = g;
            qVertices1[nodeSrcPos].tables_indexing[idx] = g;
            qVertices0[nodeSrcPos].numTablesIndexing++;
            qVertices1[nodeSrcPos].numTablesIndexing++;

            idx = qVertices0[nodeDstPos].numTablesIndexed;
            qVertices0[nodeDstPos].tables_indexed[idx] = g;
            qVertices0[nodeDstPos].vertex_indexing[idx] = nodeSrcPos;
            qVertices1[nodeDstPos].tables_indexed[idx] = g;
            qVertices1[nodeDstPos].vertex_indexing[idx] = nodeSrcPos;
            
            qVertices0[nodeDstPos].numTablesIndexed++;
            qVertices1[nodeDstPos].numTablesIndexed++;
        } else {
            unsigned char idx = qVertices0[nodeDstPos].numTablesIndexing;
            qVertices0[nodeDstPos].tables_indexing[idx] = g;
            qVertices1[nodeDstPos].tables_indexing[idx] = g;
            qVertices0[nodeDstPos].numTablesIndexing++;
            qVertices1[nodeDstPos].numTablesIndexing++;

            idx = qVertices0[nodeSrcPos].numTablesIndexed;
            qVertices0[nodeSrcPos].tables_indexed[idx] = g;
            qVertices0[nodeSrcPos].vertex_indexing[idx] = nodeDstPos;
            qVertices1[nodeSrcPos].tables_indexed[idx] = g;
            qVertices1[nodeSrcPos].vertex_indexing[idx] = nodeDstPos;
            
            qVertices0[nodeSrcPos].numTablesIndexed++;
            qVertices1[nodeSrcPos].numTablesIndexed++;
        }
    }
}


/* Translating counter addresses to ram addresses and
 * increasing of one the counter specified. */
template <typename T_DDR,
         size_t CNT_LOG,
         size_t ROW_LOG,
         size_t NODE_W,
         size_t HASH1_W,
         size_t HASH2_W>
void increaseCounter(
        T_DDR *htb_buf,
        AdjHT *hTables,
        hls::stream<ap_uint<HASH1_W>> &stream_index1,
        hls::stream<ap_uint<HASH2_W>> &stream_index2,
        hls::stream<unsigned short> &stream_ntable,
        hls::stream<bool> &stream_end)
{
    const size_t CNT_ROW = 1UL << (ROW_LOG - CNT_LOG);
    T_DDR ram_row, mask;
    ap_uint<HASH1_W> index1;
    ap_uint<HASH2_W> index2;
    ap_uint<(1UL << CNT_LOG)> counter;
    unsigned short ntb;
    ap_uint<64> addr_outrow;
    ap_uint<64> addr_inrow;
    ap_uint<64> addr_counter;

    bool last = stream_end.read();
INCREASE_COUNTER_DDR_LOOP:
    while(!last){
        index1 = stream_index1.read();
        index2 = stream_index2.read();
        ntb = stream_ntable.read();

        /* Address in the matrix [HASH1_W][HASH2_W] */
        addr_counter = index1;
        addr_counter <<= HASH2_W;
        addr_counter += index2;

        /* Compute address of row storing the counter */
        addr_outrow = hTables[ntb].start_offset + (addr_counter >> (ROW_LOG - CNT_LOG));

        /* Compute address of counter inside the row */
        addr_inrow = addr_counter.range((ROW_LOG - CNT_LOG) - 1, 0);

        /* Read, modify and write the counter */
        ram_row = htb_buf[addr_outrow];
        mask = ram_row >> (addr_inrow << CNT_LOG);
        mask += 1;
        mask <<= (addr_inrow << CNT_LOG);
        ram_row <<= ((CNT_ROW - addr_inrow) << CNT_LOG);
        ram_row >>= ((CNT_ROW - addr_inrow) << CNT_LOG);
        htb_buf[addr_outrow] = ram_row | mask;

        hTables[ntb].n_edges++;     
        last = stream_end.read();
    }
}

/* Transform counters to offsets. */
template <typename T_DDR,
         size_t CNT_LOG,
         size_t ROW_LOG,
         size_t HASH1_W,
         size_t HASH2_W>
void counterToOffset(
        unsigned short numTables,
        AdjHT *hTables,
        T_DDR *htb_buf)
{
    const size_t CNT_ROW = 1UL << (ROW_LOG - CNT_LOG);
    T_DDR row, row_new;

COUNTER_TO_OFFSET_DDR_LOOP:
    for (unsigned int ntb = 0; ntb < numTables; ntb++){
        ap_uint<(1UL << CNT_LOG)> base_addr = 0;
        ap_uint<(1UL << CNT_LOG)> counter;

COUNTER_TO_OFFSET_TABLE_LOOP:
        for(unsigned int start = 0; start < HTB_SIZE; start++){
            row = htb_buf[start + hTables[ntb].start_offset];

COUNTER_TO_OFFSET_ROW_LOOP:
            for (int g = 0; g < CNT_ROW; g++){
#pragma HLS unroll
                counter = row.range((1UL << CNT_LOG) - 1, 0);
                row_new >>= (1UL << CNT_LOG);
                row >>= (1UL << CNT_LOG);

#if DEBUG_STATS
                if (counter > debug::max_collisions) 
                    debug::max_collisions = counter;
                debug::avg_collisions += (float)counter / (1 << (HASH1_W + HASH2_W));
#endif /* DEBUG_STATS */

                row_new.range((CNT_ROW << CNT_LOG) - 1, (CNT_ROW - 1) << CNT_LOG) = base_addr;
                base_addr += counter;
            }

            htb_buf[start + hTables[ntb].start_offset] = row_new;
        }
    }       
}

/* Store edges based on offsets */
template <typename T_DDR,
         size_t EDGE_LOG,
         size_t CNT_LOG,
         size_t ROW_LOG,
         size_t NODE_W,
         size_t HASH1_W,
         size_t HASH2_W>
void storeEdges(
        AdjHT                                       *hTables,
        T_DDR                                       *htb_buf,
        hls::stream<ap_uint<(1UL << EDGE_LOG)>>     &stream_edge,
        hls::stream<ap_uint<HASH1_W>>               &stream_index1,
        hls::stream<ap_uint<HASH2_W>>               &stream_index2,
        hls::stream<unsigned short>                 &stream_ntable,
        hls::stream<bool>                           &stream_end)
{
    const size_t CNT_ROW = 1UL << (ROW_LOG - CNT_LOG);
    T_DDR row, mask;
    ap_uint<HASH1_W> index1;
    ap_uint<HASH2_W> index2;
    ap_uint<(1UL << EDGE_LOG)> edge;
    ap_uint<(1UL << CNT_LOG)> offset;
    unsigned short ntb;
    ap_uint<64> addr_outrow;
    ap_uint<64> addr_inrow;
    ap_uint<64> addr_offset;

    bool last = stream_end.read();
INCREASE_COUNTER_DDR_LOOP:
    while(!last){
        index1 = stream_index1.read();
        index2 = stream_index2.read();
        ntb = stream_ntable.read();
        edge = stream_edge.read();

        /* Address in the matrix [HASH1_W][HASH2_W] */
        addr_offset = index1;
        addr_offset <<= HASH2_W;
        addr_offset += index2;

        /* Compute address of row storing the offset */
        addr_outrow = hTables[ntb].start_offset + (addr_offset >> (ROW_LOG - CNT_LOG));

        /* Compute address of offset inside the row */
        addr_inrow = addr_offset.range((ROW_LOG - CNT_LOG) - 1, 0);

        // Select the correct offset from the word by shifting
        // to the right and masking instead of using range().
        row = htb_buf[addr_outrow];
        mask = row >> (addr_inrow << CNT_LOG);
        offset = mask & ((1UL << (1UL << CNT_LOG)) - 1);

        // Add 1 to the offset to move it forward
        mask += 1;

        // Rewrite the offset in the word
        mask <<= (addr_inrow << CNT_LOG);
        row <<= ((CNT_ROW - addr_inrow) << CNT_LOG);
        row >>= ((CNT_ROW - addr_inrow) << CNT_LOG);
        htb_buf[addr_outrow] = row | mask;
        
        /* Compute address of row that will store the edge */
        addr_outrow = hTables[ntb].start_edges 
            + (offset >> (ROW_LOG - EDGE_LOG));

        /* Compute address of the edge inside the row */
        addr_inrow = offset.range((ROW_LOG - EDGE_LOG) - 1, 0);

        /* Read, modify and write the edge */
        row = htb_buf[addr_outrow];
        mask = edge;
        mask <<= (addr_inrow << EDGE_LOG);
        htb_buf[addr_outrow] = row | mask;

        last = stream_end.read();
    }
}

template <typename T_BLOOM,
        size_t BLOOM_LOG,
        size_t HASH_W,
        size_t HASH1_W,
        size_t K_FUN_LOG>
void updateBloom(
        T_BLOOM                       *bloom_p,
        hls::stream<unsigned int>     &stream_bloom_a,
        hls::stream<ap_uint<HASH_W>>  &stream_index,
        hls::stream<bool>             &stream_end)
{
    constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    
    bool last = stream_end.read();
    while (!last) {
        unsigned int bloom_a = stream_bloom_a.read();
        ap_uint<HASH_W> hash_v = stream_index.read();

#ifndef __SYNTHESIS__
        assert(((bloom_a << K_FUN_LOG ) + K_FUN) < BLOOM_SPACE);
#endif

        for (int g = 0; g < K_FUN; g++){
#pragma HLS unroll
            T_BLOOM set = bloom_p[(bloom_a << K_FUN_LOG) + g];
            ap_uint<BLOOM_LOG> idx = hash_v.range((HASH_W / K_FUN) * (g + 1) - 1,
                    (HASH_W / K_FUN) * (g + 1) - BLOOM_LOG);
            set[idx] = 1;
            bloom_p[(bloom_a << K_FUN_LOG) + g] = set;
        }
        last = stream_end.read();
    }
}

template <size_t NODE_W,
         size_t LAB_W,
         size_t HASH1_W,
         size_t HASH2_W>
void edgeToHashFirst(
        edge_t                        *edge_buf,
        TableDescriptor               *tDescriptors,
        unsigned short                numTables,
        unsigned int                  numDataEdges,

        hls::stream<unsigned int>     &stream_b_addr,
        hls::stream<ap_uint<64>>      &stream_b_index,
        hls::stream<bool>             &stream_b_end,
        hls::stream<ap_uint<HASH1_W>> &stream_c_index1,
        hls::stream<ap_uint<HASH2_W>> &stream_c_index2,
        hls::stream<unsigned short>   &stream_c_ntable,
        hls::stream<bool>             &stream_c_end)
{
    
COUNT_OCCURENCIES_TOP_LOOP:
    for (int s = 0; s < numDataEdges; s++) {
        edge_t edge = edge_buf[s];
        ap_uint<LAB_W> labeldst = edge.dst.range(LAB_W - 1, 0);
		ap_uint<LAB_W> labelsrc = edge.src.range(LAB_W - 1, 0);
		ap_uint<NODE_W> nodedst = edge.dst >> LAB_W;
		ap_uint<NODE_W> nodesrc = edge.src >> LAB_W;

        /* Finding correct table */
        unsigned int g = 0;
COUNT_OCCURENCIES_FIND_TABLE_LOOP:
        for (; g < numTables; g++){
            if (tDescriptors[g].src_label == labelsrc &&
                    tDescriptors[g].dst_label == labeldst){

                ap_uint<NODE_W> vertexIndexing, vertexIndexed; 
                if(tDescriptors[g].dir){
                    vertexIndexing = nodesrc;
                    vertexIndexed = nodedst;
                } else {
                    vertexIndexing = nodedst;
                    vertexIndexed = nodesrc;
                }

                /* Compute indices for hash table */
                ap_uint<64> hash_out;
                xf::database::details::hashlookup3_core<NODE_W>(
                        vertexIndexing,
                        hash_out);

                ap_uint<HASH1_W> indexAdj = hash_out.range(HASH1_W - 1, 0); 

                xf::database::details::hashlookup3_core<NODE_W>(
                        vertexIndexed,
                        hash_out);

                ap_uint<HASH2_W> indexEdge = hash_out.range(HASH2_W - 1, 0); 

                stream_b_addr.write(g * (1 << HASH1_W) + indexAdj);
                stream_b_index.write(hash_out);
                stream_b_end.write(false);

                stream_c_index1.write(indexAdj);
                stream_c_index2.write(indexEdge);
                stream_c_ntable.write(g);
                stream_c_end.write(false);
            }

#if UNDIRECTED   
            if (tDescriptors[g].src_label == labeldst &&
                    tDescriptors[g].dst_label == labelsrc){

                ap_uint<NODE_W> vertexIndexing, vertexIndexed; 
                if(!tDescriptors[g].dir){
                    vertexIndexing = nodesrc;
                    vertexIndexed = nodedst;
                } else {
                    vertexIndexing = nodedst;
                    vertexIndexed = nodesrc;
                }
                
                /* Compute indices for hash table */
                ap_uint<64> hash_out;
                xf::database::details::hashlookup3_core<NODE_W>(
                        vertexIndexing,
                        hash_out);

                ap_uint<HASH1_W> indexAdj = hash_out.range(HASH1_W - 1, 0); 

                xf::database::details::hashlookup3_core<NODE_W>(
                        vertexIndexed,
                        hash_out);

                ap_uint<HASH2_W> indexEdge = hash_out.range(HASH2_W - 1, 0); 
                
                stream_b_addr.write(g * (1 << HASH1_W) + indexAdj);
                stream_b_index.write(hash_out);
                stream_b_end.write(false);

                stream_c_index1.write(indexAdj);
                stream_c_index2.write(indexEdge);
                stream_c_ntable.write(g);
                stream_c_end.write(false);
            }
#endif

        }
    }
    stream_b_end.write(true);
    stream_c_end.write(true);
}

template <typename T_DDR,
         typename T_BLOOM,
         size_t CNT_LOG,
         size_t ROW_LOG,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t NODE_W,
         size_t LAB_W,
         size_t HASH1_W,
         size_t HASH2_W,
         size_t STREAM_D>
void countEdges(
        edge_t          *edge_p,
        T_DDR           *htb_p,
        T_DDR           *bloom_p,
        TableDescriptor *tDescriptors,
        AdjHT           *hTables,
        unsigned short  numTables,
        unsigned long   numDataEdges)
{

#pragma HLS DATAFLOW
    hls::stream<ap_uint<HASH1_W>, STREAM_D> stream_c_index1("Counters index 1");
    hls::stream<ap_uint<HASH2_W>, STREAM_D> stream_c_index2("Counters index 2");
    hls::stream<unsigned short, STREAM_D> stream_c_ntable("Counters table");
    hls::stream<bool> stream_c_end("Counters end");

    hls::stream<unsigned int, STREAM_D> stream_b_addr("Bloom filter address");
    hls::stream<ap_uint<64>, STREAM_D> stream_b_index("Bloom index");
    hls::stream<bool> stream_b_end("Bloom end");

    // Compute hash values of vertices and find the correct table //
    edgeToHashFirst<
        NODE_W,
        LAB_W,
        HASH1_W,
        HASH2_W>(
            edge_p,
            tDescriptors,
            numTables,
            numDataEdges,
            stream_b_addr,
            stream_b_index,
            stream_b_end,
            stream_c_index1,
            stream_c_index2,
            stream_c_ntable,
            stream_c_end);

    // Update specific counter //
    increaseCounter<
        T_DDR,
        CNT_LOG,
        ROW_LOG,
        NODE_W,
        HASH1_W,
        HASH2_W>(
            htb_p,
            hTables,
            stream_c_index1,
            stream_c_index2,
            stream_c_ntable,
            stream_c_end);

    // Update indexed bloom filters //
    updateBloom<
        T_BLOOM,
        BLOOM_LOG,
        64,
        HASH1_W,
        K_FUN_LOG>(
        bloom_p,
        stream_b_addr,
        stream_b_index,
        stream_b_end);
    
}

template <size_t EDGE_LOG,
         size_t NODE_W,
         size_t LAB_W,
         size_t HASH1_W,
         size_t HASH2_W>
void edgeToHashSecond(
        edge_t                                  *edge_buf,
        TableDescriptor                         *tDescriptors,
        unsigned short                          numTables,
        unsigned int                            numDataEdges,

        hls::stream<ap_uint<(1UL << EDGE_LOG)>> &stream_c_edge,
        hls::stream<ap_uint<HASH1_W>>           &stream_c_index1,
        hls::stream<ap_uint<HASH2_W>>           &stream_c_index2,
        hls::stream<unsigned short>             &stream_c_ntable,
        hls::stream<bool>                       &stream_c_end)
{
    
COUNT_OCCURENCIES_TOP_LOOP:
    for (int s = 0; s < numDataEdges; s++) {
        edge_t edge = edge_buf[s];
        ap_uint<LAB_W> labeldst = edge.dst.range(LAB_W - 1, 0);
		ap_uint<LAB_W> labelsrc = edge.src.range(LAB_W - 1, 0);
		ap_uint<NODE_W> nodedst = edge.dst >> LAB_W;
		ap_uint<NODE_W> nodesrc = edge.src >> LAB_W;

        /* Finding correct table */
        unsigned int g = 0;
COUNT_OCCURENCIES_FIND_TABLE_LOOP:
        for (; g < numTables; g++){
            if (tDescriptors[g].src_label == labelsrc &&
                    tDescriptors[g].dst_label == labeldst){

                ap_uint<NODE_W> vertexIndexing, vertexIndexed; 
                if(tDescriptors[g].dir){
                    vertexIndexing = nodesrc;
                    vertexIndexed = nodedst;
                } else {
                    vertexIndexing = nodedst;
                    vertexIndexed = nodesrc;
                }

                /* Compute indices for hash table */
                ap_uint<64> hash_out;
                xf::database::details::hashlookup3_core<NODE_W>(
                        vertexIndexing,
                        hash_out);

                ap_uint<HASH1_W> indexAdj = hash_out.range(HASH1_W - 1, 0); 

                xf::database::details::hashlookup3_core<NODE_W>(
                        vertexIndexed,
                        hash_out);

                ap_uint<HASH2_W> indexEdge = hash_out.range(HASH2_W - 1, 0); 

                stream_c_index1.write(indexAdj);
                stream_c_index2.write(indexEdge);
                stream_c_ntable.write(g);
                stream_c_edge.write(vertexIndexing.concat(vertexIndexed));
                stream_c_end.write(false);
            }

#if UNDIRECTED   
            if (tDescriptors[g].src_label == labeldst &&
                    tDescriptors[g].dst_label == labelsrc){

                ap_uint<NODE_W> vertexIndexing, vertexIndexed; 
                if(!tDescriptors[g].dir){
                    vertexIndexing = nodesrc;
                    vertexIndexed = nodedst;
                } else {
                    vertexIndexing = nodedst;
                    vertexIndexed = nodesrc;
                }
                
                /* Compute indices for hash table */
                ap_uint<64> hash_out;
                xf::database::details::hashlookup3_core<NODE_W>(
                        vertexIndexing,
                        hash_out);

                ap_uint<HASH1_W> indexAdj = hash_out.range(HASH1_W - 1, 0); 

                xf::database::details::hashlookup3_core<NODE_W>(
                        vertexIndexed,
                        hash_out);

                ap_uint<HASH2_W> indexEdge = hash_out.range(HASH2_W - 1, 0); 
                
                stream_c_index1.write(indexAdj);
                stream_c_index2.write(indexEdge);
                stream_c_ntable.write(g);
                stream_c_edge.write(vertexIndexing.concat(vertexIndexed));
                stream_c_end.write(false);
            }
#endif

        }
    }
    stream_c_end.write(true);
}

template <typename T_DDR,
         size_t EDGE_LOG,
         size_t CNT_LOG,
         size_t ROW_LOG,
         size_t NODE_W,
         size_t LAB_W,
         size_t HASH1_W,
         size_t HASH2_W,
         size_t STREAM_D>
void writeEdges(
        edge_t          *edge_p,
        T_DDR           *htb_p,
        TableDescriptor *tDescriptors,
        AdjHT           *hTables,
        unsigned short  numTables,
        unsigned long   numDataEdges)
{

#pragma HLS DATAFLOW

    hls::stream<ap_uint<(1UL << EDGE_LOG)>, STREAM_D> stream_o_edge("Offset edge");
    hls::stream<ap_uint<HASH1_W>, STREAM_D> stream_o_index1("Offsets index 1");
    hls::stream<ap_uint<HASH2_W>, STREAM_D> stream_o_index2("Offsets index 2");
    hls::stream<unsigned short, STREAM_D> stream_o_ntable("Offset table");
    hls::stream<bool> stream_o_end("Offset end");
    
    /* Write edges based on precomputed offsets */
    edgeToHashSecond<
        EDGE_LOG,
        NODE_W,
        LAB_W,
        HASH1_W,
        HASH2_W>(
            edge_p,
            tDescriptors,
            numTables,
            numDataEdges,
            stream_o_edge,
            stream_o_index1,
            stream_o_index2,
            stream_o_ntable,
            stream_o_end);

    storeEdges<
        T_DDR,
        EDGE_LOG,
        CNT_LOG,
        ROW_LOG,
        NODE_W,
        HASH1_W,
        HASH2_W>(
            hTables,
            htb_p,
            stream_o_edge,
            stream_o_index1,
            stream_o_index2,
            stream_o_ntable,
            stream_o_end);
}

/* Reads two times the data graph and fills the data stuctures */
template <typename T_DDR,
         typename T_BLOOM,
         size_t EDGE_LOG,
         size_t CNT_LOG,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t ROW_LOG,
         size_t NODE_W,
         size_t LAB_W,
         size_t HASH1_W,
         size_t HASH2_W,
         size_t STREAM_D,
         size_t HTB_SPACE>
void fillTables(
        edge_t  *edge_buf,
        T_DDR   *htb_buf,
        T_DDR   *bloom_p,
        AdjHT   *hTables0,
        AdjHT   *hTables1,
        TableDescriptor *tDescriptors,
        unsigned long numDataEdges,
        unsigned short numTables)
{

    /* Resetting portion of memory dedicated to counters 
     * 1 << HASH1_W * HASH2_W is the number of counters needed
     * for each table, then it should be divided by the number
     * of counters stored in each row which is 1 << (ROW_LOG - CNT_LOG)*/
    unsigned long end_addr = numTables * HTB_SIZE;
    unsigned long start_addr = 0;

STORE_HASHTABLES_POINTER_LOOP:
    for (unsigned int ntb = 0; ntb < numTables; ntb++){
        hTables0[ntb].start_offset = start_addr;
        hTables0[ntb].n_edges = 0;
        start_addr += HTB_SIZE;
    }

    countEdges<
         T_DDR,
         T_BLOOM,
         CNT_LOG,
         ROW_LOG,
         BLOOM_LOG,
         K_FUN_LOG,
         NODE_W,
         LAB_W,
         HASH1_W,
         HASH2_W,
         STREAM_D>(
            edge_buf,
            htb_buf,
            bloom_p,
            tDescriptors,
            hTables0,
            numTables,
            numDataEdges);
    
    /* From counts to offsets */
    counterToOffset<T_DDR,
         CNT_LOG,
         ROW_LOG,
         HASH1_W,
         HASH2_W>(
            numTables,
            hTables0,
            htb_buf);

    start_addr = (start_addr + (1UL << CACHE_WORDS_PER_LINE)) & ~((1UL << CACHE_WORDS_PER_LINE) - 1);
STORE_EDGES_POINTER_LOOP:
    for (unsigned short ntb = 0; ntb < numTables; ntb++){
        hTables0[ntb].start_edges = start_addr;
        start_addr += (hTables0[ntb].n_edges >> (ROW_LOG - EDGE_LOG)) + 1;
        start_addr = (start_addr + (1UL << CACHE_WORDS_PER_LINE)) & ~((1UL << CACHE_WORDS_PER_LINE) - 1);
#ifndef __SYNTHESIS__
        assert(start_addr < HTB_SPACE);
#endif
        hTables1[ntb].start_offset = hTables0[ntb].start_offset;
        hTables1[ntb].start_edges = hTables0[ntb].start_edges;
        hTables1[ntb].n_edges = hTables0[ntb].n_edges;
    }

    writeEdges<
         T_DDR,
         EDGE_LOG,
         CNT_LOG,
         ROW_LOG,
         NODE_W,
         LAB_W,
         HASH1_W,
         HASH2_W,
         STREAM_D>(
            edge_buf,
            htb_buf,
            tDescriptors,
            hTables0,
            numTables,
            numDataEdges);

#ifndef __SYNTHESIS__
    end_addr = start_addr * (1UL << (ROW_LOG - 3)) + ((numTables * ((1 << HASH1_W) + 1)) << (BLOOM_LOG - 3));
    std::cout << "Occupied " << end_addr << " bytes, " << 
        end_addr / (float)(1UL << 20) << " MB. " << std::endl;
#endif

    /* 
    const size_t CNT_ROW = 1UL << (ROW_LOG - CNT_LOG);
    std::ofstream f("CHECK.txt");
    for(unsigned short tab = 0; tab < numTables; tab++){
        f << "Table " << tab << std::endl;
        ap_uint<64> start = hTables0[tab].start_offset;
        ap_uint<HASH1_W> counter = 0;
        for(ap_uint<64> addr = 0;
                addr < HTB_SIZE;
                addr++){
            T_DDR row = htb_buf[start + addr];
            for(int g = 0; g < CNT_ROW; g++){
                if (((addr << (ROW_LOG - CNT_LOG)) + g) % (1UL << HASH2_W) == 0){
                    f << counter << ": " << std::endl;
                    counter++;
                }

                ap_uint<(1UL << CNT_LOG)> edge = row.range(((g+1)<<CNT_LOG)-1, g<<CNT_LOG);
                f << "\t" << edge.range((1UL << CNT_LOG) - 1, 0) << std::endl;

            }
        }
    }
    f << std::endl;

    for(unsigned short tab = 0; tab < numTables; tab++){
        f << "Table " << tab << std::endl;
        ap_uint<64> start = hTables0[tab].start_edges;
        for(ap_uint<64> addr = 0;
                addr <= (hTables0[tab].n_edges >> (ROW_LOG-EDGE_LOG));
                addr++){
            T_DDR row = htb_buf[start + addr];
            for(int g = 0; g < EDGE_ROW; g++){
                ap_uint<(1UL << EDGE_LOG)> edge = row.range(((g+1)<<EDGE_LOG)-1, g<<EDGE_LOG);
                f << "\t" << edge.range(2*NODE_W-1, NODE_W) << " -> "
                    << edge.range(NODE_W-1, 0) << std::endl;
            }
        }
    }

    constexpr size_t BLOOM_W = (1UL << BLOOM_LOG);
    constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    float one_mean = 0;
    for (unsigned int tab = 0; tab < numTables; tab++){
        f << "Table " << tab << " indexed set" ;
        for(unsigned long addr = 0; addr < (1UL << HASH1_W); addr++){
            T_BLOOM row = bloom_p[addr + tab * (1UL << HASH1_W)];
            f << std::endl << "\t";
            for(int g = 0; g < K_FUN; g++){
                ap_uint<(1UL << (BLOOM_LOG - K_FUN_LOG)) > bloom_f = 
                    row.range(((g + 1) << (BLOOM_LOG - K_FUN_LOG)) - 1, g << (BLOOM_LOG - K_FUN_LOG));
                f << std::hex << (unsigned long)bloom_f << " ";
                unsigned int count = 0;
                while( bloom_f > 0) {
                    count++;
                    bloom_f = bloom_f & (bloom_f - 1);
                }
                one_mean += (float)count / ((1UL << HASH1_W + 1)*numTables);
            }
        }
        f << std::endl;
    }
    std::cout << one_mean << std::endl;

    for (unsigned int tab = 0; tab < numTables; tab++){
        f << "Table " << tab << " indexing set";
        T_BLOOM row = bloom_p[numTables * (1UL << HASH1_W) + tab];
        f << std::endl << "\t";
        for(int g = 0; g < K_FUN; g++){
            ap_uint<(1UL << (BLOOM_LOG - K_FUN_LOG)) > bloom_f = 
                row.range(((g + 1) << (BLOOM_LOG - K_FUN_LOG)) - 1, g << (BLOOM_LOG - K_FUN_LOG));
            f << std::hex << (unsigned long)bloom_f << " ";
        }
        f << std::endl;
    }
    f.close();
    */      
    
    /* for (unsigned int tab = 0; tab < numTables; tab++){ */
/* std::cout << "Table " << tab << " indexing set"; */
/* ap_uint<1024> row = bloom_indexing[tab]; */
/* std::cout << std::endl << "\t"; */
/* for(int g = 0; g < 16; g++){ */
/* unsigned long bloom_f = */
/* row.range(((g + 1)*64) - 1, g*64); */
/* std::cout << std::hex << (unsigned long)bloom_f << " "; */
/* } */
/* std::cout << std::endl; */
/* } */

#if DEBUG_STATS
        constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    for (unsigned int tab = 0; tab < numTables; tab++){
        for(unsigned long addr = 0; addr < (1UL << HASH1_W); addr++){
            T_BLOOM row = bloom_p[addr + tab * (1UL << HASH1_W)];
            for(int g = 0; g < K_FUN; g++){
                ap_uint<(1UL << (BLOOM_LOG - K_FUN_LOG)) > bloom_f = 
                    row.range(((g + 1) << (BLOOM_LOG - K_FUN_LOG)) - 1, g << (BLOOM_LOG - K_FUN_LOG));
                unsigned int count = 0;
                while( bloom_f > 0) {
                    count++;
                    bloom_f = bloom_f & (bloom_f - 1);
                }
                debug::bloom_fullness += (float)count / ((1UL << HASH1_W + 1)*numTables);
            }
        }
    }
#endif /* DEBUG_STATS */

}

template <typename T_DDR,
         typename T_BLOOM,
         size_t EDGE_LOG,
         size_t CNT_LOG,
         size_t BLOOM_LOG,
         size_t K_FUN_LOG,
         size_t ROW_LOG,
         size_t NODE_W,
         size_t LAB_W,
         size_t HASH1_W,
         size_t HASH2_W,
         size_t STREAM_D,
         size_t HTB_SPACE,
         size_t MAX_QV,
         size_t MAX_TB>
void preprocess(
        edge_t *edge_buf,
        T_DDR *htb_buf,
        T_DDR *bloom_p,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        AdjHT *hTables0,
        AdjHT *hTables1,
        unsigned short numQueryVert,
        unsigned short numQueryEdges,
        unsigned long numDataEdges)
{
    TableDescriptor tDescriptors[MAX_TB];
    unsigned short numTables = 0;

    buildTableDescriptors<
        MAX_QV,
        MAX_TB,
        NODE_W,
        LAB_W>(
            &edge_buf[numDataEdges],
            qVertices0,
            qVertices1,
            tDescriptors,
            numTables,
            numQueryVert,
            numQueryEdges);

    fillTables<
         T_DDR,
         T_BLOOM,
         EDGE_LOG,
         CNT_LOG,
         BLOOM_LOG,
         K_FUN_LOG,
         ROW_LOG,
         NODE_W,
         LAB_W,
         HASH1_W,
         HASH2_W,
         STREAM_D,
         HTB_SPACE>(
            edge_buf,
            htb_buf,
            bloom_p,
            hTables0,
            hTables1,
            tDescriptors,
            numDataEdges,
            numTables);
}
