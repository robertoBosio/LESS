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

#ifndef __SYNTHESIS__
#define DEBUG_STATS
#include "debug.hpp"
#endif /*__SYNTHESIS__*/


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
            qVertices0[nodeSrcPos].addTableIndexing(g);
            qVertices1[nodeSrcPos].addTableIndexing(g);
            qVertices0[nodeDstPos].addTableIndexed(g, nodeSrcPos);
            qVertices1[nodeDstPos].addTableIndexed(g, nodeSrcPos);
        } else {
            qVertices0[nodeSrcPos].addTableIndexed(g, nodeDstPos);
            qVertices0[nodeDstPos].addTableIndexing(g);
            qVertices1[nodeSrcPos].addTableIndexed(g, nodeDstPos);
            qVertices1[nodeDstPos].addTableIndexing(g);
        }
    }
}


/* Translating counter addresses to ram addresses and
 * increasing of one the counter specified. */
template <typename T_DDR,
         size_t EDGE_LOG,
         size_t CNT_LOG,
         size_t ROW_LOG,
         size_t NODE_W,
         size_t HASH1_W,
         size_t HASH2_W>
void increaseCounter(
        hls::stream<ap_uint<(1UL << EDGE_LOG)>> &stream_edge,
        hls::stream<ap_uint<HASH1_W>> &stream_index1,
        hls::stream<ap_uint<HASH2_W>> &stream_index2,
        hls::stream<unsigned short> &stream_ntable,
        hls::stream<bool> &stream_end,
        AdjHT *hTables,
        T_DDR *htb_buf)
{
    const size_t CNT_ROW = 1UL << (ROW_LOG - CNT_LOG);
    T_DDR ram_row, mask;
    ap_uint<HASH1_W> index1;
    ap_uint<HASH2_W> index2;
    ap_uint<(1UL << EDGE_LOG)> edge;
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
        edge = stream_edge.read();

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
        ap_uint<(1UL << CNT_LOG)> prev_addr = 0;
        ap_uint<(1UL << CNT_LOG)> counter;
        unsigned int hash_used = 0;

COUNTER_TO_OFFSET_TABLE_LOOP:
        for(unsigned int start = 0; start < HTB_SIZE; start++){
            row = htb_buf[start + hTables[ntb].start_offset];

COUNTER_TO_OFFSET_ROW_LOOP:
            for (int g = 0; g < CNT_ROW; g++){
                counter = row.range((1UL << CNT_LOG) - 1, 0);
                row_new >>= (1UL << CNT_LOG);
                row >>= (1UL << CNT_LOG);

#ifdef DEBUG_STATS
                if (counter > debug::max_collisions) 
                    debug::max_collisions = counter;
                debug::avg_collisions += (float)counter / (1 << (HASH1_W + HASH2_W));
#endif /* DEBUG_STATS */

                row_new.range((CNT_ROW << CNT_LOG) - 1, (CNT_ROW - 1) << CNT_LOG) = base_addr;
                base_addr += counter;

                /* Check if an hash is used */
                if (((start << (ROW_LOG - CNT_LOG)) + g) % (1UL << HASH2_W) == 0){
                    if(prev_addr < base_addr){
                        hash_used++;
                        prev_addr = base_addr;
                    }
                }
            }
            htb_buf[start + hTables[ntb].start_offset] = row_new;
        }
        hTables[ntb].hash_set = hash_used;
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

template <size_t EDGE_LOG,
         size_t NODE_W,
         size_t LAB_W,
         size_t HASH1_W,
         size_t HASH2_W>
void edgeToHash(
        edge_t *edge_buf,
        TableDescriptor *tDescriptors,
        unsigned short numTables,
        unsigned int numDataEdges,

        hls::stream<ap_uint<(1UL << EDGE_LOG)>> &stream_edge,
        hls::stream<ap_uint<HASH1_W>> &stream_index1,
        hls::stream<ap_uint<HASH2_W>> &stream_index2,
        hls::stream<unsigned short> &stream_ntable,
        hls::stream<bool> &stream_end_out)
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

                stream_index1.write(indexAdj);
                stream_index2.write(indexEdge);
                stream_ntable.write(g);
                stream_edge.write(vertexIndexing.concat(vertexIndexed));
                stream_end_out.write(false);
            }

#ifdef UNDIRECTED   
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

                stream_index1.write(indexAdj);
                stream_index2.write(indexEdge);
                stream_ntable.write(g);
                stream_edge.write(vertexIndexing.concat(vertexIndexed));
                stream_end_out.write(false);
            }
#endif

        }
    }
    stream_end_out.write(true);
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
void countEdges(
        edge_t *edge_buf,
        TableDescriptor *tDescriptors,
        AdjHT *hTables,
        T_DDR *htb_buf,
        unsigned short numTables,
        unsigned long numDataEdges)
{

#pragma HLS DATAFLOW

    hls::stream<ap_uint<(1UL << EDGE_LOG)>, STREAM_D> stream_edge("Offset edge");
    hls::stream<ap_uint<HASH1_W>, STREAM_D> stream_c_index1("Counters index 1");
    hls::stream<ap_uint<HASH2_W>, STREAM_D> stream_c_index2("Counters index 2");
    hls::stream<unsigned short, STREAM_D> stream_c_ntable("Counters table");
    hls::stream<bool> stream_end_count("Counters end");

    /* Count edges per vertex source */
    edgeToHash<
        EDGE_LOG,
        NODE_W,
        LAB_W,
        HASH1_W,
        HASH2_W>(
            edge_buf,
            tDescriptors,
            numTables,
            numDataEdges,
            stream_edge,
            stream_c_index1,
            stream_c_index2,
            stream_c_ntable,
            stream_end_count);

    /* Update specific counter */
    increaseCounter<
        T_DDR,
        EDGE_LOG,
        CNT_LOG,
        ROW_LOG,
        NODE_W,
        HASH1_W,
        HASH2_W>(
            stream_edge,
            stream_c_index1,
            stream_c_index2,
            stream_c_ntable,
            stream_end_count,
            hTables,
            htb_buf);

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
        edge_t *edge_buf,
        TableDescriptor *tDescriptors,
        AdjHT *hTables,
        T_DDR *htb_buf,
        unsigned short numTables,
        unsigned long numDataEdges)
{

#pragma HLS DATAFLOW

    hls::stream<ap_uint<(1UL << EDGE_LOG)>, STREAM_D> stream_edge("Offset edge");
    hls::stream<ap_uint<HASH1_W>, STREAM_D> stream_o_index1("Offsets index 1");
    hls::stream<ap_uint<HASH2_W>, STREAM_D> stream_o_index2("Offsets index 2");
    hls::stream<unsigned short, STREAM_D> stream_o_ntable("Offset table");
    hls::stream<bool> stream_end_offset("Offset end");

    /* Write edges based on precomputed offsets */
    edgeToHash<
        EDGE_LOG,
        NODE_W,
        LAB_W,
        HASH1_W,
        HASH2_W>(
            edge_buf,
            tDescriptors,
            numTables,
            numDataEdges,
            stream_edge,
            stream_o_index1,
            stream_o_index2,
            stream_o_ntable,
            stream_end_offset);

    storeEdges<
        T_DDR,
        EDGE_LOG,
        CNT_LOG,
        ROW_LOG,
        NODE_W,
        HASH1_W,
        HASH2_W>(
            hTables,
            htb_buf,
            stream_edge,
            stream_o_index1,
            stream_o_index2,
            stream_o_ntable,
            stream_end_offset);
}

/* Reads two times the data graph and fills the data stuctures */
template <typename T_DDR,
         size_t EDGE_LOG,
         size_t CNT_LOG,
         size_t ROW_LOG,
         size_t NODE_W,
         size_t LAB_W,
         size_t HASH1_W,
         size_t HASH2_W,
         size_t STREAM_D,
         size_t HTB_SPACE>
void fillTables(
        edge_t *edge_buf,
        AdjHT *hTables0,
        AdjHT *hTables1,
        T_DDR *htb_buf,
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
         EDGE_LOG,
         CNT_LOG,
         ROW_LOG,
         NODE_W,
         LAB_W,
         HASH1_W,
         HASH2_W,
         STREAM_D>(
            edge_buf,
            tDescriptors,
            hTables0,
            htb_buf,
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
    
    start_addr = end_addr;
STORE_EDGES_POINTER_LOOP:
    for (unsigned short ntb = 0; ntb < numTables; ntb++){
        hTables0[ntb].start_edges = start_addr;
        start_addr += (hTables0[ntb].n_edges >> (ROW_LOG - EDGE_LOG)) + 1;
#ifndef __SYNTHESIS__
        assert(start_addr < HTB_SPACE);
#endif
        /* Fill the last row of tables's edges of 1s.
         * Because while reading the minimum set it is
         * read the entire line, and must be possible
         * to recognize which are real edges and which no */
        hTables1[ntb].start_offset = hTables0[ntb].start_offset;
        hTables1[ntb].start_edges = hTables0[ntb].start_edges;
        hTables1[ntb].n_edges = hTables0[ntb].n_edges;
        hTables1[ntb].hash_set = hTables0[ntb].hash_set;
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
            tDescriptors,
            hTables0,
            htb_buf,
            numTables,
            numDataEdges);

#ifndef __SYNTHESIS__
    end_addr = start_addr;
    std::cout << "Occupied " << end_addr * 64 << " bytes, " << end_addr <<
        " words." << std::endl;
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
    f.close();
    */      
}

template <typename T_DDR,
         size_t EDGE_LOG,
         size_t CNT_LOG,
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
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        TableDescriptor *tDescriptors,
        AdjHT *hTables0,
        AdjHT *hTables1,
        unsigned short &numTables,
        unsigned short numQueryVert,
        unsigned short numQueryEdges,
        unsigned long numDataEdges
        ){

    buildTableDescriptors<
        MAX_QV,
        MAX_TB,
        NODE_W,
        LAB_W>(
            edge_buf,
            qVertices0,
            qVertices1,
            tDescriptors,
            numTables,
            numQueryVert,
            numQueryEdges);

    fillTables<
         T_DDR,
         EDGE_LOG,
         CNT_LOG,
         ROW_LOG,
         NODE_W,
         LAB_W,
         HASH1_W,
         HASH2_W,
         STREAM_D,
         HTB_SPACE>(
            &edge_buf[numQueryVert + numQueryEdges],
            hTables0,
            hTables1,
            htb_buf,
            tDescriptors,
            numDataEdges,
            numTables);
}
