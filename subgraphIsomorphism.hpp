#include <ap_int.h>
#include <hls_stream.h>

#include "hls_print.h"
#include "Parameters.hpp"
#include "QueryVertex.hpp"
#include "Trie.hpp"
#include "hash_lookup3.hpp"

#include <unordered_map>

#ifndef __SYNTHESIS__
#include <cassert>
#include <fstream>
#endif

#define HTB_SIZE    (1UL << (H_W_1 + H_W_2 - (9 - C_W)))
#define CNT_ROW     (1UL << (9 - C_W))
#define EDGE_ROW    (1UL << (9 - E_W))

#define V_ID_W      VERTEX_WIDTH_BIT
#define V_L_W       LABEL_WIDTH
#define MAX_QV      MAX_QUERY_VERTICES
#define MAX_CL      MAX_COLLISIONS
#define MAX_TB      MAX_TABLES
#define H_W_1       HASH_WIDTH_FIRST
#define H_W_2       HASH_WIDTH_SECOND
#define C_W         COUNTER_WIDTH
#define E_W         EDGE_WIDTH
#define S_D         STREAM_DEPTH    
#define S_D_RES     STREAM_DEPTH_RES    

/* Partial solutions still to be checked.
 * Set to 1 to start the process. */
long int nPartialSol = 1;

#ifndef __SYNTHESIS__
/* statistic purpose */
long int nStreamPar = 0;
long int nStreamParMax = 0;
#endif

/* Builds the table descriptors based on the information
 * from the query graph. */
void buildTableDescriptors(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        TableDescriptor *tDescriptors,
        ap_uint<8> &numTables,
        ap_uint<8> &numQueryVert)
{
    /* Translate from id of vertex to position in the order */
    ap_uint<8> fromNumToPos[MAX_QV];

    /* Filling information about query vertices and coping
     * the vertex order needed by multiway join */
    bool last = stream_end.read();
FILL_ORDER_LOOP:
    while(!last){
        ap_uint<V_ID_W> node = stream_src.read();
        fromNumToPos[node] = numQueryVert;
        numQueryVert++;
        last = stream_end.read();
    }
    numQueryVert--;

    /* Creating table descriptors */
    last = stream_end.read();
CREATE_TABDESC_LOOP:
    while(!last){
        bool dirEdge = false;
        ap_uint<V_ID_W> nodesrc = stream_src.read();
        ap_uint<V_ID_W> nodedst = stream_dst.read();
        ap_uint<V_L_W> labelsrc = stream_src_l.read();
        ap_uint<V_L_W> labeldst = stream_dst_l.read();
        ap_uint<8> nodeSrcPos = fromNumToPos[nodesrc];
        ap_uint<8> nodeDstPos = fromNumToPos[nodedst];

        if (nodeSrcPos < nodeDstPos)
            dirEdge = true;

#ifndef __SYNTHESIS__
        std::cout << (unsigned int)nodesrc << "(" << (char)labelsrc << ")"
            << " -> " <<  (unsigned int)nodedst << "(" << (char)labeldst << ")" 
            << std::endl;
#endif

        /* Understanding if table already exists */
        uint8_t g = 0;
FIND_CORRECT_TABLE_LOOP:
        for (; g < numTables; g++){
            if (tDescriptors[g].src_label == labelsrc &&
                    tDescriptors[g].dst_label == labeldst &&
                    tDescriptors[g].dir == dirEdge)
                break;
        }
     
        /* Add new table descriptor */   
        if (g == numTables){
            tDescriptors[g].src_label = labelsrc;
            tDescriptors[g].dst_label = labeldst;
            tDescriptors[g].dir = dirEdge;
            numTables++;
        }

#ifndef __SYNTHESIS__
        if (dirEdge){
        std::cout << "Table " << (int)g << ": " << (char)labelsrc << 
             " -> " << (char)labeldst << std::endl;
        } else {
        std::cout << "Table " << (int)g << ": " << (char)labeldst << 
             " <- " << (char)labelsrc << std::endl;
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

        last = stream_end.read();
    }
}

/* Translating counter addresses to ram addresses and
 * increasing of one the counter specified.
 * The counter is on 16 bits, so shift by 4. */
void increaseCounter(
        hls::stream<ap_uint<(1UL << E_W)>> &stream_edge,
        hls::stream<ap_uint<H_W_1>> &stream_index1,
        hls::stream<ap_uint<H_W_2>> &stream_index2,
        hls::stream<ap_uint<8>> &stream_ntable,
        hls::stream<bool> &stream_end,
        AdjHT *hTables,
        ap_uint<512> *htb_buf)
{
    ap_uint<512> ram_row;
    ap_uint<H_W_1> index1;
    ap_uint<H_W_2> index2;
    ap_uint<(1UL << E_W)> edge;
    ap_uint<8> ntb;
    ap_uint<64> addr_outrow;
    ap_uint<64> addr_inrow;
    ap_uint<64> addr_counter;
    ap_uint<(1UL << C_W)> counter;

    bool last = stream_end.read();
INCREASE_COUNTER_DDR_LOOP:
    while(!last){
        index1 = stream_index1.read();
        index2 = stream_index2.read();
        ntb = stream_ntable.read();
        edge = stream_edge.read();

        /* Address in the matrix [H_W_1][H_W_2] */
        addr_counter = index1;
        addr_counter <<= H_W_2;
        addr_counter += index2;

        /* Compute address of row storing the counter */
        addr_outrow = hTables[ntb].start_offset + (addr_counter >> (9 - C_W));

        /* Compute address of counter inside the row */
        addr_inrow = addr_counter.range((9 - C_W) - 1, 0);

        /* Read, modify and write the counter */
        ram_row = htb_buf[addr_outrow];
        counter = ram_row.range(((addr_inrow + 1) << C_W)-1, addr_inrow << C_W); 
        counter += 1;
        ram_row.range(((addr_inrow + 1) << C_W)-1, addr_inrow << C_W) = counter;
        htb_buf[addr_outrow] = ram_row;
 
        hTables[ntb].n_edges++;     
        last = stream_end.read();
    }
}

/* Transform counters to offsets. */
void counterToOffset(
        ap_uint<8> numTables,
        AdjHT *hTables,
        ap_uint<512> *htb_buf)
{
    ap_uint<512> ram_row;

COUNTER_TO_OFFSET_DDR_LOOP:
    for (ap_uint<8> ntb = 0; ntb < numTables; ntb++){
        ap_uint<(1 << C_W)> base_addr = 0;
        ap_uint<(1UL << C_W)> prev_addr = 0;
        ap_uint<(1UL << C_W)> counter;
        ap_uint<32> hash_used = 0;

COUNTER_TO_OFFSET_TABLE_LOOP:
        for(ap_uint<64> start = 0; start < HTB_SIZE; start++){
            ram_row = htb_buf[start + hTables[ntb].start_offset];
            for (int g = 0; g < CNT_ROW; g++){
                counter = ram_row.range(((g + 1) << C_W) - 1, g << C_W);
                ram_row.range(((g + 1) << C_W) - 1, g << C_W) = base_addr;
                base_addr += counter;

                /* Store used bit */
                if (counter > 0){
                    ram_row.set(((g + 1) << C_W) - 1);
                }

                /* Check if an hash is used */
                if (((start << (9 - C_W)) + g) % (1UL << H_W_2) == 0){
                    if(prev_addr < base_addr){
                        hash_used++;
                        prev_addr = base_addr;
                    }
                }
            }
            htb_buf[start + hTables[ntb].start_offset] = ram_row;
        }
        hTables[ntb].hash_set = hash_used;
    }       
}

/* Store edges based on offsets */
void storeEdges(
        hls::stream<ap_uint<(1UL << E_W)>> &stream_edge,
        hls::stream<ap_uint<H_W_1>> &stream_index1,
        hls::stream<ap_uint<H_W_2>> &stream_index2,
        hls::stream<ap_uint<8>> &stream_ntable,
        hls::stream<bool> &stream_end,
        AdjHT *hTables,
        ap_uint<512> *htb_buf)
{
    ap_uint<512> ram_row;
    ap_uint<H_W_1> index1;
    ap_uint<H_W_2> index2;
    ap_uint<(1UL << E_W)> edge;
    ap_uint<8> ntb;
    ap_uint<64> addr_outrow;
    ap_uint<64> addr_inrow;
    ap_uint<64> addr_offset;
    ap_uint<(1UL << C_W) - 1> offset;

    bool last = stream_end.read();
INCREASE_COUNTER_DDR_LOOP:
    while(!last){
        index1 = stream_index1.read();
        index2 = stream_index2.read();
        ntb = stream_ntable.read();
        edge = stream_edge.read();

        /* Address in the matrix [H_W_1][H_W_2] */
        addr_offset = index1;
        addr_offset <<= H_W_2;
        addr_offset += index2;

        /* Compute address of row storing the offset */
        addr_outrow = hTables[ntb].start_offset + (addr_offset >> (9 - C_W));

        /* Compute address of offset inside the row */
        addr_inrow = addr_offset.range((9 - C_W) - 1, 0);

        /* Read, modify and write the offset */
        ram_row = htb_buf[addr_outrow];
        offset = ram_row.range(((addr_inrow + 1) << C_W) - 2, addr_inrow << C_W); 
        ram_row.range(((addr_inrow + 1) << C_W) - 2, 
                addr_inrow << C_W) = offset + 1;
        htb_buf[addr_outrow] = ram_row;
 
        /* Compute address of row that will store the edge */
        addr_outrow = hTables[ntb].start_edges 
            + (offset >> (9 - E_W));

        /* Compute address of the edge inside the row */
        addr_inrow = offset.range((9 - E_W) - 1, 0);

        /* Read, modify and write the edge */
        ram_row = htb_buf[addr_outrow];
        ram_row.range(((addr_inrow + 1) << E_W) - 1,
            addr_inrow << E_W) = edge;
        htb_buf[addr_outrow] = ram_row;

        last = stream_end.read();
    }
}

void edgeToHash(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end_in,
        TableDescriptor *tDescriptors,
        ap_uint<8> numTables,

        hls::stream<ap_uint<(1UL << E_W)>> &stream_edge,
        hls::stream<ap_uint<H_W_1>> &stream_index1,
        hls::stream<ap_uint<H_W_2>> &stream_index2,
        hls::stream<ap_uint<8>> &stream_ntable,
        hls::stream<bool> &stream_end_out)
{
    
    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    
    bool last = stream_end_in.read();
COUNT_OCCURENCIES_LOOP:
    while(!last){
        ap_uint<V_ID_W> nodesrc = stream_src.read();
        ap_uint<V_ID_W> nodedst = stream_dst.read();
        ap_uint<V_L_W> labelsrc = stream_src_l.read();
        ap_uint<V_L_W> labeldst = stream_dst_l.read();
        
        /* Finding correct table */
        ap_uint<8> g = 0;
        for (; g < numTables; g++){
            if (tDescriptors[g].src_label == labelsrc &&
                    tDescriptors[g].dst_label == labeldst){

                ap_uint<V_ID_W> vertexIndexing, vertexIndexed;               
                if(tDescriptors[g].dir){
                    vertexIndexing = nodesrc;
                    vertexIndexed = nodedst;
                } else {
                    vertexIndexing = nodedst;
                    vertexIndexed = nodesrc;
                }
                
                /* Compute indices for hash table */
                stream_hash_in.write(vertexIndexing);
                xf::database::hashLookup3<V_ID_W>(
                        stream_hash_in,
                        stream_hash_out);

                ap_uint<H_W_1> indexAdj = 
                    stream_hash_out.read().range(H_W_1 - 1, 0); 

                stream_hash_in.write(vertexIndexed);
                xf::database::hashLookup3<V_ID_W>(
                        stream_hash_in, 
                        stream_hash_out);

                ap_uint<H_W_2> indexEdge = 
                    stream_hash_out.read().range(H_W_2 - 1, 0); 

                stream_index1.write(indexAdj);
                stream_index2.write(indexEdge);
                stream_ntable.write(g);
                stream_edge.write(vertexIndexing.concat(vertexIndexed));
                stream_end_out.write(false);
            }
        }
        last = stream_end_in.read();
    }
    stream_end_out.write(true);
}

void countEdges(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end_in,
        TableDescriptor *tDescriptors,
        AdjHT *hTables,
        ap_uint<512> *htb_buf,
        ap_uint<8> numTables)
{

#pragma HLS DATAFLOW

    hls::stream<ap_uint<2*V_ID_W>, S_D> stream_edge("Offset edge");
    hls::stream<ap_uint<H_W_1>, S_D> stream_c_index1("Counters index 1");
    hls::stream<ap_uint<H_W_2>, S_D> stream_c_index2("Counters index 2");
    hls::stream<ap_uint<8>, S_D> stream_c_ntable("Counters table");
    hls::stream<bool> stream_end_count("Counters end");

    /* Count edges per vertex source */
    edgeToHash(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end_in,
            tDescriptors,
            numTables,
            stream_edge,
            stream_c_index1,
            stream_c_index2,
            stream_c_ntable,
            stream_end_count);

    /* Update specific counter */
    increaseCounter(
            stream_edge,
            stream_c_index1,
            stream_c_index2,
            stream_c_ntable,
            stream_end_count,
            hTables,
            htb_buf);

}

void writeEdges(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end_in,
        TableDescriptor *tDescriptors,
        AdjHT *hTables,
        ap_uint<512> *htb_buf,
        ap_uint<8> numTables)
{

#pragma HLS DATAFLOW

    hls::stream<ap_uint<2*V_ID_W>, S_D> stream_edge("Offset edge");
    hls::stream<ap_uint<H_W_1>, S_D> stream_o_index1("Offsets index 1");
    hls::stream<ap_uint<H_W_2>, S_D> stream_o_index2("Offsets index 2");
    hls::stream<ap_uint<8>, S_D> stream_o_ntable("Offset table");
    hls::stream<bool> stream_end_offset("Offset end");

    /* Write edges based on precomputed offsets */
    edgeToHash(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end_in,
            tDescriptors,
            numTables,
            stream_edge,
            stream_o_index1,
            stream_o_index2,
            stream_o_ntable,
            stream_end_offset);

    storeEdges(
            stream_edge,
            stream_o_index1,
            stream_o_index2,
            stream_o_ntable,
            stream_end_offset,
            hTables,
            htb_buf);
}

/* Reads two times the data graph and fills the data stuctures */
void fillTables(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end,
        AdjHT *hTables0,
        AdjHT *hTables1,
        ap_uint<512> *htb_buf,
        TableDescriptor *tDescriptors,
        ap_uint<8> numTables)
{

    /* Resetting portion of memory dedicated to counters 
     * 1 << H_W_1 * H_W_2 is the number of counters needed
     * for each table, then it should be divided by the number
     * of counters stored in each row which is 1 << (9 - C_W)*/
    ap_uint<64> end_addr = numTables * HTB_SIZE;
/*
RESET_HASHTABLES_LOOP:
    for (ap_uint<64> addr = 0; addr < end_addr; addr++){
        htb_buf[addr] = 0;
    }
*/
    ap_uint<64> start_addr = 0;
STORE_HASHTABLES_POINTER_LOOP:
    for (ap_uint<8> ntb = 0; ntb < numTables; ntb++){
        hTables0[ntb].start_offset = start_addr;
        hTables0[ntb].n_edges = 0;
        start_addr += HTB_SIZE;
    }

    countEdges(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end,
            tDescriptors,
            hTables0,
            htb_buf,
            numTables);
    
    /* From counts to offsets */
    counterToOffset(
            numTables,
            hTables0,
            htb_buf);
    
    start_addr = end_addr;
STORE_EDGES_POINTER_LOOP:
    for (ap_uint<8> ntb = 0; ntb < numTables; ntb++){
        hTables0[ntb].start_edges = start_addr;
        start_addr += (hTables0[ntb].n_edges >> (9 - E_W)) + 1;
        /* Fill the last row of tables's edges of 1s.
         * Because while reading the minimum set it is
         * read the entire line, and must be possible
         * to recognize which are real edges and which no */
        hTables1[ntb].start_offset = hTables0[ntb].start_offset;
        hTables1[ntb].start_edges = hTables0[ntb].start_edges;
        hTables1[ntb].n_edges = hTables0[ntb].n_edges;
        hTables1[ntb].hash_set = hTables0[ntb].hash_set;
        //htb_buf[start_addr - 1] = ~((ap_uint<512>)0);
    }

    writeEdges(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end,
            tDescriptors,
            hTables0,
            htb_buf,
            numTables);

#ifndef __SYNTHESIS__
    end_addr = start_addr;
    std::cout << "Occupied " << end_addr * 64 << " bytes, " << end_addr <<
        " words." << std::endl;
#endif

//   std::ofstream f("CHECK.txt");
//   for(ap_uint<8> tab = 0; tab < numTables; tab++){
//       f << "Table " << tab << std::endl;
//       ap_uint<64> start = hTables0[tab].start_offset;
//       ap_uint<H_W_1> counter = 0;
//       for(ap_uint<64> addr = 0; 
//               addr < HTB_SIZE; 
//               addr++){
//           ap_uint<512> row = htb_buf[start + addr];
//           for(int g = 0; g < CNT_ROW; g++){
//               if (((addr << (9 - C_W)) + g) % (1UL << H_W_2) == 0){
//                   f << counter << ": " << std::endl;
//                   counter++;
//               }
//
//               ap_uint<(1UL << C_W)> edge = row.range(((g+1)<<C_W)-1, g<<C_W);
//               f << "\t" << edge.range((1UL << C_W)-2, 0) << " -> "
//                  << edge.test((1UL << C_W)-1) << std::endl;
//               
//           }
//       }
//   } 
//   f << std::endl;
//   for(ap_uint<8> tab = 0; tab < numTables; tab++){
//       f << "Table " << tab << std::endl;
//       ap_uint<64> start = hTables0[tab].start_edges;
//       for(ap_uint<64> addr = 0; 
//               addr <= (hTables0[tab].n_edges >> (9-E_W)); 
//               addr++){
//           ap_uint<512> row = htb_buf[start + addr];
//           for(int g = 0; g < EDGE_ROW; g++){
//               ap_uint<(1UL << E_W)> edge = row.range(((g+1)<<E_W)-1, g<<E_W);
//               f << "\t" << edge.range(2*V_ID_W-1, V_ID_W) << " -> "
//                  << edge.range(V_ID_W-1, 0) << std::endl;
//           }
//       }
//   } 
//   f.close();
}

template <uint8_t W_1,
        uint8_t W_2,
        uint8_t SHF,
        uint8_t T>
ap_uint<(1UL << T)> read_table(
        ap_uint<W_1> index1,
        ap_uint<W_2> index2,
        ap_uint<512> *htb_buf,
        ap_uint<32> start_addr)
{
    //static ap_uint<512> ram_row;
    ap_uint<512> ram_row;
    /* Impossible that the first word read is the all 1s
     * since before edges are read offsets which are on 
     * the lower address space of the memory */
    //static ap_uint<32> prev_addr_row = ~((ap_uint<32>)0);
    ap_uint<32> addr_row;
    ap_uint<32> addr_inrow;
    ap_uint<32> addr_counter;

    addr_counter = index1;
    addr_counter <<= SHF;
    addr_counter += index2;

    /* Compute address of row storing the counter */
    addr_row = start_addr + (addr_counter >> (9 - T));

    /* Compute address of data inside the row */
    addr_inrow = addr_counter.range((9 - T) - 1, 0);

    /* Read the data */
    //if (prev_addr_row != addr_row){
        ram_row = htb_buf[addr_row];
        //prev_addr_row = addr_row;
    //}

    return ram_row.range(((addr_inrow + 1) << T) - 1,
            addr_inrow << T);

}

/* Multiway join propose: retrieves the tables in which
 * the current query vertex is involved, computes the sizes
 * and keep track of the smallest one */
void mwj_propose_findmin(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<512> *htb_buf,
        
        hls::stream<ap_uint<32>> &stream_minread_out,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_out,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    ap_uint<8> tableIndex = 0;
    ap_uint<8> curQV = 0;
    ap_uint<32> minSize = (1UL << 32) - 1;
    ap_uint<32> minStart, minEnd, minOff;
    ap_uint<8 + 1 + V_ID_W> minData;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    static bool start = true;
    bool last;

    if (start){
    	start = false;
    	last = true;
    } else {
    	last = stream_end_embed_in.read();
    }
    
PROPOSE_COPYING_EMBEDDING_LOOP:
    while(!last){
        curEmb[curQV] = stream_embed_in.read();
        stream_embed_out.write(curEmb[curQV]);
        stream_end_embed_out.write(false);
        curQV++;
        last = stream_end_embed_in.read();
    }
    stream_end_embed_out.write(true);

    //hls::print("Proposing for %d\n", (int)curQV);
    /* Find sizes of sets indexed by the current query vertex */
PROPOSE_TBINDEXING_LOOP:
    for(int g = 0; g < qVertices[curQV].numTablesIndexing; g++){
        tableIndex = qVertices[curQV].tables_indexing[g];
        
        //std::cout << "\tEvaluating " << tableIndex << " : size " <<
        //   hTables[tableIndex].hash_set << std::endl; 
        if (hTables[tableIndex].hash_set < minSize){
            minSize = hTables[tableIndex].hash_set;
            minOff = hTables[tableIndex].start_offset;
            minStart = 0;
            minEnd = HTB_SIZE;
            minData.range(7, 0) = tableIndex;
            minData.clear(8);
            minData.range(V_ID_W + 8, 9) = 0;
        }
    }

    /* Find sizes of sets in which the current query vertex
     * is indexed by an other query vertex */
PROPOSE_TBINDEXED_LOOP:
    for(int g = 0; g < qVertices[curQV].numTablesIndexed; g++){
        ap_uint<(1UL << C_W)> start_off = 0;
        ap_uint<(1UL << C_W)> end_off;
        tableIndex = qVertices[curQV].tables_indexed[g];
        uint8_t ivPos = qVertices[curQV].vertex_indexing[g];

        stream_hash_in.write(curEmb[ivPos]);
        xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
        ap_uint<H_W_1> index = stream_hash_out.read().range(H_W_1 - 1, 0);

        if (index != 0){ 
            start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                    index - 1,
                    (1UL << H_W_2) - 1,
                    htb_buf,
                    hTables[tableIndex].start_offset);
            start_off = start_off.range((1UL << C_W) - 2, 0);
        }

        end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                index,
                (1UL << H_W_2) - 1,
                htb_buf,
                hTables[tableIndex].start_offset);
        end_off = end_off.range((1UL << C_W) - 2, 0);

       // std::cout << "\tEvaluating " << tableIndex << " : size " <<
       //    end_off << " - " << start_off << std::endl; 
        if ((end_off - start_off) < minSize) {
            minSize = end_off - start_off;
            minOff = hTables[tableIndex].start_edges;
            minStart = start_off;
            minEnd = end_off;
            minData.range(7, 0) = tableIndex;
            minData.set(8);
            minData.range(V_ID_W + 8, 9) = curEmb[ivPos];
        }
    }
    stream_minread_out.write(minStart);
    stream_minread_out.write(minEnd);
    stream_minread_out.write(minOff);
    stream_min_out.write(minData);
}

/* Read the minimum set from memory */
void mwj_propose_readmin(
        hls::stream<ap_uint<32>> &stream_minread_in,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_in,
        ap_uint<512> *htb_buf,

        hls::stream<ap_uint<H_W_1>> &stream_min_out,
        hls::stream<bool> &stream_end_out,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_set)
{

    ap_uint<32> minStart = stream_minread_in.read();
    ap_uint<32> minEnd = stream_minread_in.read();
    ap_uint<32> minOff = stream_minread_in.read();
    ap_uint<V_ID_W + 1 + 8> minData = stream_min_in.read();
    stream_min_set.write(minData);

    //std::cout << "\tMinimum set is " << minData.range(7, 0) << std::endl;
    //std::cout << "\tReading from " << minStart << " to " << minEnd 
    //    << " with offset " << minOff << std::endl;
    if (minData.test(8)){
        hls::stream<ap_uint<V_ID_W>> stream_hash_in;
        hls::stream<ap_uint<64>> stream_hash_out;
        ap_uint<V_ID_W*2> edge;
        ap_uint<V_ID_W> vertex, vertexCheck;

PROPOSE_READ_MIN_INDEXED_LOOP:
        for (ap_uint<32> i = minStart; i < minEnd; i++){
            edge = read_table<32, 1, 0, E_W>(
                    i,
                    0,
                    htb_buf,
                    minOff);

            vertexCheck = edge.range(2*V_ID_W - 1, V_ID_W);
            vertex = edge.range(V_ID_W - 1, 0);

            if (minData.range(V_ID_W + 8, 9) == vertexCheck){
                stream_hash_in.write(vertex);
                xf::database::hashLookup3<V_ID_W>(
                        stream_hash_in, 
                        stream_hash_out);
                ap_uint<H_W_1> vertexHash =
                    stream_hash_out.read().range(H_W_1 - 1, 0);
                stream_min_out.write(vertexHash);
                stream_end_out.write(false);
            }
        }
    } else {
        ap_uint<(1UL << C_W)> start_off = 0;
        ap_uint<(1UL << C_W)> end_off;

PROPOSE_READ_MIN_INDEXING_LOOP:
        for (int i = 0; i < (1UL << H_W_1); i++){
            end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                    i,
                    (1UL << H_W_2) - 1,
                    htb_buf,
                    minOff);
            end_off = end_off.range((1UL << C_W) - 2, 0);
            
            /* Compare offsets and stream if they are different, 
             * in other words if some stored edge has that hash */
            if (end_off > start_off){
                stream_min_out.write((ap_uint<H_W_1>)i);
                stream_end_out.write(false);
                start_off = end_off;
            }
        }
    }
    stream_end_out.write(true);
}

void mwj_propose(
        AdjHT *hTables,
        QueryVertex *qVertices,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        ap_uint<512> *htb_buf0,
        ap_uint<512> *htb_buf1,
        
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<H_W_1>> &stream_min_out,
        hls::stream<bool> &stream_end_min_out,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_set)
{

#pragma HLS DATAFLOW

    //std::cout << "Propose: " << std::endl;
    hls::stream<ap_uint<32>, S_D> stream_minread;
    hls::stream<ap_uint<V_ID_W + 1 + 8>, S_D> stream_min;
    
    mwj_propose_findmin(
            stream_embed_in,
            stream_end_embed_in,
            hTables,
            qVertices,
            htb_buf0,
            stream_minread,
            stream_min,
            stream_embed_out,
            stream_end_embed_out);

    mwj_propose_readmin(
            stream_minread,
            stream_min,
            htb_buf1,
            stream_min_out,
            stream_end_min_out,
            stream_min_set);

}

void mwj_intersect(
        AdjHT *hTables,
        QueryVertex *qVertices,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<H_W_1>> &stream_min_in,
        hls::stream<bool> &stream_end_min_in,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_desc_in,
        ap_uint<512> *htb_buf,

        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<H_W_1>> &stream_inter_out,
        hls::stream<bool> &stream_end_out,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_desc_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    ap_uint<H_W_1> candidate;
    ap_uint<8> tableIndex, curQV;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    curQV = 0;

    bool last = stream_end_embed_in.read();
INTERSECT_COPYING_EMBEDDING_LOOP:
    while(!last){
        curEmb[curQV] = stream_embed_in.read();
        stream_embed_out.write(curEmb[curQV]);
        stream_end_embed_out.write(false);
        curQV++;
        last = stream_end_embed_in.read();
    }
    stream_end_embed_out.write(true);

    /* bypass data for extract phase */
    stream_min_desc_out.write(stream_min_desc_in.read());

    //std::cout << "Intersection: probing" << std::endl;
    last = stream_end_min_in.read();
INTERSECT_LOOP:
    while(!last){
        candidate = stream_min_in.read();
        //std::cout << "\t" << candidate << std::endl;
        bool inter = true;
INTERSECT_TBINDEXED_LOOP:
        for(int g = 0; 
                g < qVertices[curQV].numTablesIndexed && inter;
                g++)
        {
            ap_uint<(1UL << C_W)> offset;
            tableIndex = qVertices[curQV].tables_indexed[g];
            uint8_t ivPos = qVertices[curQV].vertex_indexing[g];

            stream_hash_in.write(curEmb[ivPos]);
            xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
            ap_uint<H_W_1> index1 = stream_hash_out.read().range(H_W_1 - 1, 0);
            ap_uint<H_W_2> index2 = candidate.range(H_W_2 - 1, 0);
            
            offset = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                    index1,
                    index2,
                    htb_buf,
                    hTables[tableIndex].start_offset);
            inter = offset.test((1UL << C_W) - 1);
        }

INTERSECT_TBINDEXING_LOOP:
        for(int g = 0; 
                g < qVertices[curQV].numTablesIndexing && inter; 
                g++)
        {
            ap_uint<(1UL << C_W)> start_off, end_off;
            tableIndex = qVertices[curQV].tables_indexing[g];
            
            if (candidate != 0){ 
                start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                        candidate - 1,
                        (1UL << H_W_2) - 1,
                        htb_buf,
                        hTables[tableIndex].start_offset);
                start_off = start_off.range((1UL << C_W) - 2, 0);
            } else {
                start_off = 0;
            }

            end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                    candidate,
                    (1UL << H_W_2) - 1,
                    htb_buf,
                    hTables[tableIndex].start_offset);
            end_off = end_off.range((1UL << C_W) - 2, 0);
            inter = (end_off > start_off);
        }

        if (inter){
            stream_inter_out.write(candidate);
            stream_end_out.write(false);
        }

        last = stream_end_min_in.read();
    }
    stream_end_out.write(true);
}

void mwj_extract_hashtovid(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<H_W_1>> &stream_inter_in,
        hls::stream<bool> &stream_end_in,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_set,
        AdjHT *hTables,
        ap_uint<512> *htb_buf,

        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_out,
        hls::stream<bool> &stream_set_end_out,
        hls::stream<bool> &stream_end_out)
{

    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    ap_uint<V_ID_W + 1 + 8> minData = stream_min_set.read();
    ap_uint<H_W_1> hashinter;
    ap_uint<V_ID_W> vertex;
    ap_uint<V_ID_W*2> edge;
    
    /* Bypassing solution for verify */
    bool last = stream_end_embed_in.read();
EXTRACT_COPYING_EMBEDDING_LOOP:
    while(!last){
        stream_embed_out.write(stream_embed_in.read());
        stream_end_embed_out.write(last);
        last = stream_end_embed_in.read();
    }
    stream_end_embed_out.write(last);

    /* filter to remove hash collisions when passing from H_W_1 to H_W_2 */
    ap_uint<16> filter[(1UL << (H_W_2 - 4))];
EXTRACT_RESET_FILTER_LOOP:
    for (int g = 0; g < (1UL << (H_W_2 - 4)); filter[g++] = 0);

    //std::cout << "\tHash to vid: " << std::endl;
    last = stream_end_in.read();
EXTRACT_HASHTOVID_LOOP:
    while(!last){
        stream_end_out.write(false);
        hashinter = stream_inter_in.read();
        ap_uint<(1UL << C_W)> start_off = 0;
        ap_uint<(1UL << C_W)> end_off = 0;
        ap_uint<H_W_1> index1, index1_f;
        ap_uint<H_W_2> index2, index2_f;

        if (minData.test(8)) {

            /* Retriving hash for vertex indexing the table */
            stream_hash_in.write(minData.range(V_ID_W + 8, 9));
            xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
            index1 = stream_hash_out.read().range(H_W_1 - 1, 0);
            index2 = hashinter.range(H_W_2 - 1, 0);
        } else {
            index1 = hashinter;
            index2 = 0;
        }

        if (index2 != 0){ 
            index1_f = index1;
            index2_f = index2 - 1;
        } else if (index2 == 0 && index1 != 0){
            index1_f = index1 - 1;
            index2_f = (1UL << H_W_2) - 1;
        }

        if (!(index2 == 0 && index1 == 0)){
            start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                    index1_f,
                    index2_f,
                    htb_buf,
                    hTables[minData.range(7,0)].start_offset);
            start_off = start_off.range((1UL << C_W) - 2, 0);
        }

        if (!minData.test(8)){
            index2 = (1UL << H_W_2) - 1;
        }

        end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                index1,
                index2,
                htb_buf,
                hTables[minData.range(7,0)].start_offset);
        end_off = end_off.range((1UL << C_W) - 2, 0);

        /* If hash bit already present do no stream anything */
        if (minData.test(8)){
            ap_uint<4> bitindex = index2.range(3, 0);
            ap_uint<H_W_2 - 4> arrayindex = index2.range(H_W_2 - 1, 4);
            ap_uint<16> filterword = filter[arrayindex];
            if (filterword[bitindex] == 1){
                end_off = start_off;
            } else {
                filterword[bitindex] = 1;
                filter[arrayindex] = filterword;
            }
        }

        //std::cout << "\t" << hashinter << ": " <<
        //    start_off << " to " << end_off << 
        //    "from tab " << minData.range(7,0) <<
        //    "which starts at " << hTables[minData.range(7,0)].start_offset <<
        //    std::endl;
EXTRACT_HASHTOVID_READ_LOOP:
        for (; start_off < end_off; start_off++){

            edge = read_table<32, 1, 0, E_W>(
                    start_off,
                    0,
                    htb_buf,
                    hTables[minData.range(7, 0)].start_edges);
            
            if (minData.test(8)){
                vertex = edge.range(V_ID_W - 1, 0);
            } else {
                vertex = edge.range(2*V_ID_W - 1, V_ID_W);
            }

            //std::cout << "\t\t" << vertex << std::endl;
            stream_inter_out.write(vertex);
            stream_set_end_out.write(false);            
        }

        stream_set_end_out.write(true);
        last = stream_end_in.read();
    }
    stream_end_out.write(true);
}

void mwj_extract_bagtoset(
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_set_end_in,
        hls::stream<bool> &stream_end_in,

        hls::stream<ap_uint<V_ID_W>> &stream_inter_out,
        hls::stream<bool> &stream_end_out)
{
   
    ap_uint<V_ID_W> set[MAX_CL];
    uint8_t counter = 0;
    bool last = stream_end_in.read();
EXTRACT_BAGTOSET_LOOP:
    while(!last){
        bool lastSet = stream_set_end_in.read();
        counter = 0;
EXTRACT_BAGTOSET_SET_LOOP:
        while(!lastSet){
            ap_uint<V_ID_W> vertex = stream_inter_in.read();
            uint8_t nSet = 0;
EXTRACT_BAGTOSET_SETCHECKER_LOOP:
            for(; nSet < counter; nSet++){
                if (vertex == set[nSet]){
                    break;
                }
            }
            if (nSet == counter){
                set[counter++] = vertex;
                stream_inter_out.write(vertex);
                stream_end_out.write(false);
            }
            lastSet = stream_set_end_in.read();    
        }
        last = stream_end_in.read();
    }
    stream_end_out.write(true);
}

void mwj_extract(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<H_W_1>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_in,
        AdjHT *hTables,
        ap_uint<512> *htb_buf,

        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<V_ID_W>> &stream_candidates_out,
        hls::stream<bool> &stream_end_candidates_out)
{

#pragma HLS DATAFLOW

    //std::cout << "Extract: " << std::endl;
    hls::stream<ap_uint<V_ID_W>, S_D> stream_htv_inter("Bags of vertices");
    hls::stream<bool, S_D> stream_htv_bag_end("Bag delimeter");
    hls::stream<bool, S_D> stream_htv_end("Bags delimeter");

    mwj_extract_hashtovid(
            stream_embed_in,
            stream_end_embed_in,
            stream_inter_in,
            stream_end_inter_in,
            stream_min_in,
            hTables,
            htb_buf,
            stream_embed_out,
            stream_end_embed_out,
            stream_htv_inter,
            stream_htv_bag_end,
            stream_htv_end);

    mwj_extract_bagtoset(
            stream_htv_inter,
            stream_htv_bag_end,
            stream_htv_end,
            stream_candidates_out,
            stream_end_candidates_out);
}

void mwj_verify_nothomomorphism(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<V_ID_W>> &stream_checked_out,
        hls::stream<bool> &stream_checked_end_out)
{
    ap_uint<8> curQV = 0;
    ap_uint<V_ID_W> curEmb[MAX_QV];

    bool last = stream_end_embed_in.read();
VERIFY_HOMOMO_COPYING_EMBEDDING_LOOP:
    while(!last){
        curEmb[curQV] = stream_embed_in.read();
        stream_embed_out.write(curEmb[curQV]);
        stream_end_embed_out.write(false);
        curQV++;
        last = stream_end_embed_in.read();
    }
    stream_end_embed_out.write(true);

    last = stream_end_inter_in.read();
VERIFY_CHECK_LOOP:
    while(!last){
        ap_uint<V_ID_W> vToVerify = stream_inter_in.read();
        bool homomorphism = false;
        for (int g = 0; g < curQV; g++){
            if (vToVerify == curEmb[g])
                homomorphism = true;
        }
        if (!homomorphism){
            stream_checked_out.write(vToVerify);
            stream_checked_end_out.write(false);
        }

        last = stream_end_inter_in.read();
    }
    stream_checked_end_out.write(true);
}

void mwj_verify_edge(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<512> *htb_buf,
        
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<V_ID_W>> &stream_checked_out,
        hls::stream<bool> &stream_checked_end_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    ap_uint<8> curQV = 0;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    //std::cout << "Verify: " << std::endl;

    bool last = stream_end_embed_in.read();
VERIFY_EDGE_COPYING_EMBEDDING_LOOP:
    while(!last){
        curEmb[curQV] = stream_embed_in.read();
        stream_embed_out.write(curEmb[curQV]);
        stream_end_embed_out.write(false);
        curQV++;
        last = stream_end_embed_in.read();
    }
    stream_end_embed_out.write(true);
    
    last = stream_end_inter_in.read();
VERIFY_CHECK_LOOP:
    while(!last){
        ap_uint<V_ID_W> vToVerify = stream_inter_in.read();
        bool checked = true;
        //hls::print("Ingresso: %d\n", (int)vToVerify);

VERIFY_CHECK_EDGE_LOOP:
        for(int g = 0; g < qVertices[curQV].numTablesIndexed && checked; g++){
            checked = false;
            ap_uint<H_W_1> index1, index1_f;
            ap_uint<H_W_2> index2, index2_f;
            uint8_t tableIndex = qVertices[curQV].tables_indexed[g];
            uint8_t ivPos = qVertices[curQV].vertex_indexing[g];
            stream_hash_in.write(curEmb[ivPos]);
            stream_hash_in.write(vToVerify);
            xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
            index1 = stream_hash_out.read().range(H_W_1 - 1, 0);
            xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
            index2 = stream_hash_out.read().range(H_W_2 - 1, 0);

            ap_uint<(1UL << C_W)> start_off = 0;
            ap_uint<(1UL << C_W)> end_off = 0;

            if (index2 != 0){
                index1_f = index1;
                index2_f = index2 - 1;
            } else if (index2 == 0 && index1 != 0){
                index1_f = index1 - 1;
                index2_f = (1UL << H_W_2)-1;
            }

            if (!(index2 == 0 && index1 == 0)){
                start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                        index1_f,
                        index2_f,
                        htb_buf,
                        hTables[tableIndex].start_offset);
                start_off = start_off.range((1UL << C_W) - 2, 0);
            }

            end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                    index1,
                    index2,
                    htb_buf,
                    hTables[tableIndex].start_offset);
            end_off = end_off.range((1UL << C_W) - 2, 0);

VERIFY_READ_MEMORY_LOOP:
            for (; start_off < end_off && checked == false; start_off++){
                ap_uint<2 * V_ID_W> edge = read_table<32, 1, 0, E_W>(
                        start_off,
                        0,
                        htb_buf,
                        hTables[tableIndex].start_edges);

                ap_uint<V_ID_W> vertexIndexed, vertexIndexing;
                vertexIndexed = edge.range(V_ID_W - 1, 0);
                vertexIndexing = edge.range(2*V_ID_W - 1, V_ID_W);
                if (vertexIndexing == curEmb[ivPos] &&
                        vertexIndexed == vToVerify){
                    checked = true;
                }
            }
        }

        if (checked){
            stream_checked_out.write(vToVerify);
            stream_checked_end_out.write(false);
        }

        last = stream_end_inter_in.read();
    }
    stream_checked_end_out.write(true);

}

void mwj_verify_add(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        ap_uint<8> nQueryVer,
        
        hls::stream<ap_uint<V_ID_W>> &stream_partial_out,
        hls::stream<bool> &stream_partial_end_out,
        hls::stream<ap_uint<V_ID_W>> &stream_final_out,
        hls::stream<bool> &stream_final_end_out)
{
    ap_uint<8> curQV = 0;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    //std::cout << "Verify: " << std::endl;

    bool last = stream_end_embed_in.read();
VERIFY_ADD_COPYING_EMBEDDING_LOOP:
    while(!last){
        curEmb[curQV] = stream_embed_in.read();
        curQV++;
        last = stream_end_embed_in.read();
    }

    last = stream_end_inter_in.read();
VERIFY_CHECK_LOOP:
    while(!last){
        ap_uint<V_ID_W> vToVerify = stream_inter_in.read();
        
        /* Write in the correct stream */
        if (curQV == nQueryVer){
            stream_final_end_out.write(false);
VERIFY_WRITE_FINAL_LOOP:
            for (int g = 0; g < curQV; g++){
                stream_final_out.write(curEmb[g]);
            }
            stream_final_out.write(vToVerify);
        } else {
VERIFY_WRITE_PARTIAL_LOOP:
            for (int g = 0; g < curQV; g++){
                stream_partial_out.write(curEmb[g]);
                stream_partial_end_out.write(false);
            }
            stream_partial_out.write(vToVerify);
            stream_partial_end_out.write(false);
            stream_partial_end_out.write(true);
            nPartialSol++;

#ifndef __SYNTHESIS__
            nStreamPar += (curQV + 1);
            if (nStreamPar > nStreamParMax){
                nStreamParMax = nStreamPar;
            }
#endif
        }

        last = stream_end_inter_in.read();
    }

    nPartialSol--;

#ifndef __SYNTHESIS__
    nStreamPar -= curQV;
#endif

    if (nPartialSol <= 0){
        stream_final_end_out.write(true);
    }
}

void mwj_verify(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<512> *htb_buf,
        ap_uint<8> nQueryVer,
        
        hls::stream<ap_uint<V_ID_W>> &stream_partial_out,
        hls::stream<bool> &stream_partial_end_out,
        hls::stream<ap_uint<V_ID_W>> &stream_final_out,
        hls::stream<bool> &stream_final_end_out)
{
#pragma HLS DATAFLOW

        hls::stream<ap_uint<V_ID_W>, S_D> stream_checked1;
        hls::stream<bool, S_D> stream_checked1_end;
        hls::stream<ap_uint<V_ID_W>, MAX_QV> stream_embed1;
        hls::stream<bool, MAX_QV> stream_embed1_end;
        
        hls::stream<ap_uint<V_ID_W>, S_D> stream_checked2;
        hls::stream<bool, S_D> stream_checked2_end;
        hls::stream<ap_uint<V_ID_W>, MAX_QV> stream_embed2;
        hls::stream<bool, MAX_QV> stream_embed2_end;

        mwj_verify_nothomomorphism(
                stream_embed_in,
                stream_end_embed_in,
                stream_inter_in,
                stream_end_inter_in,
                stream_embed1,
                stream_embed1_end,
                stream_checked1,
                stream_checked1_end);
        
        mwj_verify_edge(
                stream_embed1,
                stream_embed1_end,
                stream_checked1,
                stream_checked1_end,
                hTables,
                qVertices,
                htb_buf,
                stream_embed2,
                stream_embed2_end,
                stream_checked2,
                stream_checked2_end);

        mwj_verify_add(
                stream_embed2,
                stream_embed2_end,
                stream_checked2,
                stream_checked2_end,
                nQueryVer,
                stream_partial_out,
                stream_partial_end_out,
                stream_final_out,
                stream_final_end_out);
}

void mwj_proposeWrap(
        int nCycles,
        AdjHT *hTables,
        QueryVertex *qVertices,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        ap_uint<512> *htb_buf0,
        ap_uint<512> *htb_buf1,
        
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<H_W_1>> &stream_min_out,
        hls::stream<bool> &stream_end_min_out,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_set)
{
    for (int g=0; g < nCycles; g++){
        mwj_propose(
                hTables,
                qVertices,
                stream_embed_in,
                stream_end_embed_in,
                htb_buf0,
                htb_buf1,
                stream_embed_out,
                stream_end_embed_out,
                stream_min_out,
                stream_end_min_out,
                stream_min_set);
    }
}

void mwj_intersectWrap(
        int nCycles,
        AdjHT *hTables,
        QueryVertex *qVertices,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<H_W_1>> &stream_min_in,
        hls::stream<bool> &stream_end_min_in,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_set_in,
        ap_uint<512> *htb_buf,

        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<H_W_1>> &stream_inter_out,
        hls::stream<bool> &stream_end_out,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_set_out)
{
    for (int g = 0; g < nCycles; g++){
        mwj_intersect(
                hTables,
                qVertices,
                stream_embed_in,
                stream_end_embed_in,
                stream_min_in,
                stream_end_min_in,
                stream_min_set_in,
                htb_buf,
                stream_embed_out,
                stream_end_embed_out,
                stream_inter_out,
                stream_end_out,
                stream_min_set_out);
    }
}

void mwj_extractWrap(
        int nCycles,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<H_W_1>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_in,
        AdjHT *hTables,
        ap_uint<512> *htb_buf,

        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<V_ID_W>> &stream_candidates_out,
        hls::stream<bool> &stream_end_candidates_out)
{
    
    for (int g = 0; g < nCycles; g++){
        mwj_extract(
                stream_embed_in,
                stream_end_embed_in,
                stream_inter_in,
                stream_end_inter_in,
                stream_min_in,
                hTables,
                htb_buf,
                stream_embed_out,
                stream_end_embed_out,
                stream_candidates_out,
                stream_end_candidates_out);
    }
}

void mwj_verifyWrap(
        int nCycles,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<512> *htb_buf,
        ap_uint<8> nQueryVer,
        
        hls::stream<ap_uint<V_ID_W>> &stream_partial_out,
        hls::stream<bool> &stream_partial_end_out,
        hls::stream<ap_uint<V_ID_W>> &stream_final_out,
        hls::stream<bool> &stream_final_end_out)
{
    for (int g = 0; g < nCycles; g++){
        mwj_verify(
                stream_embed_in,
                stream_end_embed_in,
                stream_inter_in,
                stream_end_inter_in,
                hTables,
                qVertices,
                htb_buf,
                nQueryVer,
                stream_partial_out,
                stream_partial_end_out,
                stream_final_out,
                stream_final_end_out);
    }
}

void multiwayJoin(
        ap_uint<512> *htb_buf0,
        ap_uint<512> *htb_buf1,
        ap_uint<512> *htb_buf2,
        ap_uint<512> *htb_buf3,
        ap_uint<512> *htb_buf4,
        AdjHT *hTables0,
        AdjHT *hTables1,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        ap_uint<8> nQueryVer,
        int nPartEmbed,

        hls::stream<ap_uint<V_ID_W>> &stream_final_out,
        hls::stream<bool> &stream_final_out_end)
{
#pragma HLS STABLE variable=htb_buf0
#pragma HLS STABLE variable=htb_buf1
#pragma HLS STABLE variable=htb_buf2
#pragma HLS STABLE variable=htb_buf3
#pragma HLS STABLE variable=htb_buf4
#pragma HLS STABLE variable=hTables0
#pragma HLS STABLE variable=hTables1
#pragma HLS STABLE variable=qVertices0
#pragma HLS STABLE variable=qVertices1
#pragma HLS STABLE variable=nQueryVer
#pragma HLS DATAFLOW

    /* Propose data out*/
    hls::stream<ap_uint<V_ID_W>, MAX_QV> p_stream_embed
        ("Partial result propose");
    hls::stream<bool, MAX_QV> p_stream_embed_end("Partial result del. propose");
    
    hls::stream<ap_uint<H_W_1>, S_D> p_stream_min
        ("Min. set vertices propose");
    hls::stream<bool, S_D> p_stream_min_end("Min. set verices del. propose");
    
    hls::stream<ap_uint<V_ID_W + 1 + 8>, S_D> p_stream_min_desc
        ("Min. set description propose");

    /* Intersect data out */    
    hls::stream<ap_uint<H_W_1>, S_D> i_stream_hash_set("Set intersection hashes");
    hls::stream<bool, S_D> i_stream_hash_set_end("Set intersection hashes del.");
    hls::stream<ap_uint<V_ID_W>, MAX_QV> i_stream_embed
        ("Partial result intersect.");
    hls::stream<bool, MAX_QV> i_stream_embed_end
        ("Partial result del. intersect.");
    hls::stream<ap_uint<V_ID_W + 1 + 8>, S_D> i_stream_min_desc
        ("Min. set description propose");
    
    /* Extract data out */
    hls::stream<ap_uint<V_ID_W>, MAX_QV> e_stream_embed
        ("Partial result intersect.");
    hls::stream<bool, MAX_QV> e_stream_embed_end
        ("Partial result del. intersect.");
    hls::stream<ap_uint<V_ID_W>, S_D> e_stream_cand("Set candidates");
    hls::stream<bool, S_D> e_stream_cand_end("Set candidates delimeter");

    /* Verify data out */
	static hls::stream<ap_uint<V_ID_W>, S_D_RES> stream_embed
        ("Partial result stream");
	static hls::stream<bool, S_D_RES> stream_embed_end
        ("Partial result delimeter");
    
        
    mwj_proposeWrap(
            nPartEmbed,
            hTables0,
            qVertices0,
            stream_embed,
            stream_embed_end,
            htb_buf0,
            htb_buf1,
            p_stream_embed,
            p_stream_embed_end,
            p_stream_min,
            p_stream_min_end,
            p_stream_min_desc);

    mwj_intersectWrap(
            nPartEmbed,
            hTables1,
            qVertices1,
            p_stream_embed,
            p_stream_embed_end,
            p_stream_min,
            p_stream_min_end,
            p_stream_min_desc,
            htb_buf2,
            i_stream_embed,
            i_stream_embed_end,
            i_stream_hash_set,
            i_stream_hash_set_end,
            i_stream_min_desc);

    mwj_extractWrap(
            nPartEmbed,
            i_stream_embed,
            i_stream_embed_end,
            i_stream_hash_set,
            i_stream_hash_set_end,
            i_stream_min_desc,
            hTables0,
            htb_buf3,
            e_stream_embed,
            e_stream_embed_end,
            e_stream_cand,
            e_stream_cand_end);

    mwj_verifyWrap(
            nPartEmbed,
            e_stream_embed,
            e_stream_embed_end,
            e_stream_cand,
            e_stream_cand_end,
            hTables1,
            qVertices0,
            htb_buf4,
            nQueryVer,
            stream_embed,
            stream_embed_end,
            stream_final_out,
            stream_final_out_end);

}

void multiwayJoinWrap(
        ap_uint<512> *htb_buf0,
        ap_uint<512> *htb_buf1,
        ap_uint<512> *htb_buf2,
        ap_uint<512> *htb_buf3,
        ap_uint<512> *htb_buf4,
        AdjHT *hTables0,
        AdjHT *hTables1,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        ap_uint<8> nQueryVer,

        hls::stream<ap_uint<V_ID_W>> &stream_final_out,
        hls::stream<bool> &stream_final_out_end)
{

MULTIWAYJOIN_LOOP:
	do {
        multiwayJoin(
                htb_buf0,
                htb_buf1,
                htb_buf2,
                htb_buf3,
                htb_buf4,
                hTables0,
                hTables1,
                qVertices0,
                qVertices1,
                nQueryVer,
                nPartialSol,
                stream_final_out,
                stream_final_out_end);

    } while(nPartialSol);

}

void subgraphIsomorphism(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end,
        ap_uint<512> *htb_buf0,
        ap_uint<512> *htb_buf1,
        ap_uint<512> *htb_buf2,
        ap_uint<512> *htb_buf3,
        ap_uint<512> *htb_buf4,

        hls::stream<ap_uint<V_ID_W>> &stream_out,
        hls::stream<bool> &stream_end_out)
{

    QueryVertex qVertices0[MAX_QV], qVertices1[MAX_QV];
    TableDescriptor tDescriptors[MAX_TB];
    AdjHT hTables0[MAX_TB], hTables1[MAX_TB];
    ap_uint<8> numTables = 0;
    ap_uint<8> numQueryVert = 0;

    buildTableDescriptors(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end,
            qVertices0,
            qVertices1,
            tDescriptors,
            numTables,
            numQueryVert);

    fillTables(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end,
            hTables0,
            hTables1,
            htb_buf0,
            tDescriptors,
            numTables);

    multiwayJoinWrap(
            htb_buf0,
            htb_buf1,
            htb_buf2,
            htb_buf3,
            htb_buf4,
            hTables0,
            hTables1,
            qVertices0,
            qVertices1,
            numQueryVert,
            stream_out,
            stream_end_out);

#ifndef __SYNTHESIS__
    std::cout << "Max depth stream partial solutions: " << 
        nStreamParMax << std::endl;
#endif
}
