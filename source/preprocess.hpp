#define HLS_STREAM_THREAD_SAFE
#ifndef __SYNTHESIS__
#include <cassert>
#include <fstream>
#include <thread>
#include <unordered_map>
#endif

/* #include "hls_stream_sizeup.h" */
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
#include "hls_burst_maxi.h"
//#include "ap_utils.h"

#define HTB_SIZE    (1UL << (H_W_1 + H_W_2 - (9 - C_W)))
#define CNT_ROW     (1UL << (9 - C_W))
#define EDGE_ROW    (1UL << (9 - E_W))
#define NODE_ROW    (1UL << (9 - (E_W - 1)))
#define STOP_S      7

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
#define DDR_WORDS   RES_WIDTH
#define DDR_W       DDR_WORD

#ifndef __SYNTHESIS__
#define DEBUG_STATS
#include "debug.hpp"
#endif /*__SYNTHESIS__*/

typedef struct __attribute__((packed)) {
    ap_uint<V_ID_W> src, dst;
    ap_uint<V_L_W> lsrc, ldst;
} edge_t;

/* Builds the table descriptors based on the information
 * from the query graph. */
void buildTableDescriptors(
        edge_t htb_buf[DDR_WIDTH],
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        TableDescriptor *tDescriptors,
        unsigned short &numTables,
        unsigned short numQueryVert,
        unsigned short numQueryEdge)
{
    /* Translate from id of vertex to position in the order */
    ap_uint<8> fromNumToPos[MAX_QV];
    ap_uint<V_ID_W> *node_buf = (ap_uint<V_ID_W>*)htb_buf0;

    /* Filling information about query vertices and coping
     * the vertex order needed by multiway join */
FILL_ORDER_LOOP:
    for (int g = 0; g < numQueryVert; g++){
        ap_uint<V_ID_W> node = node_buf[g];
#ifndef __SYNTHESIS__
        assert(numQueryVert < (MAX_QV));
#endif
        fromNumToPos[node] = g;
    }
  
    ap_uint<V_L_W * 2> *lab_buf = (ap_uint<V_L_W * 2>*)htb_buf1;
    ap_uint<V_ID_W * 2> *edge_buf = (ap_uint<V_ID_W * 2>*)&htb_buf0[numQueryVert];

    /* Creating table descriptors */
CREATE_TABDESC_LOOP:
    for (int g = 0; g < numQueryEdge; g++) {
        bool dirEdge = false;
        ap_uint<V_ID_W * 2> edge = edge_buf[g];
        ap_uint<V_ID_W
        ap_uint<V_L_W> labeldst = labeldstif.data;
		ap_uint<V_L_W> labelsrc = labelsrcif.data;
		ap_uint<V_ID_W> nodedst = nodedstif.data;
		ap_uint<V_ID_W> nodesrc = nodesrcif.data;
        ap_uint<8> nodeSrcPos = fromNumToPos[nodesrc];
        ap_uint<8> nodeDstPos = fromNumToPos[nodedst];

        if (nodeSrcPos < nodeDstPos)
            dirEdge = true;

#ifndef __SYNTHESIS__
        std::cout << (unsigned int)nodesrc << "(" << (int)labelsrc << ")"
            << " -> " <<  (unsigned int)nodedst << "(" << (int)labeldst << ")" 
            << std::endl;
#endif

        /* Understanding if table already exists */
        uint8_t g = 0;
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

template <size_t WORD,
         size_t COUNTER>
ap_uint<(1UL << WORD)> utils_increasecounter(
        ap_uint<(1UL << WORD)> row,
        unsigned int pos)
{
#pragma HLS inline
    const unsigned int CR = (1UL << (WORD - COUNTER));
    ap_uint <(1UL << WORD)> mask = row >> (pos << COUNTER);
    mask += 1;
    mask <<= (pos << COUNTER);
    row <<= ((CR - pos) << COUNTER);
    row >>= ((CR - pos) << COUNTER);
    return row |= mask;
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
        ap_uint<DDR_W> *htb_buf)
{
    ap_uint<DDR_W> ram_row;
    ap_uint<H_W_1> index1;
    ap_uint<H_W_2> index2;
    ap_uint<(1UL << E_W)> edge;
    ap_uint<8> ntb;
    ap_uint<64> addr_outrow;
    unsigned int addr_inrow;
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
        htb_buf[addr_outrow] = utils_increasecounter<9, C_W>(htb_buf[addr_outrow], addr_inrow);
 
        hTables[ntb].n_edges++;     
        last = stream_end.read();
    }
}

/* Transform counters to offsets. */
void counterToOffset(
        ap_uint<8> numTables,
        AdjHT *hTables,
        ap_uint<DDR_W> *htb_buf)
{
    ap_uint<DDR_W> ram_row, ram_row_new;

COUNTER_TO_OFFSET_DDR_LOOP:
    for (unsigned int ntb = 0; ntb < numTables; ntb++){
        ap_uint<(1 << C_W)> base_addr = 0;
        ap_uint<(1UL << C_W)> prev_addr = 0;
        ap_uint<(1UL << C_W)> counter;
        ap_uint<32> hash_used = 0;

COUNTER_TO_OFFSET_TABLE_LOOP:
        for(unsigned int start = 0; start < HTB_SIZE; start++){
            ram_row = htb_buf[start + hTables[ntb].start_offset];
            for (int g = 0; g < CNT_ROW; g++){
                counter = ram_row.range((1UL << C_W) - 1, 0);
                ram_row_new >>= (1UL << C_W);
                ram_row >>= (1UL << C_W);

#ifndef __SYNTHESIS__
                assert(base_addr < (1UL << (1UL << (C_W))-1));
#endif /* __SYNTHESIS__ */
#ifdef DEBUG_STATS
                if (counter > debug::max_collisions) 
                    debug::max_collisions = counter;
                debug::avg_collisions += (float)counter / (1 << (H_W_1 + H_W_2));
#endif /* DEBUG_STATS */

                ram_row_new.range((CNT_ROW << C_W) - 1, (CNT_ROW - 1) << C_W) = base_addr;
                base_addr += counter;

                /* Store used bit */
                if (counter > 0){
                    ram_row_new.set((CNT_ROW << C_W) - 1);
                }

                /* Check if an hash is used */
                if (((start << (9 - C_W)) + g) % (1UL << H_W_2) == 0){
                    if(prev_addr < base_addr){
                        hash_used++;
                        prev_addr = base_addr;
                    }
                }
            }
            htb_buf[start + hTables[ntb].start_offset] = ram_row_new;
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
        ap_uint<DDR_W> *htb_buf)
{
    ap_uint<DDR_W> ram_row;
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
        /* ram_row = htb_buf[addr_outrow]; */
        /* offset = ram_row.range(((addr_inrow + 1) << C_W) - 2, addr_inrow << C_W); */
        /* ram_row.range(((addr_inrow + 1) << C_W) - 2, */
        /* addr_inrow << C_W) = offset + 1; */

        ram_row = htb_buf[addr_outrow];
        ap_uint <DDR_W> mask = ram_row >> (addr_inrow << C_W);
        offset = mask & ((1UL << (1UL << C_W)) - 1);
        mask += 1;
        mask <<= (addr_inrow << C_W);
        ram_row <<= ((CNT_ROW - addr_inrow) << C_W);
        ram_row >>= ((CNT_ROW - addr_inrow) << C_W);
        ram_row |= mask;
        htb_buf[addr_outrow] = ram_row;
        
        /* Compute address of row that will store the edge */
        addr_outrow = hTables[ntb].start_edges 
            + (offset >> (9 - E_W));

        /* Compute address of the edge inside the row */
        addr_inrow = offset.range((9 - E_W) - 1, 0);

        /* Read, modify and write the edge */
        ram_row = htb_buf[addr_outrow];
/* ram_row.range(((addr_inrow + 1) << E_W) - 1, */
/* addr_inrow << E_W) = edge; */
        mask = edge;
        mask <<= (addr_inrow << E_W);
        ram_row |= mask;
        htb_buf[addr_outrow] = ram_row;

        last = stream_end.read();
    }
}

void edgeToHash(
        hls::stream<T_NODE> &stream_src,
        hls::stream<T_NODE> &stream_dst,
        hls::stream<T_LABEL> &stream_src_l,
        hls::stream<T_LABEL> &stream_dst_l,
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
    
    bool last;
COUNT_OCCURENCIES_LOOP:
    do {
        T_LABEL labeldstif = stream_dst_l.read();
        T_LABEL labelsrcif = stream_src_l.read();
        T_NODE nodedstif = stream_dst.read();
        T_NODE nodesrcif = stream_src.read();
        ap_uint<V_L_W> labeldst = labeldstif.data;
		ap_uint<V_L_W> labelsrc = labelsrcif.data;
		ap_uint<V_ID_W> nodedst = nodedstif.data;
		ap_uint<V_ID_W> nodesrc = nodesrcif.data;

        /* Finding correct table */
        unsigned int g = 0;
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

#ifdef UNDIRECTED   
            if (tDescriptors[g].src_label == labeldst &&
                    tDescriptors[g].dst_label == labelsrc){

                ap_uint<V_ID_W> vertexIndexing, vertexIndexed; 
                if(!tDescriptors[g].dir){
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
#endif

        }
        last = nodesrcif.last;
    } while(!last);
    stream_end_out.write(true);
}

void countEdges(
        hls::stream<T_NODE> &stream_src,
        hls::stream<T_NODE> &stream_dst,
        hls::stream<T_LABEL> &stream_src_l,
        hls::stream<T_LABEL> &stream_dst_l,
        TableDescriptor *tDescriptors,
        AdjHT *hTables,
        ap_uint<DDR_W> *htb_buf,
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
        hls::stream<T_NODE> &stream_src,
        hls::stream<T_NODE> &stream_dst,
        hls::stream<T_LABEL> &stream_src_l,
        hls::stream<T_LABEL> &stream_dst_l,
        TableDescriptor *tDescriptors,
        AdjHT *hTables,
        ap_uint<DDR_W> *htb_buf,
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

#ifdef DEBUG_INTERFACE
        unsigned int &debif_endpreprocess,
#endif /* DEBUG_INTERFACE */
        
        hls::stream<T_NODE> &stream_src,
        hls::stream<T_NODE> &stream_dst,
        hls::stream<T_LABEL> &stream_src_l,
        hls::stream<T_LABEL> &stream_dst_l,
        AdjHT *hTables0,
        AdjHT *hTables1,
        ap_uint<DDR_W> *htb_buf,
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
    for (unsigned int ntb = 0; ntb < numTables; ntb++){
        hTables0[ntb].start_offset = start_addr;
        hTables0[ntb].n_edges = 0;
        start_addr += HTB_SIZE;
    }

    countEdges(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
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
#ifndef __SYNTHESIS__
        assert(start_addr < DDR_WIDTH);
#endif
        /* Fill the last row of tables's edges of 1s.
         * Because while reading the minimum set it is
         * read the entire line, and must be possible
         * to recognize which are real edges and which no */
        hTables1[ntb].start_offset = hTables0[ntb].start_offset;
        hTables1[ntb].start_edges = hTables0[ntb].start_edges;
        hTables1[ntb].n_edges = hTables0[ntb].n_edges;
        hTables1[ntb].hash_set = hTables0[ntb].hash_set;
        //htb_buf[start_addr - 1] = ~((ap_uint<DDR_W>)0);
    }

    writeEdges(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            tDescriptors,
            hTables0,
            htb_buf,
            numTables);

#ifdef DEBUG_INTERFACE
    debif_endpreprocess = hTables0[0].n_edges;
#endif /* DEBUG_INTERFACE */

#ifndef __SYNTHESIS__
    end_addr = start_addr;
    std::cout << "Occupied " << end_addr * 64 << " bytes, " << end_addr <<
        " words." << std::endl;
#endif

    /* std::ofstream f("CHECK.txt"); */
    /* for(ap_uint<8> tab = 0; tab < numTables; tab++){ */
    /* f << "Table " << tab << std::endl; */
    /* ap_uint<64> start = hTables0[tab].start_offset; */
    /* ap_uint<H_W_1> counter = 0; */
    /* for(ap_uint<64> addr = 0; */
    /* addr < HTB_SIZE; */
    /* addr++){ */
    /* ap_uint<DDR_W> row = htb_buf[start + addr]; */
    /* for(int g = 0; g < CNT_ROW; g++){ */
    /* if (((addr << (9 - C_W)) + g) % (1UL << H_W_2) == 0){ */
    /* f << counter << ": " << std::endl; */
    /* counter++; */
    /* } */

    /* ap_uint<(1UL << C_W)> edge = row.range(((g+1)<<C_W)-1, g<<C_W); */
    /* f << "\t" << edge.range((1UL << C_W)-2, 0) << " -> " */
    /* << edge.test((1UL << C_W)-1) << std::endl; */

    /* } */
    /* } */
    /* } */
    /* f << std::endl; */
    /* for(ap_uint<8> tab = 0; tab < numTables; tab++){ */
    /* f << "Table " << tab << std::endl; */
    /* ap_uint<64> start = hTables0[tab].start_edges; */
    /* for(ap_uint<64> addr = 0; */
    /* addr <= (hTables0[tab].n_edges >> (9-E_W)); */
    /* addr++){ */
    /* ap_uint<DDR_W> row = htb_buf[start + addr]; */
    /* for(int g = 0; g < EDGE_ROW; g++){ */
    /* ap_uint<(1UL << E_W)> edge = row.range(((g+1)<<E_W)-1, g<<E_W); */
    /* f << "\t" << edge.range(2*V_ID_W-1, V_ID_W) << " -> " */
    /* << edge.range(V_ID_W-1, 0) << std::endl; */
    /* } */
    /* } */
    /* } */
    /* f.close(); */
}
