#define HLS_STREAM_THREAD_SAFE
#include "cache.h"
#ifndef __SYNTHESIS__
#include <cassert>
#include <fstream>
#include <unordered_map>
#endif
#include <thread>

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
#include "preprocess.hpp"
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

#if VERIFY_CACHE
typedef cache< ap_uint<DDR_W>, true, false, 2,
        HASHTABLES_SPACE, 1, 1, 8, false, 1, 1,
        false, 7> cache_type;
#endif /* VERIFY_CACHE */

#ifndef __SYNTHESIS__
#define DEBUG_STATS
#include "debug.hpp"
#endif /*__SYNTHESIS__*/

template <size_t W_1,
        size_t W_2,
        size_t SHF,
        size_t T>
ap_uint<(1UL << T)> read_table(
        ap_uint<W_1> index1,
        ap_uint<W_2> index2,
        ap_uint<DDR_W> *htb_buf,
        ap_uint<32> start_addr)
{
    //static ap_uint<DDR_W> ram_row;
    ap_uint<DDR_W> ram_row;
    
    /* It is impossible that the first word read is the all 1s
     * since before edges are read offsets which are on 
     * the lower address space of the memory */
    //static ap_uint<32> prev_addr_row = ~((ap_uint<32>)0);
    unsigned long addr_row;
    unsigned long addr_inrow;
    ap_uint<64> addr_counter;

    addr_counter = index1;
    addr_counter <<= SHF;
    addr_counter += index2;

    /* Compute address of row storing the counter */
    addr_row = start_addr + (addr_counter >> (DDR_BIT - T));

    /* Compute address of data inside the row */
    addr_inrow = addr_counter.range((DDR_BIT - T) - 1, 0);

    /* Read the data */
    ram_row = htb_buf[addr_row];
    ram_row >>= (addr_inrow << T);

/* return ram_row.range(((addr_inrow + 1) << T) - 1, */
/* addr_inrow << T); */
    return ram_row;

}

#if VERIFY_CACHE
template <size_t W_1,
        size_t W_2,
        size_t SHF,
        size_t T>
ap_uint<(1UL << T)> read_table(
        unsigned int cache_port,
        ap_uint<W_1> index1,
        ap_uint<W_2> index2,
        cache_type &htb_buf,
        ap_uint<32> start_addr)
{
    ap_uint<DDR_W> ram_row;
    
    /* It is impossible that the first word read is the all 1s
     * since before edges are read offsets which are on 
     * the lower address space of the memory */
    //static ap_uint<32> prev_addr_row = ~((ap_uint<32>)0);
    unsigned long addr_row;
    unsigned long addr_inrow;
    ap_uint<64> addr_counter;

    addr_counter = index1;
    addr_counter <<= SHF;
    addr_counter += index2;

    /* Compute address of row storing the counter */
    addr_row = start_addr + (addr_counter >> (DDR_BIT - T));

    /* Compute address of data inside the row */
    addr_inrow = addr_counter.range((DDR_BIT - T) - 1, 0);

    /* Read the data */
    ram_row = htb_buf.get(addr_row, cache_port);
    ram_row >>= (addr_inrow << T);

    return ram_row;

}
#endif /* VERIFY_CACHE */
void mwj_propose(
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<DDR_W> *htb_buf,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool, 1> &stream_stop,

        hls::stream<ap_uint<V_ID_W>> &stream_min_out,
        hls::stream<bool> &stream_end_out,
        hls::stream<ap_uint<SET_INFO_WIDTH>> &stream_setinfo,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out)
{
    
    ap_uint<8> tableIndex;
    ap_uint<V_ID_W> curQV;
    ap_uint<V_ID_W> readv;
    ap_uint<32> minSize;
    ap_uint<32> minStart, minEnd, minOff;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<V_ID_W * 2> edge;
    ap_uint<V_ID_W> vertex, vertexCheck;
    ap_uint<SET_INFO_WIDTH> setinfo;
    bool stop;

    while (1) {
        if (stream_embed_in.read_nb(readv)){
            tableIndex = 0;
            minSize = (1UL << 32) - 1;

            if (readv.test(V_ID_W - 1)){
                // Delimiter, new solution
                curQV = readv.range(V_ID_W - 2, 0);
PROPOSE_COPYING_EMBEDDING_LOOP:
                for (int g = 0; g < curQV; g++){
#pragma HLS pipeline II=1
                    curEmb[g] = stream_embed_in.read();
                }
            } else {
                curEmb[curQV - 1] = readv;
            }

            for (int g = 0; g < curQV; g++) {
#pragma HLS pipeline II=1
                stream_embed_out.write(curEmb[g]);
                stream_end_embed_out.write(false);
            }
            stream_end_embed_out.write(true);

            /* Find sizes of sets in which the current query vertex
             * is indexed by an other query vertex */
PROPOSE_TBINDEXED_LOOP:
            for(int g = 0; g < qVertices[curQV].numTablesIndexed; g++){
#pragma HLS pipeline II=1
                ap_uint<(1UL << C_W)> start_off = 0;
                ap_uint<(1UL << C_W)> end_off;
                tableIndex = qVertices[curQV].tables_indexed[g];
                uint8_t ivPos = qVertices[curQV].vertex_indexing[g];

                ap_uint<64> hash_out;
                xf::database::details::hashlookup3_core<V_ID_W>(curEmb[ivPos], hash_out);
                ap_uint<H_W_1> index = hash_out.range(H_W_1 - 1, 0);

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

#ifdef DEBUG_STATS
                debug::findmin_reads += 2;
#endif
                if ((end_off - start_off) < minSize) {
                    minSize = end_off - start_off;
                    minOff = hTables[tableIndex].start_edges;
                    minStart = start_off;
                    minEnd = end_off;
                    setinfo.range(7, 0) = tableIndex;
                    setinfo.range(15, 8) = ivPos;
                }
            }
            
            stream_setinfo.write(setinfo);
            unsigned int rowstart = minOff + (minStart >> (DDR_BIT - E_W));
            unsigned int rowend = minOff + (minEnd >> (DDR_BIT - E_W));
            unsigned int window_left = minStart.range((DDR_BIT - E_W) - 1, 0);
            unsigned int window_right = minEnd.range((DDR_BIT - E_W) - 1, 0) + 
                (rowend - rowstart) * EDGE_ROW;
            unsigned int cnt = 0;
            ap_uint<V_ID_W> vIndexing = curEmb[setinfo.range(15, 8)];
            /* htb_buf.read_request(rowstart, rowend-rowstart+1); */

PROPOSE_READ_MIN_INDEXED_LOOP:
            for (int g = rowstart; g <= rowend; g++){
                row_t row = htb_buf[g];
                /* row_t row = htb_buf.read(); */
                for (int i = 0; i < EDGE_ROW; i++, cnt++){
                    if (cnt >= window_left && cnt < window_right){
                        edge = row.range((1UL << E_W) - 1, 0);
                        vertexCheck = edge.range(V_ID_W * 2 - 1, V_ID_W);
                        vertex = edge.range(V_ID_W - 1, 0);
                        if (vIndexing == vertexCheck){
                            stream_min_out.write(vertex);
                            stream_end_out.write(false);
                        }
                    }
                    row >>= (1UL << E_W);
                }
#ifdef DEBUG_STATS
                debug::readmin_reads++;
#endif
            }
            stream_end_out.write(true);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

void mwj_batchbuild(
        hls::stream<ap_uint<V_ID_W>> &stream_set_in,
        hls::stream<bool> &stream_set_end_in,
        hls::stream<ap_uint<SET_INFO_WIDTH>> &stream_setinfo_in,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<bool, 1> &stream_stop,
        
        hls::stream<ap_uint<V_ID_W>> &stream_batch_out,
        hls::stream<ap_uint<2>> &stream_batch_end_out,
        hls::stream<ap_uint<SET_INFO_WIDTH>> &stream_setinfo_out,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out)
{
    ap_uint<8> curQV;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<V_ID_W> setinfo;
    unsigned int batch_counter;
    bool stop, last;
    
    while(1){
        if (stream_end_embed_in.read_nb(last)){
            curQV = 0;
            batch_counter = 0;

BATCHBUILD_COPYING_EMBEDDING_LOOP:
            while(!last){
                curEmb[curQV] = stream_embed_in.read();
                stream_embed_out.write(curEmb[curQV]);
                stream_end_embed_out.write(false);
                curQV++;
                last = stream_end_embed_in.read();
            }
            stream_end_embed_out.write(true);

            setinfo = stream_setinfo_in.read();
            stream_setinfo_out.write(setinfo);
            last = stream_set_end_in.read();

BATCHBUILD_MAIN_LOOP:
            while(!last){

BATCHBUILD_MOVING_SET_LOOP:
                while(!last && batch_counter != (PROPOSE_BATCH_SIZE - 1)){
#pragma HLS pipeline II=1
                    ap_uint<V_ID_W> node = stream_set_in.read();
                    stream_batch_out.write(node);
                    stream_batch_end_out.write(0);
                    batch_counter++;
                    last = stream_set_end_in.read();
                }

                // Stream again partial solution if max batch size is reached
                if (!last && batch_counter == (PROPOSE_BATCH_SIZE - 1)){
                    stream_batch_end_out.write(1);
                    batch_counter = 0;

                    for (int g = 0; g < curQV; g++){
                        stream_embed_out.write(curEmb[g]);
                        stream_end_embed_out.write(false);
                    }
                    stream_end_embed_out.write(true);
                    stream_setinfo_out.write(setinfo);
                }
            }
            stream_batch_end_out.write(3);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

void mwj_intersect(
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<DDR_W> *htb_buf,
        hls::stream<ap_uint<V_ID_W>> &stream_min_in,
        hls::stream<ap_uint<2>> &stream_end_min_in,
        hls::stream<ap_uint<SET_INFO_WIDTH>> &stream_setinfo_in,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<bool, 1> &stream_stop,

        hls::stream<ap_uint<V_ID_W>> &stream_inter_out,
        hls::stream<ap_uint<2>> &stream_end_out,
        hls::stream<ap_uint<SET_INFO_WIDTH>> &stream_setinfo_out,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out)
{
    ap_uint<64> candidate_hash;
    ap_uint<V_ID_W> candidate_v;
    ap_uint<8> tableIndex, curQV;
    ap_uint<2> last_set;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<SET_INFO_WIDTH> setinfo;
    bool stop, last_sol;

    while (1) {
        if (stream_end_embed_in.read_nb(last_sol)){
            curQV = 0;

INTERSECT_COPYING_EMBEDDING_LOOP:
            while(!last_sol){
                curEmb[curQV] = stream_embed_in.read();
                stream_embed_out.write(curEmb[curQV]);
                stream_end_embed_out.write(false);
                curQV++;
                last_sol = stream_end_embed_in.read();
            }
            stream_end_embed_out.write(true);

            setinfo = stream_setinfo_in.read();
            stream_setinfo_out.write(setinfo);
            
            last_set = stream_end_min_in.read();
INTERSECT_LOOP:
            while(!last_set.test(0)){
                candidate_v = stream_min_in.read();
                xf::database::details::hashlookup3_core<V_ID_W>(candidate_v, candidate_hash);
                
                //std::cout << "\t" << candidate << std::endl;
                bool inter = true;
#ifdef DEBUG_STATS
                bool real_inter = true;
#endif
INTERSECT_TBINDEXED_LOOP:
                for(int g = 0; 
                        g < qVertices[curQV].numTablesIndexed;
                        /* g < qVertices[curQV].numTablesIndexed && inter; */
                        g++)
                {
                    ap_uint<(1UL << C_W)> offset;
                    tableIndex = qVertices[curQV].tables_indexed[g];
                    uint8_t ivPos = qVertices[curQV].vertex_indexing[g];
                    if (setinfo.range(7, 0) != tableIndex || setinfo.range(15, 8) != ivPos) {
                        ap_uint<64> hash_out;
                        xf::database::details::hashlookup3_core<V_ID_W>(curEmb[ivPos], hash_out);
                        ap_uint<H_W_1> index1 = hash_out.range(H_W_1 - 1, 0);
                        ap_uint<H_W_2> index2 = candidate_hash.range(H_W_2 - 1, 0);

                        offset = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                                index1,
                                index2,
                                htb_buf,
                                hTables[tableIndex].start_offset);
                        inter = inter && offset.test((1UL << C_W) - 1);
#ifdef DEBUG_STATS
                        debug::intersect_reads++;

                        /* Computing when there is an alias */
                        ap_uint<H_W_1> index1_f;
                        ap_uint<H_W_2> index2_f;
                        ap_uint<(1UL << C_W)> start_off;

                        bool checked = false;
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

                        for (; start_off < offset.range((1UL << C_W) - 2,0) && checked == false; start_off++){
                            ap_uint<2 * V_ID_W> edge = read_table<32, 1, 0, E_W>(
                                    start_off,
                                    0,
                                    htb_buf,
                                    hTables[tableIndex].start_edges);

                            ap_uint<V_ID_W> vertexIndexed, vertexIndexing;
                            vertexIndexed = edge.range(V_ID_W - 1, 0);
                            vertexIndexing = edge.range(2*V_ID_W - 1, V_ID_W);
                            if (vertexIndexing == curEmb[ivPos] &&
                                    vertexIndexed == candidate_v){
                                checked = true;
                            }
                        }

                        if (offset.test((1UL << C_W) - 1) && checked) {
                            /* True positive */
                            debug::intersect_bit_truepositive++;
                        } else if (offset.test((1UL << C_W) - 1) && !checked) {
                            /* False positive */
                            real_inter = false;
                            debug::intersect_bit_falsepositive++;
                        } else if (!offset.test((1UL << C_W) - 1) && !checked) {
                            /* True negative */
                            real_inter = false;
                            debug::intersect_bit_truenegative++;
                        } else if (!offset.test((1UL << C_W) - 1) && checked) {
                            /* False negative */
                            debug::intersect_bit_falsenegative++;
                        }
#endif
                    }
                }

#if INTERESCT_INDEXING_LOOP
INTERSECT_TBINDEXING_LOOP:
                for(int g = 0;
                        g < qVertices[curQV].numTablesIndexing;
/* g < qVertices[curQV].numTablesIndexing && inter; */
                        g++)
                {
                    ap_uint<(1UL << C_W)> start_off, end_off;
                    tableIndex = qVertices[curQV].tables_indexing[g];
                    ap_uint<H_W_1> candidate = candidate_hash.range(H_W_1 - 1, 0);
                    if (candidate != 0){
                        start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                                candidate - 1,
                                (1UL << H_W_2) - 1,
                                htb_buf,
                                hTables[tableIndex].start_offset);
                        start_off = start_off.range((1UL << C_W) - 2, 0);
#ifdef DEBUG_STATS
                        debug::intersect_reads++;
#endif
                    } else {
                        start_off = 0;
                    }

                    end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                            candidate,
                            (1UL << H_W_2) - 1,
                            htb_buf,
                            hTables[tableIndex].start_offset);
                    end_off = end_off.range((1UL << C_W) - 2, 0);
                    inter = inter && (end_off > start_off);

#ifdef DEBUG_STATS
                    bool temp = (end_off > start_off);
                    debug::intersect_reads++;

                    // Computing when there is an alias
                    bool checked = false;
                    for (; start_off < end_off && checked == false; start_off++){
                        ap_uint<2 * V_ID_W> edge = read_table<32, 1, 0, E_W>(
                                start_off,
                                0,
                                htb_buf,
                                hTables[tableIndex].start_edges);

                        ap_uint<V_ID_W> vertexIndexing;
                        vertexIndexing = edge.range(2*V_ID_W - 1, V_ID_W);
                        if (vertexIndexing == candidate_v){
                            checked = true;
                        }
                    }

                    if (temp && checked) {
                        /* //True positive */
                        debug::intersect_bit_truepositive++;
                    } else if (temp && !checked) {
                        /* //False positive */
                        real_inter = false;
                        debug::intersect_bit_falsepositive++;
                    } else if (!temp && !checked) {
                        /* //True negative */
                        real_inter = false;
                        debug::intersect_bit_truenegative++;
                    } else if (!temp && checked) {
                        /* //False negative */
                        debug::intersect_bit_falsenegative++;
                    }

#endif
                }
#endif /* INTERSECT_INDEXING_LOOP */
#ifdef DEBUG_STATS
                if (inter && !real_inter) {
                    /* //False positive */
                    debug::intersect_sol_falsepositive++;
                } else if (inter && real_inter) {
                    /* //True positive */
                    debug::intersect_sol_truepositive++;
                } else if (!inter && !real_inter) {
                    /* //True negative */
                    debug::intersect_sol_truenegative++;
                } else if (!inter && real_inter) {
                    /* //False negative */
                    debug::intersect_sol_falsenegative++;
                }
#endif

                if (inter){
                    stream_inter_out.write(candidate_v);
                    stream_end_out.write(last_set);
                }

                last_set = stream_end_min_in.read();
            }
            stream_end_out.write(last_set);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

void mwj_homomorphism(
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        hls::stream<ap_uint<SET_INFO_WIDTH>> &stream_setinfo_in,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<bool, 1> &stream_stop,
        
        hls::stream<ap_uint<V_ID_W>> &stream_batch_out,
        hls::stream<bool> &stream_batch_end_out,
        hls::stream<ap_uint<SET_INFO_WIDTH>> &stream_setinfo_out,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out)
{
    ap_uint<8> curQV;
    ap_uint<2> last_set;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool stop, last_sol;

    while(1){
        if (stream_end_embed_in.read_nb(last_sol)){
            curQV = 0;

VERIFY_HOMOMO_COPYING_EMBEDDING_LOOP:
            while(!last_sol){
                curEmb[curQV] = stream_embed_in.read();
                stream_embed_out.write(curEmb[curQV]);
                stream_end_embed_out.write(false);
                curQV++;
                last_sol = stream_end_embed_in.read();
            }
            stream_end_embed_out.write(true);

            stream_setinfo_out.write(stream_setinfo_in.read());
            last_set = stream_end_inter_in.read();

VERIFY_CHECK_LOOP:
            while(!last_set.test(0)){
                ap_uint<V_ID_W> vToVerify = stream_inter_in.read();
                bool homomorphism = false;
         
                for (int g = 0; g < curQV; g++){
                    if (vToVerify == curEmb[g])
                        homomorphism = true;
                }

                if (!homomorphism){
                    stream_batch_out.write(vToVerify);
                    stream_batch_end_out.write(last_set);
                }

                last_set = stream_end_inter_in.read();
            }
            stream_batch_end_out.write(last_set);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

template <typename T>
void mwj_verify(
        AdjHT *hTables,
        QueryVertex *qVertices,
        T htb_buf,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<ap_uint<2>> &stream_end_inter_in,
        hls::stream<ap_uint<SET_INFO_WIDTH>> &stream_setinfo_in,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<bool, 1> &stream_stop,

        hls::stream<ap_uint<V_ID_W>> &stream_checked_out,
        hls::stream<ap_uint<2>> &stream_checked_end_out,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out)
{
    ap_uint<8> curQV;
    ap_uint<2> last_set;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool stop, last_sol;
    ap_uint<SET_INFO_WIDTH> setinfo;
    ap_uint<V_ID_W> vToVerify;
    unsigned short nEdgesToVerify;
    unsigned short buffer_pr, buffer_pw;
    unsigned short buffer_size;
    bool last_batch;
    ap_uint<V_ID_W> buffer[PROPOSE_BATCH_SIZE];
#pragma HLS bind_storage variable=buffer type=ram_1p impl=bram

    while(1){
        if (stream_end_embed_in.read_nb(last_sol)){
            curQV = 0;
            buffer_size = 0;

VERIFY_EDGE_COPYING_EMBEDDING_LOOP:
            while(!last_sol){
                curEmb[curQV] = stream_embed_in.read();
                stream_embed_out.write(curEmb[curQV]);
                stream_end_embed_out.write(false);
                curQV++;
                last_sol = stream_end_embed_in.read();
            }
            stream_end_embed_out.write(true);

            setinfo = stream_setinfo_in.read();
            nEdgesToVerify = qVertices[curQV].numTablesIndexed;

            for (int g = 0; g < nEdgesToVerify; g++){
                buffer_pw = buffer_pr = 0; 
                uint8_t tableIndex = qVertices[curQV].tables_indexed[g];
                uint8_t ivPos = qVertices[curQV].vertex_indexing[g];
                if (g == 0) {
                    last_set = stream_end_inter_in.read();
                    last_batch = last_set[1];
                } else {
                        last_set[1] = last_batch;
                        last_set[0] = (buffer_pr == buffer_size);
                }

VERIFY_CHECK_LOOP:
                while(!last_set.test(0)){
                     
                    if (g == 0){
                        vToVerify = stream_inter_in.read();
                    } else {
                        vToVerify = buffer[buffer_pr++];
                    }

                    bool checked = true;
                    ap_uint<H_W_1> index1, index1_f;
                    ap_uint<H_W_2> index2, index2_f;

                    if (setinfo.range(7, 0) != tableIndex || 
                            setinfo.range(15, 8) != ivPos) {
                        checked = false;
                        ap_uint<64> hash_out;
                        xf::database::details::hashlookup3_core<V_ID_W>(curEmb[ivPos], hash_out);
                        index1 = hash_out.range(H_W_1 - 1, 0);
                        xf::database::details::hashlookup3_core<V_ID_W>(vToVerify, hash_out);
                        index2 = hash_out.range(H_W_2 - 1, 0);

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
                                    0,
                                    index1_f,
                                    index2_f,
                                    htb_buf,
                                    hTables[tableIndex].start_offset);
                            start_off = start_off.range((1UL << C_W) - 2, 0);
#ifdef DEBUG_STATS
                            debug::verify_reads++;
#endif
                        }

                        end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                                0,
                                index1,
                                index2,
                                htb_buf,
                                hTables[tableIndex].start_offset);
                        end_off = end_off.range((1UL << C_W) - 2, 0);
#ifdef DEBUG_STATS
                        debug::verify_reads++;
#endif

VERIFY_READ_MEMORY_LOOP:
                        for (; start_off < end_off; start_off++){
                            ap_uint<2 * V_ID_W> edge = read_table<32, 1, 0, E_W>(
                                    1,
                                    start_off,
                                    0,
                                    htb_buf,
                                    hTables[tableIndex].start_edges);

                            ap_uint<V_ID_W> vertexIndexed, vertexIndexing;
                            vertexIndexed = edge.range(V_ID_W - 1, 0);
                            vertexIndexing = edge.range(2 * V_ID_W - 1, V_ID_W);
                            if (vertexIndexing == curEmb[ivPos] &&
                                    vertexIndexed == vToVerify){
                                checked = true;
                            }
#ifdef DEBUG_STATS
                            debug::verify_reads++;
#endif
                        }
                    }

                    if (checked){
                        if (g == (nEdgesToVerify - 1)){
/* std::cout << g << " Out " << vToVerify << std::endl; */
                            stream_checked_out.write(vToVerify);
                            stream_checked_end_out.write(last_set);
/* std::cout << last_set << std::endl; */
#ifdef DEBUG_STATS
                            debug::solution_correct++;
#endif
                        } else {
/* std::cout << g << " Writing " << vToVerify << " in " << buffer_pw << std::endl; */
                            buffer[buffer_pw++] = vToVerify;
                        } 
                    }

#ifdef DEBUG_STATS
                    else {
                        debug::solution_wrong++;
                    }
#endif
                    if (g == 0) {
                        last_set = stream_end_inter_in.read();
                        last_batch = last_set[1];
                    } else {
                        last_set[1] = last_batch;
                        last_set[0] = (buffer_pr == buffer_size);
                    }

                }

                if (g == (nEdgesToVerify - 1)){
/* std::cout << last_set << std::endl; */
                    stream_checked_end_out.write(last_set);
                } else {
/* std::cout << g << " Set size " << buffer_pw << std::endl; */
                    buffer_size = buffer_pw;
                } 
            }
        }

        if (stream_stop.read_nb(stop)){
            break;
        }
    }
}

void mwj_assembly(
        unsigned short nQueryVer,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<ap_uint<2>> &stream_end_inter_in,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<V_ID_W>> &stream_batch,
        hls::stream<bool> &stream_batch_end,
        hls::stream<bool, 1> streams_stop[STOP_S],
        
        hls::stream<ap_uint<V_ID_W>> &stream_partial_out,
#ifdef COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif
        )
{
    ap_uint<V_ID_W> curQV;
    ap_uint<2> last_set;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<SET_INFO_WIDTH> setinfo;
    unsigned long int nPartSol {0};
    bool last_sol;
    T_NODE node;

#ifdef COUNT_ONLY
    unsigned long int counter {0};
#endif

    last_sol = stream_batch_end.read();
    curQV = (1UL << (V_ID_W - 1)) + 1;
    stream_partial_out.write(curQV);
    while(!last_sol){
        stream_partial_out.write(stream_batch.read());
        nPartSol++;
        last_sol = stream_batch_end.read();
    }
    
    while(1) {
        if (stream_end_embed_in.read_nb(last_sol)){
            
            curQV = 0;
VERIFY_ADD_COPYING_EMBEDDING_LOOP:
            while(!last_sol){
                curEmb[curQV] = stream_embed_in.read();
                curQV++;
                last_sol = stream_end_embed_in.read();
            }

            last_set = stream_end_inter_in.read();
            
            // Write delimiter for new solutions
            if (!last_set.test(0) && (curQV != nQueryVer - 1)){
                stream_partial_out.write((curQV + 1) | (1UL << (V_ID_W - 1)));
                for (int g = 0; g < curQV; g++){
                    stream_partial_out.write(curEmb[g]);
                }
            }

VERIFY_CHECK_LOOP:
            while(!last_set.test(0)){
                ap_uint<V_ID_W> vToVerify = stream_inter_in.read();
                
                /* Write in the correct stream */
                if (curQV == nQueryVer - 1){
#ifdef COUNT_ONLY
                    counter++;
#else
VERIFY_WRITE_FINAL_LOOP:
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
VERIFY_WRITE_PARTIAL_LOOP:
                    stream_partial_out.write(vToVerify);
                    nPartSol++;
                }
#ifdef DEBUG_STATS
                debug::embeddings++;
#endif
                last_set = stream_end_inter_in.read();
            }

            // Last batch of a set 
            if (last_set.test(1))
                nPartSol--;
        }

        if (nPartSol == 0){
            for (int g = 0; g < STOP_S; g++){
#pragma HLS unroll
                streams_stop[g].write(true);
            }

#ifdef COUNT_ONLY
            /* Write in output number of results */
            result += counter;
#else
            /* Write last node */
            node.data = 0;
            node.last = true;
            node.keep = ~0;
            result.write(node);
#endif
            break;
        }
    }
}

#if VERIFY_CACHE
void mwj_verifyWrapper(
        AdjHT *hTables,
        QueryVertex *qVertices,
        cache_type &htb_buf,
        hls::stream<ap_uint<V_ID_W>> &stream_set_in,
        hls::stream<ap_uint<2>> &stream_set_end_in,
        hls::stream<ap_uint<SET_INFO_WIDTH>> &stream_set_info_in,
        hls::stream<ap_uint<V_ID_W>> &stream_sol_in,
        hls::stream<bool> &stream_sol_end_in,
        hls::stream<bool, 1> &stream_stop,
        
        hls::stream<ap_uint<V_ID_W>> &stream_set_out,
        hls::stream<ap_uint<2>> &stream_set_end_out,
        hls::stream<ap_uint<V_ID_W>> &stream_sol_out,
        hls::stream<bool> &stream_sol_end_out)
{
    htb_buf.init();
    std::thread mwj_verify_t(
            mwj_verify<cache_type&>,
            hTables,
            qVertices,
            std::ref(htb_buf),
            std::ref(stream_set_in),
            std::ref(stream_set_end_in),
            std::ref(stream_set_info_in),
            std::ref(stream_sol_in),
            std::ref(stream_sol_end_in),
            std::ref(stream_stop),
            std::ref(stream_set_out),
            std::ref(stream_set_end_out),
            std::ref(stream_sol_out),
            std::ref(stream_sol_end_out));
    mwj_verify_t.join();
    htb_buf.stop();
}
#endif /* VERIFY_CACHE */

template <size_t BATCH_SIZE>
void mwj_batch(
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<DDR_W> *htb_buf,
        
        hls::stream<bool> &stream_batch_end,
        hls::stream<bool> &stream_batches_end,
        hls::stream<ap_uint<V_ID_W>> &stream_batch)
{
    ap_uint<8> tableIndex {0};
    ap_uint<32> minSize = (1UL << 32) - 1;
    ap_uint<32> minStart, minEnd, minOff;
    ap_uint<8 + 1 + V_ID_W> minData;
    ap_uint<V_ID_W*2> edge;
    ap_uint<V_ID_W> vertex, vertexCheck;
    ap_uint<V_ID_W> set[MAX_CL];
    ap_uint<H_W_1> hash_buff = 0;
    ap_uint<5> set_counter = 0;
    bool flag_buff = false;
    bool flag_new = true;
    unsigned int batch_counter = BATCH_SIZE + 1;

PROPOSE_TBINDEXING_LOOP:
    for(int g = 0; g < qVertices[0].numTablesIndexing; g++){
        tableIndex = qVertices[0].tables_indexing[g];

        if (hTables[tableIndex].hash_set < minSize){
            minSize = hTables[tableIndex].hash_set;
            minOff = hTables[tableIndex].start_edges;
            minStart = 0;
            minEnd = hTables[tableIndex].n_edges;
            minData.range(7, 0) = tableIndex;
            minData.clear(8);
            minData.range(V_ID_W + 8, 9) = 0;
        }
    }

    unsigned int rowstart = minOff + (minStart >> (DDR_BIT - E_W));
    unsigned int rowend = minOff + (minEnd >> (DDR_BIT - E_W));
    unsigned int window_left = minStart.range((DDR_BIT - E_W) - 1, 0);
    unsigned int window_right = minEnd.range((DDR_BIT - E_W) - 1, 0) + 
        (rowend - rowstart) * EDGE_ROW;
    unsigned int cnt = 0;

PROPOSE_READ_MIN_INDEXING_LOOP:
    for (int g = rowstart; g <= rowend; g++){
        row_t row = htb_buf[g];
        /* row_t row = htb_buf.read(); */
        for (int i = 0; i < EDGE_ROW; i++, cnt++){
            if (cnt < window_right){
                edge = row.range((1UL << E_W) - 1, 0);
                vertex = edge.range(V_ID_W * 2 - 1, V_ID_W);
                ap_uint<64> hash_out;
                xf::database::details::hashlookup3_core<V_ID_W>(
                        vertex,
                        hash_out);
                ap_uint<H_W_1> vertexHash = hash_out.range(H_W_1 - 1, 0);

                if (flag_buff && hash_buff == vertexHash){
                    int nSet = 0;
                    flag_new = true;
EXTRACT_BAGTOSET_SETCHECKER_LOOP:
                    for(; nSet < set_counter; nSet++){
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
                    if(batch_counter == BATCH_SIZE){
                        stream_batches_end.write(false);
                        stream_batch_end.write(true);
                    }

                    set[set_counter++] = vertex;
                    batch_counter = (batch_counter >= BATCH_SIZE)? 1 : batch_counter + 1;
                    stream_batch_end.write(false);
                    stream_batch.write(vertex);
                }

                hash_buff = vertexHash;
                flag_buff = true;
            }
            row >>= (1UL << E_W);
        }
    }
    stream_batch_end.write(true);
    stream_batches_end.write(true);
}

void multiwayJoin(
        ap_uint<DDR_W> *htb_buf0,
        ap_uint<DDR_W> *htb_buf1,
/* hls::burst_maxi<row_t> htb_buf1, */
        ap_uint<DDR_W> *htb_buf2,
        row_t *res_buf,
        AdjHT *hTables0,
        AdjHT *hTables1,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        unsigned short nQueryVer,
        hls::stream<ap_uint<V_ID_W>> &stream_batch,
        hls::stream<bool> &stream_batch_end,
        
#ifdef COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif
        )
{
#pragma HLS STABLE variable=htb_buf0
#pragma HLS STABLE variable=htb_buf1
#pragma HLS STABLE variable=htb_buf2
#pragma HLS STABLE variable=hTables0
#pragma HLS STABLE variable=hTables1
#pragma HLS STABLE variable=qVertices0
#pragma HLS STABLE variable=qVertices1
#pragma HLS STABLE variable=nQueryVer
#pragma HLS DATAFLOW

    /* Propose data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> p_stream_sol
        ("Propose - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> p_stream_sol_end
        ("Propose - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> p_stream_set
        ("Propose - set nodes");
    hls_thread_local hls::stream<bool, S_D> p_stream_set_end
        ("Propose - set nodes end flag");
    hls_thread_local hls::stream<ap_uint<SET_INFO_WIDTH>, S_D> p_stream_set_info
        ("Propose - set info");

    /* Batchbuilder data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> b_stream_sol
        ("Batchbuild - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> b_stream_sol_end
        ("Batchbuild - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> b_stream_set
        ("Batchbuild - set nodes");
    hls_thread_local hls::stream<ap_uint<2>, S_D> b_stream_set_end
        ("Batchbuild - set nodes end flag");
    hls_thread_local hls::stream<ap_uint<SET_INFO_WIDTH>, S_D> b_stream_set_info
        ("Batchbuild - set info");
    
    /* Intersect data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> i_stream_sol
        ("Intersect - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> i_stream_sol_end
        ("Intersect - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> i_stream_set
        ("Intersect - set nodes");
    hls_thread_local hls::stream<ap_uint<2>, S_D> i_stream_set_end
        ("Intersect - set nodes end flag");
    hls_thread_local hls::stream<ap_uint<SET_INFO_WIDTH>, S_D> i_stream_set_info
        ("Intersect - set info");
    
    /* Homomorphism data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> h_stream_sol
        ("Homomorphism - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> h_stream_sol_end
        ("Homomorphism - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> h_stream_set
        ("Homomorphism - set nodes");
    hls_thread_local hls::stream<bool, S_D> h_stream_set_end
        ("Homomorphism - set nodes end flag");
    hls_thread_local hls::stream<ap_uint<SET_INFO_WIDTH>, S_D> h_stream_set_info
        ("Homomorphism - set info");
    
    /* Verify data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> v_stream_sol
        ("Verify - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> v_stream_sol_end
        ("Verify - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> v_stream_set
        ("Verify - set nodes");
    hls_thread_local hls::stream<ap_uint<2>, S_D> v_stream_set_end
        ("Verify - set nodes end flag");

    /* Assembly data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, 32> a_stream_sol
        ("Assembly - partial solution");

    /* Dynamic fifo data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, 32> dyn_stream_sol
        ("Dynamic fifo - partial solution");
 
    /* Stop signals */
    hls_thread_local hls::stream<bool, 1> streams_stop[STOP_S];

#ifndef __SYNTHESIS__
    for (int g = 0; g < STOP_S; g++) {
        char stream_name[10];
        sprintf(stream_name, "stop_%d", g);
        streams_stop[g].set_name(stream_name);
    } 
#endif /* __SYNTHESIS__ */

#if VERIFY_CACHE
    cache_type a_cache(htb_buf2);
#endif /* VERIFY_CACHE */

    dynfifo_init<
        ap_uint<V_ID_W>,        /* fifo data type */
        row_t,                  /* fifo data type */
        S_DEPTH,                /* in/out stream size */
        BURST_S*2,              /* load/store stream size */
        DDR_WORD,               /* bitwidth ddr word */
        BURST_S,                /* burst transaction size */
        RESULTS_SPACE>          /* memory words available */
            (a_stream_sol,
             dyn_stream_sol,
             streams_stop[STOP_S - 2],
             streams_stop[STOP_S - 1],
             res_buf);

#ifdef __SYNTHESIS__ 

    mwj_propose(
            hTables0,
            qVertices0,
            htb_buf0,
            dyn_stream_sol,
            streams_stop[0],
            p_stream_set,
            p_stream_set_end,
            p_stream_set_info,
            p_stream_sol,
            p_stream_sol_end);

    mwj_homomorphism(
            p_stream_set,
            p_stream_set_end,
            p_stream_set_info,
            p_stream_sol,
            p_stream_sol_end,
            streams_stop[3],
            h_stream_set,
            h_stream_set_end,
            h_stream_set_info,
            h_stream_sol,
            h_stream_sol_end);
    
    mwj_batchbuild(
            h_stream_set,
            h_stream_set_end,
            h_stream_set_info,
            h_stream_sol,
            h_stream_sol_end,
            streams_stop[1],
            b_stream_set,
            b_stream_set_end,
            b_stream_set_info,
            b_stream_sol,
            b_stream_sol_end);
    
    mwj_intersect(
            hTables1,
            qVertices1,
            htb_buf1,
            b_stream_set,
            b_stream_set_end,
            b_stream_set_info,
            b_stream_sol,
            b_stream_sol_end,
            streams_stop[2],
            i_stream_set,
            i_stream_set_end,
            i_stream_set_info,
            i_stream_sol,
            i_stream_sol_end);


#if VERIFY_CACHE

    cache_wrapper (
            mwj_verify<cache_type&>,
            hTables1,
            qVertices0,
            a_cache,
            i_stream_set,
            i_stream_set_end,
            i_stream_set_info,
            i_stream_sol,
            i_stream_sol_end,
            streams_stop[4],
            v_stream_set,
            v_stream_set_end,
            v_stream_sol,
            v_stream_sol_end);
#else

    mwj_verify<ap_uint<DDR_W>*>(
            hTables1,
            qVertices0,
            htb_buf2,
            i_stream_set,
            i_stream_set_end,
            i_stream_set_info,
            i_stream_sol,
            i_stream_sol_end,
            streams_stop[4],
            v_stream_set,
            v_stream_set_end,
            v_stream_sol,
            v_stream_sol_end);

#endif /* VERIFY_CACHE */

    mwj_assembly(
            nQueryVer,
            v_stream_set,
            v_stream_set_end,
            v_stream_sol,
            v_stream_sol_end,
            stream_batch,
            stream_batch_end,
            streams_stop,
            a_stream_sol,
            result);
    
#else
  
    for (int g = 0; g < STOP_S - 2; g++) 
        hls::stream_globals::incr_task_counter();
    
    std::thread mwj_propose_t(
            mwj_propose,
            hTables0,
            qVertices0,
            htb_buf0,
            std::ref(dyn_stream_sol),
            std::ref(streams_stop[0]),
            std::ref(p_stream_set),
            std::ref(p_stream_set_end),
            std::ref(p_stream_set_info),
            std::ref(p_stream_sol),
            std::ref(p_stream_sol_end));

    std::thread mwj_homomorphism_t(
            mwj_homomorphism,
            std::ref(p_stream_set),
            std::ref(p_stream_set_end),
            std::ref(p_stream_set_info),
            std::ref(p_stream_sol),
            std::ref(p_stream_sol_end),
            std::ref(streams_stop[3]),
            std::ref(h_stream_set),
            std::ref(h_stream_set_end),
            std::ref(h_stream_set_info),
            std::ref(h_stream_sol),
            std::ref(h_stream_sol_end));
    
    std::thread mwj_batchbuild_t(
            mwj_batchbuild,
            std::ref(h_stream_set),
            std::ref(h_stream_set_end),
            std::ref(h_stream_set_info),
            std::ref(h_stream_sol),
            std::ref(h_stream_sol_end),
            std::ref(streams_stop[1]),
            std::ref(b_stream_set),
            std::ref(b_stream_set_end),
            std::ref(b_stream_set_info),
            std::ref(b_stream_sol),
            std::ref(b_stream_sol_end));
    
    std::thread mwj_intersect_t(
            mwj_intersect,
            hTables1,
            qVertices1,
            htb_buf1,
            std::ref(b_stream_set),
            std::ref(b_stream_set_end),
            std::ref(b_stream_set_info),
            std::ref(b_stream_sol),
            std::ref(b_stream_sol_end),
            std::ref(streams_stop[2]),
            std::ref(i_stream_set),
            std::ref(i_stream_set_end),
            std::ref(i_stream_set_info),
            std::ref(i_stream_sol),
            std::ref(i_stream_sol_end));

    
    std::thread mwj_assembly_t(
            mwj_assembly,
            nQueryVer,
            std::ref(v_stream_set),
            std::ref(v_stream_set_end),
            std::ref(v_stream_sol),
            std::ref(v_stream_sol_end),
            std::ref(stream_batch),
            std::ref(stream_batch_end),
            std::ref(streams_stop),
            std::ref(a_stream_sol),
            std::ref(result));

#if VERIFY_CACHE

    mwj_verifyWrapper(
            hTables1,
            qVertices0,
            a_cache,
            i_stream_set,
            i_stream_set_end,
            i_stream_set_info,
            i_stream_sol,
            i_stream_sol_end,
            streams_stop[4],
            v_stream_set,
            v_stream_set_end,
            v_stream_sol,
            v_stream_sol_end);

#else

    std::thread mwj_verify_t(
            mwj_verify<ap_uint<DDR_W>>,
            hTables1,
            qVertices0,
            htb_buf2,
            std::ref(i_stream_set),
            std::ref(i_stream_set_end),
            std::ref(i_stream_set_info),
            std::ref(i_stream_sol),
            std::ref(i_stream_sol_end),
            std::ref(streams_stop[4]),
            std::ref(v_stream_set),
            std::ref(v_stream_set_end),
            std::ref(v_stream_sol),
            std::ref(v_stream_sol_end));

#endif /* VERIFY_CACHE */

    mwj_assembly_t.join();
    mwj_batchbuild_t.join();
    mwj_propose_t.join();
    mwj_intersect_t.join();
    mwj_homomorphism_t.join(); 

#if !VERIFY_CACHE
    mwj_verify_t.join();
#else
#ifdef DEBUG_STATS
    debug::cache_hit_0 += a_cache.get_hit_ratio(0);
    debug::cache_hit_1 += a_cache.get_hit_ratio(1);
    debug::cache_req_0 += a_cache.get_n_l1_reqs(0);
    debug::cache_req_1 += a_cache.get_n_l1_reqs(1);
#endif /* DEBUG_STATS */
#endif /* VERIFY_CACHE */
    

#endif /* __SYNTHESIS__ */

}

void multiwayJoinWrap(
        ap_uint<DDR_W> *htb_buf0,
        ap_uint<DDR_W> *htb_buf1,
/* hls::burst_maxi<row_t> htb_buf1, */
        ap_uint<DDR_W> *htb_buf2,
        ap_uint<DDR_W> *htb_buf3,
        row_t *res_buf,
        AdjHT *hTables0,
        AdjHT *hTables1,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        unsigned short nQueryVer,

#ifdef COUNT_ONLY
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
#pragma HLS dataflow

    hls::stream<bool, START_BATCH_SIZE> stream_batch_end("Stream batch end");
    hls::stream<bool, S_D> stream_batches_end("Stream batches end");
    hls::stream<ap_uint<V_ID_W>,  START_BATCH_SIZE> stream_batch("Stream batch");

    mwj_batch<START_BATCH_SIZE>(
        hTables0,
        qVertices0,
        htb_buf0,
        stream_batch_end,
        stream_batches_end,
        stream_batch);
  
    do {    
        multiwayJoin(
                htb_buf1,
                htb_buf2,
                htb_buf3,
                res_buf,
                hTables0,
                hTables1,
                qVertices0,
                qVertices1,
                nQueryVer,
                stream_batch,
                stream_batch_end,
                result);
   
#ifdef DEBUG_STATS
        debug::batches++;
#endif

    } while (!stream_batches_end.read());
}

void subgraphIsomorphism(
        edge_t edge_buf[GRAPHS_SPACE],
        row_t htb_buf0[HASHTABLES_SPACE],
        row_t htb_buf1[HASHTABLES_SPACE],
/* hls::burst_maxi<row_t> htb_buf1, */
        row_t htb_buf2[HASHTABLES_SPACE],
        row_t htb_buf3[HASHTABLES_SPACE],
        row_t res_buf[RESULTS_SPACE],
        unsigned short numQueryVert,
        unsigned short numQueryEdges,
        unsigned long numDataEdges,

#ifdef DEBUG_INTERFACE
        volatile unsigned int &debif_endpreprocess,
#endif /* DEBUG_INTERFACE */

#ifdef COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif /* COUNT_ONLY */

        )
{

/* #pragma HLS INTERFACE mode=m_axi port=htb_buf1 bundle=gmem1 depth=DDR_WIDTH */
#pragma HLS INTERFACE mode=m_axi port=htb_buf0 bundle=batch_prep
#pragma HLS INTERFACE mode=m_axi port=htb_buf1 bundle=prop
#pragma HLS INTERFACE mode=m_axi port=htb_buf2 bundle=inter
#pragma HLS INTERFACE mode=m_axi port=htb_buf3 bundle=verif
#pragma HLS INTERFACE mode=m_axi port=res_buf bundle=fifo
#pragma HLS INTERFACE mode=m_axi port=edge_buf bundle=graph
/* #pragma HLS alias ports=htb_buf0,htb_buf2,htb_buf3 distance=0 */
#pragma HLS alias ports=htb_buf0,htb_buf1,htb_buf2,htb_buf3 distance=0

/* #pragma HLS INTERFACE mode=axis port=stream_src */
/* #pragma HLS INTERFACE mode=axis port=stream_dst */
/* #pragma HLS INTERFACE mode=axis port=stream_src_l */
/* #pragma HLS INTERFACE mode=axis port=stream_dst_l */
#pragma HLS INTERFACE mode=s_axilite port=numQueryVert
#pragma HLS INTERFACE mode=s_axilite port=numQueryEdges
#pragma HLS INTERFACE mode=s_axilite port=numDataEdges
#pragma HLS INTERFACE mode=s_axilite port=return

#ifdef DEBUG_INTERFACE
#pragma HLS INTERFACE mode=s_axilite port=debif_endpreprocess
#endif /* DEBUG_INTERFACE */

#ifdef COUNT_ONLY
#pragma HLS INTERFACE mode=s_axilite port=result
#else
#pragma HLS INTERFACE mode=axis port=result
#endif /* COUNT_ONLY */

#ifdef DEBUG_STATS
/* statistic purposes */

    debug::init();

#endif /* DEBUG_STATS */

#ifdef DEBUG_INTERFACE
    debif_endpreprocess = 0;
#endif /* DEBUG_INTERFACE */

    QueryVertex qVertices0[MAX_QV], qVertices1[MAX_QV];
    TableDescriptor tDescriptors[MAX_TB];
    AdjHT hTables0[MAX_TB], hTables1[MAX_TB];
    unsigned short numTables = 0;
    unsigned long localResult = 0;

    preprocess<row_t,
        EDGE_WIDTH,
        COUNTER_WIDTH,
        DDR_BIT,
        VERTEX_WIDTH_BIT,
        LABEL_WIDTH,
        HASH_WIDTH_FIRST,
        HASH_WIDTH_SECOND,
        STREAM_DEPTH,
        HASHTABLES_SPACE,
        MAX_QUERY_VERTICES,
        MAX_TABLES>(
                edge_buf,
                htb_buf0,
                qVertices0,
                qVertices1,
                tDescriptors,
                hTables0,
                hTables1,
                numTables,
                numQueryVert,
                numQueryEdges,
                numDataEdges);

#ifdef DEBUG_INTERFACE
    ap_wait();
    debif_endpreprocess = 1;
    ap_wait();
#endif /* DEBUG_INTERFACE */

    multiwayJoinWrap(
            htb_buf0,
            htb_buf1,
            htb_buf2,
            htb_buf3,
            res_buf,
            hTables0,
            hTables1,
            qVertices0,
            qVertices1,
            numQueryVert,
            localResult);

    result = localResult;

#ifdef DEBUG_STATS
    {
        using namespace debug;
        std::ofstream debof("../../../../stats.txt", std::ofstream::app);

        unsigned long debug_total_reads = findmin_reads +
            readmin_reads + intersect_reads +
            hashtovid_reads + verify_reads;

        unsigned long debug_verify = solution_wrong + 
            solution_correct;

        unsigned long debug_probe = intersect_bit_falsepositive + intersect_bit_falsenegative +
            intersect_bit_truepositive + intersect_bit_truenegative;

        unsigned long debug_sol = intersect_sol_falsepositive + intersect_sol_falsenegative +
            intersect_sol_truepositive + intersect_sol_truenegative;
        
        unsigned int hw1, hw2, cnt;
        hw1 = H_W_1; 
        hw2 = H_W_2; 
        cnt = C_W;

        debof << "DEBUG STATISTICS HW1: " << hw1 << " HW2: " << hw2
            << " CNT: " << cnt << std::endl << std::endl;

        debof << "\tfindmin reads:     " << findmin_reads << "\t" <<
            findmin_reads * 100 / debug_total_reads << "%" << std::endl;
        debof << "\treadmin reads:     " << readmin_reads << "\t" <<
            readmin_reads * 100 / debug_total_reads << "%"<< std::endl;
        debof << "\tintersect reads:   " << intersect_reads << "\t" <<
            intersect_reads * 100 / debug_total_reads << "%"<< std::endl;
        debof << "\tverify reads:      " << verify_reads << "\t" << 
            verify_reads * 100 / debug_total_reads << "%"<< std::endl;
        debof << "\tTOTAL:             " <<  debug_total_reads << "\t100%\n" << 
            std::endl;
        debof << "\tread per embedding: " << debug_total_reads / solution_correct << std::endl;
        debof << "\tindexed tables:    " << indexed_tables << std::endl;
        debof << "\tindexing tables:   " << indexing_tables << std::endl;
        debof << "\n\tsolution wrong:   " << solution_wrong << "\t" <<
            solution_wrong * 100 / debug_verify << "%" << std::endl;
        debof << "\tsolution correct:   " << solution_correct << "\t" <<
            solution_correct * 100 / debug_verify << "%" << std::endl;
        debof << "\tmax collisions:   " << max_collisions << std::endl;
        debof << "\tavg collisions:   " << avg_collisions / numTables << std::endl;
#if VERIFY_CACHE
        debof << "\tmean cache hit 0:   " << cache_hit_0 / batches  << std::endl;
        debof << "\tmean cache reqs 0:   " << cache_req_0  << std::endl;
        debof << "\tmean cache hit 1:   " << cache_hit_1 / batches  << std::endl;
        debof << "\tmean cache reqs 1:   " << cache_req_1 << std::endl << std::endl;
#endif
        debof << "\tdata struct prob: " << std::setprecision(3) << std::endl <<
            "\t\tFP " << (float)intersect_bit_falsepositive / debug_probe << "\t" <<
            intersect_bit_falsepositive << std::endl <<
            "\t\tTP " << (float)intersect_bit_truepositive / debug_probe << "\t" <<
            intersect_bit_truepositive << std::endl <<
            "\t\tFN " << (float)intersect_bit_falsenegative / debug_probe << "\t" <<
            intersect_bit_falsenegative << std::endl <<
            "\t\tTN " << (float)intersect_bit_truenegative / debug_probe << "\t" <<
            intersect_bit_truenegative << std::endl << std::endl;
        debof << "\tsolution prob: " << std::setprecision(3) << std::endl <<
            "\t\tFP " << (float)intersect_sol_falsepositive / debug_sol << "\t" <<
            intersect_sol_falsepositive << std::endl <<
            "\t\tTP " << (float)intersect_sol_truepositive / debug_sol << "\t" <<
            intersect_sol_truepositive << std::endl <<
            "\t\tFN " << (float)intersect_sol_falsenegative / debug_sol << "\t" <<
            intersect_sol_falsenegative << std::endl <<
            "\t\tTN " << (float)intersect_sol_truenegative / debug_sol << "\t" <<
            intersect_sol_truenegative << std::endl << std::endl;
        debof.close();
    }
#endif
}
