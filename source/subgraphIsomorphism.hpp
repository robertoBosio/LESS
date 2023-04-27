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
#define S_D         SUB_STREAM_DEPTH    
#define DDR_WORDS   RES_WIDTH
#define DDR_W       DDR_WORD

#if VERIFY_CACHE
typedef cache< ap_uint<DDR_W>, true, false, 2,
        HASHTABLES_SPACE, 1, 1, 16, false, 32, 1,
        false, 2> cache_type;
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
    ap_uint<DDR_W> ram_row;
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

    return ram_row;

}

#if VERIFY_CACHE
template <size_t W_1,
        size_t W_2,
        size_t SHF,
        unsigned int CACHE_PORT,
        size_t T>
ap_uint<(1UL << T)> read_table(
        ap_uint<W_1> index1,
        ap_uint<W_2> index2,
        cache_type &htb_buf,
        ap_uint<32> start_addr)
{
#pragma HLS inline
/* #pragma HLS function_instantiate variable=cache_port */
    ap_uint<DDR_W> ram_row;
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
    ram_row = htb_buf.get(addr_row, CACHE_PORT);
    ram_row >>= (addr_inrow << T);

    return ram_row;

}
#endif /* VERIFY_CACHE */

void mwj_propose(
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<DDR_W> *htb_buf,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_stop,

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
                debug::set++;
                debug::readmin_reads++;
#endif
            }
            stream_end_out.write(true);
#ifdef DEBUG_STATS
            debug::n_sets++;
#endif
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

void mwj_homomorphism(
        hls::stream<ap_uint<V_ID_W>>            &stream_set_in,
        hls::stream<bool>                       &stream_set_end_in,
        hls::stream<ap_uint<SET_INFO_WIDTH>>    &stream_set_info_in,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_in,
        hls::stream<bool>                       &stream_sol_end_in,
        
        hls::stream<ap_uint<V_ID_W>>            &stream_set_out,
        hls::stream<bool>                       &stream_set_end_out,
        hls::stream<ap_uint<SET_INFO_WIDTH>>    &stream_set_info_out,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_out,
        hls::stream<bool>                       &stream_sol_end_out)
{
    ap_uint<8> curQV {0};
    ap_uint<2> last_set;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool stop, last_sol;

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

    stream_set_info_out.write(stream_set_info_in.read());
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

        last_set = stream_set_end_in.read();
    }
    stream_set_end_out.write(last_set);
}

template <size_t MAX_BATCH_SIZE>
void mwj_batchbuild(
        hls::stream<ap_uint<V_ID_W>>            &stream_set_in,
        hls::stream<bool>                       &stream_set_end_in,
        hls::stream<ap_uint<SET_INFO_WIDTH>>    &stream_set_info_in,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_in,
        hls::stream<bool>                       &stream_sol_end_in,
       
        hls::stream<bool>                       &stream_req,
        hls::stream<ap_uint<V_ID_W>>            &stream_batch_out,
        hls::stream<bool>                       &stream_batch_end_out,
        hls::stream<ap_uint<SET_INFO_WIDTH>>    &stream_set_info_out,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_out,
        hls::stream<bool>                       &stream_sol_end_out)
{
    ap_uint<8> curQV {0};
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<V_ID_W> set_info;
    unsigned int batch_counter {0};
    bool stop, last;
    
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

    set_info = stream_set_info_in.read();
    stream_set_info_out.write(set_info);
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

            for (int g = 0; g < curQV; g++){
                stream_sol_out.write(curEmb[g]);
                stream_sol_end_out.write(false);
            }
            stream_sol_end_out.write(true);
            stream_set_info_out.write(set_info);
        }
    }
    stream_batch_end_out.write(true);

}

template <size_t TUPLE_I,
         size_t BATCH_SIZE_LOG>
void mwj_tuplebuild(
        QueryVertex                             *qVertices,
        hls::stream<ap_uint<V_ID_W>>            &stream_set_in,
        hls::stream<bool>                       &stream_set_end_in,
        hls::stream<ap_uint<SET_INFO_WIDTH>>    &stream_set_info_in,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_in,
        hls::stream<bool>                       &stream_sol_end_in,
        hls::stream<bool>                       &stream_stop,
       
        hls::stream<ap_uint<TUPLE_I>>           &stream_tuple_out,
        hls::stream<bool>                       &stream_tuple_end_out,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_out,
        hls::stream<bool>                       &stream_sol_end_out)
{
    ap_uint<V_ID_W> buffer[(1UL << BATCH_SIZE_LOG)];
    ap_uint<BATCH_SIZE_LOG> buffer_size {0};
    ap_uint<BATCH_SIZE_LOG> buffer_p;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<V_ID_W> set_info;
    ap_uint<V_ID_W> vToVerify;
    ap_uint<TUPLE_I> tuple;
    unsigned char curQV {0};
    bool stop, last;
#pragma HLS bind_storage variable=buffer type=ram_1p impl=bram
   
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
            set_info = stream_set_info_in.read();
            unsigned char cycles = qVertices[curQV].numTablesIndexed;

            if (cycles > 0){
                uint8_t tableIndex = qVertices[curQV].tables_indexed[0];
                uint8_t ivPos = qVertices[curQV].vertex_indexing[0];
                bool bit_min = (set_info.range(7, 0) != tableIndex || 
                        set_info.range(15, 8) != ivPos);

                last = stream_set_end_in.read();

TUPLEBUILD_MAIN_LOOP_FIRST_IT:
                while(!last){
#pragma HLS pipeline II=1 
                    vToVerify = stream_set_in.read();

                    // Building the tuple as 
                    // (VERTEX_DST, VERTEX_SRC, TABLE, POSITION, BIT_LAST, BIT_MIN)
                    tuple.range(V_ID_W - 1, 0) = vToVerify;
                    tuple.range((V_ID_W * 2) - 1, V_ID_W) = curEmb[ivPos];
                    tuple.range((V_ID_W * 2) + 7, V_ID_W * 2) = tableIndex;
                    tuple.range(TUPLE_I - 3, V_ID_W * 2 + 8) = buffer_size;
                    tuple[TUPLE_I - 2] = (qVertices[curQV].numTablesIndexed == 1);
                    tuple[TUPLE_I - 1] = bit_min;
                    stream_tuple_out.write(tuple);
                    stream_tuple_end_out.write(false);

                    buffer[buffer_size++] = vToVerify;
                    last = stream_set_end_in.read();
                }
            }

TUPLEBUILD_EDGE_LOOP_AFTER_IT:
            for (int g = 1; g < cycles; g++){
                uint8_t tableIndex = qVertices[curQV].tables_indexed[g];
                uint8_t ivPos = qVertices[curQV].vertex_indexing[g];
                bool bit_last = (g == cycles - 1);
                bool bit_min = (set_info.range(7, 0) != tableIndex || 
                        set_info.range(15, 8) != ivPos);

TUPLEBUILD_MAIN_LOOP:
                for (int buffer_p = 0; buffer_p < buffer_size; buffer_p++){
#pragma HLS pipeline II=1 
                    vToVerify = buffer[buffer_p];

                    // Building the tuple as 
                    // (VERTEX_DST, VERTEX_SRC, TABLE, POSITION, BIT_LAST, BIT_MIN)
                    tuple.range(V_ID_W - 1, 0) = vToVerify;
                    tuple.range((V_ID_W * 2) - 1, V_ID_W) = curEmb[ivPos];
                    tuple.range((V_ID_W * 2) + 7, V_ID_W * 2) = tableIndex;
                    tuple.range(TUPLE_I - 3, V_ID_W * 2 + 8) = buffer_p;
                    tuple[TUPLE_I - 2] = bit_last;
                    tuple[TUPLE_I - 1] = bit_min;
                    stream_tuple_out.write(tuple);
                    stream_tuple_end_out.write(false);
                }
            }
            stream_tuple_end_out.write(true);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

template <typename T,
            size_t TUPLE_I,
            size_t TUPLE_V,
            size_t BATCH_SIZE_LOG>
void mwj_intersect(
        AdjHT                         *hTables,
        T                             htb_buf,
        hls::stream<ap_uint<TUPLE_I>> &stream_tuple_in,
        hls::stream<bool>             &stream_tuple_end_in,
        hls::stream<ap_uint<V_ID_W>>  &stream_sol_in,
        hls::stream<bool>             &stream_sol_end_in,
        hls::stream<bool>             &stream_stop,

        hls::stream<ap_uint<TUPLE_V>> &stream_tuple_out,
        hls::stream<bool>             &stream_tuple_end_out,
        hls::stream<ap_uint<V_ID_W>>  &stream_sol_out,
        hls::stream<bool>             &stream_sol_end_out)
{
    ap_uint<64> candidate_hash;
    ap_uint<(1UL << BATCH_SIZE_LOG)> bits;
    ap_uint<V_ID_W> candidate_v;
    ap_uint<V_ID_W> indexing_v;
    ap_uint<TUPLE_I> tuple_in;
    ap_uint<TUPLE_V> tuple_out;
    ap_uint<SET_INFO_WIDTH> set_info;
    unsigned char tableIndex;
    ap_uint<DDR_W> ram_row;
    unsigned long addr_row;
    unsigned long addr_inrow;
    ap_uint<64> addr_counter;
    bool stop, last;

    while (1) {
        if (stream_sol_end_in.read_nb(last)){
INTERSECT_COPYING_EMBEDDING_LOOP:
            while(!last){
                stream_sol_out.write(stream_sol_in.read());
                stream_sol_end_out.write(false);
                last = stream_sol_end_in.read();
            }
            stream_sol_end_out.write(true);

            bits = ~0;
            last = stream_tuple_end_in.read();

#ifdef DEBUG_STATS
            if (last)
                debug::empty_sol++;
#endif

INTERSECT_LOOP:
            while(!last){
                tuple_in = stream_tuple_in.read();
                ap_uint<64> hash_out;
                ap_uint<(1UL << C_W)> start_off = 0;
                ap_uint<(1UL << C_W)> end_off = 1;
                ap_uint<BATCH_SIZE_LOG> pos = tuple_in.range(TUPLE_I - 3, V_ID_W * 2 + 8);

                if (tuple_in.test(TUPLE_I - 1) && bits[pos]){

                    ap_uint<H_W_1> index1_f;
                    ap_uint<H_W_2> index2_f;
                    candidate_v = tuple_in.range(V_ID_W - 1, 0);
                    indexing_v = tuple_in.range((V_ID_W * 2) - 1, V_ID_W);
                    tableIndex = tuple_in.range((V_ID_W * 2) + 7, V_ID_W * 2);

                    xf::database::details::hashlookup3_core<V_ID_W>(
                            candidate_v, candidate_hash);
                    xf::database::details::hashlookup3_core<V_ID_W>(
                            indexing_v, hash_out);
                    ap_uint<H_W_1> index1 = hash_out.range(H_W_1 - 1, 0);
                    ap_uint<H_W_2> index2 = candidate_hash.range(H_W_2 - 1, 0);

                    addr_counter = index1;
                    addr_counter <<= H_W_2;
                    addr_counter += index2;

                    /* Compute address of row storing the counter */
                    addr_row = hTables[tableIndex].start_offset + 
                        (addr_counter >> (DDR_BIT - C_W));

                    /* Compute address of data inside the row */
                    addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);
                    ram_row = htb_buf.get(addr_row, 0);
                    end_off = ram_row.range(((addr_inrow + 1) << C_W) - 1, 
                            addr_inrow << C_W);
                    
                    if (addr_counter != 0){
                        addr_counter--;

                        /* Compute address of row storing the counter */
                        addr_row = hTables[tableIndex].start_offset + 
                            (addr_counter >> (DDR_BIT - C_W));

                        /* Compute address of data inside the row */
                        addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);
                        ram_row = htb_buf.get(addr_row, 0);
                        start_off = ram_row.range(((addr_inrow + 1) << C_W) - 1, 
                                addr_inrow << C_W);
#ifdef DEBUG_STATS
                        debug::intersect_reads++;
#endif
                    }

#ifdef DEBUG_STATS
                    debug::intersect_reads++;
                    debug::intersect_filter += (start_off < end_off)? 1 : 0;
#endif
                }

                // Building the tuple as 
                // ((TUPLE_I), START, END, INTER_BIT)
                tuple_out[TUPLE_V - 1] = (start_off < end_off) & bits[pos];
                tuple_out.range(TUPLE_I - 1, 0) = tuple_in;
                tuple_out.range(TUPLE_I + (1UL << C_W) - 1, TUPLE_I) = start_off;
                tuple_out.range(TUPLE_I + 2 * (1UL << C_W) - 1,
                       TUPLE_I + (1UL << C_W)) = end_off;

                bits[pos] = bits[pos] & (start_off < end_off);
                stream_tuple_out.write(tuple_out);
                stream_tuple_end_out.write(false);

/* std::cout << "(( " */
/* << tuple_out.range(V_ID_W - 1, 0) << ", " */
/* << tuple_out.range((V_ID_W * 2) - 1, V_ID_W) << ", " */
/* << tuple_out.range((V_ID_W * 2) + 7, V_ID_W * 2) << ", " */
/* << tuple_out.range(TUPLE_I - 3, V_ID_W * 2 + 8) << ", " */
/* << tuple_out[TUPLE_I - 2] << ", " */
/* << tuple_out[TUPLE_I - 1] << "), " */
/* << tuple_out.range(TUPLE_I + (1UL << C_W) - 1, TUPLE_I) << ", " */
/* << tuple_out.range(TUPLE_I + 2 * (1UL << C_W) - 1, */
/* TUPLE_I + (1UL << C_W)) << ", " */
/* << tuple_out[TUPLE_V - 1] << ") " */
/* << bits[pos] << std::endl; */

                last = stream_tuple_end_in.read();
            }
            stream_tuple_end_out.write(true);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

template <typename T,
            size_t TUPLE_I,
            size_t TUPLE_V,
            size_t TUPLE_A,
            size_t BATCH_SIZE_LOG>
void mwj_verify(
        AdjHT *hTables,
        T htb_buf,
        hls::stream<ap_uint<TUPLE_V>>           &stream_tuple_in,
        hls::stream<bool>                       &stream_tuple_end_in,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_in,
        hls::stream<bool>                       &stream_sol_end_in,
        hls::stream<bool>                       &stream_stop,

        hls::stream<ap_uint<TUPLE_A>>           &stream_tuple_out,
        hls::stream<bool>                       &stream_tuple_end_out,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_out,
        hls::stream<bool>                       &stream_sol_end_out)
{
    ap_uint<(1UL << BATCH_SIZE_LOG)> bits;
    ap_uint<(1UL << C_W)> start_off;
    ap_uint<(1UL << C_W)> end_off;
    ap_uint<TUPLE_V> tuple_in;
    ap_uint<TUPLE_A> tuple_out;
    ap_uint<V_ID_W> candidate_v;
    ap_uint<V_ID_W> indexing_v;
    unsigned char tableIndex;
    bool stop, last;

    while(1){
        if (stream_sol_end_in.read_nb(last)){

VERIFY_EDGE_COPYING_EMBEDDING_LOOP:
            while(!last){
                stream_sol_out.write(stream_sol_in.read());
                stream_sol_end_out.write(false);
                last = stream_sol_end_in.read();
            }
            stream_sol_end_out.write(true);

            bits = ~0;
            last = stream_tuple_end_in.read();
VERIFY_MAIN_LOOP:
            while(!last){
                tuple_in = stream_tuple_in.read();
                ap_uint<(1UL << C_W)> start_off = 0;
                ap_uint<(1UL << C_W)> end_off = 0;
                ap_uint<BATCH_SIZE_LOG> pos = tuple_in.range(TUPLE_I - 3, V_ID_W * 2 + 8);
                candidate_v = tuple_in.range(V_ID_W - 1, 0);
                bool checked = !tuple_in.test(TUPLE_I - 1);
                if (tuple_in.test(TUPLE_V - 1) && tuple_in.test(TUPLE_I - 1) && bits[pos]){
                    indexing_v = tuple_in.range((V_ID_W * 2) - 1, V_ID_W);
                    tableIndex = tuple_in.range((V_ID_W * 2) + 7, V_ID_W * 2);
                    start_off = tuple_in.range(TUPLE_I + (1UL << C_W) - 1, TUPLE_I);
                    end_off = tuple_in.range(TUPLE_I + 2 * (1UL << C_W) - 1,
                            TUPLE_I + (1UL << C_W));
        
VERIFY_READ_MEMORY_LOOP:
                    for (; start_off < end_off; start_off++){
                        ap_uint<2 * V_ID_W> edge = read_table<32, 1, 0, 1, E_W>(
                                start_off,
                                0,
                                htb_buf,
                                hTables[tableIndex].start_edges);

                        ap_uint<V_ID_W> vertexIndexed, vertexIndexing;
                        vertexIndexed = edge.range(V_ID_W - 1, 0);
                        vertexIndexing = edge.range(2 * V_ID_W - 1, V_ID_W);
                        if (vertexIndexing == indexing_v &&
                                vertexIndexed == candidate_v){
                            checked = true;
                        }
#ifdef DEBUG_STATS
                        debug::verify_reads++;
                        /* Computing when there is an alias */

                        if (tuple_in.test(TUPLE_V - 1) && checked) {
                            /* // True positive */
                            debug::intersect_bit_truepositive++;
                        } else if (tuple_in.test(TUPLE_V - 1) && !checked) {
                            /* // False positive */
                            debug::intersect_bit_falsepositive++;
                        }
#endif
                    }
                }

                // Building the tuple as 
                // (VERTEX_DST, POSITION, LAST_BIT, EDGE_BIT)
                // example: (v1, 107, 1)
                tuple_out[TUPLE_A - 1] = checked & bits[pos];
                tuple_out[TUPLE_A - 2] = tuple_in[TUPLE_I - 2];
                tuple_out.range(V_ID_W - 1, 0) = candidate_v;
                tuple_out.range(TUPLE_A - 3, V_ID_W) = 
                    tuple_in(TUPLE_I - 3, V_ID_W * 2 + 8);
                bits[pos] = bits[pos] & checked; 
/* std::cout << "( " */
/* << tuple_out.range(V_ID_W - 1, 0) << ", " */
/* << tuple_out.range(TUPLE_A - 3, V_ID_W) << ", " */
/* << tuple_out[TUPLE_A - 2] << ", " */
/* << tuple_out[TUPLE_A - 1] << ")" << std::endl; */
                
                stream_tuple_out.write(tuple_out);
                stream_tuple_end_out.write(false);
                last = stream_tuple_end_in.read();
            }
            stream_tuple_end_out.write(true);
        }

        if (stream_stop.read_nb(stop)){
            break;
        }
    }
}

template<size_t TUPLE_A,
        size_t MAX_BATCH_SIZE>
void mwj_filter(
        hls::stream<ap_uint<TUPLE_A>>           &stream_tuple_in,
        hls::stream<bool>                       &stream_tuple_end_in,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_in,
        hls::stream<bool>                       &stream_sol_end_in,
        
        hls::stream<ap_uint<V_ID_W>>            &stream_set_out,
        hls::stream<bool>                       &stream_set_end_out,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_out,
        hls::stream<bool>                       &stream_sol_end_out)
{
    ap_uint<MAX_BATCH_SIZE> bits;
    ap_uint<TUPLE_A> tuple_in;
    bool last;

    last = stream_sol_end_in.read();
    while(!last){
        stream_sol_out.write(stream_sol_in.read());
        stream_sol_end_out.write(false);
        last = stream_sol_end_in.read();
    }
    stream_sol_end_out.write(true);
    
    bits = ~0;
    last = stream_tuple_end_in.read();
FILTER_LOOP:
    while(!last){
        tuple_in = stream_tuple_in.read();
        unsigned short p = tuple_in.range(TUPLE_A - 3, V_ID_W);
        bits[p] = bits[p] & tuple_in[TUPLE_A - 1];
        if (tuple_in.test(TUPLE_A - 2) && bits.test(p)){
/* std::cout << tuple_in.range(V_ID_W - 1, 0) << std::endl; */
            stream_set_out.write(tuple_in.range(V_ID_W - 1, 0));
            stream_set_end_out.write(false);
        }
        last = stream_tuple_end_in.read();
    }
    stream_set_end_out.write(true);
}

template<typename FLAG_T>
void mwj_garbagecollector(
        hls::stream<ap_uint<V_ID_W>>            &stream_set_in,
        hls::stream<FLAG_T>                     &stream_set_end_in,
        hls::stream<ap_uint<SET_INFO_WIDTH>>    &stream_set_info_in,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_in,
        hls::stream<bool>                       &stream_sol_end_in,
        
        hls::stream<bool>                       &stream_req, 
        hls::stream<ap_uint<V_ID_W>>            &stream_set_out,
        hls::stream<FLAG_T>                     &stream_set_end_out,
        hls::stream<ap_uint<SET_INFO_WIDTH>>    &stream_set_info_out,
        hls::stream<ap_uint<V_ID_W>>            &stream_sol_out,
        hls::stream<bool>                       &stream_sol_end_out)
{
    ap_uint<8> curQV {0};
    ap_uint<SET_INFO_WIDTH> set_info;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool last;
    FLAG_T last_set;

    last = stream_sol_end_in.read();

GARBAGECOLLECTOR_READING_SOL_LOOP:
    while(!last){
        curEmb[curQV++] = stream_sol_in.read();
        last = stream_sol_end_in.read();
    }

    set_info = stream_set_info_in.read();
    last_set = stream_set_end_in.read();

    if (!last_set){
GARBAGECOLLECTOR_WRITING_SOL_LOOP:
        for (int g = 0; g < curQV; g++){
            stream_sol_out.write(curEmb[g]);
            stream_sol_end_out.write(false);
        }
        stream_sol_end_out.write(true);
        stream_set_info_out.write(set_info);
        stream_set_end_out.write(last_set);
    } else {
        stream_req.write(false);
    }

GARBAGECOLLECTOR_WRITING_SET_LOOP:
    while(!last_set){
        stream_set_out.write(stream_set_in.read());
        last_set = stream_set_end_in.read();
        stream_set_end_out.write(last_set);
    }
}

void mwj_assembly(
        unsigned short nQueryVer,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<V_ID_W>> &stream_batch,
        hls::stream<bool> &stream_batch_end,
       
        hls::stream<bool> &stream_stop,
        hls::stream<bool> &stream_req, 
        hls::stream<ap_uint<V_ID_W>> &stream_partial_out,
#ifdef COUNT_ONLY
        long unsigned int &result
#else  
        hls::stream<T_NODE> &result
#endif
        )
{
    ap_uint<V_ID_W> curQV;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<SET_INFO_WIDTH> setinfo;
    unsigned long int nPartSol {0};
    bool last_sol, last_set, stop;
    T_NODE node;

#ifdef COUNT_ONLY
    unsigned long int counter {0};
#endif

    last_sol = stream_batch_end.read();
    curQV = (1UL << (V_ID_W - 1)) + 1;
    stream_partial_out.write(curQV);
    while(!last_sol){
        stream_partial_out.write(stream_batch.read());
        stream_req.write(true);
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

            if (!last_set && (curQV != nQueryVer - 1)){
                stream_partial_out.write((curQV + 1) | (1UL << (V_ID_W - 1)));
                for (int g = 0; g < curQV; g++){
                    stream_partial_out.write(curEmb[g]);
                }
            }

VERIFY_CHECK_LOOP:
            while(!last_set){
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
                    stream_req.write(true);
                }
#ifdef DEBUG_STATS
                debug::embeddings++;
#endif
                last_set = stream_end_inter_in.read();
            }

            // Last batch of a set 
            stream_req.write(false);

        }

        if (stream_stop.read_nb(stop))
            break;
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
        hls::stream<bool> streams_stop[STOP_S])
{
    static unsigned long sol {0};
    bool req = stream_req.read();
    
    if (req) {
        sol++;
    } else {
        sol--;
    }

    if (sol == 0){
        for (int g = 0; g < STOP_S; g++){
#pragma HLS unroll
            streams_stop[g].write(true);
        }
    }
}

template <size_t TUPLE_I,
            size_t TUPLE_V,
            size_t TUPLE_A,
            size_t MAX_BATCH_LOG>
void mwj_Wrapper(
        AdjHT *hTables,
        cache_type &htb_buf,
        hls::stream<ap_uint<TUPLE_I>> &stream_tuple_in,
        hls::stream<bool>             &stream_tuple_end_in,
        hls::stream<ap_uint<V_ID_W>>  &stream_sol_in,
        hls::stream<bool>             &stream_sol_end_in,
        hls::stream<bool>             &stream_stop_0,
        hls::stream<bool>             &stream_stop_1,
        
        hls::stream<ap_uint<TUPLE_A>> &stream_tuple_out,
        hls::stream<bool>             &stream_tuple_end_out,
        hls::stream<ap_uint<V_ID_W>>  &stream_sol_out,
        hls::stream<bool>             &stream_sol_end_out)
{
#pragma HLS dataflow
    /* Intersect data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> i_stream_sol
        ("Intersect - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> i_stream_sol_end
        ("Intersect - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<TUPLE_V>, S_D> i_stream_tuple
        ("Intersect - tuples");
    hls_thread_local hls::stream<bool, S_D> i_stream_tuple_end
        ("Intersect - tuples end flag");

#ifdef __SYNTHESIS__
    mwj_intersect<cache_type&, 
        TUPLE_I,
        TUPLE_V,
        MAX_BATCH_LOG>(
                hTables,
                htb_buf,
                stream_tuple_in,
                stream_tuple_end_in,
                stream_sol_in,
                stream_sol_end_in,
                stream_stop_0,
                i_stream_tuple,
                i_stream_tuple_end,
                i_stream_sol,
                i_stream_sol_end);
            
    mwj_verify<cache_type&, 
        TUPLE_I,
        TUPLE_V,
        TUPLE_A,
        MAX_BATCH_LOG>(
                hTables,
                htb_buf,
                i_stream_tuple,
                i_stream_tuple_end,
                i_stream_sol,
                i_stream_sol_end,
                stream_stop_1,
                stream_tuple_out,
                stream_tuple_end_out,
                stream_sol_out,
                stream_sol_end_out);
#else

    htb_buf.init();
    std::thread mwj_intersect_t(
            mwj_intersect<cache_type&, 
        TUPLE_I,
        TUPLE_V,
        PROPOSE_BATCH_LOG>,
            hTables,
            std::ref(htb_buf),
            std::ref(stream_tuple_in),
            std::ref(stream_tuple_end_in),
            std::ref(stream_sol_in),
            std::ref(stream_sol_end_in),
            std::ref(stream_stop_0),
            std::ref(i_stream_tuple),
            std::ref(i_stream_tuple_end),
            std::ref(i_stream_sol),
            std::ref(i_stream_sol_end));
    
    std::thread mwj_verify_t(
            mwj_verify<cache_type&, 
            TUPLE_I,
            TUPLE_V,
            TUPLE_A,
            MAX_BATCH_LOG>,
            hTables,
            std::ref(htb_buf),
            std::ref(i_stream_tuple),
            std::ref(i_stream_tuple_end),
            std::ref(i_stream_sol),
            std::ref(i_stream_sol_end),
            std::ref(stream_stop_1),
            std::ref(stream_tuple_out),
            std::ref(stream_tuple_end_out),
            std::ref(stream_sol_out),
            std::ref(stream_sol_end_out));
    
    mwj_intersect_t.join();
    mwj_verify_t.join();
    htb_buf.stop();

#endif
}

void mwj_batch(
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<DDR_W> *htb_buf,
        unsigned short numBatchSize,
        
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
    unsigned int batch_counter = numBatchSize + 1;

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
                    if(batch_counter == numBatchSize){
                        stream_batches_end.write(false);
                        stream_batch_end.write(true);
                    }

                    set[set_counter++] = vertex;
                    batch_counter = (batch_counter >= numBatchSize)? 1 : batch_counter + 1;
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
/* hls::burst_maxi<row_t> htb_buf0, */
        ap_uint<DDR_W> *htb_buf0,
        ap_uint<DDR_W> *htb_buf1,
        row_t *res_buf,
        AdjHT *hTables0,
        AdjHT *hTables1,
        QueryVertex *qVertices0,
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
#pragma HLS STABLE variable=hTables0
#pragma HLS STABLE variable=hTables1
#pragma HLS STABLE variable=qVertices0
#pragma HLS STABLE variable=nQueryVer
#pragma HLS DATAFLOW

    constexpr size_t TUPLE_I = (V_ID_W * 2) + 8 + PROPOSE_BATCH_LOG + 2;
    constexpr size_t TUPLE_V = TUPLE_I + ((1UL << COUNTER_WIDTH) * 2) + 1;
    constexpr size_t TUPLE_A = V_ID_W + PROPOSE_BATCH_LOG + 2;

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
    
    /* Batchbuild data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> b_stream_sol
        ("Batchbuild - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> b_stream_sol_end
        ("Batchbuild - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> b_stream_set
        ("Batchbuild - set nodes");
    hls_thread_local hls::stream<bool, S_D> b_stream_set_end
        ("Batchbuild - set nodes end flag");
    hls_thread_local hls::stream<ap_uint<SET_INFO_WIDTH>, S_D> b_stream_set_info
        ("Batchbuild - set info");
    
    /* Tuplebuild data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> t_stream_sol
        ("Tuplebuild - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> t_stream_sol_end
        ("Tuplebuild - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<TUPLE_I>, S_D> t_stream_tuple
        ("Tuplebuild - tuples");
    hls_thread_local hls::stream<bool, S_D> t_stream_tuple_end
        ("Tuplebuild - tuples end flag");
    
    /* Verify data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> v_stream_sol
        ("Verify - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> v_stream_sol_end
        ("Verify - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<TUPLE_A>, S_D> v_stream_tuple
        ("Verify - tuples");
    hls_thread_local hls::stream<bool, S_D> v_stream_tuple_end
        ("Verify - tuples end flag");

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
    hls_thread_local hls::stream<ap_uint<V_ID_W>, 32> a_stream_sol
        ("Assembly - partial solution");

    /* Dynamic fifo data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, 32> dyn_stream_sol
        ("Dynamic fifo - partial solution");
 
    /* Stop signals */
    hls_thread_local hls::stream<bool, 1> streams_stop[STOP_S];
    hls_thread_local hls::stream<bool, S_D> merge_out;
    hls_thread_local hls::stream<bool, S_D> merge_in[MERGE_IN_STREAMS];
#pragma HLS array_partition variable=merge_in type=complete    

#ifndef __SYNTHESIS__
    for (int g = 0; g < STOP_S; g++) {
        char stream_name[10];
        sprintf(stream_name, "stop_%d", g);
        streams_stop[g].set_name(stream_name);
    } 
#endif /* __SYNTHESIS__ */

#if VERIFY_CACHE
    cache_type a_cache(htb_buf1);
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

    hls_thread_local hls::task mwj_homomorphism_t(
            mwj_homomorphism,
            p_stream_set,
            p_stream_set_end,
            p_stream_set_info,
            p_stream_sol,
            p_stream_sol_end,
            h_stream_set,
            h_stream_set_end,
            h_stream_set_info,
            h_stream_sol,
            h_stream_sol_end);

    hls_thread_local hls::task mwj_batchbuild_t(
            mwj_batchbuild<(1UL << PROPOSE_BATCH_LOG)>,
            h_stream_set,
            h_stream_set_end,
            h_stream_set_info,
            h_stream_sol,
            h_stream_sol_end,
            merge_in[0],
            b_stream_set,
            b_stream_set_end,
            b_stream_set_info,
            b_stream_sol,
            b_stream_sol_end);


    hls_thread_local hls::task mwj_filter_t(
            mwj_filter<TUPLE_A,
            (1UL << PROPOSE_BATCH_LOG)>,
            v_stream_tuple,
            v_stream_tuple_end,
            v_stream_sol,
            v_stream_sol_end,
            f_stream_set,
            f_stream_set_end,
            f_stream_sol,
            f_stream_sol_end);

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

    mwj_tuplebuild<TUPLE_I,
        PROPOSE_BATCH_LOG>(
                qVertices0,
                b_stream_set,
                b_stream_set_end,
                b_stream_set_info,
                b_stream_sol,
                b_stream_sol_end,
                streams_stop[4],
                t_stream_tuple,
                t_stream_tuple_end,
                t_stream_sol,
                t_stream_sol_end);

    cache_wrapper (
            mwj_Wrapper<TUPLE_I,
            TUPLE_V,
            TUPLE_A,
            PROPOSE_BATCH_LOG>,
            hTables1,
            a_cache,
            t_stream_tuple,
            t_stream_tuple_end,
            t_stream_sol,
            t_stream_sol_end,
            streams_stop[1],
            streams_stop[2],
            v_stream_tuple,
            v_stream_tuple_end,
            v_stream_sol,
            v_stream_sol_end);

    mwj_assembly(
            nQueryVer,
            f_stream_set,
            f_stream_set_end,
            f_stream_sol,
            f_stream_sol_end,
            stream_batch,
            stream_batch_end,
            streams_stop[3],
            merge_in[1],
            a_stream_sol,
            result);
    
    hls_thread_local hls::task mwj_merge_t(
            mwj_merge<bool>,
            merge_in,
            merge_out);

    hls_thread_local hls::task mwj_stop_t(
            mwj_stop,
            merge_out,
            streams_stop);

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
    
    std::thread mwj_tuplebuild_t(
            mwj_tuplebuild<TUPLE_I,
            PROPOSE_BATCH_LOG>,
            qVertices0,
            std::ref(b_stream_set),
            std::ref(b_stream_set_end),
            std::ref(b_stream_set_info),
            std::ref(b_stream_sol),
            std::ref(b_stream_sol_end),
            std::ref(streams_stop[4]),
            std::ref(t_stream_tuple),
            std::ref(t_stream_tuple_end),
            std::ref(t_stream_sol),
            std::ref(t_stream_sol_end));

    std::thread mwj_assembly_t(
            mwj_assembly,
            nQueryVer,
            std::ref(f_stream_set),
            std::ref(f_stream_set_end),
            std::ref(f_stream_sol),
            std::ref(f_stream_sol_end),
            std::ref(stream_batch),
            std::ref(stream_batch_end),
            std::ref(streams_stop[3]),
            std::ref(merge_in[0]),
            std::ref(a_stream_sol),
            std::ref(result));

    hls_thread_local hls::task mwj_merge_t(
            mwj_merge<bool>,
            merge_in,
            merge_out);

    hls_thread_local hls::task mwj_stop_t(
            mwj_stop,
            merge_out,
            streams_stop);

    mwj_Wrapper<TUPLE_I,
        TUPLE_V,
        TUPLE_A,
        PROPOSE_BATCH_LOG>(
            hTables1,
            a_cache,
            t_stream_tuple,
            t_stream_tuple_end,
            t_stream_sol,
            t_stream_sol_end,
            streams_stop[1],
            streams_stop[2],
            v_stream_tuple,
            v_stream_tuple_end,
            v_stream_sol,
            v_stream_sol_end);

    mwj_assembly_t.join();
    mwj_propose_t.join();
    mwj_tuplebuild_t.join();

#ifdef DEBUG_STATS
    debug::cache_hit_0 += a_cache.get_hit_ratio(0);
    debug::cache_hit_1 += a_cache.get_hit_ratio(1);
    debug::cache_req_0 += a_cache.get_n_l1_reqs(0);
    debug::cache_req_1 += a_cache.get_n_l1_reqs(1);
#endif /* DEBUG_STATS */

#endif /* __SYNTHESIS__ */

}

template <size_t MAX_BATCH_SIZE>
void multiwayJoinWrap(
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
        unsigned short numBatchSize,

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
#pragma HLS dataflow

    hls::stream<bool, MAX_BATCH_SIZE> stream_batch_end("Stream batch end");
    hls::stream<bool, S_D> stream_batches_end("Stream batches end");
    hls::stream<ap_uint<V_ID_W>,  MAX_BATCH_SIZE> stream_batch("Stream batch");

    mwj_batch(
        hTables0,
        qVertices0,
        htb_buf0,
        numBatchSize,
        stream_batch_end,
        stream_batches_end,
        stream_batch);
 
MULTIWAYJOIN_LOOP: 
    do {    
        multiwayJoin(
                htb_buf1,
                htb_buf2,
                res_buf,
                hTables0,
                hTables1,
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
        row_t res_buf[RESULTS_SPACE],
        unsigned short numQueryVert,
        unsigned short numQueryEdges,
        unsigned long numDataEdges,
        unsigned short numBatchSize,

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

#pragma HLS INTERFACE mode=m_axi port=htb_buf0 bundle=batch_prep
#pragma HLS INTERFACE mode=m_axi port=htb_buf1 bundle=prop
#pragma HLS INTERFACE mode=m_axi port=htb_buf2 bundle=cache
#pragma HLS INTERFACE mode=m_axi port=res_buf bundle=fifo
#pragma HLS INTERFACE mode=m_axi port=edge_buf bundle=graph
#pragma HLS alias ports=htb_buf0,htb_buf1,htb_buf2 distance=0

#pragma HLS INTERFACE mode=s_axilite port=numQueryVert
#pragma HLS INTERFACE mode=s_axilite port=numQueryEdges
#pragma HLS INTERFACE mode=s_axilite port=numDataEdges
#pragma HLS INTERFACE mode=s_axilite port=numBatchSize
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
        SUB_STREAM_DEPTH,
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

    multiwayJoinWrap<
        MAX_START_BATCH_SIZE>(
            htb_buf0,
            htb_buf1,
            htb_buf2,
            res_buf,
            hTables0,
            hTables1,
            qVertices0,
            qVertices1,
            numQueryVert,
            numBatchSize,
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

        unsigned long debug_probe = intersect_bit_falsepositive +
            intersect_bit_truepositive;
        
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
        debof << "\tread per embedding: " << debug_total_reads / result << std::endl;
        debof << "\tindexed tables:    " << indexed_tables << std::endl;
        debof << "\tindexing tables:   " << indexing_tables << std::endl;
/* debof << "\n\tsolution wrong:   " << solution_wrong << "\t" << */
/* solution_wrong * 100 / debug_verify << "%" << std::endl; */
/* debof << "\tsolution correct:   " << solution_correct << "\t" << */
/* solution_correct * 100 / debug_verify << "%" << std::endl; */
        debof << "\tmax collisions:   " << max_collisions << std::endl;
        debof << "\tavg collisions:   " << avg_collisions / numTables << std::endl;
        debof << "\tintersect filter:   " << intersect_filter << std::endl;
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
            intersect_bit_truepositive << std::endl << std::endl;
        debof.close();
    }
    std::cout << debug::set / (float)debug::n_sets << std::endl;
#endif
}
