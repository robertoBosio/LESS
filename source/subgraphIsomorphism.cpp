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

#define STOP_S      7    
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

#if CACHE_ENABLE
typedef cache< ap_uint<DDR_W>, true, false, 2,
        HASHTABLES_SPACE, 1, 1, (1UL << CACHE_WORDS_PER_LINE), false, 128, 1,
        false, 1, AUTO, BRAM> htb_cache_t;

typedef cache< ap_uint<DDR_W>, true, false, 1,
        HASHTABLES_SPACE, 128, 1, (1UL << K_FUNCTIONS), false, 0, 0,
        false, 4, BRAM, AUTO> bloom_cache_t;
#endif /* CACHE_ENABLE */


/******** Tuple definition ********/
typedef struct {
    ap_uint<V_ID_W>         indexing_v;
    unsigned char           tb_index;
    unsigned char           iv_pos;     //Query indexing vertex position.
} readmin_tuple_t;

typedef struct {
    unsigned char           tb_index;
    unsigned char           iv_pos;     //Query indexing vertex position.
} tuplebuild_tuple_t;

typedef struct {
    ap_uint<V_ID_W>         indexing_v;
    ap_uint<V_ID_W>         indexed_v;
    ap_uint<64>             addr_counter;
    unsigned char           tb_index;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    bool                    bit_last_edge;
    bool                    bit_min_set;
} intersect_tuple_t;

typedef struct {
    ap_uint<V_ID_W>         indexing_v;
    ap_uint<V_ID_W>         indexed_v;
    unsigned char           tb_index;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    ap_uint<(1UL << C_W)>   start_off;
    unsigned short          n_edges;
    bool                    bit_last_edge;
    bool                    bit_min_set;
    bool                    bit_no_edge;
} split_tuple_t;

typedef struct {
    ap_uint<V_ID_W>         indexing_v;
    ap_uint<V_ID_W>         indexed_v;
    unsigned char           tb_index;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    ap_uint<(1UL << C_W)>   address;
    bool                    bit_last_address;
    bool                    bit_last_edge;
    bool                    bit_min_set;
    bool                    bit_no_edge;
} verify_tuple_t;

typedef struct {
    ap_uint<V_ID_W>         indexed_v;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    bool                    bit_last_address;
    bool                    bit_last_edge;
    bool                    bit_equal;
} compact_tuple_t;

typedef struct {
    ap_uint<V_ID_W>         indexed_v;
    ap_uint<PROPOSE_BATCH_LOG> pos;
    bool                    bit_last_edge;
    bool                    bit_checked;
} assembly_tuple_t;

/******** End tuple definition ********/

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

template <typename T_MEM,
        size_t W_1,
        size_t W_2,
        size_t SHF,
        size_t T,
        size_t PORT>
ap_uint<(1UL << T)> read_table(
        ap_uint<W_1> index1,
        ap_uint<W_2> index2,
        T_MEM htb_buf,
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
    ram_row = htb_buf.get(addr_row, PORT);
    ram_row >>= (addr_inrow << T);

    return ram_row;

}

template <typename T_BLOOM,
         size_t BLOOM_LOG,
        size_t K_FUN_LOG>
unsigned int bloom_bitset(T_BLOOM filter)                         
{
#pragma HLS inline
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

template <typename T_BLOOM,
        size_t BLOOM_LOG, 
        size_t K_FUN_LOG>
void mwj_findmin(
        const unsigned char             hash1_w,
        QueryVertex                     *qVertices,
        bloom_cache_t                   &bloom_p,
        hls::stream< ap_uint<V_ID_W> >  &stream_embed_in,
        hls::stream<bool>               &stream_stop,

        hls::stream<readmin_tuple_t>    &stream_tuple_out,
        hls::stream<T_BLOOM>            &stream_filter_out,
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_out,
        hls::stream<bool>               &stream_sol_end_out)
{
    constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    T_BLOOM filter[K_FUN];
#pragma HLS array_partition variable=filter type=complete dim=1

    readmin_tuple_t tuple;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<V_ID_W> curQV;
    ap_uint<V_ID_W> readv;
    unsigned int minSize;
    unsigned int minTab; 
    unsigned int minIvPos;
    bool stop;

FINDMIN_TASK_LOOP:
    while (1) {

        if (stream_embed_in.read_nb(readv)){

            //Initializing variables for a new partial solution
            for (int g = 0; g < K_FUN; filter[g++] = ~0){
#pragma HLS unroll
            }
            minSize = ~0;
            
            if (readv.test(V_ID_W - 1)){
                // Delimiter, new solution
                curQV = readv.range(V_ID_W - 2, 0);
FINDMIN_NEW_SOLUTION_LOOP:
                for (int g = 0; g < curQV; g++){
#pragma HLS pipeline II=1
                    curEmb[g] = stream_embed_in.read();
                }
            } else {
                //Updating only the last vertex
                curEmb[curQV - 1] = readv;
            }

FINDMIN_COPYING_EMBEDDING_LOOP:
            for (int g = 0; g < curQV; g++) {
#pragma HLS pipeline II=1
                stream_sol_out.write(curEmb[g]);
                stream_sol_end_out.write(false);
            }
            stream_sol_end_out.write(true);

            /* Find sizes of sets in which the current query vertex
             * is indexed by an other query vertex */
FINDMIN_TBINDEXED_LOOP:
            for(int g = 0; g < qVertices[curQV].numTablesIndexed; g++){
#pragma HLS pipeline II=2
                unsigned char tableIndex = qVertices[curQV].tables_indexed[g];
                unsigned char ivPos = qVertices[curQV].vertex_indexing[g];

                //Computing intersection between bloom filters
                ap_uint<64> hash_out;
                xf::database::details::hashlookup3_core<V_ID_W>(curEmb[ivPos], hash_out);
                hash_out = hash_out.range(hash1_w - 1, 0);
                unsigned int address = (tableIndex * (1UL << hash1_w)) + hash_out;
                address <<= K_FUN_LOG;
                unsigned short bloom_s = 0;
                T_BLOOM set_bloom[K_FUN]; 
                bloom_p.get_line(address, 0, set_bloom);
                for (int g = 0; g < K_FUN; g++){
#pragma HLS unroll
                    bloom_s += bloom_bitset<T_BLOOM, BLOOM_LOG, 0>(set_bloom[g]);
                    filter[g] = filter[g] & set_bloom[g];
                }
                bloom_s >>= K_FUN_LOG;

#if DEBUG_STATS
                debug::findmin_reads++;
#endif
                if (bloom_s < minSize){
                    minSize = bloom_s;
                    minTab = tableIndex;
                    minIvPos = ivPos;
                }
            }
          
            tuple.indexing_v = curEmb[minIvPos]; 
            tuple.tb_index = minTab;
            tuple.iv_pos = minIvPos; 
            stream_tuple_out.write(tuple);
            
            for (int g = 0; g < K_FUN; g++){
#pragma HLS unroll
                stream_filter_out.write(filter[g]);
            }
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

template <typename T_BLOOM,
        size_t BLOOM_LOG, 
        size_t K_FUN_LOG>
void mwj_readmin(
        const unsigned char             hash1_w,
        const unsigned char             hash2_w,
        AdjHT                           *hTables,
        row_t                           *m_axi,
        hls::stream<readmin_tuple_t>    &stream_tuple_in,
        hls::stream<T_BLOOM>            &stream_filter_in,
        hls::stream<bool>               &stream_stop,

        hls::stream< ap_uint<V_ID_W> >  &stream_set_out,
        hls::stream<bool>               &stream_set_end_out,
        hls::stream<tuplebuild_tuple_t> &stream_tuple_out)
{
    constexpr size_t K_FUN = (1UL << K_FUN_LOG);
    T_BLOOM filter[K_FUN];
#pragma HLS array_partition variable=filter type=complete dim=1

    // const unsigned long hash_mask = (1UL << hash1_w) - 1; 
    readmin_tuple_t tuple_in;
    tuplebuild_tuple_t tuple_out;
    ap_uint<V_ID_W> vertexCheck;
    ap_uint<V_ID_W> curQV;
    ap_uint<V_ID_W> vertex; 
    ap_uint<V_ID_W * 2> edge;
    bool stop;

READMIN_TASK_LOOP:
    while (true) {
        if (stream_tuple_in.read_nb(tuple_in)){ 

            tuple_out.tb_index = tuple_in.tb_index;
            tuple_out.iv_pos = tuple_in.iv_pos;
            stream_tuple_out.write(tuple_out);
            for (int g = 0; g < K_FUN; g++){
#pragma HLS unroll
                filter[g] = stream_filter_in.read();
            }

            ap_uint<64> hash_out;
            xf::database::details::hashlookup3_core<V_ID_W>(tuple_in.indexing_v, hash_out);
            volatile unsigned int start_off = 0;
            volatile unsigned int end_off;
            unsigned int addr_inrow;
            ap_uint<DDR_W> ram_row;
            unsigned long addr_row;
            ap_uint<64> addr_counter;
            hash_out = hash_out.range(hash1_w - 1, 0);
            
            if (hash_out != 0){ 
                addr_counter = hash_out - 1;
                addr_counter <<= hash2_w;
                addr_counter += (1UL << hash2_w) - 1;

                /* Compute address of row storing the counter */
                addr_row = hTables[tuple_in.tb_index].start_offset 
                    + (addr_counter >> (DDR_BIT - C_W));

                /* Compute address of data inside the row */
                addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1,0);

                /* Read the data */
                ram_row = m_axi[addr_row];
                start_off = ram_row.range(((addr_inrow + 1) << C_W) - 1, addr_inrow << C_W);
            }


            addr_counter = hash_out;
            addr_counter <<= hash2_w;
            addr_counter += (1UL << hash2_w) - 1;

            /* Compute address of row storing the counter */
            addr_row = hTables[tuple_in.tb_index].start_offset 
                + (addr_counter >> (DDR_BIT - C_W));

            /* Compute address of data inside the row */
            addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1,0);

            /* Read the data */
            ram_row = m_axi[addr_row];
            end_off = ram_row.range(((addr_inrow + 1) << C_W) - 1, addr_inrow << C_W);

#if DEBUG_STATS
            debug::readmin_reads += 2;
#endif
           
            unsigned int rowstart = hTables[tuple_in.tb_index].start_edges + 
                (start_off >> (DDR_BIT - E_W));
            unsigned int rowend = hTables[tuple_in.tb_index].start_edges + 
                (end_off >> (DDR_BIT - E_W));
            
READMIN_EDGES_LOOP:
            for (int g = rowstart; g <= rowend; g++){
                row_t row = m_axi[g];

                for (int i = 0; i < EDGE_ROW; i++){
#pragma HLS unroll
                    edge = row.range(((i + 1) << E_W) - 1, i << E_W);
                    vertexCheck = edge.range(V_ID_W * 2 - 1, V_ID_W);
                    vertex = edge.range(V_ID_W - 1, 0);

                    if (tuple_in.indexing_v == vertexCheck){

                        xf::database::details::hashlookup3_core<V_ID_W>(vertex,
                                hash_out);
                        bool test = true;
                        for (int g = 0; g < K_FUN; g++){
#pragma HLS unroll
                            ap_uint<BLOOM_LOG> idx = hash_out.range((64 / K_FUN) * (g + 1) - 1,
                                    (64 / K_FUN) * (g + 1) - BLOOM_LOG);
/* std::cout << "Testing " << idx << " in " << std::hex << */
/* filter[g] << std::dec << std::endl; */
                            test = test && (filter[g][idx] == 1);
                        }

                        if (test){
                            stream_set_out.write(vertex);
                            stream_set_end_out.write(false);
#if DEBUG_STATS
                            debug::readmin_vstream++;
                        }
                        else {
                            debug::bloom_filter++;
#endif
                        }
                    }
                }
#if DEBUG_STATS
                debug::readmin_reads++;
#endif
            }
            stream_set_end_out.write(true);
#if DEBUG_STATS
            debug::readmin_n_sets++;
#endif
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

void mwj_homomorphism(
        hls::stream< ap_uint<V_ID_W> >  &stream_set_in,
        hls::stream<bool>               &stream_set_end_in,
        hls::stream<tuplebuild_tuple_t> &stream_tuple_in,
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_in,
        hls::stream<bool>               &stream_sol_end_in,
        
        hls::stream< ap_uint<V_ID_W> >  &stream_set_out,
        hls::stream<bool>               &stream_set_end_out,
        hls::stream<tuplebuild_tuple_t> &stream_tuple_out,
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_out,
        hls::stream<bool>               &stream_sol_end_out)
{
    ap_uint<8> curQV {0};
    ap_uint<2> last_set;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool last_sol;

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

    stream_tuple_out.write(stream_tuple_in.read());
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
#if DEBUG_STATS
        else {
            debug::homomo_trashed++;
        }
#endif
        last_set = stream_set_end_in.read();
    }
    stream_set_end_out.write(last_set);
}

template <size_t MAX_BATCH_SIZE>
void mwj_batchbuild(
        hls::stream< ap_uint<V_ID_W> >  &stream_set_in,
        hls::stream<bool>               &stream_set_end_in,
        hls::stream<tuplebuild_tuple_t> &stream_tuple_in,
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_in,
        hls::stream<bool>               &stream_sol_end_in,
       
        hls::stream<bool>               &stream_req,
        hls::stream< ap_uint<V_ID_W> >  &stream_batch_out,
        hls::stream<bool>               &stream_batch_end_out,
        hls::stream<tuplebuild_tuple_t> &stream_tuple_out,
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_out,
        hls::stream<bool>               &stream_sol_end_out)
{
    ap_uint<8> curQV {0};
    ap_uint<V_ID_W> curEmb[MAX_QV];
    unsigned int batch_counter {0};
    bool last;
    
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

    tuplebuild_tuple_t tuple_in = stream_tuple_in.read();
    stream_tuple_out.write(tuple_in);
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

BATCHBUILD_ADD_SOLUTION_LOOP:
            for (int g = 0; g < curQV; g++){
                stream_sol_out.write(curEmb[g]);
                stream_sol_end_out.write(false);
            }
            stream_sol_end_out.write(true);
            stream_tuple_out.write(tuple_in);
        }
    }
    stream_batch_end_out.write(true);
}

template <size_t BATCH_SIZE_LOG>
void mwj_tuplebuild(
        const unsigned char             hash1_w,
        const unsigned char             hash2_w,
        QueryVertex                     *qVertices,
        hls::stream< ap_uint<V_ID_W> >  &stream_set_in,
        hls::stream<bool>               &stream_set_end_in,
        hls::stream<tuplebuild_tuple_t> &stream_tuple_in,
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_in,
        hls::stream<bool>               &stream_sol_end_in,
        hls::stream<bool>               &stream_stop,
       
        hls::stream<intersect_tuple_t>  &stream_tuple_out,
        hls::stream<bool>               &stream_tuple_end_out,
        hls::stream< ap_uint<V_ID_W> >  &stream_sol_out,
        hls::stream<bool>               &stream_sol_end_out)
{
    const unsigned int hash1_mask = (1UL << hash1_w) - 1;
    const unsigned int hash2_mask = (1UL << hash2_w) - 1;
    ap_uint<V_ID_W> buffer[(1UL << BATCH_SIZE_LOG)];
    ap_uint<BATCH_SIZE_LOG> buffer_size {0};
    ap_uint<BATCH_SIZE_LOG> buffer_p;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<V_ID_W> vToVerify;
    hls::stream<ap_uint<V_ID_W>, 4> hash_in0, hash_in1;
    hls::stream<ap_uint<64>, 4> hash_out0, hash_out1;
    unsigned long addr_counter;
    intersect_tuple_t tuple_out;
    unsigned char curQV {0};
    bool stop, last;
   
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
            tuplebuild_tuple_t tuple_in = stream_tuple_in.read();
            unsigned char cycles = qVertices[curQV].numTablesIndexed;

            if (cycles > 0){
                uint8_t tableIndex = qVertices[curQV].tables_indexed[0];
                uint8_t ivPos = qVertices[curQV].vertex_indexing[0];
                bool bit_min = (tuple_in.tb_index != tableIndex || 
                        tuple_in.iv_pos != ivPos);

                last = stream_set_end_in.read();

TUPLEBUILD_MAIN_LOOP_FIRST_IT:
                while(!last){
#pragma HLS pipeline II=1 
                    vToVerify = stream_set_in.read();
                     
                    hash_in0.write(vToVerify);
                    hash_in1.write(curEmb[ivPos]);
                    xf::database::hashLookup3<V_ID_W>(hash_in0, hash_out0);
                    xf::database::hashLookup3<V_ID_W>(hash_in1, hash_out1);
                    
                    addr_counter = hash_out1.read() & hash1_mask;
                    addr_counter <<= hash2_w;
                    addr_counter += hash_out0.read() & hash2_mask;

                    tuple_out.indexed_v = vToVerify;
                    tuple_out.indexing_v = curEmb[ivPos];
                    tuple_out.addr_counter = addr_counter;
                    tuple_out.tb_index = tableIndex;
                    tuple_out.pos = buffer_size;
                    tuple_out.bit_last_edge = (qVertices[curQV].numTablesIndexed == 1);
                    tuple_out.bit_min_set = bit_min;

                    stream_tuple_out.write(tuple_out);
                    stream_tuple_end_out.write(false);
/* std::cout << "( " */
/* << tuple_out.indexing_v */
/* << ", " << tuple_out.indexed_v */
/* << ", " << tuple_out.indexing_h */
/* << ", " << tuple_out.indexed_h */
/* << ", " << (int)tuple_out.tb_index */
/* << ", " << tuple_out.pos */
/* << ", " << tuple_out.bit_last_edge */
/* << ", " << tuple_out.bit_min_set << ")" << std::endl; */

                    buffer[buffer_size++] = vToVerify;
                    last = stream_set_end_in.read();
                }
            }

TUPLEBUILD_EDGE_LOOP_AFTER_IT:
            for (int g = 0; g < cycles - 1; g++){
                uint8_t tableIndex = qVertices[curQV].tables_indexed[g + 1];
                uint8_t ivPos = qVertices[curQV].vertex_indexing[g + 1];
                bool bit_last = (g == cycles - 2);
                bool bit_min = (tuple_in.tb_index != tableIndex || 
                        tuple_in.iv_pos != ivPos);

TUPLEBUILD_MAIN_LOOP:
                for (int buffer_p = 0; buffer_p < buffer_size; buffer_p++){
#pragma HLS pipeline II=1 
                    vToVerify = buffer[buffer_p];
                    
                    hash_in0.write(vToVerify);
                    hash_in1.write(curEmb[ivPos]);
                    xf::database::hashLookup3<V_ID_W>(hash_in0, hash_out0);
                    xf::database::hashLookup3<V_ID_W>(hash_in1, hash_out1);
                    
                    addr_counter = hash_out1.read() & hash1_mask;
                    addr_counter <<= hash2_w;
                    addr_counter += hash_out0.read() & hash2_mask;
                    
                    tuple_out.indexed_v = vToVerify;
                    tuple_out.indexing_v = curEmb[ivPos];
                    tuple_out.addr_counter = addr_counter;
                    tuple_out.tb_index = tableIndex;
                    tuple_out.pos = buffer_p;
                    tuple_out.bit_last_edge = bit_last;
                    tuple_out.bit_min_set = bit_min;
                    
                    stream_tuple_out.write(tuple_out);
                    stream_tuple_end_out.write(false);
                
/* std::cout << "( " */
/* << tuple_out.indexing_v */
/* << ", " << tuple_out.indexed_v */
/* << ", " << tuple_out.indexing_h */
/* << ", " << tuple_out.indexed_h */
/* << ", " << (int)tuple_out.tb_index */
/* << ", " << tuple_out.pos */
/* << ", " << tuple_out.bit_last_edge */
/* << ", " << tuple_out.bit_min_set << ")" << std::endl; */
                }
            }
            stream_tuple_end_out.write(true);
        }

        if (stream_stop.read_nb(stop))
            break;
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
        AdjHT                          *hTables,
        htb_cache_t                    &htb_buf,
        hls::stream<intersect_tuple_t> &stream_tuple_in,
        hls::stream<bool>              &stream_tuple_end_in,
        hls::stream<bool>              &stream_stop,

        hls::stream<split_tuple_t>     &stream_tuple_out,
        hls::stream<bool>              &stream_tuple_end_out)
{
    ap_uint<V_ID_W> indexing_v;
    intersect_tuple_t tuple_in;
    split_tuple_t tuple_out;
    unsigned char tableIndex;
    ap_uint<DDR_W> ram_row;
    unsigned long addr_row;
    unsigned long addr_inrow;
    ap_uint<64> addr_counter;
    bool stop, last;

INTERSECT_TASK_LOOP:
    while (1) {
#pragma HLS pipeline II=2
        if (stream_tuple_end_in.read_nb(last)){

            if (!last){
                tuple_in = stream_tuple_in.read();
                ap_uint<(1UL << C_W)> start_off = 0;
                ap_uint<(1UL << C_W)> end_off = 1;
                tableIndex = tuple_in.tb_index;
                addr_counter = tuple_in.addr_counter;

                if (tuple_in.bit_min_set){
                    
                    /* Compute address of row storing the counter */
                    addr_row = hTables[tableIndex].start_offset + 
                        (addr_counter >> (DDR_BIT - C_W));

                    /* Compute address of data inside the row */
                    addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);
                    ram_row = htb_buf.get(addr_row, 1);
                    end_off = ram_row.range(((addr_inrow + 1) << C_W) - 1, 
                            addr_inrow << C_W);

                    if (addr_counter != 0){
                        addr_counter--;

                        //Compute address of row storing the counter
                        addr_row = hTables[tableIndex].start_offset +
                            (addr_counter >> (DDR_BIT - C_W));

                        //Compute address of data inside the row
                        addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);
                        ram_row = htb_buf.get(addr_row, 1);
                        start_off = ram_row.range(((addr_inrow + 1) << C_W) - 1,
                                addr_inrow << C_W);
                    }
#if DEBUG_STATS
                    debug::intersect_reads += 2;
#endif /* DEBUG_STATS */
                }

                tuple_out.indexed_v     = tuple_in.indexed_v;
                tuple_out.indexing_v    = tuple_in.indexing_v;
                tuple_out.tb_index      = tuple_in.tb_index;
                tuple_out.pos           = tuple_in.pos;
                tuple_out.start_off     = start_off;
                tuple_out.n_edges       = end_off - start_off;
                tuple_out.bit_last_edge = tuple_in.bit_last_edge;
                tuple_out.bit_min_set   = tuple_in.bit_min_set;
                tuple_out.bit_no_edge   = (start_off == end_off);

#if DEBUG_STATS
                debug::intersect_filter += (start_off == end_off)? 1 : 0;
#endif /* DEBUG_STATS */
/* std::cout << "( " */
/* << tuple_out.indexing_v */
/* << ", " << tuple_out.indexed_v */
/* << ", " << (int)tuple_out.tb_index */
/* << ", " << tuple_out.pos */
/* << ", " << tuple_out.start_off */
/* << ", " << tuple_out.n_edges */
/* << ", " << tuple_out.bit_last_edge */
/* << ", " << tuple_out.bit_min_set */
/* << ", " << tuple_out.bit_no_edge << ")" << std::endl; */
                /* bits[pos] = bits[pos] & (start_off < end_off); */
                stream_tuple_out.write(tuple_out);
            }
            stream_tuple_end_out.write(last);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

void mwj_split(
        hls::stream<split_tuple_t>  &stream_tuple_in,
        hls::stream<bool>           &stream_tuple_end_in,

        hls::stream<verify_tuple_t> &stream_tuple_out,
        hls::stream<bool>           &stream_tuple_end_out)
{
    constexpr size_t EDGE_BLOCK = (CACHE_WORDS_PER_LINE + DDR_BIT - E_W);
    bool last;
    split_tuple_t tuple_in;
    verify_tuple_t tuple_out;
    unsigned short n_edges;

    last = stream_tuple_end_in.read();
    if (!last){
        tuple_in = stream_tuple_in.read();
        n_edges = (tuple_in.n_edges == 0) ? 1 : tuple_in.n_edges;
        
        tuple_out.indexing_v    = tuple_in.indexing_v;
        tuple_out.indexed_v     = tuple_in.indexed_v;
        tuple_out.tb_index      = tuple_in.tb_index;
        tuple_out.pos           = tuple_in.pos;
        tuple_out.bit_last_edge = tuple_in.bit_last_edge;
        tuple_out.bit_min_set   = tuple_in.bit_min_set;
        tuple_out.bit_no_edge   = tuple_in.bit_no_edge;

        int first_block = tuple_in.start_off >> EDGE_BLOCK ;
        int end_block = (tuple_in.start_off + n_edges - 1) >> EDGE_BLOCK;
SPLIT_MAIN_LOOP:
        for (int s = 0; first_block <= end_block; first_block++, s += (1UL << EDGE_BLOCK)){
#pragma HLS pipeline II=1
            tuple_out.address = tuple_in.start_off + s;
            tuple_out.bit_last_address = (first_block == end_block);
            stream_tuple_out.write(tuple_out);
            stream_tuple_end_out.write(false);

/* std::cout << "( " */
/* << tuple_out.indexing_v */
/* << ", " << tuple_out.indexed_v */
/* << ", " << (int)tuple_out.tb_index */
/* << ", " << tuple_out.pos */
/* << ", " << tuple_out.address */
/* << ", " << tuple_out.bit_last_address */
/* << ", " << tuple_out.bit_last_edge */
/* << ", " << tuple_out.bit_no_edge */
/* << ", " << tuple_out.bit_min_set << ")" << std::endl; */
        }
    } else {
        stream_tuple_end_out.write(true);
    }

}

template <size_t BATCH_SIZE_LOG>
void mwj_verify(
        AdjHT *hTables,
        htb_cache_t                  &htb_buf,
        hls::stream<verify_tuple_t>  &stream_tuple_in,
        hls::stream<bool>            &stream_tuple_end_in,
        hls::stream<bool>            &stream_stop,

        hls::stream<compact_tuple_t> &stream_tuple_out,
        hls::stream<bool>            &stream_tuple_end_out)
{
    constexpr size_t EDGE_PER_WORD = (DDR_BIT - E_W);
    
    verify_tuple_t tuple_in;
    compact_tuple_t tuple_out;
    ap_uint<V_ID_W> candidate_v;
    ap_uint<V_ID_W> indexing_v;
    unsigned char tableIndex;
    bool stop, last;
    ap_uint<128> edge_block[(1UL << CACHE_WORDS_PER_LINE)];
#pragma HLS array_partition variable=edge_block type=complete

VERIFY_TASK_LOOP:
    while(1){
#pragma HLS pipeline II=1
        if (stream_tuple_end_in.read_nb(last)){

            if(!last){
                tuple_in = stream_tuple_in.read();
                candidate_v = tuple_in.indexed_v;
                tuple_out.bit_equal = !tuple_in.bit_min_set && !tuple_in.bit_no_edge;

                if (tuple_in.bit_min_set && !tuple_in.bit_no_edge){
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

                tuple_out.bit_last_edge     = tuple_in.bit_last_edge;
                tuple_out.bit_last_address  = tuple_in.bit_last_address;
                tuple_out.indexed_v         = candidate_v;
                tuple_out.pos               = tuple_in.pos;
                stream_tuple_out.write(tuple_out);
            }
            stream_tuple_end_out.write(last);
        }

        if (stream_stop.read_nb(stop)){
            break;
        }
    }
}

void mwj_compact(
        hls::stream<compact_tuple_t>  &stream_tuple_in,
        hls::stream<bool>             &stream_tuple_end_in,
        
        hls::stream<assembly_tuple_t> &stream_tuple_out,
        hls::stream<bool>             &stream_tuple_end_out)
{
#pragma HLS pipeline II=1
    static bool checked = false;
    
    compact_tuple_t tuple_in;
    assembly_tuple_t tuple_out;
    bool last;

    last = stream_tuple_end_in.read();
    if (!last) {
        tuple_in = stream_tuple_in.read();
        checked |= tuple_in.bit_equal;
        if (tuple_in.bit_last_address){
            tuple_out.bit_checked = checked;
            tuple_out.indexed_v = tuple_in.indexed_v;
            tuple_out.pos = tuple_in.pos;
            tuple_out.bit_last_edge = tuple_in.bit_last_edge;
            stream_tuple_out.write(tuple_out); 
            stream_tuple_end_out.write(false);
#if DEBUG_STATS
            if (!checked)
                debug::verify_filter++;
#endif /* DEBUG_STATS */
            checked = false;
        }
    } else {
        stream_tuple_end_out.write(true);
        checked = false; 
    }
}

template<size_t MAX_BATCH_SIZE>
void mwj_filter(
        hls::stream<assembly_tuple_t>   &stream_tuple_in,
        hls::stream<bool>               &stream_tuple_end_in,

        hls::stream< ap_uint<V_ID_W> >  &stream_set_out,
        hls::stream<bool>               &stream_set_end_out)
{
#pragma HLS pipeline II=1
    static ap_uint<MAX_BATCH_SIZE> bits = ~0;
    assembly_tuple_t tuple_in;
    bool last;
    
    last = stream_tuple_end_in.read();
    if (!last){
        tuple_in = stream_tuple_in.read();
        unsigned short p = tuple_in.pos;
        bits[p] = bits[p] && tuple_in.bit_checked;
        if (tuple_in.bit_last_edge && bits.test(p)){
            stream_set_out.write(tuple_in.indexed_v);
            stream_set_end_out.write(false);
        }
    } else {
        bits = ~0;
        stream_set_end_out.write(true);
    }
}

void mwj_assembly(
        unsigned short nQueryVer,
        hls::stream< ap_uint<V_ID_W> >    &stream_inter_in,
        hls::stream<bool>               &stream_end_inter_in,
        hls::stream< ap_uint<V_ID_W> >    &stream_embed_in,
        hls::stream<bool>               &stream_end_embed_in,
        hls::stream< ap_uint<V_ID_W> >    &stream_batch,
        hls::stream<bool>               &stream_batch_end,
       
        hls::stream<bool>               streams_stop[STOP_S],
        hls::stream<bool>               &stream_sol0,
        hls::stream<bool>               &stream_req, 
        hls::stream< ap_uint<V_ID_W> >    &stream_partial_out,
#if COUNT_ONLY
        long unsigned int               &result
#else  
        hls::stream<T_NODE>             &result
#endif
        )
{
    ap_uint<V_ID_W> curQV;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool last_sol, last_set, stop;
    bool last_start, token_new_start;
    T_NODE node;

#if COUNT_ONLY
    unsigned long int counter {0};
#endif

    last_start = stream_batch_end.read();
    last_start = stream_batch_end.read();
    token_new_start = false;
    stream_partial_out.write((1UL << (V_ID_W - 1)) + 1);
    stream_partial_out.write(stream_batch.read());
    stream_req.write(true);

ASSEMBLY_TASK_LOOP:
    while(1) {
        if (stream_end_embed_in.read_nb(last_sol)){

            curQV = 0;
ASSEMBLY_COPYING_EMBEDDING_LOOP:
            while(!last_sol){
                curEmb[curQV] = stream_embed_in.read();
                curQV++;
                last_sol = stream_end_embed_in.read();
            }

            last_set = stream_end_inter_in.read();

            //Last bit used to signal new solution
            if (!last_set && (curQV != nQueryVer - 1)){
                stream_partial_out.write((curQV + 1) | (1UL << (V_ID_W - 1)));
ASSEMBLY_RADIX_LOOP:
                for (int g = 0; g < curQV; g++){
                    stream_partial_out.write(curEmb[g]);
                }
            }

ASSEMBLY_SET_LOOP:
            while(!last_set){
                ap_uint<V_ID_W> vToVerify = stream_inter_in.read();
                /* Write in the correct stream */
                if (curQV == nQueryVer - 1){
#if COUNT_ONLY
                    if (!last_start && token_new_start){
                        last_start = stream_batch_end.read();
                        token_new_start = false;
                        stream_partial_out.write((1UL << (V_ID_W - 1)) + 1);
                        stream_partial_out.write(stream_batch.read());
                        stream_req.write(true);
                    }
                    counter++;
#else
ASSEMBLY_WRITE_FINAL_LOOP:
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
                    token_new_start = true;
                    stream_partial_out.write(vToVerify);
                    stream_req.write(true);
                }
                last_set = stream_end_inter_in.read();
            }

            // Last batch of a set 
            stream_req.write(false);

        }

        if (stream_sol0.read_nb(stop)){

#if DEBUG_STATS
                debug::miss_indexing++;
#endif

            // Test if there are some node from start batch 
            if (!last_start){
                last_start = stream_batch_end.read();
                token_new_start = true;
                stream_partial_out.write((1UL << (V_ID_W - 1)) + 1);
                stream_partial_out.write(stream_batch.read());
                stream_req.write(true);
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
        hls::stream<bool> &dynfifo_overflow,
        hls::stream<bool> &stream_sol0)
{
    static unsigned long sol {0};
    bool ovf {false};
    bool req = stream_req.read();
    
    if (req) {
        sol++;
    } else {
        sol--;
    }

    dynfifo_overflow.read_nb(ovf);

    if (sol == 0 || ovf == true){
        sol = 0;
        stream_sol0.write(true);
    }
}

void mwj_batch(
        const unsigned char hash1_w,
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<DDR_W> *htb_buf,
        
        hls::stream<bool> &stream_batch_end,
        hls::stream< ap_uint<V_ID_W> > &stream_batch)
{
    ap_uint<8> tableIndex {0};
    ap_uint<32> minSize = (1UL << 32) - 1;
    ap_uint<32> minStart, minOff;
    ap_uint<8 + 1 + V_ID_W> minData;
    ap_uint<V_ID_W*2> edge;
    ap_uint<V_ID_W> vertex, vertexCheck;
    ap_uint<V_ID_W> set[MAX_CL];
    ap_uint<64> hash_buff = 0;
    unsigned char set_counter = 0;
    bool flag_buff = false;
    bool flag_new = true;

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
    for (unsigned int g = rowstart; g <= rowend; g++){
        row_t row = htb_buf[g];
        for (unsigned int i = 0; i < EDGE_ROW; i++, cnt++){
            if (cnt < window_right){
                edge = row.range((1UL << E_W) - 1, 0);
                vertex = edge.range(V_ID_W * 2 - 1, V_ID_W);
                ap_uint<64> hash_out;
                xf::database::details::hashlookup3_core<V_ID_W>(
                        vertex,
                        hash_out);
                hash_out = hash_out.range(hash1_w - 1, 0);

                if (flag_buff && hash_buff == hash_out){
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
#if DEBUG_STATS
                    debug::start_set++;
#endif
                    set[set_counter++] = vertex;
                    stream_batch_end.write(false);
                    stream_batch.write(vertex);
                }

                hash_buff = hash_out;
                flag_buff = true;
            }
#if DEBUG_STATS
            debug::batch_reads++;
#endif
            row >>= (1UL << E_W);
        }
    }
    stream_batch_end.write(true);
}

template <typename T_BLOOM,
          size_t BLOOM_LOG,
          size_t K_FUN_LOG>
void readmincache_wrapper(
    AdjHT *hTables,
    htb_cache_t &htb_buf,
    hls::stream<readmin_tuple_t> &stream_tuple_in,
    hls::stream<T_BLOOM> &stream_filter_in,
    hls::stream<bool> &stream_stop,

    hls::stream< ap_uint<V_ID_W> > &stream_set_out,
    hls::stream<bool> &stream_set_end_out,
    hls::stream<tuplebuild_tuple_t> &stream_tuple_out)
{
    htb_buf.init(2);
    mwj_readmin<T_BLOOM,
                BLOOM_LOG,
                K_FUN_LOG>
                (
        hTables,
        htb_buf,
        stream_tuple_in,
        stream_filter_in,
        stream_stop,
        stream_set_out,
        stream_set_end_out,
        stream_tuple_out);
}

template <size_t BATCH_SIZE_LOG>
void intersectcache_wrapper(
    AdjHT *hTables,
    htb_cache_t &htb_buf,
    hls::stream<intersect_tuple_t> &stream_tuple_in,
    hls::stream<bool> &stream_tuple_end_in,
    hls::stream<bool> &stream_stop,

    hls::stream<split_tuple_t> &stream_tuple_out,
    hls::stream<bool> &stream_tuple_end_out)
{
    htb_buf.init(1);
    mwj_intersect<BATCH_SIZE_LOG>(
        hTables,
        htb_buf,
        stream_tuple_in,
        stream_tuple_end_in,
        stream_stop,
        stream_tuple_out,
        stream_tuple_end_out);
}

template <size_t BATCH_SIZE_LOG>
void verifycache_wrapper(
        AdjHT *hTables,
        htb_cache_t                  &htb_buf,
        hls::stream<verify_tuple_t>  &stream_tuple_in,
        hls::stream<bool>            &stream_tuple_end_in,
        hls::stream<bool>            &stream_stop,

        hls::stream<compact_tuple_t> &stream_tuple_out,
        hls::stream<bool>            &stream_tuple_end_out)
{
    htb_buf.init(0);
    mwj_verify<BATCH_SIZE_LOG>(
            hTables,
            htb_buf,
            stream_tuple_in,
            stream_tuple_end_in,          
            stream_stop,                  
            stream_tuple_out,             
            stream_tuple_end_out);
    htb_buf.stop();
}

template <typename T_BLOOM,
        size_t BLOOM_LOG, 
        size_t K_FUN_LOG>
void multiwayJoin(
        ap_uint<DDR_W>                  *htb_buf0,
        ap_uint<DDR_W>                  *htb_buf1,
        T_BLOOM                         *bloom_p,
        row_t                           *res_buf,
        AdjHT                           *hTables0,
        AdjHT                           *hTables1,
        QueryVertex                     *qVertices0,
        const unsigned short            nQueryVer,
        const unsigned char             hash1_w,
        const unsigned char             hash2_w,
        hls::stream< ap_uint<V_ID_W> >  &stream_batch,
        hls::stream<bool>               &stream_batch_end,
        unsigned long                   &dynfifo_diagnostic,
        
#if COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif
        )
{
#pragma HLS STABLE variable=htb_buf0
#pragma HLS STABLE variable=htb_buf1
#pragma HLS STABLE variable=bloom_p
#pragma HLS STABLE variable=hTables0
#pragma HLS STABLE variable=hTables1
#pragma HLS STABLE variable=qVertices0
#pragma HLS STABLE variable=nQueryVer
#pragma HLS STABLE variable=hash1_w
#pragma HLS STABLE variable=hash2_w

#pragma HLS DATAFLOW

    /* Findmin data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> p_stream_sol
        ("Propose - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> p_stream_sol_end
        ("Propose - partial solution end flag");
    hls_thread_local hls::stream<T_BLOOM, 3> p_stream_filter
        ("Propose - filter");
    hls_thread_local hls::stream<readmin_tuple_t, S_D> p_stream_tuple
        ("Propose - tuples");

    /* Readmin data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> r_stream_sol
        ("Readmin - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> r_stream_sol_end
        ("Readmin - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> r_stream_set
        ("Readmin - set nodes");
    hls_thread_local hls::stream<bool, S_D> r_stream_set_end
        ("Readmin - set nodes end flag");
    hls_thread_local hls::stream<tuplebuild_tuple_t, S_D> r_stream_tuple
        ("Readmin - tuples");
    
    /* Homomorphism data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> h_stream_sol
        ("Homomorphism - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> h_stream_sol_end
        ("Homomorphism - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> h_stream_set
        ("Homomorphism - set nodes");
    hls_thread_local hls::stream<bool, S_D> h_stream_set_end
        ("Homomorphism - set nodes end flag");
    hls_thread_local hls::stream<tuplebuild_tuple_t, S_D> h_stream_tuple
        ("Homomorphism - tuples");
    
    /* Batchbuild data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> b_stream_sol
        ("Batchbuild - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> b_stream_sol_end
        ("Batchbuild - partial solution end flag");
    hls_thread_local hls::stream<ap_uint<V_ID_W>, S_D> b_stream_set
        ("Batchbuild - set nodes");
    hls_thread_local hls::stream<bool, S_D> b_stream_set_end
        ("Batchbuild - set nodes end flag");
    hls_thread_local hls::stream<tuplebuild_tuple_t, S_D> b_stream_tuple
        ("Batchbuild - tuples");
    
    /* Tuplebuild data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> t_stream_sol
        ("Tuplebuild - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> t_stream_sol_end
        ("Tuplebuild - partial solution end flag");
    hls_thread_local hls::stream<intersect_tuple_t, S_D> t_stream_tuple
        ("Tuplebuild - tuples");
    hls_thread_local hls::stream<bool, S_D> t_stream_tuple_end
        ("Tuplebuild - tuples end flag");
    
    /* Intersect data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> i_stream_sol
        ("Intersect - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> i_stream_sol_end
        ("Intersect - partial solution end flag");
    hls_thread_local hls::stream<split_tuple_t, S_D> i_stream_tuple
        ("Intersect - tuples");
    hls_thread_local hls::stream<bool, S_D> i_stream_tuple_end
        ("Intersect - tuples end flag");
    
    /* Split data out */    
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> s_stream_sol
        ("Split - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> s_stream_sol_end
        ("Split - partial solution end flag");
    hls_thread_local hls::stream<verify_tuple_t, S_D> s_stream_tuple
        ("Split - tuples");
    hls_thread_local hls::stream<bool, S_D> s_stream_tuple_end
        ("Split - tuples end flag");
    
    /* Verify data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> v_stream_sol
        ("Verify - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> v_stream_sol_end
        ("Verify - partial solution end flag");
    hls_thread_local hls::stream<compact_tuple_t, S_D> v_stream_tuple
        ("Verify - tuples");
    hls_thread_local hls::stream<bool, S_D> v_stream_tuple_end
        ("Verify - tuples end flag");
    
    /* Compact data out */
    hls_thread_local hls::stream<ap_uint<V_ID_W>, MAX_QV> c_stream_sol
        ("Compact - partial solution");
    hls_thread_local hls::stream<bool, MAX_QV> c_stream_sol_end
        ("Compact - partial solution end flag");
    hls_thread_local hls::stream<assembly_tuple_t, S_D> c_stream_tuple
        ("Compact - tuples");
    hls_thread_local hls::stream<bool, S_D> c_stream_tuple_end
        ("Compact - tuples end flag");

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
    hls_thread_local hls::stream<bool, 1> dyn_stream_ovf
        ("Dynamic fifo - overflow");
 
    /* Stop signals */
    hls_thread_local hls::stream<bool, 1> streams_stop[STOP_S];
    hls_thread_local hls::stream<bool, S_D> merge_out;
    hls_thread_local hls::stream<bool, S_D> merge_in[MERGE_IN_STREAMS];
    hls_thread_local hls::stream<bool, 1> stream_sol0;
#pragma HLS array_partition variable=merge_in type=complete    

#ifndef __SYNTHESIS__
    for (int g = 0; g < STOP_S; g++) {
        char stream_name[10];
        sprintf(stream_name, "stop_%d", g);
        streams_stop[g].set_name(stream_name);
    } 
#endif /* __SYNTHESIS__ */

    htb_cache_t htb_cache(htb_buf0);
    bloom_cache_t bloom_cache(bloom_p);
    
    dynfifo_init<
        ap_uint<V_ID_W>,        /* fifo data type */
        row_t,                  /* fifo data type */
        DYN_FIFO_DEPTH,         /* in/out stream size */
        DYN_FIFO_BURST * 2,     /* load/store stream size */
        DDR_WORD,               /* bitwidth ddr word */
        DYN_FIFO_BURST,         /* burst transaction size */
        RESULTS_SPACE>          /* memory words available */
            (res_buf,   
             dynfifo_diagnostic,
             a_stream_sol,
             dyn_stream_sol,
             streams_stop[STOP_S - 2],
             streams_stop[STOP_S - 1],
             dyn_stream_ovf);

    hls_thread_local hls::task mwj_homomorphism_t(
            mwj_homomorphism,
            r_stream_set,
            r_stream_set_end,
            r_stream_tuple,
            r_stream_sol,
            r_stream_sol_end,
            h_stream_set,
            h_stream_set_end,
            h_stream_tuple,
            h_stream_sol,
            h_stream_sol_end);

    hls_thread_local hls::task mwj_batchbuild_t(
            mwj_batchbuild<(1UL << PROPOSE_BATCH_LOG)>,
            h_stream_set,
            h_stream_set_end,
            h_stream_tuple,
            h_stream_sol,
            h_stream_sol_end,
            merge_in[0],
            b_stream_set,
            b_stream_set_end,
            b_stream_tuple,
            b_stream_sol,
            b_stream_sol_end);
    
    hls_thread_local hls::task mwj_bypass_sol_t5(
            mwj_bypass_sol,
            p_stream_sol,
            p_stream_sol_end,
            r_stream_sol,
            r_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_t0(
            mwj_bypass_sol,
            t_stream_sol,
            t_stream_sol_end,
            i_stream_sol,
            i_stream_sol_end);

    hls_thread_local hls::task mwj_bypass_sol_t1(
            mwj_bypass_sol,
            i_stream_sol,
            i_stream_sol_end,
            s_stream_sol,
            s_stream_sol_end);
    
    hls_thread_local hls::task mwj_bypass_sol_t2(
            mwj_bypass_sol,
            s_stream_sol,
            s_stream_sol_end,
            v_stream_sol,
            v_stream_sol_end);
    
    hls_thread_local hls::task mwj_bypass_sol_t3(
            mwj_bypass_sol,
            v_stream_sol,
            v_stream_sol_end,
            c_stream_sol,
            c_stream_sol_end);
    
    hls_thread_local hls::task mwj_bypass_sol_t4(
            mwj_bypass_sol,
            c_stream_sol,
            c_stream_sol_end,
            f_stream_sol,
            f_stream_sol_end);
    
    hls_thread_local hls::task mwj_split_t(
            mwj_split,
            i_stream_tuple,
            i_stream_tuple_end,
            s_stream_tuple,
            s_stream_tuple_end);

    hls_thread_local hls::task mwj_compact_t(
            mwj_compact,
            v_stream_tuple,
            v_stream_tuple_end,
            c_stream_tuple,
            c_stream_tuple_end);

    hls_thread_local hls::task mwj_filter_t(
            mwj_filter<(1UL << PROPOSE_BATCH_LOG)>,
            c_stream_tuple,
            c_stream_tuple_end,
            f_stream_set,
            f_stream_set_end);

#ifdef __SYNTHESIS__

    cache_wrapper(
        mwj_findmin<
            T_BLOOM,
            BLOOM_LOG,
            K_FUN_LOG>,
        hash1_w,
        qVertices0,
        bloom_cache,
        dyn_stream_sol,
        streams_stop[0],
        p_stream_tuple,
        p_stream_filter,
        p_stream_sol,
        p_stream_sol_end);

    mwj_readmin<
        T_BLOOM,
        BLOOM_LOG,
        K_FUN_LOG>(
        hash1_w,
        hash2_w,
        hTables0,
        htb_buf1,
        p_stream_tuple,
        p_stream_filter,
        streams_stop[1],
        r_stream_set,
        r_stream_set_end,
        r_stream_tuple);

    mwj_tuplebuild<PROPOSE_BATCH_LOG>(
        hash1_w,
        hash2_w,
        qVertices0,
        b_stream_set,
        b_stream_set_end,
        b_stream_tuple,
        b_stream_sol,
        b_stream_sol_end,
        streams_stop[2],
        t_stream_tuple,
        t_stream_tuple_end,
        t_stream_sol,
        t_stream_sol_end);

    intersectcache_wrapper<PROPOSE_BATCH_LOG>(
        hTables1,
        htb_cache,
        t_stream_tuple,
        t_stream_tuple_end,
        streams_stop[3],
        i_stream_tuple,
        i_stream_tuple_end);

    verifycache_wrapper<PROPOSE_BATCH_LOG>(
            hTables1,
            htb_cache,
            s_stream_tuple,
            s_stream_tuple_end,
            streams_stop[4],
            v_stream_tuple,
            v_stream_tuple_end);
    
    mwj_assembly(
            nQueryVer,
            f_stream_set,
            f_stream_set_end,
            f_stream_sol,
            f_stream_sol_end,
            stream_batch,
            stream_batch_end,
            streams_stop,
            stream_sol0,
            merge_in[1],
            a_stream_sol,
            result);

#else
  
    for (int g = 0; g < STOP_S - 2; g++) 
        hls::stream_globals::incr_task_counter();
    
    htb_cache.init();
    bloom_cache.init();

    std::thread mwj_findmin_t(
            mwj_findmin<T_BLOOM,
            BLOOM_LOG,
            K_FUN_LOG>,
            hash1_w,
            qVertices0,
            std::ref(bloom_cache),
            std::ref(dyn_stream_sol),
            std::ref(streams_stop[0]),
            std::ref(p_stream_tuple),
            std::ref(p_stream_filter),
            std::ref(p_stream_sol),
            std::ref(p_stream_sol_end));
        
    std::thread mwj_tuplebuild_t(
            mwj_tuplebuild<PROPOSE_BATCH_LOG>,
            hash1_w,
            hash2_w,
            qVertices0,
            std::ref(b_stream_set),
            std::ref(b_stream_set_end),
            std::ref(b_stream_tuple),
            std::ref(b_stream_sol),
            std::ref(b_stream_sol_end),
            std::ref(streams_stop[2]),
            std::ref(t_stream_tuple),
            std::ref(t_stream_tuple_end),
            std::ref(t_stream_sol),
            std::ref(t_stream_sol_end));

    std::thread mwj_readmin_t(
            mwj_readmin<T_BLOOM,
            BLOOM_LOG,
            K_FUN_LOG>,
            hash1_w,
            hash2_w,
            hTables0,
            htb_buf1,
            std::ref(p_stream_tuple),
            std::ref(p_stream_filter),
            std::ref(streams_stop[1]),
            std::ref(r_stream_set),
            std::ref(r_stream_set_end),
            std::ref(r_stream_tuple));
    
    std::thread mwj_intersect_t(
            mwj_intersect<PROPOSE_BATCH_LOG>,
            hTables1,
            std::ref(htb_cache),
            std::ref(t_stream_tuple),
            std::ref(t_stream_tuple_end),
            std::ref(streams_stop[3]),
            std::ref(i_stream_tuple),
            std::ref(i_stream_tuple_end));

    std::thread mwj_verify_t(
            mwj_verify<PROPOSE_BATCH_LOG>,
            hTables1,
            std::ref(htb_cache),
            std::ref(s_stream_tuple),
            std::ref(s_stream_tuple_end),
            std::ref(streams_stop[4]),
            std::ref(v_stream_tuple),
            std::ref(v_stream_tuple_end));
    
    std::thread mwj_assembly_t(
            mwj_assembly,
            nQueryVer,
            std::ref(f_stream_set),
            std::ref(f_stream_set_end),
            std::ref(f_stream_sol),
            std::ref(f_stream_sol_end),
            std::ref(stream_batch),
            std::ref(stream_batch_end),
            std::ref(streams_stop),
            std::ref(stream_sol0),
            std::ref(merge_in[0]),
            std::ref(a_stream_sol),
            std::ref(result));

#endif /* __SYNTHESIS__ */
    
    hls_thread_local hls::task mwj_merge_t(
            mwj_merge<bool>,
            merge_in,
            merge_out);

    hls_thread_local hls::task mwj_stop_t(
            mwj_stop,
            merge_out,
            dyn_stream_ovf,
            stream_sol0);

#ifndef __SYNTHESIS__

    mwj_findmin_t.join();
    mwj_readmin_t.join();
    mwj_tuplebuild_t.join();
    mwj_intersect_t.join(); 
    mwj_verify_t.join();
    mwj_assembly_t.join();

#if DEBUG_STATS
    debug::cache_hit_prop   = bloom_cache.get_n_hits(0);
    debug::cache_hit_inter  = htb_cache.get_n_l1_hits(1);
    debug::cache_hit_verify = htb_cache.get_n_l1_hits(0);
    debug::cache_req_prop   = bloom_cache.get_n_reqs(0);
    debug::cache_req_inter  = htb_cache.get_n_l1_reqs(1);
    debug::cache_req_verify = htb_cache.get_n_l1_reqs(0);
#endif /* DEBUG_STATS */

    htb_cache.stop();
    bloom_cache.stop();

#endif /* __SYNTHESIS__ */
}

template <typename T_BLOOM,
        size_t BLOOM_LOG, 
        size_t K_FUN_LOG>
void multiwayJoinWrap(
        ap_uint<DDR_W> *htb_buf0,
        ap_uint<DDR_W> *htb_buf1,
        ap_uint<DDR_W> *htb_buf2,
        T_BLOOM *bloom_p,
        row_t *res_buf,
        AdjHT *hTables0,
        AdjHT *hTables1,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        const unsigned short nQueryVer,
        const unsigned char hash1_w,
        const unsigned char hash2_w,
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
#pragma HLS STABLE variable=hTables0
#pragma HLS STABLE variable=hTables1
#pragma HLS STABLE variable=qVertices0
#pragma HLS STABLE variable=qVertices1
#pragma HLS STABLE variable=nQueryVer
#pragma HLS STABLE variable=hash1_w
#pragma HLS STABLE variable=hash2_w
#pragma HLS dataflow

    hls::stream<bool, S_D> stream_batch_end("Stream batch end");
    hls::stream<ap_uint<V_ID_W>, S_D> stream_batch("Stream batch");

    mwj_batch(
        hash1_w,
        hTables0,
        qVertices0,
        htb_buf0,
        stream_batch_end,
        stream_batch);
    
    multiwayJoin<
        T_BLOOM,
        BLOOM_LOG,
        K_FUN_LOG>(
                htb_buf1,
                htb_buf2,
                bloom_p,
                res_buf,
                hTables0,
                hTables1,
                qVertices1,
                nQueryVer,
                hash1_w,
                hash2_w,
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
        row_t bloom_p[BLOOM_SPACE],
        row_t res_buf[RESULTS_SPACE],
        const unsigned short numQueryVert,
        const unsigned char hash1_w,
        const unsigned char hash2_w,
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
    max_widen_bitwidth=128
#pragma HLS INTERFACE mode=m_axi port=htb_buf1 bundle=cache \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2
#pragma HLS INTERFACE mode=m_axi port=htb_buf2 bundle=readmin \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2
#pragma HLS INTERFACE mode=m_axi port=res_buf bundle=fifo \
    max_widen_bitwidth=128 max_read_burst_length=32 max_write_burst_length=32 
#pragma HLS INTERFACE mode=m_axi port=bloom_p bundle=bloom \
    max_widen_bitwidth=128

#pragma HLS alias ports=htb_buf0,htb_buf1,htb_buf2 distance=0

#pragma HLS INTERFACE mode=s_axilite port=numQueryVert
#pragma HLS INTERFACE mode=s_axilite port=hash1_w
#pragma HLS INTERFACE mode=s_axilite port=hash2_w
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

    multiwayJoinWrap<
        bloom_t,
        BLOOM_FILTER_WIDTH,
        K_FUNCTIONS>(
            htb_buf0,
            htb_buf1,
            htb_buf2,
            bloom_p,
            res_buf,
            hTables0,
            hTables1,
            qVertices0,
            qVertices1,
            numQueryVert,
            hash1_w,
            hash2_w,
            dynfifo_diagnostic,
            localResult);

    result = localResult;

#if DEBUG_STATS
    debug::print(hash1_w, hash2_w);
#endif /* DEBUG_STATS */
}

#else

void subgraphIsomorphism(
        edge_t edge_buf[GRAPHS_SPACE],
        row_t htb_buf0[HASHTABLES_SPACE],
        row_t htb_buf1[HASHTABLES_SPACE],
        row_t htb_buf2[HASHTABLES_SPACE],
        row_t bloom_p[BLOOM_SPACE],
        row_t res_buf[RESULTS_SPACE],
        const unsigned short numQueryVert,
        const unsigned short numQueryEdges,
        const unsigned long numDataEdges,
        const unsigned char hash1_w,
        const unsigned char hash2_w,
        unsigned long &dynfifo_diagnostic,

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
    max_widen_bitwidth=128
#pragma HLS INTERFACE mode=m_axi port=htb_buf1 bundle=cache \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2
#pragma HLS INTERFACE mode=m_axi port=htb_buf2 bundle=readmin \
    max_widen_bitwidth=128 num_write_outstanding=1 max_write_burst_length=2
#pragma HLS INTERFACE mode=m_axi port=res_buf bundle=fifo \
    max_widen_bitwidth=128 max_read_burst_length=32 max_write_burst_length=32
#pragma HLS INTERFACE mode=m_axi port=bloom_p bundle=bloom \
    max_widen_bitwidth=128
#pragma HLS INTERFACE mode=m_axi port=edge_buf bundle=graph \
    max_widen_bitwidth=128

#pragma HLS alias ports=htb_buf0,htb_buf1,htb_buf2 distance=0

#pragma HLS INTERFACE mode=s_axilite port=numQueryVert
#pragma HLS INTERFACE mode=s_axilite port=numQueryEdges
#pragma HLS INTERFACE mode=s_axilite port=numDataEdges
#pragma HLS INTERFACE mode=s_axilite port=hash1_w
#pragma HLS INTERFACE mode=s_axilite port=hash2_w
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
        64,
        LABEL_WIDTH,
        DEFAULT_STREAM_DEPTH,
        HASHTABLES_SPACE,
        MAX_QUERY_VERTICES,
        MAX_TABLES>(
                edge_buf,
                htb_buf0,
                bloom_p,
                qVertices0,
                qVertices1,
                hTables0,
                hTables1,
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

    multiwayJoinWrap<
        bloom_t,
        BLOOM_FILTER_WIDTH,
        K_FUNCTIONS>(
            htb_buf0,
            htb_buf1,
            htb_buf2,
            bloom_p,
            res_buf,
            hTables0,
            hTables1,
            qVertices0,
            qVertices1,
            numQueryVert,
            hash1_w,
            hash2_w,
            dynfifo_diagnostic,
            localResult);

    result = localResult;

#if DEBUG_STATS
    debug::print(hash1_w, hash2_w);
#endif /* DEBUG_STATS */

}
#endif /* SOFTWARE_PREPROC */

#pragma GCC diagnostic pop
