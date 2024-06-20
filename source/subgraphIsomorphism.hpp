#pragma once
#define HLS_STREAM_THREAD_SAFE 1

#include <hls_stream.h>
#include <ap_int.h>
#include <ap_axi_sdata.h>
#include "Parameters.hpp"
#include "QueryVertex.hpp"
#include "Trie.hpp"

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
                         QueryVertex qVertices[MAX_QUERY_VERTICES],
                         AdjHT hTables0[MAX_TABLES],
                         AdjHT hTables1[MAX_TABLES],

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

#if COUNT_ONLY
                         long unsigned int &result
#else
                         hls::stream<T_NODE> &result
#endif /* COUNT_ONLY */

);

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
                         unsigned long &p_hmsb0,
                         unsigned long &p_hmsb1,
                         unsigned long &p_hmsb2,
                         unsigned long &p_hmsb3,
#endif /* DEBUG_INTERFACE */

#if COUNT_ONLY
                         long unsigned int &result
#else
                         hls::stream<T_NODE> &result
#endif /* COUNT_ONLY */

);
#endif
