#pragma once

/* ------ Definitions regarding data graph ------ */
#define MAX_DATA_VERTICES	(1UL << VERTEX_WIDTH)
#define MAX_DATA_EDGES		(1UL << 12)
#define MAX_OUT_DEGREE      (1UL << COUNTER_WIDTH)

/* ------ Definition regarding query graph ------ */

/* Cannot be bigger than 255 since in the code it is
 * always contained in 8 bits */
#define MAX_QUERY_VERTICES	(1UL << 4)
#define MAX_TABLES 		    32

/* Max space given to query data in fifo for graph read*/
#define MAX_QUERYDATA       300

/* 5 bits are dedicated to labels, to compact as
 * much as possible input data.
 * In the multiway join instead, labels are no more
 * packed with vertex but 1 bit is used to signal
 * the end of a solution */
#define VERTEX_WIDTH        5   /* 4 bytes */
#define VERTEX_WIDTH_BIT    (1UL << VERTEX_WIDTH)
#define EDGE_WIDTH          (VERTEX_WIDTH + 1)
#define LABEL_WIDTH         5

/* Count collision and offsets, it should have enough
 * to contain the number of edges in every table.
 * The final bitwidth is equal to (2^COUNTER_WIDTH-1).
 * The MSB bit is used to check if that hash is used. */
#define COUNTER_WIDTH       5

/* Bloom filter parameters */
#define BLOOM_FILTER_WIDTH  7
#define K_FUNCTIONS         3

/* bitwidth of the hash to index source vertices, 1st level. */
// #define HASH_WIDTH_FIRST    11
// #define HASH_WIDTH_SECOND   7

#define MAX_COLLISIONS      (1UL << 5)

/* Depth of streams connecting tasks */
#define DEFAULT_STREAM_DEPTH 32

/* Dynamic fifo parameters */
#define DYN_FIFO_DEPTH      128
#define DYN_FIFO_BURST      32

#define DDR_BIT             7
#define DDR_WORD            (1UL << DDR_BIT)

#define HASHTABLES_SPACE    ((1UL << 26) / (DDR_WORD / 8))  //~ 67 MB
#define GRAPHS_SPACE        (1UL << 23)                     //~ 8 million edges
#define BLOOM_SPACE         ((1UL << 26) / (DDR_WORD / 8))  //~ 67 MB
#define RESULTS_SPACE		(DYN_FIFO_BURST * (1UL << 17))  //~ 67 MB
// #define HTB_SIZE            (1UL << (HASH_WIDTH_FIRST + HASH_WIDTH_SECOND - (DDR_BIT - COUNTER_WIDTH)))

#define EDGE_ROW            (1UL << (DDR_BIT - EDGE_WIDTH))

#define PROPOSE_BATCH_LOG   6
#define MERGE_IN_STREAMS    2
#define CACHE_WORDS_PER_LINE 3
#include <ap_axi_sdata.h>

/* Functionality definition */
#define COUNT_ONLY          1
#define UNDIRECTED          1
#define DEBUG_INTERFACE     1
#define SOFTWARE_PREPROC    0
#define CACHE_ENABLE        1

#ifndef __SYNTHESIS__
#define DEBUG_STATS         1
#else
#define DEBUG_STATS         0
#endif /*__SYNTHESIS__*/

struct alignas(16) edge_struct{
    ap_uint<VERTEX_WIDTH_BIT> src, dst, labelsrc, labeldst;
};
typedef struct edge_struct edge_t;

typedef ap_axiu<VERTEX_WIDTH_BIT, 0, 0, 0> T_NODE;
typedef ap_axiu<LABEL_WIDTH, 0, 0, 0> T_LABEL;
typedef ap_uint<DDR_WORD> row_t;
typedef ap_uint<(1UL << BLOOM_FILTER_WIDTH)> bloom_t;
typedef ap_uint<(DDR_WORD << CACHE_WORDS_PER_LINE)> edge_block_t;
typedef ap_uint<VERTEX_WIDTH_BIT> vertex_t;

// Useful typename for streams
typedef ap_uint<(1UL << BLOOM_FILTER_WIDTH) * (1UL << K_FUNCTIONS)> bloom_flat_t;
