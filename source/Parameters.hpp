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
#define BLOOM_FILTER_WIDTH  8
#define K_FUNCTIONS         1

/* bitwidth of the hash to index source vertices, 1st level. */
#define HASH_WIDTH_FIRST    11

/* bitwidth of the hash to index a specific edge, 2nd level. */
#define HASH_WIDTH_SECOND   7

#define MAX_COLLISIONS      (1UL << 5)

#define SUB_STREAM_DEPTH    10

#define S_DEPTH             64
#define BURST_S             32
#define DDR_BIT             7
#define DDR_WORD            (1UL << DDR_BIT)
#define HASHTABLES_SPACE    ((1UL << 25) / (DDR_WORD / 8))
#define GRAPHS_SPACE        5000000
#define BLOOM_SPACE         (1UL << HASH_WIDTH_FIRST) * 32
#define RESULTS_SPACE		(BURST_S * (1UL << 17))

#define HTB_SIZE            (1UL << (HASH_WIDTH_FIRST + HASH_WIDTH_SECOND - (DDR_BIT - COUNTER_WIDTH)))
#define EDGE_ROW            (1UL << (DDR_BIT - EDGE_WIDTH))

#define COUNT_ONLY
#define UNDIRECTED
#define DEBUG_INTERFACE
#define INTERSECT_INDEXING_LOOP 0
#define VERIFY_CACHE 1
#define MAX_START_BATCH_SIZE 512
#define PROPOSE_BATCH_LOG 6
#define SET_INFO_WIDTH 16
#define MERGE_IN_STREAMS 2
#define EDGE_BLOCK 4
#include <ap_axi_sdata.h>

typedef struct {
    ap_uint<VERTEX_WIDTH_BIT> src, dst;
} edge_t;

typedef ap_axiu<VERTEX_WIDTH_BIT, 0, 0, 0> T_NODE;
typedef ap_axiu<LABEL_WIDTH, 0, 0, 0> T_LABEL;
typedef ap_uint<DDR_WORD> row_t;
typedef ap_uint<(1UL << BLOOM_FILTER_WIDTH)> bloom_t;
typedef ap_uint<(DDR_WORD << EDGE_BLOCK)> edge_block_t;
