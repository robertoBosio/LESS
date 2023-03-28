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

#define VERTEX_WIDTH        5   /* 4 bytes */
#define VERTEX_WIDTH_BIT    (1UL << VERTEX_WIDTH)
#define EDGE_WIDTH          (VERTEX_WIDTH + 1)
#define LABEL_WIDTH         5

/* Count collision and offsets, it should have enough
 * to contain the number of edges in every table.
 * The final bitwidth is equal to (2^COUNTER_WIDTH-1).
 * The MSB bit is used to check if that hash is used. */
#define COUNTER_WIDTH       4

/* bitwidth of the hash to index source vertices, 1st level.
 * Must be greater or equal to HASH_WIDTH_SECOND. */
#define HASH_WIDTH_FIRST    9

/* bitwidth of the hash to index a specific edge, 2nd level.
 * Must be greater or equal to 5 due to the filter in
 * hashtovid function. */
#define HASH_WIDTH_SECOND   6

#define MAX_COLLISIONS      (1UL << 5)

#define STREAM_DEPTH        10

#define S_DEPTH             32
#define BURST_S             32
#define DDR_BIT             7
#define DDR_WORD            (1UL << DDR_BIT)
#define HASHTABLES_SPACE    (15000000 / (DDR_WORD / 8))
#define GRAPHS_SPACE        500000
#define RESULTS_SPACE		(BURST_S * (1UL << 17))

#define HTB_SIZE            (1UL << (HASH_WIDTH_FIRST + HASH_WIDTH_SECOND - (DDR_BIT - COUNTER_WIDTH)))
#define EDGE_ROW            (1UL << (DDR_BIT - EDGE_WIDTH))

#define COUNT_ONLY
#define UNDIRECTED
#define DEBUG_INTERFACE
#define INTERSECT_INDEXING_LOOP 0

#include <ap_axi_sdata.h>
typedef struct {
    ap_uint<VERTEX_WIDTH_BIT> src, dst;
} edge_t;

typedef ap_axiu<VERTEX_WIDTH_BIT, 0, 0, 0> T_NODE;
typedef ap_axiu<LABEL_WIDTH, 0, 0, 0> T_LABEL;
typedef ap_uint<DDR_WORD> row_t;
