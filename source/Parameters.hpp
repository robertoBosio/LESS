#include <ap_axi_sdata.h>

/* Definitions regarding data graph */
#define MAX_DATA_VERTICES	(1UL << VERTEX_WIDTH)
#define MAX_DATA_EDGES		(1UL << 12)
#define MAX_OUT_DEGREE      (1UL << COUNTER_WIDTH)

/* ------ Definition regarding query graph ------ */

/* Cannot be bigger than 255 since in the code it is
 * always contained in 8 bits */
#define MAX_QUERY_VERTICES	(1UL << 4)
#define MAX_TABLES 		    10

#define VERTEX_WIDTH        4
#define VERTEX_WIDTH_BIT    (1UL << VERTEX_WIDTH)
#define EDGE_WIDTH          (VERTEX_WIDTH + 1)
#define LABEL_WIDTH         8

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
#define HASH_WIDTH_SECOND   5

#define MAX_COLLISIONS      (1UL << 5)

#define STREAM_DEPTH        10
#define STREAM_DEPTH_RES    150000
#define DDR_WIDTH           10000
#define RES_WIDTH			60000

typedef ap_axiu<VERTEX_WIDTH_BIT, 0, 0, 0> T_NODE;
typedef ap_axiu<LABEL_WIDTH, 0, 0, 0> T_LABEL;
