/* Definitions regarding data graph */
#define MAX_DATA_VERTICES	(1UL << VERTEX_WIDTH)
#define MAX_DATA_EDGES		(1UL << 12)
#define MAX_OUT_DEGREE      (1UL << COUNTER_WIDTH)

/* ------ Definition regarding query graph ------ */

/* Cannot be bigger than 255 since in the code it is
 * always contained in 8 bits */
#define MAX_QUERY_VERTICES	(1UL << 4)

#define MAX_QUERY_EDGES		(1UL << 5)
#define MAX_QUERY_DEGREE    (1UL << 3)
#define MAX_LABELS		    (1UL << 2)
#define MAX_TABLES 		    (MAX_LABELS*(MAX_LABELS-1))

/* Cannot contains vertex id (2^VERTEX_WIDTH)-1 since the
 * sequence with all the 1s is used to detect unused 
 * position for edges in tables */
#define VERTEX_WIDTH        4
#define VERTEX_WIDTH_BIT    (1UL << VERTEX_WIDTH)
#define EDGE_WIDTH          (VERTEX_WIDTH + 1)
#define LABEL_WIDTH         8

/* bitwidth of the counter inside the hash table,
 * will count collisions, in other words edges which
 * have same hash source and same hash destination.
 * The final bitwidth is equal to (2^COUNTER_WIDTH)-1.
 * The MSB bit is used to check if that hash is used. */
#define COUNTER_WIDTH       4

/* bitwidth of the hash to index source vertices, 1st level.
 * Must be greater or equal to HASH_WIDTH_SECOND. */
#define HASH_WIDTH_FIRST    5

/* bitwidth of the hash to index a specific edge, 2nd level.
 * Must be greater or equal to 5 due to the filter in
 * hashtovid function. */
#define HASH_WIDTH_SECOND   5

#define MAX_COLLISIONS      (1UL << 4)


#define DDR_WIDTH           150
