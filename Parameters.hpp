/* Definitions regarding data graph */
#define MAX_DATA_VERTICES	(1UL << VERTEX_WIDTH)
#define MAX_DATA_EDGES		2^20
#define MAX_OUT_DEGREE      (1UL << COUNTER_WIDTH)

/* Definition regarding query grap */
#define MAX_QUERY_VERTICES	2^4
#define MAX_QUERY_EDGES		2^5
#define MAX_QUERY_DEGREE    2^3
#define MAX_LABELS		    2^4
#define MAX_TABLES 		    (MAX_LABELS*(MAX_LABELS-1))

#define VERTEX_WIDTH        16
#define LABEL_WIDTH         8

/* bitwidth of the counter inside the hash table,
 * will count collisions, in other words edges which
 * have same hash source and same hash destination */
#define COUNTER_WIDTH       10

/* bitwidth of the hash to index source vertices, 1st level */
#define HASH_WIDTH_FIRST    1

/* bitwidth of the hash to index a specific edge, 2nd level */
#define HASH_WIDTH_SECOND   1

#define MAX_COLLISIONS      2^3
