/* Definitions regarding data graph */
#define MAX_DATA_VERTICES	(1UL << VERTEX_WIDTH)
#define MAX_DATA_EDGES		(1UL << 12)
#define MAX_OUT_DEGREE      (1UL << COUNTER_WIDTH)

/* Definition regarding query grap */
#define MAX_QUERY_VERTICES	(1UL << 4)
#define MAX_QUERY_EDGES		(1UL << 5)
#define MAX_QUERY_DEGREE    (1UL << 3)
#define MAX_LABELS		    (1UL << 2)
#define MAX_TABLES 		    (MAX_LABELS*(MAX_LABELS-1))

#define VERTEX_WIDTH        16
#define LABEL_WIDTH         8

/* bitwidth of the counter inside the hash table,
 * will count collisions, in other words edges which
 * have same hash source and same hash destination */
#define COUNTER_WIDTH       18

/* bitwidth of the hash to index source vertices, 1st level */
#define HASH_WIDTH_FIRST    12

/* bitwidth of the hash to index a specific edge, 2nd level */
#define HASH_WIDTH_SECOND   5

#define MAX_COLLISIONS      (1UL << 4)
