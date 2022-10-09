/* Definitions regarding data graph */
#define MAX_DATA_VERTICES	(1UL << VERTEX_WIDTH)
#define MAX_DATA_EDGES		2^20
#define MAX_OUT_DEGREE      (1UL << COUNTER_WIDTH)

/* Definition regarding query grap */
#define MAX_QUERY_VERTICES	2^4
#define MAX_QUERY_EDGES		2^5
#define MAX_LABELS		    2^4
#define MAX_TABLES 		    (MAX_LABELS*(MAX_LABELS-1))

#define VERTEX_WIDTH        16
#define LABEL_WIDTH         8
#define COUNTER_WIDTH       18
#define HASH_WIDTH          16

