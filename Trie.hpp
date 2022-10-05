#include <ap_int.h>
#include "Parameters.hpp"

class Trie{

	bool dir;
	ap_uint<32> src_label, adj_label;
	ap_uint<64> source[MAX_DATA_VERTICES];    
	ap_uint<64> offset[MAX_DATA_VERTICES];
	ap_uint<64> adjacency[MAX_DATA_EDGES];
}
