#include <ap_int.h>
#include "Parameters.hpp"
#include "types.hpp"

struct AdjHT{
    ap_uint<64> start_offset;
    ap_uint<64> start_edges;
    ap_uint<32> n_edges;
    ap_uint<32> hash_set;
};

struct TableDescriptor{
    
    /* True: src -> dst, False: dst <- src */
    bool dir;
    ap_uint<LABEL_WIDTH> src_label, dst_label;
};
