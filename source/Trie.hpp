#pragma once
#include <ap_int.h>
#include "Parameters.hpp"
#include "types.hpp"

struct AdjHT{
    unsigned int start_offset;
    unsigned int start_edges;
    unsigned int n_edges;
    unsigned int hash_set;
};

struct TableDescriptor{

    /* True: src -> dst, False: dst <- src */
    bool dir;
    ap_uint<LABEL_WIDTH> src_label, dst_label;
};
