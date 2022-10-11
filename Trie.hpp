#include <ap_int.h>
#include "Parameters.hpp"

class Trie{

    public:
        ap_uint<COUNTER_WIDTH> offset[1UL << HASH_WIDTH];
        ap_uint<2*VERTEX_WIDTH> edges[MAX_DATA_EDGES];

#if HASH_SET_VERSION
        ap_uint<HASH_WIDTH> source[1UL << HASH_WIDTH];
        ap_uint<HASH_WIDTH> sCounter;
        void addSourceVertex(ap_uint<HASH_WIDTH> hash){
            source[sCounter++] = hash;
        }
#endif

};

class TrieDescriptor{
    
    public:
        /* True: src -> dst, False: dst <- src */
        bool dir;
        ap_uint<LABEL_WIDTH> src_label, dst_label;
};

class SetDescriptor{
    public:
        uint8_t tIndex;
        uint32_t sSize;
        uint32_t sStart;
        bool indexed;

};
