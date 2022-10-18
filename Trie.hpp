#include <ap_int.h>
#include "Parameters.hpp"

struct hashTable{
    ap_uint<COUNTER_WIDTH> offset[1UL << HASH_WIDTH_SECOND];

};

struct Trie{
    hashTable adjHashTable[1UL << HASH_WIDTH_FIRST];
    ap_uint<2*VERTEX_WIDTH> edges[MAX_DATA_EDGES];
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
        bool indexed;
        ap_uint<VERTEX_WIDTH> vertexIndexing;
};
