#include <ap_int.h>

#include "Parameters.hpp"
#include "queryVertex.hpp"
#include "Trie.hpp"

#define DEBUG 1
void subgraphIsomorphism(
        ap_uint<8> numQueryVertices,
        ap_uint<8> *queryOffset,
        ap_uint<8> *queryAdjacency,
        ap_uint<8> *queryOrder,
        ap_uint<32> *queryLabels){

#if DEBUG
    for(int g = 0; g < numQueryVertices; g++){
        std::cout << (unsigned char)queryLabels[g] << std::endl;
    }
#endif

}
