#include <fstream>
#include <ap_int.h>
#include "subgraphIsomorphism.hpp"

int main(){

    /* Query data structure */
    ap_uint<8> numQueryVertices = 0;       
    ap_uint<8> *qAdjacency; 
    ap_uint<8> *qOffsets;
    ap_uint<8> *qOrder;
    ap_uint<32> *qLabels;

    std::string fLine{};
    std::cout << "Start" << std::endl;

    /* Max 254 query nodes */
    uint16_t numQueryAdj = 0;

    /* Query files */
    std::ifstream fQueryLab("data/queryLables.txt");
    std::ifstream fQueryOff("data/queryOffsets.txt");
    std::ifstream fQueryCol("data/queryColumns.txt");
    std::ifstream fQueryOrd("data/queryOrder.txt");

    std::getline(fQueryLab, fLine);
    sscanf(fLine.c_str(), "%hhu\n", &numQueryVertices);
    
    qOrder = (ap_uint<8>*)malloc(numQueryVertices);
    qLabels = (ap_uint<8>*)malloc(numQueryVertices);
    qOffsets = (ap_uint<8>*)malloc(numQueryVertices + 1);

    /* Store labels */
    for(int count = 0; count < numQueryVertices; count++){    
        std::getline(fQueryLab, fLine);
        sscanf(fLine.c_str(), "%u %u", &node, &label);
        qLabels[node] = label;
    }

    /* Store offsets */
    getline(fQueryOff, fLine);
    for(int node = 0; node < numQueryVertex + 1; node++){    
        getline (fQueryOff, fLine);
        sscanf(fLine.c_str(), "%hu", &qOffsets[node]);
    }

    /* Store adjancy array */
    getline(fQueryCol, fLine);
    sscanf(fLine.c_str(), "%hu", &numQueryAdj);
    qAdjacency = (ap_uint<8>*)malloc(numQueryAdj);
    for(int count = 0; count < numQueryAdj; count++){    
        getline (fQueryCol, fLine);
        sscanf(fLine.c_str(), "%hhu", &qAdjacency[count]);
    }

    /* Store matching order */
    for(int count = 0; count < numQueryVertices; count++){    
        getline (fQueryOrd, fLine);
        sscanf(fLine.c_str(), "%u", &qOrder[count]);
    }

    fQueryLab.close();
    fQueryOff.close();
    fQueryCol.close();
    fQueryOrd.close();

    subgraphIsomorphism(
            numQueryVertices,
            qOffsets,
            qAdjacency,
            qOrder,
            qLabels);

    std::cout << "End" << std::endl;
    return 0;	
}
