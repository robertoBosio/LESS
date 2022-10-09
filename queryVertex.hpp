#include <ap_int.h>
#include "Parameters.hpp"

class queryVertex {
    public:    
        ap_uint<8> tables_indexed[MAX_TABLES]; //tables in which the vertex is indexed by other
        ap_uint<8> tables_indexing[MAX_TABLES];   //tables in which the vertex is used as index
        ap_uint<8> numTablesIndexed = 0;
        ap_uint<8> numTablesIndexing = 0;
        ap_uint<LABEL_WIDTH> label;            //label assigned to the vertex 
        ap_uint<8> pos;                        //position in the matching order 

        bool operator<(const queryVertex &r){
            return (this->pos < r.pos);
        }

        void addTableIndexed(ap_uint<8> t){
            tables_indexed[numTablesIndexed++] = t;
        }

        void addTableIndexing(ap_uint<8> t){
            tables_indexing[numTablesIndexing++] = t;
        }
};
