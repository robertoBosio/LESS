#include <ap_int.h>
#include "Parameters.hpp"

class QueryVertex {
    public:    
        ap_uint<LABEL_WIDTH> label;            
        
        /* tables in which the vertex is indexed by other */
        ap_uint<8> tables_indexed[MAX_TABLES];
        ap_uint<8> vertex_indexing[MAX_TABLES];
        ap_uint<8> numTablesIndexed = 0;

        /* Tables in which the vertex is used as index */
        ap_uint<8> tables_indexing[MAX_TABLES];   
        ap_uint<8> numTablesIndexing = 0;

        void addTableIndexed(ap_uint<8> t, ap_uint<8> v){
            tables_indexed[numTablesIndexed] = t;
            vertex_indexing[numTablesIndexed] = v;
            numTablesIndexed++;
        }

        void addTableIndexing(ap_uint<8> t){
            tables_indexing[numTablesIndexing++] = t;
        }
};
