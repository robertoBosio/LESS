#pragma once
#include <ap_int.h>
#include "Parameters.hpp"

class QueryVertex {
    public:    
        //ap_uint<LABEL_WIDTH> label;            
        
        /* tables in which the vertex is indexed by other */
        unsigned char tables_indexed[MAX_TABLES];
        unsigned char vertex_indexing[MAX_TABLES];
        unsigned char numTablesIndexed = 0;

        /* Tables in which the vertex is used as index */
        unsigned char tables_indexing[MAX_TABLES];   
        unsigned char numTablesIndexing = 0;

        void addTableIndexed(unsigned char t, unsigned char v){
            tables_indexed[numTablesIndexed] = t;
            vertex_indexing[numTablesIndexed] = v;
            numTablesIndexed++;
        }

        void addTableIndexing(unsigned char t){
            tables_indexing[numTablesIndexing++] = t;
        }
};
