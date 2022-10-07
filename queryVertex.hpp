#include <ap_int.h>
#include "Parameters.hpp"

class queryVertex {                    
	ap_uint<8> tables_indexed[MAX_TABLES]; //tables in which the vertex is indexed by other
	ap_uint<8> tables_index[MAX_TABLES];   //tables in which the vertex is used as index
	ap_uint<32> label;                     //label assigned to the vertex 
	ap_uint<8> pos;                        //position in the matching order 

	bool operator<(const queryVertex &r){
		return (this->pos < r.pos);
	}
};
