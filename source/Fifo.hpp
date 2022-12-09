#include <ap_int.h>
#include "Parameters.hpp"
#include <hls_stream.h>

template<T, unsigned int MAX>
class Fifo {
    ap_uint<512> *start;
    hls::stream<T, MAX> int_stream;
    unsigned int count;
    unsigned int count_mem
    public:    
    
        void Fifo(ap_uint<512> *mem){
            start = mem;
            count = 0;
        }

        void write(T data){
            if (count == MAX){
                count_mem++;
            } else {
                int_stream.write(data);
                count++;
            }
        }

        T read(){
            if (count_mem != 0){
            } else {
                count--;
                return int_stream.read();
            }
        }
};
