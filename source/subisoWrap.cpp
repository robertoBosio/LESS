#include "subisoWrap.hpp"
#include "subgraphIsomorphism.hpp"

void subisoWrapper(
		hls::stream<T_NODE> &stream_src,
		hls::stream<T_NODE> &stream_dst,
		hls::stream<T_LABEL> &stream_src_l,
		hls::stream<T_LABEL> &stream_dst_l,
        ap_uint<512> htb_buf[DDR_WIDTH],
        T_DDR res_buf[RES_WIDTH],

#ifdef DEBUG_INTERFACE
        unsigned int &debug_endpreprocess_s,
#endif

#ifdef COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif

        )
{

    subgraphIsomorphism
            (stream_src,
             stream_dst,
             stream_src_l,
             stream_dst_l,
             htb_buf,
             htb_buf,
             htb_buf,
             htb_buf,
             htb_buf,
             res_buf,
#ifdef DEBUG_INTERFACE
             debug_endpreprocess_s,
#endif
             result);
}

