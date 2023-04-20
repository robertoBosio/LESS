#include "subisoWrap.hpp"
#include "subgraphIsomorphism.hpp"

void subisoWrapper(
        edge_t edge_buf[GRAPHS_SPACE],
        row_t htb_buf[HASHTABLES_SPACE],
        /* hls::burst_maxi<row_t> htb_buf1, */
        row_t res_buf[RESULTS_SPACE],
        unsigned short numQueryVert,
        unsigned short numQueryEdges,
        unsigned long numDataEdges,
        unsigned short numBatchSize,

#ifdef DEBUG_INTERFACE
        volatile unsigned int &debif_endpreprocess,
#endif

#ifdef COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif

        )
{

    subgraphIsomorphism
            (edge_buf,
             htb_buf,
             htb_buf,
             htb_buf,
             /* htb_buf1, */
             res_buf,
             numQueryVert,
             numQueryEdges,
             numDataEdges,
             numBatchSize,
#ifdef DEBUG_INTERFACE
             debif_endpreprocess,
#endif
             result);
}

