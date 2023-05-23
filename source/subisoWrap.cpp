#include "Parameters.hpp"
#include "subisoWrap.hpp"
#include "subgraphIsomorphism.hpp"

#if SOFTWARE_PREPROC
void subisoWrapper(
        row_t htb_buf[HASHTABLES_SPACE],
        bloom_t bloom_p[BLOOM_SPACE],
        row_t res_buf[RESULTS_SPACE],
        unsigned short numQueryVert,
        unsigned short numBatchSize,
        unsigned long &diagnostic,
        QueryVertex qVertices0[MAX_QUERY_VERTICES], 
        QueryVertex qVertices1[MAX_QUERY_VERTICES],
        AdjHT hTables0[MAX_TABLES],
        AdjHT hTables1[MAX_TABLES],

#if DEBUG_INTERFACE
        volatile unsigned int &debif_endpreprocess,
#endif

#if COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif
        )
{

    subgraphIsomorphism
            (htb_buf,
             htb_buf,
             htb_buf,
             htb_buf,
             bloom_p,
             res_buf,
             numQueryVert,
             numBatchSize,
             diagnostic,
             qVertices0,
             qVertices1,
             hTables0,
             hTables1,
#if DEBUG_INTERFACE
             debif_endpreprocess,
#endif
             result);
}

#else

void subisoWrapper(
        edge_t edge_buf[GRAPHS_SPACE],
        row_t htb_buf[HASHTABLES_SPACE],
        bloom_t bloom_p[BLOOM_SPACE],
        row_t res_buf[RESULTS_SPACE],
        unsigned short numQueryVert,
        unsigned short numQueryEdges,
        unsigned long numDataEdges,
        unsigned short numBatchSize,
        unsigned long &diagnostic,

#if DEBUG_INTERFACE
        volatile unsigned int &debif_endpreprocess,
#endif

#if COUNT_ONLY
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
             htb_buf,
             bloom_p,
             res_buf,
             numQueryVert,
             numQueryEdges,
             numDataEdges,
             numBatchSize,
             diagnostic,
#ifdef DEBUG_INTERFACE
             debif_endpreprocess,
#endif
             result);
}
#endif /* SOFTWARE_PREPROC */
