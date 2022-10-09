#include <ap_int.h>
#include <hls_stream.h>

#include "Parameters.hpp"
#include "queryVertex.hpp"
#include "Trie.hpp"
#include "hash_lookup3.hpp"

#define DEBUG 1
#define HASH_SET_VERSION 0

template <uint8_t V_ID_W, uint8_t V_L_W>
void buildTableDescriptors(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end,
        queryVertex *qVertices,
        TrieDescriptor *tDescriptors,
        ap_uint<8> *numTables)
{

    ap_uint<8> pos = 0;

    /* Filling information about query vertices */
    bool last = stream_end.read();
    while(!last){
        ap_uint<V_ID_W> node = stream_src.read();
        last = stream_end.read();
        qVertices[(unsigned int)node].pos = pos++;
    }

    last = stream_end.read();
    while(!last){
        bool dirEdge = false;
        ap_uint<V_ID_W> nodesrc = stream_src.read();
        ap_uint<V_ID_W> nodedst = stream_dst.read();
        ap_uint<V_L_W> labelsrc = stream_src_l.read();
        ap_uint<V_L_W> labeldst = stream_dst_l.read();
        
        if (qVertices[nodesrc] < qVertices[nodedst])
            dirEdge = true;

#if DEBUG
        std::cout << (unsigned int)nodesrc << "(" << (char)labelsrc << ")"
            << " -> " <<  (unsigned int)nodedst << "(" << (char)labeldst << ")" 
            << std::endl;
#endif
        /* Understanding if table already exists */
        uint8_t g = 0;
        for (; g < *numTables; g++){
            if (tDescriptors[g].src_label == labelsrc &&
                    tDescriptors[g].dst_label == labeldst &&
                    tDescriptors[g].dir == dirEdge)
                break;
        }
     
        /* Add new table descriptor */   
        if (g == *numTables){
            tDescriptors[g].src_label = labelsrc;
            tDescriptors[g].dst_label = labeldst;
            tDescriptors[g].dir = dirEdge;
            (*numTables)++;
        }

#if DEBUG
        if (dirEdge){
        std::cout << "Table " << (int)g << ": " << (char)labelsrc << 
             " -> " << (char)labeldst << std::endl;
        } else {
        std::cout << "Table " << (int)g << ": " << (char)labeldst << 
             " <- " << (char)labelsrc << std::endl;
        }
#endif

        /* Linking vertices to tables */
        if (dirEdge){
            qVertices[nodesrc].addTableIndexing(g);
            qVertices[nodedst].addTableIndexed(g);
        } else {
            qVertices[nodesrc].addTableIndexed(g);
            qVertices[nodedst].addTableIndexing(g);
        }

        last = stream_end.read();
    }
}

template <uint8_t V_ID_W, uint8_t V_L_W, uint8_t H_W>
void fillTables(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end,
        Trie *hTables,
        TrieDescriptor *tDescriptors,
        ap_uint<8> numTables)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;

    bool last = stream_end.read();
    
    while(!last){
        ap_uint<V_ID_W> nodesrc = stream_src.read();
        ap_uint<V_ID_W> nodedst = stream_dst.read();
        ap_uint<V_L_W> labelsrc = stream_src_l.read();
        ap_uint<V_L_W> labeldst = stream_dst_l.read();
        
        /* Finding correct table */
        uint8_t g = 0;
        for (; g < numTables; g++){
            if (tDescriptors[g].src_label == labelsrc &&
                    tDescriptors[g].dst_label == labeldst){
                
                if(tDescriptors[g].dir){
                    stream_hash_in.write(nodesrc);
                } else {
                    stream_hash_in.write(nodedst);
                }

                /* Compute index in hash table */
                xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
                ap_uint<H_W> index = stream_hash_out.read().range(H_W - 1, 0); 
#if DEBUG
                std::cout <<  ((tDescriptors[g].dir)?nodesrc:nodedst) << ": " << index << std::endl;
#endif

#if HASH_SET_VERSION
                if (hTables[g].offset[index] == 0)
                    hTables[g].addSourceVertex(index);
#endif
                hTables[g].offset[index]++;
            }
        }

        last = stream_end.read();
    }



}
template <uint8_t V_ID_W, uint8_t V_L_W, uint8_t MAX_QV, uint64_t MAX_TB, uint8_t H_W>
void subgraphIsomorphism(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end
        )
{

    queryVertex qVertices[MAX_QV];
    Trie hTables[MAX_TB];
    TrieDescriptor tDescriptors[MAX_TB];
    ap_uint<8> numTables = 0;

#if DEBUG
    std::cout << "Allocating " << sizeof(Trie)*MAX_TB << " bytes." << std::endl;
#endif

    buildTableDescriptors<V_ID_W, V_L_W>(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end,
            qVertices,
            tDescriptors,
            &numTables);


    fillTables<V_ID_W, V_L_W, H_W>(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end,
            hTables,
            tDescriptors,
            numTables);

}
