#include <ap_int.h>
#include <hls_stream.h>

#include "Parameters.hpp"
#include "queryVertex.hpp"
#include "Trie.hpp"
#include "hash_lookup3.hpp"

#define DEBUG 1

template <uint8_t V_ID_W, uint8_t V_L_W>
void buildTableDescriptors(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end,
        queryVertex &qVertices,
        TrieDescriptor &tDescriptors,
        ap_uint<8> &numTables,
        
        hls::stream<ap_uint<V_ID_W>> &stream_ord,
        hls::stream<bool> &stream_ord_end)
{

    ap_uint<8> pos = 0;

    /* Filling information about query vertices and recoping
     * the vertex order needed by multiway join */
    bool last = stream_end.read();
    while(!last){
        ap_uint<V_ID_W> node = stream_src.read();
        last = stream_end.read();
        stream_ord.write(node);
        stream_ord_end.write(false);
        qVertices[(unsigned int)node].pos = pos++;
    }
    stream_ord_end.write(true);

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
        for (; g < numTables; g++){
            if (tDescriptors[g].src_label == labelsrc &&
                    tDescriptors[g].dst_label == labeldst &&
                    tDescriptors[g].dir == dirEdge)
                break;
        }
     
        /* Add new table descriptor */   
        if (g == numTables){
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
            qVertices[nodedst].addTableIndexed(g, nodesrc);
        } else {
            qVertices[nodesrc].addTableIndexed(g, nodedst);
            qVertices[nodedst].addTableIndexing(g);
        }

        last = stream_end.read();
    }
}

template <uint8_t V_ID_W, uint8_t V_L_W, uint8_t H_W, uint8_t C_W>
void fillTables(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end,
        Trie &hTables,
        TrieDescriptor &tDescriptors,
        ap_uint<8> numTables)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;

    /* Reset memory */
    for (uint8_t g = 0; g < numTables; g++){
        for(int i = 0; i < (1UL << H_W); i++){
            hTables[g].offset[i] = 0;
        }
    }

    /* Count edges per vertex source */
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

    for (uint8_t g = 0; g < numTables; g++){
        ap_uint<C_W> base_address = 0;
        ap_uint<C_W> step = 0;
        for(int i = 0; i < (1UL << H_W); i++){
            step = hTables[g].offset[i];
            hTables[g].offset[i] = base_address;
            base_address += step;   
        }
    }

    /* From counts to offsets */
    last = stream_end.read();
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
                
                /* Store edge based on offset and update offset */
                hTables[g].edges[hTables[g].offset[index]].write(nodesrc.concat(nodedst));
                hTables[g].offset[index]++;
            }
        }

        last = stream_end.read();
    }

}

template <uint8_t V_ID_W, uint8_t H_W, uint8_t C_W>
void mwj_propose_addrgen(
        Trie &hTables,
        TrieDescriptor &tDescriptors,
        queryVertex &qVertices,
        ap_uint<V_ID_W> current_qv,
        ap_uint<V_ID_W> &current_em,

        hls::stream<SetDescriptor> &stream_set_desc_out,
        hls::stream<bool> &stream_end
        )
{

    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    SetDescriptor curSetDesc;   
    uint8_t tableIndex = 0;
    uint32_t size = 0;

    for(int g = 0; g < qVertices[current_qv].numTablesIndexing; g++){
        tableIndex = qVertices[current_qv].tables_indexing[g];

#if HASH_SET_VERSION
        size = hTables[tableIndex].sCounter;
#else
        size = hTables[tableIndex].offset[(1UL << H_W)-1];
#endif

        curSetDesc.tIndex = tableIndex;
        curSetDesc.sSize = size;
        curSetDesc.sStart = 0;
        curSetDesc.indexed = false;
        stream_set_desc_out.write(curSetDesc);
        stream_end.write(false);
    }

    for(int g = 0; g < qVertices[current_qv].numTablesIndexed; g++){
        tableIndex = qVertices[current_qv].tables_indexed[g];
        uint8_t indexingVertex = qVertices[current_qv].vertex_indexing[g];
        uint8_t ivPos = qVertices[indexingVertex].pos;
        stream_hash_in.write(current_em[ivPos]);
        xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
        ap_uint<H_W> index = stream_hash_out.read().range(H_W - 1, 0);
        ap_uint<C_W> start_off = 0;
        if (index != 0){
            start_off = hTables[tableIndex].offset[index-1];
        }
        size = hTables[tableIndex].offset[index] - start_off;
        
        curSetDesc.tIndex = tableIndex;
        curSetDesc.sSize = size;
        curSetDesc.sStart = start_off;
        curSetDesc.indexed = true;
        stream_set_desc_out.write(curSetDesc);
        stream_end.write(false);
    }
    stream_end.write(true);
}

template <uint8_t V_ID_W, uint8_t H_W, uint8_t C_W>
void mwj_propose_findmin(
        hls::stream<SetDescriptor> &stream_set_desc_in,
        hls::stream<bool> &stream_end_in,

        hls::stream<SetDescriptor> &stream_set_desc_out,
        hls::stream<bool> &stream_end_out,
        SetDescriptor &min_set_desc,
        )
{

    SetDescriptor curSet;

    bool last = stream_end_in.read();
    if (!last){
        min_set_desc = stream_set_desc_in.read();
        last = stream_end_in.read();
    }

    while(!last){
        curSet = stream_set_desc_in.read();

        /* Write the old set descriptor with minimum size and keep
         * as minimum the new one */
        if (curSet.sSize < min_set_desc.sSize){
            SetDescriptor t = min_set_desc;
            min_set_desc = curSet;
            curSet = t;
        }

        stream_set_desc_out.write(curSet);
        stream_end_out.write(false); 
    }
    stream_end_out.write(true);
}

template<uint8_t V_ID_W, uint8_t C_W, uint8_t H_W>
void mvj_propose_readmem(
        hls::stream<SetDescriptor> &stream_set_desc;
        hls::stream<bool> &stream_end_in,
        SetDescriptor &min_set_desc;

#if HASH_SET_VERSION
        hls::stream<ap_uint<H_W>> stream_sets_out,
#else
        hls::stream<ap_uint<V_ID_W>> stream_sets_out,
#endif
        hls::stream<bool> stream_set_ends_out,
        hls::stream<bool> stream_end_out
        )
{

    SetDescriptor curr_desc = min_set_desc;
    bool last = false;
    ap_uint<2*V_ID_W> edge;
    ap_uint<V_ID_W> vertex;
 
#if HASH_SET_VERSION
    ap_uint<H_W> vertex_hash; 
    while(!last){
     
        /* Read one set at a time, starting from the smallest */
        for(int i = curr_desc.sStart; i < (curr_desc.sStart + curr_desc.sSize); i++){
            if (curr_desc.indexed){
                edge = hTables[curr_desc.tIndex].edges[i];
                if (tDescriptors[curr_desc.tIndex].dir){
                    vertex = edge.range(2*V_ID_W-1, V_ID_W);
                } else {
                    vertex = edge.range(V_ID_W-1, 0);
                }
                //manca hashing
            } else {
                vertex_hash = hTables[curr_desc.tIndex].source[i];
            }
            stream_sets_out.write(vertex_hash);
            stream_set_ends_out.write(false);
        }
        stream_set_ends_out.write(true);
        stream_end_out.write(false);

        /* Read next set descriptor */
        last = stream_end_in.read();
        if(!last){
            curr_desc = stream_set_desc.read();
        }
    }
    stream_end_out.write(true);
#else

    while(!last){
     
        /* Read one set at a time, starting from the smallest */
        for(int i = curr_desc.sStart; i < (curr_desc.sStart + curr_desc.sSize); i++){
            edge = hTables[curr_desc.tIndex].edges[i];
            if (tDescriptors[curr_desc.tIndex].dir ^ curr_desc.indexed){
                vertex = edge.range(2*V_ID_W-1, V_ID_W);
            } else {
                vertex = edge.range(V_ID_W-1, 0);
            }
            stream_sets_out.write(vertex);
            stream_set_ends_out.write(false);
        }
        stream_set_ends_out.write(true);
        stream_end_out.write(false);

        /* Read next set descriptor */
        last = stream_end_in.read();
        if(!last){
            curr_desc = stream_set_desc.read();
        }
    }
    stream_end_out.write(true);
#endif
}

mwj_propose(){
}

template <uint8_t V_ID_W>
void mwj_extract(
        hls::stream<ap_uint<V_ID_W>> &buffer_res,
        int debug_counter)
{
    switch(debug_counter){
        case 0:
            {
                buffer_res.write(ap_uint<V_ID_W>(1));
                buffer_res.write(ap_uint<V_ID_W>(2));
                buffer_res.write(ap_uint<V_ID_W>(1));
                buffer_res.write(ap_uint<V_ID_W>(6));
                break;
            }
        case 1:
            {
                buffer_res.write(ap_uint<V_ID_W>(2));
                buffer_res.write(ap_uint<V_ID_W>(2));
                buffer_res.write(ap_uint<V_ID_W>(8));
                break;
            }
        case 2:
            {
                buffer_res.write(ap_uint<V_ID_W>(2));
                buffer_res.write(ap_uint<V_ID_W>(6));
                buffer_res.write(ap_uint<V_ID_W>(5));
                break;
            }
        case 3:
            {
                buffer_res.write(ap_uint<V_ID_W>(3));
                buffer_res.write(ap_uint<V_ID_W>(2));
                buffer_res.write(ap_uint<V_ID_W>(8));
                buffer_res.write(ap_uint<V_ID_W>(0));
                break;
            }
        case 4:
            {
                buffer_res.write(ap_uint<V_ID_W>(3));
                buffer_res.write(ap_uint<V_ID_W>(6));
                buffer_res.write(ap_uint<V_ID_W>(5));
                buffer_res.write(ap_uint<V_ID_W>(7));
                break;
            }
        case 5:
            {
                buffer_res.write(ap_uint<V_ID_W>(4));
                buffer_res.write(ap_uint<V_ID_W>(2));
                buffer_res.write(ap_uint<V_ID_W>(8));
                buffer_res.write(ap_uint<V_ID_W>(0));
                buffer_res.write(ap_uint<V_ID_W>(3));
                break;
            }
        case 6:
            {
                buffer_res.write(ap_uint<V_ID_W>(4));
                buffer_res.write(ap_uint<V_ID_W>(6));
                buffer_res.write(ap_uint<V_ID_W>(5));
                buffer_res.write(ap_uint<V_ID_W>(7));
                buffer_res.write(ap_uint<V_ID_W>(3));
                break;
            }
        default:
            {
                std::cout << "Error: " << debug_counter << std::endl;
            }
    }   
}

template <uint8_t V_ID_W, uint8_t MAX_QV>
void multiwayJoin(
        hls::stream<ap_uint<V_ID_W>> &stream_ord,
        hls::stream<bool> &stream_end,
        ap_uint<8> numTables,
        
        hls::stream<ap_uint<V_ID_W>> &out_embeddings){

    /* Embeddings stored in the buffer in sequential way
     * buffer  <---- 2, v0, v3, 3, v2, v3, v4 ---- */ 
    hls::stream<ap_uint<V_ID_W>> buffer("buffer stream");

    ap_uint<V_ID_W> current_embedding_v[MAX_QV];
    ap_uint<8> current_embedding_c = 0;
    ap_uint<V_ID_W> current_qv = 0;
    ap_uint<8> counter = 0;
    bool no_sol = false;

    int debug_counter = 0;

    bool last_qv = stream_end.read();
    while(!last_qv && !no_sol){
        current_qv = stream_ord.read();
        
        /* Compute until buffer of result is empty or
         * new query vertex is needed */
        while(current_embedding_c == counter && !no_sol){
            //propose
            //intersect
            mwj_extract<V_ID_W>(buffer, debug_counter);
            debug_counter++;

            /* Check if other embedding exist */
            if(!buffer.empty()){    
                current_embedding_c = buffer.read();
                for (int g=0; g < current_embedding_c; g++){
                    current_embedding_v[g] = buffer.read();
                }
            } else {
                no_sol = true;
            }
        }

        counter++;
        last_qv = stream_end.read();
    }
     
#if DEBUG
    std::cout << "count: " << counter << std::endl;
#endif

    for(int g=0; g < counter; g++){
        out_embeddings.write(current_embedding_v[g]);
    }
    while(!buffer.empty()){
        buffer.read();
        for(int g=0; g < counter; g++){
            out_embeddings.write(buffer.read());
        }
    }
}

template <uint8_t V_ID_W, uint8_t V_L_W, uint8_t MAX_QV, uint64_t MAX_TB, uint8_t H_W, uint8_t C_W>
void subgraphIsomorphism(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end,

        hls::stream<ap_uint<V_ID_W>> &stream_out
        )
{

    hls::stream<ap_uint<V_ID_W>> stream_ord("stream order");
    hls::stream<bool> stream_ord_end("stream order end");

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
            &numTables,
            stream_ord,
            stream_ord_end);


    fillTables<V_ID_W, V_L_W, H_W, C_W>(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end,
            hTables,
            tDescriptors,
            numTables);

#if DEBUG

    for (uint8_t g = 0; g < numTables; g++){
        if (tDescriptors[g].dir){
        std::cout << "Table " << (int)g << ": " << (char)tDescriptors[g].src_label << 
             " -> " << (char)tDescriptors[g].dst_label << std::endl;
        } else {
        std::cout << "Table " << (int)g << ": " << (char)tDescriptors[g].dst_label << 
             " <- " << (char)tDescriptors[g].src_label << std::endl;
        }
        for (int s = 0; s < hTables[g].offset[(1UL << H_W)-1]; s++){
            std::cout << hTables[g].edges[s].range(2*V_ID_W-1, V_ID_W) <<
                " " << hTables[g].edges[s].range(V_ID_W-1, 0);
#if HASH_SET_VERSION
            if (s < hTables[g].sCounter)
                std::cout << "\t" << hTables[g].source[s];
#endif
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
#endif

    multiwayJoin<V_ID_W, MAX_QV>(
            stream_ord,
            stream_ord_end,
            numTables,
            stream_out);
}
