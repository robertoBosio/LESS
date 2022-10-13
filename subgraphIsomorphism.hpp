#include <ap_int.h>
#include <hls_stream.h>

#include "Parameters.hpp"
#include "queryVertex.hpp"
#include "Trie.hpp"
#include "hash_lookup3.hpp"

#include <unordered_map>

#define DEBUG 1

template <uint8_t V_ID_W, uint8_t V_L_W>
void buildTableDescriptors(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end,
        queryVertex *qVertices,
        TrieDescriptor *tDescriptors,
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
            numTables++;
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
        Trie *hTables,
        TrieDescriptor *tDescriptors,
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
        Trie *hTables,
        queryVertex *qVertices,
        ap_uint<V_ID_W> current_qv,
        ap_uint<V_ID_W> *current_em,

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

void mwj_propose_findmin(
        hls::stream<SetDescriptor> &stream_set_desc_in,
        hls::stream<bool> &stream_end_in,

        hls::stream<SetDescriptor> &stream_set_desc_out,
        hls::stream<bool> &stream_end_out,
        SetDescriptor &min_set_desc)
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
        last = stream_end_in.read();
    }
    stream_end_out.write(true);
}

template<uint8_t V_ID_W, uint8_t H_W>
void mwj_propose_readmem(
        hls::stream<SetDescriptor> &stream_set_desc,
        hls::stream<bool> &stream_end_in,
        SetDescriptor &min_set_desc,
        Trie *hTables,
        TrieDescriptor *tDescriptors,

#if HASH_SET_VERSION
        hls::stream<ap_uint<H_W>> &stream_sets_out,
#else
        hls::stream<ap_uint<V_ID_W>> &stream_sets_out,
#endif
        hls::stream<bool> &stream_set_ends_out,
        hls::stream<bool> &stream_end_out)
{

    SetDescriptor curr_desc = min_set_desc;
    bool last = false;
    ap_uint<2*V_ID_W> edge;
    ap_uint<V_ID_W> vertex;
 
#if HASH_SET_VERSION
    ap_uint<H_W> vertex_hash;
    hls::stream<ap_uint<V_ID_W>> stream_hash_in; 
    hls::stream<ap_uint<64>> stream_hash_out; 
#endif

    while(!last){
     
        /* Read one set at a time, starting from the smallest */
        for(int i = curr_desc.sStart; i < (curr_desc.sStart + curr_desc.sSize); i++){
     
#if HASH_SET_VERSION
            
            if (curr_desc.indexed){
                edge = hTables[curr_desc.tIndex].edges[i];
                if (tDescriptors[curr_desc.tIndex].dir){
                    vertex = edge.range(V_ID_W-1, 0);
                } else {
                    vertex = edge.range(2*V_ID_W-1, V_ID_W);
                }
                stream_hash_in.write(vertex);
                xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
                vertex_hash = stream_hash_out.read().range(H_W - 1, 0);
            } else {
                vertex_hash = hTables[curr_desc.tIndex].source[i];
            }
            stream_sets_out.write(vertex_hash);

#else

            edge = hTables[curr_desc.tIndex].edges[i];
            // NEED TO COMPARE FOR COLLISION
            /* Could implement cache for remove already read vertex
             * since in some case we are not reading sets but arrays */
            if (tDescriptors[curr_desc.tIndex].dir ^ curr_desc.indexed){
                vertex = edge.range(2*V_ID_W-1, V_ID_W);
            } else {
                vertex = edge.range(V_ID_W-1, 0);
            }
            stream_sets_out.write(vertex);

#endif

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
}

template <uint8_t V_ID_W, uint8_t H_W, uint8_t C_W>
void mwj_propose(
        Trie *hTables,
        queryVertex *qVertices,
        TrieDescriptor *tDescriptors,
        ap_uint<V_ID_W> curQV,
        ap_uint<V_ID_W> *curEmb,
#if HASH_SET_VERSION
        hls::stream<ap_uint<H_W>> &stream_sets_out,
#else
        hls::stream<ap_uint<V_ID_W>> &stream_sets_out,
#endif
        hls::stream<bool> &stream_set_ends_out,
        hls::stream<bool> &stream_end_out)
{

    hls::stream<SetDescriptor> stream_addrGen("Set descriptors with min");
    hls::stream<bool> stream_addrGen_end("Set descriptors delimiter for addrGen");

    hls::stream<SetDescriptor> stream_findMin("Set descriptors without min");
    hls::stream<bool> stream_findMin_end("Set descriptors delimiter for findMin");
    SetDescriptor minSet;


    mwj_propose_addrgen<V_ID_W, H_W, C_W>(
            hTables,
            qVertices,
            curQV,
            curEmb,
            stream_addrGen,
            stream_addrGen_end);

    mwj_propose_findmin(
            stream_addrGen,
            stream_addrGen_end,
            stream_findMin,
            stream_findMin_end,
            minSet);

    mwj_propose_readmem<V_ID_W, H_W>(
            stream_findMin,
            stream_findMin_end,
            minSet,
            hTables,
            tDescriptors,
            stream_sets_out,
            stream_set_ends_out,
            stream_end_out);
}


template<uint8_t V_ID_W, uint8_t H_W>
void mwj_intersect(
#if HASH_SET_VERSION
    hls::stream<ap_uint<H_W>> &stream_sets_in,
#else
    hls::stream<ap_uint<V_ID_W>> &stream_sets_in,
#endif
    hls::stream<bool> &stream_set_ends_in,
    hls::stream<bool> &stream_end_in,
    
#if HASH_SET_VERSION
    hls::stream<ap_uint<H_W>> &stream_intersection_out,
#else
    hls::stream<ap_uint<V_ID_W>> &stream_intersection_out,
#endif
    hls::stream<bool> &stream_end_out)
{
    bool first = true;
    std::unordered_map<int, bool> hashT;
    bool last = stream_end_in.read();
    while(!last){
        
        bool last_s = stream_set_ends_in.read();
        while(!last_s){

#if HASH_SET_VERSION
            ap_uint<H_W> v = stream_sets_in.read();
#else
            ap_uint<V_ID_W> v = stream_sets_in.read();
#endif
            auto check = hashT.find(v.to_int());
            if (first){
                if (check == hashT.end())
                    hashT.insert(std::make_pair(v.to_int(), false));
            } else {
                if (check != hashT.end())
                    hashT[v.to_int()] = false;
            }
            last_s = stream_set_ends_in.read();
        }
        
        for(auto iter = hashT.begin(); iter != hashT.end();){
            if (iter->second == true){
                iter = hashT.erase(iter);
            } else {
                iter->second = true;
                ++iter;
            }
        }

        first = false;
        last = stream_end_in.read();
    }

    for(auto iter = hashT.begin(); iter != hashT.end(); iter++){
        stream_intersection_out.write(iter->first);
        stream_end_out.write(false);       
    }
    stream_end_out.write(true);
}

template <uint8_t V_ID_W>
void mwj_extract(
        int debug_counter,
        ap_uint<8> nCurEmb,
        ap_uint<V_ID_W> *curEmb,
        hls::stream<ap_uint<V_ID_W>> &stream_intersection_in,
        hls::stream<bool> &stream_end_in,

        hls::stream<ap_uint<V_ID_W>> &buffer_res)
{

    bool last = stream_end_in.read();
    while(!last){
        buffer_res.write(ap_uint<V_ID_W>(nCurEmb+1));
        for (int g = 0; g < nCurEmb; g++){
            buffer_res.write(curEmb[g]);
        }
        buffer_res.write(stream_intersection_in.read());
        last = stream_end_in.read();
    }
}

template <uint8_t V_ID_W, uint8_t MAX_QV, uint8_t H_W, uint8_t C_W>
void multiwayJoin(
        hls::stream<ap_uint<V_ID_W>> &stream_ord,
        hls::stream<bool> &stream_end,
        ap_uint<8> numTables,
        Trie *hTables,
        TrieDescriptor *tDescriptors,
        queryVertex *qVertices,

        hls::stream<ap_uint<V_ID_W>> &out_embeddings,
        hls::stream<bool> &stream_set_end_out,
        hls::stream<bool> &stream_end_out)
{

    /* Embeddings stored in the buffer in sequential way
     * buffer  <---- 2, v0, v3, 3, v2, v3, v4 ---- */ 
    hls::stream<ap_uint<V_ID_W>> buffer("buffer stream");

    /* Propose data out*/
#if HASH_SET_VERSION
    hls::stream<ap_uint<H_W>> p_stream_sets("Sets propose");
#else
    hls::stream<ap_uint<V_ID_W>> p_stream_sets("Sets propose");
#endif
    hls::stream<bool> p_stream_set_ends("Set delimeter propose");
    hls::stream<bool> p_stream_end("End stream propose");

    /* Intersect data out */    
#if HASH_SET_VERSION
    hls::stream<ap_uint<H_W>> i_stream_set;
#else
    hls::stream<ap_uint<V_ID_W>> i_stream_set;
#endif
    hls::stream<bool> i_stream_end;

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
            mwj_propose<V_ID_W, H_W, C_W>(
                    hTables,
                    qVertices,
                    tDescriptors,
                    current_qv,
                    current_embedding_v,
                    p_stream_sets,
                    p_stream_set_ends,
                    p_stream_end);

            mwj_intersect<V_ID_W, H_W>(
                    p_stream_sets,
                    p_stream_set_ends,
                    p_stream_end,
                    i_stream_set,
                    i_stream_end);

            mwj_extract<V_ID_W>(
                    debug_counter,
                    counter,
                    current_embedding_v,
                    i_stream_set,
                    i_stream_end,
                    buffer);

            /* Check if other embedding exist */
            if(!buffer.empty()){    
                current_embedding_c = buffer.read();
                std::cout << "{ ";
                for (int g=0; g < current_embedding_c; g++){
                    current_embedding_v[g] = buffer.read();
                    std::cout << current_embedding_v[g] << " ";
                }
                std::cout << "}" << std::endl;
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
        stream_set_end_out.write(false);
    }
    stream_set_end_out.write(true);
    stream_end_out.write(false);

    while(!buffer.empty()){
        buffer.read();
        for(int g=0; g < counter; g++){
            out_embeddings.write(buffer.read());
            stream_set_end_out.write(false);
        }
        stream_set_end_out.write(true);
        stream_end_out.write(false);
    }
    stream_end_out.write(true);
}

template <uint8_t V_ID_W, uint8_t V_L_W, uint8_t MAX_QV, uint64_t MAX_TB, uint8_t H_W, uint8_t C_W>
void subgraphIsomorphism(
        hls::stream<ap_uint<V_ID_W>> &stream_src,
        hls::stream<ap_uint<V_ID_W>> &stream_dst,
        hls::stream<ap_uint<V_L_W>> &stream_src_l,
        hls::stream<ap_uint<V_L_W>> &stream_dst_l,
        hls::stream<bool> &stream_end,

        hls::stream<ap_uint<V_ID_W>> &stream_out,
        hls::stream<bool> &stream_set_end_out,
        hls::stream<bool> &stream_end_out)
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
            numTables,
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

    multiwayJoin<V_ID_W, MAX_QV, H_W, C_W>(
            stream_ord,
            stream_ord_end,
            numTables,
            hTables,
            tDescriptors,
            qVertices,
            stream_out,
            stream_set_end_out,
            stream_end_out);
}
