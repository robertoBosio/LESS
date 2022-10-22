#include <ap_int.h>
#include <hls_stream.h>

#include "Parameters.hpp"
#include "queryVertex.hpp"
#include "Trie.hpp"
#include "hash_lookup3.hpp"

#include <unordered_map>

#define DEBUG 1

bool deb = false;

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

template <uint8_t V_ID_W, uint8_t V_L_W, uint8_t H_W_1, uint8_t H_W_2, uint8_t C_W>
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

    int max_coll = 0;

    /* Reset memory */
    for (int c1 = 0; c1 < numTables; c1++){
        for(int c2 = 0; c2 < (1UL << H_W_1); c2++){
            for (int c3 = 0; c3 < (1UL << H_W_2); c3++){
                hTables[c1].adjHashTable[c2].offset[c3] = 0;
            }
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
                
                ap_uint<V_ID_W> vertexIndexing, vertexIndexed;               
                if(tDescriptors[g].dir){
                    vertexIndexing = nodesrc;
                    vertexIndexed = nodedst;
                } else {
                    vertexIndexing = nodedst;
                    vertexIndexed = nodesrc;
                }
                
                /* Compute index in hash table */
                stream_hash_in.write(vertexIndexing);
                xf::database::hashLookup3<V_ID_W>(
                        stream_hash_in,
                        stream_hash_out);

                ap_uint<H_W_1> indexAdj = 
                    stream_hash_out.read().range(H_W_1 - 1, 0); 

                stream_hash_in.write(vertexIndexed);
                xf::database::hashLookup3<V_ID_W>(
                        stream_hash_in, 
                        stream_hash_out);

                ap_uint<H_W_2> indexEdge = 
                    stream_hash_out.read().range(H_W_2 - 1, 0); 

#if DEBUG_HASH
                std::cout << (int)((tDescriptors[g].dir)?nodesrc:nodedst) << ": " << indexAdj << std::endl;
#endif

                hTables[g].adjHashTable[indexAdj].offset[indexEdge]++;
                if (hTables[g].adjHashTable[indexAdj].offset[indexEdge] > max_coll){
                    max_coll = hTables[g].adjHashTable[indexAdj].offset[indexEdge];
                }
            }
        }
        
        last = stream_end.read();
    }
    
    for (uint8_t iTab = 0; iTab < numTables; iTab++){
        ap_uint<C_W> base_address = 0;
        ap_uint<C_W> step = 0;
        for(unsigned int iAdj = 0; iAdj < (1UL << H_W_1); iAdj++){
            for (unsigned int iEdge = 0; iEdge < (1UL << H_W_2); iEdge++){
                step = hTables[iTab].adjHashTable[iAdj].offset[iEdge];
                hTables[iTab].adjHashTable[iAdj].offset[iEdge] = base_address;
                base_address += step;
            }
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
                
                ap_uint<V_ID_W> vertexIndexing, vertexIndexed;               
                if(tDescriptors[g].dir){
                    vertexIndexing = nodesrc;
                    vertexIndexed = nodedst;
                } else {
                    vertexIndexing = nodedst;
                    vertexIndexed = nodesrc;
                }

                /* Compute index in hash table */
                stream_hash_in.write(vertexIndexing);
                xf::database::hashLookup3<V_ID_W>(
                        stream_hash_in,
                        stream_hash_out);

                ap_uint<H_W_1> indexAdj = 
                    stream_hash_out.read().range(H_W_1 - 1, 0); 

                stream_hash_in.write(vertexIndexed);
                xf::database::hashLookup3<V_ID_W>(
                        stream_hash_in, 
                        stream_hash_out);

                ap_uint<H_W_2> indexEdge = 
                    stream_hash_out.read().range(H_W_2 - 1, 0); 

                /* Store edge based on offset and update offset */
                ap_uint<C_W> off = hTables[g].
                    adjHashTable[indexAdj].
                    offset[indexEdge];

                hTables[g].edges[off] = vertexIndexing.concat(vertexIndexed);
                hTables[g].adjHashTable[(int)indexAdj].offset[(int)indexEdge]++;

            }
        }

        last = stream_end.read();
    }
}

template <uint8_t V_ID_W, uint8_t H_W_1, uint8_t H_W_2, uint8_t C_W>
void mwj_propose(
        Trie *hTables,
        queryVertex *qVertices,
        ap_uint<V_ID_W> curQV,
        ap_uint<V_ID_W> *curEmb,
        
        hls::stream<SetDescriptor> &stream_setdescs_out,
        hls::stream<bool> &stream_end_out)
{

    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    SetDescriptor curSetDesc;   
    uint8_t tableIndex = 0;
    uint32_t size = 0;

    for(int g = 0; g < qVertices[curQV].numTablesIndexing; g++){
        tableIndex = qVertices[curQV].tables_indexing[g];
        size = hTables[tableIndex].adjHashTable[(1UL << H_W_1)-1].offset[(1UL << H_W_2)-1];
        curSetDesc.tIndex = tableIndex;
        curSetDesc.sSize = size;
        curSetDesc.indexed = false;
        curSetDesc.vertexIndexing = (ap_uint<V_ID_W>)0;
        stream_setdescs_out.write(curSetDesc);
        stream_end_out.write(false);
    }

    for(int g = 0; g < qVertices[curQV].numTablesIndexed; g++){
        tableIndex = qVertices[curQV].tables_indexed[g];
        uint8_t indexingVertex = qVertices[curQV].vertex_indexing[g];
        uint8_t ivPos = qVertices[indexingVertex].pos;
        stream_hash_in.write(curEmb[ivPos]);
        xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
        ap_uint<H_W_1> index = stream_hash_out.read().range(H_W_1 - 1, 0);
        ap_uint<C_W> start_off = 0;
        if (index != 0){
            start_off = hTables[tableIndex].adjHashTable[index-1].offset[(1UL << H_W_2)-1];
        }
        size = hTables[tableIndex].adjHashTable[index].offset[(1UL << H_W_2)-1] - start_off;
        
        curSetDesc.tIndex = tableIndex;
        curSetDesc.sSize = size;
        curSetDesc.indexed = true;
        curSetDesc.vertexIndexing = curEmb[ivPos];
        stream_setdescs_out.write(curSetDesc);
        stream_end_out.write(false);
    }
    stream_end_out.write(true);
}


template<uint8_t V_ID_W, uint8_t H_W_1, uint8_t H_W_2, uint8_t C_W>
void mwj_intersect_streammin(
    hls::stream<SetDescriptor> &stream_set_descs_in,
    hls::stream<bool> &stream_end_in,
    SetDescriptor *setDesc,
    uint8_t &counter,
    Trie *hTables,
    TrieDescriptor *tDescriptors,

    hls::stream<ap_uint<H_W_1>> &stream_min_out,
    hls::stream<bool> &stream_end_out)
{
    counter = 0;
    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    ap_uint<C_W> start_off = 0;
    ap_uint<C_W> end_off;
    ap_uint<V_ID_W*2> edge;
    ap_uint<V_ID_W> vertex, vertexCheck;

    bool last = stream_end_in.read();
    if (!last){
        setDesc[counter++] = stream_set_descs_in.read();
        last = stream_end_in.read();
    }

    while(!last){
        setDesc[counter] = stream_set_descs_in.read();

        /* Write the old set descriptor with minimum size and keep
         * as minimum the new one */
        if (setDesc[counter].sSize < setDesc[0].sSize){
            SetDescriptor t = setDesc[0];
            setDesc[0] = setDesc[counter];
            setDesc[counter] = t;
        }
        counter++;
        last = stream_end_in.read();
    }
  
    /* Read the smallest set */
    if (setDesc[0].indexed){
        
        stream_hash_in.write(setDesc[0].vertexIndexing);
        xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
        ap_uint<H_W_1> index = stream_hash_out.read().range(H_W_1-1, 0);

        if (index != 0){
            start_off = hTables[setDesc[0].tIndex].
                adjHashTable[index-1].
                offset[(1UL << H_W_2) - 1];
        }
        end_off = hTables[setDesc[0].tIndex].
                adjHashTable[index].
                offset[(1UL << H_W_2) - 1];
     
        /* Read bag of indexed vertices and write hashes */
        for(ap_uint<C_W> i = start_off; i < end_off; i++){
            edge = hTables[setDesc[0].tIndex].edges[i];
            
            /* Vertex indexing the table */
            vertexCheck = edge.range(2*V_ID_W-1, V_ID_W);

            vertex = edge.range(V_ID_W-1, 0);

            /* Check vertex indexing to remove collisions and thus transform
             * a bag into a set */
            if (setDesc[0].vertexIndexing == vertexCheck){
                stream_hash_in.write(vertex);
                xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
                ap_uint<H_W_1> vertexHash = stream_hash_out.read().range(H_W_1-1, 0);
                stream_min_out.write(vertexHash);
                stream_end_out.write(false);
            }
            
        }

    } else {
     
        /* Read set of source hashes */
        for(uint32_t i = 0; i < (1UL << H_W_1); i++){
            end_off = hTables[setDesc[0].tIndex].
                adjHashTable[(ap_uint<H_W_1>)i].
                offset[(1UL << H_W_2)-1];
            
            if (end_off > start_off){
                stream_min_out.write((ap_uint<H_W_1>)i);
                stream_end_out.write(false);
                start_off = end_off;
            }
        }
    }
    stream_end_out.write(true);
}


template <uint8_t V_ID_W, uint8_t H_W_1, uint8_t H_W_2, uint8_t C_W>
void mwj_intersect_probe(
    hls::stream<ap_uint<H_W_1>> &stream_min_in,
    hls::stream<bool> &stream_end_in,
    SetDescriptor *setDesc,
    uint8_t &counter,
    Trie *hTables,

    hls::stream<ap_uint<H_W_1>> &stream_inter_out,
    hls::stream<bool> &stream_end_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    bool last = stream_end_in.read();
    ap_uint<H_W_1> candidate;

    while(!last){
        candidate = stream_min_in.read();
        bool inter = true;
        for (int i = 1; i < counter && inter == true; i++){
            ap_uint<C_W> start_off = 0;
            ap_uint<C_W> end_off = 0;

            if (setDesc[i].indexed) {

                /* Retriving hash for vertex indexing the table */
                stream_hash_in.write(setDesc[i].vertexIndexing);
                xf::database::hashLookup3<V_ID_W>(
                        stream_hash_in, 
                        stream_hash_out);

                ap_uint<H_W_1> index1 = 
                    stream_hash_out.read().range(H_W_1 - 1, 0);
                
                ap_uint<H_W_2> index2 = candidate.range(H_W_2-1, 0);

                /* Handling corner cases */
                if (index2 != 0){
                    start_off = hTables[setDesc[i].tIndex].
                        adjHashTable[index1].
                        offset[index2-1];
                } else {
                    if (index1 != 0){
                        start_off = hTables[setDesc[i].tIndex].
                            adjHashTable[index1-1].
                            offset[(1UL << H_W_2)-1];
                    }
                }

                /* Retriving end address of edges with
                 * vertex indexing -> candidate */
                end_off = hTables[setDesc[i].tIndex].
                    adjHashTable[index1].
                    offset[index2];
 
            } else {

                /* Handling corner cases */
                if (candidate != 0){
                    start_off = hTables[setDesc[i].tIndex].
                        adjHashTable[candidate-1].
                        offset[(1UL << H_W_2)-1];
                }
                end_off = hTables[setDesc[i].tIndex].
                    adjHashTable[candidate].
                    offset[(1UL << H_W_2)-1];
            }                
           
            /* if start address and end address are equal
             * the candidate is not present in the set */
            if (end_off == start_off){
                inter = false;
            }

        }

        if (inter){
            stream_inter_out.write(candidate);
            stream_end_out.write(false);
        }

        last = stream_end_in.read();
    }
    stream_end_out.write(true);
}

template <uint8_t V_ID_W, uint8_t H_W_1, uint8_t H_W_2, uint8_t C_W, uint8_t MAX_QD>
void mwj_intersect(
        hls::stream<SetDescriptor> &stream_setdescs_in,
        hls::stream<bool> &stream_end_in,
        Trie *hTables,
        TrieDescriptor *tDescriptors,
        SetDescriptor &minSet,

        hls::stream<ap_uint<H_W_1>> &stream_inter_out,
        hls::stream<bool> &stream_end_out)
{
    SetDescriptor setDescs[MAX_QD];
    uint8_t counter;
    hls::stream<ap_uint<H_W_1>> stream_min;
    hls::stream<bool> stream_end_min;
   
    mwj_intersect_streammin<V_ID_W, H_W_1, H_W_2, C_W>(
            stream_setdescs_in,
            stream_end_in,
            setDescs,
            counter,
            hTables,
            tDescriptors,
            stream_min,
            stream_end_min);
    
    mwj_intersect_probe<V_ID_W, H_W_1, H_W_2, C_W>(
            stream_min,
            stream_end_min,
            setDescs,
            counter,
            hTables,
            stream_inter_out,
            stream_end_out);
    
    minSet = setDescs[0];
}

template <uint8_t V_ID_W, uint8_t H_W_1, uint8_t H_W_2, uint8_t C_W>
void mwj_extract_hashtovid(
        hls::stream<ap_uint<H_W_1>> &stream_inter_in,
        hls::stream<bool> &stream_end_in,
        SetDescriptor minSet,
        Trie *hTables,

        hls::stream<ap_uint<V_ID_W>> &stream_inter_out,
        hls::stream<bool> &stream_set_end_out,
        hls::stream<bool> &stream_end_out)
{

    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    ap_uint<V_ID_W*2> edge;
    bool last = stream_end_in.read();
    ap_uint<H_W_1> hashinter;
    ap_uint<V_ID_W> vertex;

    /* filter to remove hash collission when passing from H_W_1 to H_W_2 */
    ap_uint<16> filter[(1UL << (H_W_2 - 4))];
    for (int g = 0; g < (1UL << (H_W_2 - 4)); filter[g++] = 0);

    while(!last){
        hashinter = stream_inter_in.read();
        if (deb) {
            std::cout << hashinter << std::endl;
        }
        
        ap_uint<C_W> start_off = 0;
        ap_uint<C_W> end_off = 0;
        if (minSet.indexed) {

            /* Retriving hash for vertex indexing the table */
            stream_hash_in.write(minSet.vertexIndexing);
            xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
            ap_uint<H_W_1> index1 = stream_hash_out.read().range(H_W_1 - 1, 0);
            ap_uint<H_W_2> index2 = hashinter.range(H_W_2 - 1, 0);
            
            /* Handling corner cases */
            if (index2 != 0){
                start_off = hTables[minSet.tIndex].
                    adjHashTable[index1].
                    offset[index2-1];
            } else {
                if (index1 != 0){
                    start_off = hTables[minSet.tIndex].
                        adjHashTable[index1-1].
                        offset[(1UL << H_W_2)-1];
                }
            }

            /* Retriving end address of edges with
             * vertex indexing -> candidate */
            end_off = hTables[minSet.tIndex].
                adjHashTable[index1].
                offset[index2];
            
            /* If hash bit already present do no stream anything */
            ap_uint<4> bitindex = index2.range(3, 0);
            ap_uint<H_W_2 - 4> arrayindex = index2.range(H_W_2-1, 4);
            ap_uint<16> filterword = filter[arrayindex];
            if (filterword[bitindex] == 1){
                end_off = start_off;
            } else {
                filterword[bitindex] = 1;
                filter[arrayindex] = filterword;
            }

        } else {

            /* Handling corner cases */
            if (hashinter != 0){
                start_off = hTables[minSet.tIndex].
                    adjHashTable[hashinter-1].
                    offset[(1UL << H_W_2)-1];
            }
            end_off = hTables[minSet.tIndex].
                adjHashTable[hashinter].
                offset[(1UL << H_W_2)-1];
        }

        for (; start_off < end_off; start_off++){
            edge = hTables[minSet.tIndex].edges[start_off];
            if (minSet.indexed){
                vertex = edge.range(V_ID_W-1, 0);
            } else {
                vertex = edge.range(2*V_ID_W-1, V_ID_W);
            }
            stream_inter_out.write(vertex);
            stream_set_end_out.write(false);            
        }
        stream_set_end_out.write(true);
        stream_end_out.write(false);
        last = stream_end_in.read();
    }
    stream_end_out.write(true);
}

template <uint8_t V_ID_W, uint8_t H_W_1, uint8_t H_W_2, uint8_t C_W, uint8_t MAX_CL>
void mwj_extract_bagtoset(
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_set_end_in,
        hls::stream<bool> &stream_end_in,

        hls::stream<ap_uint<V_ID_W>> &stream_inter_out,
        hls::stream<bool> &stream_end_out)
{
    
    ap_uint<V_ID_W> set[MAX_CL];
    uint8_t counter = 0;
    bool last = stream_end_in.read();
    while(!last){
        bool lastSet = stream_set_end_in.read();
        counter = 0;
        while(!lastSet){
            ap_uint<V_ID_W> vertex = stream_inter_in.read();
            if (deb) {
                std::cout << vertex << std::endl;
            }
            uint8_t nSet = 0;
            for(; nSet < counter; nSet++){
                if (vertex == set[nSet]){
                    break;
                }
            }
            if (nSet == counter){
                set[counter++] = vertex;
                stream_inter_out.write(vertex);
                stream_end_out.write(false);
            }
            lastSet = stream_set_end_in.read();    
        }
        if (deb) {
            std::cout << std::endl;
        }
        last = stream_end_in.read();
    }
    if (deb) {
        std::cout << std::endl;
    }
    stream_end_out.write(true);
}

template <uint8_t V_ID_W>
void mwj_extract_product(
        ap_uint<8> nCurEmb,
        ap_uint<V_ID_W> *curEmb,
        hls::stream<ap_uint<V_ID_W>> &stream_intersection_in,
        hls::stream<bool> &stream_end_in,

        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_embed_end_out,
        hls::stream<bool> &stream_end_out)
{
    ap_uint<V_ID_W> vertex;
    bool last = stream_end_in.read();
    while(!last){
        for (int g = 0; g < nCurEmb; g++){
            stream_embed_out.write(curEmb[g]);
            stream_embed_end_out.write(false);
        }
        vertex = stream_intersection_in.read();
        stream_embed_out.write(vertex);
        if (deb) {
            std::cout << vertex << std::endl;
        }
        stream_embed_end_out.write(false);
        stream_embed_end_out.write(true);
        stream_end_out.write(false);
        last = stream_end_in.read();
    }
    stream_end_out.write(true);
}

template <uint8_t V_ID_W, uint8_t H_W_1, uint8_t H_W_2, uint8_t C_W, uint8_t MAX_CL>
void mwj_extract(
        hls::stream<ap_uint<H_W_1>> &stream_inter_in,
        hls::stream<bool> &stream_end_in,
        Trie *hTables,
        SetDescriptor &minSet,
        ap_uint<8> nCurEmb,
        ap_uint<V_ID_W> *curEmb,

        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_embed_end_out,
        hls::stream<bool> &stream_end_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_htv_inter("Bags of vertices");
    hls::stream<bool> stream_htv_bag_end("Bag delimeter");
    hls::stream<bool> stream_htv_end("Bags delimeter");

    hls::stream<ap_uint<V_ID_W>> stream_bts("Set of embeddings");
    hls::stream<bool> stream_bts_end;

    mwj_extract_hashtovid<V_ID_W, H_W_1, H_W_2, C_W>(
            stream_inter_in,
            stream_end_in,
            minSet,
            hTables,
            stream_htv_inter,
            stream_htv_bag_end,
            stream_htv_end);

    mwj_extract_bagtoset<V_ID_W, H_W_1, H_W_2, C_W, MAX_CL>(
            stream_htv_inter,
            stream_htv_bag_end,
            stream_htv_end,
            stream_bts,
            stream_bts_end);

    mwj_extract_product<V_ID_W>(
            nCurEmb,
            curEmb,
            stream_bts,
            stream_bts_end,
            stream_embed_out,
            stream_embed_end_out,
            stream_end_out);
}

template <uint8_t V_ID_W, uint8_t H_W_1, uint8_t H_W_2, uint8_t C_W, uint8_t MAX_QV>
void mwj_verify(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_embed_end_in,
        hls::stream<bool> &stream_end_in,
        Trie *hTables,
        queryVertex *qVertices,
        ap_uint<V_ID_W> curQV,
        ap_uint<V_ID_W> *curEmb,
        
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_embed_end_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    bool last = stream_end_in.read();
    while(!last){
        ap_uint<V_ID_W> curEmb[MAX_QV];
        int counter = 0;
        bool checked = true;

        bool lastEmb = stream_embed_end_in.read();
        while(!lastEmb){    
            curEmb[counter++] = stream_embed_in.read();
            lastEmb = stream_embed_end_in.read();
        }

        for(int g = 0; g < qVertices[curQV].numTablesIndexed && checked; g++){
            checked = false;
            uint8_t tableIndex = qVertices[curQV].tables_indexed[g];
            uint8_t indexingVertex = qVertices[curQV].vertex_indexing[g];
            uint8_t ivPos = qVertices[indexingVertex].pos;
            stream_hash_in.write(curEmb[ivPos]);
            xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
            ap_uint<H_W_1> index1 = stream_hash_out.read().range(H_W_1 - 1, 0);

            stream_hash_in.write(curEmb[counter-1]);
            xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
            ap_uint<H_W_2> index2 = stream_hash_out.read().range(H_W_2 - 1, 0);
            
            ap_uint<C_W> start_off = 0;

            /* Handling corner cases */
            if (index2 != 0){
                start_off = hTables[tableIndex].
                    adjHashTable[index1].
                    offset[index2-1];
            } else {
                if (index1 != 0){
                    start_off = hTables[tableIndex].
                        adjHashTable[index1-1].
                        offset[(1UL << H_W_2)-1];
                }
            }

            /* Retriving end address of edges with
             * vertex indexing -> candidate */
            ap_uint<C_W> end_off = hTables[tableIndex].
                adjHashTable[index1].
                offset[index2];

            if (deb){
                std::cout << "Table " << (int)tableIndex << ": searching for " <<
                    curEmb[ivPos] << " (" << index1 << ") -> " <<
                    curEmb[counter -1] << "(" << index2 << ")" << 
                    "starting from " << start_off << " to " <<
                    end_off << std::endl;
            }
            for(; start_off < end_off; start_off++){
                ap_uint<2*V_ID_W> edge =
                    hTables[tableIndex].edges[start_off];
                ap_uint<V_ID_W> vertexIndexed, vertexIndexing;
                vertexIndexed = edge.range(V_ID_W-1, 0);
                vertexIndexing = edge.range(2*V_ID_W-1, V_ID_W);
                if (deb){
                    std::cout << "\t" << vertexIndexing << " -> " <<
                        vertexIndexed << std::endl;
                }
                if (vertexIndexing == curEmb[ivPos] &&
                        vertexIndexed == curEmb[counter-1]){
                    checked = true;
                    break;
                }
            }
        }

        if (checked){
            for (int g = 0; g < counter; g++){
                stream_embed_out.write(curEmb[g]);
                stream_embed_end_out.write(false);
            }
            if (deb) {
                std::cout << curEmb[counter-1] << ": OK" << std::endl; 
            }
            stream_embed_end_out.write(true);
        } else {
            if (deb) {
                std::cout << curEmb[counter-1] << ": NO" << std::endl; 
            }
        
        }
        last = stream_end_in.read();
    }
}

template <uint8_t V_ID_W, uint8_t MAX_QV, uint8_t H_W_1, uint8_t H_W_2, uint8_t C_W, uint8_t MAX_QD, uint8_t MAX_CL>
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
    hls::stream<ap_uint<V_ID_W>> stream_embed("Embeddings stream");
    hls::stream<bool> stream_embed_end("Embedding delimeter");

    /* Propose data out*/
    hls::stream<SetDescriptor> p_stream_setdescs("Set descs propose");
    hls::stream<bool> p_stream_end("End stream propose");

    /* Intersect data out */    
    hls::stream<ap_uint<H_W_1>> i_stream_set;
    hls::stream<bool> i_stream_end;
    
    /* Extract data out */
    hls::stream<ap_uint<V_ID_W>> e_stream_embed("Embeddings extract");
    hls::stream<bool> e_stream_embed_end("Embedding extract delimeter");
    hls::stream<bool> e_stream_end("Extract delimeter");

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
            SetDescriptor minSet;

            mwj_propose<V_ID_W, H_W_1, H_W_2, C_W>(
                    hTables,
                    qVertices,
                    current_qv,
                    current_embedding_v,
                    p_stream_setdescs,
                    p_stream_end);

            mwj_intersect<V_ID_W, H_W_1, H_W_2, C_W, MAX_QD>(
                    p_stream_setdescs,
                    p_stream_end,
                    hTables,
                    tDescriptors,
                    minSet,
                    i_stream_set,
                    i_stream_end);

            mwj_extract<V_ID_W, H_W_1, H_W_2, C_W, MAX_CL>(
                    i_stream_set,
                    i_stream_end,
                    hTables,
                    minSet,
                    counter,
                    current_embedding_v,
                    e_stream_embed,
                    e_stream_embed_end,
                    e_stream_end);

            mwj_verify<V_ID_W, H_W_1, H_W_2, C_W, MAX_QV>(
                    e_stream_embed,
                    e_stream_embed_end,
                    e_stream_end,
                    hTables,
                    qVertices,
                    current_qv,
                    current_embedding_v,
                    stream_embed,
                    stream_embed_end);
            
           
            /* Check if other embedding exist */
            if(!stream_embed_end.empty()){    
                current_embedding_c = 0;
                bool last = stream_embed_end.read();
                while(!last){
                    current_embedding_v[current_embedding_c++] = 
                        stream_embed.read();
                    last = stream_embed_end.read();
                }
            } else {
                no_sol = true;
            }
            debug_counter++;
        }
        std::cout << "counter :" << counter << std::endl;
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

    while(!stream_embed_end.empty()){
        bool last = stream_embed_end.read();
        while(!last){
            out_embeddings.write(stream_embed.read());
            stream_set_end_out.write(false);
            last = stream_embed_end.read();
        }
        stream_set_end_out.write(true);
        stream_end_out.write(false);
    }
    stream_end_out.write(true);
}

template <uint8_t V_ID_W, uint8_t V_L_W, uint8_t MAX_QV, uint64_t MAX_TB, uint8_t H_W_1, uint8_t H_W_2, uint8_t C_W, uint8_t MAX_QD, uint8_t MAX_CL>
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
    TrieDescriptor tDescriptors[MAX_TB];
    ap_uint<8> numTables = 0;

#ifndef __SYNTHESIS__ 
    Trie *hTables;
    hTables = (Trie*)malloc(sizeof(Trie)*MAX_TB);
#else
    Trie hTables[MAX_TB];
#endif

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


    fillTables<V_ID_W, V_L_W, H_W_1, H_W_2, C_W>(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            stream_end,
            hTables,
            tDescriptors,
            numTables);

#ifdef DEBUG_TAB
    for(int g = 0; g < numTables; g++){
        if (tDescriptors[g].dir){
            std::cout << "Table " << (int)g << ": " << (char)tDescriptors[g].src_label << 
                " -> " << (char)tDescriptors[g].dst_label << std::endl;
        } else {
            std::cout << "Table " << (int)g << ": " << (char)tDescriptors[g].dst_label << 
                " <- " << (char)tDescriptors[g].src_label << std::endl;
        }
        int start = 0;
        for (int s = 0; s < (1UL << H_W_1); s++){
            std::cout << s << ": " << std::endl;
            for (int ss = 0; ss < (1UL << H_W_2); ss++){
                std::cout << "\t" << ss << ":" << std::endl;
                for(; start < hTables[g].adjHashTable[s].offset[ss]; start++){
                    std::cout << "\t\t" << 
                        hTables[g].edges[start].range(2*V_ID_W-1, V_ID_W) <<
                        hTables[g].edges[start].range(V_ID_W-1, 0) << std::endl;
                }
            }
        }
        std::cout << std::endl;
    }
#endif

    multiwayJoin<V_ID_W, MAX_QV, H_W_1, H_W_2, C_W, MAX_QD, MAX_CL>(
            stream_ord,
            stream_ord_end,
            numTables,
            hTables,
            tDescriptors,
            qVertices,
            stream_out,
            stream_set_end_out,
            stream_end_out);

#ifndef __SYNTHESIS__ 
    free(hTables);
#endif
}
