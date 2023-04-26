
template <typename T,
            size_t TUPLE_I,
            size_t TUPLE_V,
            size_t BATCH_SIZE_LOG>
void mwj_intersect(
        AdjHT                         *hTables,
        T                             htb_buf,
        hls::stream<ap_uint<TUPLE_I>> &stream_tuple_in,
        hls::stream<bool>             &stream_tuple_end_in,
        hls::stream<ap_uint<V_ID_W>>  &stream_sol_in,
        hls::stream<bool>             &stream_sol_end_in,
        hls::stream<bool>             &stream_stop,

        hls::stream<ap_uint<TUPLE_V>> &stream_tuple_out,
        hls::stream<bool>             &stream_tuple_end_out,
        hls::stream<ap_uint<V_ID_W>>  &stream_sol_out,
        hls::stream<bool>             &stream_sol_end_out)
{
    ap_uint<64> candidate_hash;
    ap_uint<(1UL << BATCH_SIZE_LOG)> bits;
    ap_uint<V_ID_W> candidate_v;
    ap_uint<V_ID_W> indexing_v;
    ap_uint<TUPLE_I> tuple_in;
    ap_uint<TUPLE_V> tuple_out;
    ap_uint<SET_INFO_WIDTH> set_info;
    unsigned char tableIndex;
    ap_uint<DDR_W> ram_row;
    unsigned long addr_row;
    unsigned long addr_inrow;
    ap_uint<64> addr_counter;
    bool stop, last;

    while (1) {
        if (stream_sol_end_in.read_nb(last)){
INTERSECT_COPYING_EMBEDDING_LOOP:
            while(!last){
                stream_sol_out.write(stream_sol_in.read());
                stream_sol_end_out.write(false);
                last = stream_sol_end_in.read();
            }
            stream_sol_end_out.write(true);

            bits = ~0;
            last = stream_tuple_end_in.read();
INTERSECT_LOOP:
            while(!last){
                tuple_in = stream_tuple_in.read();
                ap_uint<64> hash_out;
                ap_uint<(1UL << C_W)> start_off = 0;
                ap_uint<(1UL << C_W)> end_off = 1;
                ap_uint<BATCH_SIZE_LOG> pos = tuple_in.range(TUPLE_I - 3, V_ID_W * 2 + 8);

                if (tuple_in.test(TUPLE_I - 1) && bits[pos]){

                    ap_uint<H_W_1> index1_f;
                    ap_uint<H_W_2> index2_f;
                    candidate_v = tuple_in.range(V_ID_W - 1, 0);
                    indexing_v = tuple_in.range((V_ID_W * 2) - 1, V_ID_W);
                    tableIndex = tuple_in.range((V_ID_W * 2) + 7, V_ID_W * 2);

                    xf::database::details::hashlookup3_core<V_ID_W>(
                            candidate_v, candidate_hash);
                    xf::database::details::hashlookup3_core<V_ID_W>(
                            indexing_v, hash_out);
                    ap_uint<H_W_1> index1 = hash_out.range(H_W_1 - 1, 0);
                    ap_uint<H_W_2> index2 = candidate_hash.range(H_W_2 - 1, 0);

                    addr_counter = index1;
                    addr_counter <<= H_W_2;
                    addr_counter += index2;

                    /* Compute address of row storing the counter */
                    addr_row = hTables[tableIndex].start_offset + 
                        (addr_counter >> (DDR_BIT - C_W));

                    /* Compute address of data inside the row */
                    addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);
                    ram_row = htb_buf.get(addr_row, 0);
                    end_off = ram_row.range(((addr_inrow + 1) << C_W) - 1, 
                            addr_inrow << C_W);
                    
                    if (addr_counter != 0){
                        addr_counter--;

                        /* Compute address of row storing the counter */
                        addr_row = hTables[tableIndex].start_offset + 
                            (addr_counter >> (DDR_BIT - C_W));

                        /* Compute address of data inside the row */
                        addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);
                        ram_row = htb_buf.get(addr_row, 0);
                        start_off = ram_row.range(((addr_inrow + 1) << C_W) - 1, 
                                addr_inrow << C_W);
#ifdef DEBUG_STATS
                        debug::intersect_reads++;
#endif
                    }

#ifdef DEBUG_STATS
                    debug::intersect_reads++;
                    debug::intersect_filter += (start_off < end_off)? 1 : 0;
#endif
                }

                // Building the tuple as 
                // ((TUPLE_I), START, END, INTER_BIT)
                tuple_out[TUPLE_V - 1] = (start_off < end_off) & bits[pos];
                tuple_out.range(TUPLE_I - 1, 0) = tuple_in;
                tuple_out.range(TUPLE_I + (1UL << C_W) - 1, TUPLE_I) = start_off;
                tuple_out.range(TUPLE_I + 2 * (1UL << C_W) - 1,
                       TUPLE_I + (1UL << C_W)) = end_off;

                bits[pos] = bits[pos] & (start_off < end_off);
                stream_tuple_out.write(tuple_out);
                stream_tuple_end_out.write(false);

/* std::cout << "(( " */
/* << tuple_out.range(V_ID_W - 1, 0) << ", " */
/* << tuple_out.range((V_ID_W * 2) - 1, V_ID_W) << ", " */
/* << tuple_out.range((V_ID_W * 2) + 7, V_ID_W * 2) << ", " */
/* << tuple_out.range(TUPLE_I - 3, V_ID_W * 2 + 8) << ", " */
/* << tuple_out[TUPLE_I - 2] << ", " */
/* << tuple_out[TUPLE_I - 1] << "), " */
/* << tuple_out.range(TUPLE_I + (1UL << C_W) - 1, TUPLE_I) << ", " */
/* << tuple_out.range(TUPLE_I + 2 * (1UL << C_W) - 1, */
/* TUPLE_I + (1UL << C_W)) << ", " */
/* << tuple_out[TUPLE_V - 1] << ") " */
/* << bits[pos] << std::endl; */

                last = stream_tuple_end_in.read();
            }
            stream_tuple_end_out.write(true);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

template <size_t LOCK_WORD_LOG>
bool checknlock(
        ap_uint<(1UL << LOCK_WORD_LOG)> &locks,
        ap_uint<LOCK_WORD_LOG> idx_lock){
#pragma HLS inline
    bool test = locks.test(idx_lock);
    if (!test){
        locks[idx_lock] = 1;
    }
    return test;
}

template <size_t LOCK_WORD_LOG>
void slowMutex(
        hls::stream<ap_uint<LOCK_WORD_LOG>> &checknlock_stream,
        hls::stream<bool> &sem_stream,
        hls::stream<ap_uint<LOCK_WORD_LOG>> &unlock_stream
        ){
    static ap_uint<LOCK_WORD_LOG> checknlock_mutex;
    static ap_uint<(1UL << LOCK_WORD_LOG)> locks = 0;
    static bool checknlock_flag {false};

    ap_uint<LOCK_WORD_LOG> unlock_mutex;
    bool test {true};

#pragma HLS pipeline II=1
    bool unlock_flag = unlock_stream.read_nb(unlock_mutex);
    if (!checknlock_flag){
        checknlock_flag = checknlock_stream.read_nb(checknlock_mutex);
    }

    if (checknlock_flag)
        test = checknlock<LOCK_WORD_LOG>(locks, checknlock_mutex);
    
    if (unlock_flag)
        locks[unlock_mutex] = 0;

    ap_wait();

    if (!test){
        sem_stream.write(true);
        checknlock_flag = false;
    }
}

/* Builds the table descriptors based on the information
 * from the query graph. */
void buildTableDescriptors(
        hls::stream<T_NODE> &stream_src,
        hls::stream<T_NODE> &stream_dst,
        hls::stream<T_LABEL> &stream_src_l,
        hls::stream<T_LABEL> &stream_dst_l,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        TableDescriptor *tDescriptors,
        ap_uint<8> &numTables,
        ap_uint<8> &numQueryVert)
{
    /* Translate from id of vertex to position in the order */
    ap_uint<8> fromNumToPos[MAX_QV];

    /* Filling information about query vertices and coping
     * the vertex order needed by multiway join */
    bool last;
FILL_ORDER_LOOP:
    do{
        T_NODE nodeif = stream_src.read();
        ap_uint<V_ID_W> node = nodeif.data;
#ifndef __SYNTHESIS__
        assert(numQueryVert < (MAX_QV));
#endif
        fromNumToPos[node] = numQueryVert;
        numQueryVert++;
        last = nodeif.last;
    } while(!last);

    numQueryVert--;
    
    /* Creating table descriptors */
CREATE_TABDESC_LOOP:
    do{
        bool dirEdge = false;
        T_LABEL labeldstif = stream_dst_l.read();
        T_LABEL labelsrcif = stream_src_l.read();
        T_NODE nodedstif = stream_dst.read();
        T_NODE nodesrcif = stream_src.read();
        ap_uint<V_L_W> labeldst = labeldstif.data;
		ap_uint<V_L_W> labelsrc = labelsrcif.data;
		ap_uint<V_ID_W> nodedst = nodedstif.data;
		ap_uint<V_ID_W> nodesrc = nodesrcif.data;
        ap_uint<8> nodeSrcPos = fromNumToPos[nodesrc];
        ap_uint<8> nodeDstPos = fromNumToPos[nodedst];

        if (nodeSrcPos < nodeDstPos)
            dirEdge = true;

#ifndef __SYNTHESIS__
        std::cout << (unsigned int)nodesrc << "(" << (int)labelsrc << ")"
            << " -> " <<  (unsigned int)nodedst << "(" << (int)labeldst << ")" 
            << std::endl;
#endif

        /* Understanding if table already exists */
        uint8_t g = 0;
FIND_CORRECT_TABLE_LOOP:
        for (; g < numTables; g++){
            if (tDescriptors[g].src_label == labelsrc 
                    && tDescriptors[g].dst_label == labeldst 
                    && tDescriptors[g].dir == dirEdge)
                break;
        }
     
        /* Add new table descriptor */   
        if (g == numTables){
#ifndef __SYNTHESIS__
            assert(numTables < MAX_TB);
#endif
            tDescriptors[g].src_label = labelsrc;
            tDescriptors[g].dst_label = labeldst;
            tDescriptors[g].dir = dirEdge;
            numTables++;
        }

#ifndef __SYNTHESIS__
        if (dirEdge){
        std::cout << "Table " << (int)g << ": " << (int)labelsrc << 
             " -> " << (int)labeldst << std::endl;
        } else {
        std::cout << "Table " << (int)g << ": " << (int)labeldst << 
             " <- " << (int)labelsrc << std::endl;
        }
#endif

        /* Linking vertices to tables */
        if (dirEdge){
            qVertices0[nodeSrcPos].addTableIndexing(g);
            qVertices1[nodeSrcPos].addTableIndexing(g);
            qVertices0[nodeDstPos].addTableIndexed(g, nodeSrcPos);
            qVertices1[nodeDstPos].addTableIndexed(g, nodeSrcPos);
        } else {
            qVertices0[nodeSrcPos].addTableIndexed(g, nodeDstPos);
            qVertices0[nodeDstPos].addTableIndexing(g);
            qVertices1[nodeSrcPos].addTableIndexed(g, nodeDstPos);
            qVertices1[nodeDstPos].addTableIndexing(g);
        }

        last = nodesrcif.last;
    } while(!last);
}

template <size_t WORD,
         size_t COUNTER>
ap_uint<(1UL << WORD)> utils_increasecounter(
        ap_uint<(1UL << WORD)> row,
        unsigned int pos)
{
#pragma HLS inline
    const unsigned int CR = (1UL << (WORD - COUNTER));
    ap_uint <(1UL << WORD)> mask = row >> (pos << COUNTER);
    mask += 1;
    mask <<= (pos << COUNTER);
    row <<= ((CR - pos) << COUNTER);
    row >>= ((CR - pos) << COUNTER);
    return row |= mask;
}


/* Translating counter addresses to ram addresses and
 * increasing of one the counter specified.
 * The counter is on 16 bits, so shift by 4. */
void increaseCounter(
        hls::stream<ap_uint<(1UL << E_W)>> &stream_edge,
        hls::stream<ap_uint<H_W_1>> &stream_index1,
        hls::stream<ap_uint<H_W_2>> &stream_index2,
        hls::stream<ap_uint<8>> &stream_ntable,
        hls::stream<bool> &stream_end,
        AdjHT *hTables,
        ap_uint<DDR_W> *htb_buf)
{
    ap_uint<DDR_W> ram_row;
    ap_uint<H_W_1> index1;
    ap_uint<H_W_2> index2;
    ap_uint<(1UL << E_W)> edge;
    ap_uint<8> ntb;
    ap_uint<64> addr_outrow;
    unsigned int addr_inrow;
    ap_uint<64> addr_counter;
    ap_uint<(1UL << C_W)> counter;

    bool last = stream_end.read();
INCREASE_COUNTER_DDR_LOOP:
    while(!last){
        index1 = stream_index1.read();
        index2 = stream_index2.read();
        ntb = stream_ntable.read();
        edge = stream_edge.read();

        /* Address in the matrix [H_W_1][H_W_2] */
        addr_counter = index1;
        addr_counter <<= H_W_2;
        addr_counter += index2;

        /* Compute address of row storing the counter */
        addr_outrow = hTables[ntb].start_offset + (addr_counter >> (DDR_BIT - C_W));

        /* Compute address of counter inside the row */
        addr_inrow = addr_counter.range((DDR_BIT - C_W) - 1, 0);

        /* Read, modify and write the counter */
        htb_buf[addr_outrow] = utils_increasecounter<DDR_BIT, C_W>(htb_buf[addr_outrow], addr_inrow);
 
        hTables[ntb].n_edges++;     
        last = stream_end.read();
    }
}

/* Transform counters to offsets. */
void counterToOffset(
        ap_uint<8> numTables,
        AdjHT *hTables,
        ap_uint<DDR_W> *htb_buf)
{
    ap_uint<DDR_W> ram_row, ram_row_new;

COUNTER_TO_OFFSET_DDR_LOOP:
    for (unsigned int ntb = 0; ntb < numTables; ntb++){
        ap_uint<(1 << C_W)> base_addr = 0;
        ap_uint<(1UL << C_W)> prev_addr = 0;
        ap_uint<(1UL << C_W)> counter;
        ap_uint<32> hash_used = 0;

COUNTER_TO_OFFSET_TABLE_LOOP:
        for(unsigned int start = 0; start < HTB_SIZE; start++){
            ram_row = htb_buf[start + hTables[ntb].start_offset];
            for (int g = 0; g < CNT_ROW; g++){
                counter = ram_row.range((1UL << C_W) - 1, 0);
                ram_row_new >>= (1UL << C_W);
                ram_row >>= (1UL << C_W);

#ifndef __SYNTHESIS__
                assert(base_addr < (1UL << (1UL << (C_W))-1));
#endif /* __SYNTHESIS__ */
#ifdef DEBUG_STATS
                if (counter > debug::max_collisions) 
                    debug::max_collisions = counter;
                debug::avg_collisions += (float)counter / (1 << (H_W_1 + H_W_2));
#endif /* DEBUG_STATS */

                ram_row_new.range((CNT_ROW << C_W) - 1, (CNT_ROW - 1) << C_W) = base_addr;
                base_addr += counter;

                /* Store used bit */
                if (counter > 0){
                    ram_row_new.set((CNT_ROW << C_W) - 1);
                }

                /* Check if an hash is used */
                if (((start << (DDR_BIT - C_W)) + g) % (1UL << H_W_2) == 0){
                    if(prev_addr < base_addr){
                        hash_used++;
                        prev_addr = base_addr;
                    }
                }
            }
            htb_buf[start + hTables[ntb].start_offset] = ram_row_new;
        }
        hTables[ntb].hash_set = hash_used;
    }       
}

/* Store edges based on offsets */
void storeEdges(
        hls::stream<ap_uint<(1UL << E_W)>> &stream_edge,
        hls::stream<ap_uint<H_W_1>> &stream_index1,
        hls::stream<ap_uint<H_W_2>> &stream_index2,
        hls::stream<ap_uint<8>> &stream_ntable,
        hls::stream<bool> &stream_end,
        AdjHT *hTables,
        ap_uint<DDR_W> *htb_buf)
{
    ap_uint<DDR_W> ram_row;
    ap_uint<H_W_1> index1;
    ap_uint<H_W_2> index2;
    ap_uint<(1UL << E_W)> edge;
    ap_uint<8> ntb;
    ap_uint<64> addr_outrow;
    ap_uint<64> addr_inrow;
    ap_uint<64> addr_offset;
    ap_uint<(1UL << C_W) - 1> offset;

    bool last = stream_end.read();
INCREASE_COUNTER_DDR_LOOP:
    while(!last){
        index1 = stream_index1.read();
        index2 = stream_index2.read();
        ntb = stream_ntable.read();
        edge = stream_edge.read();

        /* Address in the matrix [H_W_1][H_W_2] */
        addr_offset = index1;
        addr_offset <<= H_W_2;
        addr_offset += index2;

        /* Compute address of row storing the offset */
        addr_outrow = hTables[ntb].start_offset + (addr_offset >> (DDR_BIT - C_W));

        /* Compute address of offset inside the row */
        addr_inrow = addr_offset.range((DDR_BIT - C_W) - 1, 0);

        /* Read, modify and write the offset */
        /* ram_row = htb_buf[addr_outrow]; */
        /* offset = ram_row.range(((addr_inrow + 1) << C_W) - 2, addr_inrow << C_W); */
        /* ram_row.range(((addr_inrow + 1) << C_W) - 2, */
        /* addr_inrow << C_W) = offset + 1; */

        ram_row = htb_buf[addr_outrow];
        ap_uint <DDR_W> mask = ram_row >> (addr_inrow << C_W);
        offset = mask & ((1UL << (1UL << C_W)) - 1);
        mask += 1;
        mask <<= (addr_inrow << C_W);
        ram_row <<= ((CNT_ROW - addr_inrow) << C_W);
        ram_row >>= ((CNT_ROW - addr_inrow) << C_W);
        ram_row |= mask;
        htb_buf[addr_outrow] = ram_row;
        
        /* Compute address of row that will store the edge */
        addr_outrow = hTables[ntb].start_edges 
            + (offset >> (DDR_BIT - E_W));

        /* Compute address of the edge inside the row */
        addr_inrow = offset.range((DDR_BIT - E_W) - 1, 0);

        /* Read, modify and write the edge */
        ram_row = htb_buf[addr_outrow];
/* ram_row.range(((addr_inrow + 1) << E_W) - 1, */
/* addr_inrow << E_W) = edge; */
        mask = edge;
        mask <<= (addr_inrow << E_W);
        ram_row |= mask;
        htb_buf[addr_outrow] = ram_row;

        last = stream_end.read();
    }
}

void edgeToHash(
        hls::stream<T_NODE> &stream_src,
        hls::stream<T_NODE> &stream_dst,
        hls::stream<T_LABEL> &stream_src_l,
        hls::stream<T_LABEL> &stream_dst_l,
        TableDescriptor *tDescriptors,
        ap_uint<8> numTables,

        hls::stream<ap_uint<(1UL << E_W)>> &stream_edge,
        hls::stream<ap_uint<H_W_1>> &stream_index1,
        hls::stream<ap_uint<H_W_2>> &stream_index2,
        hls::stream<ap_uint<8>> &stream_ntable,
        hls::stream<bool> &stream_end_out)
{
    
    hls::stream<ap_uint<V_ID_W>> stream_hash_in;
    hls::stream<ap_uint<64>> stream_hash_out;
    
    bool last;
COUNT_OCCURENCIES_LOOP:
    do {
        T_LABEL labeldstif = stream_dst_l.read();
        T_LABEL labelsrcif = stream_src_l.read();
        T_NODE nodedstif = stream_dst.read();
        T_NODE nodesrcif = stream_src.read();
        ap_uint<V_L_W> labeldst = labeldstif.data;
		ap_uint<V_L_W> labelsrc = labelsrcif.data;
		ap_uint<V_ID_W> nodedst = nodedstif.data;
		ap_uint<V_ID_W> nodesrc = nodesrcif.data;

        /* Finding correct table */
        unsigned int g = 0;
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
                
                /* Compute indices for hash table */
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

                stream_index1.write(indexAdj);
                stream_index2.write(indexEdge);
                stream_ntable.write(g);
                stream_edge.write(vertexIndexing.concat(vertexIndexed));
                stream_end_out.write(false);
            }

#ifdef UNDIRECTED   
            if (tDescriptors[g].src_label == labeldst &&
                    tDescriptors[g].dst_label == labelsrc){

                ap_uint<V_ID_W> vertexIndexing, vertexIndexed; 
                if(!tDescriptors[g].dir){
                    vertexIndexing = nodesrc;
                    vertexIndexed = nodedst;
                } else {
                    vertexIndexing = nodedst;
                    vertexIndexed = nodesrc;
                }
                
                /* Compute indices for hash table */
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

                stream_index1.write(indexAdj);
                stream_index2.write(indexEdge);
                stream_ntable.write(g);
                stream_edge.write(vertexIndexing.concat(vertexIndexed));
                stream_end_out.write(false);
            }
#endif

        }
        last = nodesrcif.last;
    } while(!last);
    stream_end_out.write(true);
}

void countEdges(
        hls::stream<T_NODE> &stream_src,
        hls::stream<T_NODE> &stream_dst,
        hls::stream<T_LABEL> &stream_src_l,
        hls::stream<T_LABEL> &stream_dst_l,
        TableDescriptor *tDescriptors,
        AdjHT *hTables,
        ap_uint<DDR_W> *htb_buf,
        ap_uint<8> numTables)
{

#pragma HLS DATAFLOW

    hls::stream<ap_uint<2*V_ID_W>, S_D> stream_edge("Offset edge");
    hls::stream<ap_uint<H_W_1>, S_D> stream_c_index1("Counters index 1");
    hls::stream<ap_uint<H_W_2>, S_D> stream_c_index2("Counters index 2");
    hls::stream<ap_uint<8>, S_D> stream_c_ntable("Counters table");
    hls::stream<bool> stream_end_count("Counters end");

    /* Count edges per vertex source */
    edgeToHash(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            tDescriptors,
            numTables,
            stream_edge,
            stream_c_index1,
            stream_c_index2,
            stream_c_ntable,
            stream_end_count);

    /* Update specific counter */
    increaseCounter(
            stream_edge,
            stream_c_index1,
            stream_c_index2,
            stream_c_ntable,
            stream_end_count,
            hTables,
            htb_buf);

}

void writeEdges(
        hls::stream<T_NODE> &stream_src,
        hls::stream<T_NODE> &stream_dst,
        hls::stream<T_LABEL> &stream_src_l,
        hls::stream<T_LABEL> &stream_dst_l,
        TableDescriptor *tDescriptors,
        AdjHT *hTables,
        ap_uint<DDR_W> *htb_buf,
        ap_uint<8> numTables)
{

#pragma HLS DATAFLOW

    hls::stream<ap_uint<2*V_ID_W>, S_D> stream_edge("Offset edge");
    hls::stream<ap_uint<H_W_1>, S_D> stream_o_index1("Offsets index 1");
    hls::stream<ap_uint<H_W_2>, S_D> stream_o_index2("Offsets index 2");
    hls::stream<ap_uint<8>, S_D> stream_o_ntable("Offset table");
    hls::stream<bool> stream_end_offset("Offset end");

    /* Write edges based on precomputed offsets */
    edgeToHash(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            tDescriptors,
            numTables,
            stream_edge,
            stream_o_index1,
            stream_o_index2,
            stream_o_ntable,
            stream_end_offset);

    storeEdges(
            stream_edge,
            stream_o_index1,
            stream_o_index2,
            stream_o_ntable,
            stream_end_offset,
            hTables,
            htb_buf);
}

/* Reads two times the data graph and fills the data stuctures */
void fillTables(

#ifdef DEBUG_INTERFACE
        unsigned int &debif_endpreprocess,
#endif /* DEBUG_INTERFACE */
        
        hls::stream<T_NODE> &stream_src,
        hls::stream<T_NODE> &stream_dst,
        hls::stream<T_LABEL> &stream_src_l,
        hls::stream<T_LABEL> &stream_dst_l,
        AdjHT *hTables0,
        AdjHT *hTables1,
        ap_uint<DDR_W> *htb_buf,
        TableDescriptor *tDescriptors,
        ap_uint<8> numTables)
{

    /* Resetting portion of memory dedicated to counters 
     * 1 << H_W_1 * H_W_2 is the number of counters needed
     * for each table, then it should be divided by the number
     * of counters stored in each row which is 1 << (DDR_BIT - C_W)*/
    ap_uint<64> end_addr = numTables * HTB_SIZE;
/*
RESET_HASHTABLES_LOOP:
    for (ap_uint<64> addr = 0; addr < end_addr; addr++){
        htb_buf[addr] = 0;
    }
*/
    ap_uint<64> start_addr = 0;
STORE_HASHTABLES_POINTER_LOOP:
    for (unsigned int ntb = 0; ntb < numTables; ntb++){
        hTables0[ntb].start_offset = start_addr;
        hTables0[ntb].n_edges = 0;
        start_addr += HTB_SIZE;
    }

    countEdges(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            tDescriptors,
            hTables0,
            htb_buf,
            numTables);
    
    /* From counts to offsets */
    counterToOffset(
            numTables,
            hTables0,
            htb_buf);
    
    start_addr = end_addr;
STORE_EDGES_POINTER_LOOP:
    for (ap_uint<8> ntb = 0; ntb < numTables; ntb++){
        hTables0[ntb].start_edges = start_addr;
        start_addr += (hTables0[ntb].n_edges >> (DDR_BIT - E_W)) + 1;
#ifndef __SYNTHESIS__
        assert(start_addr < DDR_WIDTH);
#endif
        /* Fill the last row of tables's edges of 1s.
         * Because while reading the minimum set it is
         * read the entire line, and must be possible
         * to recognize which are real edges and which no */
        hTables1[ntb].start_offset = hTables0[ntb].start_offset;
        hTables1[ntb].start_edges = hTables0[ntb].start_edges;
        hTables1[ntb].n_edges = hTables0[ntb].n_edges;
        hTables1[ntb].hash_set = hTables0[ntb].hash_set;
        //htb_buf[start_addr - 1] = ~((ap_uint<DDR_W>)0);
    }

    writeEdges(
            stream_src,
            stream_dst,
            stream_src_l,
            stream_dst_l,
            tDescriptors,
            hTables0,
            htb_buf,
            numTables);

#ifdef DEBUG_INTERFACE
    debif_endpreprocess = hTables0[0].n_edges;
#endif /* DEBUG_INTERFACE */

#ifndef __SYNTHESIS__
    end_addr = start_addr;
    std::cout << "Occupied " << end_addr * (1 << (DDR_BIT-3)) << " bytes, " << end_addr <<
        " words." << std::endl;
#endif

/* std::ofstream f("CHECK.txt"); */
/* for(ap_uint<8> tab = 0; tab < numTables; tab++){ */
/* f << "Table " << tab << std::endl; */
/* ap_uint<64> start = hTables0[tab].start_offset; */
/* ap_uint<H_W_1> counter = 0; */
/* for(ap_uint<64> addr = 0; */
/* addr < HTB_SIZE; */
/* addr++){ */
/* ap_uint<DDR_W> row = htb_buf[start + addr]; */
/* for(int g = 0; g < CNT_ROW; g++){ */
/* if (((addr << (DDR_BIT - C_W)) + g) % (1UL << H_W_2) == 0){ */
/* f << counter << ": " << std::endl; */
/* counter++; */
/* } */

/* ap_uint<(1UL << C_W)> edge = row.range(((g+1)<<C_W)-1, g<<C_W); */
/* f << "\t" << edge.range((1UL << C_W)-2, 0) << " -> " */
/* << edge.test((1UL << C_W)-1) << std::endl; */

/* } */
/* } */
/* } */
/* f << std::endl; */
/* for(ap_uint<8> tab = 0; tab < numTables; tab++){ */
/* f << "Table " << tab << std::endl; */
/* ap_uint<64> start = hTables0[tab].start_edges; */
/* for(ap_uint<64> addr = 0; */
/* addr <= (hTables0[tab].n_edges >> (DDR_BIT-E_W)); */
/* addr++){ */
/* ap_uint<DDR_W> row = htb_buf[start + addr]; */
/* for(int g = 0; g < EDGE_ROW; g++){ */
/* ap_uint<(1UL << E_W)> edge = row.range(((g+1)<<E_W)-1, g<<E_W); */
/* f << "\t" << edge.range(2*V_ID_W-1, V_ID_W) << " -> " */
/* << edge.range(V_ID_W-1, 0) << std::endl; */
/* } */
/* } */
/* } */
/* f.close(); */
}

void mwj_extract_hashtovid(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<H_W_1>> &stream_inter_in,
        hls::stream<bool> &stream_end_in,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_set,
        AdjHT *hTables,
        ap_uint<DDR_W> *htb_buf,
        hls::stream<bool, 1> &stream_stop,

        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_out,
        hls::stream<bool> &stream_set_end_out,
        hls::stream<bool> &stream_end_out)
{

    hls::stream<ap_uint<V_ID_W>> stream_hash_in("Extract hashtovid hash_in");
    hls::stream<ap_uint<64>> stream_hash_out("Extract hashtovid hash_out");
    ap_uint<V_ID_W + 1 + 8> minData;
    ap_uint<H_W_1> hashinter;
    ap_uint<V_ID_W> vertex;
    ap_uint<V_ID_W*2> edge;
    bool stop, last;

    while(1){
        /* Bypassing solution for verify */
        if (stream_end_embed_in.read_nb(last)){

EXTRACT_COPYING_EMBEDDING_LOOP:
            while(!last){
                stream_embed_out.write(stream_embed_in.read());
                stream_end_embed_out.write(last);
                last = stream_end_embed_in.read();
            }
            stream_end_embed_out.write(last);

            minData = stream_min_set.read();

            /* filter to remove hash collisions when passing from H_W_1 to H_W_2 */
            ap_uint<16> filter[(1UL << (H_W_2 - 4))];
EXTRACT_RESET_FILTER_LOOP:
            for (int g = 0; g < (1UL << (H_W_2 - 4)); filter[g++] = 0);

            //std::cout << "\tHash to vid: " << std::endl;
            last = stream_end_in.read();
EXTRACT_HASHTOVID_LOOP:
            while(!last){
                stream_end_out.write(false);
                hashinter = stream_inter_in.read();
                ap_uint<(1UL << C_W)> start_off = 0;
                ap_uint<(1UL << C_W)> end_off = 0;
                ap_uint<H_W_1> index1, index1_f;
                ap_uint<H_W_2> index2, index2_f;

                if (minData.test(8)) {

                    /* Retriving hash for vertex indexing the table */
                    stream_hash_in.write(minData.range(V_ID_W + 8, 9));
                    xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
                    index1 = stream_hash_out.read().range(H_W_1 - 1, 0);
                    index2 = hashinter.range(H_W_2 - 1, 0);
                } else {
                    index1 = hashinter;
                    index2 = 0;
                }

                if (index2 != 0){ 
                    index1_f = index1;
                    index2_f = index2 - 1;
                } else if (index2 == 0 && index1 != 0){
                    index1_f = index1 - 1;
                    index2_f = (1UL << H_W_2) - 1;
                }

                if (!(index2 == 0 && index1 == 0)){
                    start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                            index1_f,
                            index2_f,
                            htb_buf,
                            hTables[minData.range(7,0)].start_offset);
                    start_off = start_off.range((1UL << C_W) - 2, 0);

#ifdef DEBUG_STATS
                    debug::hashtovid_reads++;
#endif
                }

                if (!minData.test(8)){
                    index2 = (1UL << H_W_2) - 1;
                }

                end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                        index1,
                        index2,
                        htb_buf,
                        hTables[minData.range(7,0)].start_offset);
                end_off = end_off.range((1UL << C_W) - 2, 0);

#ifdef DEBUG_STATS
                debug::hashtovid_reads++;
#endif

                /* If hash bit already present do no stream anything */
                if (minData.test(8)){
                    ap_uint<4> bitindex = index2.range(3, 0);
                    ap_uint<H_W_2 - 4> arrayindex = index2.range(H_W_2 - 1, 4);
                    ap_uint<16> filterword = filter[arrayindex];
                    if (filterword[bitindex] == 1){
                        end_off = start_off;
                    } else {
                        filterword[bitindex] = 1;
                        filter[arrayindex] = filterword;
                    }
                }

/* std::cout << "\t" << hashinter << ": " << */
/* start_off << " to " << end_off << */
/* "from tab " << minData.range(7,0) << */
/* "which starts at " << hTables[minData.range(7,0)].start_offset << */
/* std::endl; */
EXTRACT_HASHTOVID_READ_LOOP:
                for (; start_off < end_off; start_off++){

                    edge = read_table<32, 1, 0, E_W>(
                            start_off,
                            0,
                            htb_buf,
                            hTables[minData.range(7, 0)].start_edges);

                    if (minData.test(8)){
                        vertex = edge.range(V_ID_W - 1, 0);
                    } else {
                        vertex = edge.range(2*V_ID_W - 1, V_ID_W);
                    }

                    stream_inter_out.write(vertex);
                    stream_set_end_out.write(false);            

#ifdef DEBUG_STATS
                    debug::hashtovid_reads++;
#endif
                }

                stream_set_end_out.write(true);
                last = stream_end_in.read();
            }
            stream_end_out.write(true);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

void mwj_extract_bagtoset(
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_set_end_in,
        hls::stream<bool> &stream_end_in,
        hls::stream<bool, 1> &stream_stop,

        hls::stream<ap_uint<V_ID_W>> &stream_inter_out,
        hls::stream<bool> &stream_end_out)
{
   
    ap_uint<V_ID_W> set[MAX_CL];
    uint8_t counter;
    bool stop, last;

    while(1){
        if (stream_end_in.read_nb(last)){
            counter = 0;

EXTRACT_BAGTOSET_LOOP:
            while(!last){
                bool lastSet = stream_set_end_in.read();
                counter = 0;
EXTRACT_BAGTOSET_SET_LOOP:
                while(!lastSet){
                    ap_uint<V_ID_W> vertex = stream_inter_in.read();
                    uint8_t nSet = 0;
EXTRACT_BAGTOSET_SETCHECKER_LOOP:
                    for(; nSet < counter; nSet++){
                        if (vertex == set[nSet]){
                            break;
                        }
                    }
                    if (nSet == counter){
#ifndef __SYNTHESIS__
                        assert(counter < MAX_CL);
#endif
                        set[counter++] = vertex;
                        stream_inter_out.write(vertex);
                        stream_end_out.write(false);
                    }
                    lastSet = stream_set_end_in.read();    
                }
                last = stream_end_in.read();
            }
            stream_end_out.write(true);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

/* Multiway join propose: retrieves the tables in which
 * the current query vertex is involved, computes the sizes
 * and keep track of the smallest one */
void mwj_propose_findmin(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<DDR_W> *htb_buf,
        long unsigned cycles,

        hls::stream<ap_uint<32>> &stream_minread_out,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_out,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in("Propose findmin hash_in");
    hls::stream<ap_uint<64>> stream_hash_out("Propose readmin hash_in");
    ap_uint<8> tableIndex;
    ap_uint<V_ID_W> curQV;
    ap_uint<32> minSize;
    ap_uint<32> minStart, minEnd, minOff;
    ap_uint<8 + 1 + V_ID_W> minData;
    ap_uint<V_ID_W> curEmb[MAX_QV];

    for(long unsigned c = 0; c < cycles; c++){
        curQV = stream_embed_in.read();
        tableIndex = 0;
        minSize = (1UL << 32) - 1;

        /* std::cout << "Current solution (" << curQV << "): "; */
PROPOSE_COPYING_EMBEDDING_LOOP:
        for (int g = 0; g < curQV; g++){
#pragma HLS pipeline II=1
            curEmb[g] = stream_embed_in.read();
            stream_embed_out.write(curEmb[g]);
            stream_end_embed_out.write(false);
        }
        stream_end_embed_out.write(true);

        //hls::print("Proposing for %d\n", (int)curQV);
        /* Find sizes of sets indexed by the current query vertex */
PROPOSE_TBINDEXING_LOOP:
        for(int g = 0; g < qVertices[curQV].numTablesIndexing; g++){
#pragma HLS pipeline II=1
            tableIndex = qVertices[curQV].tables_indexing[g];

            //std::cout << "\tEvaluating " << tableIndex << " : size " <<
            //   hTables[tableIndex].hash_set << std::endl; 
            if (hTables[tableIndex].hash_set < minSize){
                minSize = hTables[tableIndex].hash_set;
                minOff = hTables[tableIndex].start_edges;
                minStart = 0;
                minEnd = hTables[tableIndex].n_edges;
                minData.range(7, 0) = tableIndex;
                minData.clear(8);
                minData.range(V_ID_W + 8, 9) = 0;
            }
        }

        /* Find sizes of sets in which the current query vertex
         * is indexed by an other query vertex */
PROPOSE_TBINDEXED_LOOP:
        for(int g = 0; g < qVertices[curQV].numTablesIndexed; g++){
#pragma HLS pipeline II=1
            ap_uint<(1UL << C_W)> start_off = 0;
            ap_uint<(1UL << C_W)> end_off;
            tableIndex = qVertices[curQV].tables_indexed[g];
            uint8_t ivPos = qVertices[curQV].vertex_indexing[g];

            stream_hash_in.write(curEmb[ivPos]);
            xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
            ap_uint<H_W_1> index = stream_hash_out.read().range(H_W_1 - 1, 0);

            if (index != 0){ 
                start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                        index - 1,
                        (1UL << H_W_2) - 1,
                        htb_buf,
                        hTables[tableIndex].start_offset);
                start_off = start_off.range((1UL << C_W) - 2, 0);
            }

            end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                    index,
                    (1UL << H_W_2) - 1,
                    htb_buf,
                    hTables[tableIndex].start_offset);
            end_off = end_off.range((1UL << C_W) - 2, 0);

#ifdef DEBUG_STATS
            debug::findmin_reads += 2;
#endif
            // std::cout << "\tEvaluating " << tableIndex << " : size " <<
            //    end_off << " - " << start_off << std::endl; 
            if ((end_off - start_off) < minSize) {
                minSize = end_off - start_off;
                minOff = hTables[tableIndex].start_edges;
                minStart = start_off;
                minEnd = end_off;
                minData.range(7, 0) = tableIndex;
                minData.set(8);
                minData.range(V_ID_W + 8, 9) = curEmb[ivPos];
            }
        }
        stream_minread_out.write(minStart);
        stream_minread_out.write(minEnd);
        stream_minread_out.write(minOff);
        stream_min_out.write(minData);
    }
}

/* Read the minimum set from memory */
void mwj_propose_readmin(
        hls::stream<ap_uint<32>> &stream_minread_in,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_in,
        ap_uint<DDR_W> *htb_buf,
        /* hls::burst_maxi<row_t> htb_buf, */
        long unsigned cycles,

        hls::stream<ap_uint<V_ID_W>> &stream_min_out,
        hls::stream<bool> &stream_end_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in("Read min hash_in");
    hls::stream<ap_uint<64>> stream_hash_out("Read min hash_out");
    ap_uint<32> minStart; 
    ap_uint<32> minEnd;
    ap_uint<32> minOff;
    ap_uint<V_ID_W + 1 + 8> minData;
    ap_uint<V_ID_W*2> edge;
    ap_uint<V_ID_W> vertex, vertexCheck;
    ap_uint<V_ID_W> set[MAX_CL];
    
    for(long unsigned c = 0; c < cycles; c++){

        minStart = stream_minread_in.read();
        minEnd = stream_minread_in.read();
        minOff = stream_minread_in.read();
        minData = stream_min_in.read();
        unsigned int rowstart = minOff + (minStart >> (DDR_BIT - E_W));
        unsigned int rowend = minOff + (minEnd >> (DDR_BIT - E_W));
        unsigned int window_left = minStart.range((DDR_BIT - E_W) - 1, 0);
        unsigned int window_right = minEnd.range((DDR_BIT - E_W) - 1, 0) + 
            (rowend - rowstart) * EDGE_ROW;
        unsigned int cnt = 0;

        /* std::cout << "\tMinimum set is " << minData.range(7, 0) << std::endl; */
        /* std::cout << "\tReading from " << rowstart << " to " << rowend */
        /* << " window " << window_left << "-" << window_right << std::endl; */
        /* htb_buf.read_request(rowstart, rowend-rowstart+1); */

        if (minData.test(8)){
PROPOSE_READ_MIN_INDEXED_LOOP:
            for (int g = rowstart; g <= rowend; g++){
                row_t row = htb_buf[g];
                /* row_t row = htb_buf.read(); */
                for (int i = 0; i < EDGE_ROW; i++, cnt++){
                    if (cnt >= window_left && cnt < window_right){
                        /* edge = row.range(((i + 1) * V_ID_W*2) - 1, i * V_ID_W*2); */
                        edge = row.range((1UL << E_W) - 1, 0);
                        vertexCheck = edge.range(V_ID_W*2 - 1, V_ID_W);
                        vertex = edge.range(V_ID_W - 1, 0);
                        if (minData.range(V_ID_W + 8, 9) == vertexCheck){
                            stream_min_out.write(vertex);
                            stream_end_out.write(false);
                        }
                    }
                    row >>= (1UL << E_W);
                }
#ifdef DEBUG_STATS
                debug::readmin_reads++;
#endif
            }
#ifdef DEBUG_STATS
            debug::indexed_tables++;
#endif
        } else {
            ap_uint<H_W_1> hash_buff = 0;
            ap_uint<5> set_counter = 0;
            bool flag_buff = false;
            bool flag_new = true;

PROPOSE_READ_MIN_INDEXING_LOOP:
            for (int g = rowstart; g <= rowend; g++){
                row_t row = htb_buf[g];
                /* row_t row = htb_buf.read(); */
                for (int i = 0; i < EDGE_ROW; i++, cnt++){
                    if (cnt < window_right){
                        edge = row.range((1UL << E_W) - 1, 0);
                        /* edge = row.range(((i + 1) * V_ID_W*2) - 1, i * V_ID_W*2); */
                        vertex = edge.range(V_ID_W*2 - 1, V_ID_W);
                        stream_hash_in.write(vertex);
                        xf::database::hashLookup3<V_ID_W>(
                                stream_hash_in,
                                stream_hash_out);
                        ap_uint<H_W_1> vertexHash = stream_hash_out.read();

                        if (flag_buff && hash_buff == vertexHash){
                            int nSet = 0;
                            flag_new = true;
EXTRACT_BAGTOSET_SETCHECKER_LOOP:
                            for(; nSet < set_counter; nSet++){
                                if (vertex == set[nSet]){
                                    flag_new = false;
                                    break;
                                }
                            }
                        } else {
                            flag_new = true;
                            set_counter = 0;
                        }

#ifndef __SYNTHESIS__
                        assert(set_counter < MAX_CL);
#endif
                        if (flag_new) {
                            set[set_counter++] = vertex;
                            stream_min_out.write(vertex);
                            stream_end_out.write(false);
                        }
                        hash_buff = vertexHash;
                        flag_buff = true;
                    }
                    row >>= (1UL << E_W);
                }
#ifdef DEBUG_STATS
                debug::readmin_reads++;
#endif
            }
#ifdef DEBUG_STATS
            debug::indexing_tables++;
#endif
        }
        stream_end_out.write(true);
    }
}

void mwj_intersect(
        AdjHT *hTables,
        QueryVertex *qVertices,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<V_ID_W>> &stream_min_in,
        hls::stream<bool> &stream_end_min_in,
        ap_uint<DDR_W> *htb_buf,
        long unsigned cycles,

        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_out,
        hls::stream<bool> &stream_end_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in("Intersect hash_in");
    hls::stream<ap_uint<64>> stream_hash_out("Intersect hash_out");
    ap_uint<64> candidate_hash;
    ap_uint<V_ID_W> candidate_v;
    ap_uint<8> tableIndex, curQV;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool last;

    for(long unsigned c = 0; c < cycles; c++){
        last = stream_end_embed_in.read();
        curQV = 0;

INTERSECT_COPYING_EMBEDDING_LOOP:
        while(!last){
            curEmb[curQV] = stream_embed_in.read();
            stream_embed_out.write(curEmb[curQV]);
            stream_end_embed_out.write(false);
            curQV++;
            last = stream_end_embed_in.read();
        }
        stream_end_embed_out.write(true);

        /* bypass data for extract phase */
        /* stream_min_desc_out.write(stream_min_desc_in.read()); */

        /* std::cout << "Intersection: probing" << std::endl; */
        last = stream_end_min_in.read();
INTERSECT_LOOP:
        while(!last){
            candidate_v = stream_min_in.read();
            stream_hash_in.write(candidate_v);
            xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
            candidate_hash = stream_hash_out.read();

            //std::cout << "\t" << candidate << std::endl;
            bool inter = true;
#ifdef DEBUG_STATS
            bool real_inter = true;
#endif
INTERSECT_TBINDEXED_LOOP:
            for(int g = 0; 
                    g < qVertices[curQV].numTablesIndexed;
                    /* g < qVertices[curQV].numTablesIndexed && inter; */
                    g++)
            {
                ap_uint<(1UL << C_W)> offset;
                tableIndex = qVertices[curQV].tables_indexed[g];
                uint8_t ivPos = qVertices[curQV].vertex_indexing[g];

                stream_hash_in.write(curEmb[ivPos]);
                xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
                ap_uint<H_W_1> index1 = stream_hash_out.read().range(H_W_1 - 1, 0);
                ap_uint<H_W_2> index2 = candidate_hash.range(H_W_2 - 1, 0);

                offset = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                        index1,
                        index2,
                        htb_buf,
                        hTables[tableIndex].start_offset);
                inter = inter && offset.test((1UL << C_W) - 1);
#ifdef DEBUG_STATS
                debug::intersect_reads++;

                /* Computing when there is an alias */
                ap_uint<H_W_1> index1_f;
                ap_uint<H_W_2> index2_f;
                ap_uint<(1UL << C_W)> start_off;

                bool checked = false;
                if (index2 != 0){
                    index1_f = index1;
                    index2_f = index2 - 1;
                } else if (index2 == 0 && index1 != 0){
                    index1_f = index1 - 1;
                    index2_f = (1UL << H_W_2)-1;
                }

                if (!(index2 == 0 && index1 == 0)){
                    start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                            index1_f,
                            index2_f,
                            htb_buf,
                            hTables[tableIndex].start_offset);
                    start_off = start_off.range((1UL << C_W) - 2, 0);
                }

                for (; start_off < offset.range((1UL << C_W) - 2,0) && checked == false; start_off++){
                    ap_uint<2 * V_ID_W> edge = read_table<32, 1, 0, E_W>(
                            start_off,
                            0,
                            htb_buf,
                            hTables[tableIndex].start_edges);

                    ap_uint<V_ID_W> vertexIndexed, vertexIndexing;
                    vertexIndexed = edge.range(V_ID_W - 1, 0);
                    vertexIndexing = edge.range(2*V_ID_W - 1, V_ID_W);
                    if (vertexIndexing == curEmb[ivPos] &&
                            vertexIndexed == candidate_v){
                        checked = true;
                    }
                }

                if (offset.test((1UL << C_W) - 1) && checked) {
                    /* True positive */
                    debug::intersect_bit_truepositive++;
                } else if (offset.test((1UL << C_W) - 1) && !checked) {
                    /* False positive */
                    real_inter = false;
                    debug::intersect_bit_falsepositive++;
                } else if (!offset.test((1UL << C_W) - 1) && !checked) {
                    /* True negative */
                    real_inter = false;
                    debug::intersect_bit_truenegative++;
                } else if (!offset.test((1UL << C_W) - 1) && checked) {
                    /* False negative */
                    debug::intersect_bit_falsenegative++;
                }
#endif
            }

#if INTERESCT_INDEXING_LOOP
INTERSECT_TBINDEXING_LOOP:
            for(int g = 0;
                    g < qVertices[curQV].numTablesIndexing;
                    /* g < qVertices[curQV].numTablesIndexing && inter; */
                    g++)
            {
                ap_uint<(1UL << C_W)> start_off, end_off;
                tableIndex = qVertices[curQV].tables_indexing[g];
                ap_uint<H_W_1> candidate = candidate_hash.range(H_W_1 - 1, 0);
                if (candidate != 0){
                    start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                            candidate - 1,
                            (1UL << H_W_2) - 1,
                            htb_buf,
                            hTables[tableIndex].start_offset);
                    start_off = start_off.range((1UL << C_W) - 2, 0);
#ifdef DEBUG_STATS
                    debug::intersect_reads++;
#endif
                } else {
                    start_off = 0;
                }

                end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                        candidate,
                        (1UL << H_W_2) - 1,
                        htb_buf,
                        hTables[tableIndex].start_offset);
                end_off = end_off.range((1UL << C_W) - 2, 0);
                inter = inter && (end_off > start_off);

#ifdef DEBUG_STATS
                bool temp = (end_off > start_off);
                debug::intersect_reads++;

                // Computing when there is an alias
                bool checked = false;
                for (; start_off < end_off && checked == false; start_off++){
                    ap_uint<2 * V_ID_W> edge = read_table<32, 1, 0, E_W>(
                            start_off,
                            0,
                            htb_buf,
                            hTables[tableIndex].start_edges);

                    ap_uint<V_ID_W> vertexIndexing;
                    vertexIndexing = edge.range(2*V_ID_W - 1, V_ID_W);
                    if (vertexIndexing == candidate_v){
                        checked = true;
                    }
                }

                if (temp && checked) {
                    /* //True positive */
                    debug::intersect_bit_truepositive++;
                } else if (temp && !checked) {
                    /* //False positive */
                    real_inter = false;
                    debug::intersect_bit_falsepositive++;
                } else if (!temp && !checked) {
                    /* //True negative */
                    real_inter = false;
                    debug::intersect_bit_truenegative++;
                } else if (!temp && checked) {
                    /* //False negative */
                    debug::intersect_bit_falsenegative++;
                }

#endif
            }
#endif /* INTERSECT_INDEXING_LOOP */
#ifdef DEBUG_STATS
            if (inter && !real_inter) {
                /* //False positive */
                debug::intersect_sol_falsepositive++;
            } else if (inter && real_inter) {
                /* //True positive */
                debug::intersect_sol_truepositive++;
            } else if (!inter && !real_inter) {
                /* //True negative */
                debug::intersect_sol_truenegative++;
            } else if (!inter && real_inter) {
                /* //False negative */
                debug::intersect_sol_falsenegative++;
            }
#endif

            if (inter){
                stream_inter_out.write(candidate_v);
                stream_end_out.write(false);
            }

            last = stream_end_min_in.read();
        }
        stream_end_out.write(true);
    }
}

void mwj_verify_nothomomorphism(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        long unsigned cycles,
        
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<V_ID_W>> &stream_checked_out,
        hls::stream<bool> &stream_checked_end_out)
{
    ap_uint<8> curQV;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool last;

    for(long unsigned c = 0; c < cycles; c++){
        last = stream_end_embed_in.read();
        curQV = 0;

VERIFY_HOMOMO_COPYING_EMBEDDING_LOOP:
        while(!last){
            curEmb[curQV] = stream_embed_in.read();
            stream_embed_out.write(curEmb[curQV]);
            stream_end_embed_out.write(false);
            curQV++;
            last = stream_end_embed_in.read();
        }
        stream_end_embed_out.write(true);


        last = stream_end_inter_in.read();
VERIFY_CHECK_LOOP:
        while(!last){
            ap_uint<V_ID_W> vToVerify = stream_inter_in.read();
            bool homomorphism = false;
            for (int g = 0; g < curQV; g++){
                if (vToVerify == curEmb[g])
                    homomorphism = true;
            }
            if (!homomorphism){
                stream_checked_out.write(vToVerify);
                stream_checked_end_out.write(false);
            }

            last = stream_end_inter_in.read();
        }
        stream_checked_end_out.write(true);
    }
}

void mwj_verify_edge(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<DDR_W> *htb_buf,
        long unsigned cycles,
        
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<V_ID_W>> &stream_checked_out,
        hls::stream<bool> &stream_checked_end_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in("Verify edge hash_in");
    hls::stream<ap_uint<64>> stream_hash_out("Verify edge hash_out");
    ap_uint<8> curQV;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    ap_uint<DDR_W> cache;
    ap_uint<32> addr_prev_row = 0;
    bool last;

#ifdef DEBUG_STATS    
    std::ofstream of("../../../../verstats.txt", std::ofstream::app);
#endif

    for(long unsigned c = 0; c < cycles; c++){
        last = stream_end_embed_in.read();
        curQV = 0;

VERIFY_EDGE_COPYING_EMBEDDING_LOOP:
        while(!last){
            curEmb[curQV] = stream_embed_in.read();
            stream_embed_out.write(curEmb[curQV]);
            stream_end_embed_out.write(false);
            curQV++;
            last = stream_end_embed_in.read();
        }
        stream_end_embed_out.write(true);

        last = stream_end_inter_in.read();
VERIFY_CHECK_LOOP:
        while(!last){
            ap_uint<V_ID_W> vToVerify = stream_inter_in.read();
            bool checked = true;
            int g = 0;

VERIFY_CHECK_EDGE_LOOP:
            for(; g < qVertices[curQV].numTablesIndexed && checked; g++){
                checked = false;
                ap_uint<H_W_1> index1, index1_f;
                ap_uint<H_W_2> index2, index2_f;
                uint8_t tableIndex = qVertices[curQV].tables_indexed[g];
                uint8_t ivPos = qVertices[curQV].vertex_indexing[g];
                stream_hash_in.write(curEmb[ivPos]);
                stream_hash_in.write(vToVerify);
                xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
                index1 = stream_hash_out.read().range(H_W_1 - 1, 0);
                xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
                index2 = stream_hash_out.read().range(H_W_2 - 1, 0);

                ap_uint<(1UL << C_W)> start_off = 0;
                ap_uint<(1UL << C_W)> end_off = 0;

                if (index2 != 0){
                    index1_f = index1;
                    index2_f = index2 - 1;
                } else if (index2 == 0 && index1 != 0){
                    index1_f = index1 - 1;
                    index2_f = (1UL << H_W_2)-1;
                }

                if (!(index2 == 0 && index1 == 0)){
                    start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                            index1_f,
                            index2_f,
                            htb_buf,
                            hTables[tableIndex].start_offset);
                    start_off = start_off.range((1UL << C_W) - 2, 0);

#ifdef DEBUG_STATS
                    debug::verify_reads++;
#endif
                }

                end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                        index1,
                        index2,
                        htb_buf,
                        hTables[tableIndex].start_offset);
                end_off = end_off.range((1UL << C_W) - 2, 0);

#ifdef DEBUG_STATS
                debug::verify_reads++;
                /* of << "(" << vToVerify << ") " << start_off + hTables[tableIndex].start_edges << */
                /* " - " << end_off - start_off << std::endl; */
#endif
VERIFY_READ_MEMORY_LOOP:
                for (; start_off < end_off; start_off++){
#pragma HLS pipeline II=1
                    ap_uint<DDR_W> ram_row;
                    ap_uint<32> addr_row;
                    ap_uint<32> addr_inrow;

                    /* Compute address of row storing the counter */
                    addr_row = hTables[tableIndex].start_edges + (start_off >> (DDR_BIT - E_W));

                    /* Compute address of data inside the row */
                    addr_inrow = start_off.range((DDR_BIT - E_W) - 1, 0);

                    /* Read the data */
                    if (addr_prev_row == addr_row){
                        ram_row = cache;
#ifdef DEBUG_STATS
                        debug::verify_reusage++;
#endif
                    } else {
                        ram_row = htb_buf[addr_row];
                        cache = ram_row;
                        addr_prev_row = addr_row;
                    }

                    ap_uint<2 * V_ID_W> edge = ram_row >> (addr_inrow << E_W);
                    /* ap_uint<2 * V_ID_W> edge = read_table<32, 1, 0, E_W>( */
                    /* start_off, */
                    /* 0, */
                    /* htb_buf, */
                    /* hTables[tableIndex].start_edges); */

                    ap_uint<V_ID_W> vertexIndexed, vertexIndexing;
                    vertexIndexed = edge.range(V_ID_W - 1, 0);
                    vertexIndexing = edge.range(2*V_ID_W - 1, V_ID_W);
                    if (vertexIndexing == curEmb[ivPos] &&
                            vertexIndexed == vToVerify){
                        checked = true;
                    }
#ifdef DEBUG_STATS
                    debug::verify_reads++;
#endif
                }
            }

            /* #ifdef DEBUG_STATS */
            /* of << checked << std::endl << std::endl; */
            /* #endif */

            if (checked){
                stream_checked_out.write(vToVerify);
                stream_checked_end_out.write(false);
#ifdef DEBUG_STATS
                /* if (g > 0) */
                debug::solution_correct++;
            } else { 
                debug::solution_wrong++; 
#endif
            }

            last = stream_end_inter_in.read();
        }
        stream_checked_end_out.write(true);
    }
}

void mwj_verify_add(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        long unsigned cycles,
        long unsigned &nPartSol,
        unsigned short nQueryVer,
        hls::stream<bool> streams_stop[2],

        hls::stream<ap_uint<V_ID_W>> &stream_partial_out,
#ifdef COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif
        )
{
    ap_uint<8> curQV;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool last;
    T_NODE node;

#ifdef COUNT_ONLY
    unsigned long int counter {0};
#endif
    
    /* start the pipeline with empty solution */
    for(long unsigned c = 0; c < cycles; c++){
        last = stream_end_embed_in.read();
        curQV = 0;

VERIFY_ADD_COPYING_EMBEDDING_LOOP:
        while(!last){
            curEmb[curQV] = stream_embed_in.read();
            curQV++;
            last = stream_end_embed_in.read();
        }

        last = stream_end_inter_in.read();
VERIFY_CHECK_LOOP:
        while(!last){
            ap_uint<V_ID_W> vToVerify = stream_inter_in.read();

            /* Write in the correct stream */
            if (curQV == nQueryVer - 1){
#ifdef COUNT_ONLY
                counter++;
#else
VERIFY_WRITE_FINAL_LOOP:
                for (int g = 0; g < curQV; g++){
                    node.data = curEmb[g];
                    node.last = false;
                    node.keep = ~0;
                    result.write(node);
                }
                node.data = vToVerify;
                node.last = false;
                node.keep = ~0;
                result.write(node);
#endif
            } else {
VERIFY_WRITE_PARTIAL_LOOP:
                stream_partial_out.write((ap_uint<V_ID_W>)curQV + 1);
                for (int g = 0; g < curQV; g++){
                    stream_partial_out.write(curEmb[g]);
                }
                stream_partial_out.write(vToVerify);
                nPartSol++;
            }
#ifdef DEBUG_STATS
            debug::embeddings++;
#endif
            last = stream_end_inter_in.read();
        }
        nPartSol--;
    }
        
    if (nPartSol == 0){
#ifdef COUNT_ONLY
        /* Write in output number of results */
        result += counter;
#else
        /* Write last node */
        node.data = 0;
        node.last = true;
        node.keep = ~0;
        result.write(node);
#endif
        for (int g = 0; g < STOP_S; g++)
            streams_stop[g].write(true);
    }
}

void starter(
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<DDR_W> *htb_buf,

        hls::stream<long unsigned> &nset,
        hls::stream<bool> &endset,
        hls::stream<ap_uint<V_ID_W>> &start_nodes)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in("Starter hash_in");
    hls::stream<ap_uint<64>> stream_hash_out("Starter hash_out");
    ap_uint<8> tableIndex = 0;
    ap_uint<32> minSize = (1UL << 32) - 1;
    ap_uint<H_W_1> hash_buff = 0;
    ap_uint<V_ID_W> set[MAX_CL];
    ap_uint<V_ID_W*2> edge;
    ap_uint<V_ID_W> vertex;
    ap_uint<5> set_counter = 0;
    bool flag_buff = false;
    bool flag_new = true;
    unsigned int rowstart; 
    unsigned int rowend;
    unsigned int window_right;
    unsigned int cnt = 0;
    long unsigned nPartSol = 0;

PROPOSE_TBINDEXING_LOOP:
    for(int g = 0; g < qVertices[0].numTablesIndexing; g++){
#pragma HLS pipeline II=1
        tableIndex = qVertices[0].tables_indexing[g];

        if (hTables[tableIndex].hash_set < minSize){
            minSize = hTables[tableIndex].hash_set;
            rowstart = hTables[tableIndex].start_edges;
            rowend = hTables[tableIndex].start_edges + 
                (hTables[tableIndex].n_edges >> (DDR_BIT - E_W));
        }
    }
    window_right = hTables[tableIndex].n_edges.range((DDR_BIT - E_W) - 1, 0) + 
        (rowend - rowstart) * EDGE_ROW;

PROPOSE_READ_MIN_INDEXING_LOOP:
    for (int g = rowstart; g <= rowend; g++){
        row_t row = htb_buf[g];
        /* row_t row = htb_buf.read(); */
        for (int i = 0; i < EDGE_ROW; i++, cnt++){
            if (cnt < window_right){
                edge = row.range((1UL << E_W) - 1, 0);
                /* edge = row.range(((i + 1) * V_ID_W*2) - 1, i * V_ID_W*2); */
                vertex = edge.range(V_ID_W*2 - 1, V_ID_W);
                stream_hash_in.write(vertex);
                xf::database::hashLookup3<V_ID_W>(
                        stream_hash_in,
                        stream_hash_out);
                ap_uint<H_W_1> vertexHash = stream_hash_out.read();

                if (flag_buff && hash_buff == vertexHash){
                    int nSet = 0;
                    flag_new = true;
EXTRACT_BAGTOSET_SETCHECKER_LOOP:
                    for(; nSet < set_counter; nSet++){
                        if (vertex == set[nSet]){
                            flag_new = false;
                            break;
                        }
                    }
                } else {
                    flag_new = true;
                    set_counter = 0;
                }

#ifndef __SYNTHESIS__
                assert(set_counter < MAX_CL);
#endif
                if (flag_new) {
                    set[set_counter++] = vertex;
                    nPartSol++;
                    start_nodes.write(1);
                    start_nodes.write(vertex);
                }
                 
                if (nPartSol == 1000){
                    nset.write(1000);
                    endset.write(false);
                    nPartSol = 0;
                }
                hash_buff = vertexHash;
                flag_buff = true;
            }
            row >>= (1UL << E_W);
        }
    }
    nset.write(nPartSol);
    endset.write(true);
}

void mwj_2to1(
        hls::stream<ap_uint<V_ID_W>> &from_loop,
        hls::stream<ap_uint<V_ID_W>> &from_ext,
        hls::stream<ap_uint<V_ID_W>> &output,
        long unsigned cycles,
        bool &sel) 
{
    ap_uint<V_ID_W> curQV;
    for(long unsigned c = 0; c < cycles; c++){
        if (sel){
            curQV = from_ext.read();
            output.write(curQV);
            for (int g = 0; g < curQV; g++){
                output.write(from_ext.read());
            }
        } else {
            curQV = from_loop.read();
            output.write(curQV);
            for (int g = 0; g < curQV; g++){
                output.write(from_loop.read());
            }
        }
    }
    sel = false;
}

void multiwayJoin(
        ap_uint<DDR_W> *htb_buf0,
        ap_uint<DDR_W> *htb_buf1,
/* hls::burst_maxi<row_t> htb_buf1, */
        ap_uint<DDR_W> *htb_buf2,
        ap_uint<DDR_W> *htb_buf3,
        row_t *res_buf,
        AdjHT *hTables0,
        AdjHT *hTables1,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        unsigned short nQueryVer,
        hls::stream<ap_uint<V_ID_W>> &start_sol,
        hls::stream<unsigned long> &nstartset,
#ifdef COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif
        )
{

#pragma HLS STABLE variable=htb_buf0
#pragma HLS STABLE variable=htb_buf1
#pragma HLS STABLE variable=htb_buf2
#pragma HLS STABLE variable=htb_buf3
#pragma HLS STABLE variable=hTables0
#pragma HLS STABLE variable=hTables1
#pragma HLS STABLE variable=qVertices0
#pragma HLS STABLE variable=qVertices1
#pragma HLS STABLE variable=nQueryVer
#pragma HLS dataflow 
   
    /* dynfifo signals */ 
    hls::stream<ap_uint<V_ID_W>, 32> dyn_in_stream("in stream");
    hls::stream<ap_uint<V_ID_W>, 32> dyn_out_stream("out stream");
    hls::stream<bool, 1> streams_stop[STOP_S]; 

    hls::stream<ap_uint<V_ID_W>, S_D> p_in_stream("in stream propose");
    
    /* Propose data out */
    hls::stream<ap_uint<V_ID_W>, MAX_QV> p_stream_embed
        ("Partial result propose");
    hls::stream<bool, MAX_QV> p_stream_embed_end("Partial result del. propose");
    
    hls::stream<ap_uint<V_ID_W>, S_D> p_stream_min
        ("Min. set vertices propose");
    hls::stream<bool, S_D> p_stream_min_end("Min. set verices del. propose");
    
    hls::stream<ap_uint<V_ID_W + 1 + 8>, S_D> p_stream_min_desc
        ("Min. set description propose");
    
    /* Propose internal streams */
    hls::stream<ap_uint<32>, S_D> p_find_stream_minread("Min read");
    hls::stream<ap_uint<V_ID_W + 1 + 8>, S_D> p_find_stream_min("Min desc");

    /* Intersect data out */    
    hls::stream<ap_uint<V_ID_W>, S_D> i_stream_hash_set("Set intersection hashes");
    hls::stream<bool, S_D> i_stream_hash_set_end("Set intersection hashes del.");
    hls::stream<ap_uint<V_ID_W>, MAX_QV> i_stream_embed
        ("Partial result intersect.");
    hls::stream<bool, MAX_QV> i_stream_embed_end
        ("Partial result del. intersect.");
    hls::stream<ap_uint<V_ID_W + 1 + 8>, S_D> i_stream_min_desc
        ("Min. set description propose");

    /* Verify internal streams */
    hls::stream<ap_uint<V_ID_W>, S_D> stream_checked1("Stream checked1");
    hls::stream<bool, S_D> stream_checked1_end;
    hls::stream<ap_uint<V_ID_W>, MAX_QV> stream_embed1("Stream embed1");
    hls::stream<bool, MAX_QV> stream_embed1_end;

    hls::stream<ap_uint<V_ID_W>, S_D> stream_checked2("Stream checked2");
    hls::stream<bool, S_D> stream_checked2_end;
    hls::stream<ap_uint<V_ID_W>, MAX_QV> stream_embed2("Stream embed2");
    hls::stream<bool, MAX_QV> stream_embed2_end;
    
    bool select = true;
    unsigned long nSolPart = nstartset.read();
    /* unsigned long nSolPart = 1; */
    
    dynfifo_init<
        ap_uint<V_ID_W>,         /* fifo data type */
        row_t,                   /* fifo data type */
        S_DEPTH,                 /* in/out stream size */
        BURST_S*2,               /* load/store stream size */
        DDR_WORD,                /* bitwidth ddr word */
        BURST_S,                 /* burst transaction size */
        RESULTS_SPACE>           /* memory words available */
            (dyn_in_stream,
             dyn_out_stream,
             streams_stop,
             res_buf);

    do {
/* std::cout << nSolPart << std::endl; */
        mwj_2to1(
                dyn_out_stream,
                start_sol,
                p_in_stream,
                nSolPart,
                select);

        mwj_propose_findmin(
                p_in_stream,
                hTables0,
                qVertices0,
                htb_buf0,
                nSolPart,
                p_find_stream_minread,
                p_find_stream_min,
                p_stream_embed,
                p_stream_embed_end);

        mwj_propose_readmin(
                p_find_stream_minread,
                p_find_stream_min,
                htb_buf1,
                nSolPart,
                p_stream_min,
                p_stream_min_end);

        mwj_intersect(
                hTables1,
                qVertices1,
                p_stream_embed,
                p_stream_embed_end,
                p_stream_min,
                p_stream_min_end,
                htb_buf2,
                nSolPart,
                i_stream_embed,
                i_stream_embed_end,
                i_stream_hash_set,
                i_stream_hash_set_end);

        mwj_verify_nothomomorphism(
                i_stream_embed,
                i_stream_embed_end,
                i_stream_hash_set,
                i_stream_hash_set_end,
                nSolPart,
                stream_embed1,
                stream_embed1_end,
                stream_checked1,
                stream_checked1_end);

        mwj_verify_edge(
                stream_embed1,
                stream_embed1_end,
                stream_checked1,
                stream_checked1_end,
                hTables1,
                qVertices0,
                htb_buf3,
                nSolPart,
                stream_embed2,
                stream_embed2_end,
                stream_checked2,
                stream_checked2_end);

        mwj_verify_add(
                stream_embed2,
                stream_embed2_end,
                stream_checked2,
                stream_checked2_end,
                nSolPart,
                nSolPart,
                nQueryVer,
                streams_stop,
                dyn_in_stream,
                result);
    
    } while(nSolPart);
    
}

void multiwayJoinWrap(
        ap_uint<DDR_W> *htb_buf0,
        ap_uint<DDR_W> *htb_buf1,
/* hls::burst_maxi<row_t> htb_buf1, */
        ap_uint<DDR_W> *htb_buf2,
        ap_uint<DDR_W> *htb_buf3,
        row_t *res_buf,
        AdjHT *hTables0,
        AdjHT *hTables1,
        QueryVertex *qVertices0,
        QueryVertex *qVertices1,
        unsigned short nQueryVer,

#ifdef COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif
        )
{

/* #pragma HLS dataflow */
    hls::stream<ap_uint<V_ID_W>, S_D> start_stream("start nodes");
    hls::stream<unsigned long, S_D> nstartsol("start nodes");
    hls::stream<bool, 1> endsol("start nodes");

    starter(hTables1,
            qVertices0,
            htb_buf0,
            nstartsol,
            endsol,
            start_stream);

    bool last;
    do { 
        multiwayJoin(
                htb_buf0,
                htb_buf1,
                htb_buf2,
                htb_buf3,
                res_buf,
                hTables0,
                hTables1,
                qVertices0,
                qVertices1,
                nQueryVer,
                start_stream,
                nstartsol,
                result);

        last = endsol.read();

    } while (!last);      
}

void subgraphIsomorphism(
        edge_t edge_buf[GRAPHS_SPACE],
        row_t htb_buf0[HASHTABLES_SPACE],
        row_t htb_buf2[HASHTABLES_SPACE],
        row_t htb_buf3[HASHTABLES_SPACE],
        row_t htb_buf1[HASHTABLES_SPACE],
/* hls::burst_maxi<row_t> htb_buf1, */
        row_t res_buf[RESULTS_SPACE],
        unsigned short numQueryVert,
        unsigned short numQueryEdges,
        unsigned long numDataEdges,

#ifdef DEBUG_INTERFACE
        unsigned int &debif_endpreprocess,
#endif /* DEBUG_INTERFACE */

#ifdef COUNT_ONLY
        long unsigned int &result
#else
        hls::stream<T_NODE> &result
#endif /* COUNT_ONLY */

        )
{

/* #pragma HLS INTERFACE mode=m_axi port=htb_buf1 bundle=gmem1 depth=DDR_WIDTH */
#pragma HLS INTERFACE mode=m_axi port=htb_buf0 bundle=gmem0
#pragma HLS INTERFACE mode=m_axi port=htb_buf1 bundle=gmem1 
#pragma HLS INTERFACE mode=m_axi port=htb_buf2 bundle=gmem2
#pragma HLS INTERFACE mode=m_axi port=htb_buf3 bundle=gmem3
#pragma HLS INTERFACE mode=m_axi port=res_buf bundle=gmem4
#pragma HLS INTERFACE mode=m_axi port=edge_buf bundle=gmem5
/* #pragma HLS alias ports=htb_buf0,htb_buf2,htb_buf3 distance=0 */
#pragma HLS alias ports=htb_buf0,htb_buf1,htb_buf2,htb_buf3 distance=0

/* #pragma HLS INTERFACE mode=axis port=stream_src */
/* #pragma HLS INTERFACE mode=axis port=stream_dst */
/* #pragma HLS INTERFACE mode=axis port=stream_src_l */
/* #pragma HLS INTERFACE mode=axis port=stream_dst_l */
#pragma HLS INTERFACE mode=s_axilite port=numQueryVert
#pragma HLS INTERFACE mode=s_axilite port=numQueryEdges
#pragma HLS INTERFACE mode=s_axilite port=numDataEdges
#pragma HLS INTERFACE mode=s_axilite port=return

#ifdef DEBUG_INTERFACE
#pragma HLS INTERFACE mode=s_axilite port=debif_endpreprocess
#endif /* DEBUG_INTERFACE */

#ifdef COUNT_ONLY
#pragma HLS INTERFACE mode=s_axilite port=result
#else
#pragma HLS INTERFACE mode=axis port=result
#endif /* COUNT_ONLY */

#ifdef DEBUG_STATS
/* statistic purposes */

    debug::init();

#endif /* DEBUG_STATS */

    QueryVertex qVertices0[MAX_QV], qVertices1[MAX_QV];
    TableDescriptor tDescriptors[MAX_TB];
    AdjHT hTables0[MAX_TB], hTables1[MAX_TB];
    unsigned short numTables = 0;
    unsigned long localResult = 0;

    preprocess<row_t,
        EDGE_WIDTH,
        COUNTER_WIDTH,
        DDR_BIT,
        VERTEX_WIDTH_BIT,
        LABEL_WIDTH,
        HASH_WIDTH_FIRST,
        HASH_WIDTH_SECOND,
        STREAM_DEPTH,
        HASHTABLES_SPACE,
        MAX_QUERY_VERTICES,
        MAX_TABLES>(
#ifdef DEBUG_INTERFACE
                debif_endpreprocess,
#endif /* DEBUG_INTERFACE */
                edge_buf,
                htb_buf0,
                qVertices0,
                qVertices1,
                tDescriptors,
                hTables0,
                hTables1,
                numTables,
                numQueryVert,
                numQueryEdges,
                numDataEdges);


    multiwayJoinWrap(
            htb_buf0,
            htb_buf1,
            htb_buf2,
            htb_buf3,
            res_buf,
            hTables0,
            hTables1,
            qVertices0,
            qVertices1,
            numQueryVert,
            localResult);

    result = localResult;

#ifdef DEBUG_STATS
    {
        using namespace debug;
        std::ofstream of("../../../../stats.txt", std::ofstream::app);

        unsigned long debug_total_reads = findmin_reads +
            readmin_reads + intersect_reads +
            hashtovid_reads + verify_reads;

        unsigned long debug_verify = solution_wrong + 
            solution_correct;

        unsigned long debug_probe = intersect_bit_falsepositive + 
            intersect_bit_falsenegative +
            intersect_bit_truepositive + 
            intersect_bit_truenegative;

        unsigned long debug_sol = intersect_sol_falsepositive + 
            intersect_sol_falsenegative +
            intersect_sol_truepositive + 
            intersect_sol_truenegative;
        
        unsigned int hw1, hw2, cnt;
        hw1 = H_W_1; 
        hw2 = H_W_2; 
        cnt = C_W;

        of << "DEBUG STATISTICS HW1: " << hw1 << " HW2: " << hw2
            << " CNT: " << cnt << std::endl << std::endl;

        of << "\tfindmin reads:     " << findmin_reads << "\t" <<
            findmin_reads * 100 / debug_total_reads << "%" << std::endl;
        of << "\treadmin reads:     " << readmin_reads << "\t" <<
            readmin_reads * 100 / debug_total_reads << "%"<< std::endl;
        of << "\tintersect reads:   " << intersect_reads << "\t" <<
            intersect_reads * 100 / debug_total_reads << "%"<< std::endl;
        of << "\tverify reads:      " << verify_reads << "\t" << 
            verify_reads * 100 / debug_total_reads << "%"<< std::endl;
        of << "\tTOTAL:             " <<  debug_total_reads << "\t100%\n" << 
            std::endl;
        of << "\tread per embedding: " << debug_total_reads / solution_correct << 
            std::endl;
        of << "\tindexed tables:    " << indexed_tables << std::endl;
        of << "\tindexing tables:   " << indexing_tables << std::endl;
        of << "\n\tsolution wrong:   " << solution_wrong << "\t" <<
            solution_wrong * 100 / debug_verify << "%" << std::endl;
        of << "\tsolution correct:   " << solution_correct << "\t" <<
            solution_correct * 100 / debug_verify << "%" << std::endl;
        of << "\tmax collisions:   " << max_collisions << std::endl;
        of << "\tavg collisions:   " << avg_collisions / numTables << std::endl;
        of << "\tverify reusage:   " << verify_reusage << std::endl << std::endl;
        of << "\tdata struct prob: " << std::setprecision(3) << std::endl <<
            "\t\tFP " << (float)intersect_bit_falsepositive / debug_probe << "\t" <<
            intersect_bit_falsepositive << std::endl <<
            "\t\tTP " << (float)intersect_bit_truepositive / debug_probe << "\t" <<
            intersect_bit_truepositive << std::endl <<
            "\t\tFN " << (float)intersect_bit_falsenegative / debug_probe << "\t" <<
            intersect_bit_falsenegative << std::endl <<
            "\t\tTN " << (float)intersect_bit_truenegative / debug_probe << "\t" <<
            intersect_bit_truenegative << std::endl << std::endl;
        of << "\tsolution prob: " << std::setprecision(3) << std::endl <<
            "\t\tFP " << (float)intersect_sol_falsepositive / debug_sol << "\t" <<
            intersect_sol_falsepositive << std::endl <<
            "\t\tTP " << (float)intersect_sol_truepositive / debug_sol << "\t" <<
            intersect_sol_truepositive << std::endl <<
            "\t\tFN " << (float)intersect_sol_falsenegative / debug_sol << "\t" <<
            intersect_sol_falsenegative << std::endl <<
            "\t\tTN " << (float)intersect_sol_truenegative / debug_sol << "\t" <<
            intersect_sol_truenegative << std::endl << std::endl;
        of.close();
    }

#endif
}

template <typename T>
void mwj_verify(
        AdjHT *hTables,
        QueryVertex *qVertices,
        T htb_buf,
        hls::stream<ap_uint<V_ID_W>> &stream_inter_in,
        hls::stream<bool> &stream_end_inter_in,
        hls::stream<ap_uint<16>> &stream_setinfo_in,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        hls::stream<bool> &stream_end_embed_in,
        hls::stream<bool, 1> &stream_stop,
        
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out,
        hls::stream<ap_uint<V_ID_W>> &stream_checked_out,
        hls::stream<bool> &stream_checked_end_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in("Verify edge hash_in");
    hls::stream<ap_uint<64>> stream_hash_out("Verify edge hash_out");
    ap_uint<8> curQV;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool stop, last;
    ap_uint<16> setinfo;

#ifdef DEBUG_STATS
    std::ofstream of("../../../../verstats.txt", std::ofstream::app);
#endif

    while(1){
        if (stream_end_embed_in.read_nb(last)){
            curQV = 0;

VERIFY_EDGE_COPYING_EMBEDDING_LOOP:
            while(!last){
                curEmb[curQV] = stream_embed_in.read();
                stream_embed_out.write(curEmb[curQV]);
                stream_end_embed_out.write(false);
                curQV++;
                last = stream_end_embed_in.read();
            }
            stream_end_embed_out.write(true);

            setinfo = stream_setinfo_in.read();
            last = stream_end_inter_in.read();
VERIFY_CHECK_LOOP:
            while(!last){
                ap_uint<V_ID_W> vToVerify = stream_inter_in.read();
                bool checked = true;
                int g = 0;

VERIFY_CHECK_EDGE_LOOP:
                for(; g < qVertices[curQV].numTablesIndexed && checked; g++){
                    ap_uint<H_W_1> index1, index1_f;
                    ap_uint<H_W_2> index2, index2_f;
                    uint8_t tableIndex = qVertices[curQV].tables_indexed[g];
                    uint8_t ivPos = qVertices[curQV].vertex_indexing[g];
    
                    if (setinfo.range(7, 0) != tableIndex || setinfo.range(15, 8) != ivPos) {
                        checked = false;
                        stream_hash_in.write(curEmb[ivPos]);
                        stream_hash_in.write(vToVerify);
                        xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
                        index1 = stream_hash_out.read().range(H_W_1 - 1, 0);
                        xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
                        index2 = stream_hash_out.read().range(H_W_2 - 1, 0);

                        ap_uint<(1UL << C_W)> start_off = 0;
                        ap_uint<(1UL << C_W)> end_off = 0;

                        if (index2 != 0){
                            index1_f = index1;
                            index2_f = index2 - 1;
                        } else if (index2 == 0 && index1 != 0){
                            index1_f = index1 - 1;
                            index2_f = (1UL << H_W_2)-1;
                        }

                        if (!(index2 == 0 && index1 == 0)){
                            start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                                    0,
                                    index1_f,
                                    index2_f,
                                    htb_buf,
                                    hTables[tableIndex].start_offset);
                            start_off = start_off.range((1UL << C_W) - 2, 0);

#ifdef DEBUG_STATS
                            debug::verify_reads++;
#endif
                        }

                        end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                                0,
                                index1,
                                index2,
                                htb_buf,
                                hTables[tableIndex].start_offset);
                        end_off = end_off.range((1UL << C_W) - 2, 0);

#ifdef DEBUG_STATS
                        debug::verify_reads++;
                        /* of << "(" << vToVerify << ") " << start_off + hTables[tableIndex].start_edges << */
                        /* " - " << end_off - start_off << std::endl; */
#endif
VERIFY_READ_MEMORY_LOOP:
                        for (; start_off < end_off; start_off++){
                            ap_uint<2 * V_ID_W> edge = read_table<32, 1, 0, E_W>(
                                    1,
                                    start_off,
                                    0,
                                    htb_buf,
                                    hTables[tableIndex].start_edges);

                            ap_uint<V_ID_W> vertexIndexed, vertexIndexing;
                            vertexIndexed = edge.range(V_ID_W - 1, 0);
                            vertexIndexing = edge.range(2*V_ID_W - 1, V_ID_W);
                            if (vertexIndexing == curEmb[ivPos] &&
                                    vertexIndexed == vToVerify){
                                checked = true;
                            }
#ifdef DEBUG_STATS
                            debug::verify_reads++;
#endif
                        }
                    }
                }

/* #ifdef DEBUG_STATS */
/* of << checked << std::endl << std::endl; */
/* #endif */

                if (checked){
                    stream_checked_out.write(vToVerify);
                    stream_checked_end_out.write(false);
#ifdef DEBUG_STATS
/* if (g > 0) */
                        debug::solution_correct++;
                } else { 
                    debug::solution_wrong++; 
#endif
                }

                last = stream_end_inter_in.read();
            }
            stream_checked_end_out.write(true);
        }
        if (stream_stop.read_nb(stop))
            break;
    }
}

/* Multiway join propose: retrieves the tables in which
 * the current query vertex is involved, computes the sizes
 * and keep track of the smallest one */
void mwj_propose_findmin(
        hls::stream<ap_uint<V_ID_W>> &stream_embed_in,
        AdjHT *hTables,
        QueryVertex *qVertices,
        ap_uint<DDR_W> *htb_buf,
        hls::stream<bool, 1> &stream_stop,

        hls::stream<ap_uint<32>> &stream_minread_out,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_out,
        hls::stream<ap_uint<V_ID_W>> &stream_embed_out,
        hls::stream<bool> &stream_end_embed_out)
{
    
    hls::stream<ap_uint<V_ID_W>> stream_hash_in("Propose findmin hash_in");
    hls::stream<ap_uint<64>> stream_hash_out("Propose readmin hash_in");
    ap_uint<8> tableIndex;
    ap_uint<V_ID_W> curQV;
    ap_uint<32> minSize;
    ap_uint<32> minStart, minEnd, minOff;
    ap_uint<8 + 1 + V_ID_W> minData;
    ap_uint<V_ID_W> curEmb[MAX_QV];
    bool stop;

    while (1) {
        if (stream_embed_in.read_nb(curQV)){
            tableIndex = 0;
            minSize = (1UL << 32) - 1;

PROPOSE_COPYING_EMBEDDING_LOOP:
            for (int g = 0; g < curQV; g++){
#pragma HLS pipeline II=1
                curEmb[g] = stream_embed_in.read();
                stream_embed_out.write(curEmb[g]);
                stream_end_embed_out.write(false);
            }
            stream_end_embed_out.write(true);

            /* Find sizes of sets in which the current query vertex
             * is indexed by an other query vertex */
PROPOSE_TBINDEXED_LOOP:
            for(int g = 0; g < qVertices[curQV].numTablesIndexed; g++){
#pragma HLS pipeline II=1
                ap_uint<(1UL << C_W)> start_off = 0;
                ap_uint<(1UL << C_W)> end_off;
                tableIndex = qVertices[curQV].tables_indexed[g];
                uint8_t ivPos = qVertices[curQV].vertex_indexing[g];

                stream_hash_in.write(curEmb[ivPos]);
                xf::database::hashLookup3<V_ID_W>(stream_hash_in, stream_hash_out);
                ap_uint<H_W_1> index = stream_hash_out.read().range(H_W_1 - 1, 0);

                if (index != 0){ 
                    start_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                            index - 1,
                            (1UL << H_W_2) - 1,
                            htb_buf,
                            hTables[tableIndex].start_offset);
                    start_off = start_off.range((1UL << C_W) - 2, 0);
                }

                end_off = read_table<H_W_1, H_W_2, H_W_2, C_W>(
                        index,
                        (1UL << H_W_2) - 1,
                        htb_buf,
                        hTables[tableIndex].start_offset);
                end_off = end_off.range((1UL << C_W) - 2, 0);

#ifdef DEBUG_STATS
                debug::findmin_reads += 2;
#endif
                if ((end_off - start_off) < minSize) {
                    minSize = end_off - start_off;
                    minOff = hTables[tableIndex].start_edges;
                    minStart = start_off;
                    minEnd = end_off;
                    minData.range(7, 0) = tableIndex;
                    minData.set(8);
                    minData.range(V_ID_W + 8, 9) = curEmb[ivPos];
                }
            }
            stream_minread_out.write(minStart);
            stream_minread_out.write(minEnd);
            stream_minread_out.write(minOff);
            stream_min_out.write(minData);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

/* Read the minimum set from memory */
void mwj_propose_readmin(
        hls::stream<ap_uint<32>> &stream_minread_in,
        hls::stream<ap_uint<V_ID_W + 1 + 8>> &stream_min_in,
        ap_uint<DDR_W> *htb_buf,
        /* hls::burst_maxi<row_t> htb_buf, */
        hls::stream<bool, 1> &stream_stop,

        hls::stream<ap_uint<V_ID_W>> &stream_min_out,
        hls::stream<bool> &stream_end_out)
{
    hls::stream<ap_uint<V_ID_W>> stream_hash_in("Read min hash_in");
    hls::stream<ap_uint<64>> stream_hash_out("Read min hash_out");
    ap_uint<32> minStart; 
    ap_uint<32> minEnd;
    ap_uint<32> minOff;
    ap_uint<V_ID_W + 1 + 8> minData;
    ap_uint<V_ID_W * 2> edge;
    ap_uint<V_ID_W> vertex, vertexCheck;
    bool stop;

    while(1) {

        if (stream_minread_in.read_nb(minStart)){
            minEnd = stream_minread_in.read();
            minOff = stream_minread_in.read();
            minData = stream_min_in.read();
            unsigned int rowstart = minOff + (minStart >> (DDR_BIT - E_W));
            unsigned int rowend = minOff + (minEnd >> (DDR_BIT - E_W));
            unsigned int window_left = minStart.range((DDR_BIT - E_W) - 1, 0);
            unsigned int window_right = minEnd.range((DDR_BIT - E_W) - 1, 0) + 
                (rowend - rowstart) * EDGE_ROW;
            unsigned int cnt = 0;

            /* std::cout << "\tMinimum set is " << minData.range(7, 0) << std::endl; */
            /* std::cout << "\tReading from " << rowstart << " to " << rowend */
            /* << " window " << window_left << "-" << window_right << std::endl; */
            /* htb_buf.read_request(rowstart, rowend-rowstart+1); */

PROPOSE_READ_MIN_INDEXED_LOOP:
            for (int g = rowstart; g <= rowend; g++){
                row_t row = htb_buf[g];
                /* row_t row = htb_buf.read(); */
                for (int i = 0; i < EDGE_ROW; i++, cnt++){
                    if (cnt >= window_left && cnt < window_right){
                        edge = row.range((1UL << E_W) - 1, 0);
                        vertexCheck = edge.range(V_ID_W * 2 - 1, V_ID_W);
                        vertex = edge.range(V_ID_W - 1, 0);
                        if (minData.range(V_ID_W + 8, 9) == vertexCheck){
                            stream_min_out.write(vertex);
                            stream_end_out.write(false);
                        }
                    }
                    row >>= (1UL << E_W);
                }
#ifdef DEBUG_STATS
                debug::readmin_reads++;
#endif
            }
            stream_end_out.write(true);
        }

        if (stream_stop.read_nb(stop))
            break;
    }
}

