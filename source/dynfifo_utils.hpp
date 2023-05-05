#ifndef DYNFIFO_UTILS_HPP
#define DYNFIFO_UTILS_HPP

#ifndef __SYNTHESIS__
#define HLS_STREAM_THREAD_SAFE
#include <thread>
#endif /* __SYNTHESIS__ */

#include <ap_int.h>
#include "hls_task.h"
#include <hls_stream.h>

/**
 * Send stop signals to the tasks.
 */
void MMU_stop(
        hls::stream<bool> &stop_req,
        hls::stream<bool> *stop_streams)
{
    stop_req.read();
    for(int g = 0; g < 3; g++){
        stop_streams[g].write(true);
    }
}

/**
 * Pack data in DDR word. 
 */
template<typename DATA_T,
    typename DDR_WORD_T,
    size_t DATA_WIDTH,
    size_t DATA_PER_WORD,
    size_t MEM_WIDTH>
void MMU_pack(
        hls::stream<DATA_T> &store_stream,
        hls::stream<DDR_WORD_T> &packed_stream)
{
#pragma HLS pipeline II=1

    static unsigned int counter {DATA_PER_WORD - 1};
    static DDR_WORD_T word;

    word <<= DATA_WIDTH;
    word(DATA_WIDTH - 1, 0) = store_stream.read();

    if (counter == 0){
        packed_stream.write(word);
        counter = DATA_PER_WORD - 1;
    } else {
        counter--;
    }
}

/**
 * Unpack a DDR word in data. 
 */
template<typename DATA_T,
    typename DDR_WORD_T,
    size_t DATA_WIDTH,
    size_t DATA_PER_WORD,
    size_t MEM_WIDTH>
void MMU_unpack(
        hls::stream<DATA_T> &load_stream,
        hls::stream<DDR_WORD_T> &packed_stream)
{
#pragma HLS pipeline II=1

    static unsigned char counter {0};
    static DDR_WORD_T word;

    if (counter == 0){
        word = packed_stream.read();
        counter = DATA_PER_WORD - 1;
    } else {
        counter--;
    }

    load_stream.write(word(MEM_WIDTH - 1, (DATA_PER_WORD - 1) * DATA_WIDTH));
    word <<= DATA_WIDTH;

}

template <typename DDR_WORD_T,
            size_t BURST_SIZE>
void MMU_burst(
        DDR_WORD_T* mem,
        hls::stream<DDR_WORD_T>     &in_stream,
        hls::stream<DDR_WORD_T>     &out_stream,
        hls::stream<bool>           &stop_stream,
        hls::stream<bool>           &req_stream,
        hls::stream<bool>           &resp_stream,
        hls::stream<unsigned int>   &p_stream)
{
    unsigned int p;
    bool stop, req;

    while (true){
        if (req_stream.read_nb(req)){
            p = p_stream.read();
            if (req) {
                for (int g = 0; g < BURST_SIZE; g++){
#pragma HLS unroll
                    mem[p + g] = in_stream.read();
                }
            } else {
                for (int g = 0; g < BURST_SIZE; g++){
#pragma HLS unroll
                    out_stream.write(mem[p + g]);
                }
            }
            ap_wait();
            ap_wait();
            resp_stream.write(req);
        }
        
        if (stop_stream.read_nb(stop))
            break;
    }
}

/** 
 * Main process with II=1 in charge of simulate a normal hls::stream,
 * in bypass mode shift data from in_stream to out_stream,
 * otherwise move data from in_stream to store_stream and from
 * load_stream to out_stream.
 *
 * DATA_T: data type,
 * DATA_PER_WORD: number of data inside a packet
 */
template<typename DATA_T,
    size_t DATA_PER_WORD>
void MMU_fast(
        hls::stream<DATA_T> &in_stream,
        hls::stream<DATA_T> &out_stream,
        hls::stream<DATA_T> &store_stream,
        hls::stream<DATA_T> &load_stream,
        hls::stream<bool> &stop_req_stream)
{

    enum State_fast { bypass, stall, shift, stall_ddr };

    unsigned int p_counter {DATA_PER_WORD};
    unsigned int ddr_data {0};
    DATA_T buff_r;
    DATA_T buff_w;
    bool flagbuff_r {false};
    bool flagbuff_w {false};
    State_fast state {bypass};

/* #ifndef __SYNTHESIS__ */
/* std::this_thread::sleep_for(std::chrono::milliseconds(100)); */
/* #endif */

    while(true) {
#pragma HLS pipeline II=1 
        switch (state) {
            
            case bypass: 

                if (in_stream.read_nb(buff_w)){
                    if (!out_stream.write_nb(buff_w)){
                        /* Out stream full and input stream not empty */
                        state = stall;
                        flagbuff_w = true;
                        ddr_data = 1;
                    }
                }

                break;

            case stall:
                if (in_stream.size() < DATA_PER_WORD){
                    if (out_stream.write_nb(buff_w)){
                        /* Out stream not full and packet not ready,
                         * return to bypass state */
                        flagbuff_w = false;
                        state = bypass;
                        ddr_data = 0;
                    }
                } else {
                    /* Out stream full and packet ready */
                    state = shift;
                    p_counter = DATA_PER_WORD;
                }
                break;

            case shift:
                if (!flagbuff_r){
                    buff_r = in_stream.read();
                    flagbuff_r = true;
                }

                if (store_stream.write_nb(buff_r)) {
                    p_counter--;
                    ddr_data++;
                    flagbuff_r = false;
                    if (p_counter == 0) {
                        state = stall_ddr;
                    }
                }

                if (flagbuff_w){
                    if (out_stream.write_nb(buff_w)){
                        ddr_data--;
                        if (!load_stream.read_nb(buff_w)){
                            flagbuff_w = false;
                        }   
                    }
                } else {
                    if (load_stream.read_nb(buff_w)){
                        flagbuff_w = true;
                    }   
                }

                break;

            case stall_ddr:
                if (in_stream.size() >= DATA_PER_WORD){
                    state = shift;
                    p_counter = DATA_PER_WORD;
                }

                if (flagbuff_w){
                    if (out_stream.write_nb(buff_w)){
                        ddr_data--;
                        if (ddr_data == 0){
                            state = bypass;
                        }
                        if (!load_stream.read_nb(buff_w)){
                            flagbuff_w = false;
                        }   
                    }
                } else {
                    if (load_stream.read_nb(buff_w)){
                        flagbuff_w = true;
                    }   
                }
        }

        bool req_s;
        if (stop_req_stream.read_nb(req_s)){
            break;
        }
    }
}
        
/**
 * Process handeling the communication with the ddr.
 * Gathers packet to do burst transactions.
 *
 * DDR_WORD_T: ddr word data type,
 * FIFO_SIZE: internal fifo size storing packets,
 * DDR_WORDS: number of words allocated for the dynfifo
 * BURST_SIZE: size of a burst transaction,
 */
template <typename DDR_WORD_T,
            size_t FIFO_SIZE,
            size_t DDR_WORDS,
            size_t BURST_SIZE>
void MMU_slow(
        hls::stream<DDR_WORD_T> &in_stream,
        hls::stream<DDR_WORD_T> &out_stream,
        hls::stream<unsigned int> &ptr_stream,
        hls::stream<bool> &req_stream,
        hls::stream<bool> &resp_stream,
        hls::stream<bool> &stop_req_stream)
{

    enum State_slow { bypass, stall, burst_read, burst_write, stall_ddr };
    unsigned int mem_tail {0};
    unsigned int mem_head {0};
    unsigned long max_space {0};
    unsigned long space_used {0};
    State_slow state {bypass};

MMU_SLOW_TASK_LOOP:
    while(true) {
        switch (state) {

            case bypass:
                {
                    if (!in_stream.size()){
                        if (out_stream.full()){
                            state = stall;
                        } else {
                            out_stream.write(in_stream.read());
                        }
                    }
                    break;
                }

            case stall:
                {
                    if (in_stream.size() < BURST_SIZE){
                        if (!out_stream.full()){
                            state = bypass;
                        }
                    } else {
                        state = burst_write;
                    }
                    break;
                }

            case burst_write:
                {
                    bool dep, rsp;
                    req_stream.write(true);
                    dep = ptr_stream.write_dep(mem_head, false);
                    mem_head = (mem_head + BURST_SIZE) % DDR_WORDS;
                    space_used += BURST_SIZE;
                    if (space_used > max_space)
                        max_space = space_used;
                    state = stall_ddr;
                    resp_stream.read_dep(rsp, dep);
                    break;
                }

            case burst_read:
                {
                    bool dep, rsp;
                    req_stream.write(false);
                    dep = ptr_stream.write_dep(mem_tail, false);
                    mem_tail = (mem_tail + BURST_SIZE) % DDR_WORDS;
                    space_used -= BURST_SIZE;
                    state = (mem_head == mem_tail)? bypass: stall_ddr;
                    resp_stream.read_dep(rsp, dep);
                    break;
                }

            case stall_ddr:
                {
                    if (in_stream.size() >= BURST_SIZE){
                        state = burst_write;
                    } else {
                        if (out_stream.size() <= (FIFO_SIZE - BURST_SIZE)){
                            state = burst_read;
                        }
                    }
                }
        }

        bool req_s;
        if (stop_req_stream.read_nb(req_s)){
            break;
        }
    }

}

/**
 * DATA_T: data type,
 * FIFO_SIZE_T: input/output fifo size,
 * FIFO_SIZE_P: internal fifo size storing packets,
 * MEM_WIDTH: number of bit for a word,
 * BURST_SIZE: size of a burst transaction,
 * DDR_WORDS: number of words allocated for the dynfifo
 *
 * Starts the processes needed by the dynamic fifo.
 */
template<typename DATA_T,
    typename DDR_WORD_T, 
    size_t FIFO_SIZE_T,
    size_t FIFO_SIZE_P,
    size_t MEM_WIDTH,
    size_t BURST_SIZE,
    size_t DDR_WORDS>
void dynfifo_init(
        hls::stream<DATA_T> &in_stream,
        hls::stream<DATA_T> &out_stream,
        hls::stream<bool> &stop_req_stream,
        DDR_WORD_T mem[DDR_WORDS]){

#pragma HLS inline 
    
    static_assert((DDR_WORDS % BURST_SIZE) == 0,
            "The number of ddr words must be a multiple of a burst transaction");
    static_assert(FIFO_SIZE_P >= BURST_SIZE,
            "The packet fifo size must be equal or greater than a burst transaction");
    static_assert(FIFO_SIZE_T >= (MEM_WIDTH / (8*sizeof(DATA_T))),
            "The input/output fifo size must be equal or greater than a packet");
    
    const size_t DATA_WIDTH = sizeof(DATA_T)*8;
    const size_t DATA_PER_WORD = MEM_WIDTH / DATA_WIDTH;

    hls_thread_local hls::stream<DATA_T, FIFO_SIZE_T> store_stream;
    hls_thread_local hls::stream<DATA_T, FIFO_SIZE_T> load_stream;
    hls_thread_local hls::stream<DDR_WORD_T, FIFO_SIZE_P> store_p_stream;
    hls_thread_local hls::stream<DDR_WORD_T, FIFO_SIZE_P> load_p_stream;

    hls_thread_local hls::stream<unsigned int, 2> ptr_stream;
    hls_thread_local hls::stream<bool, 2> req_stream;
    hls_thread_local hls::stream<bool, 2> resp_stream;
    
    hls_thread_local hls::stream<bool, 1> streams_stop[3];

#ifdef __SYNTHESIS__    

    MMU_fast<DATA_T,
            DATA_PER_WORD>(
            in_stream,
            out_stream,
            store_stream,
            load_stream,
            streams_stop[0]); 

    MMU_slow<DDR_WORD_T,
            FIFO_SIZE_P,
            DDR_WORDS,
            BURST_SIZE>(
            store_p_stream,
            load_p_stream,
            ptr_stream,
            req_stream,
            resp_stream,
            streams_stop[1]); 

    MMU_burst<DDR_WORD_T,
            BURST_SIZE>(
            mem,
            store_p_stream,
            load_p_stream,
            streams_stop[2], 
            req_stream,
            resp_stream,
            ptr_stream);

#else

    for (int g = 0; g < 3; g++) 
        hls::stream_globals::incr_task_counter();

    std::thread mmu_fast_t(
            MMU_fast<DATA_T,
            DATA_PER_WORD>,
            std::ref(in_stream),
            std::ref(out_stream),
            std::ref(store_stream),
            std::ref(load_stream),
            std::ref(streams_stop[0])); 

    std::thread mmu_slow_t(
            MMU_slow<DDR_WORD_T,
            FIFO_SIZE_P,
            DDR_WORDS,
            BURST_SIZE>,
            std::ref(store_p_stream),
            std::ref(load_p_stream),
            std::ref(ptr_stream),
            std::ref(req_stream),
            std::ref(resp_stream),
            std::ref(streams_stop[1]));

    std::thread mmu_burst_t(
            MMU_burst<DDR_WORD_T,
            BURST_SIZE>,
            mem,
            std::ref(store_p_stream),
            std::ref(load_p_stream),
            std::ref(streams_stop[2]),
            std::ref(req_stream),
            std::ref(resp_stream),
            std::ref(ptr_stream));
    
    mmu_fast_t.detach();
    mmu_slow_t.detach();
    mmu_burst_t.detach();

#endif /* __SYNTHESIS__ */
    
    hls_thread_local hls::task mmu_pack(
            MMU_pack<DATA_T,
            DDR_WORD_T,
            DATA_WIDTH,
            DATA_PER_WORD,
            MEM_WIDTH>,
            store_stream,
            store_p_stream);

    hls_thread_local hls::task mmu_unpack(
            MMU_unpack<DATA_T,
            DDR_WORD_T,
            DATA_WIDTH,
            DATA_PER_WORD,
            MEM_WIDTH>,
            load_stream,
            load_p_stream);

    hls_thread_local hls::task mmu_stop(
            MMU_stop,
            stop_req_stream,
            streams_stop);
}

#endif /* DYNFIFO_UTILS_HPP */
