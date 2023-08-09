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

    enum State_slow
    {
        bypass = 0,
        stall = 1,
        shift = 2,
        stall_ddr = 3
    };
    typedef ap_uint<2> state_type;
    state_type state = bypass;

    unsigned int p_counter {DATA_PER_WORD};
    unsigned int ddr_data {0};
    DATA_T buff_r;
    DATA_T buff_w;
    bool flagbuff_r {false};
    bool flagbuff_w {false};

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
        DDR_WORD_T* mem,
        const unsigned long space,
        unsigned long &diagnostic,
        hls::stream<DDR_WORD_T> &in_stream,
        hls::stream<DDR_WORD_T> &out_stream,
        hls::stream<bool> &stop_req_stream,
        hls::stream<bool> &overflow_stream)
{

    enum State_slow
    {
        bypass = 0,
        stall = 1,
        burst_read = 2,
        burst_write = 3,
        stall_ddr = 4
    };
    typedef ap_uint<3> state_type;
    state_type state = bypass;

    unsigned int mem_tail {0};
    unsigned int mem_head {0};
    unsigned long max_space {0};
    unsigned long space_used {0};

MMU_SLOW_TASK_LOOP:
    while(true) {
#pragma HLS pipeline II=40

        switch (state) {

            case bypass: 
                if (!in_stream.empty()){
                    if (out_stream.full()){
                        state = stall;
                    } else {
                        out_stream.write(in_stream.read());
                    }
                }
                break;

            case stall:
                if (in_stream.size() < BURST_SIZE){
                    if (!out_stream.full()){
                        state = bypass;
                    }
                } else {
                    state = burst_write;
                }
                break;

            case burst_write:
                for (int g = 0; g < BURST_SIZE; g++){
#pragma HLS unroll
                    mem[mem_head + g] = in_stream.read();
                }
                mem_head = (mem_head + BURST_SIZE);
                if (mem_head == space)
                    mem_head = 0;

                space_used += BURST_SIZE;
                if (space_used > max_space){
                    max_space = space_used;
                }

                if (space_used > space){
                    overflow_stream.write(true);
                }
                state = stall_ddr;
                break;

            case burst_read:
                for (int g = 0; g < BURST_SIZE; g++){
#pragma HLS unroll
                    out_stream.write(mem[mem_tail + g]);
                }
                mem_tail = (mem_tail + BURST_SIZE);
                if (mem_tail == space)
                    mem_tail = 0;
                space_used -= BURST_SIZE;
                state = (mem_head == mem_tail)? bypass: stall_ddr;
                break;

            case stall_ddr:
                if (in_stream.size() >= BURST_SIZE){
                    state = burst_write;
                } else {
                    if (out_stream.size() <= (FIFO_SIZE - BURST_SIZE)){
                        state = burst_read;
                    }
                }
        }

        bool req_s;
        if (stop_req_stream.read_nb(req_s)){
            diagnostic = max_space;
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
        DDR_WORD_T mem[DDR_WORDS],
        unsigned long &diagnostic,
        const unsigned long space,
        hls::stream<DATA_T> &in_stream,
        hls::stream<DATA_T> &out_stream,
        hls::stream<bool> &stop_req_stream_fast,
        hls::stream<bool> &stop_req_stream_slow,
        hls::stream<bool> &overflow)
{

#pragma HLS inline 
#pragma HLS dataflow
    
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

#ifdef __SYNTHESIS__    

    MMU_fast<DATA_T,
            DATA_PER_WORD>(
            in_stream,
            out_stream,
            store_stream,
            load_stream,
            stop_req_stream_fast); 

    MMU_slow<DDR_WORD_T,
            FIFO_SIZE_P,
            DDR_WORDS,
            BURST_SIZE>(
            mem,
            space,
            diagnostic,
            store_p_stream,
            load_p_stream,
            stop_req_stream_slow,
            overflow);

#else

    for (int g = 0; g < 2; g++) 
        hls::stream_globals::incr_task_counter();

    std::thread mmu_fast_t(
            MMU_fast<DATA_T,
            DATA_PER_WORD>,
            std::ref(in_stream),
            std::ref(out_stream),
            std::ref(store_stream),
            std::ref(load_stream),
            std::ref(stop_req_stream_fast)); 

    std::thread mmu_slow_t(
            MMU_slow<DDR_WORD_T,
            FIFO_SIZE_P,
            DDR_WORDS,
            BURST_SIZE>,
            mem,
            space,
            std::ref(diagnostic),
            std::ref(store_p_stream),
            std::ref(load_p_stream),
            std::ref(stop_req_stream_slow),
            std::ref(overflow));

    mmu_fast_t.detach();
    mmu_slow_t.detach();

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

}

#endif /* DYNFIFO_UTILS_HPP */
