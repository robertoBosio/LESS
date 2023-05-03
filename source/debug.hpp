#pragma once

namespace debug {

    unsigned long findmin_reads {0};
    unsigned long readmin_reads {0};
    unsigned long intersect_reads {0};
    unsigned long verify_reads {0};
    
    unsigned long solution_correct {0};
    unsigned long solution_wrong {0};

    unsigned long intersect_filter {0};
    unsigned long intersect_bit_truepositive {0};
    unsigned long intersect_bit_falsepositive {0};
    unsigned long intersect_bit_missindexing {0};
    
    unsigned long empty_sol {0};
    unsigned long embeddings {0};
    unsigned long max_collisions {0};
    unsigned long n_sets = {0};
    unsigned long set = {0};
    float avg_collisions {0};
    float bloom_fullness {0};
    unsigned long bloom_filter {0};

    float cache_hit_0 {0};
    float cache_hit_1 {0};
    unsigned long cache_req_0 {0};
    unsigned long cache_req_1 {0};
    unsigned long batches {0};

    static void init(){
        findmin_reads  = 0;
        readmin_reads  = 0;
        intersect_reads  = 0;
        verify_reads  = 0;
        solution_correct  = 0;
        solution_wrong  = 0;
        embeddings = 0;
        max_collisions = 0;
        avg_collisions = 0;
        intersect_filter = 0;
        intersect_bit_truepositive = 0;
        intersect_bit_falsepositive = 0;
        intersect_bit_missindexing = 0;
        empty_sol = 0;
        bloom_filter = 0;
        bloom_fullness = 0;
        cache_hit_0 = 0;
        cache_hit_1 = 0;
        cache_req_0 = 0;
        cache_req_1 = 0;
        batches = 0;
        n_sets = 0;
        set = 0;
    }
};
