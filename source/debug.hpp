#pragma once

namespace debug {

    static unsigned long hashtable_reads {0};
    static unsigned long edge_reads {0};
    static unsigned long bloom_reads {0};
    static unsigned long tuple_ver {0};
    static unsigned long solution_correct {0};
    static unsigned long solution_wrong {0};
    static unsigned long intersect_filter {0};
    static unsigned long intersect_bit_truepositive {0};
    static unsigned long intersect_bit_falsepositive {0};
    static unsigned long intersect_bit_missindexing {0};
    static unsigned long empty_sol {0};
    static unsigned long embeddings {0};
    static unsigned long max_collisions {0};
    static unsigned long n_sets = {0};
    static unsigned long set = {0};
    static float avg_collisions {0};
    static float bloom_fullness {0};
    static unsigned long bloom_filter {0};
    static float cache_hit_0 {0};
    static float cache_hit_1 {0};
    static unsigned long cache_req_0 {0};
    static unsigned long cache_req_1 {0};
    static unsigned long batches {0};

    static void init(){
        hashtable_reads = 0;
        bloom_reads = 0;
        edge_reads = 0;
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
        tuple_ver = 0;
    }
};
