#pragma once

namespace debug {

    static unsigned long batch_reads;
    static unsigned long findmin_reads;
    static unsigned long readmin_reads;
    static unsigned long readmin_n_sets;
    static unsigned long readmin_vstream;
    static unsigned long bloom_filter;
    static unsigned long intersect_reads;
    static unsigned long intersect_filter;
    static unsigned long verify_filter;
    static unsigned long verify_reads;
    static unsigned long homomo_trashed;
    static unsigned long miss_indexing;
    static unsigned long start_set;
    static unsigned long max_collisions;
    static unsigned long cache_req_prop;
    static unsigned long cache_req_inter;
    static unsigned long cache_req_verify;
    static unsigned long hash_collisions;
    static float bloom_fullness;
    static float cache_hit_prop;
    static float cache_hit_inter;
    static float cache_hit_verify;
    
    static void init(){
        batch_reads      = 0;
        findmin_reads    = 0;
        readmin_reads    = 0;
        readmin_n_sets   = 0;
        readmin_vstream  = 0;
        bloom_filter     = 0;
        bloom_fullness   = 0;
        intersect_reads  = 0;
        intersect_filter = 0;
        verify_filter    = 0;
        verify_reads     = 0;
        homomo_trashed   = 0;
        start_set        = 0;
        miss_indexing    = 0;
        max_collisions   = 0;
        hash_collisions  = 0;
        cache_hit_prop   = 0;
        cache_hit_inter  = 0;
        cache_hit_verify = 0;
        cache_req_prop   = 0;
        cache_req_inter  = 0;
        cache_req_verify = 0;
    }

    static void print(
        const unsigned char hash1_w,
        const unsigned char hash2_w
    ){
        std::ofstream debof("../../../../stats.txt", std::ofstream::app);

        unsigned long mem_reads =  batch_reads +
        findmin_reads    +
        readmin_reads    +
        intersect_reads  +
        verify_reads;     
        
        unsigned int cnt, k_fun;
        cnt = COUNTER_WIDTH;
        k_fun = (1UL << K_FUNCTIONS);
        verify_filter -= intersect_filter;

        debof << "DEBUG STATISTICS HW1: " << hash1_w << " HW2: " << hash2_w
            << " CNT: " << cnt << " K_FUN: " << k_fun << std::endl << std::endl;

        debof << "\tbatch reads:     " << batch_reads << "\t" << 
            batch_reads * 100 / mem_reads << "%"<< std::endl;
        debof << "\tfindmin reads:   " << findmin_reads << "\t" <<
            findmin_reads * 100 / mem_reads << "%" << std::endl;
        debof << "\treadmin reads:   " << readmin_reads << "\t" <<
            readmin_reads * 100 / mem_reads << "%"<< std::endl;
        debof << "\tintersect reads: " << intersect_reads << "\t" << 
            intersect_reads * 100 / mem_reads << "%"<< std::endl;
        debof << "\tverify reads:    " << verify_reads << "\t" << 
            verify_reads * 100 / mem_reads << "%"<< std::endl;
        debof << "\tTOTAL:           " <<  mem_reads << "\t100%\n" << 
            std::endl;

#if CACHE_ENABLE
        mem_reads = mem_reads - cache_hit_verify - cache_hit_inter - cache_hit_prop;
        intersect_reads -= cache_hit_inter;
        verify_reads -= cache_hit_verify;
        findmin_reads -= cache_hit_prop;
        debof << "With Caching on:" << std::endl;
        debof << "\tbatch reads:     " << batch_reads << "\t" << 
            batch_reads * 100 / mem_reads << "%"<< std::endl;
        debof << "\tfindmin reads:   " << findmin_reads << "\t" <<
            findmin_reads * 100 / mem_reads << "%" << std::endl;
        debof << "\treadmin reads:   " << readmin_reads << "\t" <<
            readmin_reads * 100 / mem_reads << "%"<< std::endl;
        debof << "\tintersect reads: " << intersect_reads << "\t" << 
            intersect_reads * 100 / mem_reads << "%"<< std::endl;
        debof << "\tverify reads:    " << verify_reads << "\t" << 
            verify_reads * 100 / mem_reads << "%"<< std::endl;
        debof << "\tTOTAL:           " <<  mem_reads << "\t100%\n" << 
            std::endl;

        debof << "\tmean cache hit propose: " << cache_hit_prop / cache_req_prop
            << std::endl;
        debof << "\tcache reqs propose:     " << cache_req_prop << std::endl; 
        debof << "\tmean cache hit inter:   " << cache_hit_inter / cache_req_inter  
            << std::endl;
        debof << "\tcache reqs inter:       " << cache_req_inter  << std::endl;
        debof << "\tmean cache hit verify:  " << cache_hit_verify / cache_req_verify  
            << std::endl;
        debof << "\tcache reqs verify:      " << cache_req_verify << std::endl; 
#endif

        debof << "\tmax collisions:      " << max_collisions << std::endl;
//        debof << "\tavg collisions:      " << hash_collisions / (float)counters << std::endl;
        debof << "\tbloom fullness:      " << bloom_fullness << std::endl;
        debof << "\tbloom filter:        " << bloom_filter << std::endl;
        debof << "\tintersect filter:    " << intersect_filter << std::endl;
        debof << "\tverify filter:       " << verify_filter << std::endl;
        debof << "\thomomorphism filter: " << homomo_trashed << std::endl;
        debof << "\tstart_set:           " << start_set << std::endl;
        debof << "\tmiss indexing        " << miss_indexing << std::endl;
        debof << "\taverage set:         " << readmin_vstream / (float)readmin_n_sets 
            << std::endl << std::endl;

        debof.close();
    }
};
