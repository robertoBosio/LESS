#pragma once

namespace debug
{

      static unsigned long batch_reads;
      static unsigned long findmin_reads;
      static unsigned long readmin_counter_reads;
      static unsigned long readmin_edge_reads;
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
      static unsigned long cache_req_findmin;
      static unsigned long cache_req_read_c;
      static unsigned long cache_req_read_e;
      static unsigned long cache_req_inter;
      static unsigned long cache_req_verify;
      static unsigned long hash_collisions;
      static float bloom_fullness;
      static float cache_hit_findmin;
      static float cache_hit_read_c;
      static float cache_hit_read_e;
      static float cache_hit_inter;
      static float cache_hit_verify;

      static void init()
      {
            batch_reads = 0;
            findmin_reads = 0;
            readmin_counter_reads = 0;
            readmin_edge_reads = 0;
            readmin_n_sets = 0;
            readmin_vstream = 0;
            bloom_filter = 0;
            bloom_fullness = 0;
            intersect_reads = 0;
            intersect_filter = 0;
            verify_filter = 0;
            verify_reads = 0;
            homomo_trashed = 0;
            start_set = 0;
            miss_indexing = 0;
            max_collisions = 0;
            hash_collisions = 0;
            cache_hit_findmin = 0;
            cache_hit_read_c = 0;
            cache_hit_read_e = 0;
            cache_hit_inter = 0;
            cache_hit_verify = 0;
            cache_req_findmin = 0;
            cache_req_read_c = 0;
            cache_req_read_e = 0;
            cache_req_inter = 0;
            cache_req_verify = 0;
      }

      static void print(
          const unsigned char hash1_w,
          const unsigned char hash2_w)
      {
            std::ofstream debof("../../../../stats.txt", std::ofstream::app);

            unsigned long mem_accesses = batch_reads + findmin_reads +
                                         readmin_counter_reads + readmin_edge_reads +
                                         intersect_reads + verify_reads;

            unsigned int cnt, k_fun;
            cnt = COUNTER_WIDTH;
            k_fun = (1UL << K_FUNCTIONS);
            //   verify_filter -= intersect_filter;

            debof << "DEBUG STATISTICS HW1: " << (int)hash1_w << " HW2: " << (int)hash2_w
                  << " CNT: " << cnt << " K_FUN: " << k_fun << std::endl
                  << std::endl;

            debof << "\tbatch accesses:              " << std::setw(15)
                  << batch_reads << "\t" << std::setw(4)
                  << batch_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\tfindmin accesses:            " << std::setw(15)
                  << findmin_reads << "\t" << std::setw(4)
                  << findmin_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\treadmin counters accesses:   " << std::setw(15)
                  << readmin_counter_reads << "\t" << std::setw(4)
                  << readmin_counter_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\treadmin edge accesses:       " << std::setw(15)
                  << readmin_edge_reads << "\t" << std::setw(4)
                  << readmin_edge_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\tintersect accesses:          " << std::setw(15)
                  << intersect_reads << "\t" << std::setw(4)
                  << intersect_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\tverify accesses:             " << std::setw(15)
                  << verify_reads << "\t" << std::setw(4)
                  << verify_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\tTOTAL:                       " << std::setw(15)
                  << mem_accesses << "\t" << std::setw(4) << "100" << "%" << std::endl
                  << std::endl;

#if CACHE_ENABLE
            findmin_reads -= cache_hit_findmin;
            intersect_reads -= cache_hit_inter;
            verify_reads -= cache_hit_verify;
            readmin_edge_reads -= cache_hit_read_e;
            readmin_counter_reads -= cache_hit_read_c;
            mem_accesses = batch_reads + findmin_reads +
                           readmin_counter_reads + readmin_edge_reads +
                           intersect_reads + verify_reads;

            debof << "\tbatch accesses:              " << std::setw(15)
                  << batch_reads << "\t" << std::setw(4)
                  << batch_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\tfindmin accesses:            " << std::setw(15)
                  << findmin_reads << "\t" << std::setw(4)
                  << findmin_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\treadmin counters accesses:   " << std::setw(15)
                  << readmin_counter_reads << "\t" << std::setw(4)
                  << readmin_counter_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\treadmin edge accesses:       " << std::setw(15)
                  << readmin_edge_reads << "\t" << std::setw(4)
                  << readmin_edge_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\tintersect accesses:          " << std::setw(15)
                  << intersect_reads << "\t" << std::setw(4)
                  << intersect_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\tverify accesses:             " << std::setw(15)
                  << verify_reads << "\t" << std::setw(4)
                  << verify_reads * 100 / mem_accesses << "%" << std::endl;
            debof << "\tTOTAL:                       " << std::setw(15)
                  << mem_accesses << "\t" << std::setw(4) << "100" << "%" << std::endl
                  << std::endl;

            if (cache_req_findmin != 0)
            {
                  debof << "\tmean cache hit findmin:  " << cache_hit_findmin / cache_req_findmin
                        << std::endl;
                  debof << "\tcache reqs findmin:      " << cache_req_findmin << std::endl;
            }
            if (cache_req_read_e != 0)
            {
                  debof << "\tmean cache hit readmin edges: "
                        << cache_hit_read_e / cache_req_read_e << std::endl;
                  debof << "\tcache reqs readmin edges:     " << cache_req_read_e << std::endl;
            }
            if (cache_req_read_c != 0)
            {
                  debof << "\tmean cache hit readmin counters: "
                        << cache_hit_read_c / cache_req_read_c << std::endl;
                  debof << "\tcache reqs readmin counters:     " << cache_req_read_c << std::endl;
            }
            debof << "\tmean cache hit inter:   " << cache_hit_inter / cache_req_inter
                  << std::endl;
            debof << "\tcache reqs inter:       " << cache_req_inter << std::endl;
            debof << "\tmean cache hit verify:  " << cache_hit_verify / cache_req_verify
                  << std::endl;
            debof << "\tcache reqs verify:      " << cache_req_verify << std::endl;
#endif

            debof << "\tmax collisions:         " << max_collisions << std::endl;
            //        debof << "\tavg collisions:         " << hash_collisions / (float)counters << std::endl;
            debof << "\tbloom fullness:         " << bloom_fullness << std::endl;
            debof << "\tbloom filter:           " << bloom_filter << std::endl;
            debof << "\tintersect filter:       " << intersect_filter << std::endl;
            debof << "\tverify filter:          " << verify_filter << std::endl;
            debof << "\thomomorphism filter:    " << homomo_trashed << std::endl;
            debof << "\tstart_set:              " << start_set << std::endl;
            debof << "\tmiss indexing           " << miss_indexing << std::endl;
            debof << "\taverage set:            " << readmin_vstream / (float)readmin_n_sets
                  << std::endl
                  << std::endl;

            debof.close();
      }
};
