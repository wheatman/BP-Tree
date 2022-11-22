#include "ParallelTools/reducer.h"
#include <algorithm>
#include <functional>
#include <random>
#include <sys/time.h>
#include <vector>
#include <set>
#include <map>
#include <ParallelTools/parallel.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

#include <tlx/container/btree_set.hpp>
#include "tlx/container/btree_map.hpp"

#if CILK != 1
#define cilk_for for
#endif

#define CORRECTNESS 1

static long get_usecs() {
  struct timeval st;
  gettimeofday(&st, NULL);
  return st.tv_sec * 1000000 + st.tv_usec;
}

template <class T>
std::vector<T> create_random_data(size_t n, size_t max_val,
                                  std::seed_seq &seed) {

  std::mt19937_64 eng(seed); // a source of random data

  std::uniform_int_distribution<T> dist(0, max_val);
  std::vector<T> v(n);

  generate(begin(v), end(v), bind(dist, eng));
  return v;
}

void test_workload_a(int trials) {
  std::ifstream infile("/home/ubuntu/xvdf_mounted/ycsb/index-microbench/workloads/loada_unif_int.dat");
  std::string a;
  uint64_t b;
  std::vector<uint64_t> vect;
  
  while (infile >> a >> b)
  {
    vect.push_back(b);
  }
  printf("Loaded %lu elts from YCSB\n", vect.size());

  std::this_thread::sleep_for(std::chrono::nanoseconds(3000000000));
  
  printf("Stopped waiting, started timing\n");

  tlx::btree_map<uint64_t, uint64_t, std::less<uint64_t>, tlx::btree_default_traits<uint64_t, uint64_t>,
                 std::allocator<uint64_t>, true> concurrent_map;
  uint64_t start_time, end_time;
  
  for (int trial = 0; trial <= trials; trial++) {
    tlx::btree_map<uint64_t, uint64_t, std::less<uint64_t>, tlx::btree_default_traits<uint64_t, uint64_t>,
                 std::allocator<uint64_t>, true> concurrent_map;
    start_time = get_usecs();
    uint64_t sum_key_range = 0;
    uint64_t sum_val_range = 0;
    cilk_for(uint32_t i = 0; i < vect.size(); i++) {
      concurrent_map.insert({vect[i], vect[i]});
    }
    end_time = get_usecs();

    concurrent_map.map_range_length(1, vect.size(), [&sum_key_range, &sum_val_range]([[maybe_unused]] auto key, auto val) {
                  sum_key_range += key;
                  sum_val_range += val;
                });
    
    printf("did %lu ops in %lu microsecs, throughput of %f, sum key = %lu, sum_val = %lu\n", vect.size(), end_time - start_time, (float)vect.size() / (float)(end_time - start_time), sum_key_range, sum_val_range);
  }
  printf("FINISHED YCSB VERSION \n\n");
  return;

/*
  for (int trial = 0; trial < 6; trial++) {
    std::seed_seq seed{trial};
    std::vector<uint64_t> data =
        create_random_data<uint64_t>(100000000, std::numeric_limits<uint64_t>::max(), seed);

    tlx::btree_map<uint64_t, uint64_t, std::less<uint64_t>, tlx::btree_default_traits<uint64_t, uint64_t>,
                 std::allocator<uint64_t>, true> concurrent_map;
    start_time = get_usecs();
    uint64_t sum_key_range = 0;
    uint64_t sum_val_range = 0;
    cilk_for(uint32_t i = 0; i < data.size(); i++) {
      concurrent_map.insert({data[i], data[i]});
    }
    end_time = get_usecs();

    concurrent_map.map_range_length(1, vect.size(), [&sum_key_range, &sum_val_range]([[maybe_unused]] auto key, auto val) {
                  sum_key_range += key;
                  sum_val_range += val;
                });
    
    printf("did %lu ops in %lu microsecs, throughput of %f, sum key = %lu, sum_val = %lu\n", data.size(), end_time - start_time, (float)data.size() / (float)(end_time - start_time), sum_key_range, sum_val_range);
  }
  printf("FINISHED MICRO VERSION \n\n");
*/
}

int main(int argc, char *argv[]) {
  int trials = 5;
  test_workload_a(trials);
  return 0;
}