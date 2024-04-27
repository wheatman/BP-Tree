#include <algorithm>
#include <functional>
#include <random>
#include <sys/time.h>
#include <vector>
#include <set>
#include <map>
#include <ParallelTools/parallel.h>
#include "leafDS/cxxopts.hpp"
#include <fstream>
#include <thread>

#include <container/btree_set.hpp>
#include <container/btree_map.hpp>
#include "timers.hpp"

// #define CORRECTNESS 0
#define TRIALS 5
#define SLOTS 0

static long get_usecs() {
  struct timeval st;
  gettimeofday(&st, NULL);
  return st.tv_sec * 1000000 + st.tv_usec;
}

struct ThreadArgs {
    std::function<void(int, int)> func;
    int start;
    int end;
};


void* threadFunction(void* arg) {
    ThreadArgs* args = static_cast<ThreadArgs*>(arg);
    args->func(args->start, args->end);
    pthread_exit(NULL);
}


template <typename F> inline void parallel_for(size_t start, size_t end, F f) {
    const int numThreads = 48;
    pthread_t threads[numThreads];
    ThreadArgs threadArgs[numThreads];
    int per_thread = (end - start)/numThreads;

    // Create the threads and start executing the lambda function
    for (int i = 0; i < numThreads; i++) {
        threadArgs[i].func = [&f](int arg1, int arg2) {
            for (int k = arg1 ; k < arg2; k++) {
                f(k);
            }
        };

        threadArgs[i].start = start + (i * per_thread);
        if (i == numThreads - 1) {
          threadArgs[i].end = end;
        } else {
          threadArgs[i].end = start + ((i+1) * per_thread);
        }
        int result = pthread_create(&threads[i], NULL, threadFunction, &threadArgs[i]);

        if (result != 0) {
            std::cerr << "Failed to create thread " << i << std::endl;
            exit(-1);
        }
    }

    // Wait for the threads to finish
    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }

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

template <class T>
std::tuple<bool, uint64_t, uint64_t, uint64_t, uint64_t>
test_concurrent_btreeset(uint64_t max_size, std::seed_seq &seed) {
  // std::vector<T> data =
  //     create_random_data<T>(max_size, 10000, seed);
  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);

  for(uint32_t i = 0; i < max_size; i++) {
    data[i]++; // no zeroes
  }
  uint64_t start, end;
#if CORRECTNESS
  std::set<T> serial_set;
  std::vector<T> inserted_elts;
  std::set<T> duplicated_elts;

  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    // if (std::count(serial_set.begin(), serial_set.end(), data[i])) {
    //   printf("inserting twice %lu\n", data[i]);
    //   duplicated_elts.insert(data[i]);
    // }
    serial_set.insert(data[i]);
  }
  printf("%lu unique elements\n", serial_set.size());
  end = get_usecs();

  for (auto e: serial_set) {
    inserted_elts.push_back(e);
  }
  // int64_t serial_time = end - start;
  // printf("inserted all the data serially in %lu\n", end - start);
#endif

  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, false>
      serial_test_set;
  start = get_usecs();
  for(uint32_t i = 0; i < max_size; i++) {
    serial_test_set.insert(data[i]);
  }
  end = get_usecs();
  uint64_t serial_time = end - start;
  printf("\tinserted %d elts serially in %lu\n", max_size, end - start);

  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, true>
      concurrent_set;
  start = get_usecs();
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_set.insert(data[i]);
  }
  end = get_usecs();
  uint64_t parallel_time = end - start;
  printf("\tinserted %d elts concurrently in %lu\n", max_size, end - start);

#if CORRECTNESS
  bool wrong = false;
  uint64_t correct_sum = 0;
  for(auto e : serial_set) {
    // might need to change if you get values because exists takes in a key
    if(!concurrent_set.exists(e)) {
      printf("insertion, didn't find %lu\n", e);
      wrong = true;
    }
    correct_sum += e;
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  }

  /*
  // test sum
  uint64_t serial_test_sum = serial_test_set.psum();
  ASSERT("serial sum got %lu, should be %lu\n", serial_test_sum, correct_sum);

  uint64_t concurrent_test_sum = concurrent_set.psum();
  ASSERT("concurrent sum got %lu, should be %lu\n", concurrent_sum, correct_sum);
  */
#endif

  std::vector<uint64_t> indxs_to_remove =
    create_random_data<uint64_t>(max_size / 2, data.size(), seed);

#if CORRECTNESS
  std::vector<T> elems_removed;
  for (uint32_t i = 0; i < indxs_to_remove.size(); i++) {
    serial_set.erase(data[indxs_to_remove[i]]);
    elems_removed.push_back(data[indxs_to_remove[i]]);
    if (std::count(duplicated_elts.begin(), duplicated_elts.end(), data[indxs_to_remove[i]])) {
      printf("trying to delete a duplicated elt %lu\n", data[indxs_to_remove[i]]);
    }
  }
#endif

  uint64_t serial_remove_start = get_usecs();
  for (uint32_t i = 0; i < indxs_to_remove.size(); i++) {
    serial_test_set.erase(data[indxs_to_remove[i]]);
  }
  uint64_t serial_remove_end = get_usecs();
  uint64_t serial_remove_time = serial_remove_end - serial_remove_start;

  printf("removed half the data serially in %lu\n",
         serial_remove_end - serial_remove_start);
         
  uint64_t parallel_remove_start = get_usecs();
  cilk_for(uint32_t i = 0; i < indxs_to_remove.size(); i++) {
    concurrent_set.erase(data[indxs_to_remove[i]]);
  }

  uint64_t parallel_remove_end = get_usecs();
  uint64_t parallel_remove_time = parallel_remove_end - parallel_remove_start;


  printf("removed half the data in parallel in %lu\n",
         parallel_remove_end - parallel_remove_start);

#if CORRECTNESS
  // // Deletion correctness check

  // check that everything in serial set (not erased) is still in btree
  for (auto e : serial_set) {
    if (!serial_test_set.exists(e)) {
      printf("serial btree, incorrectly deleted %lu \n", e);
      wrong = true;
    }
  }
  // check that everything erased from serial set is not in btree
  for (int i = 0; i < elems_removed.size(); i++) {
    if (serial_test_set.exists(elems_removed[i])) {
      printf("serial btree, didn't delete %lu, i = %u \n", elems_removed[i], i);
      wrong = true;
    }
  }

  if (wrong) {
    return {false, 0, 0, 0, 0};
  }

  for (auto e : serial_set) {
    if (!concurrent_set.exists(e)) {
      printf("concurrent btree, incorrectly deleted %lu \n", e);
      wrong = true;
    }
  }
  // check that everything erased from serial set is not in btree
  for (int i = 0; i < elems_removed.size(); i++) {
    if (concurrent_set.exists(elems_removed[i])) {
      printf("concurrent btree, didn't delete %lu, i = %u \n", elems_removed[i], i);
      wrong = true;
    }
  }

  if (wrong) {
    return {false, 0, 0, 0, 0};
  }
#endif

  // return {true, 0, 0, 0, 0};

  // if (serial_set.size() != concurrent_set.size()) {
  //   printf("the sizes of concurrent set and correct set don't match, got %lu, expected %lu\n",
  //          concurrent_set.size(), serial_set.size());
  //   return {false, 0, 0, 0, 0};
  // }

  // for(auto e : indxs_to_remove) {
  //   // might need to change if you get values because exists takes in a key
  //   if(concurrent_set.exists(data[e])) {
  //     printf("concurrent, didn't delete %lu\n", data[e]);
  //     wrong = true;
  //   }
  // }
  // if (wrong) {
  //   return {false, 0, 0, 0, 0};
  // }
// #endif

  return {true, serial_time, parallel_time, serial_remove_time, parallel_remove_time};
}

template <class T, uint32_t internal_bytes>
void test_concurrent_sum_correctness(uint64_t max_size, std::seed_seq &seed) {
  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);
  std::set<T> serial_set;

  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
                    std::allocator<T>, false> serial_test_set;

  for (uint32_t i = 0; i < max_size; i++) {
    serial_set.insert(data[i]);
    serial_test_set.insert(data[i]);
  }

  uint64_t correct_sum = 0;
  for(auto e : serial_set) {
          correct_sum += e;
  }
  auto serial_sum = serial_test_set.psum_with_map();
  printf("serial btree sum with map got %lu, should be %lu\n", serial_sum, correct_sum);
  assert(serial_sum == correct_sum);
  serial_sum = serial_test_set.psum_with_subtract();
  printf("serial btree sum with subtract got %lu, should be %lu\n", serial_sum, correct_sum);
  assert(serial_sum == correct_sum);

  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, true>
      concurrent_set;
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_set.insert(data[i]);
  }
  auto concurrent_sum = concurrent_set.psum_with_map();
  printf("\tconcurrent btree sum with map got %lu, should be %lu\n", concurrent_sum, correct_sum);

  assert(concurrent_sum == correct_sum);
  concurrent_sum = concurrent_set.psum_with_subtract();
  assert(concurrent_sum == correct_sum);
  printf("\tconcurrent btree sum with subtract got %lu, should be %lu\n", concurrent_sum, correct_sum);

  //printf("concurrent btree sum got %lu, should be %lu\n", concurrent_sum, correct_sum);
}


template <class T, uint32_t internal_bytes>
void test_concurrent_sum_time(uint64_t max_size, std::seed_seq &seed, int trials) {
  std::vector<T> data =
      create_random_data<T>(2*max_size, std::numeric_limits<T>::max(), seed);

  std::vector<uint64_t> insert_times(trials);
  std::vector<uint64_t> with_map_times(trials);
  std::vector<uint64_t> with_subtract_times(trials);
  uint64_t start, end, map_time, subtract_time, insert_time;
  for(int i = 0; i < trials + 1; i++) {
	  tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
                  std::allocator<T>, true> concurrent_map;

    // TIME INSERTS
    start = get_usecs();
    parallel_for(0, 2*max_size, [&](const uint32_t &i) {
      concurrent_map.insert({data[i], 2*data[i]});
    });
    end = get_usecs();
    printf("\tDone inserting %lu elts in %lu\n",max_size, end - start);

	  start = get_usecs();
	  auto concurrent_sum = concurrent_map.psum_with_map();
	  end = get_usecs();
	  map_time = end - start;
	  
	  printf("\tconcurrent btree sum with map got %lu in %lu\n", concurrent_sum, map_time);

		  if(i > 0) {
			  insert_times[i-1] = insert_time;
			  with_map_times[i-1] = map_time;
		  }

  }
  for(int i = 0; i < trials + 1; i++) {
	  tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
                  std::allocator<T>, true> concurrent_map;

    // TIME INSERTS
    start = get_usecs();
    parallel_for(0, 2*max_size, [&](const uint32_t &i) {
      concurrent_map.insert({data[i], 2*data[i]});
    });
    end = get_usecs();
    printf("\tDone inserting %lu elts in %lu\n",max_size, end - start);

	  start = get_usecs();
	  auto concurrent_sum = concurrent_map.psum_with_subtract();
	  end = get_usecs();
	  subtract_time = end - start;
	  printf("\tconcurrent btree sum with subtract got %lu in %lu\n", concurrent_sum, subtract_time);
		  if(i > 0) {
			  with_subtract_times[i-1] = subtract_time;
		  }
  }
  std::sort(insert_times.begin(), insert_times.end());
  std::sort(with_map_times.begin(), with_map_times.end());
  std::sort(with_subtract_times.begin(), with_subtract_times.end());
  printf("insert time = %lu\n", insert_times[trials / 2]);
  printf("with map time = %lu, with subtract time = %lu\n", with_map_times[trials / 2], with_subtract_times[trials / 2]);

}


template <class T, uint32_t internal_bytes>
void test_concurrent_btreeset_scalability(uint64_t max_size, uint64_t NUM_QUERIES,
                                          std::seed_seq &seed, int trials) {
  std::mt19937_64 eng(seed);
  std::vector<T> data = create_random_data<T>(max_size, std::numeric_limits<T>::max() - (trials + 1), seed);
  
  std::seed_seq query_seed{1};
  std::seed_seq query_seed_2{2};
  std::vector<uint64_t> range_query_start_idxs =
        create_random_data<uint64_t>(NUM_QUERIES, data.size() - 1, query_seed);

  uint64_t MAX_QUERY_SIZE = 100000;

  printf("done generating the data\n");
  uint64_t start_time, end_time;
  uint64_t parallel_time = 0;
#if ENABLE_TRACE_TIMER == 1
  std::vector<uint64_t> concurrent_latencies(trials * max_size);
  std::vector<uint64_t> concurrent_find_latencies(trials * max_size);
  std::vector<uint64_t> concurrent_sorted_range_latencies(trials * max_size);
  std::vector<uint64_t> concurrent_unsorted_range_latencies(trials * max_size);
#endif
  for (int j = 0; j < trials; j++) {
    tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
                    std::allocator<T>, true>
        concurrent_map;
    start_time = get_usecs();
    cilk_for(uint32_t i = 0; i < max_size; i++) {
      uint64_t start_indv = timer::get_time();
      concurrent_map.insert({data[i] + j, 2*(data[i] + j)});
      uint64_t end_indv = timer::get_time();
  #if ENABLE_TRACE_TIMER == 1
      concurrent_latencies[j * max_size + i] = end_indv - start_indv;
  #endif
    }
    end_time = get_usecs();
    parallel_time += end_time - start_time;
    printf("concurrent insert trial took  %lu\n", end_time - start_time);
    
    // TIME POINT FINDS
    std::vector<bool> found_count(NUM_QUERIES);
    start_time = get_usecs();
    cilk_for(uint32_t i = 0; i < NUM_QUERIES; i++) {
      uint64_t start_indv = timer::get_time();
      found_count[i] = concurrent_map.exists(data[range_query_start_idxs[i]]);
      uint64_t end_indv = timer::get_time();
    #if ENABLE_TRACE_TIMER == 1
      concurrent_find_latencies[j * max_size + i] = end_indv - start_indv;
    #endif
    }
    end_time = get_usecs();
    int count_found = 0;
    for (auto e : found_count) {
      count_found += e ? 1 : 0;
    }
    printf("\tDone finding %lu elts in %lu, count = %d \n",NUM_QUERIES, end_time - start_time, count_found);

    // TIME RANGE QUERIES
    std::vector<uint64_t> range_query_lengths =
      create_random_data<uint64_t>(NUM_QUERIES, MAX_QUERY_SIZE - 1, query_seed_2);

    std::vector<T> concurrent_range_query_length_maxs(NUM_QUERIES);
    std::vector<uint64_t> concurrent_range_query_length_counts(NUM_QUERIES);
    std::vector<T> concurrent_range_query_length_key_sums(NUM_QUERIES);
    std::vector<T> concurrent_range_query_length_val_sums(NUM_QUERIES);

    // SORTED RANGE QUERY
    cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
      T start;
      start = data[range_query_start_idxs[i]];

      uint64_t num_in_range = 0;
      T max_key_in_range = start;
      concurrent_map.map_range_length(start, range_query_lengths[i], [&num_in_range, &max_key_in_range]([[maybe_unused]] auto key, auto val) {
                num_in_range++;
                if (key > max_key_in_range) {
                  max_key_in_range = key;
                }
              });
      concurrent_range_query_length_counts[i] = num_in_range;
      concurrent_range_query_length_maxs[i] = max_key_in_range;
    }

    start_time = get_usecs();
    cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
      T start;
      start = data[range_query_start_idxs[i]];
      T sum_key_range = 0;
      T sum_val_range = 0;
      uint64_t start_indv = timer::get_time();
      concurrent_map.map_range_length(start, range_query_lengths[i], [&sum_key_range, &sum_val_range]([[maybe_unused]] auto key, auto val) {
                sum_key_range += key;
                sum_val_range += val;
              });
      concurrent_range_query_length_key_sums[i] = sum_key_range;
      concurrent_range_query_length_val_sums[i] = sum_val_range;
      uint64_t end_indv = timer::get_time();
    #if ENABLE_TRACE_TIMER == 1
      concurrent_sorted_range_latencies[j * max_size + i] = end_indv - start_indv;
    #endif
    }
    end_time = get_usecs();
    printf("\t did sorted range queries with max len %lu concurrently in %lu\n", MAX_QUERY_SIZE, end_time - start_time);
    T sum_all_keys = 0;
    for (auto e: concurrent_range_query_length_key_sums) {
      sum_all_keys += e;
    }
    T sum_all_vals = 0;
    for (auto e: concurrent_range_query_length_val_sums) {
      sum_all_vals += e;
    }
    if (sum_all_keys * 2 != sum_all_vals) {
      printf("\t\t wrong, sum keys * 2 not equal to sum vals\n");
    }
    printf("\t\t sum keys %lu sum vals %lu\n", sum_all_keys, sum_all_vals);

    std::vector<T> concurrent_range_query_key_sums(NUM_QUERIES);
    std::vector<T> concurrent_range_query_val_sums(NUM_QUERIES);

    start_time = get_usecs();
    cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
      T start, end;
      start = data[range_query_start_idxs[i]];
      end = concurrent_range_query_length_maxs[i];
      if (range_query_lengths[i] != 0) {
        end++;
      }

      T sum_key_range = 0;
      T sum_val_range = 0;
      uint64_t start_indv = timer::get_time();
      concurrent_map.map_range(start, end, [&sum_key_range, &sum_val_range]([[maybe_unused]] auto key, auto val) {
                sum_key_range += key;
                sum_val_range += val;
              });
      concurrent_range_query_key_sums[i] = sum_key_range;
      concurrent_range_query_val_sums[i] = sum_val_range;
      uint64_t end_indv = timer::get_time();
    #if ENABLE_TRACE_TIMER == 1
      concurrent_unsorted_range_latencies[j * max_size + i] = end_indv - start_indv;
    #endif
    }
    end_time = get_usecs();
    printf("\t did unsorted range queries with max len %lu concurrently in %lu\n", MAX_QUERY_SIZE, end_time - start_time);
    sum_all_keys = 0;
    for (auto e: concurrent_range_query_key_sums) {
      sum_all_keys += e;
    }
    sum_all_vals = 0;
    for (auto e: concurrent_range_query_val_sums) {
      sum_all_vals += e;
    }
    if (sum_all_keys * 2 != sum_all_vals) {
      printf("\t\t wrong, sum keys * 2 not equal to sum vals\n");
    }
    printf("\t\t sum keys %lu sum vals %lu\n", sum_all_keys, sum_all_vals);

  }

  // printf("Parallel insert time / trials %lu\n", parallel_time / trials);
#if ENABLE_TRACE_TIMER == 1
  std::sort(concurrent_latencies.begin(), concurrent_latencies.end());
  std::sort(concurrent_find_latencies.begin(), concurrent_find_latencies.end());
  std::sort(concurrent_sorted_range_latencies.begin(), concurrent_sorted_range_latencies.end());
  std::sort(concurrent_unsorted_range_latencies.begin(), concurrent_unsorted_range_latencies.end());
  printf("concurrent inserts: 50%% = %lu, 90%% = %lu, 99%% = %lu, 99.9%% = %lu, max = "
         "%lu\n",
         concurrent_latencies[concurrent_latencies.size() / 2],
         concurrent_latencies[concurrent_latencies.size() * 9 / 10],
         concurrent_latencies[concurrent_latencies.size() * 99 / 100],
         concurrent_latencies[concurrent_latencies.size() * 999 / 1000],
         concurrent_latencies[concurrent_latencies.size() - 1]);
  printf("concurrent finds: 50%% = %lu, 90%% = %lu, 99%% = %lu, 99.9%% = %lu, max = "
         "%lu\n",
         concurrent_find_latencies[concurrent_find_latencies.size() / 2],
         concurrent_find_latencies[concurrent_find_latencies.size() * 9 / 10],
         concurrent_find_latencies[concurrent_find_latencies.size() * 99 / 100],
         concurrent_find_latencies[concurrent_find_latencies.size() * 999 / 1000],
         concurrent_find_latencies[concurrent_find_latencies.size() - 1]);
  printf("concurrent sorted range: 50%% = %lu, 90%% = %lu, 99%% = %lu, 99.9%% = %lu, max = "
         "%lu\n",
         concurrent_sorted_range_latencies[concurrent_sorted_range_latencies.size() / 2],
         concurrent_sorted_range_latencies[concurrent_sorted_range_latencies.size() * 9 / 10],
         concurrent_sorted_range_latencies[concurrent_sorted_range_latencies.size() * 99 / 100],
         concurrent_sorted_range_latencies[concurrent_sorted_range_latencies.size() * 999 / 1000],
         concurrent_sorted_range_latencies[concurrent_sorted_range_latencies.size() - 1]);
  printf("concurrent unsorted range: 50%% = %lu, 90%% = %lu, 99%% = %lu, 99.9%% = %lu, max = "
         "%lu\n",
         concurrent_unsorted_range_latencies[concurrent_unsorted_range_latencies.size() / 2],
         concurrent_unsorted_range_latencies[concurrent_unsorted_range_latencies.size() * 9 / 10],
         concurrent_unsorted_range_latencies[concurrent_unsorted_range_latencies.size() * 99 / 100],
         concurrent_unsorted_range_latencies[concurrent_unsorted_range_latencies.size() * 999 / 1000],
         concurrent_unsorted_range_latencies[concurrent_unsorted_range_latencies.size() - 1]);
#endif
}

template <class T>
std::tuple<bool, uint64_t, uint64_t, uint64_t, uint64_t>
test_concurrent_btreemap(uint64_t max_size, std::seed_seq &seed) {
  // std::vector<T> data =
  //     create_random_data<T>(max_size, 10000, seed);
  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);

  for(uint32_t i = 0; i < max_size; i++) {
    data[i]++; // no zeroes
  }
  uint64_t start, end;
#if CORRECTNESS
  std::map<T, T> serial_map;
  std::vector<std::pair<T, T>> inserted_elts;

  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    serial_map.insert({data[i], 2*data[i]});
  }
  printf("%lu unique elements\n", serial_map.size());
  end = get_usecs();

  for (auto e: serial_map) {
    inserted_elts.push_back(e);
  }
  // int64_t serial_time = end - start;
  // printf("inserted all the data serially in %lu\n", end - start);
#endif

  tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, false>
      serial_test_map;
  start = get_usecs();
  for(uint32_t i = 0; i < max_size; i++) {
    serial_test_map.insert({data[i], 2*data[i]});
  }
  end = get_usecs();
  uint64_t serial_time = end - start;
  printf("\tinserted %lu elts serially in %lu\n", max_size, end - start);

  tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, true>
      concurrent_map;
  start = get_usecs();
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_map.insert({data[i], 2*data[i]});
  }
  end = get_usecs();
  uint64_t parallel_time = end - start;
  printf("\tinserted %lu elts concurrently in %lu\n", max_size, end - start);

#if CORRECTNESS
  bool wrong = false;
  uint64_t correct_sum = 0;
  // check serial
  for(auto e : serial_map) {
    // might need to change if you get values because exists takes in a key
    if(!serial_test_map.exists(e.first)) {
      printf("insertion, didn't find key %lu\n", e.first);
      wrong = true;
    }
    if (serial_test_map.value(e.first) != e.second) {
      printf("got the wrong value for key %lu, got %lu, expected %lu ", e.first, serial_test_map.value(e.first), e.second);
    }
    correct_sum += e.first;
  }
  // check concurrent
  for(auto e : serial_map) {
    // might need to change if you get values because exists takes in a key
    if(!concurrent_map.exists(e.first)) {
      printf("insertion, didn't find key %lu\n", e.first);
      wrong = true;
    }
    if (concurrent_map.value(e.first) != e.second) {
      printf("got the wrong value for key %lu, got %lu, expected %lu ", e.first, concurrent_map.value(e.first), e.second);
    }
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  }
#endif

  std::vector<uint64_t> indxs_to_remove =
    create_random_data<uint64_t>(max_size / 2, data.size(), seed);

#if CORRECTNESS
  std::vector<T> elems_removed;
  for (uint32_t i = 0; i < indxs_to_remove.size(); i++) {
    serial_map.erase(data[indxs_to_remove[i]]);
    elems_removed.push_back(data[indxs_to_remove[i]]);
  }
#endif

  uint64_t serial_remove_start = get_usecs();
  for (uint32_t i = 0; i < indxs_to_remove.size(); i++) {
    serial_test_map.erase(data[indxs_to_remove[i]]);
  }
  uint64_t serial_remove_end = get_usecs();
  uint64_t serial_remove_time = serial_remove_end - serial_remove_start;

  printf("removed half the data serially in %lu\n",
         serial_remove_end - serial_remove_start);
         
  uint64_t parallel_remove_start = get_usecs();
  cilk_for(uint32_t i = 0; i < indxs_to_remove.size(); i++) {
    concurrent_map.erase(data[indxs_to_remove[i]]);
  }

  uint64_t parallel_remove_end = get_usecs();
  uint64_t parallel_remove_time = parallel_remove_end - parallel_remove_start;


  printf("removed half the data in parallel in %lu\n",
         parallel_remove_end - parallel_remove_start);

#if CORRECTNESS
  // // Deletion correctness check

  // check that everything in serial set (not erased) is still in btree
  for (auto e : serial_map) {
    if (!serial_test_map.exists(e.first)) {
      printf("serial btree, incorrectly deleted %lu \n", e);
      wrong = true;
    }
  }
  // check that everything erased from serial set is not in btree
  for (int i = 0; i < elems_removed.size(); i++) {
    if (serial_test_map.exists(elems_removed[i])) {
      printf("serial btree, didn't delete %lu, i = %u \n", elems_removed[i], i);
      wrong = true;
    }
  }

  if (wrong) {
    return {false, 0, 0, 0, 0};
  }

  for (auto e : serial_map) {
    if (!concurrent_map.exists(e.first)) {
      printf("concurrent btree, incorrectly deleted %lu \n", e);
      wrong = true;
    }
  }
  // check that everything erased from serial set is not in btree
  for (int i = 0; i < elems_removed.size(); i++) {
    if (concurrent_map.exists(elems_removed[i])) {
      printf("concurrent btree, didn't delete %lu, i = %u \n", elems_removed[i], i);
      wrong = true;
    }
  }

  if (wrong) {
    return {false, 0, 0, 0, 0};
  }
#endif

  return {true, serial_time, parallel_time, serial_remove_time, parallel_remove_time};
}

template <class T>
std::tuple<bool, uint64_t, uint64_t, uint64_t, uint64_t>
test_concurrent_range_query(uint64_t max_size, std::seed_seq &seed) {
  uint64_t NUM_QUERIES = 100000;
  uint64_t MAX_QUERY_SIZE = 100;

  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);
  // std::vector<T> data =
  //     create_random_data<T>(max_size, 10000, seed);

  // // do 1/2 * max_size range queries
  // std::vector<T> range_queries = 
  //     create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);

  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, false> serial_set;
  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, true> concurrent_set;
  std::set<T> checker_set;

  for (uint32_t i = 0; i < max_size; i++) {
    // serial_set.insert(data[i]);
    checker_set.insert(data[i]);
  }

  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_set.insert(data[i]);
  }

  std::vector<T> checker_sorted;
  for (auto e: checker_set) {
    checker_sorted.push_back(e);
  }
  std::sort(checker_sorted.begin(), checker_sorted.end());

  printf("Done inserting %lu elts\n",max_size);

  uint64_t start_time, end_time, serial_time, concurrent_time;

  std::seed_seq query_seed{1};
  std::seed_seq query_seed_2{2};

  std::vector<uint64_t> range_query_start_idxs =
      create_random_data<uint64_t>(NUM_QUERIES, checker_sorted.size(), query_seed);
  
  std::vector<uint64_t> range_query_lengths =
      create_random_data<uint64_t>(NUM_QUERIES, MAX_QUERY_SIZE, query_seed_2);
  
  // std::vector<uint64_t> range_query_end_idxs =
  //     create_random_data<uint64_t>(checker_sorted.size() / 2, checker_sorted.size(), seed);

  std::vector<T> correct_range_query_sums(NUM_QUERIES);
  std::vector<uint64_t> correct_range_query_counts(NUM_QUERIES);

  std::vector<T> serial_range_query_sums(NUM_QUERIES);
  std::vector<uint64_t> serial_range_query_counts(NUM_QUERIES);

  std::vector<T> concurrent_range_query_sums(NUM_QUERIES);
  std::vector<uint64_t> concurrent_range_query_counts(NUM_QUERIES);

  T start, end;
  bool wrong = false;

#if CORRECTNESS
  // get correct range sums
  cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
    start = checker_sorted[range_query_start_idxs[i]];
    end = checker_sorted[std::min(range_query_start_idxs[i] + range_query_lengths[i], checker_sorted.size() - 1)];

    correct_range_query_counts[i] = 0;
    correct_range_query_sums[i] = 0;

    for (int j = 0; checker_sorted[j] < end; j++) {
      if (checker_sorted[j] < start) {
        continue;
      }
      correct_range_query_counts[i] += 1;
      correct_range_query_sums[i] += checker_sorted[j];
    }
  }
  printf("\t did %lu correct range queries with max query size %lu\n", 
        NUM_QUERIES,
        MAX_QUERY_SIZE);
#endif

  /*
  start_time = get_usecs();
  // serial btree range sums
  for (uint32_t i = 0; i < NUM_QUERIES; i++) {
    start = checker_sorted[range_query_start_idxs[i]];
    end = checker_sorted[std::min(range_query_start_idxs[i] + range_query_lengths[i], checker_sorted.size() - 1)];

    uint64_t num_in_range = 0;
    uint64_t sum_in_range = 0;
    serial_set.map_range(start, end, [&num_in_range, &sum_in_range]([[maybe_unused]] auto el) {
              num_in_range += 1;
              sum_in_range += el;
            });
    serial_range_query_counts[i] = num_in_range;
    serial_range_query_sums[i] = sum_in_range;
  }
  end_time = get_usecs();
  serial_time = end_time - start_time;
  printf("\t did %lu range queries serially in %lu\n", NUM_QUERIES, serial_time);

  // correctness check of serial 
  for (size_t i = 0; i < correct_range_query_sums.size(); i++) {
    if (correct_range_query_sums[i] != serial_range_query_sums[i]) {
      printf("wrong serial range query sum, expected %lu, got %lu\n", correct_range_query_sums[i], serial_range_query_sums[i]);
      wrong = true;
    }
    if (correct_range_query_counts[i] != serial_range_query_counts[i]) {
      printf("wrong serial range query count, expected %lu, got %lu\n", correct_range_query_counts[i], serial_range_query_counts[i]);
      wrong = true;
    }
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  } 
  */

  start_time = get_usecs();
  // concurrent btree range sums
  cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
    start = checker_sorted[range_query_start_idxs[i]];
    end = checker_sorted[std::min(range_query_start_idxs[i] + range_query_lengths[i], checker_sorted.size() - 1)];

    uint64_t num_in_range = 0;
    uint64_t sum_in_range = 0;
    concurrent_set.map_range(start, end, [&num_in_range, &sum_in_range]([[maybe_unused]] auto el) {
              num_in_range += 1;
              sum_in_range += el;
            });
    concurrent_range_query_counts[i] = num_in_range;
    concurrent_range_query_sums[i] = sum_in_range;
  }
  end_time = get_usecs();
  concurrent_time = end_time - start_time;
  printf("\t did %lu range queries concurrently in %lu\n", NUM_QUERIES, concurrent_time);

#if CORRECTNESS
  // correctness check of concurrent 
  wrong = false;
  for (size_t i = 0; i < correct_range_query_sums.size(); i++) {
    if (correct_range_query_sums[i] != concurrent_range_query_sums[i]) {
      printf("wrong concurrent range query sum, expected %lu, got %lu\n", correct_range_query_sums[i], concurrent_range_query_sums[i]);
      wrong = true;
    }
    if (correct_range_query_counts[i] != concurrent_range_query_counts[i]) {
      printf("wrong concurrent range query count, expected %lu, got %lu\n", correct_range_query_counts[i], concurrent_range_query_counts[i]);
      wrong = true;
    }
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  } 
#endif
  // return {true, 0, 0, 0, 0};
  // */

  // *** testing map_range_length ***

  std::vector<T> serial_range_query_length_sums(NUM_QUERIES);
  std::vector<uint64_t> serial_range_query_length_counts(NUM_QUERIES);

  std::vector<T> concurrent_range_query_length_sums(NUM_QUERIES);
  std::vector<uint64_t> concurrent_range_query_length_counts(NUM_QUERIES);

  /*
  start_time = get_usecs();
  // serial btree range sums
  for (uint32_t i = 0; i < NUM_QUERIES; i++) {
    start = checker_sorted[range_query_start_idxs[i]];
    // end = checker_sorted[std::min(range_query_start_idxs[i] + range_query_lengths[i], checker_sorted.size() - 1)];

    uint64_t num_in_range = 0;
    uint64_t sum_in_range = 0;
    serial_set.map_range_length(start, correct_range_query_counts[i], [&num_in_range, &sum_in_range]([[maybe_unused]] auto el) {
              num_in_range += 1;
              sum_in_range += el;
            });
    // printf("hi ");
    serial_range_query_length_counts[i] = num_in_range;
    serial_range_query_length_sums[i] = sum_in_range;
  }
  end_time = get_usecs();
  serial_time = end_time - start_time;
  printf("\t did %lu range queries by length serially in %lu\n", NUM_QUERIES, serial_time);

  // correctness check of serial 
  wrong = false;
  for (size_t i = 0; i < correct_range_query_sums.size(); i++) {
    if (correct_range_query_sums[i] != serial_range_query_length_sums[i]) {
      printf("wrong serial range query length sum, expected %lu, got %lu\n", correct_range_query_sums[i], serial_range_query_length_sums[i]);
      wrong = true;
    }
    if (correct_range_query_counts[i] != serial_range_query_length_counts[i]) {
      printf("wrong serial range query length count, expected %lu, got %lu\n", correct_range_query_counts[i], serial_range_query_length_counts[i]);
      wrong = true;
    }
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  } 
  */

  start_time = get_usecs();
  // concurrent btree range sums
  cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
    start = checker_sorted[range_query_start_idxs[i]];
    // end = checker_sorted[std::min(range_query_start_idxs[i] + range_query_lengths[i], checker_sorted.size() - 1)];

    uint64_t num_in_range = 0;
    uint64_t sum_in_range = 0;

#if CORRECTNESS
    concurrent_set.map_range_length(start, correct_range_query_counts[i], [&num_in_range, &sum_in_range]([[maybe_unused]] auto el) {
              num_in_range += 1;
              sum_in_range += el;
            });
#else
    concurrent_set.map_range_length(start, concurrent_range_query_counts[i], [&num_in_range, &sum_in_range]([[maybe_unused]] auto el) {
              num_in_range += 1;
              sum_in_range += el;
            });
#endif
    concurrent_range_query_length_counts[i] = num_in_range;
    concurrent_range_query_length_sums[i] = sum_in_range;
  }
  end_time = get_usecs();
  concurrent_time = end_time - start_time;
  printf("\t did %lu range queries by length concurrently in %lu\n", NUM_QUERIES, concurrent_time);

#if CORRECTNESS
  // correctness check of concurrent 
  wrong = false;
  for (size_t i = 0; i < correct_range_query_sums.size(); i++) {
    if (correct_range_query_sums[i] != concurrent_range_query_length_sums[i]) {
      printf("wrong concurrent range query length sum, expected %lu, got %lu\n", correct_range_query_sums[i], concurrent_range_query_length_sums[i]);
      wrong = true;
    }
    if (correct_range_query_counts[i] != concurrent_range_query_length_counts[i]) {
      printf("wrong concurrent range query length count, expected %lu, got %lu\n", correct_range_query_counts[i], concurrent_range_query_length_counts[i]);
      wrong = true;
    }
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  } 
#endif
  return {true, 0, 0, 0, 0};
}

template <class T>
std::tuple<bool, uint64_t, uint64_t, uint64_t, uint64_t>
test_concurrent_range_query_map(uint64_t max_size, std::seed_seq &seed) {
  uint64_t NUM_QUERIES = 100000;
  uint64_t MAX_QUERY_SIZE = 100;

  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);
  // std::vector<T> data =
  //     create_random_data<T>(max_size, 10000, seed);

  // // do 1/2 * max_size range queries
  // std::vector<T> range_queries = 
  //     create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);

  tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, false> serial_map;
  tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, true> concurrent_map;
  std::map<T, T> checker_map;
  uint64_t start_time, end_time, serial_time, concurrent_time;

  for (uint32_t i = 0; i < max_size; i++) {
    // serial_set.insert({data[i], 2*data[i]});
    checker_map.insert({data[i], 2*data[i]});
  }

  start_time = get_usecs();
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_map.insert({data[i], 2*data[i]});
  }
  end_time = get_usecs();

  std::vector<std::tuple<T, T>> checker_sorted;
  for (auto e: checker_map) {
    checker_sorted.push_back(e);
  }
  std::sort(checker_sorted.begin(), checker_sorted.end());

  printf("Done inserting %lu elts in %lu\n",max_size, end_time - start_time);


  std::seed_seq query_seed{1};
  std::seed_seq query_seed_2{2};

  std::vector<uint64_t> range_query_start_idxs =
      create_random_data<uint64_t>(NUM_QUERIES, checker_sorted.size(), query_seed);
  
  std::vector<uint64_t> range_query_lengths =
      create_random_data<uint64_t>(NUM_QUERIES, MAX_QUERY_SIZE, query_seed_2);
  
  // std::vector<uint64_t> range_query_end_idxs =
  //     create_random_data<uint64_t>(checker_sorted.size() / 2, checker_sorted.size(), seed);

  std::vector<T> correct_range_query_sums(NUM_QUERIES);
  std::vector<uint64_t> correct_range_query_counts(NUM_QUERIES);

  std::vector<T> serial_range_query_sums(NUM_QUERIES);
  std::vector<uint64_t> serial_range_query_counts(NUM_QUERIES);

  std::vector<T> concurrent_range_query_sums(NUM_QUERIES);
  std::vector<uint64_t> concurrent_range_query_counts(NUM_QUERIES);

  T start, end;
  bool wrong = false;

#if CORRECTNESS
  // get correct range sums
  cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
    start = std::get<0>(checker_sorted[range_query_start_idxs[i]]);
    end = std::get<0>(checker_sorted[std::min(range_query_start_idxs[i] + range_query_lengths[i], checker_sorted.size() - 1)]);

    correct_range_query_counts[i] = 0;
    correct_range_query_sums[i] = 0;

    for (int j = 0; std::get<0>(checker_sorted[j]) < end; j++) {
      if (std::get<0>(checker_sorted[j]) < start) {
        continue;
      }
      correct_range_query_counts[i] += 1;
      correct_range_query_sums[i] += std::get<1>(checker_sorted[j]);
    }
  }
  printf("\t did %lu correct range queries with max query size %lu\n", 
        NUM_QUERIES,
        MAX_QUERY_SIZE);
#endif

  /*
  start_time = get_usecs();
  // serial btree range sums
  for (uint32_t i = 0; i < NUM_QUERIES; i++) {
    start = std::get<0>(checker_sorted[range_query_start_idxs[i]]);
    end = std::get<0>(checker_sorted[std::min(range_query_start_idxs[i] + range_query_lengths[i], checker_sorted.size() - 1)]);

    uint64_t num_in_range = 0;
    uint64_t sum_in_range = 0;
    serial_map.map_range(start, end, [&num_in_range, &sum_in_range]([[maybe_unused]] auto key, auto val) {
              num_in_range += 1;
              sum_in_range += val;
            });
    serial_range_query_counts[i] = num_in_range;
    serial_range_query_sums[i] = sum_in_range;
  }
  end_time = get_usecs();
  serial_time = end_time - start_time;
  printf("\t did %lu range queries serially in %lu\n", NUM_QUERIES, serial_time);

  // correctness check of serial 
  for (size_t i = 0; i < correct_range_query_sums.size(); i++) {
    if (correct_range_query_sums[i] != serial_range_query_sums[i]) {
      printf("wrong serial range query sum, expected %lu, got %lu\n", correct_range_query_sums[i], serial_range_query_sums[i]);
      wrong = true;
    }
    if (correct_range_query_counts[i] != serial_range_query_counts[i]) {
      printf("wrong serial range query count, expected %lu, got %lu\n", correct_range_query_counts[i], serial_range_query_counts[i]);
      wrong = true;
    }
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  } 
  */

  start_time = get_usecs();
  // concurrent btree range sums
  cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
    start = std::get<0>(checker_sorted[range_query_start_idxs[i]]);
    end = std::get<0>(checker_sorted[std::min(range_query_start_idxs[i] + range_query_lengths[i], checker_sorted.size() - 1)]);

    uint64_t num_in_range = 0;
    uint64_t sum_in_range = 0;
    concurrent_map.map_range(start, end, [&num_in_range, &sum_in_range]([[maybe_unused]] auto key, auto val) {
              num_in_range += 1;
              sum_in_range += val;
            });
    concurrent_range_query_counts[i] = num_in_range;
    concurrent_range_query_sums[i] = sum_in_range;
  }
  end_time = get_usecs();
  concurrent_time = end_time - start_time;
  printf("\t did %lu range queries concurrently in %lu\n", NUM_QUERIES, concurrent_time);

#if CORRECTNESS
  // correctness check of concurrent 
  wrong = false;
  for (size_t i = 0; i < correct_range_query_sums.size(); i++) {
    if (correct_range_query_sums[i] != concurrent_range_query_sums[i]) {
      printf("wrong concurrent range query sum, expected %lu, got %lu\n", correct_range_query_sums[i], concurrent_range_query_sums[i]);
      wrong = true;
    }
    if (correct_range_query_counts[i] != concurrent_range_query_counts[i]) {
      printf("wrong concurrent range query count, expected %lu, got %lu\n", correct_range_query_counts[i], concurrent_range_query_counts[i]);
      wrong = true;
    }
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  } 
#endif
  // return {true, 0, 0, 0, 0};
  // */

  // *** testing map_range_length ***

  std::vector<T> serial_range_query_length_sums(NUM_QUERIES);
  std::vector<uint64_t> serial_range_query_length_counts(NUM_QUERIES);

  std::vector<T> concurrent_range_query_length_sums(NUM_QUERIES);
  std::vector<uint64_t> concurrent_range_query_length_counts(NUM_QUERIES);

  /*
  start_time = get_usecs();
  // serial btree range sums
  for (uint32_t i = 0; i < NUM_QUERIES; i++) {
    start = std::get<0>(checker_sorted[range_query_start_idxs[i]]);

    uint64_t num_in_range = 0;
    uint64_t sum_in_range = 0;
    serial_map.map_range_length(start, correct_range_query_counts[i], [&num_in_range, &sum_in_range]([[maybe_unused]] auto key, auto val) {
              num_in_range += 1;
              sum_in_range += val;
            });
    serial_range_query_length_counts[i] = num_in_range;
    serial_range_query_length_sums[i] = sum_in_range;
  }
  end_time = get_usecs();
  serial_time = end_time - start_time;
  printf("\t did %lu range queries by length serially in %lu\n", NUM_QUERIES, serial_time);

  // correctness check of serial 
  wrong = false;
  for (size_t i = 0; i < correct_range_query_sums.size(); i++) {
    if (correct_range_query_sums[i] != serial_range_query_length_sums[i]) {
      printf("wrong serial range query length sum, expected %lu, got %lu\n", correct_range_query_sums[i], serial_range_query_length_sums[i]);
      wrong = true;
    }
    if (correct_range_query_counts[i] != serial_range_query_length_counts[i]) {
      printf("wrong serial range query length count, expected %lu, got %lu\n", correct_range_query_counts[i], serial_range_query_length_counts[i]);
      wrong = true;
    }
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  } 
  */

  start_time = get_usecs();
  // concurrent btree range sums
  cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
    start = std::get<0>(checker_sorted[range_query_start_idxs[i]]);
    // end = checker_sorted[std::min(range_query_start_idxs[i] + range_query_lengths[i], checker_sorted.size() - 1)];

    uint64_t num_in_range = 0;
    uint64_t sum_in_range = 0;

#if CORRECTNESS
    concurrent_map.map_range_length(start, correct_range_query_counts[i], [&num_in_range, &sum_in_range]([[maybe_unused]] auto key, auto val) {
              num_in_range += 1;
              sum_in_range += val;
            });
#else
    concurrent_map.map_range_length(start, concurrent_range_query_counts[i], [&num_in_range, &sum_in_range]([[maybe_unused]] auto key, auto val) {
              num_in_range += 1;
              sum_in_range += val;
            });
#endif
    concurrent_range_query_length_counts[i] = num_in_range;
    concurrent_range_query_length_sums[i] = sum_in_range;
  }
  end_time = get_usecs();
  concurrent_time = end_time - start_time;
  printf("\t did %lu range queries by length concurrently in %lu\n", NUM_QUERIES, concurrent_time);

#if CORRECTNESS
  // correctness check of concurrent 
  wrong = false;
  for (size_t i = 0; i < correct_range_query_sums.size(); i++) {
    if (correct_range_query_sums[i] != concurrent_range_query_length_sums[i]) {
      printf("wrong concurrent range query length sum, expected %lu, got %lu\n", correct_range_query_sums[i], concurrent_range_query_length_sums[i]);
      wrong = true;
    }
    if (correct_range_query_counts[i] != concurrent_range_query_length_counts[i]) {
      printf("wrong concurrent range query length count, expected %lu, got %lu\n", correct_range_query_counts[i], concurrent_range_query_length_counts[i]);
      wrong = true;
    }
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  } 
#endif
  return {true, 0, 0, 0, 0};
}

template <class T, uint32_t internal_bytes>
bool
test_sequential_inserts(uint64_t max_size, int trials) {
  printf("*** testing sequential inserts ***\n");
  uint64_t start_time, end_time;
  std::vector<uint64_t> insert_times;

  for (int cur_trial = 0; cur_trial <= trials; cur_trial++) {
    printf("trials %d, inserts %lu\n", cur_trial, max_size);
    tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
                  std::allocator<T>, true> concurrent_map;

    uint32_t num_threads = 48;
    uint32_t inserts_per_thread = max_size / num_threads;
    // TIME INSERTS

    start_time = get_usecs();
		parallel_for(0, num_threads, [&](const uint32_t &thread_id) {
      uint32_t num_to_insert = inserts_per_thread;
      if (thread_id == num_threads - 1) {
        num_to_insert += (max_size % num_threads);
      }
      for(uint64_t i = 0; i < num_to_insert; i++) {
        uint64_t val_to_insert = thread_id * 1000000000 + i + 1;
        concurrent_map.insert({val_to_insert, val_to_insert});
      }
    });
    end_time = get_usecs();
    if (cur_trial > 0) {insert_times.push_back(end_time - start_time);}
    printf("\tDone inserting %lu elts in %lu\n",max_size, end_time - start_time);
    
   auto concurrent_sum = concurrent_map.psum_with_subtract();
   printf("concurrent sum = %lu\n", concurrent_sum);
  }
  std::sort(insert_times.begin(), insert_times.end());
  printf("median insert time = %lu\n", insert_times[trials / 2]);
  return true;
}

template <class T, uint32_t internal_bytes, uint32_t header_size, uint32_t block_size>
bool
test_map_cache_misses(uint64_t max_size, uint64_t NUM_QUERIES, uint32_t MAX_QUERY_SIZE, std::seed_seq &seed, int trials) {
  printf("num trials = %d\n", trials);
  for (int cur_trial = 0; cur_trial <= trials; cur_trial++) {
    printf("\nRunning leafds btree with internal bytes = %u with leafds slots %d, trial = %d, total trials %d\n",internal_bytes, SLOTS, cur_trial, trials);
    printf("inserting %lu elts, num queries = %lu, max query size = %u\n", max_size * 2, NUM_QUERIES, MAX_QUERY_SIZE);

    std::vector<T> data =
        create_random_data<T>(2*max_size, std::numeric_limits<T>::max(), seed);
    
    std::set<T> checker_set;
    bool wrong;

    std::seed_seq query_seed{1};
    std::seed_seq query_seed_2{2};

  #if CORRECTNESS
    std::vector<uint64_t> range_query_start_idxs =
        create_random_data<uint64_t>(NUM_QUERIES, checker_sorted.size() - 1, query_seed);
  #else
    std::vector<uint64_t> range_query_start_idxs =
        create_random_data<uint64_t>(NUM_QUERIES, data.size() - 1, query_seed);
  #endif

    tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes, header_size, block_size>,
                  std::allocator<T>, true> concurrent_map;

    // TIME INSERTS
    parallel_for(0, 2*max_size, [&](const uint32_t &i) {
      concurrent_map.insert({data[i], 2*data[i]});
    });

    printf("trial %d: inserted %lu\n", cur_trial, 2*max_size);
    // do range queries for given length 
    std::vector<uint64_t> range_query_lengths =
      create_random_data<uint64_t>(NUM_QUERIES, MAX_QUERY_SIZE - 1, query_seed_2);

    // sorted query first to get end key
    std::vector<T> concurrent_range_query_length_maxs(NUM_QUERIES);
    std::vector<T> concurrent_range_query_length_val_sums(NUM_QUERIES);

    printf("*** MAP BY LENGTH ***\n");
    parallel_for(0, NUM_QUERIES, [&](const uint32_t &i) {
      T start = data[range_query_start_idxs[i]];
      T sum_val_range = 0;
      T max_key_in_range = 0;
      concurrent_map.map_range_length(start, range_query_lengths[i], [&max_key_in_range, &sum_val_range]([[maybe_unused]] auto key, auto val) {
                sum_val_range += val;
                if (key > max_key_in_range) {
                  max_key_in_range = key;
                }
              });
      concurrent_range_query_length_val_sums[i] = sum_val_range;
      concurrent_range_query_length_maxs[i] = max_key_in_range;
    });
    T sum_all_vals = 0;
    for (auto e: concurrent_range_query_length_val_sums) {
      sum_all_vals += e;
    }

    printf("\t\t\tsum vals %lu\n", sum_all_vals);

/*
    printf("*** MAP UNSORTED ***\n");
    std::vector<T> concurrent_range_query_key_sums(NUM_QUERIES);
    std::vector<T> concurrent_range_query_val_sums(NUM_QUERIES);

    parallel_for(0, NUM_QUERIES, [&](const uint32_t &i) {
      T start, end;
#if CORRECTNESS
      start = checker_sorted[range_query_start_idxs[i]];
#else 
      start = data[range_query_start_idxs[i]];
#endif
      end = concurrent_range_query_length_maxs[i];
      if (range_query_lengths[i] != 0) {
        end++;
      }

      T sum_key_range = 0;
      T sum_val_range = 0;

      concurrent_map.map_range(start, end, [&sum_key_range, &sum_val_range]([[maybe_unused]] auto key, auto val) {
                sum_key_range += key;
                sum_val_range += val;
              });
      concurrent_range_query_key_sums[i] = sum_key_range;
      concurrent_range_query_val_sums[i] = sum_val_range;
    });
    T sum_all_keys = 0;
    for (auto e: concurrent_range_query_key_sums) {
      sum_all_keys += e;
    }
    sum_all_vals = 0;
    for (auto e: concurrent_range_query_val_sums) {
      sum_all_vals += e;
    }
    if (sum_all_keys * 2 != sum_all_vals) {
      printf("\t\t\t wrong, sum keys * 2 not equal to sum vals\n");
      // return false;
    }
    printf("\t\t\tmap range sum keys %lu sum vals %lu\n", sum_all_keys, sum_all_vals);
    printf("END: trial %d, total trials %d\n", cur_trial, trials);

*/
  }
  return true;
}

template <class T, uint32_t internal_bytes>
bool
test_serial_microbenchmarks_map(uint64_t max_size, uint64_t NUM_QUERIES, std::seed_seq &seed, bool write_csv, int trials) {
  // std::vector<uint32_t> num_query_sizes{1000};
  std::vector<uint32_t> num_query_sizes{100, 1000, 10000, 100000};

  uint64_t start_time, end_time;
  std::vector<uint64_t> insert_times;
  std::vector<uint64_t> find_times;
  std::vector<uint64_t> sorted_range_query_times_by_size;
  std::vector<uint64_t> unsorted_range_query_times_by_size;

  for (int cur_trial = 0; cur_trial < trials; cur_trial++) {
    printf("\nRunning leafds btree with internal bytes = %u with leafds slots %d, trial = %d\n",internal_bytes, SLOTS, cur_trial);

    std::vector<T> data =
        create_random_data<T>(2*max_size, std::numeric_limits<T>::max(), seed);
    
    std::set<T> checker_set;
    bool wrong;

    std::seed_seq query_seed{1};
    std::seed_seq query_seed_2{2};

    std::vector<uint64_t> range_query_start_idxs =
        create_random_data<uint64_t>(NUM_QUERIES, data.size() - 1, query_seed);
    tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
                  std::allocator<T>, false> serial_map;

    // TIME INSERTS
    for(uint32_t i = 0; i < max_size; i++) {
      serial_map.insert({data[i], 2*data[i]});
    }
    printf("trial %d: inserted %lu\n", cur_trial, max_size);

    start_time = get_usecs();
    for(uint32_t i = max_size; i < 2*max_size; i++) {
      serial_map.insert({data[i], 2*data[i]});
    }
    end_time = get_usecs();
    insert_times.push_back(end_time - start_time);
    printf("\tDone inserting %lu elts in %lu\n",max_size, end_time - start_time);
    
    // TIME POINT QUERIES
    std::vector<bool> found_count(NUM_QUERIES);
    start_time = get_usecs();
    for(uint32_t i = 0; i < NUM_QUERIES; i++) {
      found_count[i] = serial_map.exists(data[range_query_start_idxs[i]]);
    }
    end_time = get_usecs();
    find_times.push_back(end_time - start_time);
    int count_found = 0;
    for (auto e : found_count) {
      count_found += e ? 1 : 0;
    }
    printf("\tDone finding %lu elts in %lu, count = %d \n",NUM_QUERIES, end_time - start_time, count_found);

    // TIME RANGE QUERIES FOR ALL LENGTHS
    for (size_t query_size_i = 0; query_size_i < num_query_sizes.size(); query_size_i++) {
      uint64_t MAX_QUERY_SIZE = num_query_sizes[query_size_i];
      std::vector<uint64_t> range_query_lengths =
        create_random_data<uint64_t>(NUM_QUERIES, MAX_QUERY_SIZE - 1, query_seed_2);

      // sorted query first to get end key
      std::vector<T> concurrent_range_query_length_maxs(NUM_QUERIES);
      std::vector<uint64_t> concurrent_range_query_length_counts(NUM_QUERIES);
      std::vector<T> concurrent_range_query_length_key_sums(NUM_QUERIES);
      std::vector<T> concurrent_range_query_length_val_sums(NUM_QUERIES);

      for (uint32_t i = 0; i < NUM_QUERIES; i++) {
        T start;
        start = data[range_query_start_idxs[i]];

        uint64_t num_in_range = 0;
        T max_key_in_range = start;

        serial_map.map_range_length(start, range_query_lengths[i], [&max_key_in_range]([[maybe_unused]] auto key, auto val) {
                  if (key > max_key_in_range) {
                    max_key_in_range = key;
                  }
                });
        serial_map.map_range_length(start, range_query_lengths[i], [&num_in_range]([[maybe_unused]] auto key, auto val) {
                  num_in_range++;
                });
        concurrent_range_query_length_counts[i] = num_in_range;
        concurrent_range_query_length_maxs[i] = max_key_in_range;
      }

      start_time = get_usecs();
      for (uint32_t i = 0; i < NUM_QUERIES; i++) {
        T start;
        start = data[range_query_start_idxs[i]];
        T sum_key_range = 0;
        T sum_val_range = 0;
        serial_map.map_range_length(start, range_query_lengths[i], [&sum_key_range, &sum_val_range]([[maybe_unused]] auto key, auto val) {
                  sum_key_range += key;
                  sum_val_range += val;
                });
        concurrent_range_query_length_key_sums[i] = sum_key_range;
        concurrent_range_query_length_val_sums[i] = sum_val_range;
      }
      end_time = get_usecs();
      sorted_range_query_times_by_size.push_back(end_time - start_time);
      printf("\t\t did sorted range queries with max len %lu serially in %lu\n", MAX_QUERY_SIZE, end_time - start_time);
      T sum_all_keys = 0;
      for (auto e: concurrent_range_query_length_key_sums) {
        sum_all_keys += e;
      }
      T sum_all_vals = 0;
      for (auto e: concurrent_range_query_length_val_sums) {
        sum_all_vals += e;
      }
      if (sum_all_keys * 2 != sum_all_vals) {
        printf("\t\t\t wrong, sum keys * 2 %lu not equal to sum vals %lu\n", sum_all_keys, sum_all_vals);
        // return false;
      }
      printf("\t\t\t sum keys %lu sum vals %lu\n", sum_all_keys, sum_all_vals);
      std::vector<T> concurrent_range_query_key_sums(NUM_QUERIES);
      std::vector<T> concurrent_range_query_val_sums(NUM_QUERIES);

      start_time = get_usecs();
      for (uint32_t i = 0; i < NUM_QUERIES; i++) {
        T start, end;

        start = data[range_query_start_idxs[i]];
        end = concurrent_range_query_length_maxs[i];
        if (range_query_lengths[i] != 0) {
          end++;
        }

        T sum_key_range = 0;
        T sum_val_range = 0;
        serial_map.map_range(start, end, [&sum_key_range, &sum_val_range]([[maybe_unused]] auto key, auto val) {
                  sum_key_range += key;
                  sum_val_range += val;
                });
        concurrent_range_query_key_sums[i] = sum_key_range;
        concurrent_range_query_val_sums[i] = sum_val_range;
      }
      end_time = get_usecs();
      unsorted_range_query_times_by_size.push_back(end_time - start_time);
      printf("\t\t did unsorted range queries with max len %lu serially in %lu\n", MAX_QUERY_SIZE, end_time - start_time);
      sum_all_keys = 0;
      for (auto e: concurrent_range_query_key_sums) {
        sum_all_keys += e;
      }
      sum_all_vals = 0;
      for (auto e: concurrent_range_query_val_sums) {
        sum_all_vals += e;
      }
      if (sum_all_keys * 2 != sum_all_vals) {
        printf("\t\t\t wrong, sum keys * 2 not equal to sum vals\n");
        // return false;
      }
    }
  }

  // DONE WITH TESTS FOR THIS SIZE
  if (!write_csv) { 
    printf("correct\n");
  } else {
    // printf("tree_type, internal bytes, leaf bytes, num_inserted, insert_time, num_finds, find_time, \n");
    // tree_type, internal bytes, leaf bytes, num_inserted,num_range_queries, max_query_size,  unsorted_query_time, sorted_query_time
    std::ofstream outfile;
    uint64_t insert_time_med, find_time_med;
    std::sort(insert_times.begin(), insert_times.end());
    std::sort(find_times.begin(), find_times.end());
    insert_time_med = insert_times[trials/2];
    find_time_med = find_times[trials/2];

    outfile.open("insert_finds.csv", std::ios_base::app); 
    outfile << "leafds, " << internal_bytes << ", " << SLOTS << ", " <<  max_size << ", " <<  insert_time_med << ", " <<  NUM_QUERIES << ", " <<  find_time_med << ", \n";
    outfile.close();
    outfile.open("range_queries.csv", std::ios_base::app); 
    printf("num sorted trials %lu, num unsorted trials %lu\n", sorted_range_query_times_by_size.size(), unsorted_range_query_times_by_size.size());
    for (size_t i = 0; i < num_query_sizes.size(); i++) {
      std::vector<uint64_t> curr_unsorted_query_times;
      std::vector<uint64_t> curr_sorted_query_times;
      for (int t = 0; t < trials; t++) {
        curr_unsorted_query_times.push_back(unsorted_range_query_times_by_size[t*num_query_sizes.size() + i]);
        curr_sorted_query_times.push_back(sorted_range_query_times_by_size[t*num_query_sizes.size() + i]);
      }
      std::sort(curr_unsorted_query_times.begin(), curr_unsorted_query_times.end());
      std::sort(curr_sorted_query_times.begin(), curr_sorted_query_times.end());
      outfile << "leafds, " << internal_bytes << ", " << SLOTS << ", " <<  max_size << ", " <<  NUM_QUERIES << ", " <<  num_query_sizes[i] << ", " <<  curr_unsorted_query_times[trials/2] << ", " <<  curr_sorted_query_times[trials/2] << ", \n";
    }
    outfile.close();
  }
  return true;
}




 
template <class T, uint32_t internal_bytes, uint32_t header_size, uint32_t block_size>
bool
test_concurrent_microbenchmarks_map(uint64_t max_size, uint64_t NUM_QUERIES, std::seed_seq &seed, bool write_csv, int trials) {
  // std::vector<uint32_t> num_query_sizes{1000};
  std::vector<uint32_t> num_query_sizes{100, 1000, 10000, 100000};

  uint64_t start_time, end_time;
  std::vector<uint64_t> insert_times;
  std::vector<uint64_t> find_times;
  std::vector<uint64_t> sorted_range_query_times_by_size;
  std::vector<uint64_t> unsorted_range_query_times_by_size;

  for (int cur_trial = 0; cur_trial <= trials; cur_trial++) {
    printf("\nRunning leafds btree with internal bytes = %u with header size %u, block size %u, trial = %d\n",internal_bytes, header_size, block_size, cur_trial);

    // std::vector<uint32_t> num_query_sizes{100, 1000, 10000, 100000};

    // std::vector<T> data =
    //     create_random_data<T>(max_size, max_size * 2, seed);
    std::vector<T> data =
        create_random_data<T>(2*max_size, std::numeric_limits<T>::max(), seed);
    
    std::set<T> checker_set;
    bool wrong;

  #if CORRECTNESS
    for (uint32_t i = 0; i < 2*max_size; i++) {
      checker_set.insert(data[i]);
    }

    std::vector<T> checker_sorted;
    for (auto key: checker_set) {
      checker_sorted.push_back(key);
    }
    std::sort(checker_sorted.begin(), checker_sorted.end());
  #endif

    std::seed_seq query_seed{1};
    std::seed_seq query_seed_2{2};

  #if CORRECTNESS
    std::vector<uint64_t> range_query_start_idxs =
        create_random_data<uint64_t>(NUM_QUERIES, checker_sorted.size() - 1, query_seed);
  #else
    std::vector<uint64_t> range_query_start_idxs =
        create_random_data<uint64_t>(NUM_QUERIES, data.size() - 1, query_seed);
  #endif
    tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes, header_size, block_size>,
                  std::allocator<T>, true> concurrent_map;

    // TIME INSERTS
    parallel_for(0, max_size, [&](const uint32_t &i) {
      concurrent_map.insert({data[i], 2*data[i]});
    });
    printf("trial %d: inserted %lu\n", cur_trial, max_size);

    start_time = get_usecs();
    parallel_for(max_size, 2*max_size, [&](const uint32_t &i) {
      concurrent_map.insert({data[i], 2*data[i]});
    });
    end_time = get_usecs();
    if (cur_trial > 0) {insert_times.push_back(end_time - start_time);}
    printf("\tDone inserting %lu elts in %lu\n",max_size, end_time - start_time);
    
    // TIME POINT QUERIES
    std::vector<bool> found_count(NUM_QUERIES);
    start_time = get_usecs();
    parallel_for(0, NUM_QUERIES, [&](const uint32_t &i) {
      found_count[i] = concurrent_map.exists(data[range_query_start_idxs[i]]);
    });
    end_time = get_usecs();
    if (cur_trial > 0) {find_times.push_back(end_time - start_time);}
    int count_found = 0;
    for (auto e : found_count) {
      count_found += e ? 1 : 0;
    }
    printf("\tDone finding %lu elts in %lu, count = %d \n",NUM_QUERIES, end_time - start_time, count_found);

    // TIME RANGE QUERIES FOR ALL LENGTHS
    for (size_t query_size_i = 0; query_size_i < num_query_sizes.size(); query_size_i++) {
      uint64_t MAX_QUERY_SIZE = num_query_sizes[query_size_i];
      std::vector<uint64_t> range_query_lengths =
        create_random_data<uint64_t>(NUM_QUERIES, MAX_QUERY_SIZE - 1, query_seed_2);

      // T start, end;
  #if CORRECTNESS
      std::vector<T> correct_range_query_maxs(NUM_QUERIES);
      std::vector<uint64_t> correct_range_query_counts(NUM_QUERIES);

      // get correct range sums
      parallel_for(0, NUM_QUERIES, [&](const uint32_t &i) {
        T start;
        start = checker_sorted[range_query_start_idxs[i]];

        correct_range_query_counts[i] = 0;
        correct_range_query_maxs[i] = start;

        size_t curr_index = range_query_start_idxs[i];

        for (uint64_t count = 0; count < range_query_lengths[i]; count++) {
          if (curr_index >= checker_sorted.size()) {
            break;
          }
          correct_range_query_maxs[i] = checker_sorted[curr_index];
          correct_range_query_counts[i] += 1;
          curr_index++;
        }
      });
      printf("\t did %lu correct range queries with max query size %lu on checker set\n", 
            NUM_QUERIES,
            MAX_QUERY_SIZE);
  #endif

      // sorted query first to get end key
      std::vector<T> concurrent_range_query_length_maxs(NUM_QUERIES);
      std::vector<uint64_t> concurrent_range_query_length_counts(NUM_QUERIES);
      std::vector<T> concurrent_range_query_length_key_sums(NUM_QUERIES);
      std::vector<T> concurrent_range_query_length_val_sums(NUM_QUERIES);

      parallel_for(0, NUM_QUERIES, [&](const uint32_t &i) {
        T start;
  #if CORRECTNESS
        start = checker_sorted[range_query_start_idxs[i]];
  #else 
        start = data[range_query_start_idxs[i]];
  #endif

        uint64_t num_in_range = 0;
        T max_key_in_range = start;

        concurrent_map.map_range_length(start, range_query_lengths[i], [&max_key_in_range]([[maybe_unused]] auto key, auto val) {
                  if (key > max_key_in_range) {
                    max_key_in_range = key;
                  }
                });
        concurrent_map.map_range_length(start, range_query_lengths[i], [&num_in_range]([[maybe_unused]] auto key, auto val) {
                  num_in_range++;
                });
        concurrent_range_query_length_counts[i] = num_in_range;
        concurrent_range_query_length_maxs[i] = max_key_in_range;
      });
      start_time = get_usecs();
      parallel_for(0, NUM_QUERIES, [&](const uint32_t &i) {
        T start;
  #if CORRECTNESS
        start = checker_sorted[range_query_start_idxs[i]];
  #else 
        start = data[range_query_start_idxs[i]];
  #endif
        T sum_key_range = 0;
        T sum_val_range = 0;
        concurrent_map.map_range_length(start, range_query_lengths[i], [&sum_key_range, &sum_val_range]([[maybe_unused]] auto key, auto val) {
                  sum_key_range += key;
                  sum_val_range += val;
                });
        concurrent_range_query_length_key_sums[i] = sum_key_range;
        concurrent_range_query_length_val_sums[i] = sum_val_range;
      });
      end_time = get_usecs();
      if (cur_trial > 0) {sorted_range_query_times_by_size.push_back(end_time - start_time);}
      printf("\t\t did sorted range queries with max len %lu concurrently in %lu\n", MAX_QUERY_SIZE, end_time - start_time);
      T sum_all_keys = 0;
      for (auto e: concurrent_range_query_length_key_sums) {
        sum_all_keys += e;
      }
      T sum_all_vals = 0;
      for (auto e: concurrent_range_query_length_val_sums) {
        sum_all_vals += e;
      }
      if (sum_all_keys * 2 != sum_all_vals) {
        printf("\t\t\t wrong, sum keys * 2 %lu not equal to sum vals %lu\n", sum_all_keys, sum_all_vals);
        // return false;
      }
      printf("\t\t\t sum keys %lu sum vals %lu\n", sum_all_keys, sum_all_vals);

  #if CORRECTNESS
      // correctness check of concurrent sorted query
      printf("*** CHECKING CONCURRENT SORTED QUERY ***\n");
      wrong = false;
      for (size_t i = 0; i < NUM_QUERIES; i++) {
        if (correct_range_query_counts[i] != concurrent_range_query_length_counts[i]) {
          printf("wrong concurrent range query length count at %lu with len %lu, expected %lu, got %lu\n", i, range_query_lengths[i], correct_range_query_counts[i], concurrent_range_query_length_counts[i]);
          T num_in_range = 0;
          printf("got: \n");

          T start = checker_sorted[range_query_start_idxs[i]];
          concurrent_map.map_range_length(start, range_query_lengths[i], [&num_in_range]([[maybe_unused]] auto key, auto val) {
                    printf("\telt %lu: %lu\n", num_in_range, key);
                    num_in_range++;
                  });
            wrong = true;
          printf("\n\nshould be:\n");
          for (uint64_t count = 0; count < range_query_lengths[i]; count++) {
            T curr_index = range_query_start_idxs[i] + count;
            if (curr_index >= checker_sorted.size()) {
              break;
            }
            printf("\telt %lu: %lu\n", count, checker_sorted[curr_index]);
          }
        }
        if (correct_range_query_maxs[i] != concurrent_range_query_length_maxs[i]) {
          printf("wrong concurrent range query length max, expected %lu, got %lu\n", correct_range_query_maxs[i], concurrent_range_query_length_maxs[i]);
          wrong = true;
        }
      }
      if (wrong) {
        return false;
      }
      printf("*** PASSED CONCURRENT SORTED QUERY ***\n");
  #endif

#if CORRECTNESS
      std::vector<T> concurrent_range_query_maxs(NUM_QUERIES);
      std::vector<uint64_t> concurrent_range_query_counts(NUM_QUERIES);
#endif
      std::vector<T> concurrent_range_query_key_sums(NUM_QUERIES);
      std::vector<T> concurrent_range_query_val_sums(NUM_QUERIES);

#if CORRECTNESS
      parallel_for(0, NUM_QUERIES, [&](const uint32_t &i) {
        T start, end;
        start = checker_sorted[range_query_start_idxs[i]];
        end = concurrent_range_query_length_maxs[i];
        if (range_query_lengths[i] != 0) {
          end++;
        }
        uint64_t num_in_range = 0;
        T max_key_in_range = start;
        concurrent_map.map_range(start, end, [&num_in_range]([[maybe_unused]] auto key, auto val) {
                  num_in_range += 1;
                });

        concurrent_map.map_range(start, end, [&max_key_in_range]([[maybe_unused]] auto key, auto val) {
                  if (key > max_key_in_range) {
                    max_key_in_range = key;
                  }
                });
        concurrent_range_query_counts[i] = num_in_range;
        concurrent_range_query_maxs[i] = max_key_in_range;
      });
      // printf("\n");
#endif

      start_time = get_usecs();
      parallel_for(0, NUM_QUERIES, [&](const uint32_t &i) {
        T start, end;
  #if CORRECTNESS
        start = checker_sorted[range_query_start_idxs[i]];
  #else 
        start = data[range_query_start_idxs[i]];
  #endif
        end = concurrent_range_query_length_maxs[i];
        if (range_query_lengths[i] != 0) {
          end++;
        }

        T sum_key_range = 0;
        T sum_val_range = 0;
        concurrent_map.map_range(start, end, [&sum_key_range, &sum_val_range]([[maybe_unused]] auto key, auto val) {
                  sum_key_range += key;
                  sum_val_range += val;
                });
        concurrent_range_query_key_sums[i] = sum_key_range;
        concurrent_range_query_val_sums[i] = sum_val_range;
      });
      end_time = get_usecs();
      if (cur_trial > 0) {unsorted_range_query_times_by_size.push_back(end_time - start_time);}
      printf("\t\t did unsorted range queries with max len %lu concurrently in %lu\n", MAX_QUERY_SIZE, end_time - start_time);
      sum_all_keys = 0;
      for (auto e: concurrent_range_query_key_sums) {
        sum_all_keys += e;
      }
      sum_all_vals = 0;
      for (auto e: concurrent_range_query_val_sums) {
        sum_all_vals += e;
      }
      if (sum_all_keys * 2 != sum_all_vals) {
        printf("\t\t\t wrong, sum keys * 2 not equal to sum vals\n");
        // return false;
      }
      printf("\t\t\t sum keys %lu sum vals %lu\n", sum_all_keys, sum_all_vals);
  #if CORRECTNESS
      // correctness check of concurrent 
      wrong = false;

      printf("*** CHECKING UNSORTED CONCURRENT QUERY ***\n");
      for (size_t i = 0; i < NUM_QUERIES; i++) {
        if (correct_range_query_counts[i] != concurrent_range_query_counts[i]) {
          printf("wrong concurrent range query count, expected %lu, got %lu\n", correct_range_query_counts[i], concurrent_range_query_counts[i]);
          wrong = true;
        }
        if (correct_range_query_maxs[i] != concurrent_range_query_maxs[i]) {
          printf("wrong concurrent range query max, expected %lu, got %lu\n", correct_range_query_maxs[i], concurrent_range_query_maxs[i]);
          wrong = true;
        }
      }
      if (wrong) {
        return false;
      }
      printf("*** PASSED UNSORTED CONCURRENT QUERY ***\n");
  #endif
 
      /* 
#if CORRECTNESS
      std::vector<T> concurrent_range_query_sorted_end_maxs(NUM_QUERIES);
      std::vector<uint64_t> concurrent_range_query_sorted_end_counts(NUM_QUERIES);
#endif
      std::vector<T> concurrent_range_query_sorted_end_key_sums(NUM_QUERIES);
      std::vector<T> concurrent_range_query_sorted_end_val_sums(NUM_QUERIES);

#if CORRECTNESS
      parallel_for(0, NUM_QUERIES, [&](const uint32_t &i) {
        T start, end;
        start = checker_sorted[range_query_start_idxs[i]];
        end = concurrent_range_query_length_maxs[i];
        if (range_query_lengths[i] != 0) {
          end++;
        }
        uint64_t num_in_range = 0;
        T max_key_in_range = start;
        concurrent_map.map_range_sorted_end(start, end, [&num_in_range, &max_key_in_range]([[maybe_unused]] auto key, auto val) {
                  num_in_range += 1;
                  if (key > max_key_in_range) {
                    max_key_in_range = key;
                  }
                });
        concurrent_range_query_sorted_end_counts[i] = num_in_range;
        concurrent_range_query_sorted_end_maxs[i] = max_key_in_range;
      });
#endif

      start_time = get_usecs();
      parallel_for(0, NUM_QUERIES, [&](const uint32_t &i) {
        T start, end;
  #if CORRECTNESS
        start = checker_sorted[range_query_start_idxs[i]];
  #else 
        start = data[range_query_start_idxs[i]];
  #endif
        end = concurrent_range_query_length_maxs[i];
        if (range_query_lengths[i] != 0) {
          end++;
        }

        T sum_key_range = 0;
        T sum_val_range = 0;
        concurrent_map.map_range_sorted_end(start, end, [&sum_key_range, &sum_val_range]([[maybe_unused]] auto key, auto val) {
                  sum_key_range += key;
                  sum_val_range += val;
                });
        concurrent_range_query_sorted_end_key_sums[i] = sum_key_range;
        concurrent_range_query_sorted_end_val_sums[i] = sum_val_range;
      });
      end_time = get_usecs();
      if (cur_trial > 0) {unsorted_range_query_times_by_size.push_back(end_time - start_time);}
      printf("\t\t did sorted with end range queries with max len %lu concurrently in %lu\n", MAX_QUERY_SIZE, end_time - start_time);
      sum_all_keys = 0;
      for (auto e: concurrent_range_query_sorted_end_key_sums) {
        sum_all_keys += e;
      }
      sum_all_vals = 0;
      for (auto e: concurrent_range_query_sorted_end_val_sums) {
        sum_all_vals += e;
      }
      if (sum_all_keys * 2 != sum_all_vals) {
        printf("\t\t\t wrong, sum keys * 2 not equal to sum vals\n");
        // return false;
      }
      printf("\t\t\t sum keys %lu sum vals %lu\n", sum_all_keys, sum_all_vals);
  #if CORRECTNESS
      // correctness check of concurrent 
      wrong = false;
      for (size_t i = 0; i < NUM_QUERIES; i++) {
        if (correct_range_query_counts[i] != concurrent_range_query_sorted_end_counts[i]) {
          printf("wrong concurrent sorted end range query count, expected %lu, got %lu\n", correct_range_query_counts[i], concurrent_range_query_sorted_end_counts[i]);
          wrong = true;
        }
        if (correct_range_query_maxs[i] != concurrent_range_query_sorted_end_maxs[i]) {
          printf("wrong concurrent sorted end range query max, expected %lu, got %lu\n", correct_range_query_maxs[i], concurrent_range_query_sorted_end_maxs[i]);
          wrong = true;
        }
      }
      if (wrong) {
        return false;
      } 
  #endif

  #if CORRECTNESS 
      // non checker_sorted correctness check
      wrong = false;
      for (size_t i = 0; i < NUM_QUERIES; i++) {
        if (concurrent_range_query_length_counts[i] != concurrent_range_query_counts[i]) {
          printf("wrong concurrent range query count, expected %lu, got %lu\n", concurrent_range_query_length_counts[i], concurrent_range_query_counts[i]);
          wrong = true;
        }
        if (concurrent_range_query_length_maxs[i] != concurrent_range_query_maxs[i]) {
          printf("wrong concurrent range query max, expected %lu, got %lu\n", concurrent_range_query_length_maxs[i], concurrent_range_query_maxs[i]);
          wrong = true;
        }
      }
      if (wrong) { return false;}
  #endif

    */
    }
  }

  // DONE WITH TESTS FOR THIS SIZE
  if (!write_csv) { 
    printf("correct\n");
  } else {
    // printf("tree_type, internal bytes, leaf bytes, num_inserted, insert_time, num_finds, find_time, \n");
    // tree_type, internal bytes, leaf bytes, num_inserted,num_range_queries, max_query_size,  unsorted_query_time, sorted_query_time
    std::ofstream outfile;
    uint64_t insert_time_med, find_time_med;
    std::sort(insert_times.begin(), insert_times.end());
    std::sort(find_times.begin(), find_times.end());
    insert_time_med = insert_times[trials/2];
    find_time_med = find_times[trials/2];

    outfile.open("insert_finds.csv", std::ios_base::app); 
    outfile << "leafds, " << internal_bytes << ", " << header_size << ", " << block_size << ", " <<  max_size << ", " <<  insert_time_med << ", " <<  NUM_QUERIES << ", " <<  find_time_med << ", \n";
    outfile.close();
    outfile.open("range_queries.csv", std::ios_base::app); 
    printf("num sorted trials %lu, num unsorted trials %lu\n", sorted_range_query_times_by_size.size(), unsorted_range_query_times_by_size.size());
    for (size_t i = 0; i < num_query_sizes.size(); i++) {
      std::vector<uint64_t> curr_unsorted_query_times;
      std::vector<uint64_t> curr_sorted_query_times;
      for (int t = 0; t < trials; t++) {
        curr_unsorted_query_times.push_back(unsorted_range_query_times_by_size[t*num_query_sizes.size() + i]);
        curr_sorted_query_times.push_back(sorted_range_query_times_by_size[t*num_query_sizes.size() + i]);
      }
      std::sort(curr_unsorted_query_times.begin(), curr_unsorted_query_times.end());
      std::sort(curr_sorted_query_times.begin(), curr_sorted_query_times.end());
      outfile << "leafds, " << internal_bytes << ", " <<  header_size << ", " << block_size << ", " <<  max_size << ", " <<  NUM_QUERIES << ", " <<  num_query_sizes[i] << ", " <<  curr_unsorted_query_times[trials/2] << ", " <<  curr_sorted_query_times[trials/2] << ", \n";
    }
    outfile.close();
  }
  return true;
}

template <class T, uint32_t internal_bytes>
bool
test_iterator_merge_map(uint64_t max_size, std::seed_seq &seed, bool write_csv, int num_trials) {
  // std::vector<uint32_t> num_query_sizes{50};
  // std::vector<uint32_t> num_query_sizes{100, 1000, 10000, 100000};

  uint64_t start_time, end_time;
  std::vector<uint64_t> insert_times;
  // std::vector<uint64_t> find_times;

  printf("\nRunning leafds btree with internal bytes = %u with leafds slots %lu\n",internal_bytes, SLOTS);

  // std::vector<uint32_t> num_query_sizes{100, 1000, 10000, 100000};

  // std::vector<T> data =
  //     create_random_data<T>(max_size, max_size * 2, seed);
  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);
  
  std::set<T> checker_set;
  bool wrong;
  tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
                std::allocator<T>, true> concurrent_map1, concurrent_map2;

  // TIME INSERTS
  start_time = get_usecs();
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_map1.insert({data[i], 2*data[i]});
  }
  end_time = get_usecs();
  printf("\tDone inserting %lu elts in %lu\n",max_size, end_time - start_time);

#if CORRECTNESS
  for (uint32_t i = 0; i < max_size; i++) {
    checker_set.insert(data[i]);
  }

  std::vector<T> elts_sorted;
  for (auto key: checker_set) {
    elts_sorted.push_back(key);
  }
  std::sort(elts_sorted.begin(), elts_sorted.end());
  auto it_correct = elts_sorted.begin();
  auto it_leafds = concurrent_map1.begin();
  int count = 0;

  // check iterator correctness
  while (it_correct != elts_sorted.end() && it_leafds != concurrent_map1.end()) {
    T correct_key = *it_correct;
    T leafds_key = it_leafds.key();
    // auto leafds_key_deref = *it_leafds;

    if (correct_key != leafds_key) {
      printf("wrong iterator value, expected %lu but got %lu on count = %lu\n", correct_key, leafds_key, count);
      return false;
    }
    ++it_correct;
    ++it_leafds;
    count++;
  }
  if (it_correct != elts_sorted.end()) {
    printf("leafds iterator counted too few elts\n");
    return false;
  }
  if (it_leafds != concurrent_map1.end()) {
    printf("leafds iterator counted too many elts\n");
    return false;
  } 
  printf("\tcorrect iterator\n");
#endif
  std::seed_seq seed2{1};
  std::vector<T> data2 =
          create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed2);

  cilk_for(uint32_t i = 0; i < max_size; i++) {
      concurrent_map2.insert({data2[i], 2*data2[i]});
  }
  
#if CORRECTNESS
  for (uint32_t i = 0; i < max_size; i++) {
    checker_set.insert(data2[i]);
  }

  std::vector<T> elts_sorted_merged;
  for (auto key: checker_set) {
    elts_sorted_merged.push_back(key);
  }
  std::sort(elts_sorted_merged.begin(), elts_sorted_merged.end());
#endif

  typedef tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T>,
          std::allocator<T>, true> btree_type;

  for(int i = 0; i < num_trials; i++) {
      tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T>,
              std::allocator<T>, true> merged_tree;
      typename btree_type::iterator iterator1 = concurrent_map1.begin();
      typename btree_type::iterator iterator2 = concurrent_map2.begin();

      uint64_t start_time, end_time;
      start_time = get_usecs();
      while (iterator1 != concurrent_map1.end() && iterator2 != concurrent_map2.end()) {
          auto val1 = *iterator1;
          auto val2 = *iterator2;
          if (std::get<0>(val1) < std::get<0>(val2)) {
              merged_tree.insert({std::get<0>(val1), std::get<1>(val1)});
              iterator1++;
          } else {
              merged_tree.insert({std::get<0>(val2), std::get<1>(val2)});
              iterator2++;
          }
      }
      if (iterator1 != concurrent_map1.end()) {
          auto val1 = *iterator1;
          while (iterator1 != concurrent_map1.end()) {
              merged_tree.insert({std::get<0>(val1), std::get<1>(val1)});
              iterator1++;
          }
      } else if (iterator2 != concurrent_map2.end()) {
          auto val2 = *iterator2;
          while (iterator2 != concurrent_map2.end()) {
              merged_tree.insert({std::get<0>(val2), std::get<1>(val2)});
              iterator2++;
          }
      } else {
          std::cout << "Invalid scenario!\n";
          exit(0);
      }
      end_time = get_usecs();
      printf("\tDone merging %lu elts in %lu\n", max_size, end_time - start_time);
  }
  return true;
}

template <class T, uint32_t internal_bytes>
bool
test_iterator_merge_range_version_map(uint64_t max_size, std::seed_seq &seed, bool write_csv, int num_trials) {

  uint64_t start_time, end_time;
  std::vector<uint64_t> insert_times;
  // std::vector<uint64_t> find_times;

  printf("\nRunning leafds btree with internal bytes = %u with leafds slots %lu\n",internal_bytes, SLOTS);

  // std::vector<uint32_t> num_query_sizes{100, 1000, 10000, 100000};

  // std::vector<T> data =
  //     create_random_data<T>(max_size, max_size * 2, seed);
  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);
  
  std::set<T> checker_set;
  bool wrong;
  tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
                std::allocator<T>, true> concurrent_map1, concurrent_map2;

  // TIME INSERTS
  start_time = get_usecs();
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_map1.insert({data[i], 2*data[i]});
  }
  end_time = get_usecs();
  printf("\tDone inserting %lu elts in %lu\n",max_size, end_time - start_time);

#if CORRECTNESS
  for (uint32_t i = 0; i < max_size; i++) {
    checker_set.insert(data[i]);
  }

  std::vector<T> elts_sorted;
  for (auto key: checker_set) {
    elts_sorted.push_back(key);
  }
  std::sort(elts_sorted.begin(), elts_sorted.end());
#endif

  std::seed_seq seed2{1};
  std::vector<T> data2 =
          create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed2);

  cilk_for(uint32_t i = 0; i < max_size; i++) {
      concurrent_map2.insert({data2[i], 2*data2[i]});
  }
  
#if CORRECTNESS
  for (uint32_t i = 0; i < max_size; i++) {
    checker_set.insert(data2[i]);
  }

  std::vector<T> elts_sorted_merged;
  for (auto key: checker_set) {
    elts_sorted_merged.push_back(key);
  }
  std::sort(elts_sorted_merged.begin(), elts_sorted_merged.end());
#endif

  // typedef tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
  //         std::allocator<T>, true> btree_type;

  for(int trial = 0; trial < num_trials; trial++) {
      tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
              std::allocator<T>, true> merged_tree, merged_unsorted_tree;

      std::vector<std::tuple<T, T>> vec1(max_size);
      std::vector<std::tuple<T, T>> vec2(max_size);

      uint64_t count_elts1 = 0;
      uint64_t count_elts2 = 0;

      uint64_t start_time, end_time, end_merge_time;
      start_time = get_usecs();

      concurrent_map1.map_range_length(1, max_size, [&vec1, &count_elts1]([[maybe_unused]] auto key, auto val) {
                  vec1[count_elts1] = {key, val};
                  count_elts1++;
                });

      concurrent_map2.map_range_length(1, max_size, [&vec2, &count_elts2]([[maybe_unused]] auto key, auto val) {
                  vec2[count_elts2] = {key, val};
                  count_elts2++;
                });
      end_time = get_usecs();
      printf("\tDone sweeping %lu elts via sorted range in %lu\n", max_size, end_time - start_time);
      
      // std::mt19937_64 eng(0);
      // std::shuffle(std::begin(vec1), std::end(vec1), eng);
      // std::shuffle(std::begin(vec2), std::end(vec2), eng);

      start_time = get_usecs();
      cilk_for (int i = 0; i < vec1.size(); i++) {
        merged_unsorted_tree.insert({std::get<0>(vec1[i]), std::get<1>(vec1[i])});
      }
      cilk_for (int i = 0; i < vec2.size(); i++) {
        merged_unsorted_tree.insert({std::get<0>(vec2[i]), std::get<1>(vec2[i])});
      }

      end_merge_time = get_usecs();
      printf("\t\tDone concurrent merging sorted %lu elts via sorted_range output in %lu\n", max_size, end_merge_time - start_time);

  }
  return true;
}


template <class T, uint32_t internal_bytes>
bool
test_parallel_iter_merge_map(uint64_t max_size, uint64_t num_chunk_multiplier, std::seed_seq &seed, bool write_csv, int num_trials) {
  uint64_t start_time, end_time;

  uint64_t num_chunks = 48 * num_chunk_multiplier;
  uint64_t chunk_size = std::numeric_limits<T>::max() / num_chunks;

  printf("\nRunning leafds btree with internal bytes = %u with leafds slots %d\n",internal_bytes, SLOTS);

  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);
  
  std::set<T> checker_set;
  bool wrong;
  tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
                std::allocator<T>, true> concurrent_map1, concurrent_map2;

  // TIME INSERTS
  start_time = get_usecs();
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_map1.insert({data[i], 2*data[i]});
  }
  end_time = get_usecs();
  printf("\tDone inserting %lu elts in %lu\n",max_size, end_time - start_time);

#if CORRECTNESS
  for (uint32_t i = 0; i < max_size; i++) {
    checker_set.insert(data[i]);
  }

  std::vector<T> elts_sorted;
  for (auto key: checker_set) {
    elts_sorted.push_back(key);
  }
  std::sort(elts_sorted.begin(), elts_sorted.end());
#endif

  std::seed_seq seed2{1};
  std::vector<T> data2 =
          create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed2);

  start_time = get_usecs();
  cilk_for(uint32_t i = 0; i < max_size; i++) {
      concurrent_map2.insert({data2[i], 2*data2[i]});
  }
  end_time = get_usecs();

  printf("\tDone inserting second %lu elts in %lu\n",max_size, end_time - start_time);
  
#if CORRECTNESS
  for (uint32_t i = 0; i < max_size; i++) {
    checker_set.insert(data2[i]);
  }

  std::vector<std::tuple<T,T>> elts_sorted_merged;
  for (auto key: checker_set) {
    elts_sorted_merged.push_back({key, 2 * key});
  }
  std::sort(elts_sorted_merged.begin(), elts_sorted_merged.end());
  printf("\tDone inserting %lu elts for correctness %lu\n",max_size);
#endif

  // typedef tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
  //        std::allocator<T>, true> btree_type;

  for(int trial = 0; trial <= num_trials; trial++) {
    // tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
    //         std::allocator<T>, true> merged_tree;
    
    std::vector<std::vector<std::tuple<T, T>>> merged_vecs(num_chunks);
    std::vector<uint64_t> merged_vecs_prefix_sums(num_chunks + 1);
    merged_vecs_prefix_sums[0] = 0;
    
    start_time = get_usecs();

    cilk_for (uint64_t chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {
      
      uint64_t start_key = 1 + chunk_idx * chunk_size;
      uint64_t end_key = (chunk_idx == num_chunks - 1) ? std::numeric_limits<T>::max() : 1 + (chunk_idx + 1) * chunk_size;
      std::vector<std::tuple<T, T>> merged_chunk;

      uint64_t start_time, end_time, end_merge_time;
      start_time = get_usecs();

      auto it_btree_1 = concurrent_map1.lower_bound(start_key);
      auto it_btree_2 = concurrent_map2.lower_bound(start_key);
      while (it_btree_1 != concurrent_map1.end() && it_btree_2 != concurrent_map2.end() && std::get<0>(*it_btree_1) < end_key && std::get<0>(*it_btree_2) < end_key) {
        // printf("\t\t lower start_key_1 = %lu, lower start_key_2 = %lu \n", (*it_btree_1).first, (*it_btree_2).first);
        if (std::get<0>(*it_btree_1) < std::get<0>(*it_btree_2)) {
          merged_chunk.push_back(*it_btree_1);
          it_btree_1++;
        }
        else if (std::get<0>(*it_btree_2) < std::get<0>(*it_btree_1)) {
          merged_chunk.push_back(*it_btree_2);
          it_btree_2++;
        }
        else {
          merged_chunk.push_back(*it_btree_1);
          it_btree_1++;
          it_btree_2++;
        }
      }
      while (it_btree_1 != concurrent_map1.end() && std::get<0>(*it_btree_1) < end_key) {
        merged_chunk.push_back(*it_btree_1);
        it_btree_1++;
      }
      while (it_btree_2 != concurrent_map2.end() && std::get<0>(*it_btree_2) < end_key) {
        merged_chunk.push_back(*it_btree_2);
        it_btree_2++;
      }

      // std::merge(concurrent_map1.lower_bound(start_key), concurrent_map1.lower_bound(end_key), concurrent_map2.lower_bound(start_key), concurrent_map2.lower_bound(end_key), std::back_inserter(merged_chunk));
      merged_vecs[chunk_idx] = merged_chunk;
      merged_vecs_prefix_sums[chunk_idx + 1] = merged_chunk.size();
    }

    end_time = get_usecs();
    printf("\tDone merging chunks for %lu elts via sorted range end in %lu\n", max_size, end_time - start_time);


    start_time = get_usecs();
    for (size_t chunk_idx = 1; chunk_idx < num_chunks + 1; chunk_idx++) {
      merged_vecs_prefix_sums[chunk_idx] += merged_vecs_prefix_sums[chunk_idx - 1];
    }
    end_time = get_usecs();
    uint64_t summing_time = end_time - start_time;
    
    std::vector<std::tuple<T, T>> concat_merged(merged_vecs_prefix_sums[num_chunks]);

    start_time = get_usecs();
    cilk_for (size_t chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {
      for (uint64_t index = merged_vecs_prefix_sums[chunk_idx]; index < merged_vecs_prefix_sums[chunk_idx + 1]; index++) {
        concat_merged[index] = merged_vecs[chunk_idx][index - merged_vecs_prefix_sums[chunk_idx]];
      }
    }
    end_time = get_usecs();
    printf("\tDone concating chunks for %lu elts via sorted range end in %lu\n", max_size, summing_time + end_time - start_time);

#if CORRECTNESS
    if (elts_sorted_merged != concat_merged) {
      printf("\tMerged vector not correct\n");
      for (uint64_t i = 0; i < elts_sorted_merged.size(); i++) {
        if (elts_sorted_merged[i] != concat_merged[i]) {
          printf("wrong at index %lu , correct = %lu, wrong = %lu\n", i, std::get<0>(elts_sorted_merged[i]), std::get<0>(concat_merged[i]));
        }
      }
      return false;
    }
    printf("\tMerged vector correct\n");
#endif
  }
  return true;
}

template <class T, uint32_t internal_bytes>
bool
test_parallel_merge_map(uint64_t max_size, uint64_t num_chunk_multiplier, std::seed_seq &seed, bool write_csv, int num_trials) {
  uint64_t start_time, end_time;

  uint64_t num_chunks = 48 * num_chunk_multiplier;
  uint64_t chunk_size = std::numeric_limits<T>::max() / num_chunks;

  printf("\nRunning leafds btree with internal bytes = %u with leafds slots %lu\n",internal_bytes, SLOTS);

  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);
  
  std::set<T> checker_set;
  bool wrong;
  tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
                std::allocator<T>, true> concurrent_map1, concurrent_map2;

  // TIME INSERTS
  start_time = get_usecs();
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_map1.insert({data[i], 2*data[i]});
  }
  end_time = get_usecs();
  printf("\tDone inserting %lu elts in %lu\n",max_size, end_time - start_time);

#if CORRECTNESS
  for (uint32_t i = 0; i < max_size; i++) {
    checker_set.insert(data[i]);
  }

  std::vector<T> elts_sorted;
  for (auto key: checker_set) {
    elts_sorted.push_back(key);
  }
  std::sort(elts_sorted.begin(), elts_sorted.end());
#endif

  std::seed_seq seed2{1};
  std::vector<T> data2 =
          create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed2);

  cilk_for(uint32_t i = 0; i < max_size; i++) {
      concurrent_map2.insert({data2[i], 2*data2[i]});
  }
  
#if CORRECTNESS
  for (uint32_t i = 0; i < max_size; i++) {
    checker_set.insert(data2[i]);
  }

  std::vector<std::tuple<T,T>> elts_sorted_merged;
  for (auto key: checker_set) {
    elts_sorted_merged.push_back({key, 2 * key});
  }
  std::sort(elts_sorted_merged.begin(), elts_sorted_merged.end());
#endif

  // typedef tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
  //         std::allocator<T>, true> btree_type;

  for(int trial = 0; trial <= num_trials; trial++) {
    // tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
    //         std::allocator<T>, true> merged_tree;
    
    // std::vector<std::vector<std::tuple<T, T>>*> merged_vecs(num_chunks);
    std::vector<std::vector<std::tuple<T, T>>> merged_vecs(num_chunks);
    std::vector<uint64_t> merged_vecs_prefix_sums(num_chunks + 1);
    merged_vecs_prefix_sums[0] = 0;
    
    start_time = get_usecs();

    cilk_for (uint64_t chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {
      std::vector<std::tuple<T, T>> vec1, vec2;
      uint64_t start_key = 1 + chunk_idx * chunk_size;
      uint64_t end_key = (chunk_idx == num_chunks - 1) ? std::numeric_limits<T>::max() : 1 + (chunk_idx + 1) * chunk_size;

      concurrent_map1.map_range_sorted_end(start_key, end_key, [&vec1](auto key, auto val) {
                  vec1.emplace_back(std::tuple<T, T>{key, val});
                });
      concurrent_map2.map_range_sorted_end(start_key, end_key, [&vec2](auto key, auto val) {
                  vec2.emplace_back(std::tuple<T, T>{key, val});
                });
      std::vector<std::tuple<T, T>> merged_chunk;

      std::merge(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), std::back_inserter(merged_chunk));
      merged_vecs[chunk_idx] = merged_chunk;
      // merged_vecs[chunk_idx] = &merged_chunk;
      merged_vecs_prefix_sums[chunk_idx + 1] = merged_chunk.size();
    }

    end_time = get_usecs();
    printf("\tDone merging chunks for %lu elts via sorted range end in %lu\n", max_size, end_time - start_time);


    start_time = get_usecs();
    for (size_t chunk_idx = 1; chunk_idx < num_chunks + 1; chunk_idx++) {
      merged_vecs_prefix_sums[chunk_idx] += merged_vecs_prefix_sums[chunk_idx - 1];
    }
    end_time = get_usecs();
    uint64_t summing_time = end_time - start_time;
    
    std::vector<std::tuple<T, T>> concat_merged(merged_vecs_prefix_sums[num_chunks]);

    start_time = get_usecs();
    cilk_for (size_t chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {
      for (uint64_t index = merged_vecs_prefix_sums[chunk_idx]; index < merged_vecs_prefix_sums[chunk_idx + 1]; index++) {
        // std::vector<std::tuple<T, T>>& vecRef = *merged_vecs[chunk_idx];
        // printf("chunk = %lu, vec size = %lu, real size = %lu\n", chunk_idx, vecRef.size(), merged_vecs_prefix_sums[chunk_idx + 1] - merged_vecs_prefix_sums[chunk_idx]);
        // concat_merged[index] = vecRef[index - merged_vecs_prefix_sums[chunk_idx]];
        concat_merged[index] = merged_vecs[chunk_idx][index - merged_vecs_prefix_sums[chunk_idx]];
      }
    }
    end_time = get_usecs();
    printf("\tDone concating chunks for %lu elts via sorted range end in %lu\n", max_size, summing_time + end_time - start_time);

#if CORRECTNESS
    if (elts_sorted_merged != concat_merged) {
      printf("\tMerged vector not correct\n");
      for (uint64_t i = 0; i < elts_sorted_merged.size(); i++) {
        if (elts_sorted_merged[i] != concat_merged[i]) {
          printf("wrong at index %lu , correct = %lu, wrong = %lu\n", i, std::get<0>(elts_sorted_merged[i]), std::get<0>(concat_merged[i]));
        }
      }
      return false;
    }
#endif
  }
  return true;
}

template <class T, uint32_t internal_bytes>
bool
test_bulk_load_map(uint64_t max_size, std::seed_seq &seed, bool write_csv, int num_trials) {
  uint64_t start_time, end_time;
  std::vector<uint64_t> insert_times;
  std::vector<uint64_t> find_times;

  std::vector<T> data = create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);
  std::vector<std::tuple<T, T>> bulk_data;

  std::set<T> checker_set;
  bool wrong;
  for (uint32_t i = 0; i < max_size; i++) {
    checker_set.insert(data[i]);
  }
  for (auto e : checker_set) {
    bulk_data.push_back({e, 2*e});
  }
  std::sort(bulk_data.begin(), bulk_data.end());

  for (int cur_trial = 0; cur_trial <= num_trials; cur_trial++) {
    printf("\nRunning leafds btree with internal bytes = %u with leafds slots %lu , trial = %lu\n",internal_bytes, SLOTS, cur_trial);

  #if CORRECTNESS

    std::vector<T> checker_sorted;
    for (auto key: checker_set) {
      checker_sorted.push_back(key);
    }
    std::sort(checker_sorted.begin(), checker_sorted.end());
  #endif

    // output to tree_type, internal bytes, leaf bytes, num_inserted, insert_time, num_finds, find_time, num_range_queries, range_time_maxlen_{}*

    tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
                  std::allocator<T>, true> concurrent_map, concurrent_map_bulk_load;

    // TIME INSERTS
    start_time = get_usecs();
    cilk_for(uint32_t i = 0; i < max_size; i++) {
      concurrent_map.insert({data[i], 2*data[i]});
    }
    end_time = get_usecs();
    // if (cur_trial > 0) {insert_times.push_back(end_time - start_time);}
    printf("\tDone inserting %lu elts in %lu\n",max_size, end_time - start_time);

    // TIME BULK LOAD
    start_time = get_usecs();
    concurrent_map_bulk_load.bulk_load(bulk_data.begin(), bulk_data.end());
    end_time = get_usecs();
    // if (cur_trial > 0) {insert_times.push_back(end_time - start_time);}
    printf("\tDone bulk loading %lu elts in %lu\n",max_size, end_time - start_time);

#if CORRECTNESS
    std::vector<int> found_count(max_size);
    std::vector<int> found_count_bulk(max_size);

    T sum_in_tree = 0;
    uint64_t num_in_tree = 0;
    concurrent_map.map_range_length(1, max_size, [&sum_in_tree, &num_in_tree]([[maybe_unused]] auto key, auto val) {
                  sum_in_tree += key;
                  num_in_tree++;
    });
    T sum_in_tree_bulk = 0;
    uint64_t num_in_tree_bulk = 0;
    concurrent_map_bulk_load.map_range_length(1, max_size, [&sum_in_tree_bulk, &num_in_tree_bulk]([[maybe_unused]] auto key, auto val) {
                  sum_in_tree_bulk += key;
                  num_in_tree_bulk++;
    });
    printf("Did map range, sum regular inserts = %lu elts, bulk load = %lu elts\n", sum_in_tree, sum_in_tree_bulk);
    printf("Did map range, num regular inserts = %lu elts, bulk load = %lu elts\n", num_in_tree, num_in_tree_bulk);

    cilk_for(uint32_t i = 0; i < max_size; i++) {
      found_count[i] = concurrent_map.exists(data[i]) ? 1 : 0;
      if (!found_count[i]) {
        printf("Elt doesn't exist: key = %lu\n", data[i]);
      }
    }
    cilk_for(uint32_t i = 0; i < max_size; i++) {
      found_count_bulk[i] = concurrent_map_bulk_load.exists(data[i]) ? 1 : 0;
      if (!found_count_bulk[i]) {
        printf("Elt doesn't exist bulk: key = %lu\n", data[i]);
      }
    }

    for (uint64_t i = 0; i < max_size; i ++) {
      if (!found_count[i]) {
        printf("Elt doesn't exist serial check: key = %lu\n", data[i]);
      }
      if (!found_count_bulk[i]) {
        printf("Elt doesn't exist bulk serial check: key = %lu\n", data[i]);
      }
    }
    int count_found = 0;
    int count_found_bulk = 0;
    for (auto e : found_count) {
      count_found += e ? 1 : 0;
    }
    for (auto e : found_count_bulk) {
      count_found_bulk += e ? 1 : 0;
    }
    if (count_found != count_found_bulk) {
      printf("Did not find some elements, regular inserts = %lu elts, bulk load = %lu elts\n", count_found, count_found_bulk);
      return false;
    } 
    printf("Correct bulk load\n");
    return true;
#endif
  }
  return true;
}

int main(int argc, char *argv[]) {

  cxxopts::Options options("BtreeTester",
                           "allows testing different attributes of the btree");

  options.positional_help("Help Text");

  // clang-format off
  options.add_options()
    ("trials", "how many values to insert", cxxopts::value<int>()->default_value( "5"))
    ("num_inserts", "number of values to insert", cxxopts::value<int>()->default_value( "100000000"))
    ("num_queries", "number of queries for query tests", cxxopts::value<int>()->default_value( "1000000"))
    ("query_size", "query size for cache test", cxxopts::value<int>()->default_value( "10000"))
    ("num_chunks", "number of chunks for merge tests", cxxopts::value<int>()->default_value( "480"))
    ("write_csv", "whether to write timings to disk")
    ("serial_microbenchmark_leafds", "run leafds 1024 byte btree microbenchmark with [trials] [num_inserts] [num_queries] [write_csv]")
    ("microbenchmark_leafds", "run leafds 1024 byte btree microbenchmark with [trials] [num_inserts] [num_queries] [write_csv]")
    ("sequential_inserts", "run parallel inserts where each threads inserts sequential elements")
    ("psum_leafds", "run leafds 1024 byte btree psum with [trials] [num_inserts]")
    ("cache_misses", "insert and query")
    ("merge_iter", "run baseline 1024 btree merge with iterators with [trials] [num_inserts per tree]");
    
  std::seed_seq seed{0};
  auto result = options.parse(argc, argv);
  uint32_t trials = result["trials"].as<int>();
  uint32_t num_inserts = result["num_inserts"].as<int>();
  uint32_t num_queries = result["num_queries"].as<int>();
  uint32_t num_chunks = result["num_chunks"].as<int>();
  uint32_t query_size = result["query_size"].as<int>();
  uint32_t write_csv = result["write_csv"].as<bool>();

  std::ofstream outfile;
  outfile.open("insert_finds.csv", std::ios_base::app); 
  outfile << "tree_type, internal bytes, leaf slots, num_inserted, insert_time, num_finds, find_time, \n";
  outfile.close();
  outfile.open("range_queries.csv", std::ios_base::app); 
  outfile << "tree_type, internal bytes, leaf slots, num_inserted,num_range_queries, max_query_size,  unsorted_query_time, sorted_query_time, \n";
  outfile.close();

  if (result["serial_microbenchmark_leafds"].as<bool>()) {
    bool correct = test_serial_microbenchmarks_map<unsigned long, 1024>(num_inserts, num_queries, seed, write_csv, trials);
    return !correct;
  }

  if (result["microbenchmark_leafds"].as<bool>()) {
    bool correct = test_concurrent_microbenchmarks_map<unsigned long, 1024, 4, 4>(num_inserts, num_queries, seed, write_csv, trials);

    correct |= test_concurrent_microbenchmarks_map<unsigned long, 1024, 4, 8>(num_inserts, num_queries, seed, write_csv, trials);
    correct |= test_concurrent_microbenchmarks_map<unsigned long, 1024, 8, 8>(num_inserts, num_queries, seed, write_csv, trials);
    correct |= test_concurrent_microbenchmarks_map<unsigned long, 1024, 8, 16>(num_inserts, num_queries, seed, write_csv, trials);
    correct |= test_concurrent_microbenchmarks_map<unsigned long, 1024, 16, 16>(num_inserts, num_queries, seed, write_csv, trials);
    correct |= test_concurrent_microbenchmarks_map<unsigned long, 1024, 16, 32>(num_inserts, num_queries, seed, write_csv, trials);
    correct |= test_concurrent_microbenchmarks_map<unsigned long, 1024, 32, 32>(num_inserts, num_queries, seed, write_csv, trials);
    correct |= test_concurrent_microbenchmarks_map<unsigned long, 1024, 32, 64>(num_inserts, num_queries, seed, write_csv, trials);
    correct |= test_concurrent_microbenchmarks_map<unsigned long, 1024, 64, 64>(num_inserts, num_queries, seed, write_csv, trials);


    return !correct;
  }
  if (result["sequential_inserts"].as<bool>()) {
    return test_sequential_inserts<unsigned long, 1024>(num_inserts, trials);
  }
  if (result["psum_leafds"].as<bool>()) {
    test_concurrent_sum_time<unsigned long, 1024>(num_inserts, seed, trials);
    return 0;
  }
  if (result["merge_iter"].as<bool>()) {
    test_parallel_iter_merge_map<unsigned long, 1024>(num_inserts, num_chunks, seed, write_csv, trials);
    return 0;
  }
  if (result["cache_misses"].as<bool>()) {
    test_map_cache_misses<unsigned long, 1024, 32, 32>(num_inserts, num_queries, query_size, seed, trials);
  }

  // bool correct = test_bulk_load_map<unsigned long, 1024>(n, seed, write_csv, trials);
  // bool correct = test_parallel_merge_map<unsigned long, 1024>(n, num_queries, seed, write_csv, trials);
  // bool correct = test_parallel_iter_merge_map<unsigned long, 1024>(n, num_queries, seed, write_csv, trials);
  // bool correct = test_iterator_merge_range_version_map<unsigned long, 1024>(n, seed, write_csv, trials);
  // bool correct = test_iterator_merge_map<unsigned long, 1024>(n, seed, write_csv, trials);
  // test_concurrent_btreeset_scalability<unsigned long, 1024>(n, num_queries, seed, trials);
  // bool correct = test_concurrent_microbenchmarks_map<unsigned long, 1024>(n, num_queries, seed, write_csv, trials);

  // if (!correct) {
  //   printf("got the wrong answer :(\n");
  //   return -1;
  // }

  return 0;
}
