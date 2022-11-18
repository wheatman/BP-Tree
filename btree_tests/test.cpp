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

#include <tlx/container/btree_set.hpp>
#include "tlx/container/btree_map.hpp"

#if CILK != 1
#define cilk_for for
#endif

#define CORRECTNESS 0

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

template <class T>
std::tuple<bool, uint64_t, uint64_t, uint64_t, uint64_t>
test_concurrent_btreeset(uint64_t max_size, std::seed_seq &seed) {
  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);

  uint64_t start, end, serial_time;

// #if DEBUG
  tlx::btree_set<T> serial_set;
  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    serial_set.insert(data[i]);
  }
  end = get_usecs();
  serial_time = end - start;
  printf("inserted all the data serially in %lu\n", end - start);
// #endif
  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, false>
      serial_test_set;
  start = get_usecs();
  for(uint32_t i = 0; i < max_size; i++) {
    serial_test_set.insert(data[i]);
  }
  end = get_usecs();
  serial_time = end - start;
  printf("\tinserted %lu elts serially in %lu\n", max_size, serial_time);

  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, true> concurrent_set;
  start = get_usecs();
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_set.insert(data[i]);
  }
  end = get_usecs();
  uint64_t parallel_time = end - start;
  printf("\tinserted %lu elts concurrently in %lu\n", max_size, parallel_time);

// #if DEBUG
  if (serial_set.size() != concurrent_set.size()) {
    printf("the sizes don't match, got %lu, expetected %lu\n",
           concurrent_set.size(), serial_set.size());
    return {false, 0, 0, 0, 0};
  }
  auto it_serial = serial_set.begin();
  auto it_concurrent = concurrent_set.begin();
  bool wrong = false;
  for (uint64_t i = 0; i < serial_set.size(); i++) {
    if (*it_serial != *it_concurrent) {
      printf("don't match in position %lu, got %lu, expected %lu\n", i,
             *it_concurrent, *it_serial);
      wrong = true;
    }
    it_serial++;
    it_concurrent++;
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  }

  std::vector<uint64_t> indxs_to_remove =
      create_random_data<uint64_t>(max_size / 2, data.size(), seed);

  uint64_t serial_remove_start = get_usecs();
  for (uint32_t i = 0; i < indxs_to_remove.size(); i++) {
    serial_set.erase(data[indxs_to_remove[i]]);
  }
  uint64_t serial_remove_end = get_usecs();
  uint64_t serial_remove_time = serial_remove_end - serial_remove_start;

  uint64_t parallel_remove_start = get_usecs();
  cilk_for(uint32_t i = 0; i < indxs_to_remove.size(); i++) {
    concurrent_set.erase(data[indxs_to_remove[i]]);
  }
  // return {true, serial_time, parallel_time, serial_remove_time, parallel_remove_time};
// #endif

  uint64_t parallel_remove_end = get_usecs();
  uint64_t parallel_remove_time = parallel_remove_end - parallel_remove_start;

  printf("removed half the data serially in %lu\n",
         serial_remove_end - serial_remove_start);

  printf("removed half the data in parallel in %lu\n",
         parallel_remove_end - parallel_remove_start);

  if (serial_set.size() != concurrent_set.size()) {
    printf("the sizes don't match, got %lu, expetected %lu\n",
           concurrent_set.size(), serial_set.size());
    return {false, 0, 0, 0, 0};
  }
  printf("inserted %lu elements\n", serial_set.size());
  it_serial = serial_set.begin();
  it_concurrent = concurrent_set.begin();
  wrong = false;
  for (uint64_t i = 0; i < serial_set.size(); i++) {
    if (*it_serial != *it_concurrent) {
      printf("don't match in position %lu, got %lu, expected %lu\n", i,
             *it_concurrent, *it_serial);
      wrong = true;
    }
    it_serial++;
    it_concurrent++;
  }
  if (wrong) {
    return {false, 0, 0, 0, 0};
  }

  return {true, serial_time, parallel_time, serial_remove_time,
          parallel_remove_time};
}

template <class T>
void test_concurrent_find(uint64_t max_size, std::seed_seq &seed) {
  std::vector<T> data =
      create_random_data<T>(max_size * 2, std::numeric_limits<T>::max(), seed);

  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, true>
      concurrent_set;
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_set.insert(data[i]);
  }
  printf("have %lu elements\n", concurrent_set.size());
  ParallelTools::Reducer_sum<uint64_t> found;
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    found += concurrent_set.exists(data[i]);
    concurrent_set.insert(data[i + max_size]);
  }
  printf("found %lu elements\n", found.get());
}

template <class T>
std::tuple<bool, uint64_t, uint64_t, uint64_t, uint64_t>
test_concurrent_range_query(uint64_t max_size, std::seed_seq &seed) {
  uint64_t NUM_QUERIES = 10000000;
  uint64_t MAX_QUERY_SIZE = 100000;

  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);


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
test_concurrent_range_query_old(uint64_t max_size, std::seed_seq &seed) {
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
    serial_set.insert(data[i]);
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

  // std::vector<T> checker_cache_sums;
  // T curr_sum = 0;
  // for (auto e: checker_sorted) {
  //   curr_sum += e;
  //   checker_cache_sums.push_back(curr_sum);
  // }
  printf("Done inserting %lu elts\n",max_size);

  uint64_t start_time, end_time, serial_time, concurrent_time;

  std::vector<uint64_t> range_query_idxs =
      create_random_data<uint64_t>(checker_sorted.size(), checker_sorted.size(), seed);
  
  // std::vector<uint64_t> range_query_end_idxs =
  //     create_random_data<uint64_t>(checker_sorted.size() / 2, checker_sorted.size(), seed);

  std::vector<T> correct_range_query_sums(checker_sorted.size()/2);
  std::vector<uint64_t> correct_range_query_counts(checker_sorted.size()/2);

  std::vector<T> serial_range_query_sums(checker_sorted.size()/2);
  std::vector<uint64_t> serial_range_query_counts(checker_sorted.size()/2);

  std::vector<T> concurrent_range_query_sums(checker_sorted.size()/2);
  std::vector<uint64_t> concurrent_range_query_counts(checker_sorted.size()/2);

  T start, end;

  // get correct range sums
  for (uint32_t i = 0; i < checker_sorted.size()/2; i++) {
    start = data[range_query_idxs[i]];
    end = data[range_query_idxs[i + checker_sorted.size()/2]];
    if (start > end) {
      std::swap(start, end);
    }
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
  printf("\t did %lu correct range queries with avg query size %lu\n", 
        checker_sorted.size()/2,
        std::accumulate(correct_range_query_counts.begin(), correct_range_query_counts.end(), 0.0) / correct_range_query_counts.size());

  start_time = get_usecs();
  // serial btree range sums
  for (uint32_t i = 0; i < checker_sorted.size()/2; i++) {
    start = data[range_query_idxs[i]];
    end = data[range_query_idxs[i + checker_sorted.size()/2]];
    if (start > end) {
      std::swap(start, end);
    }
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
  printf("\t did %lu range queries serially in %lu\n", checker_sorted.size()/2, serial_time);

  // correctness check of serial 
  bool wrong = false;
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

  start_time = get_usecs();
  // concurrent btree range sums
  cilk_for (uint32_t i = 0; i < checker_sorted.size()/2; i++) {
    start = data[range_query_idxs[i]];
    end = data[range_query_idxs[i + checker_sorted.size()/2]];
    if (start > end) {
      std::swap(start, end);
    }
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
  printf("\t did %lu range queries concurrently in %lu\n", checker_sorted.size()/2, concurrent_time);

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
  // return {true, 0, 0, 0, 0};

  // *** testing map_range_length ***

  std::vector<T> serial_range_query_length_sums(checker_sorted.size()/2);
  std::vector<uint64_t> serial_range_query_length_counts(checker_sorted.size()/2);

  std::vector<T> concurrent_range_query_length_sums(checker_sorted.size()/2);
  std::vector<uint64_t> concurrent_range_query_length_counts(checker_sorted.size()/2);

  start_time = get_usecs();
  // serial btree range sums
  for (uint32_t i = 0; i < checker_sorted.size()/2; i++) {
    start = data[range_query_idxs[i]];
    end = data[range_query_idxs[i + checker_sorted.size()/2]];
    if (start > end) {
      std::swap(start, end);
    }
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
  printf("\t did %lu range queries by length serially in %lu\n", checker_sorted.size()/2, serial_time);

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

  start_time = get_usecs();
  // concurrent btree range sums
  cilk_for (uint32_t i = 0; i < checker_sorted.size()/2; i++) {
    start = data[range_query_idxs[i]];
    end = data[range_query_idxs[i + checker_sorted.size()/2]];
    if (start > end) {
      std::swap(start, end);
    }
    uint64_t num_in_range = 0;
    uint64_t sum_in_range = 0;
    concurrent_set.map_range_length(start, correct_range_query_counts[i], [&num_in_range, &sum_in_range]([[maybe_unused]] auto el) {
              num_in_range += 1;
              sum_in_range += el;
            });
    concurrent_range_query_length_counts[i] = num_in_range;
    concurrent_range_query_length_sums[i] = sum_in_range;
  }
  end_time = get_usecs();
  concurrent_time = end_time - start_time;
  printf("\t did %lu range queries by length concurrently in %lu\n", checker_sorted.size()/2, concurrent_time);

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
  return {true, 0, 0, 0, 0};
}

template <class T>
void test_concurrent_sum(uint64_t max_size, std::seed_seq &seed) {
  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);
  std::set<T> serial_set;

  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                    std::allocator<T>, false> serial_test_set;

  for (uint32_t i = 0; i < max_size; i++) {
    serial_set.insert(data[i]);
    serial_test_set.insert(data[i]);
  }
  
  uint64_t correct_sum = 0;
  for(auto e : serial_set) {
	  correct_sum += e;
  }
  auto serial_sum = serial_test_set.psum();
  // printf("serial btree sum got %lu, should be %lu\n", serial_sum, correct_sum);
  assert(serial_sum == correct_sum);

  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, true>
      concurrent_set;
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_set.insert(data[i]);
  }
  auto concurrent_sum = concurrent_set.psum();
  assert(concurrent_sum == correct_sum);
  //printf("concurrent btree sum got %lu, should be %lu\n", concurrent_sum, correct_sum);
}


template <class T>
void test_concurrent_sum_time(uint64_t max_size, std::seed_seq &seed, int trials) {
  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);

  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, true>
      concurrent_set;
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_set.insert(data[i]);
  }
  std::vector<uint64_t> psum_times(trials);
  uint64_t start, end;
  for(int i = 0; i < trials + 1; i++) {
          start = get_usecs();
          auto concurrent_sum = concurrent_set.psum();
          end = get_usecs();

          printf("concurrent btree sum got %lu\n", concurrent_sum);

	  if(i > 0) {
		  psum_times[i-1] = end - start;
	  }
  }
  std::sort(psum_times.begin(), psum_times.end());
  printf("plain btree psum time = %lu\n", psum_times[trials / 2]);
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
    // if (std::count(serial_set.begin(), serial_set.end(), data[i])) {
    //   printf("inserting twice %lu\n", data[i]);
    //   duplicated_elts.insert(data[i]);
    // }
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

  /*
  // test sum
  uint64_t serial_test_sum = serial_test_set.psum();
  ASSERT("serial sum got %lu, should be %lu\n", serial_test_sum, correct_sum);
  uint64_t concurrent_test_sum = concurrent_set.psum();
  ASSERT("concurrent sum got %lu, should be %lu\n", concurrent_sum, correct_sum);
  */
#endif

  return {true, serial_time, parallel_time, 0, 0};
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
    concurrent_map.map_range(start, end, [&num_in_range, &sum_in_range]([[maybe_unused]] auto el) {
              num_in_range += 1;
              sum_in_range += el.second;
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
    concurrent_map.map_range_length(start, correct_range_query_counts[i], [&num_in_range, &sum_in_range]([[maybe_unused]] auto el) {
              num_in_range += 1;
              sum_in_range += el.second;
            });
#else
    concurrent_map.map_range_length(start, concurrent_range_query_counts[i], [&num_in_range, &sum_in_range]([[maybe_unused]] auto el) {
              num_in_range += 1;
              sum_in_range += el.second;
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

template <class T, uint32_t internal_bytes, uint32_t leaf_bytes>
bool
test_concurrent_microbenchmarks_map(uint64_t max_size, uint64_t NUM_QUERIES, std::seed_seq &seed, bool write_csv, int trials) {
  std::vector<uint32_t> num_query_sizes{100, 1000, 10000, 100000};

  uint64_t start_time, end_time;
  std::vector<uint64_t> insert_times;
  std::vector<uint64_t> find_times;
  std::vector<uint64_t> sorted_range_query_times_by_size;
  std::vector<uint64_t> unsorted_range_query_times_by_size;

  for (int cur_trial = 0; cur_trial <= trials; cur_trial++) {
    printf("\nRunning plain btree with internal bytes = %u, leaf bytes = %u, trial = %d\n",internal_bytes, leaf_bytes, cur_trial);

    std::vector<T> data =
        create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);
    
    std::set<T> checker_set;
    bool wrong;

  #if CORRECTNESS
    for (uint32_t i = 0; i < max_size; i++) {
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

    // output to tree_type, internal bytes, leaf bytes, num_inserted, insert_time, num_finds, find_time, num_range_queries, range_time_maxlen_{}*

    tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes, leaf_bytes>,
                  std::allocator<T>, true> concurrent_map;

    // TIME INSERTS
    start_time = get_usecs();
    cilk_for(uint32_t i = 0; i < max_size; i++) {
      concurrent_map.insert({data[i], 2*data[i]});
    }
    end_time = get_usecs();
    if (cur_trial > 0) {insert_times.push_back(end_time - start_time);}
    printf("\tDone inserting %lu elts in %lu\n",max_size, end_time - start_time);
    

    // TIME POINT QUERIES
    start_time = get_usecs();
    cilk_for(uint32_t i = 0; i < NUM_QUERIES; i++) {
      concurrent_map.exists(data[range_query_start_idxs[i]]);
    }
    end_time = get_usecs();
    if (cur_trial > 0) {find_times.push_back(end_time - start_time);}
    printf("\tDone finding %lu elts in %lu\n",NUM_QUERIES, end_time - start_time);

    // TIME RANGE QUERIES FOR ALL LENGTHS
    for (size_t query_size_i = 0; query_size_i < num_query_sizes.size(); query_size_i++) {
      uint64_t MAX_QUERY_SIZE = num_query_sizes[query_size_i];
      std::vector<uint64_t> range_query_lengths =
        create_random_data<uint64_t>(NUM_QUERIES, MAX_QUERY_SIZE - 1, query_seed_2);

      T start, end;
  #if CORRECTNESS
      std::vector<T> correct_range_query_maxs(NUM_QUERIES);
      std::vector<uint64_t> correct_range_query_counts(NUM_QUERIES);

      // get correct range sums
      cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
        start = checker_sorted[range_query_start_idxs[i]];

        correct_range_query_counts[i] = 0;
        correct_range_query_maxs[i] = start;

        size_t curr_index = range_query_start_idxs[i];
        // if (i == 12280) {
        //   printf("problem with len %lu\n", range_query_lengths[i]);
        //   printf("curr index %lu, checker size %lu \n", curr_index, checker_sorted.size());
        // }

        // printf("query len %lu", range_query_lengths[i]);
        for (uint64_t count = 0; count < range_query_lengths[i]; count++) {
          if (curr_index >= checker_sorted.size()) {
            break;
          }
          correct_range_query_maxs[i] = checker_sorted[curr_index];
          correct_range_query_counts[i] += 1;
          // printf("range count %lu", correct_range_query_counts[i] );
          curr_index++;
        }
      }
      printf("\t did %lu correct range queries with max query size %lu\n", 
            NUM_QUERIES,
            MAX_QUERY_SIZE);
  #endif

      // sorted query first to get end key
      std::vector<T> concurrent_range_query_length_maxs(NUM_QUERIES);
      std::vector<uint64_t> concurrent_range_query_length_counts(NUM_QUERIES);

      start_time = get_usecs();
      cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
  #if CORRECTNESS
        start = checker_sorted[range_query_start_idxs[i]];
  #else 
        start = data[range_query_start_idxs[i]];
  #endif

        uint64_t num_in_range = 0;
        T max_key_in_range = start;

        concurrent_map.map_range_length(start, range_query_lengths[i], [&num_in_range, &max_key_in_range]([[maybe_unused]] auto el) {
                  num_in_range++;
                  if (el.first > max_key_in_range) {
                    max_key_in_range = el.first;
                  }
                });
        concurrent_range_query_length_counts[i] = num_in_range;
        concurrent_range_query_length_maxs[i] = max_key_in_range;
      }
      end_time = get_usecs();
      if (cur_trial > 0) {sorted_range_query_times_by_size.push_back(end_time - start_time);}
      printf("\t\t did sorted range queries with max len %lu concurrently in %lu\n", MAX_QUERY_SIZE, end_time - start_time);

  #if CORRECTNESS
      // correctness check of concurrent sorted query
      wrong = false;
      for (size_t i = 0; i < NUM_QUERIES; i++) {
        if (correct_range_query_counts[i] != concurrent_range_query_length_counts[i]) {
          printf("wrong concurrent range query length count at %lu with len %lu, expected %lu, got %lu\n", i, range_query_lengths[i], correct_range_query_counts[i], concurrent_range_query_length_counts[i]);
          wrong = true;
        }
        if (correct_range_query_maxs[i] != concurrent_range_query_length_maxs[i]) {
          printf("wrong concurrent range query length max, expected %lu, got %lu\n", correct_range_query_maxs[i], concurrent_range_query_length_maxs[i]);
          wrong = true;
        }
      }
      if (wrong) {
        return false;
      } 
  #endif

      std::vector<T> concurrent_range_query_maxs(NUM_QUERIES);
      std::vector<uint64_t> concurrent_range_query_counts(NUM_QUERIES);

      start_time = get_usecs();
      cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
  #if CORRECTNESS
        start = checker_sorted[range_query_start_idxs[i]];
  #else 
        start = data[range_query_start_idxs[i]];
  #endif
        end = concurrent_range_query_length_maxs[i];
        if (range_query_lengths[i] != 0) {
          end++;
        }

        uint64_t num_in_range = 0;
        T max_key_in_range = start;
        concurrent_map.map_range(start, end, [&num_in_range, &max_key_in_range]([[maybe_unused]] auto el) {
                  num_in_range += 1;
                  if (el.first > max_key_in_range) {
                    max_key_in_range = el.first;
                  }
                });
        concurrent_range_query_counts[i] = num_in_range;
        concurrent_range_query_maxs[i] = max_key_in_range;
      }
      end_time = get_usecs();
      if (cur_trial > 0) {unsorted_range_query_times_by_size.push_back(end_time - start_time);}
      printf("\t\t did unsorted range queries with max len %lu concurrently in %lu\n", MAX_QUERY_SIZE, end_time - start_time);

  #if CORRECTNESS
      // correctness check of concurrent 
      wrong = false;
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
    outfile << "plain, " << internal_bytes << ", " <<  leaf_bytes << ", " <<  max_size << ", " <<  insert_time_med << ", " <<  NUM_QUERIES << ", " <<  find_time_med << ", \n";
    outfile.close();
    outfile.open("range_queries.csv", std::ios_base::app); 
    for (size_t i = 0; i < num_query_sizes.size(); i ++ ) {
      std::vector<uint64_t> curr_unsorted_query_times;
      std::vector<uint64_t> curr_sorted_query_times;
      for (int t = 0; t < trials; t++) {
        curr_unsorted_query_times.push_back(unsorted_range_query_times_by_size[t*num_query_sizes.size() + i]);
        curr_sorted_query_times.push_back(sorted_range_query_times_by_size[t*num_query_sizes.size() + i]);
      }
      std::sort(curr_unsorted_query_times.begin(), curr_unsorted_query_times.end());
      std::sort(curr_sorted_query_times.begin(), curr_sorted_query_times.end());
      outfile << "plain, " << internal_bytes << ", " <<  leaf_bytes << ", " <<  max_size << ", " <<  NUM_QUERIES << ", " <<  num_query_sizes[i] << ", " <<  curr_unsorted_query_times[trials/2] << ", " <<  curr_sorted_query_times[trials/2] << ", \n";
    }
    outfile.close();
  }
  return true;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("call with the number of elements to insert\n");
    return -1;
  }
  int trials = 5;
  // if (argc > 2) {
  //   trials = atoi(argv[2]);
  // }
  std::seed_seq seed{0};
  int n = atoi(argv[1]);
  int num_queries = atoi(argv[2]);
  bool write_csv = true;
  // { test_concurrent_sum_time<uint64_t>(n, seed, trials); }
  // { test_concurrent_find<uint64_t>(n, seed); }
  // { test_concurrent_sum<uint64_t>(n, seed); }

  // std::vector<uint64_t> serial_times;
  // std::vector<uint64_t> serial_remove_times;
  // std::vector<uint64_t> parallel_times;
  // std::vector<uint64_t> parallel_remove_times;
  // for (int i = 0; i < trials; i++) {
    // auto [correct, serial, parallel, serial_remove, parallel_remove] =
    //     test_concurrent_btreemap<unsigned long, 256, 1024>(n, num_queries, seed);
    std::ofstream outfile;
    outfile.open("insert_finds.csv", std::ios_base::app); 
    outfile << "tree_type, internal bytes, leaf bytes, num_inserted, insert_time, num_finds, find_time, \n";
    outfile.close();
    outfile.open("range_queries.csv", std::ios_base::app); 
    outfile << "tree_type, internal bytes, leaf bytes, num_inserted,num_range_queries, max_query_size,  unsorted_query_time, sorted_query_time, \n";
    outfile.close();

    // bool correct = test_concurrent_microbenchmarks_map<unsigned long, 256, 256>(n, num_queries, seed, write_csv, trials);
    // correct = test_concurrent_microbenchmarks_map<unsigned long, 512, 512>(n, num_queries, seed, write_csv, trials);
    // correct = test_concurrent_microbenchmarks_map<unsigned long, 1024, 1024>(n, num_queries, seed, write_csv, trials);
    // correct = test_concurrent_microbenchmarks_map<unsigned long, 2048, 2048>(n, num_queries, seed, write_csv, trials);
    // correct = test_concurrent_microbenchmarks_map<unsigned long, 4096, 4096>(n, num_queries, seed, write_csv, trials);
    // correct = test_concurrent_microbenchmarks_map<unsigned long, 8192, 8192>(n, num_queries, seed, write_csv, trials);
    // correct = test_concurrent_microbenchmarks_map<unsigned long, 16384, 16384>(n, num_queries, seed, write_csv, trials);
    // correct = test_concurrent_microbenchmarks_map<unsigned long, 32768, 32768>(n, num_queries, seed, write_csv, trials);
    // correct = test_concurrent_microbenchmarks_map<unsigned long, 65536, 65536>(n, num_queries, seed, write_csv, trials);

    bool correct = test_concurrent_microbenchmarks_map<unsigned long, 1024, 1024>(n, num_queries, seed, write_csv, trials);
    correct = test_concurrent_microbenchmarks_map<unsigned long, 1024, 2048>(n, num_queries, seed, write_csv, trials);
    correct = test_concurrent_microbenchmarks_map<unsigned long, 1024, 4096>(n, num_queries, seed, write_csv, trials);
    correct = test_concurrent_microbenchmarks_map<unsigned long, 1024, 8192>(n, num_queries, seed, write_csv, trials);
    correct = test_concurrent_microbenchmarks_map<unsigned long, 1024, 16384>(n, num_queries, seed, write_csv, trials);
    correct = test_concurrent_microbenchmarks_map<unsigned long, 1024, 32768>(n, num_queries, seed, write_csv, trials);
    correct = test_concurrent_microbenchmarks_map<unsigned long, 1024, 65536>(n, num_queries, seed, write_csv, trials);


    if (!correct) {
      printf("got the wrong answer :(\n");
      return -1;
    }

    // if (i > 0) {
	  //   serial_times.push_back(serial);
	  //   serial_remove_times.push_back(serial_remove);
	  //   parallel_times.push_back(parallel);
	  //   parallel_remove_times.push_back(parallel_remove);
    // }
  // }
  // std::sort(serial_times.begin(), serial_times.end());
  // std::sort(serial_remove_times.begin(), serial_remove_times.end());
  // std::sort(parallel_times.begin(), parallel_times.end());
  // std::sort(parallel_remove_times.begin(), parallel_remove_times.end());
  // printf("%lu, %lu, %lu, %lu\n", serial_times[trials / 2],
  //        parallel_times[trials / 2], serial_remove_times[trials / 2],
  //        parallel_remove_times[trials / 2]);
  return 0;
}
