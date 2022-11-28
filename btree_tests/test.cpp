#include <algorithm>
#include <functional>
#include <random>
#include <sys/time.h>
#include <vector>
#include <set>
#include <map>
#include <ParallelTools/parallel.h>
#include <fstream>

#include <tlx/container/btree_set.hpp>
#include <tlx/container/btree_map.hpp>

#define CORRECTNESS 0
#define TRIALS 5
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

template <class T>
void test_concurrent_sum_correctness(uint64_t max_size, std::seed_seq &seed) {
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


template <class T>
void test_concurrent_sum_time(uint64_t max_size, std::seed_seq &seed, int trials) {
  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);

  std::vector<uint64_t> insert_times(trials);
  std::vector<uint64_t> with_map_times(trials);
  std::vector<uint64_t> with_subtract_times(trials);
  uint64_t start, end, map_time, subtract_time, insert_time;
  for(int i = 0; i < trials + 1; i++) {
	  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
			 std::allocator<T>, false>
	      concurrent_set;
	  start = get_usecs();
	  for(uint32_t i = 0; i < max_size; i++) {
	    concurrent_set.insert(data[i]);
	  }
	  end = get_usecs();
	  insert_time = end - start;


	  start = get_usecs();
	  auto concurrent_sum = concurrent_set.psum_with_map();
	  end = get_usecs();
	  map_time = end - start;
	  
	  printf("\tconcurrent btree sum with map got %lu\n", concurrent_sum);

	  start = get_usecs();
	  concurrent_sum = concurrent_set.psum_with_subtract();
	  end = get_usecs();
	  subtract_time = end - start;
	  printf("\tconcurrent btree sum with subtract got %lu\n", concurrent_sum);
		  if(i > 0) {
			  insert_times[i-1] = insert_time;
			  with_map_times[i-1] = map_time;
			  with_subtract_times[i-1] = subtract_time;
		  }
  }
  std::sort(insert_times.begin(), insert_times.end());
  std::sort(with_map_times.begin(), with_map_times.end());
  std::sort(with_subtract_times.begin(), with_subtract_times.end());
  printf("insert time = %lu\n", insert_times[trials / 2]);
  printf("with map time = %lu, with subtract time = %lu\n", with_map_times[trials / 2], with_subtract_times[trials / 2]);

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
test_concurrent_microbenchmarks_map(uint64_t max_size, uint64_t NUM_QUERIES, std::seed_seq &seed, bool write_csv, int trials) {
  std::vector<uint32_t> num_query_sizes{};
  // std::vector<uint32_t> num_query_sizes{100, 1000, 10000, 100000};

  uint64_t start_time, end_time;
  std::vector<uint64_t> insert_times;
  std::vector<uint64_t> find_times;
  std::vector<uint64_t> sorted_range_query_times_by_size;
  std::vector<uint64_t> unsorted_range_query_times_by_size;

  for (int cur_trial = 0; cur_trial <= trials; cur_trial++) {
    printf("\nRunning leafds btree with internal bytes = %u with leafds slots %lu , trial = %lu\n",internal_bytes, SLOTS, cur_trial);

    // std::vector<uint32_t> num_query_sizes{100, 1000, 10000, 100000};

    // std::vector<T> data =
    //     create_random_data<T>(max_size, max_size * 2, seed);
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

    // uint64_t start_time, end_time, insert_time, find_time;

    // std::vector<uint64_t> sorted_range_query_times_by_size;
    // std::vector<uint64_t> unsorted_range_query_times_by_size;

    // output to tree_type, internal bytes, leaf bytes, num_inserted, insert_time, num_finds, find_time, num_range_queries, range_time_maxlen_{}*

    tlx::btree_map<T, T, std::less<T>, tlx::btree_default_traits<T, T, internal_bytes>,
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
    std::vector<bool> found_count(NUM_QUERIES);
    start_time = get_usecs();
    cilk_for(uint32_t i = 0; i < NUM_QUERIES; i++) {
      found_count[i] = concurrent_map.exists(data[range_query_start_idxs[i]]);
    }
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
      std::vector<T> concurrent_range_query_length_key_sums(NUM_QUERIES);
      std::vector<T> concurrent_range_query_length_val_sums(NUM_QUERIES);

      cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
        T start;
  #if CORRECTNESS
        start = checker_sorted[range_query_start_idxs[i]];
  #else 
        start = data[range_query_start_idxs[i]];
  #endif

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
      // printf("\n");

      start_time = get_usecs();
      cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
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
      }
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

#if CORRECTNESS
      std::vector<T> concurrent_range_query_maxs(NUM_QUERIES);
      std::vector<uint64_t> concurrent_range_query_counts(NUM_QUERIES);
#endif
      std::vector<T> concurrent_range_query_key_sums(NUM_QUERIES);
      std::vector<T> concurrent_range_query_val_sums(NUM_QUERIES);

#if CORRECTNESS
      cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
        T start, end;
        start = checker_sorted[range_query_start_idxs[i]];
        end = concurrent_range_query_length_maxs[i];
        if (range_query_lengths[i] != 0) {
          end++;
        }
        uint64_t num_in_range = 0;
        T max_key_in_range = start;
        concurrent_map.map_range(start, end, [&num_in_range, &max_key_in_range]([[maybe_unused]] auto key, auto val) {
                  num_in_range += 1;
                  if (key > max_key_in_range) {
                    max_key_in_range = key;
                  }
                });
        concurrent_range_query_counts[i] = num_in_range;
        concurrent_range_query_maxs[i] = max_key_in_range;
      }
      // printf("\n");
#endif

      start_time = get_usecs();
      cilk_for (uint32_t i = 0; i < NUM_QUERIES; i++) {
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
      }
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
    outfile << "leafds, " << internal_bytes << ", " << SLOTS << ", " <<  max_size << ", " <<  insert_time_med << ", " <<  NUM_QUERIES << ", " <<  find_time_med << ", \n";
    outfile.close();
    outfile.open("range_queries.csv", std::ios_base::app); 
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
      printf("wrong iterator value, expected %lu but got %lu on count = %lu, iter = %lu\n", correct_key, leafds_key, count, it_leafds);
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

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("call with the number of elements to insert\n");
    return -1;
  }
  int trials = 5;

  std::seed_seq seed{0};
  int n = atoi(argv[1]);
  int num_queries = atoi(argv[2]);
  bool write_csv = true;

  std::ofstream outfile;
  outfile.open("insert_finds.csv", std::ios_base::app); 
  outfile << "tree_type, internal bytes, leaf slots, num_inserted, insert_time, num_finds, find_time, \n";
  outfile.close();
  outfile.open("range_queries.csv", std::ios_base::app); 
  outfile << "tree_type, internal bytes, leaf slots, num_inserted,num_range_queries, max_query_size,  unsorted_query_time, sorted_query_time, \n";
  outfile.close();

  // bool correct = test_iterator_merge_map<unsigned long, 1024>(n, seed, write_csv, trials);
  bool correct = test_concurrent_microbenchmarks_map<unsigned long, 1024>(n, num_queries, seed, write_csv, trials);

  if (!correct) {
    printf("got the wrong answer :(\n");
    return -1;
  }

  return 0;
}
