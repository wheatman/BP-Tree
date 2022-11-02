#include <algorithm>
#include <functional>
#include <random>
#include <sys/time.h>
#include <vector>
#include <set>
#include <map>
#include <ParallelTools/parallel.h>

#include <tlx/container/btree_set.hpp>
#include <tlx/container/btree_map.hpp>

#define CORRECTNESS 1
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
  bool wrong = false;

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

  // /*
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
  // */

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

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("call with the number of elements to insert\n");
    return -1;
  }
  int trials = 0; // used to be 5
  if (argc > 2) {
    trials = atoi(argv[2]);
  }
  std::seed_seq seed{0};
  int n = atoi(argv[1]);

  // { test_concurrent_sum_time<uint64_t>(n, seed, trials); }

  // return 0; 
  std::vector<uint64_t> serial_times;
  std::vector<uint64_t> parallel_times;
  std::vector<uint64_t> serial_remove_times;
  std::vector<uint64_t> parallel_remove_times;
  for (int i = 0; i < trials + 1; i++) {
    auto [correct, serial_insert, parallel_insert, serial_remove, parallel_remove] =
        // test_concurrent_btreemap<uint64_t>(n, seed);
        test_concurrent_range_query<uint64_t>(n, seed);
    if (!correct) {
      printf("got the wrong answer\n");
      return -1;
    }
    // throw out the first trial
    // if (i > 0) {
    if (true) {
    	serial_times.push_back(serial_insert);
    	parallel_times.push_back(parallel_insert);
      serial_remove_times.push_back(serial_remove);
    	parallel_remove_times.push_back(parallel_remove);
    }
  }
  std::sort(serial_times.begin(), serial_times.end());
  std::sort(parallel_times.begin(), parallel_times.end());
  std::sort(serial_remove_times.begin(), serial_remove_times.end());
  std::sort(parallel_remove_times.begin(), parallel_remove_times.end());
  printf("%lu, %lu, %lu, %lu\n", serial_times[trials / 2], parallel_times[trials / 2],
                                serial_remove_times[trials / 2], parallel_remove_times[trials / 2]);
  return 0;
}
