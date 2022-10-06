#include <algorithm>
#include <functional>
#include <random>
#include <sys/time.h>
#include <vector>

#include <ParallelTools/parallel.h>

#include <tlx/container/btree_set.hpp>

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

  uint64_t start, end;

  tlx::btree_set<T> serial_set;

  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    serial_set.insert(data[i]);
  }
  end = get_usecs();
  int64_t serial_time = end - start;
  printf("inserted all the data serially in %lu\n", end - start);

  tlx::btree_set<T, std::less<T>, tlx::btree_default_traits<T, T>,
                 std::allocator<T>, true>
      concurrent_set;
  start = get_usecs();
  cilk_for(uint32_t i = 0; i < max_size; i++) {
    concurrent_set.insert(data[i]);
  }
  end = get_usecs();
  uint64_t parallel_time = end - start;
  printf("inserted all the data concurrently in %lu\n", end - start);

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

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("call with the number of elements to insert\n");
    return -1;
  }
  int trials = 1;
  if (argc > 2) {
    trials = atoi(argv[2]);
  }
  std::seed_seq seed{0};
  int n = atoi(argv[1]);
  std::vector<uint64_t> serial_times;
  std::vector<uint64_t> serial_remove_times;
  std::vector<uint64_t> parallel_times;
  std::vector<uint64_t> parallel_remove_times;
  for (int i = 0; i < trials; i++) {
    auto [correct, serial, parallel, serial_remove, parallel_remove] =
        test_concurrent_btreeset<unsigned long>(n, seed);
    if (!correct) {
      printf("got the wrong answer\n");
      return -1;
    }
    serial_times.push_back(serial);
    serial_remove_times.push_back(serial_remove);
    parallel_times.push_back(parallel);
    parallel_remove_times.push_back(parallel_remove);
  }
  std::sort(serial_times.begin(), serial_times.end());
  std::sort(serial_remove_times.begin(), serial_remove_times.end());
  std::sort(parallel_times.begin(), parallel_times.end());
  std::sort(parallel_remove_times.begin(), parallel_remove_times.end());
  printf("%lu, %lu, %lu, %lu\n", serial_times[trials / 2],
         parallel_times[trials / 2], serial_remove_times[trials / 2],
         parallel_remove_times[trials / 2]);
  return 0;
}