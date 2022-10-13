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
std::tuple<bool, uint64_t, uint64_t>
test_concurrent_btreeset(uint64_t max_size, std::seed_seq &seed) {
  std::vector<T> data =
      create_random_data<T>(max_size, std::numeric_limits<T>::max(), seed);

  uint64_t start, end;

#if DEBUG
  tlx::btree_set<T> serial_set;
  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    serial_set.insert(data[i]);
  }
  end = get_usecs();
  int64_t serial_time = end - start;
  printf("inserted all the data serially in %lu\n", end - start);
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
#if DEBUG
  if (serial_set.size() != concurrent_set.size()) {
    printf("the sizes don't match, got %lu, expetected %lu\n",
           concurrent_set.size(), serial_set.size());
    return {false, 0, 0};
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
    return {false, 0, 0};
  }
  return {true, serial_time, parallel_time};
#endif
  return {true, serial_time, parallel_time};
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
  std::vector<uint64_t> parallel_times;
  for (int i = 0; i < trials + 1; i++) {
    printf("trial %d for %d elts\n", i, n);
    auto [correct, serial, parallel] =
        test_concurrent_btreeset<unsigned long>(n, seed);
    if (!correct) {
      printf("got the wrong answer\n");
      return -1;
    }
    if (i > 0) {
	    serial_times.push_back(serial);
	    parallel_times.push_back(parallel);
    }
  }
  std::sort(serial_times.begin(), serial_times.end());
  std::sort(parallel_times.begin(), parallel_times.end());
  printf("***serial = %lu, parallel = %lu***\n", serial_times[trials / 2], parallel_times[trials / 2]);
  return 0;
}
