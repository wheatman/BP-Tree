#include <iostream>
#include <chrono>
#include <random>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sys/time.h>
// #include "tbb/tbb.h"
#include <container/btree_set.hpp>
#include <container/btree_map.hpp>
#include <../../btree_tests/BwTree/test/test_suite.h>
#include<cilk/cilk.h>
#include<thread>

static long get_usecs() {
  struct timeval st;
  gettimeofday(&st, NULL);
  return st.tv_sec * 1000000 + st.tv_usec;
}

using namespace std;

using Key=uint64_t;

using TID=uint64_t;

// index types
enum {
    TYPE_BTREE,
    TYPE_ART,
    TYPE_HOT,
    TYPE_BWTREE,
    TYPE_MASSTREE,
    TYPE_CLHT,
    TYPE_FASTFAIR,
    TYPE_LEVELHASH,
    TYPE_CCEH,
    TYPE_WOART,
};

enum {
    OP_INSERT,
    OP_UPDATE,
    OP_READ,
    OP_SCAN,
    OP_SCAN_END,
    OP_DELETE,
};

enum {
    WORKLOAD_A,
    WORKLOAD_B,
    WORKLOAD_C,
    WORKLOAD_D,
    WORKLOAD_E,
    WORKLOAD_X,
    WORKLOAD_Y,
};

enum {
    RANDINT_KEY,
    STRING_KEY,
};

enum {
    UNIFORM,
    ZIPFIAN,
};

namespace Dummy {
    inline void mfence() {asm volatile("mfence":::"memory");}

    inline void clflush(char *data, int len, bool front, bool back)
    {
        if (front)
            mfence();
        volatile char *ptr = (char *)((unsigned long)data & ~(64 - 1));
        for (; ptr < data+len; ptr += 64){
#ifdef CLFLUSH
            asm volatile("clflush %0" : "+m" (*(volatile char *)ptr));
#elif CLFLUSH_OPT
            asm volatile(".byte 0x66; clflush %0" : "+m" (*(volatile char *)(ptr)));
#elif CLWB
            asm volatile(".byte 0x66; xsaveopt %0" : "+m" (*(volatile char *)(ptr)));
#endif
        }
        if (back)
            mfence();
    }
}


static uint64_t LOAD_SIZE = 100000000;
static uint64_t RUN_SIZE = 100000000;

void loadKey(TID tid, Key &key) {
    return ;
}

void ycsb_load_run_string(int index_type, int wl, int kt, int ap, int num_thread,
        std::vector<Key *> &init_keys,
        std::vector<Key *> &keys,
        std::vector<int> &ranges,
        std::vector<int> &ops)
{
    
}

template <typename F> inline void parallel_for(size_t start, size_t end, F f) {
  cilk_for(size_t i = start; i < end; i++) f(i);
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

void ycsb_load_run_randint(int index_type, int wl, int kt, int ap, int num_thread,
        std::vector<uint64_t> &init_keys,
        std::vector<uint64_t> &keys,
        std::vector<uint64_t> &range_end,
        std::vector<int> &ranges,
        std::vector<int> &ops)
{
    std::string init_file;
    std::string txn_file;

    if (ap == UNIFORM) {
        if (kt == RANDINT_KEY && wl == WORKLOAD_A) {
            init_file = "../../ycsb/index-microbench/workloads/uniform/loada_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/uniform/txnsa_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_B) {
            init_file = "../../ycsb/index-microbench/workloads/uniform/loadb_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/uniform/txnsb_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_C) {
            init_file = "../../ycsb/index-microbench/workloads/uniform/loadc_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/uniform/txnsc_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_D) {
            init_file = "../../ycsb/index-microbench/workloads/uniform/loadd_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/uniform/txnsd_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_E) {
            init_file = "../../ycsb/index-microbench/workloads/uniform/loade_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/uniform/txnse_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_X) {
            init_file = "../../ycsb/index-microbench/workloads/uniform/loadx_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/uniform/txnsx_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_Y) {
            init_file = "../../ycsb/index-microbench/workloads/uniform/loady_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/uniform/txnsy_unif_int.dat";
        }
    } else {
        if (kt == RANDINT_KEY && wl == WORKLOAD_A) {
            init_file = "../../ycsb/index-microbench/workloads/zipfian/loada_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/zipfian/txnsa_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_B) {
            init_file = "../../ycsb/index-microbench/workloads/zipfian/loadb_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/zipfian/txnsb_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_C) {
            init_file = "../../ycsb/index-microbench/workloads/zipfian/loadc_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/zipfian/txnsc_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_D) {
            init_file = "../../ycsb/index-microbench/workloads/zipfian/loadd_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/zipfian/txnsd_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_E) {
            init_file = "../../ycsb/index-microbench/workloads/zipfian/loade_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/zipfian/txnse_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_X) {
            init_file = "../../ycsb/index-microbench/workloads/zipfian/loadx_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/zipfian/txnsx_unif_int.dat";
        } else if (kt == RANDINT_KEY && wl == WORKLOAD_Y) {
            init_file = "../../ycsb/index-microbench/workloads/zipfian/loady_unif_int.dat";
            txn_file = "../../ycsb/index-microbench/workloads/zipfian/txnsy_unif_int.dat";
        }
    }

    std::ifstream infile_load(init_file);

    std::string op;
    uint64_t key;
    uint64_t rend;
    int range;

    std::string insert("INSERT");
    std::string update("UPDATE");
    std::string read("READ");
    std::string scan("SCAN");
    std::string scanend("SCANEND");

    int count = 0;
    while ((count < LOAD_SIZE) && infile_load.good()) {
        infile_load >> op >> key;
        if (op.compare(insert) != 0) {
            std::cout << "READING LOAD FILE FAIL!\n";
            return ;
        }
        init_keys.push_back(key);
        count++;
    }

    fprintf(stderr, "Loaded %d keys\n", count);

    std::ifstream infile_txn(txn_file);

    count = 0;
    while ((count < RUN_SIZE) && infile_txn.good()) {
        infile_txn >> op >> key;
        if (op.compare(insert) == 0) {
            ops.push_back(OP_INSERT);
            keys.push_back(key);
            ranges.push_back(1);
            range_end.push_back(1);
        } else if (op.compare(update) == 0) {
            ops.push_back(OP_UPDATE);
            keys.push_back(key);
            ranges.push_back(1);
            range_end.push_back(1);
        } else if (op.compare(read) == 0) {
            ops.push_back(OP_READ);
            keys.push_back(key);
            ranges.push_back(1);
            range_end.push_back(1);
        } else if (op.compare(scan) == 0) {
            infile_txn >> range;
            ops.push_back(OP_SCAN);
            keys.push_back(key);
            ranges.push_back(range);
            range_end.push_back(1);
        } else if (op.compare(scanend) == 0) {
            infile_txn >> rend;
            ops.push_back(OP_SCAN_END);
            keys.push_back(key);
            range_end.push_back(rend);
            ranges.push_back(1);
        } else {
            std::cout << "UNRECOGNIZED CMD!\n";
            return;
        }
        count++;
    }

    std::atomic<int> range_complete, range_incomplete;
    range_complete.store(0);
    range_incomplete.store(0);

    fprintf(stderr, "Loaded %d more keys\n", count);

    std::this_thread::sleep_for(std::chrono::nanoseconds(3000000000));

    fprintf(stderr, "Slept\n");

    if (index_type == TYPE_BTREE) {


        for(int k =0; k<6; k++){
            std::vector<uint64_t> query_results_keys(RUN_SIZE);
            std::vector<uint64_t> query_results_vals(RUN_SIZE);
            tlx::btree_map<uint64_t, uint64_t, std::less<uint64_t>, tlx::btree_default_traits<uint64_t, uint64_t>,
                            std::allocator<uint64_t>, true> concurrent_map;
            {
                // Load
                auto starttime = get_usecs(); // std::chrono::system_clock::now();
                parallel_for(0, LOAD_SIZE, [&](const uint64_t &i) {
                    concurrent_map.insert({init_keys[i], init_keys[i]});
                });
                auto end = get_usecs();
                auto duration = end- starttime; //std::chrono::duration_cast<std::chrono::microseconds>(
                        //std::chrono::system_clock::now() - starttime);
                printf("\tLoad took %lu us, throughput = %f ops/us\n", duration, ((double)LOAD_SIZE)/duration);
                //printf("Throughput: load, %f ,ops/us and time %ld in us\n", (LOAD_SIZE * 1.0) / duration.count(), duration.count());
            }
        {
            // Run
            auto starttime = std::chrono::system_clock::now();
            parallel_for(0, RUN_SIZE, [&](const uint64_t &i) {
                    if (ops[i] == OP_INSERT) {
                        concurrent_map.insert({keys[i], keys[i]});
                    } else if (ops[i] == OP_READ) {
                        if(concurrent_map.exists(keys[i]) &&  concurrent_map.value(keys[i]) != keys[i]) {
                            std::cout << "[BTREE] wrong key read: " << concurrent_map.value(keys[i]) << " expected:" << keys[i] << std::endl;
                            exit(0);
                        }
                    } else if (ops[i] == OP_SCAN) {
                        uint64_t start= keys[i];
			            uint64_t key_sum = 0, val_sum = 0;
                        concurrent_map.map_range_length(keys[i], ranges[i], [&key_sum, &val_sum]([[maybe_unused]] auto el) {
                            key_sum += el.first;
                            val_sum += el.second;
                        });

                        // concurrent_map.map_range_length(keys[i], ranges[i], [&key_sum, &val_sum](auto key, auto value) {
                        //     key_sum += key;
                        //     val_sum += value;
                        // });

                        query_results_keys[i] = key_sum;
                        query_results_vals[i] = val_sum;
                    } else if (ops[i] == OP_SCAN_END) {
			            uint64_t key_sum = 0, val_sum = 0;
                        concurrent_map.map_range(keys[i], range_end[i], [&key_sum, &val_sum]([[maybe_unused]] auto el) {
                            key_sum += el.first;
                            val_sum += el.second;
                        });
                        // concurrent_map.map_range(keys[i], range_end[i], [&key_sum, &val_sum](auto key, auto value) {
                        //     key_sum += key;
                        //     val_sum += value;
                        // });
                        query_results_keys[i] = key_sum;
                        query_results_vals[i] = val_sum;
                    } else if (ops[i] == OP_UPDATE) {
                        std::cout << "NOT SUPPORTED CMD!\n";
                        exit(0);
                    }
            });
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::system_clock::now() - starttime);
            printf("\tRun, throughput: %f ,ops/us\n", (RUN_SIZE * 1.0) / duration.count());
        }
        uint64_t key_sum = 0;
        uint64_t val_sum = 0;
        for(int i = 0; i < RUN_SIZE; i++) {
            key_sum += query_results_keys[i];
            val_sum += query_results_vals[i];
        }
        printf("\ttotal key sum = %lu, total val sum = %lu\n\n", key_sum, val_sum);
        }
    }
}


int main(int argc, char **argv) {
    if (argc != 6) {
        std::cout << "Usage: ./ycsb [index type] [ycsb workload type] [key distribution] [access pattern] [number of threads]\n";
        std::cout << "1. index type: art hot bwtree masstree clht\n";
        std::cout << "               fastfair levelhash cceh woart\n";
        std::cout << "2. ycsb workload type: a, b, c, e\n";
        std::cout << "3. key distribution: randint, string\n";
        std::cout << "4. access pattern: uniform, zipfian\n";
        std::cout << "5. number of threads (integer)\n";
        return 1;
    }

    printf("%s, workload%s, %s, %s, threads %s\n", argv[1], argv[2], argv[3], argv[4], argv[5]);

    int index_type;
    if (strcmp(argv[1], "art") == 0)
        index_type = TYPE_ART;

    else if (strcmp(argv[1], "btree") == 0)
        index_type = TYPE_BTREE;
    else if (strcmp(argv[1], "hot") == 0) {
#ifdef HOT
        index_type = TYPE_HOT;
#else
        return 1;
#endif
    } else if (strcmp(argv[1], "bwtree") == 0)
        index_type = TYPE_BWTREE;
    else if (strcmp(argv[1], "masstree") == 0)
        index_type = TYPE_MASSTREE;
    else if (strcmp(argv[1], "clht") == 0)
        index_type = TYPE_CLHT;
    else if (strcmp(argv[1], "fastfair") == 0)
        index_type = TYPE_FASTFAIR;
    else if (strcmp(argv[1], "levelhash") == 0)
        index_type = TYPE_LEVELHASH;
    else if (strcmp(argv[1], "cceh") == 0)
        index_type = TYPE_CCEH;
    else if (strcmp(argv[1], "woart") == 0)
        index_type = TYPE_WOART;
    else {
        fprintf(stderr, "Unknown index type: %s\n", argv[1]);
        exit(1);
    }

    int wl;
    if (strcmp(argv[2], "a") == 0) {
        wl = WORKLOAD_A;
    } else if (strcmp(argv[2], "b") == 0) {
        wl = WORKLOAD_B;
    } else if (strcmp(argv[2], "c") == 0) {
        wl = WORKLOAD_C;
    } else if (strcmp(argv[2], "d") == 0) {
        wl = WORKLOAD_D;
    } else if (strcmp(argv[2], "e") == 0) {
        wl = WORKLOAD_E;
    } else if (strcmp(argv[2], "x") == 0) {
        wl = WORKLOAD_X;
    } else if (strcmp(argv[2], "y") == 0) {
        wl = WORKLOAD_Y;
    } else {
        fprintf(stderr, "Unknown workload: %s\n", argv[2]);
        exit(1);
    }

    int kt;
    if (strcmp(argv[3], "randint") == 0) {
        kt = RANDINT_KEY;
    } else if (strcmp(argv[3], "string") == 0) {
        kt = STRING_KEY;
    } else {
        fprintf(stderr, "Unknown key type: %s\n", argv[3]);
        exit(1);
    }

    int ap;
    if (strcmp(argv[4], "uniform") == 0) {
        ap = UNIFORM;
    } else if (strcmp(argv[4], "zipfian") == 0) {
        ap = ZIPFIAN;
    } else {
        fprintf(stderr, "Unknown access pattern: %s\n", argv[4]);
        exit(1);
    }

    int num_thread = atoi(argv[5]);
    // tbb::task_scheduler_init init(num_thread);

    if (kt != STRING_KEY) {
        std::vector<uint64_t> init_keys;
        std::vector<uint64_t> keys;
        std::vector<uint64_t> ranges_end;
        std::vector<int> ranges;
        std::vector<int> ops;

        init_keys.reserve(LOAD_SIZE);
        keys.reserve(RUN_SIZE);
        ranges_end.reserve(RUN_SIZE);
        ranges.reserve(RUN_SIZE);
        ops.reserve(RUN_SIZE);

        memset(&init_keys[0], 0x00, LOAD_SIZE * sizeof(uint64_t));
        memset(&keys[0], 0x00, RUN_SIZE * sizeof(uint64_t));
        memset(&ranges_end[0], 0x00, RUN_SIZE * sizeof(uint64_t));
        memset(&ranges[0], 0x00, RUN_SIZE * sizeof(int));
        memset(&ops[0], 0x00, RUN_SIZE * sizeof(int));

        ycsb_load_run_randint(index_type, wl, kt, ap, num_thread, init_keys, keys,ranges_end, ranges, ops);
    } else {
        // std::vector<Key *> init_keys;
        // std::vector<Key *> keys;
        // std::vector<int> ranges;
        // std::vector<int> ops;

        // init_keys.reserve(LOAD_SIZE);
        // keys.reserve(RUN_SIZE);
        // ranges.reserve(RUN_SIZE);
        // ops.reserve(RUN_SIZE);

        // memset(&init_keys[0], 0x00, LOAD_SIZE * sizeof(Key *));
        // memset(&keys[0], 0x00, RUN_SIZE * sizeof(Key *));
        // memset(&ranges[0], 0x00, RUN_SIZE * sizeof(int));
        // memset(&ops[0], 0x00, RUN_SIZE * sizeof(int));

        // ycsb_load_run_string(index_type, wl, kt, ap, num_thread, init_keys, keys, ranges, ops);
    }

    return 0;
}