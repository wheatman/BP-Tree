numactl -N 1 perf stat -a -e cache-misses,LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,LLC-prefetch-misses,dTLB-load-misses,dTLB-prefetches,iTLB-load-misses ./basic --cache_misses

notes: look at the file
first do the inserts


then rebuild with commenting back in the sorted range
  run with --query_size {100, 1000, 10000, 100000}

then unsorted range
  run with --query_size {100, 1000, 10000, 100000}

