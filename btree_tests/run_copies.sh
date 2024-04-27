: <<'END'
for i in {0..23}
  do
    taskset -c $i ./basic --microbenchmark_baseline --num_inserts=100000000 --num_queries=0 --trials=5 &
  done
for i in {48..71}
  do
    taskset -c $i ./basic --microbenchmark_baseline --num_inserts=100000000 --num_queries=0 --trials=5 &
  done
  wait
END

taskset -c 0-23 ./basic --microbenchmark_baseline --num_inserts=100000000 --num_queries=0 --trials=5 &
taskset -c 48-71 ./basic --microbenchmark_baseline --num_inserts=100000000 --num_queries=0 --trials=5 &
wait
