#!/bin/sh
numactl -N 1 -m 1 ./ycsb btree y randint uniform 48
numactl -N 1 -m 1 ./ycsb btree y randint zipfian 48
