#!/bin/bash

# 1st arg is distribution, second arg is num threads

./ycsb btree a randint $1 $2
./ycsb btree b randint $1 $2
./ycsb btree c randint $1 $2
./ycsb btree e randint $1 $2
./ycsb btree x randint $1 $2
./ycsb btree y randint $1 $2
