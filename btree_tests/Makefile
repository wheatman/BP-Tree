OPT?=fast
VALGRIND?=0
INFO?=0
ENABLE_TRACE_TIMER?=0
CYCLE_TIMER?=0
DEBUG?=0
CILK?=0
SANITIZE?=0
CORRECTNESS?=0
LEAFDS?=0

ifeq ($(LEAFDS),1)
CFLAGS := -Wall -Wextra -O$(OPT) -g  -std=c++20 -IParallelTools/ -IleafDS/ -I../tlx-leafds -gdwarf-4
else
CFLAGS := -Wall -Wextra -O$(OPT) -g  -std=c++20 -IParallelTools/ -I../tlx-plain -gdwarf-4
endif

ifeq ($(SANITIZE),1)
ifeq ($(CILK),1)
CFLAGS += -fsanitize=cilk,undefined,address -fno-omit-frame-pointer
# CFLAGS += -fsanitize=undefined,address -fno-omit-frame-pointer
else
CFLAGS += -fsanitize=undefined,address -fno-omit-frame-pointer
endif
endif

LDFLAGS := -lrt -lm -lpthread -lm -ldl -latomic -ljemalloc
# -ljemalloc -ltcmalloc

ifeq ($(VALGRIND),0)
CFLAGS += -march=native
endif

ifeq ($(LEAFDS),1)
DEFINES := -DENABLE_TRACE_TIMER=$(ENABLE_TRACE_TIMER) -DLEAFDS=$(LEAFDS) -DCYCLE_TIMER=$(CYCLE_TIMER) -DCILK=$(CILK) -DDEBUG=$(DEBUG) -DVALGRIND=$(VALGRIND) -DAVX512=$(AVX512) -DCORRECTNESS=$(CORRECTNESS)
else
DEFINES := -DENABLE_TRACE_TIMER=$(ENABLE_TRACE_TIMER) -DLEAFDS=$(LEAFDS) -DCYCLE_TIMER=$(CYCLE_TIMER) -DCILK=$(CILK) -DDEBUG=$(DEBUG) -DCORRECTNESS=$(CORRECTNESS) -DENABLE_TRACE_TIMER=$(ENABLE_TRACE_TIMER)
endif

ifeq ($(CILK),1)
CFLAGS += -fopencilk -DPARLAY_CILK
ONE_WORKER = CILK_NWORKERS=1
endif


ifeq ($(DEBUG),0)
CFLAGS += -DNDEBUG
endif


ifeq ($(INFO), 1) 
# CFLAGS +=  -Rpass-missed="(inline|loop*)" 
#CFLAGS += -Rpass="(inline|loop*)" -Rpass-missed="(inline|loop*)" -Rpass-analysis="(inline|loop*)" 
CFLAGS += -Rpass=.* -Rpass-missed=.* -Rpass-analysis=.* 
endif


all: basic

ifeq ($(LEAFDS),1)
basic: test-leafds.cpp
	$(CXX) $(CFLAGS) $(DEFINES) $(LDFLAGS) -o $@ test-leafds.cpp
else
$(info    VAR is $(LEAFDS))
basic: test-plain.cpp
	$(CXX) $(CFLAGS) $(DEFINES) $(LDFLAGS) -o $@ test-plain.cpp
endif

ycsb: ycsb.cpp
	$(CXX) $(CFLAGS) $(DEFINES) $(LDFLAGS) -o $@ ycsb.cpp

bwtree: bwtree_test.cpp
	$(CXX) $(CFLAGS) $(DEFINES) $(LDFLAGS) -o $@ bwtree_test.cpp

clean:
	rm -f *.o opt profile basic code.profdata default.profraw
