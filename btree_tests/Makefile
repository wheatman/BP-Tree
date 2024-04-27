OPT?=fast
VALGRIND?=0
INFO?=0
ENABLE_TRACE_TIMER?=0
CYCLE_TIMER?=0
DEBUG?=0
SANITIZE?=0
CORRECTNESS?=0
LEAFDS?=0
AVX512?=0

ifeq ($(LEAFDS),1)
CFLAGS := -pthread -Wall -Wextra -O$(OPT) -g  -std=c++20 -IParallelTools/ -I../tlx-leafds -fmax-errors=1
else
CFLAGS := -pthread -Wall -Wextra -O$(OPT) -g  -std=c++20 -IParallelTools/ -I../tlx-plain -fmax-errors=1
endif

ifeq ($(SANITIZE),1)
CFLAGS += -fsanitize=undefined,address -fno-omit-frame-pointer
endif

LDFLAGS := -lrt -lm -lm -ldl -latomic -ljemalloc
# -ljemalloc -ltcmalloc

ifeq ($(VALGRIND),0)
CFLAGS += -march=native
endif

ifeq ($(LEAFDS),1)
DEFINES := -DENABLE_TRACE_TIMER=$(ENABLE_TRACE_TIMER) -DLEAFDS=$(LEAFDS) -DCYCLE_TIMER=$(CYCLE_TIMER) -DDEBUG=$(DEBUG) -DVALGRIND=$(VALGRIND) -DAVX512=$(AVX512) -DCORRECTNESS=$(CORRECTNESS)
else
DEFINES := -DENABLE_TRACE_TIMER=$(ENABLE_TRACE_TIMER) -DLEAFDS=$(LEAFDS) -DCYCLE_TIMER=$(CYCLE_TIMER) -DDEBUG=$(DEBUG) -DCORRECTNESS=$(CORRECTNESS) -DENABLE_TRACE_TIMER=$(ENABLE_TRACE_TIMER)
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
	rm -f *.o opt profile basic ycsb code.profdata default.profraw
