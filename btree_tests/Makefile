OPT?=fast
VALGRIND?=0
INFO?=0
ENABLE_TRACE_TIMER?=0
CYCLE_TIMER?=0
DEBUG?=0
CILK?=0
SANITIZE?=0



CFLAGS := -Wall -Wextra -O$(OPT) -g  -std=c++20 -IParallelTools/ -IleafDS/ -I../ -gdwarf-4 -ferror-limit=1

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

DEFINES := -DENABLE_TRACE_TIMER=$(ENABLE_TRACE_TIMER) -DCYCLE_TIMER=$(CYCLE_TIMER) -DCILK=$(CILK) -DDEBUG=$(DEBUG) -DVALGRIND=$(VALGRIND) -DAVX512=$(AVX512)

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
 
basic: test.cpp
	$(CXX) $(CFLAGS) $(DEFINES) $(LDFLAGS) -o $@ test.cpp

ycsb: ycsb.cpp
	$(CXX) $(CFLAGS) $(DEFINES) $(LDFLAGS) -o $@ ycsb.cpp

ycsb_test_copy: ycsb_test_copy.cpp
	$(CXX) $(CFLAGS) $(DEFINES) $(LDFLAGS) -o $@ ycsb_test_copy.cpp

clean:
	rm -f *.o opt profile basic code.profdata default.profraw
