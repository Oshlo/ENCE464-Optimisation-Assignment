# Makefile for building poisson.c and threads.c

TARGETS = poisson

CROSS_TOOL =
CC_CPP = $(CROSS_TOOL)g++
CC_C = $(CROSS_TOOL)gcc

CPP_FLAGS = -O3
C_FLAGS = -Wall -Werror -g -std=c99 -O3

all: clean $(TARGETS)

$(TARGETS):
	$(CC_CPP) $(CPP_FLAGS) $@.c -o $@ -lpthread 

clean:
	-rm -f $(TARGETS)

run:
	./$(TARGETS)

test:
	 ./test.sh

debug:
	./$(TARGETS) --debug --test

show_iterations:
	./$(TARGETS) --show-iterations

run100: 
	./$(TARGETS) -n 101 -i 100

run200:
	./$(TARGETS) -n 201 -i 100

extra_credit:
	./$(TARGETS) -n 301 -i 500 -t 8

thread1:
	./$(TARGETS) -n 101 -i 200 -t 1

thread2:
	./$(TARGETS) -n 101 -i 200 -t 2

thread4:
	./$(TARGETS) -n 101 -i 200 -t 4








