CC = clang++ -std=c++20 -O3 -ffast-math -march=native -g -Wall -Wextra -Wpedantic
ifeq (0,1)
CC += -fsanitize=address,undefined,thread
endif
LDFLAGS = -lSDL3

run: demo
	./demo

perf: demo
	perf stat -d -d -d -r 5 ./demo --perf

ensim.o: ensim/ensim.cc ensim/*.hh Makefile
	$(CC) -c ensim/ensim.cc
	objdump -dr -C ensim.o > ensim.asm

demo.o: demo.cc Makefile
	$(CC) -c demo.cc

demo: demo.o ensim.o Makefile
	$(CC) $(LDFLAGS) demo.o ensim.o -o demo

clean:
	rm -f ensim.asm demo ensim.o
