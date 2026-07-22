CFLAGS = -O3 -ffast-math -march=native -g -fsanitize=address,undefined
WFLAGS = -Wall -Wextra -Wpedantic -Wshadow -Wpedantic
CC = clang++ -std=c++20
LDFLAGS = -lSDL3
ENSIM_PERF = #-DENSIM_PERF

all:
	$(CC) $(CFLAGS) $(WFLAGS) ensim.cc -c
	$(CC) $(CFLAGS) $(WFLAGS) $(ENSIM_PERF) $(LDFLAGS) main.cc ensim.o
	objdump -dr -C ensim.o > out.asm
	perf stat -d -d -d ./a.out 44800

clean:
	rm out.asm
	rm a.out
	rm ensim.o
