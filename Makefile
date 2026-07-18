CFLAGS = -march=native -O3 -ffast-math -g
WFLAGS = -Wconversion -Wsign-conversion -Wall -Wextra -Wpedantic -Wshadow -Wpedantic
CC = clang++ -std=c++20

all:
	$(CC) $(CFLAGS) $(WFLAGS) ensim5.cc -c
	$(CC) $(CFLAGS) $(WFLAGS) main.cc ensim5.o
	objdump -dr -C ensim5.o > out.asm
	perf stat -d -d -d ./a.out 44800

clean:
	rm out.asm
	rm a.out
	rm ensim5.o
