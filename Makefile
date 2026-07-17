CFLAGS = -march=native -O3 -ffast-math -g
WFLAGS = -Wconversion -Wsign-conversion -Wall -Wextra -Wpedantic -Wshadow -Wpedantic
CC = clang++ -std=c++20

all:
	$(CC) $(CFLAGS) $(WFLAGS) -c main.cc
	objdump -dr -C main.o > out.asm
