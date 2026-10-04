CC = clang++ -std=c++20 -O3 -ffast-math -march=native -g -Wall -Wextra -Wpedantic -Wfatal-errors

SAN = 0

ifeq ($(SAN),1)
CC += -fsanitize=undefined,thread
endif

ifeq ($(SAN),2)
CC += -fsanitize=undefined,address
endif

all: raylib

ensim.o: ensim.cc ensim.hh Makefile
	$(CC) -c ensim.cc -Wdouble-promotion
	objdump -dr -C ensim.o > ensim.asm

raylib.o: raylib.cc ensim.hh Makefile
	$(CC) -c raylib.cc

raylib: raylib.o ensim.o Makefile
	$(CC) -lraylib raylib.o ensim.o -o raylib

clean:
	rm -f ensim.asm raylib ensim.o raylib.o
