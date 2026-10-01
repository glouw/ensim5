CC = clang++ -std=c++20 -O3 -ffast-math -march=native -g -Wall -Wextra -Wpedantic -Wfatal-errors

SAN = 0

ifeq ($(SAN),1)
CC += -fsanitize=undefined,thread
endif

ifeq ($(SAN),2)
CC += -fsanitize=undefined,address
endif

all: sdl raylib

perf: sdl
	perf stat -d -d -d -r 5 ./sdl --perf

ensim.o: ensim.cc ensim.hh Makefile
	$(CC) -c ensim.cc -Wdouble-promotion
	objdump -dr -C ensim.o > ensim.asm

sdl.o: sdl.cc ensim.hh Makefile
	$(CC) -c sdl.cc

sdl: sdl.o ensim.o Makefile
	$(CC) -lSDL3 sdl.o ensim.o -o sdl

raylib.o: raylib.cc ensim.hh Makefile
	$(CC) -c raylib.cc

raylib: raylib.o ensim.o Makefile
	$(CC) -lraylib raylib.o ensim.o -o raylib

clean:
	rm -f ensim.asm raylib sdl ensim.o sdl.o raylib.o
