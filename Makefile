CC = clang++ -std=c++20 -O3 -ffast-math -march=native -g -Wall -Wextra -Wpedantic -Wfatal-errors

SAN = 0

ifeq ($(SAN),1)
CC += -fsanitize=undefined,thread
endif

ifeq ($(SAN),2)
CC += -fsanitize=undefined,address
endif

all: gui

ensim.o: ensim.cc ensim.hh Makefile
	$(CC) -c ensim.cc -Wdouble-promotion
	objdump -dr -C ensim.o > ensim.asm

gui.o: gui.cc ensim.hh Makefile
	$(CC) -c gui.cc

gui: gui.o ensim.o Makefile
	$(CC) -lraylib gui.o ensim.o -o gui

clean:
	rm -f ensim.asm gui ensim.o gui.o
