ASANFLAGS =-fsanitize=address
CFLAGS = -O3 -ffast-math -march=native -g
WFLAGS = -Wall -Wextra -Wpedantic
CC = clang++ -std=c++20
LDFLAGS = -lSDL3
OUT = ensim
ASM = ensim.asm
OBJS = ensim.o main.o
DEPS = ensim.hh Makefile
PERF = perf stat -d -d -d -r 20
DUMP = objdump -dr -C

run: all
	./$(OUT)

perf: all
	$(PERF) ./$(OUT) --perf

all: $(OBJS)
	$(CC) $(CFLAGS) $(WFLAGS) $(ASANFLAGS) $(LDFLAGS) $^ -o $(OUT)

main.o: main.cc $(DEPS)
	$(CC) $(CFLAGS) $(WFLAGS) $(ASANFLAGS) $< -c

ensim.o: ensim.cc $(DEPS)
	$(CC) $(CFLAGS) $(WFLAGS) $(ASANFLAGS) $< -c
	$(DUMP) $@ > $(ASM)

clean:
	rm -f $(ASM) $(OUT) $(OBJS)
