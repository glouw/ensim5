ASANFLAGS = #-fsanitize=address,undefined
CFLAGS = -O2 -ffast-math -march=native -g
WFLAGS = -Wall -Wextra -Wpedantic -Wshadow -Wpedantic
CC = clang++ -std=c++20
LDFLAGS = -lSDL3
OUT = ensim
ASM = ensim.asm
OBJS = ensim.o main.o
DEPS = ensim.hh Makefile
PERF = perf stat -d -d -d
DUMP = objdump -dr -C

run: all
	./$(OUT)

perf: all
	$(PERF) ./$(OUT) --perf

all: $(OBJS)
	$(CC) $(CFLAGS) $(WFLAGS) $(ASANFLAGS) $(LDFLAGS) $^ -o $(OUT)

main.o: main.cc $(DEPS)
	$(CC) $(CFLAGS) $(WFLAGS) $(ASANFLAGS) $(CPPFLAGS) $< -c

ensim.o: ensim.cc $(DEPS)
	$(CC) $(CFLAGS) $(WFLAGS) $(ASANFLAGS) $(CPPFLAGS) $< -c
	$(DUMP) $@ > $(ASM)

clean:
	rm -f $(ASM) $(OUT) $(OBJS)
