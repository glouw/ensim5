#include "ensim5.hh"
#include <iostream>
#include <chrono>

int main(int argc, char* argv[])
{
    size_t steps = atoi(argv[1]);
    auto t0 = std::chrono::high_resolution_clock::now();
    auto engine = new_engine_8_9_5();
    engine->run_engine(steps);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    std::cout << "Time: " << dt.count() << "ms\n";
}

// 44800
