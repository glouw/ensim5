#include <iostream>
#include <chrono>

#include "ensim5.hh"

int main(int argc, char* argv[])
{
    auto t0 = std::chrono::high_resolution_clock::now();
    size_t steps = atoi(argv[1]);
    auto engine = ensim5::init_engine_8_9_5();
    for(size_t step = 0; step < steps; step++)
    {
        engine.step_engine();
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    std::cout << "Time: " << dt.count() << "us\n";
}

// 44800
