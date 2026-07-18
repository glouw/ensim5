#ifndef __ENSIM5_HH__
#define __ENSIM5_HH__

#include <cstddef>
#include <memory>

class ensim5_c
{
public:
    virtual void run_engine(const size_t steps) = 0;
    virtual void reset_engine() = 0;
    virtual ~ensim5_c() = default;
};

std::unique_ptr<ensim5_c> new_engine_8_9_5();

#endif
