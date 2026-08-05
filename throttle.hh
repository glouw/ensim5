#pragma once

struct throttle
{
    static constexpr size_t size = 4;
    std::array<real, size> table = {};

    real lookup(const real open_ratio)
    {
        const size_t last = size - 1;
        const real at = last * open_ratio;
        const size_t index = at;
        const real ratio = at - index;
        const size_t next = index + 1;
        const real delta = table[next] - table[index];
        return table[index] + delta * ratio;
    }
};
