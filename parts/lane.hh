#pragma once

template<size_t N>
using lane = std::array<real, N>;

template<size_t N>
using mask = std::array<bool, N>;
