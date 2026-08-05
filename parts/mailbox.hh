#pragma once

template<size_t W, size_t H>
struct alignas(std::hardware_destructive_interference_size) mailbox
{
    /*
     * Recieve
     */

    atom throttle_open_ratio = 0.5_r;
    std::atomic<size_t> log_x = -1;
    std::atomic<size_t> log_y = -1;
    trip injection_enabled = true;

    /*
     * Send
     */

    std::atomic<size_t> swap_drops = 0;
    atom engine_angular_velocity_r_per_s = 0.0_r;
    std::array<std::array<atom, W>, H> port_open_ratios = {};
    std::array<std::array<trip, W>, H> panics = {};
};
