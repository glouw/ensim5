#pragma once

template<size_t W, size_t H>
struct alignas(std::hardware_destructive_interference_size) mailbox
{
    /*
     * Recieve
     */

    std::atomic<real> throttle_open_ratio = 0.5_r;
    std::atomic<size_t> log_x = -1;
    std::atomic<size_t> log_y = -1;
    std::atomic<bool> injection_enabled = true;

    /*
     * Send
     */

    std::atomic<size_t> swap_drops = 0;
    std::atomic<real> engine_angular_velocity_r_per_s = 0.0_r;
    std::array<std::array<std::atomic<real>, W>, H> port_open_ratios = {};
    std::array<std::array<std::atomic<bool>, W>, H> panics = {};
};
