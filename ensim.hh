#pragma once

#include <cstddef>
#include <memory>
#include <vector>
#include <array>
#include <atomic>
#include <span>
#include <string_view>

namespace ensim
{
    static constexpr size_t g_sample_rate_hz = 48000;

    struct engine
    {
        virtual void run(const size_t steps) = 0;
        virtual void reset() = 0;
        virtual size_t get_width() const = 0;
        virtual size_t get_height() const = 0;
        virtual size_t get_piston_y() const = 0;
        virtual size_t get_audio_y() const = 0;
        virtual size_t get_throttle_y() const = 0;
        virtual size_t get_bytes() const = 0;
        virtual std::string_view get_signal_name(const size_t index) const = 0;
        virtual const std::atomic<double>& get_load_angular_velocity_r_per_s() const = 0;
        virtual const std::atomic<double>& get_engine_angular_velocity_r_per_s() const = 0;
        virtual const std::atomic<double>& get_port_open_ratio(const size_t x, const size_t y) const = 0;
        virtual const std::atomic<bool>& get_panic(const size_t x, const size_t y) const = 0;
        virtual size_t get_swap_drops() const = 0;
        virtual const std::vector<double>& get_signal(const size_t index) const = 0;
        virtual const std::vector<double>& get_static_temperature_signal_k() const = 0;
        virtual const std::vector<double>& get_static_pressure_signal_pa() const = 0;
        virtual const std::vector<double>& get_volume_signal_m3() const = 0;
        virtual const std::vector<float>& get_audio_signal() const = 0;
        virtual const std::vector<float>& get_impulse_signal() const = 0;
        virtual const std::span<const std::vector<float>> get_pipe_pressure_signals() const = 0;
        virtual void set_throttle_open_ratio(const double) = 0;
        virtual void set_injection_on() = 0;
        virtual void set_injection_off() = 0;
        virtual void increment_gear() = 0;
        virtual void decrement_gear() = 0;
        virtual void engage_clutch() = 0;
        virtual void disengage_clutch() = 0;
        virtual void set_logger(const size_t x, const size_t y) = 0;
        virtual void set_swap_lock_on() = 0;
        virtual void set_swap_lock_off() = 0;

        virtual ~engine() = default;
    };

    enum class type : size_t
    {
        inline8
    };

    std::unique_ptr<engine> new_engine(const type);
}
