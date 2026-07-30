#ifndef __ENSIM_HH__
#define __ENSIM_HH__

#include <cstddef>
#include <memory>
#include <vector>
#include <array>
#include <string_view>

namespace ensim
{

static constexpr size_t sample_rate_hz = 48000;

/*
 * Configure SIMD scheme depending on your platform.
 * Float stability is not guaranteed.
 *
 */

#if 0
using real = float;
#else
using real = double;
#endif

consteval real operator""_r(const long double x)
{
    return static_cast<real>(x);
}

using atom = std::atomic<real>;
using line = std::vector<real>;
using grid = std::vector<line>;

enum class type : size_t
{
    inline8
};

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
    virtual const atom& get_angular_velocity_r_per_s() const = 0;
    virtual const atom& get_port_open_ratio(const size_t x, const size_t y) const = 0;
    virtual size_t get_swap_drops() const = 0;

    virtual const line& get_signal(const size_t index) const = 0;
    virtual const line& get_temperature_signal_k() const = 0;
    virtual const line& get_pressure_signal_pa() const = 0;
    virtual const line& get_volume_signal_m3() const = 0;
    virtual const line& get_audio_signal() const = 0;
    virtual const line& get_impulse_signal() const = 0;
    virtual const std::vector<float> get_audio_data() const = 0;
    virtual void set_throttle_open_ratio(const real open_ratio) = 0;
    virtual void set_injection_on() = 0;
    virtual void set_injection_off() = 0;
    virtual void set_logger(const size_t x, const size_t y) = 0;
    virtual void set_swap_lock_on() = 0;
    virtual void set_swap_lock_off() = 0;

    virtual ~engine() = default;
};

std::unique_ptr<engine> new_engine(const type);

} /* namespace ensim */

#endif
