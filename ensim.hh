#ifndef __ENSIM_HH__
#define __ENSIM_HH__

#include <cstddef>
#include <memory>
#include <vector>
#include <array>
#include <string_view>

namespace ensim
{

/*
 * Configure SIMD scheme depending on your platform.
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

using line = std::vector<real>;
using grid = std::vector<line>;

enum class type : size_t
{
    inline8
};

struct engine
{
    virtual void run(const size_t steps, const size_t x = -1, const size_t y = -1) = 0;
    virtual void reset() = 0;
    virtual size_t get_w() const = 0;
    virtual size_t get_h() const = 0;
    virtual size_t get_y() const = 0;
    virtual size_t bytes() const = 0;
    virtual real get_panic_status(const size_t x, const size_t y) const = 0;
    virtual real get_port_open_ratio(const size_t x, const size_t y) const = 0;
    virtual std::string_view get_signal_name(const size_t index) const = 0;
    virtual const line& get_signal(const size_t index) const = 0;
    virtual const line& get_static_temperature_signal_k() const = 0;
    virtual const line& get_static_pressure_signal_pa() const = 0;
    virtual const line& get_volume_signal_m3() const = 0;
    virtual ~engine() = default;
};

std::unique_ptr<engine> new_engine(const type);

} /* namespace ensim */

#endif
