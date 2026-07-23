#ifndef __ENSIM_HH__
#define __ENSIM_HH__

#include <cstddef>
#include <memory>
#include <vector>
#include <array>
#include <string_view>

#define ENSIM_FLUIDS_LIST(X)          \
    X(chamber_volume_m3)              \
    X(chamber_nozzle_open_ratio)      \
    X(nozzle_mass_flow_rate_kg_per_s) \
    X(chamber_mass_kg)                \
    X(chamber_static_pressure_pa)     \
    X(chamber_static_temperature_k)   \

#define ENSIM_SPARKPLUG_LIST(X)       \
    X(fired)

#define ENSIM_PISTONS_LIST(X)         \
    X(gas_torque_n_m)                 \
    X(inertia_torque_n_m)             \

#define ENSIM_DIAGS_LIST(X)           \
    ENSIM_FLUIDS_LIST(X)              \
    ENSIM_SPARKPLUG_LIST(X)           \
    ENSIM_PISTONS_LIST(X)             \

namespace ensim
{

/*
 * Configure SIMD scheme depending on your platform.
 *
 */

#if 1
using real = float;
#else
using real = double;
#endif

consteval real operator""_r(const long double x)
{
    return static_cast<real>(x);
}

struct diags
{
    enum class channel : size_t
    {
        #define X(name) name,
        ENSIM_DIAGS_LIST(X)
        #undef X
        count,
    };

    static constexpr size_t channels = static_cast<size_t>(channel::count);

    static constexpr std::array<std::string_view, channels> name = {
        #define X(name) #name,
        ENSIM_DIAGS_LIST(X)
        #undef X
    };

    std::vector<real>& operator[](const channel c)
    {
        return plots[static_cast<size_t>(c)];
    }

    const std::vector<real>& operator[](const channel c) const
    {
        return plots[static_cast<size_t>(c)];
    }

    void clear()
    {
        for(auto& plot : plots)
        {
            plot.clear();
        }
    }

private:
    std::array<std::vector<real>, channels> plots;
};

struct engine
{
    enum class type : size_t { inline8 };
    virtual void run(const size_t steps, const size_t x = -1, const size_t y = -1) = 0;
    virtual void reset() = 0;
    virtual size_t get_width() = 0;
    virtual size_t get_height() = 0;
    virtual size_t get_piston_y() = 0;
    virtual size_t get_bytes() = 0;
    virtual const diags& get_diags() const = 0;
    virtual std::vector<std::vector<real>> get_panics() = 0;
    virtual std::vector<std::vector<real>> get_port_open_ratios() = 0;
    virtual ~engine() = default;
};

std::unique_ptr<engine> new_engine(const engine::type type);

} /* namespace ensim */

#endif
