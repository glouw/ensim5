#ifndef __ENSIM_HH__
#define __ENSIM_HH__

#include <cstddef>
#include <memory>
#include <vector>
#include <array>

namespace ensim
{
    struct diags
    {
        #define ENSIM_FLUIDS_LIST(X)          \
            X(chamber_volume_m3)              \
            X(chamber_nozzle_open_ratio)      \
            X(chamber_mass_kg)                \
            X(chamber_static_pressure_pa)     \
            X(chamber_static_temperature_k)   \
            X(nozzle_mach)                    \
            X(nozzle_mass_flow_rate_kg_per_s) \

        #define ENSIM_PISTONS_LIST(X) \
            X(gas_torque_n_m)         \
            X(inertia_torque_n_m)     \

        #define ENSIM_DIAGS_LIST(X) \
            ENSIM_FLUIDS_LIST(X)    \
            ENSIM_PISTONS_LIST(X)   \

        enum class channel : int
        {
            #define X(name) name,
            ENSIM_DIAGS_LIST(X)
            #undef X
            count,
        };

        static constexpr int channels = static_cast<int>(channel::count);

        static constexpr std::array<const char*, channels> name = {
            #define X(name) #name,
            ENSIM_DIAGS_LIST(X)
            #undef X
        };

    void clear() { for(auto& plot : plots) plot.clear(); }
    std::vector<double>& operator[](const int index) { return plots[index]; }
    std::vector<double>& operator[](const channel channel) { return plots[static_cast<int>(channel)]; }
    const std::vector<double>& operator[](const int index) const { return plots[index]; }
    const std::vector<double>& operator[](const channel channel) const { return plots[static_cast<int>(channel)]; }

    private:
        std::array<std::vector<double>, channels> plots;
    };

    struct engine
    {
        enum class type : int { inline8 };
        virtual void run(const unsigned steps, const unsigned x = -1, const unsigned y = -1) = 0;
        virtual void reset() = 0;
        virtual int width() = 0;
        virtual int height() = 0;
        virtual int piston_y() = 0;
        virtual size_t bytes() = 0;
        virtual const diags& get_diags() const = 0;
        virtual std::vector<std::vector<double>> get_port_open_ratios() = 0;
        virtual ~engine() = default;
    };

    std::unique_ptr<engine> new_engine(const engine::type type);
}

#endif
