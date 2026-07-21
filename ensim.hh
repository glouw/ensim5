#ifndef __ENSIM5_HH__
#define __ENSIM5_HH__

#include <cstddef>
#include <memory>
#include <vector>
#include <array>

namespace ensim
{
    struct diags
    {
        #define ENSIM_DIAGS_LIST(X)            \
            X(chamber_volume_m3)               \
            X(chamber_nozzle_flow_area_m2)     \
            X(chamber_static_pressure_pa)      \
            X(chamber_static_temperature_k)    \
            X(chamber_mass_kg)                 \
            X(nozzle_mach)                     \
            X(nozzle_velocity_m_per_s)         \
            X(nozzle_static_density_kg_per_m3) \
            X(nozzle_mass_flow_rate_kg_per_s)

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

        void clear();
        std::vector<float>& operator[](const int index);
        std::vector<float>& operator[](const channel channel);

    private:
        std::array<std::vector<float>, channels> plots;
    };

    struct engine
    {
        enum class type : int { inline8 };
        virtual void run(const int steps, const int x=-1, const int y=-1) = 0;
        virtual void reset() = 0;
        virtual int get_w() = 0;
        virtual int get_h() = 0;
        virtual int get_y() = 0;
        virtual diags& get_diags() = 0;
        virtual size_t get_size() = 0;
        virtual ~engine() = default;
    };

    std::unique_ptr<engine> new_engine(const engine::type type);
}

#endif
