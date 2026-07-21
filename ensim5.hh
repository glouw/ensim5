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
        enum class channel : int
        {
            chamber_volume_m3,
            chamber_nozzle_flow_area_m2,
            chamber_static_pressure_pa,
            chamber_static_temperature_k,
            chamber_mass_kg,
            nozzle_mach,
            nozzle_velocity_m_per_s,
            nozzle_static_density_kg_per_m3,
            nozzle_mass_flow_rate_kg_per_s,
            count,
        };
        void clear();
        std::vector<float>& operator[](const int index);
        std::vector<float>& operator[](const channel channel);

    private:
        std::array<std::vector<float>, static_cast<int>(channel::count)> plots;
    };

    struct engine
    {
        enum class type : int
        {
            inline8
        };

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
