#ifndef __ENSIM5_HH__
#define __ENSIM5_HH__

#include <cstddef>
#include <memory>
#include <vector>
#include <array>

enum class ensim5_diag_e : int
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
};

struct ensim5_diag_s
{
    std::array<std::vector<float>, 9> plots;

    void clear()
    {
        for(auto& plot : plots)
        {
            plot.clear();
        }
    }

    std::vector<float>& operator[](const int index)
    {
        return plots[index];
    }

    std::vector<float>& operator[](const ensim5_diag_e index)
    {
        return plots[static_cast<int>(index)];
    }
};

struct ensim5_s
{
    virtual void run_engine(const int steps, const int x=-1, const int y=-1) = 0;
    virtual void reset_engine() = 0;
    virtual int get_w() = 0;
    virtual int get_h() = 0;
    virtual int get_y() = 0;
    virtual ensim5_diag_s& get_diags() = 0;
    virtual ~ensim5_s() = default;
};

std::unique_ptr<ensim5_s> new_ensim5_inline_8();

#endif
