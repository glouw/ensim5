#ifndef __ENSIM5_HH__
#define __ENSIM5_HH__

#include <cstddef>
#include <memory>
#include <vector>
#include <array>

enum class ensim5_diag_e : int
{
    chamber_volume_m3               = 0,
    chamber_static_pressure_pa      = 1,
    chamber_static_temperature_k    = 2,
    chamber_molar_mass_kg_per_mol   = 3,
    chamber_cv_j_per_mol_k          = 4,
    chamber_moles                   = 5,
    chamber_nozzle_flow_area_m2     = 6,
    nozzle_velocity_m_per_s         = 7,
    nozzle_static_density_kg_per_m3 = 8,
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
