#pragma once

#define FLUIDS(X) \
    X(chamber_volume_m3) \
    X(chamber_nozzle_real_flow_area_m2) \
    X(chamber_static_pressure_pa) \
    X(chamber_static_temperature_k) \
    X(chamber_mass_kg) \
    X(nozzle_mass_flow_rate_kg_per_s) \
    X(nozzle_velocity_m_per_s)

#define PISTONS(X) \
    X(gas_torque_n_m) \
    X(inertia_torque_n_m)

#define DIAGS(X) FLUIDS(X) PISTONS(X)

enum
{
    #define X(name) g_##name,
    DIAGS(X)
    #undef X
    g_diags_size,
};

static constexpr std::array<std::string_view, g_diags_size> g_signal_names = {
    #define X(name) #name,
    DIAGS(X)
    #undef X
};

struct diags
{
    grid front = grid(g_diags_size);
    grid back = grid(g_diags_size);
};
