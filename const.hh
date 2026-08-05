#pragma once

#define fn __attribute__((used))

static constexpr real g_dt_s = 1.0_r / g_sample_rate_hz;
static constexpr real g_pi_r = std::numbers::pi_v<real>;
static constexpr real g_otto_cycle_r = 4.0_r * g_pi_r;
static constexpr real g_otto_intake_cycle_r = 0.0_r * g_pi_r;
static constexpr real g_otto_compression_cycle_r = 1.0_r * g_pi_r;
static constexpr real g_otto_combustion_cycle_r = 2.0_r * g_pi_r;
static constexpr real g_otto_exhaust_cycle_r = 3.0_r * g_pi_r;
static constexpr real g_resevoir_volume_m3 = 1e9_r;
static constexpr real g_ambient_temperature_k = 300.0_r;
static constexpr real g_ambient_pressure_pa = 132'800.0_r;
static constexpr real g_ambient_density_kg_per_m3 = 1.225_r;
static constexpr real g_gamma = 3.0_r / 2.0_r;
static constexpr real g_universal_gas_constant_j_per_mol_k = 8.3144598_r;
static constexpr real g_cv_j_per_mol_k = g_universal_gas_constant_j_per_mol_k / (g_gamma - 1.0_r);
static constexpr real g_molar_mass_kg_per_mol = 0.023_r;
static constexpr real g_cv_j_per_kg_k = g_cv_j_per_mol_k / g_molar_mass_kg_per_mol;
static constexpr real g_specific_gas_constant_j_per_kg_k = g_universal_gas_constant_j_per_mol_k / g_molar_mass_kg_per_mol;
static constexpr real g_energy_octane_j_per_kg = 47.9e6_r;
static constexpr real g_stoich_air_fuel_ratio = 14.7_r;
