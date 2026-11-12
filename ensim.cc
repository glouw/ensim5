#include "ensim.hh"

#include <array>
#include <numbers>
#include <cmath>
#include <mutex>
#include <cassert>
#include <atomic>

#define fn __attribute__((used))

namespace ensim {

constexpr size_t g_cacheline_size = std::hardware_destructive_interference_size;
constexpr double g_dt_s = 1.0 / g_sample_rate_hz;
constexpr double g_pi_r = std::numbers::pi_v<double>;
constexpr double g_otto_cycle_r = 4.0 * g_pi_r;
constexpr double g_otto_intake_cycle_r = 0.0 * g_pi_r;
constexpr double g_otto_combustion_cycle_r = 2.0 * g_pi_r;
constexpr double g_otto_exhaust_cycle_r = 3.0 * g_pi_r;
constexpr double g_ambient_temperature_k = 300.0;
constexpr double g_ambient_pressure_pa = 132'800.0;
constexpr double g_ambient_density_kg_per_m3 = 1.22453421;
constexpr double g_ambient_speed_of_sound_m_per_s = 403.328774;
constexpr double g_gamma = 1.5;
constexpr double g_cv_j_per_kg_k = 722.996;
constexpr double g_specific_gas_constant_j_per_kg_k = 361.498;
constexpr double g_energy_octane_j_per_kg = 47.9e6;
constexpr double g_stoich_air_fuel_ratio = 14.7;
constexpr float g_pi_r_f = g_pi_r;
constexpr float g_gamma_f = g_gamma;
constexpr float g_dt_s_f = g_dt_s;

using std::sin, std::cos, std::fmax, std::fmin, std::log, std::sqrt, std::trunc, std::exp, std::fabs, std::abs, std::llround;

template<typename T>
fn inline T clamper(const T value, const T lower, const T upper)
{
    return fmax(fmin(value, upper), lower);
}

template<typename T>
fn inline T modulos(const T value, const T by)
{
    return value - trunc(value / by) * by;
}

template<typename T>
fn inline T cuberoot(const T value)
{
    return exp(log(value) / static_cast<T>(3.0));
}

template<size_t H, size_t PY>
struct flow
{
    static constexpr size_t N = H - 1;
    static_assert(N % 2 == 0);

    std::array<double, H> chamber_prev_volume_m3 = {};
    std::array<double, H> chamber_volume_m3 = {};
    std::array<double, H> chamber_nozzle_flow_area_m2 = {};
    std::array<double, H> chamber_nozzle_real_flow_area_m2 = {};
    std::array<double, H> chamber_nozzle_open_ratio = {};
    std::array<double, H> chamber_static_pressure_pa = {};
    std::array<double, H> chamber_dynamic_pressure_pa = {};
    std::array<double, H> chamber_total_pressure_pa = {};
    std::array<double, H> chamber_static_temperature_k = {};
    std::array<double, H> chamber_dynamic_temperature_k = {};
    std::array<double, H> chamber_total_temperature_k = {};
    std::array<double, H> chamber_mass_kg = {};
    std::array<double, H> chamber_bulk_momentum_kg_m_per_s = {};
    std::array<double, H> nozzle_mach = {};
    std::array<double, H> nozzle_velocity_m_per_s = {};
    std::array<double, H> nozzle_static_density_kg_per_m3 = {};
    std::array<double, H> nozzle_static_temperature_k = {};
    std::array<double, H> nozzle_mass_flow_rate_kg_per_s = {};
    std::array<double, H> parcel_mass_kg = {};
    std::array<double, H> parcel_static_temperature_k = {};
    std::array<bool, H> panic = {};
    double piston_injector_enabled = 0.0;
    double piston_chamber_flame_height_m = 0.0;
    double piston_chamber_mass_burned_m3 = 0.0;
    double piston_chamber_radius_m = 0.0;
    bool piston_chamber_on_fire = false;

    /*
     *     Ps * V
     * m = -------
     *     Rs * Ts
     */

    fn void calc_chamber_ambients()
    {
        for(size_t i = 0; i < H; i++)
        {
            chamber_static_temperature_k[i] = g_ambient_temperature_k;
            chamber_static_pressure_pa[i] = g_ambient_pressure_pa;
        }
        for(size_t i = 0; i < H; i++)
        {
            const double Ps = chamber_static_pressure_pa[i];
            const double V = chamber_volume_m3[i];
            const double Rs = g_specific_gas_constant_j_per_kg_k;
            const double Ts = chamber_static_temperature_k[i];
            chamber_mass_kg[i] = (Ps * V) / (Rs * Ts);
        }
    }

    /*
     *       m * Rs * T
     * Ps = ------------
     *           V
     */

    fn void calc_chamber_static_pressures()
    {
        for(size_t i = 0; i < N; i++)
        {
            const double m = chamber_mass_kg[i];
            const double Rs = g_specific_gas_constant_j_per_kg_k;
            const double Ts = chamber_static_temperature_k[i];
            const double V = chamber_volume_m3[i];
            chamber_static_pressure_pa[i] = m * Rs * Ts / V;
        }
    }

    /*              ____________________
     *             /
     *            /            y - 1
     *           /             -----
     *          /                y
     *         /   2         Pt
     * M = _  /  ----- * [ (----) - 1 ]
     *      \/   y - 1       Ps
     */

    fn void calc_nozzle_machs()
    {
        for(size_t i = 0; i < N; i++)
        {
            const size_t j = i + 1;
            const double X = (2.0 / (g_gamma - 1.0));
            const double Pi = chamber_total_pressure_pa[i];
            const double Pj = chamber_total_pressure_pa[j];
            const double Pt = fmax(Pi, Pj);
            const double Ps = fmin(Pi, Pj);
            const double direction = Pi > Pj ? 1.0 : -1.0;

            /*
             *      y - 1                        3
             *      ----- = 0.3333... where y = ---
             *        y                          2
             * Term
             */

            static_assert(g_gamma == 1.5);
            const double Y = cuberoot(Pt / Ps);
            const double M = direction * sqrt(X * (Y - 1.0));
            const double A = chamber_nozzle_real_flow_area_m2[i];
            const double mute = A == 0.0 ? 0.0 : 1.0;
            nozzle_mach[i] = mute * M;
        }
    }

    /*
     *
     * Pt = Ps + Pd
     *
     */

    fn void calc_chamber_total_pressures()
    {
        for(size_t i = 0; i < H; i++)
        {
            const double Ps = chamber_static_pressure_pa[i];
            const double Pd = chamber_dynamic_pressure_pa[i];
            chamber_total_pressure_pa[i] = Ps + Pd;
        }
    }

    /*
     *
     * Tt = Ts + Td
     *
     */

    fn void calc_chamber_total_temperatures()
    {
        for(size_t i = 0; i < H; i++)
        {
            const double Ts = chamber_static_temperature_k[i];
            const double Td = chamber_dynamic_temperature_k[i];
            chamber_total_temperature_k[i] = Ts + Td;
        }
    }

    /*
     *               ______________
     *              /
     *             /  y * Pt / pt
     *            /  -------------
     *           /       y - 1  2
     * u = M _  /    1 + ----- M
     *        \/           2
     *
     *        Pt
     * pt = -------
     *      Rs * Tt
     *
     */

    fn void calc_nozzle_velocities()
    {
        for(size_t i = 0; i < N; i++)
        {
            const double Rs = g_specific_gas_constant_j_per_kg_k;
            const double Tt = chamber_total_temperature_k[i];
            const double M = nozzle_mach[i];
            const double X = g_gamma * Rs * Tt;
            const double Y = 0.5 * (g_gamma - 1.0) * M * M;
            const double u = M * sqrt(X / (1.0 + Y));
            nozzle_velocity_m_per_s[i] = u;
        }
    }

    /*
     *              pt
     * ps = -------------------
     *                      1
     *                    -----
     *                    y - 1
     *           y - 1  2
     *      (1 + ----- M )
     *             2
     */

    fn void calc_nozzle_static_densities()
    {
        for(size_t i = 0; i < N; i++)
        {
            const double Pt = chamber_total_pressure_pa[i];
            const double Rs = g_specific_gas_constant_j_per_kg_k;
            const double Tt = chamber_total_temperature_k[i];
            const double M = nozzle_mach[i];
            const double pt = Pt / (Rs * Tt);

            /*
             *        1                  3
             *      ----- = 2 where y = ---
             *      y - 1                2
             * Term
             */

            static_assert(g_gamma == 1.5);
            const double C = 1.0 + 0.5 * (g_gamma - 1.0) * M * M;
            nozzle_static_density_kg_per_m3[i] = pt / (C * C);
        }
    }

    /*
     *
     *              Tt
     * Ts = ------------------
     *            (y - 1)  2
     *        1 + ------- M
     *               2
     */

    fn void calc_nozzle_static_temperatures()
    {
        for(size_t i = 0; i < N; i++)
        {
            const double Tt = chamber_total_temperature_k[i];
            const double M = nozzle_mach[i];
            const double Tns = Tt / (1.0 + 0.5 * (g_gamma - 1.0) * M * M);
            nozzle_static_temperature_k[i] = Tns;
        }
    }

    fn void calc_nozzle_real_flow_areas()
    {
        for(size_t i = 0; i < N; i++)
        {
            const double r = chamber_nozzle_open_ratio[i];
            const double A = chamber_nozzle_flow_area_m2[i];
            chamber_nozzle_real_flow_area_m2[i] = r * A;
        }
    }

    /* .
     * m = ps A u
     *
     */

    fn void calc_nozzle_mass_flow_rates()
    {
        for(size_t i = 0; i < N; i++)
        {
            const double ps = nozzle_static_density_kg_per_m3[i];
            const double A = chamber_nozzle_real_flow_area_m2[i];
            const double u = nozzle_velocity_m_per_s[i];
            const double mdot = ps * A * u;
            nozzle_mass_flow_rate_kg_per_s[i] = mdot;
        }
    }

    /*
     *      .
     * mp = m dt
     * Tsp = Ts,upstream
     *
     */

    fn void calc_nozzle_parcels()
    {
        for(size_t i = 0; i < N; i++)
        {
            const size_t j = i + 1;
            const double mdot = nozzle_mass_flow_rate_kg_per_s[i];
            const double dm = mdot * g_dt_s;
            parcel_mass_kg[i] = dm;
            const double Tsi = chamber_static_temperature_k[i];
            const double Tsj = chamber_static_temperature_k[j];
            const double Tsp = mdot > 0.0 ? Tsi : Tsj;
            parcel_static_temperature_k[i] = Tsp;
        }
    }

    /*                  y - 1
     *               V1
     * Ts2 = Ts1 * (----)
     *               V2
     */

    fn void calc_compressions()
    {
        for(size_t i = 0; i < N; i++)
        {
            const double Ts1 = chamber_static_temperature_k[i];
            const double V1 = chamber_prev_volume_m3[i];
            const double V2 = chamber_volume_m3[i];
            const double dv = V1 / V2;

            /*
             *               1             3
             *      y - 1 = --- where y = ---
             * Term          2             2
             *
             */

            static_assert(g_gamma == 1.5);
            chamber_static_temperature_k[i] = Ts1 * sqrt(dv);
        }
    }

    /*
     *       Ts m + Tp dm
     * Ts = --------------
     *          m + dm
     */

    fn void calc_forward_energy_transfers()
    {
        for(size_t i = 0; i < N; i++)
        {
            const size_t j = i + 1;
            const double dm = parcel_mass_kg[i];
            const double m = chamber_mass_kg[j];
            const double Tsp = parcel_static_temperature_k[i];
            const double Ts0 = chamber_static_temperature_k[j];
            const double Ts1 = (Ts0 * m + Tsp * dm) / (m + dm);
            chamber_static_temperature_k[j] = dm > 0.0 ? Ts1 : Ts0;
        }
    }

    fn void calc_reverse_energy_transfers()
    {
        for(size_t i = N; i > 0; i--)
        {
            const double dm = parcel_mass_kg[i];
            const double m = chamber_mass_kg[i];
            const double Tsp = parcel_static_temperature_k[i];
            const double Ts0 = chamber_static_temperature_k[i];
            const double Ts1 = (Ts0 * m - Tsp * dm) / (m - dm);
            chamber_static_temperature_k[i] = dm < 0.0 ? Ts1 : Ts0;
        }
    }

    /*                  y - 1
     *               m1
     * Ts2 = Ts1 * (----)
     *               m2
     */

    fn void calc_mass_transfers()
    {
        for(size_t i = 0; i < N; i++)
        {
            const size_t j = i + 1;
            const double mi = chamber_mass_kg[i];
            const double mj = chamber_mass_kg[j];
            const double dm = parcel_mass_kg[i];
            const double m0 = mi - dm;
            const double m1 = mj + dm;

            /*
             *               1             3
             *      y - 1 = --- where y = ---
             * Term          2             2
             *
             */

            static_assert(g_gamma == 1.5);
            chamber_static_temperature_k[i] *= sqrt(m0 / mi);
            chamber_static_temperature_k[j] *= sqrt(m1 / mj);
            chamber_mass_kg[i] = m0;
            chamber_mass_kg[j] = m1;
        }

        for(size_t i = 0; i < N; i++)
        {
            const size_t j = i + 1;
            const double dm = parcel_mass_kg[i];
            const double u = nozzle_velocity_m_per_s[i];
            const double p = dm * u;
            chamber_bulk_momentum_kg_m_per_s[i] -= p;
            chamber_bulk_momentum_kg_m_per_s[j] += p;
        }

        /*               ___________
         *              /
         * pmax = m _  / y * Rs * Ts
         *           \/
         */

        for(size_t i = 0; i < N; i++)
        {
            const double Rs = g_specific_gas_constant_j_per_kg_k;
            const double Ts = chamber_static_temperature_k[i];
            const double m = chamber_mass_kg[i];
            const double p = chamber_bulk_momentum_kg_m_per_s[i];
            const double pmax = m * sqrt(g_gamma * Rs * Ts);
            chamber_bulk_momentum_kg_m_per_s[i] = clamper(p, -pmax, pmax);
        }
    }

    /*
     *      1     2
     * q = --- p u
     *      2
     */

    fn void calc_chamber_dynamic_pressures()
    {
        for(size_t i = 0; i < N; i++)
        {
            const double u = chamber_bulk_momentum_kg_m_per_s[i] / chamber_mass_kg[i];
            const double p = chamber_mass_kg[i] / chamber_volume_m3[i];
            const double q = 0.5 * p * u * u;
            chamber_dynamic_pressure_pa[i] = q;
        }
    }

    /*         2
     *        u
     * Td = ------
     *       2 Cp
     */

    fn void calc_chamber_dynamic_temperatures()
    {
        for(size_t i = 0; i < N; i++)
        {
            const double u = chamber_bulk_momentum_kg_m_per_s[i] / chamber_mass_kg[i];
            const double Cv = g_cv_j_per_kg_k;
            const double Cp = g_gamma * Cv;
            const double Td = 0.5 * u * u / Cp;
            chamber_dynamic_temperature_k[i] = Td;
        }
    }

    fn void ignite_piston_chamber()
    {
        piston_chamber_flame_height_m = 0.0;
        piston_chamber_mass_burned_m3 = 0.0;
        piston_chamber_on_fire = true;
    }

    fn void calc_combustion()
    {
        if(piston_chamber_on_fire)
        {
            if(not piston_injector_enabled)
            {
                piston_chamber_on_fire = false;
                return;
            }
            const double M = chamber_mass_kg[PY];
            const double V = chamber_volume_m3[PY];

            /*
             *              Ts    2
             *           [------]
             *             Ts0
             * S =  0.4 ----------------
             *              Ps    0.125
             *           [------]
             *             Ps0
             */

            const double Ts = chamber_static_temperature_k[PY];
            const double Ts0 = g_ambient_temperature_k;
            const double Ps = chamber_static_pressure_pa[PY];
            const double Ps0 = g_ambient_pressure_pa;
            const double Tr = Ts / Ts0;
            const double Pr = Ps/ Ps0;
            const double S = 0.4 * Tr * Tr / sqrt(sqrt(sqrt(Pr)));

            /*
             *                 2
             * Vb = dh * pi * r
             *
             */

            const double dh = S * g_dt_s;
            const double h1 = piston_chamber_flame_height_m;
            const double h2 = h1 + dh;
            const double r = piston_chamber_radius_m;
            const double Vb = (h2 - h1) * g_pi_r * r * r;

            /*
             *      M
             * p = ---
             *      V
             *
             * Mburned = Vb * p
             *
             */

            const double p = M / V;
            const double Mb = Vb * p;

            /*
             *             Q
             * Ts = Ts + ------
             *            M Cv
             *
             */

            const double TMb = piston_chamber_mass_burned_m3 + Mb;
            if(TMb / M < 1.0)
            {
                const double MFb = Mb / (1.0 + g_stoich_air_fuel_ratio);
                const double Q = MFb * g_energy_octane_j_per_kg;
                const double Cv = g_cv_j_per_kg_k;
                const double dTs = Q / (M * Cv);
                chamber_static_temperature_k[PY] += dTs;
            }
            else
            {
                piston_chamber_on_fire = false;
            }
            piston_chamber_flame_height_m = h2;
            piston_chamber_mass_burned_m3 = TMb;
        }
    }

    fn void calc_panics()
    {
        for(size_t i = 0; i < N; i++)
        {
            panic[i] |= chamber_mass_kg[i] <= 0.0;
            panic[i] |= chamber_static_temperature_k[i] <= 0.0;
            panic[i] |= chamber_static_pressure_pa[i] <= 0.0;
        }
    }

    fn void update()
    {
        calc_chamber_dynamic_pressures();
        calc_chamber_dynamic_temperatures();
        calc_chamber_static_pressures();
        calc_chamber_total_pressures();
        calc_chamber_total_temperatures();
        calc_nozzle_real_flow_areas();
        calc_nozzle_machs();
        calc_nozzle_velocities();
        calc_nozzle_static_densities();
        calc_nozzle_mass_flow_rates();
        calc_nozzle_static_temperatures();
        calc_nozzle_parcels();
        calc_compressions();
        calc_forward_energy_transfers();
        calc_reverse_energy_transfers();
        calc_mass_transfers();
        calc_combustion();
        calc_panics();
    }
};

template<size_t W>
struct sparkplugs
{
    static constexpr double fire_delay_theta_r = 1e-1;
    std::array<double, W> engage_theta_r = {};
    std::array<bool, W> prev_fired = {};
    std::array<bool, W> fired = {};
    std::array<bool, W> rising_edge = {};
    double flywheel_theta_r = 0.0;

    fn void calc_fired()
    {
        for(size_t i = 0; i < W; i++)
        {
            prev_fired[i] = fired[i];
            const double theta0_r = modulos(flywheel_theta_r, g_otto_cycle_r);
            const double theta1_r = engage_theta_r[i] >= g_otto_cycle_r ? (engage_theta_r[i] - g_otto_cycle_r) : engage_theta_r[i];
            fired[i] = theta0_r > theta1_r + fire_delay_theta_r;
        }
    }

    fn void calc_rising_edge()
    {
        for(size_t i = 0; i < W; i++)
        {
            rising_edge[i] = fired[i] and not prev_fired[i];
        }
    }

    fn void update()
    {
        calc_fired();
        calc_rising_edge();
    }
};

template<size_t W>
struct cams
{
    std::array<double, W> engage_theta_r = {};
    std::array<double, W> ramp_theta_r = {};
    std::array<double, W> open_ratio = {};
    double flywheel_theta_r = 0.0;
    double flywheel_angular_velocity_r_per_s = 0.0;

    /*      4            1      2      3
     * r = t  [ 35 - 84 t + 70 t - 20 t ]
     *
     */

    fn void calc_open_ratios()
    {
        for(size_t i = 0; i < W; i++)
        {
            double theta0_r = modulos(engage_theta_r[i], g_otto_cycle_r);
            if(theta0_r < 0.0)
            {
                theta0_r += g_otto_cycle_r;
            }
            double theta1_r = modulos(flywheel_theta_r, g_otto_cycle_r);
            if(theta1_r < theta0_r)
            {
                theta1_r += g_otto_cycle_r;
            }
            const double open_r = theta1_r - theta0_r;
            const double t = open_r / ramp_theta_r[i];
            const double a = t * t * t * t;
            const double b = t * a;
            const double c = t * b;
            const double d = t * c;
            const double A = 35.0 * a;
            const double B = 84.0 * b;
            const double C = 70.0 * c;
            const double D = 20.0 * d;
            const double O = A - B + C - D;
            const double R = clamper(O, 0.0, 1.0);
            open_ratio[i] = theta1_r < theta0_r ? 0.0 : R;
        }
    }

    void update()
    {
        calc_open_ratios();
    }
};

template<size_t W>
struct vtec_cams: cams<W>
{
    double engage_angular_velocity_r_per_s = 0.0;
    double ramp_boost_ratio = 0.0;
    double gain_boost_ratio = 0.0;

    auto boost_ramp()
    {
        const double w0 = engage_angular_velocity_r_per_s;
        const double w1 = this->flywheel_angular_velocity_r_per_s;
        const auto backup = this->ramp_theta_r;
        for(auto& x : this->ramp_theta_r)
        {
            if(w1 > w0)
            {
                x *= ramp_boost_ratio;
            }
        }
        return backup;
    }

    void boost_gain()
    {
        const double w0 = engage_angular_velocity_r_per_s;
        const double w1 = this->flywheel_angular_velocity_r_per_s;
        for(auto& x : this->open_ratio)
        {
            if(w1 > w0)
            {
                x *= gain_boost_ratio;
                x = clamper(x, 0.0, 1.0);
            }
        }
    }

    void just_kicked_in_yo()
    {
        const auto backup = boost_ramp();
        cams<W>::calc_open_ratios();
        boost_gain();
        this->ramp_theta_r = backup;
    }

    void update()
    {
        just_kicked_in_yo();
    }
};

template<size_t N>
struct lerp_table
{
    std::array<double, N> table = {};

    double at(const double ratio) const
    {
        if(ratio < 0.0)
        {
            return table[0];
        }
        if(ratio >= 1.0)
        {
            return table[N - 1];
        }
        const double x = ratio * (N - 1);
        const size_t i = x;
        const double t = x - i;
        return table[i] + (table[i + 1] - table[i]) * t;
    }
};

template<size_t N>
struct throttle
{
    lerp_table<N> mapping = {};
    double open_ratio = 0.0;
};

struct limiter
{
    double max_angular_velocity_r_per_s = 0.0;
    double flywheel_angular_velocity_r_per_s = 0.0;
    double limit_time_s = 0.1;
    size_t cycles = 0;
    bool limiting = false;

    void update()
    {
        const double w0 = max_angular_velocity_r_per_s;
        const double w1 = flywheel_angular_velocity_r_per_s;
        if(not limiting)
        {
            if(w1 > w0)
            {
                limiting = true;
                cycles = 0;
            }
        }
        if(limiting)
        {
            const double time_s = g_dt_s * cycles;
            if(time_s > limit_time_s and w1 < w0)
            {
                limiting = false;
            }
            cycles++;
        }
    }
};

struct injector
{
    bool enabled = true;
};

struct disk
{
    double angular_velocity_r_per_s = 0.0;
    double mass_kg = 0.0;
    double radius_m = 0.0;
    double moment_of_inertia_kg_m2 = 0.0;
    double theta_r = 0.0;

    /*
     * dw = a * dt
     *
     */

    fn void accelerate(const double angular_acceleration_r_per_s_s)
    {
        const double a = angular_acceleration_r_per_s_s;
        angular_velocity_r_per_s += a * g_dt_s;
        angular_velocity_r_per_s = fmax(angular_velocity_r_per_s, 0.0);
        theta_r += angular_velocity_r_per_s * g_dt_s;
    }

    fn void calc_moment_of_inertia()
    {
        moment_of_inertia_kg_m2 = 0.5 * mass_kg * radius_m * radius_m;
    }
};

struct load: disk
{
    double rod_friction_n_m_s_per_r = 0.0;
    double brake_friction_n_m_s_per_r = 0.0;
    double friction_torque_n_m = 0.0;
    double total_torque_n_m = 0.0;
    double max_angular_velocity_r_per_s = 0.0;
    bool brake_engaged = false;

    /*
     * Tf = -K w
     *
     */

    fn void calc_friction_torque()
    {
        const double K1 = rod_friction_n_m_s_per_r;
        const double K2 = brake_engaged ? brake_friction_n_m_s_per_r : 0.0;
        const double w = angular_velocity_r_per_s;
        friction_torque_n_m = -w * (K1 + K2);
    }

    fn void calc_total_torque()
    {
        total_torque_n_m = friction_torque_n_m;
    }

    fn void update()
    {
        calc_moment_of_inertia();
        calc_friction_torque();
        calc_total_torque();
    }
};

struct crankshaft: disk
{
    fn void update()
    {
        calc_moment_of_inertia();
    }
};

struct flywheel: disk
{
    double theta_r = 0.0;
    double last_theta_r = 0.0;

    /*
     * dth = w * dt
     *
     */

    fn void turn()
    {
        last_theta_r = theta_r;
        const double w = angular_velocity_r_per_s;
        theta_r += w * g_dt_s;
    }

    fn bool otto_cycled()
    {
        const double t0 = modulos(last_theta_r, g_otto_cycle_r);
        const double t1 = modulos(theta_r, g_otto_cycle_r);
        return t0 > t1;
    }

    fn void update()
    {
        calc_moment_of_inertia();
    }
};

struct clutch
{
    double damping_coefficient_n_m_s = 0.0;
    bool engaged = false;

    double calc_torque(const double angular_velocity_flywheel_r_per_s, const double angular_velocity_transmission_r_per_s)
    {
        const double wl = angular_velocity_flywheel_r_per_s;
        const double wr = angular_velocity_transmission_r_per_s;
        const double Tk = engaged ? damping_coefficient_n_m_s : 0.0;
        const double Tc = Tk * (wl - wr);
        return Tc;
    }
};

template<size_t N>
struct gearbox
{
    static constexpr size_t last_gear = N - 1;
    std::array<double, N> ratios = {};
    size_t gear = 0;

    double lookup_gear_ratio()
    {
        return ratios[gear];
    }
};

/* ------- + block_deck_surface_m
 *         | head_clearance_height_m
 * ------- +
 * |     | | head_compression_height_m
 * |  o  | + pin_x_m, pin_y_m
 * |     | |
 * |-----| |
 *   | |   |
 *   | |   | connecting_rod_length_m
 *   | |   |
 *   | |   |
 *   | |   |
 *   |o|   + bearing_x_m, bearing_y_m
 *    |    |
 *    |    | crank_throw_length_m
 *    |    |
 *    o    + origin
 */

template<size_t W>
struct inline_pistons
{
    std::array<double, W> diameter_m = {};
    std::array<double, W> crank_throw_length_m = {};
    std::array<double, W> connecting_rod_length_m = {};
    std::array<double, W> connecting_rod_mass_kg = {};
    std::array<double, W> head_mass_density_kg_per_m3 = {};
    std::array<double, W> head_compression_height_m = {};
    std::array<double, W> head_clearance_height_m = {};
    std::array<double, W> theta0_r = {};
    std::array<double, W> theta_r = {};
    std::array<double, W> sint = {};
    std::array<double, W> cost = {};
    std::array<double, W> pin_phi_r = {};
    std::array<double, W> pin_x_m = {};
    std::array<double, W> pin_y_m = {};
    std::array<double, W> bearing_x_m = {};
    std::array<double, W> bearing_y_m = {};
    std::array<double, W> volumes_m3 = {};
    std::array<double, W> head_mass_kg = {};
    std::array<double, W> moment_of_inertia_kg_m2 = {};
    std::array<double, W> gas_torque_n_m = {};
    std::array<double, W> inertia_torque_n_m = {};
    std::array<double, W> friction_torque_n_m = {};
    std::array<double, W> total_torque_n_m = {};
    std::array<double, W> chamber_static_pressure_pa = {};
    std::array<double, W> friction_n_m_s2_per_r2 = {};
    double flywheel_angular_velocity_r_per_s = 0.0;
    double flywheel_theta_r = 0.0;

    /*
     * t = t0 + t1
     */

    fn void calc_thetas()
    {
        for(size_t i = 0; i < W; i++)
        {
            theta_r[i] = flywheel_theta_r - theta0_r[i];
        }
    }

    fn void calc_sin_cos()
    {
        for(size_t i = 0; i < W; i++)
        {
            const double t = theta_r[i];
            sint[i] = sin(t);
            cost[i] = cos(t);
        }
    }

    fn void calc_pin_phi()
    {
        for(size_t i = 0; i < W; i++)
        {
            const double r = crank_throw_length_m[i];
            const double l = connecting_rod_length_m[i];
            pin_phi_r[i] = -asin(sint[i] * r / l);
        }
    }

    /*
     * Hailemariam Nigus. Kinematics and Load Formulation of Engine Crank Mechanism. Mechanics, Materials/
     * Science & Engineering Journal, 2015, ⟨10.13140/RG.2.1.3257.1928⟩. ⟨hal-01305936⟩
     *
     *                        ________________
     *                       /
     *                      /  2    2    2
     * y = r * cos(t) + _  /  l  - r  sin (t)
     *                   \/
     */

    fn void calc_positions()
    {
        for(size_t i = 0; i < W; i++)
        {
            const double r = crank_throw_length_m[i];
            const double l = connecting_rod_length_m[i];
            const double x = r * sint[i];
            const double y = r * cost[i];
            bearing_x_m[i] = x;
            bearing_y_m[i] = y;
            pin_x_m[i] = 0.0;
            pin_y_m[i] = y + sqrt(l * l - x * x);
        }
    }

    /*
     *         2
     * v = pi r  h
     *
     */

    fn void calc_volumes()
    {
        for(size_t i = 0; i < W; i++)
        {
            const double r = crank_throw_length_m[i];
            const double l = connecting_rod_length_m[i];
            const double cm = head_compression_height_m[i];
            const double cl = head_clearance_height_m[i];
            const double block_deck_surface_m = r + l + cm + cl;
            const double y = pin_y_m[i] + cm;
            const double radius = diameter_m[i] / 2.0;
            const double h = block_deck_surface_m - y;
            volumes_m3[i] = g_pi_r * radius * radius * h;
        }
    }

    /*           2
     * M = pi * r  * h * p
     *
     */

    fn void calc_masses()
    {
        for(size_t i = 0; i < W; i++)
        {
            const double r = 0.5 * diameter_m[i];
            const double h = 2.0 * head_compression_height_m[i];
            const double p = head_mass_density_kg_per_m3[i];
            head_mass_kg[i] = g_pi_r * r * r * h * p;
        }
    }

    /*             1        2
     * I = [ mp + --- mr ] r
     *             3
     */

    fn void calc_moments_of_inertia()
    {
        for(size_t i = 0; i < W; i++)
        {
            const double r = crank_throw_length_m[i];
            const double mp = head_mass_kg[i];
            const double mr = connecting_rod_mass_kg[i];
            moment_of_inertia_kg_m2[i] = (mp + (1.0 / 3.0) * mr) * r * r;
        }
    }

    /*
     * Hailemariam Nigus. Kinematics and Load Formulation of Engine Crank Mechanism. Mechanics, Materials/
     * Science & Engineering Journal, 2015, ⟨10.13140/RG.2.1.3257.1928⟩. ⟨hal-01305936⟩
     *
     *                           r
     * Tg = Pg A r sin(t) [ 1 + --- cos(t) ]
     *                           l
     */

    fn void calc_gas_torques()
    {
        for(size_t i = 0; i < W; i++)
        {
            const double Pg = chamber_static_pressure_pa[i] - g_ambient_pressure_pa;
            const double ar = diameter_m[i] / 2.0;
            const double A = g_pi_r * ar * ar;
            const double r = crank_throw_length_m[i];
            const double l = connecting_rod_length_m[i];
            const double X = Pg * A * r * sint[i];
            const double Y = 1.0 + (r / l) * cost[i];
            gas_torque_n_m[i] = X * Y;
        }
    }

    /*
     * Hailemariam Nigus. Kinematics and Load Formulation of Engine Crank Mechanism. Mechanics, Materials/
     * Science & Engineering Journal, 2015, ⟨10.13140/RG.2.1.3257.1928⟩. ⟨hal-01305936⟩
     *
     *           2      r            1             3r
     * Ti = I * w * [ ---- sin(t) - --- sin(2t) - ---- * sin(3t) ]
     *                 4l            2             4l
     *
     * These identities free up the SIMD lanes:
     *
     *     sin(2t) = 2 sin(t) * 1 cos(t)
     *
     *                               3
     *     sin(3t) = 3 sin(t) − 4 sin (t)
     *
     */

    fn void calc_inertia_torques()
    {
        for(size_t i = 0; i < W; i++)
        {
            const double r = crank_throw_length_m[i];
            const double l = connecting_rod_length_m[i];
            const double I = moment_of_inertia_kg_m2[i];
            const double w = flywheel_angular_velocity_r_per_s;
            const double rl = r / l;
            const double s = sint[i];
            const double c = cost[i];
            const double X = 0.25 * rl * s;
            const double Y = s * c;
            const double Z = 0.75 * rl * (3.0 * s - 4.0 * s * s * s);
            inertia_torque_n_m[i] = I * w * w * (X - Y - Z);
        }
    }

    /*
     * Tf = -K w | w |
     *
     */

    fn void calc_friction_torque()
    {
        for(size_t i = 0; i < W; i++)
        {
            const double K = friction_n_m_s2_per_r2[i];
            const double w = flywheel_angular_velocity_r_per_s;
            friction_torque_n_m[i] = -K * w * fabs(w);
        }
    }

    /*
     * Tt = Tg + Ti + Tf
     *
     */

    fn void calc_total_torque()
    {
        for(size_t i = 0; i < W; i++)
        {
            const double Tg = gas_torque_n_m[i];
            const double Ti = inertia_torque_n_m[i];
            const double Tf = friction_torque_n_m[i];
            total_torque_n_m[i] = Tg + Ti + Tf;
        }
    }

    fn void calc_volumetrics()
    {
        calc_thetas();
        calc_sin_cos();
        calc_pin_phi();
        calc_positions();
        calc_volumes();
        calc_masses();
        calc_moments_of_inertia();
    }

    fn void update()
    {
        calc_volumetrics();
        calc_gas_torques();
        calc_inertia_torques();
        calc_friction_torque();
        calc_total_torque();
    }
};

template<size_t WW, size_t L, size_t S>
struct pipe
{
    static constexpr size_t M = L - 1;
    std::array<float, WW> piston_connect_m = {};
    float length_m = 0.0f;
    float mic_position0_m = 0.0f;
    float mic_position1_m = 0.0f;

    pipe()
    {
        reset();
    }

    std::array<float, WW> in_velocity_m_per_s = {};
    std::array<float, WW> in_static_density_kg_per_m3 = {};
    std::array<float, WW> in_static_temperature_k = {};

    /*
     * Cells. U is Conserved state: F is flux state: Ff is flux face state.
     *
     * +------+     +------+     +------+       +--------+
     * |  U0  |     |  U1  |     |  U2  | ..... |  UL-1  |
     * +------+     +------+     +------+       +--------+
     * +------+     +------+     +------+       +--------+
     * |  F0  |     |  F1  |     |  F2  | ..... |  FL-1  |
     * +------+     +------+     +------+       +--------+
     *        +-----+      +-----+      +-------+
     *        | Ff0 |      | Ff1 |      | FfM-1 | M = L-1
     *        +-----+      +-----+      +-------+
     *
     * |---------------------- L ------------------------|
     *
     */

    std::array<float, L> U_inject_r = {};
    std::array<float, L> U_inject_ru = {};
    std::array<float, L> U_inject_rEs = {};

    std::array<float, L> U_r = {};
    std::array<float, L> U_ru = {};
    std::array<float, L> U_rEs = {};
    std::array<float, L> F_r = {};
    std::array<float, L> F_ru = {};
    std::array<float, L> F_rEs = {};
    std::array<float, M> Ff_r = {};
    std::array<float, M> Ff_ru = {};
    std::array<float, M> Ff_rEs = {};

    std::array<float, L> speed_of_sound_m_per_s = {};
    std::array<float, L> absolute_speed_of_sound_m_per_s = {};
    std::array<float, L> static_pressure_pa = {};

    /*
     *                1   2    r Rs         1     2
     * rEs = Cv Ts + --- u  = ------- Ts + --- r u
     *                2        y - 1        2
     */

    fn float calc_rEs_from_Ts(const float r, const float u, const float Ts)
    {
        const float ru = r * u;
        const float Rs = g_specific_gas_constant_j_per_kg_k;
        const float rEs = r * Rs * Ts / (g_gamma_f - 1.0f) + ru * ru / (2.0f * r);
        return rEs;
    }

    /*
     *                     2
     *         Ps       r u
     * rEs = ------- + -----
     *        y - 1      2
     */

    fn float calc_rEs_from_Ps(const float r, const float u, const float Ps)
    {
        const float rEs = Ps / (g_gamma_f - 1.0f) + 0.5f * r * u * u;
        return rEs;
    }

    /*
     *                            1     2
     * Ps = (y - 1) * r * [ Es - --- * u ]
     *                            2
     */

    fn float calc_Ps_from_Es(const float r, const float u, const float Es)
    {
        const float Ps = (g_gamma_f - 1.0f) * r * (Es - 0.5f * u * u);
        return Ps;
    }

    fn void to_ambient(const size_t i)
    {
        const float r = g_ambient_density_kg_per_m3;
        const float u = 0.0f;
        const float Ts = g_ambient_temperature_k;
        U_r[i] = r;
        U_ru[i] = r * u;
        U_rEs[i] = calc_rEs_from_Ts(r, u, Ts);
    }

    fn void reset()
    {
        for(size_t i = 0; i < L; i++)
        {
            to_ambient(i);
        }
    }

    /*
     *                     Ghost
     *                     Patm
     * +---+     +-----+ +-----+
     * |   |     |     | |     |
     * | 0 | ... | L-2 | | L-1 |
     * |   |     |  Y  | |  Z  |
     * +---+     +-----+ +-----+
     */

    fn void calc_pipe_open_right()
    {
        const size_t Y = L - 2;
        const size_t Z = L - 1;
        const float u = U_ru[Y] / U_r[Y];
        const float r = u >= 0.0f ? U_r[Y] : g_ambient_density_kg_per_m3;
        U_r[Z] = r;
        U_ru[Z] = r * u;
        U_rEs[Z] = calc_rEs_from_Ps(r, u, g_ambient_pressure_pa);
    }

    fn float calc_audio_sample()
    {
        const size_t Z = L - 2;
        const size_t x = Z * mic_position0_m / length_m;
        const size_t y = Z * mic_position1_m / length_m;
        return 0.1f * static_pressure_pa[x] + static_pressure_pa[y];
    }

    fn void inject()
    {
        U_inject_r = {};
        U_inject_ru = {};
        U_inject_rEs = {};
        for(size_t i = 0; i < WW; i++)
        {
            const float r = in_static_density_kg_per_m3[i];
            const float u = in_velocity_m_per_s[i];
            const float ru = r * u;
            const float Ts = in_static_temperature_k[i];
            const float ratio = piston_connect_m[i] / length_m;
            const size_t j = ratio * (L - 1);
            U_inject_r[j] += r;
            U_inject_ru[j] += ru;
            U_inject_rEs[j] += calc_rEs_from_Ts(r, u, Ts);
        }
        for(size_t i = 0; i < WW; i++)
        {
            const float ratio = piston_connect_m[i] / length_m;
            const size_t j = ratio * (L - 1);
            U_r[j] = U_inject_r[j];
            U_ru[j] = U_inject_ru[j];
            U_rEs[j] = U_inject_rEs[j];
        }
    }

    /*
     *  F = [ rr \ ruu + Ps \ u(rEs + Ps) ]
     */

    fn void calc_fluxes()
    {
        for(size_t i = 0; i < L; i++)
        {
            const float r = U_r[i];
            const float ru = U_ru[i];
            const float rEs = U_rEs[i];
            const float u = ru / r;
            const float Es = rEs / r;
            const float Ps = calc_Ps_from_Es(r, u, Es);
            F_r[i] = ru;
            F_ru[i] = ru * u + Ps;
            F_rEs[i] = u * (rEs + Ps);
        }
    }

    /*
     *          dt
     * U = U - ---- * [ Ffr - Ffl ]
     *          dx
     */

    fn void calc_conserved()
    {
        const float dx_m = length_m / L;
        const float dt_s = g_dt_s / S;
        const float dt_dx = dt_s / dx_m;
        for(size_t i = 1; i < L - 1; i++)
        {
            const size_t j = i - 1;
            U_r  [i] -= dt_dx * (Ff_r  [i] - Ff_r  [j]);
            U_ru [i] -= dt_dx * (Ff_ru [i] - Ff_ru [j]);
            U_rEs[i] -= dt_dx * (Ff_rEs[i] - Ff_rEs[j]);
        }
    }

    /*
     *                           1   2
     * Ps = (y - 1) * p * [ E - --- u ]
     *                           2
     */

    fn float calc_static_pressure(const size_t i)
    {
        const float r = U_r[i];
        const float ru = U_ru[i];
        const float rEs = U_rEs[i];
        const float u = ru / r;
        const float Es = rEs / r;
        const float Ps = calc_Ps_from_Es(r, u, Es);
        return Ps;
    }

    fn void calc_static_pressures()
    {
        for(size_t i = 0; i < L; i++)
        {
            const float Ps = calc_static_pressure(i);
            static_pressure_pa[i] = Ps;
        }
    }

    /*
     *       1                 1
     * Ff = --- [ Fl + Fr ] - --- alpha * [ Ur - Ul ]
     *       2                 2
     */

    fn void calc_flux_faces()
    {
        for(size_t i = 0; i < M; i++)
        {
            const size_t j = i + 1;
            const float Al = absolute_speed_of_sound_m_per_s[i];
            const float Ar = absolute_speed_of_sound_m_per_s[j];
            const float alpha = fmax(Al, Ar);
            Ff_r  [i] = 0.5f * ((F_r  [i] + F_r  [j]) - alpha * (U_r  [j] - U_r  [i]));
            Ff_ru [i] = 0.5f * ((F_ru [i] + F_ru [j]) - alpha * (U_ru [j] - U_ru [i]));
            Ff_rEs[i] = 0.5f * ((F_rEs[i] + F_rEs[j]) - alpha * (U_rEs[j] - U_rEs[i]));
        }
    }

    fn void update()
    {
        inject();
        for(size_t i = 0; i < S; i++)
        {
            calc_static_pressures();
            calc_speed_of_sounds();
            calc_absolute_speed_of_sounds();
            calc_pipe_open_right();
            calc_fluxes();
            calc_flux_faces();
            calc_conserved();
        }
    }

    fn void calc_speed_of_sounds()
    {
        for(size_t i = 0; i < L; i++)
        {
            const float Ps = static_pressure_pa[i];
            const float C = sqrt(g_gamma_f * Ps / U_r[i]);
            speed_of_sound_m_per_s[i] = C;
        }
    }

    fn void calc_absolute_speed_of_sounds()
    {
        for(size_t i = 0; i < L; i++)
        {
            const float a = U_ru[i] / U_r[i];
            const float c = speed_of_sound_m_per_s[i];
            absolute_speed_of_sound_m_per_s[i] = fabs(a) + c;
        }
    }
};

#define FLUIDS(X)                       \
    X(chamber_volume_m3)                \
    X(chamber_nozzle_real_flow_area_m2) \
    X(chamber_static_pressure_pa)       \
    X(chamber_static_temperature_k)     \
    X(chamber_mass_kg)                  \
    X(nozzle_static_temperature_k)      \
    X(nozzle_static_density_kg_per_m3)  \
    X(nozzle_velocity_m_per_s)          \
    X(nozzle_mach)                      \

#define PISTONS(X) X(total_torque_n_m)

#define DIAGS(X) FLUIDS(X) PISTONS(X)

enum
{
    #define X(name) g_##name,
    DIAGS(X)
    #undef X
    g_diags_size,
};

constexpr std::array<std::string_view, g_diags_size> g_signal_names = {
    #define X(name) #name,
    DIAGS(X)
    #undef X
};

struct diags
{
    std::vector<std::vector<double>> front = {};
    std::vector<std::vector<double>> back = {};

    diags()
    {
        front.resize(g_diags_size);
        back.resize(g_diags_size);
    }
};

struct dc_filter
{
    float x_prev = 0.0f;
    float y_prev = 0.0f;
    float alpha = 0.0f;

    dc_filter()
    {
        set_cutoff_frequency(5.0f);
    }

    void set_cutoff_frequency(const float cutoff_freq_hz)
    {
        const float rc = 1.0f / (2.0f * g_pi_r_f * cutoff_freq_hz);
        alpha = rc / (rc + g_dt_s_f);
    }

    float filter(const float x)
    {
        const float y = alpha * (y_prev + x - x_prev);
        x_prev = x;
        y_prev = y;
        return y;
    }
};

struct gain_filter
{
    float ratio = 0.0f;

    float filter(const float x)
    {
        return x * ratio;
    }
};

struct clamp_filter
{
    float filter(const float x)
    {
        return clamper(x, -1.0f, 1.0f);
    }
};

struct fade_in_filter
{
    gain_filter gain;
    clamp_filter clamp;

    float time_s = 1.0;

    float filter(const float x)
    {
        const float step_s = g_dt_s_f / time_s;
        gain.ratio += step_s;
        gain.ratio = clamp.filter(gain.ratio);
        return x * gain.ratio;
    }
};

struct vec3
{
    double x = {};
    double y = {};
    double z = {};
};

struct room
{
    static constexpr float dc_cutoff_hz = 100.0f;
    vec3 dim = {};
    vec3 imp = {};
    vec3 mic = {};
    double beta = 0.0;
    int order = 0;
    size_t impulse_size = 0;
    std::vector<float> impulse = {};
    dc_filter filter = {};

    room()
    {
        filter.set_cutoff_frequency(dc_cutoff_hz);
    }

    void generate()
    {
        impulse.resize(impulse_size);
        std::fill(impulse.begin(), impulse.end(), 0.0f);
        const int bound = order / 2 + 1;
        for(int ix = -bound; ix <= bound; ix++)
        for(int iy = -bound; iy <= bound; iy++)
        for(int iz = -bound; iz <= bound; iz++)
        for(int tx = 0; tx < 2; tx++)
        for(int ty = 0; ty < 2; ty++)
        for(int tz = 0; tz < 2; tz++)
        {
            const int ox = 2 * abs(ix) + tx;
            const int oy = 2 * abs(iy) + ty;
            const int oz = 2 * abs(iz) + tz;
            const int n  = ox + oy + oz;
            if(n > order)
            {
                continue;
            }
            const double px = 2.0 * ix * dim.x + (tx ? -imp.x : imp.x);
            const double py = 2.0 * iy * dim.y + (ty ? -imp.y : imp.y);
            const double pz = 2.0 * iz * dim.z + (tz ? -imp.z : imp.z);
            const double dx = px - mic.x;
            const double dy = py - mic.y;
            const double dz = pz - mic.z;
            const double d = sqrt(dx*dx + dy*dy + dz*dz);
            if(d < 1e-6)
            {
                continue;
            }
            const size_t delay = llround(d / g_ambient_speed_of_sound_m_per_s * g_sample_rate_hz);
            if(delay < 0 || delay >= impulse_size)
            {
                continue;
            }
            impulse[delay] += pow(beta, n) / d;
        }
        float max = 0.0f;
        for(float s: impulse)
        {
            max = std::max(max, s);
        }
        for(float& s: impulse)
        {
            s /= max;
        }
        for(float& s: impulse)
        {
            s = filter.filter(s);
        }
    }
};

struct convolution_filter
{
    std::vector<float> impulse = {};
    std::vector<float> history = {};
    size_t head = 0;
    size_t size = 0;
    size_t mask = 0;

    bool is_power_of_two(const size_t n)
    {
        return (n & (n - 1)) == 0;
    }

    void reset(std::vector<float>&& other)
    {
        size = other.size();
        assert(is_power_of_two(size));
        mask = size - 1;
        history.resize(2 * size);
        impulse = std::move(other);
        head = 0;
    }

    /*
     *         N - 1
     *         _____
     *         \
     *  y[n] = /____ h[i] * x[n - i]
     *         i = 0
     */

    fn float filter(const float x)
    {
        history[head] = x;
        history[head + size] = x;
        float y = 0.0f;
        for(size_t i = 0; i < size; i++)
        {
            y += impulse[i] * history[i + head];
        }
        head = (head - 1) & mask;
        return y;
    }
};

template<size_t W, size_t H>
struct mailbox
{
    struct alignas(g_cacheline_size)
    {
        std::atomic<double> throttle_open_ratio = 0.0;
        std::atomic<size_t> log_x = -1;
        std::atomic<size_t> log_y = -1;
        std::atomic<bool> injector_enabled = true;
        std::atomic<bool> clutch_engaged = true;
        std::atomic<bool> load_brake_engaged = false;
        std::atomic<size_t> gear = 0;
    }
    in;

    struct alignas(g_cacheline_size)
    {
        std::atomic<size_t> swap_drops = 0;
        std::atomic<double> load_angular_velocity_r_per_s = 0.0;
        std::atomic<double> engine_angular_velocity_r_per_s = 0.0;
        std::atomic<double> limiter_angular_velocity_r_per_s = 0.0;
        std::atomic<double> load_max_angular_velocity_r_per_s = 0.0;
        std::atomic<double> flywheel_mass_kg = 0.0;
        std::atomic<double> flywheel_radius_m = 0.0;
        std::atomic<double> flywheel_theta_r = 0.0;
        std::atomic<double> load_mass_kg = 0.0;
        std::atomic<double> load_radius_m = 0.0;
        std::atomic<double> load_theta_r = 0.0;
        std::atomic<bool> clutch_engaged = false;
        std::atomic<bool> clutch_disengaged = false;
        std::array<std::array<std::atomic<double>, W>, H> port_flow_area_m2 = {};
        std::array<std::array<std::atomic<double>, W>, H> port_open_ratios = {};
        std::array<std::array<std::atomic<bool>, W>, H> panics = {};
        std::array<std::array<std::atomic<double>, W>, H> chamber_volume_m3 = {};
        std::array<std::atomic<double>, W> pistons_head_radius_m = {};
        std::array<std::atomic<double>, W> pistons_head_height_m = {};
        std::array<std::atomic<double>, W> pistons_pin_y_m = {};
        std::array<std::atomic<double>, W> pistons_connecting_rod_length_m = {};
        std::array<std::atomic<double>, W> pistons_crank_diameter_m = {};
        std::array<std::atomic<double>, W> pistons_pin_phi_r = {};
        std::array<std::atomic<double>, W> pistons_top_dead_center_m = {};
        std::array<std::atomic<double>, W> pistons_crank_theta_r = {};
    }
    out;
};

template<size_t W, size_t H>
struct logger
{
    size_t x = -1;
    size_t y = -1;

    bool can_log()
    {
        return x < W and y < H;
    }
};

/*
 * [ ]  ...  [ ] | <- Source
 * [ ]  ...  [ | | <- Intake
 * [ ]  ...  [ ] | <- Intake Manifold
 * [ ]  ...  [ ] | <- Intake Runner
 * [ ]  ...  [ ] H <- Piston (PISTON_Y)
 * [ ]  ...  [ ] | <- Chamber0 --+
 * [ ]  ...  [ ] | <- Chamber1   |---- These three form an audio sampling division.
 * [ ]  ...  [ ] | <- Chamber2 --+     Use AUDIO_Y to sample the best sounding one.
 * [ ]  ...  [ ] | <- Sink
 * +---- W ----+ +
 *
 */

template<
    size_t W,
    size_t H,
    size_t THROTTLE_Y,
    size_t PISTON_Y,
    size_t AUDIO_Y,
    size_t PIPE_CELLS,
    size_t PIPE_SUBSTEPS,
    size_t PIPE_COUNT,
    size_t THROTTLE_STEPS,
    size_t GEAR_COUNT,
    template<size_t> class PISTONS,
    template<size_t> class CAMS,
    template<size_t> class SPARKPLUGS>
struct implements_engine: engine
{
    static constexpr size_t WW = W / PIPE_COUNT;
    double lumped_parasitic_torque_n_m = {};
    struct PISTONS<W> pistons = {};
    struct CAMS<W> inlet_cam = {};
    struct CAMS<W> outlet_cam = {};
    struct SPARKPLUGS<W> sparkplugs = {};
    std::array<struct flow<H, PISTON_Y>, W> flows = {};
    struct injector injector = {};
    struct limiter limiter = {};
    struct throttle<THROTTLE_STEPS> throttle = {};
    struct flywheel flywheel = {};
    struct crankshaft crankshaft = {};
    struct dc_filter dc = {};
    struct gain_filter gain = {};
    struct clamp_filter clamp = {};
    struct fade_in_filter fade_in = {};
    struct room room = {};
    struct convolution_filter convolution = {};
    struct gearbox<GEAR_COUNT> gearbox = {};
    struct clutch clutch = {};
    struct load load = {};
    struct diags diags = {};
    struct logger<W, H> logger = {};
    std::array<struct pipe<WW, PIPE_CELLS, PIPE_SUBSTEPS>, PIPE_COUNT> pipes = {};
    std::vector<float> audio_signal = {};
    struct mailbox<W, H> mailbox = {};
    std::array<std::vector<float>, PIPE_COUNT> pipe_static_pressures_pa = {};
    std::mutex swap_mutex = {};

    void log_logger()
    {
        if(logger.can_log())
        {
            #define X(name) diags.back[g_##name].push_back(flows[logger.x].name[logger.y]);
            FLUIDS(X)
            #undef X
            if(logger.y == PISTON_Y)
            {
                #define X(name) diags.back[g_##name].push_back(pistons.name[logger.x]);
                PISTONS(X)
                #undef X
            }
        }
    }

    fn void broadcast_state()
    {
        inlet_cam.flywheel_theta_r = flywheel.theta_r;
        outlet_cam.flywheel_theta_r = flywheel.theta_r;
        pistons.flywheel_theta_r = flywheel.theta_r;
        sparkplugs.flywheel_theta_r = flywheel.theta_r;
        pistons.flywheel_angular_velocity_r_per_s = flywheel.angular_velocity_r_per_s;
        limiter.flywheel_angular_velocity_r_per_s = flywheel.angular_velocity_r_per_s;
        inlet_cam.flywheel_angular_velocity_r_per_s = flywheel.angular_velocity_r_per_s;
        outlet_cam.flywheel_angular_velocity_r_per_s = flywheel.angular_velocity_r_per_s;
        for(size_t x = 0; x < W; x++)
        {
            flows[x].chamber_nozzle_open_ratio[PISTON_Y - 1] = inlet_cam.open_ratio[x];
            flows[x].chamber_nozzle_open_ratio[PISTON_Y + 0] = outlet_cam.open_ratio[x];
            flows[x].chamber_nozzle_open_ratio[THROTTLE_Y] = throttle.mapping.at(throttle.open_ratio);
            flows[x].piston_injector_enabled = injector.enabled && not limiter.limiting;
            flows[x].piston_chamber_radius_m = pistons.diameter_m[x] / 2.0;
            flows[x].chamber_volume_m3[PISTON_Y] = pistons.volumes_m3[x];
            pistons.chamber_static_pressure_pa[x] = flows[x].chamber_static_pressure_pa[PISTON_Y];
        }
    }

    void log_volumes()
    {
        for(size_t x = 0; x < W; x++)
        {
            flows[x].chamber_prev_volume_m3 = flows[x].chamber_volume_m3;
        }
    }

    fn void reset_chambers()
    {
        for(size_t x = 0; x < W; x++)
        {
            flows[x].calc_chamber_ambients();
        }
    }

    void reset() override
    {
        pistons.calc_volumetrics();
        broadcast_state();
        log_volumes();
        reset_chambers();
        post_mail(0);
        room.generate();
        convolution.reset(std::move(room.impulse));
    }

    bool diags_swap()
    {
        if(logger.can_log())
        {
            if(swap_mutex.try_lock())
            {
                for(size_t i = 0; i < diags.front.size(); i++)
                {
                    std::swap(diags.front[i], diags.back[i]);
                }
                for(auto& line: diags.back)
                {
                    line.clear();
                }
                size_t i = 0;
                for(auto& pipe: pipes)
                {
                    const auto& signal = pipe.static_pressure_pa;
                    pipe_static_pressures_pa[i++] = { signal.begin(), signal.end() };
                }
                swap_mutex.unlock();
                return true;
            }
            else
            {
                for(auto& line: diags.back)
                {
                    line.clear();
                }
                return false;
            }
        }
        return false;
    }

    void update_pipe()
    {
        for(size_t x = 0; x < W; x++)
        {
            const size_t i = x / WW;
            const size_t j = x % WW;
            pipes[i].in_velocity_m_per_s[j] = flows[x].nozzle_velocity_m_per_s[AUDIO_Y];
            pipes[i].in_static_temperature_k[j] = flows[x].nozzle_static_temperature_k[AUDIO_Y];
            pipes[i].in_static_density_kg_per_m3[j] = flows[x].nozzle_static_density_kg_per_m3[AUDIO_Y];
        }
        for(auto& pipe: pipes)
        {
            pipe.update();
        }
    }

    void update_limiter()
    {
        limiter.update();
    }

    double calc_engine_moment_of_inertia()
    {
        double I = 0.0;
        for(size_t x = 0; x < W; x++)
        {
            I += pistons.moment_of_inertia_kg_m2[x];
        }
        I += crankshaft.moment_of_inertia_kg_m2;
        I += flywheel.moment_of_inertia_kg_m2;
        return I;
    }

    double calc_engine_torque()
    {
        double T = 0.0;
        for(size_t x = 0; x < W; x++)
        {
            T += pistons.total_torque_n_m[x];
        }
        T -= lumped_parasitic_torque_n_m;
        return T;
    }

    void update_drivetrain()
    {
        crankshaft.update();
        flywheel.update();
        load.update();
        const double N = gearbox.lookup_gear_ratio();
        if(N == 0.0)
        {
            clutch.engaged = false;
        }
        const double Il = calc_engine_moment_of_inertia();
        const double Tl = calc_engine_torque();
        const double wl = flywheel.angular_velocity_r_per_s;
        const double Ir = load.moment_of_inertia_kg_m2;
        const double Tr = load.total_torque_n_m;
        const double wr = load.angular_velocity_r_per_s;
        const double Tc = clutch.calc_torque(wl, N * wr);
        const double Tt = clutch.engaged ? Tc : 0.0;
        const double al = (Tl - Tt) / Il;
        const double ar = (Tr + N * Tt) / Ir;
        crankshaft.accelerate(al);
        flywheel.accelerate(al);
        load.accelerate(ar);
        flywheel.turn();
    }

    void update_cams()
    {
        inlet_cam.update();
        outlet_cam.update();
    }

    void update_sparkplugs()
    {
        sparkplugs.update();
    }

    void update_pistons()
    {
        pistons.update();
    }

    void update_ignition()
    {
        for(size_t x = 0; x < W; x++)
        {
            if(sparkplugs.rising_edge[x])
            {
                flows[x].ignite_piston_chamber();
            }
        }
    }

    void update_flows()
    {
        for(size_t x = 0; x < W; x++)
        {
            flows[x].update();
        }
    }

    float calc_audio_sample()
    {
        float x = 0.0f;
        for(auto& pipe: pipes)
        {
            x += pipe.calc_audio_sample();
        }
        x = dc.filter(x);
        x = convolution.filter(x);
        x = fade_in.filter(x);
        x = gain.filter(x);
        x = clamp.filter(x);
        return x;
    }

    void collect_mail()
    {
        throttle.open_ratio = mailbox.in.throttle_open_ratio;
        logger.x = mailbox.in.log_x;
        logger.y = mailbox.in.log_y;
        injector.enabled = mailbox.in.injector_enabled;
        gearbox.gear = mailbox.in.gear;
        clutch.engaged = mailbox.in.clutch_engaged;
        load.brake_engaged = mailbox.in.load_brake_engaged;
        audio_signal.clear();
    }

    void post_mail(const size_t swap_drops)
    {
        mailbox.out.limiter_angular_velocity_r_per_s = limiter.max_angular_velocity_r_per_s;
        mailbox.out.load_max_angular_velocity_r_per_s = load.max_angular_velocity_r_per_s;
        mailbox.out.load_angular_velocity_r_per_s = load.angular_velocity_r_per_s;
        mailbox.out.engine_angular_velocity_r_per_s = flywheel.angular_velocity_r_per_s;
        mailbox.out.flywheel_mass_kg = flywheel.mass_kg;
        mailbox.out.flywheel_radius_m = flywheel.radius_m;
        mailbox.out.flywheel_theta_r = flywheel.theta_r;
        mailbox.out.load_mass_kg = load.mass_kg;
        mailbox.out.load_radius_m = load.radius_m;
        mailbox.out.load_theta_r = load.theta_r;
        mailbox.out.clutch_engaged = clutch.engaged;
        mailbox.out.clutch_disengaged = not clutch.engaged;
        for(size_t y = 0; y < H; y++)
        for(size_t x = 0; x < W; x++)
        {
            mailbox.out.port_flow_area_m2[y][x] = flows[x].chamber_nozzle_flow_area_m2[y];
            mailbox.out.port_open_ratios[y][x] = flows[x].chamber_nozzle_open_ratio[y];
            mailbox.out.panics[y][x] = flows[x].panic[y];
            mailbox.out.chamber_volume_m3[y][x] = flows[x].chamber_volume_m3[y];
        }
        for(size_t x = 0; x < W; x++)
        {
            mailbox.out.pistons_head_radius_m[x] = pistons.diameter_m[x] / 2.0;
            mailbox.out.pistons_head_height_m[x] = 2.0 * pistons.head_compression_height_m[x];
            mailbox.out.pistons_pin_y_m[x] = pistons.pin_y_m[x];
            mailbox.out.pistons_connecting_rod_length_m[x] = pistons.connecting_rod_length_m[x];
            mailbox.out.pistons_crank_diameter_m[x] = 2.0 * pistons.crank_throw_length_m[x];
            mailbox.out.pistons_pin_phi_r[x] = pistons.pin_phi_r[x];
            mailbox.out.pistons_crank_theta_r[x] = pistons.theta_r[x];
            mailbox.out.pistons_top_dead_center_m[x] = pistons.connecting_rod_length_m[x] + pistons.crank_throw_length_m[x] + pistons.head_compression_height_m[x];
        }
        mailbox.out.swap_drops += swap_drops;
    }

    void run(size_t steps) override
    {
        size_t swap_drops = 0;
        collect_mail();
        while(steps--)
        {
            update_cams();
            if(flywheel.otto_cycled())
            {
                if(not diags_swap())
                {
                    swap_drops++;
                }
            }
            update_sparkplugs();
            update_ignition();
            update_flows();
            update_pistons();
            update_limiter();
            update_drivetrain();
            update_pipe();
            log_logger();
            log_volumes();
            broadcast_state();
            const float sample = calc_audio_sample();
            audio_signal.push_back(sample);
        }
        post_mail(swap_drops);
    }

    void                          set_logger                           (const size_t x, const size_t y)       override { mailbox.in.log_x = x; mailbox.in.log_y = y;                }
    void                          set_swap_lock_on                     ()                                     override { swap_mutex.lock();                                         }
    void                          set_swap_lock_off                    ()                                     override { swap_mutex.unlock();                                       }
    void                          set_throttle                         (const double open_ratio)              override { mailbox.in.throttle_open_ratio = open_ratio;               }
    void                          set_injection_on                     ()                                     override { mailbox.in.injector_enabled = true;                        }
    void                          set_injection_off                    ()                                     override { mailbox.in.injector_enabled = false;                       }
    void                          engage_clutch                        ()                                     override { mailbox.in.clutch_engaged = true;                          }
    void                          disengage_clutch                     ()                                     override { mailbox.in.clutch_engaged = false;                         }
    void                          engage_load_brake                    ()                                     override { mailbox.in.load_brake_engaged = true;                      }
    void                          disengage_load_brake                 ()                                     override { mailbox.in.load_brake_engaged = false;                     }
    void                          increment_gear                       ()                                     override { if(mailbox.in.gear < gearbox.last_gear) mailbox.in.gear++; }
    void                          decrement_gear                       ()                                     override { if(mailbox.in.gear != 0) mailbox.in.gear--;                }
    size_t                        get_width                            ()                               const override { return W;                                                  }
    size_t                        get_height                           ()                               const override { return H;                                                  }
    size_t                        get_source_y                         ()                               const override { return 0;                                                  }
    size_t                        get_sink_y                           ()                               const override { return H - 1;                                              }
    size_t                        get_pipe_count                       ()                               const override { return PIPE_COUNT;                                         }
    size_t                        get_pipe_cells                       ()                               const override { return PIPE_CELLS;                                         }
    float                         get_pipe_length_m                    (const size_t pipe)              const override { return pipes[pipe].length_m;                               }
    size_t                        get_piston_y                         ()                               const override { return PISTON_Y;                                           }
    size_t                        get_audio_y                          ()                               const override { return AUDIO_Y;                                            }
    size_t                        get_throttle_y                       ()                               const override { return THROTTLE_Y;                                         }
    size_t                        get_bytes                            ()                               const override { return sizeof *this;                                       }
    size_t                        get_signal_count                     ()                               const override { return std::size(g_signal_names);                          }
    size_t                        get_max_gear                         ()                               const override { return GEAR_COUNT - 1;                                     }
    size_t                        get_gear                             ()                               const override { return mailbox.in.gear;                                    }
    double                        get_load_angular_velocity_r_per_s    ()                               const override { return mailbox.out.load_angular_velocity_r_per_s;          }
    double                        get_engine_angular_velocity_r_per_s  ()                               const override { return mailbox.out.engine_angular_velocity_r_per_s;        }
    size_t                        get_swap_drops                       ()                               const override { return mailbox.out.swap_drops;                             }
    double                        get_port_flow_area_m2                (const size_t x, const size_t y) const override { return mailbox.out.port_flow_area_m2[y][x];                }
    double                        get_port_open_ratio                  (const size_t x, const size_t y) const override { return mailbox.out.port_open_ratios[y][x];                 }
    bool                          get_panic                            (const size_t x, const size_t y) const override { return mailbox.out.panics[y][x];                           }
    double                        get_chamber_volume_m3                (const size_t x, const size_t y) const override { return mailbox.out.chamber_volume_m3[y][x];                }
    double                        get_piston_head_radius_m             (const size_t x)                 const override { return mailbox.out.pistons_head_radius_m[x];               }
    double                        get_piston_head_height_m             (const size_t x)                 const override { return mailbox.out.pistons_head_height_m[x];               }
    double                        get_piston_pin_y_m                   (const size_t x)                 const override { return mailbox.out.pistons_pin_y_m[x];                     }
    double                        get_piston_connecting_rod_length_m   (const size_t x)                 const override { return mailbox.out.pistons_connecting_rod_length_m[x];     }
    double                        get_piston_crank_diameter_m          (const size_t x)                 const override { return mailbox.out.pistons_crank_diameter_m[x];            }
    double                        get_piston_top_dead_center           (const size_t x)                 const override { return mailbox.out.pistons_top_dead_center_m[x];           }
    double                        get_piston_pin_phi_r                 (const size_t x)                 const override { return mailbox.out.pistons_pin_phi_r[x];                   }
    double                        get_piston_crank_theta_r             (const size_t x)                 const override { return mailbox.out.pistons_crank_theta_r[x];               }
    double                        get_flywheel_mass_kg                 ()                               const override { return mailbox.out.flywheel_mass_kg;                       }
    double                        get_flywheel_radius_m                ()                               const override { return mailbox.out.flywheel_radius_m;                      }
    double                        get_flywheel_theta_r                 ()                               const override { return mailbox.out.flywheel_theta_r;                       }
    double                        get_load_mass_kg                     ()                               const override { return mailbox.out.load_mass_kg;                           }
    double                        get_load_radius_m                    ()                               const override { return mailbox.out.load_radius_m;                          }
    double                        get_load_theta_r                     ()                               const override { return mailbox.out.load_theta_r;                           }
    double                        get_load_max_angular_velocity_r_per_s()                               const override { return mailbox.out.load_max_angular_velocity_r_per_s;      }
    double                        get_limiter_angular_velocity_r_per_s ()                               const override { return mailbox.out.limiter_angular_velocity_r_per_s;       }
    bool                          get_clutch_engaged                   ()                               const override { return mailbox.out.clutch_engaged;                         }
    bool                          get_clutch_disengaged                ()                               const override { return mailbox.out.clutch_disengaged;                      }
    std::string_view              get_signal_name                      (const size_t index)             const override { return g_signal_names[index];                              }
    const std::vector<double>&    get_signal                           (const size_t index)             const override { return diags.front[index];                                 }
    const std::vector<double>&    get_static_temperature_signal_k      ()                               const override { return get_signal(g_chamber_static_temperature_k);         }
    const std::vector<double>&    get_static_pressure_signal_pa        ()                               const override { return get_signal(g_chamber_static_pressure_pa);           }
    const std::vector<double>&    get_volume_signal_m3                 ()                               const override { return get_signal(g_chamber_volume_m3);                    }
    const std::vector<float>&     get_audio_signal                     ()                               const override { return audio_signal;                                       }
    const std::vector<float>&     get_pipe_pressure_signal             (const size_t pipe)              const override { return pipe_static_pressures_pa[pipe];                     }
    const std::vector<float>&     get_impulse_signal                   ()                               const override { return convolution.impulse;                                }
    std::tuple<float,float,float> get_room_dimensions_m                ()                               const override { return { room.dim.x, room.dim.y, room.dim.z };             }
    std::tuple<float,float,float> get_mic_position                     ()                               const override { return { room.mic.x, room.mic.y, room.mic.z };             }
    std::tuple<float,float,float> get_impulse_position                 ()                               const override { return { room.imp.x, room.imp.y, room.imp.z };             }
};

struct trx450r: implements_engine<1, 9, 2, 4, 5, 200, 20, 1, 4, 6, inline_pistons, cams, sparkplugs>
{
    trx450r()
    {
        this->room.dim = { 9.19, 3.11, 9.08 };
        this->room.imp = { 3.85, 0.21, 3.85 };
        this->room.mic = { 3.01, 1.22, 3.12 };
        this->room.beta = 0.94;
        this->room.order = 256;
        this->room.impulse_size = 8192;
        this->gearbox.ratios = { 0.0, 2.3, 1.8, 1.4, 1.2, 0.95 };
        this->clutch.damping_coefficient_n_m_s = 0.6;
        this->lumped_parasitic_torque_n_m = 0.2;
        this->limiter.max_angular_velocity_r_per_s = 800.0;
        this->limiter.limit_time_s = 0.05;
        this->load.mass_kg = 100.0;
        this->load.radius_m = 0.15;
        this->load.rod_friction_n_m_s_per_r = 0.05;
        this->load.brake_friction_n_m_s_per_r = 1.0;
        this->load.max_angular_velocity_r_per_s = 900.0;
        this->load.angular_velocity_r_per_s = 200.0;
        this->flywheel.angular_velocity_r_per_s = 200.0;
        this->flywheel.mass_kg = 1.45;
        this->flywheel.radius_m = 0.079;
        this->crankshaft.mass_kg = 8.4;
        this->crankshaft.radius_m = 0.047;
        this->pistons.friction_n_m_s2_per_r2.fill(0.00004);
        this->pistons.diameter_m.fill(0.0960);
        this->pistons.crank_throw_length_m.fill(0.03105);
        this->pistons.connecting_rod_length_m.fill(0.1150);
        this->pistons.connecting_rod_mass_kg.fill(0.42);
        this->pistons.head_mass_density_kg_per_m3.fill(2700.0);
        this->pistons.head_compression_height_m.fill(0.01);
        this->pistons.head_clearance_height_m.fill(0.0056);
        this->inlet_cam.ramp_theta_r.fill(g_pi_r * 0.7);
        this->outlet_cam.ramp_theta_r.fill(g_pi_r * 0.6);
        this->inlet_cam.engage_theta_r[0] = g_otto_intake_cycle_r - 0.7;
        this->sparkplugs.engage_theta_r[0] = g_otto_combustion_cycle_r;
        this->outlet_cam.engage_theta_r[0] = g_otto_exhaust_cycle_r + 0.6;
        flows[0].chamber_nozzle_open_ratio.fill(1.0);
        flows[0].chamber_nozzle_flow_area_m2 = {
            0.00200,
            0.00150,
            0.00030,
            0.00081,
            0.00105,
            0.00150,
            0.00175,
            0.00180,
        };
        flows[0].chamber_volume_m3 = {
            9999.0,
            0.0003,
            0.0001,
            0.0002,
            0.0000,
            0.0001,
            0.0003,
            0.0004,
            9999.9,
        };
        this->throttle.mapping.table = { 0.005, 0.05, 0.30, 1.00 };
        pipes[0].piston_connect_m = { 0.0 };
        pipes[0].length_m = 0.48;
        pipes[0].mic_position0_m = 0.11;
        pipes[0].mic_position1_m = 0.33;
        this->dc.set_cutoff_frequency(5.0);
        this->gain.ratio = 0.00001;
    }
};

struct i3_1500cc: implements_engine<3, 9, 2, 4, 5, 256, 5, 1, 4, 7, inline_pistons, vtec_cams, sparkplugs>
{
    i3_1500cc()
    {
        this->room.dim = { 7.99, 2.81, 7.98 };
        this->room.imp = { 6.51, 0.21, 5.31 };
        this->room.mic = { 5.91, 0.75, 4.92 };
        this->room.beta = 0.95;
        this->room.order = 256;
        this->room.impulse_size = 16384;
        this->gearbox.ratios = { 0.0, 3.50, 2.45, 1.80, 1.35, 1.00, 0.75 };
        this->clutch.damping_coefficient_n_m_s = 2.4;
        this->lumped_parasitic_torque_n_m = 0.0;
        this->limiter.max_angular_velocity_r_per_s = 800.0;
        this->limiter.limit_time_s = 0.05;
        this->load.mass_kg = 400.0;
        this->load.radius_m = 0.24;
        this->load.rod_friction_n_m_s_per_r = 0.15;
        this->load.brake_friction_n_m_s_per_r = 3.2;
        this->load.max_angular_velocity_r_per_s = 700.0;
        this->flywheel.mass_kg = 4.5;
        this->flywheel.radius_m = 0.105;
        this->flywheel.angular_velocity_r_per_s = 90.0;
        this->load.angular_velocity_r_per_s = 100.0;
        this->crankshaft.mass_kg = 12.0;
        this->crankshaft.radius_m = 0.075;
        this->pistons.diameter_m.fill(0.0820);
        this->pistons.crank_throw_length_m.fill(0.0473);
        this->pistons.connecting_rod_length_m.fill(0.1482);
        this->pistons.connecting_rod_mass_kg.fill(0.300);
        this->pistons.head_clearance_height_m.fill(0.008);
        this->pistons.head_mass_density_kg_per_m3.fill(2700.0);
        this->pistons.friction_n_m_s2_per_r2.fill(0.00005);
        this->pistons.head_compression_height_m.fill(0.020);
        this->inlet_cam.ramp_theta_r.fill(g_pi_r * 0.95);
        this->inlet_cam.ramp_boost_ratio = 1.0;
        this->inlet_cam.gain_boost_ratio = 1.5;
        this->inlet_cam.engage_angular_velocity_r_per_s = 800.0;
        this->outlet_cam.ramp_theta_r.fill(g_pi_r * 0.90);
        this->outlet_cam.ramp_boost_ratio = 1.05;
        this->outlet_cam.gain_boost_ratio = 1.5;
        this->outlet_cam.engage_angular_velocity_r_per_s = 800.0;
        double theta0_r = 0.0;
        for(size_t i = 0; i < get_width(); i++)
        {
            this->pistons.theta0_r[i] = theta0_r;
            this->inlet_cam.engage_theta_r[i]  = theta0_r + g_otto_intake_cycle_r - 0.75;
            this->sparkplugs.engage_theta_r[i] = theta0_r + g_otto_combustion_cycle_r - 0.1;
            this->outlet_cam.engage_theta_r[i] = theta0_r + g_otto_exhaust_cycle_r - 0.48;
            theta0_r += g_otto_cycle_r / get_width();
        }
        for(auto& flow: this->flows)
        {
            flow.chamber_nozzle_open_ratio.fill(1.0);
            flow.chamber_nozzle_flow_area_m2 = {
                0.005,
                0.004,
                0.0005,
                0.003,
                0.0025,
                0.004,
                0.005,
                0.006,
            };
        }
        flows[0].chamber_volume_m3 = { 9999.0, 0.001, 0.001, 0.001, 0.000, 0.0010, 0.001, 0.001, 9999.0, };
        flows[1].chamber_volume_m3 = { 9999.0, 0.001, 0.001, 0.001, 0.000, 0.0012, 0.001, 0.001, 9999.0, };
        flows[2].chamber_volume_m3 = { 9999.0, 0.001, 0.001, 0.001, 0.000, 0.0013, 0.001, 0.001, 9999.0, };
        this->throttle.mapping.table = { 0.00001, 0.02, 0.1, 1.0 };
        pipes[0].piston_connect_m = { 0.00, 0.05, 0.02, };
        pipes[0].length_m = 1.7;
        pipes[0].mic_position0_m = 1.61;
        pipes[0].mic_position1_m = 1.7;
        this->dc.set_cutoff_frequency(10.0);
        this->gain.ratio = 0.00004;
    }
};

struct i8_2400cc: implements_engine<8, 9, 2, 4, 5, 200, 10, 2, 4, 9, inline_pistons, cams, sparkplugs>
{
    i8_2400cc()
    {
        this->room.dim = { 5.59, 2.51, 5.58 };
        this->room.imp = { 3.51, 0.12, 3.31 };
        this->room.mic = { 4.21, 0.74, 4.42 };
        this->room.beta = 0.95;
        this->room.order = 256;
        this->room.impulse_size = 16384;
        this->gearbox.ratios = { 0.0, 3.20, 2.35, 1.85, 1.50, 1.25, 1.08, 0.94, 0.84 };
        this->clutch.damping_coefficient_n_m_s = 5.0;
        this->lumped_parasitic_torque_n_m = 0.0;
        this->limiter.max_angular_velocity_r_per_s = 1900.0;
        this->limiter.limit_time_s = 0.025;
        this->load.mass_kg = 100.0;
        this->load.radius_m = 0.2;
        this->load.rod_friction_n_m_s_per_r = 0.1;
        this->load.brake_friction_n_m_s_per_r = 1.5;
        this->load.max_angular_velocity_r_per_s = 3000.0;
        this->flywheel.mass_kg = 10.5;
        this->flywheel.radius_m = 0.15;
        this->flywheel.angular_velocity_r_per_s = 100.0;
        this->load.angular_velocity_r_per_s = 100.0;
        this->crankshaft.mass_kg = 20.0;
        this->crankshaft.radius_m = 0.10;
        this->pistons.diameter_m.fill(0.098);
        this->pistons.crank_throw_length_m.fill(0.0199);
        this->pistons.connecting_rod_length_m.fill(0.105);
        this->pistons.connecting_rod_mass_kg.fill(0.30);
        this->pistons.head_mass_density_kg_per_m3.fill(7800.0);
        this->pistons.head_compression_height_m.fill(0.012);
        this->pistons.head_clearance_height_m.fill(0.004);
        this->pistons.friction_n_m_s2_per_r2.fill(0.000005);
        this->inlet_cam.ramp_theta_r.fill(g_pi_r * 1.0);
        this->outlet_cam.ramp_theta_r.fill(g_pi_r * 0.85);
        double theta0_r = 0.0;
        for(size_t i = 0; i < get_width(); i++)
        {
            this->pistons.theta0_r[i] = theta0_r;
            this->inlet_cam.engage_theta_r[i]  = theta0_r + g_otto_intake_cycle_r - 0.8;
            this->sparkplugs.engage_theta_r[i] = theta0_r + g_otto_combustion_cycle_r - 0.4;
            this->outlet_cam.engage_theta_r[i] = theta0_r + g_otto_exhaust_cycle_r - 0.5;
            theta0_r += g_otto_cycle_r / get_width();
        }
        for(auto& flow: this->flows)
        {
            flow.chamber_nozzle_open_ratio.fill(1.0);
            flow.chamber_nozzle_flow_area_m2 = {
                0.00350,
                0.00300,
                0.00100,
                0.00250,
                0.00190,
                0.00300,
                0.00350,
                0.00350,
            };
        }
        this->flows[0].chamber_volume_m3 = { 9999.0, 0.0005, 0.0005, 0.001105, 0.0000, 0.000665, 0.0005, 0.0005, 9999.0 };
        this->flows[1].chamber_volume_m3 = { 9999.0, 0.0005, 0.0005, 0.001042, 0.0000, 0.000728, 0.0005, 0.0005, 9999.0 };
        this->flows[2].chamber_volume_m3 = { 9999.0, 0.0005, 0.0005, 0.000980, 0.0000, 0.000891, 0.0005, 0.0005, 9999.0 };
        this->flows[3].chamber_volume_m3 = { 9999.0, 0.0005, 0.0005, 0.000917, 0.0000, 0.000954, 0.0005, 0.0005, 9999.0 };
        this->flows[4].chamber_volume_m3 = { 9999.0, 0.0005, 0.0005, 0.000854, 0.0000, 0.001017, 0.0005, 0.0005, 9999.0 };
        this->flows[5].chamber_volume_m3 = { 9999.0, 0.0005, 0.0005, 0.000791, 0.0000, 0.001180, 0.0005, 0.0005, 9999.0 };
        this->flows[6].chamber_volume_m3 = { 9999.0, 0.0005, 0.0005, 0.000628, 0.0000, 0.001242, 0.0005, 0.0005, 9999.0 };
        this->flows[7].chamber_volume_m3 = { 9999.0, 0.0005, 0.0005, 0.000565, 0.0000, 0.001305, 0.0005, 0.0005, 9999.0 };
        this->throttle.mapping.table = { 0.001, 0.01, 0.05, 0.20 };
        pipes[0].piston_connect_m = { 0.00, 0.21, 0.11, 0.28, };
        pipes[1].piston_connect_m = { 0.00, 0.23, 0.18, 0.33, };
        for(auto& pipe: this->pipes)
        {
            pipe.length_m = 0.69;
            pipe.mic_position0_m = 0.36;
            pipe.mic_position1_m = 0.69;
        }
        this->dc.set_cutoff_frequency(5.0);
        this->gain.ratio = 0.00005;
    }
};

std::unique_ptr<engine> new_engine(const type type)
{
    std::unique_ptr<engine> engine;
    switch(type)
    {
    default:
    case type::trx450r:
        engine = std::make_unique<ensim::trx450r>();
        break;
    case type::i3_1500cc:
        engine = std::make_unique<ensim::i3_1500cc>();
        break;
    case type::i8_2400cc:
        engine = std::make_unique<ensim::i8_2400cc>();
        break;
    }
    engine->reset();
    return engine;
}

}
