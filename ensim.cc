#include "ensim.hh"

#include <array>
#include <numbers>
#include <cmath>

#define fn __attribute__((used))

namespace ensim {

using std::sin;
using std::cos;
using std::fmax;
using std::fmin;
using std::log;
using std::sqrt;
using std::trunc;
using std::exp;

static constexpr size_t sample_rate_hz = 48000;
static constexpr real dt_s = 1.0_r / sample_rate_hz;
static constexpr real pi_r = std::numbers::pi_v<real>;
static constexpr real otto_cycle_r = 4.0_r * pi_r;
static constexpr real otto_intake_cycle_r = 0.0_r * pi_r;
static constexpr real otto_compression_cycle_r = 1.0_r * pi_r;
static constexpr real otto_combustion_cycle_r = 2.0_r * pi_r;
static constexpr real otto_exhaust_cycle_r = 3.0_r * pi_r;
static constexpr real resevoir_volume_m3 = 1e9_r;

fn constexpr real clamper(const real value, const real lower, const real upper)
{
    return fmax(fmin(value, upper), lower);
}

fn constexpr real modulos(const real value, const real by)
{
    return value - trunc(value / by) * by;
}

fn real cuberoot(const real x)
{
    return exp(log(x) / 3.0_r);
}

template<size_t N> using lane = std::array<real, N>;
template<size_t N> using mask = std::array<bool, N>;

template<size_t H, size_t PY>
requires(H % 2 == 1)
struct flow
{
    static constexpr size_t N = H - 1;
    static_assert(N % 2 == 0);

    static constexpr real gamma = 3.0_r / 2.0_r;
    static constexpr real universal_gas_constant_j_per_mol_k = 8.3144598_r;
    static constexpr real cv_j_per_mol_k = universal_gas_constant_j_per_mol_k / (gamma - 1.0_r);
    static constexpr real molar_mass_kg_per_mol = 0.023_r;
    static constexpr real cv_j_per_kg_k = cv_j_per_mol_k / molar_mass_kg_per_mol;
    static constexpr real specific_gas_constant_j_per_kg_k = universal_gas_constant_j_per_mol_k / molar_mass_kg_per_mol;
    static constexpr real ambient_temperature_k = 300.0_r;
    static constexpr real ambient_pressure_pa = 101325.0_r;
    static constexpr real energy_octane_j_per_kg = 47.9e6_r;
    static constexpr real stoich_air_fuel_ratio = 14.7_r;

    lane<H> chamber_prev_volume_m3;
    lane<H> chamber_volume_m3;
    lane<H> chamber_nozzle_flow_area_m2;
    lane<H> chamber_nozzle_real_flow_area_m2;
    lane<H> chamber_nozzle_open_ratio;
    lane<H> chamber_pressure_pa;
    lane<H> chamber_temperature_k;
    lane<H> chamber_mass_kg;
    lane<H> nozzle_mach;
    lane<H> nozzle_velocity_m_per_s;
    lane<H> nozzle_static_density_kg_per_m3;
    lane<H> nozzle_mass_flow_rate_kg_per_s;
    lane<H> parcel_mass_kg;
    lane<H> parcel_temperature_k;
    mask<H> panic;

    real piston_injection_enabled;
    real piston_chamber_flame_height_m;
    real piston_chamber_mass_burned_m3;
    real piston_chamber_radius_m;
    bool piston_chamber_on_fire;

    /*
     *     Ps * V
     * m = -------
     *     Rs * Ts
     */

    fn void calc_chamber_ambients()
    {
        for(size_t i = 0; i < H; i++)
        {
            chamber_temperature_k[i] = ambient_temperature_k;
            chamber_pressure_pa[i] = ambient_pressure_pa;
        }
        for(size_t i = 0; i < H; i++)
        {
            const real Ps = chamber_pressure_pa[i];
            const real V = chamber_volume_m3[i];
            const real Rs = specific_gas_constant_j_per_kg_k;
            const real Ts = chamber_temperature_k[i];
            chamber_mass_kg[i] = (Ps * V) / (Rs * Ts);
        }
    }

    /*
     *       m * Rs * T
     * Ps = ------------
     *           V
     */

    fn void calc_chamber_pressures()
    {
        for(size_t i = 0; i < N; i++)
        {
            const real m = chamber_mass_kg[i];
            const real Rs = specific_gas_constant_j_per_kg_k;
            const real Ts = chamber_temperature_k[i];
            const real V = chamber_volume_m3[i];
            chamber_pressure_pa[i] = m * Rs * Ts / V;
        }
    }

    /*
     * Assume upstream chamber total pressure (Pt) is equal to upstream
     * chamber static pressure. Downstream chamber static pressure (Ps)
     * remains as-is.
     *              ____________________
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
            const real X = (2.0_r / (gamma - 1.0_r));
            const real Pi = chamber_pressure_pa[i];
            const real Pj = chamber_pressure_pa[j];
            const real Pt = fmax(Pi, Pj);
            const real Ps = fmin(Pi, Pj);
            const real direction = Pi > Pj ? 1.0_r : -1.0_r;

            /*
             *      y - 1                        3
             *      ----- = 0.3333... where y = ---
             *        y                          2
             * Term
             */

            static_assert(gamma == 3.0_r / 2.0_r);
            const real Y = cuberoot(Pt / Ps);
            const real M = direction * sqrt(X * (Y - 1.0_r));
            nozzle_mach[i] = clamper(M, -1.0_r, 1.0_r);
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
            const real Rs = specific_gas_constant_j_per_kg_k;
            const real Tt = chamber_temperature_k[i];
            const real M = nozzle_mach[i];
            const real X = gamma * Rs * Tt;
            const real Y = 0.5_r * (gamma - 1.0_r) * M * M;
            const real u = M * sqrt(X / (1.0_r + Y));
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
            const real Pt = chamber_pressure_pa[i];
            const real Rs = specific_gas_constant_j_per_kg_k;
            const real Tt = chamber_temperature_k[i];
            const real M = nozzle_mach[i];
            const real pt = Pt / (Rs * Tt);

            /*
             *        1                  3
             *      ----- = 2 where y = ---
             *      y - 1                2
             * Term
             */

            static_assert(gamma == 3.0_r / 2.0_r);
            const real C = 1.0_r + 0.5_r * (gamma - 1.0_r) * M * M;
            nozzle_static_density_kg_per_m3[i] = pt / (C * C);
        }
    }

    fn void calc_nozzle_real_flow_areas()
    {
        for(size_t i = 0; i < N; i++)
        {
            chamber_nozzle_real_flow_area_m2[i] = chamber_nozzle_open_ratio[i] * chamber_nozzle_flow_area_m2[i];
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
            const real ps = nozzle_static_density_kg_per_m3[i];
            const real A = chamber_nozzle_real_flow_area_m2[i];
            const real u = nozzle_velocity_m_per_s[i];
            const real mdot = ps * A * u;
            nozzle_mass_flow_rate_kg_per_s[i] = mdot;
        }
    }

    /*
     *      .
     * mp = m dt
     * Tp = Ts,upstream
     *
     */

    fn void calc_nozzle_parcels()
    {
        for(size_t i = 0; i < N; i++)
        {
            const size_t j = i + 1;
            const real mdot = nozzle_mass_flow_rate_kg_per_s[i];
            const real dm = mdot * dt_s;
            parcel_mass_kg[i] = dm;
            const real Ti = chamber_temperature_k[i];
            const real Tj = chamber_temperature_k[j];
            parcel_temperature_k[i] = mdot > 0.0_r ? Ti : Tj;
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
            const real Ts1 = chamber_temperature_k[i];
            const real V1 = chamber_prev_volume_m3[i];
            const real V2 = chamber_volume_m3[i];
            const real dv = V1 / V2;

            /*
             *               1             3
             *      y - 1 = --- where y = ---
             * Term          2             2
             *
             */

            static_assert(gamma == 3.0_r / 2.0_r);
            chamber_temperature_k[i] = Ts1 * sqrt(dv);
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
            const real dm = parcel_mass_kg[i];
            const real m = chamber_mass_kg[j];
            const real Tp = parcel_temperature_k[i];
            const real Ts0 = chamber_temperature_k[j];
            const real Ts1 = (Ts0 * m + Tp * dm) / (m + dm);
            chamber_temperature_k[j] = dm > 0.0_r ? Ts1 : Ts0;
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
            const real mi = chamber_mass_kg[i];
            const real mj = chamber_mass_kg[j];
            const real dm = parcel_mass_kg[i];
            const real m0 = mi - dm;
            const real m1 = mj + dm;

            /*
             *               1             3
             *      y - 1 = --- where y = ---
             * Term          2             2
             *
             */

            static_assert(gamma == 3.0_r / 2.0_r);
            chamber_temperature_k[i] *= sqrt(m0 / mi);
            chamber_temperature_k[j] *= sqrt(m1 / mj);
            chamber_mass_kg[i] = m0;
            chamber_mass_kg[j] = m1;
        }
    }

    fn void ignite_piston_chamber()
    {
        piston_chamber_flame_height_m = 0.0_r;
        piston_chamber_mass_burned_m3 = 0.0_r;
        piston_chamber_on_fire = true;
    }

    fn void calc_combustion()
    {
#if 1
        if(piston_chamber_on_fire)
        {
            if(piston_injection_enabled)
            {
                const real M = chamber_mass_kg[PY];
                const real Mf = M / (1.0 + stoich_air_fuel_ratio);
                const real Q = Mf * energy_octane_j_per_kg;
                const real Cv = cv_j_per_kg_k;
                const real dTs = Q / (M * Cv);
                chamber_temperature_k[PY] += dTs;
            }
            else
            {
            }
            piston_chamber_on_fire = false;
        }
#else
        if(piston_chamber_on_fire)
        {
            const real M = chamber_mass_kg[PY];
            const real V = chamber_volume_m3[PY];

            /*
             *              Ts    2
             *           [------]
             *             Ts0
             * S =  0.4 ----------------
             *              Ps    0.125
             *           [------]
             *             Ps0
             */

            const real Ts = chamber_temperature_k[PY];
            const real Ts0 = ambient_temperature_k;
            const real Ps = chamber_pressure_pa[PY];
            const real Ps0 = ambient_pressure_pa;
            const real Tr = Ts / Ts0;
            const real S = 0.4_r * Tr * Tr / sqrt(sqrt(sqrt(Ps/ Ps0)));

            /*
             *                 2
             * Vb = dh * pi * r
             *
             */

            const real dh = S * dt_s;
            const real h1 = piston_chamber_flame_height_m;
            const real h2 = h1 + dh;
            const real r = piston_chamber_radius_m;
            const real Vb = (h2 - h1) * pi_r * r * r;

            const real p = M / V;
            const real dm = Vb * p;
            const real Mb = piston_chamber_mass_burned_m3 + dm;

            /*
             *             M Q
             * Ts = Ts + --------
             *            AFR Cv
             *
             */

            if(Mb / M < 1.0_r)
            {
                const real Q = energy_octane_j_per_kg;
                const real Cv = cv_j_per_kg_k;
                const real X = (dm * Q);
                const real Y = (1.0_r + stoich_air_fuel_ratio) * M * Cv;
                const real dTs = X / Y;
                chamber_temperature_k[PY] += dTs;
            }
            else
            {
                piston_chamber_on_fire = false;
            }

            piston_chamber_flame_height_m = h2;
            piston_chamber_mass_burned_m3 = Mb;
        }
#endif
    }

    fn void calc_panics()
    {
        for(size_t i = 0; i < N; i++)
        {
            panic[i] |= chamber_mass_kg[i] <= 0.0_r;
            panic[i] |= chamber_temperature_k[i] <= 0.0_r;
            panic[i] |= chamber_pressure_pa[i] <= 0.0_r;
        }
    }

    fn void update()
    {
        calc_chamber_pressures();
        calc_nozzle_machs();
        calc_nozzle_velocities();
        calc_nozzle_static_densities();
        calc_nozzle_real_flow_areas();
        calc_nozzle_mass_flow_rates();
        calc_nozzle_parcels();
        calc_compressions();
        calc_forward_energy_transfers();
        calc_mass_transfers();
        calc_combustion();
        calc_panics();
    }
};

struct crankshaft
{
    real mass_kg;
    real radius_m;
    real theta_r;
    real last_theta_r;
    real angular_velocity_r_per_s;
    real angular_acceleration_r_per_s2;
    real moment_of_inertia_kg_m2;

    /*
     * dw = a * dt
     *
     */

    fn void accelerate()
    {
        const real a = angular_acceleration_r_per_s2;
        angular_velocity_r_per_s += a * dt_s;
    }

    /*
     * dth = w * dt
     *
     */

    fn void turn()
    {
        last_theta_r = theta_r;
        const real w = angular_velocity_r_per_s;
        theta_r += w * dt_s;
    }

    fn bool otto_cycled()
    {
        const real t0 = modulos(last_theta_r, otto_cycle_r);
        const real t1 = modulos(theta_r, otto_cycle_r);
        return t0 > t1;
    }

    /*
     *      1     2
     * I = --- m r
     *      2
     *
     */

    fn void calc_moment_of_inertia()
    {
        moment_of_inertia_kg_m2 = 0.5_r * mass_kg * radius_m * radius_m;
    }

    fn bool update()
    {
        calc_moment_of_inertia();
        accelerate();
        turn();
        return otto_cycled();
    }
};

template<size_t W>
struct sparkplugs
{
    lane<W> engage_theta_r;
    mask<W> prev_fired;
    mask<W> fired;
    mask<W> rising_edge;

    static constexpr real fire_delay_theta_r = 1e-1_r;
    real crankshaft_theta_r;

    fn void calc_fired()
    {
        for(size_t i = 0; i < W; i++)
        {
            prev_fired[i] = fired[i];
            const real theta0_r = modulos(crankshaft_theta_r, otto_cycle_r);
            const real theta1_r = engage_theta_r[i] >= otto_cycle_r ? (engage_theta_r[i] - otto_cycle_r) : engage_theta_r[i];
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
struct simple_cam
{
    lane<W> engage_theta_r;
    lane<W> ramp_theta_r;
    lane<W> open_ratio;

    real crankshaft_theta_r;

    /*      4            1      2      3
     * r = t  [ 35 - 84 t + 70 t - 20 t ]
     *
     */

    fn void calc_open_ratios()
    {
        for(size_t i = 0; i < W; i++)
        {
            real theta0_r = modulos(engage_theta_r[i], otto_cycle_r);
            if(theta0_r < 0.0_r)
            {
                theta0_r += otto_cycle_r;
            }
            real theta1_r = modulos(crankshaft_theta_r, otto_cycle_r);
            if(theta1_r < theta0_r)
            {
                theta1_r += otto_cycle_r;
            }
            const real open_r = theta1_r - theta0_r;
            const real t = open_r / ramp_theta_r[i];
            const real a = t * t * t * t;
            const real b = t * a;
            const real c = t * b;
            const real d = t * c;
            const real A = 35.0_r * a;
            const real B = 84.0_r * b;
            const real C = 70.0_r * c;
            const real D = 20.0_r * d;
            const real R = clamper(A - B + C - D, 0.0_r, 1.0_r);
            open_ratio[i] = theta1_r < theta0_r ? 0.0_r : R;
        }
    }

    void update()
    {
        calc_open_ratios();
    }
};

struct flywheel
{
    real mass_kg;
    real radius_m;
    real moment_of_inertia_kg_m2;

    /*
     *      1     2
     * I = --- m r
     *      2
     *
     */

    fn void calc_moment_of_inertia()
    {
        moment_of_inertia_kg_m2 = 0.5_r * mass_kg * radius_m * radius_m;
    }

    void update()
    {
        calc_moment_of_inertia();
    }
};

template<size_t W>
struct inline_pistons
{
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

    lane<W> diameter_m;
    lane<W> crank_throw_length_m;
    lane<W> connecting_rod_length_m;
    lane<W> connecting_rod_mass_kg;
    lane<W> head_mass_density_kg_per_m3;
    lane<W> head_compression_height_m;
    lane<W> head_clearance_height_m;
    lane<W> theta0_r;
    lane<W> theta_r;
    lane<W> sint;
    lane<W> cost;
    lane<W> pin_x_m;
    lane<W> pin_y_m;
    lane<W> bearing_x_m;
    lane<W> bearing_y_m;
    lane<W> volumes_m3;
    lane<W> head_mass_kg;
    lane<W> moment_of_inertia_kg_m2;
    lane<W> gas_torque_n_m;
    lane<W> inertia_torque_n_m;
    lane<W> friction_torque_n_m;
    lane<W> total_torque_n_m;
    lane<W> chamber_pressure_pa;

    real crankshaft_angular_velocity_r_per_s;
    real crankshaft_theta_r;

    static constexpr real friction_n_m_s2_per_r2 = 0.00005_r;

    /*
     * t = t0 + t1
     */

    fn void calc_thetas()
    {
        for(size_t i = 0; i < W; i++)
        {
            theta_r[i] = crankshaft_theta_r - theta0_r[i];
        }
    }

    fn void calc_sin_cos()
    {
        for(size_t i = 0; i < W; i++)
        {
            const real t = theta_r[i];
            sint[i] = sin(t);
            cost[i] = cos(t);
        }
    }

    /*
     * Hailemariam Nigus. Kinematics and Load Formulation of Engine Crank Mechanism. Mechanics, Materials/
     * Science & Engineering Journal, 2015, ⟨10.13140/RG.2.1.3257.1928⟩. ⟨hal-01305936⟩
     *
     *                        ________________
     *                       /
     *                      /  2    2    2
     * y = r * cos(t) + _  /  l  + r  sin (t)
     *                   \/
     */

    fn void calc_positions()
    {
        for(size_t i = 0; i < W; i++)
        {
            const real r = crank_throw_length_m[i];
            const real l = connecting_rod_length_m[i];
            const real x = r * sint[i];
            const real y = r * cost[i];
            bearing_x_m[i] = x;
            bearing_y_m[i] = y;
            pin_x_m[i] = 0.0_r;
            pin_y_m[i] = y + sqrt(l * l + x * x);
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
            const real r = crank_throw_length_m[i];
            const real l = connecting_rod_length_m[i];
            const real cm = head_compression_height_m[i];
            const real cl = head_clearance_height_m[i];
            const real block_deck_surface_m = r + l + cm + cl;
            const real y = pin_y_m[i] + cm;
            const real radius = diameter_m[i] / 2.0_r;
            const real h = block_deck_surface_m - y;
            volumes_m3[i] = pi_r * radius * radius * h;
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
            const real r = 0.5_r * diameter_m[i];
            const real h = 2.0_r * head_compression_height_m[i];
            const real p = head_mass_density_kg_per_m3[i];
            head_mass_kg[i] = pi_r * r * r * h * p;
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
            const real r = crank_throw_length_m[i];
            const real mp = head_mass_kg[i];
            const real mr = connecting_rod_mass_kg[i];
            moment_of_inertia_kg_m2[i] = (mp + (1.0_r / 3.0_r) * mr) * r * r;
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
            const real Pg = chamber_pressure_pa[i];
            const real A = pi_r * diameter_m[i] * diameter_m[i];
            const real r = crank_throw_length_m[i];
            const real l = connecting_rod_length_m[i];
            const real X = Pg * A * r * sint[i];
            const real Y = 1.0_r + (r / l) * cost[i];
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
            const real r = crank_throw_length_m[i];
            const real l = connecting_rod_length_m[i];
            const real I = moment_of_inertia_kg_m2[i];
            const real w = crankshaft_angular_velocity_r_per_s;
            const real rl = r / l;
            const real s = sint[i];
            const real c = cost[i];
            const real X = 0.25_r * rl * s;
            const real Y = s * c;
            const real Z = 0.75_r * rl * (3.0_r * s - 4.0_r * s * s * s);
            inertia_torque_n_m[i] = I * w * w * (X - Y - Z);
        }
    }

    /*
     *           2
     * Tf = - K w
     *
     */

    fn void calc_friction_torque()
    {
        for(size_t i = 0; i < W; i++)
        {
            const real K = friction_n_m_s2_per_r2;
            const real w = crankshaft_angular_velocity_r_per_s;
            friction_torque_n_m[i] = -K * w * w;
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
            const real Tg = gas_torque_n_m[i];
            const real Ti = inertia_torque_n_m[i];
            const real Tf = friction_torque_n_m[i];
            total_torque_n_m[i] = Tg + Ti + Tf;
        }
    }

    fn void calc_volumetrics()
    {
        calc_thetas();
        calc_sin_cos();
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

#define FLUIDS(X) \
    X(chamber_volume_m3) \
    X(chamber_nozzle_real_flow_area_m2) \
    X(nozzle_mass_flow_rate_kg_per_s) \
    X(chamber_mass_kg) \
    X(chamber_pressure_pa) \
    X(chamber_temperature_k)

#define SPARKPLUGS(X) \
    X(fired)

#define PISTONS(X) \
    X(gas_torque_n_m) \
    X(inertia_torque_n_m)

#define DIAGS(X) FLUIDS(X) SPARKPLUGS(X) PISTONS(X)

enum
{
    #define X(name) g_##name,
    DIAGS(X)
    #undef X
    g_diags_size,
};

static constexpr std::array<std::string_view, g_diags_size> signal_names = {
    #define X(name) #name,
    DIAGS(X)
    #undef X
};

struct highpass_filter
{
    real x0 = 0.0_r;
    real y0 = 0.0_r;
    const real rc;
    const real alpha;

    highpass_filter(const real cutoff_freq_hz)
        : rc(1.0_r / (2.0_r * pi_r * cutoff_freq_hz))
        , alpha(rc / (rc + dt_s))
    {
    }

    real filter(const real x)
    {
        const real y = alpha * (y0 + x - x0);
        x0 = x;
        y0 = y;
        return y;
    }
};

struct dc_filter : highpass_filter
{
    dc_filter(): highpass_filter(5.0_r) {}
};

struct volume_filter
{
    const real max = 1000.0_r;
    real ratio = 1.0_r;

    real filter(const real x)
    {
        return ratio * x / max;
    }
};

struct throttle
{
    real open_ratio = 0.1_r;
};

struct injector
{
    bool enabled = true;
};

/*
 * [ ]  ...  [ ] | <- Source
 * [ ]  ...  [ | | <- Intake
 * [ ]  ...  [ ] | <- Intake Manifold
 * [ ]  ...  [ ] | <- Intake Runner
 * [ ]  ...  [ ] H <- Piston (PY)
 * [ ]  ...  [ ] | <- Exhaust Runner
 * [ ]  ...  [ ] | <- Exhaust Manifold
 * [ ]  ...  [ ] | <- Exhaust
 * [ ]  ...  [ ] | <- Sink
 * +---- W ----+ +
 *
 */

template<typename... T>
concept aggregate_of = (std::is_aggregate_v<T> and ...);

template<
    size_t W,
    size_t H,
    size_t TY,
    size_t PY,
    size_t AY,
    template<size_t> class P,
    template<size_t> class C,
    template<size_t> class S>
requires(aggregate_of<crankshaft, flywheel, P<W>, flow<H, PY>, C<W>, S<W>>)
struct as_engine : engine
{
    crankshaft crankshaft = {};
    flywheel flywheel = {};
    throttle throttle = {};
    injector injector = {};
    P<W> pistons = {};
    C<W> inlet_cam = {};
    C<W> outlet_cam = {};
    S<W> sparkplugs = {};
    std::array<flow<H, PY>, W> flows = {};
    grid front;
    grid back;
    dc_filter dc;
    volume_filter volume;
    line audio_signal;

    as_engine()
    {
        front.resize(g_diags_size);
        back.resize(g_diags_size);
    }

    void log_at(const size_t x, const size_t y)
    {
        if(x < W and y < H)
        {
            #define X(name) back[g_##name].push_back(flows[x].name[y]);
            FLUIDS(X)
            #undef X
            if(y == PY)
            {
                #define X(name) back[g_##name].push_back(sparkplugs.name[x]);
                SPARKPLUGS(X)
                #undef X
                #define X(name) back[g_##name].push_back(pistons.name[x]);
                PISTONS(X)
                #undef X
            }
        }
    }

    fn void broadcast()
    {
        /*
         * Crankshaft theta -> inlet/outlet cams + pistons + sparkplugs thetas.
         *
         */

        inlet_cam.crankshaft_theta_r = crankshaft.theta_r;
        outlet_cam.crankshaft_theta_r = crankshaft.theta_r;
        pistons.crankshaft_theta_r = crankshaft.theta_r;
        sparkplugs.crankshaft_theta_r = crankshaft.theta_r;
        pistons.crankshaft_angular_velocity_r_per_s = crankshaft.angular_velocity_r_per_s;

        /*
         * Cam open ratios -> chamber open ratios.
         *
         */

        for(size_t x = 0; x < W; x++)
        {
            flows[x].chamber_nozzle_open_ratio[PY - 1] = inlet_cam.open_ratio[x];
            flows[x].chamber_nozzle_open_ratio[PY + 0] = outlet_cam.open_ratio[x];
        }

        /*
         * Throttle open ratios -> chamber open ratios.
         *
         */

        for(size_t x = 0; x < W; x++)
        {
            flows[x].chamber_nozzle_open_ratio[TY] = throttle.open_ratio;
        }

        /*
         * Piston shapes <-> flow shapes.
         *
         */

        for(size_t x = 0; x < W; x++)
        {
            flows[x].piston_injection_enabled = injector.enabled;
            flows[x].piston_chamber_radius_m = pistons.diameter_m[x] / 2.0_r;
            flows[x].chamber_volume_m3[PY] = pistons.volumes_m3[x];
            pistons.chamber_pressure_pa[x] = flows[x].chamber_pressure_pa[PY];
        }

        /*
         * System torque + inertia to system angular acceleration.
         *
         *      t
         * a = ---
         *      I
         */

        real I = 0.0_r;
        for(size_t x = 0; x < W; x++)
        {
            I += pistons.moment_of_inertia_kg_m2[x];
        }
        I += flywheel.moment_of_inertia_kg_m2;
        I += crankshaft.moment_of_inertia_kg_m2;
        real t = 0.0_r;
        for(size_t x = 0; x < W; x++)
        {
            t += pistons.total_torque_n_m[x];
        }
        crankshaft.angular_acceleration_r_per_s2 = t / I;
    }

    void remember_volumes()
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
        crankshaft.update();
        flywheel.update();
        pistons.calc_volumetrics();
        broadcast();
        remember_volumes();
        reset_chambers();
    }

    void diags_swap()
    {
        for(size_t i = 0; i < front.size(); i++)
        {
            std::swap(front[i], back[i]);
        }
        for(auto& line : back)
        {
            line.clear();
        }
    }

    void sample_audio()
    {
        real x0 = 0.0_r;
        for(size_t x = 0; x < W; x++)
        {
            x0 += flows[x].chamber_pressure_pa[AY];
        }
        const real x1 = dc.filter(x0);
        const real x2 = volume.filter(x1);
        audio_signal.push_back(x2);
    }

    void run(const size_t steps, const size_t xx, const size_t yy) override
    {
        audio_signal.clear();
        audio_signal.reserve(steps);
        for(size_t i = 0; i < steps; i++)
        {
            flywheel.update();
            const bool otto_cycled = crankshaft.update();
            if(otto_cycled)
            {
                diags_swap();
            }
            inlet_cam.update();
            outlet_cam.update();
            sparkplugs.update();
            pistons.update();
            for(size_t x = 0; x < W; x++)
            {
                if(sparkplugs.rising_edge[x])
                {
                    flows[x].ignite_piston_chamber();
                }
            }
            for(size_t x = 0; x < W; x++)
            {
                flows[x].update();
            }
            log_at(xx, yy);
            remember_volumes();
            broadcast();
            sample_audio();
        }
    }

    std::string_view get_signal_name(const size_t index) const override
    {
        return signal_names[index];
    }

    const line& get_signal(const size_t index) const override
    {
        return front[index];
    }

    const line& get_temperature_signal_k() const override
    {
        return get_signal(g_chamber_temperature_k);
    }

    const line& get_pressure_signal_pa() const override
    {
        return get_signal(g_chamber_pressure_pa);
    }

    const line& get_volume_signal_m3() const override
    {
        return get_signal(g_chamber_volume_m3);
    }

    const real& get_crankshaft_angular_velocity_r_per_s() const override
    {
        return crankshaft.angular_velocity_r_per_s;
    }

    size_t get_width() const override
    {
        return W;
    }

    size_t get_height() const override
    {
        return H;
    }

    size_t bytes() const override
    {
        return sizeof *this;
    }

    size_t get_piston_y() const override
    {
        return PY;
    }

    size_t get_audio_y() const override
    {
        return AY;
    }

    size_t get_throttle_y() const override
    {
        return TY;
    }

    real get_port_open_ratio(const size_t x, const size_t y) const override
    {
        return flows[x].chamber_nozzle_open_ratio[y];
    }

    real get_panic_status(const size_t x, const size_t y) const override
    {
        return flows[x].panic[y];
    }

    const line& get_audio_signal() const override
    {
        return audio_signal;
    }

    void set_throttle_open_ratio(const real open_ratio) override
    {
        throttle.open_ratio = open_ratio;
    }

    void set_injection_on() override
    {
        injector.enabled = true;
    }

    void set_injection_off() override
    {
        injector.enabled = false;
    }
};

struct inline8 : as_engine<8, 9, 2, 4, 5, inline_pistons, simple_cam, sparkplugs>
{
    inline8()
    {
        this->flywheel.mass_kg = 17.0_r;
        this->flywheel.radius_m = 0.25_r;
        this->crankshaft.mass_kg = 25.3_r;
        this->crankshaft.radius_m = 0.031_r;
        this->crankshaft.angular_velocity_r_per_s = 500.0_r;
        this->pistons.diameter_m.fill(0.085_r);
        this->pistons.crank_throw_length_m.fill(0.038_r);
        this->pistons.connecting_rod_length_m.fill(0.1_r);
        this->pistons.connecting_rod_mass_kg.fill(0.4_r);
        this->pistons.head_mass_density_kg_per_m3.fill(7800.0_r);
        this->pistons.head_compression_height_m.fill(0.025_r);
        this->pistons.head_clearance_height_m.fill(0.007_r);
        this->inlet_cam.ramp_theta_r.fill(pi_r * 0.9_r);
        this->outlet_cam.ramp_theta_r.fill(pi_r * 0.9_r);
        real theta0_r = 0.0_r;
        const size_t width = get_width();
        for(size_t i = 0; i < width; i++)
        {
            this->pistons.theta0_r[i] = theta0_r;
            this->inlet_cam.engage_theta_r[i] = theta0_r + otto_intake_cycle_r - pi_r / 8.0_r;
            this->sparkplugs.engage_theta_r[i] = theta0_r + otto_combustion_cycle_r + pi_r / 8.0_r;
            this->outlet_cam.engage_theta_r[i] = theta0_r + otto_exhaust_cycle_r - pi_r / 8.0_r;
            theta0_r += otto_cycle_r / static_cast<real>(width);
        }
        for(auto& flow : this->flows)
        {
            flow.chamber_nozzle_open_ratio.fill(1.0_r);
            flow.chamber_volume_m3 = {
                resevoir_volume_m3,
                0.0015_r,
                0.0015_r,
                0.0015_r,

                0.000_r,

                0.005_r,
                0.010_r,
                0.010_r,
                resevoir_volume_m3,
            };
            flow.chamber_nozzle_flow_area_m2 = {
                0.001_r,
                0.001_r,
                0.001_r,

                0.001_r,
                0.001_r,

                0.001_r,
                0.001_r,
                0.001_r,
            };
        }
    }
};

std::unique_ptr<engine> new_engine(const type type)
{
    std::unique_ptr<engine> engine;
    if(type == type::inline8)
    {
        engine = std::make_unique<inline8>();
    }
    engine->reset();
    return engine;
}

} /* namespace ensim */
