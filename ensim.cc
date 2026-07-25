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
static constexpr real pi = std::numbers::pi_v<real>;
static constexpr real dt_s = 1.0_r / sample_rate_hz;
static constexpr real otto_cycle_r = 4.0_r * pi;
static constexpr real otto_intake_cycle_r = 0.0_r * pi;
static constexpr real otto_compression_cycle_r = 1.0_r * pi;
static constexpr real otto_combustion_cycle_r = 2.0_r * pi;
static constexpr real otto_exhaust_cycle_r = 3.0_r * pi;
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
    static constexpr real ambient_static_temperature_k = 300.0_r;
    static constexpr real ambient_static_pressure_pa = 101325.0_r;
    static constexpr real energy_octane_j_per_kg = 47.9e6_r;
    static constexpr real stoich_air_fuel_ratio = 14.7_r;

    lane<H> chamber_prev_volume_m3;
    lane<H> chamber_volume_m3;
    lane<H> chamber_nozzle_flow_area_m2;
    lane<H> chamber_nozzle_open_ratio;
    lane<H> chamber_static_pressure_pa;
    lane<H> chamber_static_temperature_k;
    lane<H> chamber_mass_kg;
    lane<H> nozzle_mach;
    lane<H> nozzle_velocity_m_per_s;
    lane<H> nozzle_static_density_kg_per_m3;
    lane<H> nozzle_mass_flow_rate_kg_per_s;
    lane<H> parcel_mass_kg;
    lane<H> parcel_static_temperature_k;
    mask<H> panic;

    /*
     *     Ps * V
     * m = -------
     *     Rs * Ts
     */

    fn void calc_chamber_ambients()
    {
        for(size_t i = 0; i < H; i++)
        {
            chamber_static_temperature_k[i] = ambient_static_temperature_k;
            chamber_static_pressure_pa[i] = ambient_static_pressure_pa;
        }
        for(size_t i = 0; i < H; i++)
        {
            const real Ps = chamber_static_pressure_pa[i];
            const real V = chamber_volume_m3[i];
            const real Rs = specific_gas_constant_j_per_kg_k;
            const real Ts = chamber_static_temperature_k[i];
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
            const real m = chamber_mass_kg[i];
            const real Rs = specific_gas_constant_j_per_kg_k;
            const real Ts = chamber_static_temperature_k[i];
            const real V = chamber_volume_m3[i];
            chamber_static_pressure_pa[i] = m * Rs * Ts / V;
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
            const real Pi = chamber_static_pressure_pa[i];
            const real Pj = chamber_static_pressure_pa[j];
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
            const real Tt = chamber_static_temperature_k[i];
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
            const real Pt = chamber_static_pressure_pa[i];
            const real Rs = specific_gas_constant_j_per_kg_k;
            const real Tt = chamber_static_temperature_k[i];
            const real M = nozzle_mach[i];
            const real X = Pt / (Rs * Tt);

            /*
             *        1                  3
             *      ----- = 2 where y = ---
             *      y - 1                2
             * Term
             */

            static_assert(gamma == 3.0_r / 2.0_r);
            const real C = 1.0_r + 0.5_r * (gamma - 1.0_r) * M * M;
            nozzle_static_density_kg_per_m3[i] = X / (C * C);
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
            const real A = chamber_nozzle_open_ratio[i] * chamber_nozzle_flow_area_m2[i];
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
            const real Ti = chamber_static_temperature_k[i];
            const real Tj = chamber_static_temperature_k[j];
            parcel_static_temperature_k[i] = mdot > 0.0_r ? Ti : Tj;
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
            const real Ts1 = chamber_static_temperature_k[i];
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
            const real dm = parcel_mass_kg[i];
            const real m = chamber_mass_kg[j];
            const real Tp = parcel_static_temperature_k[i];
            const real Ts0 = chamber_static_temperature_k[j];
            const real Ts1 = (Ts0 * m + Tp * dm) / (m + dm);
            chamber_static_temperature_k[j] = dm > 0.0_r ? Ts1 : Ts0;
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
            chamber_static_temperature_k[i] *= sqrt(m0 / mi);
            chamber_static_temperature_k[j] *= sqrt(m1 / mj);
            chamber_mass_kg[i] = m0;
            chamber_mass_kg[j] = m1;
        }
    }

    fn void calc_panics()
    {
        for(size_t i = 0; i < N; i++)
        {
            panic[i] |= chamber_mass_kg[i] <= 0.0_r;
            panic[i] |= chamber_static_temperature_k[i] <= 0.0_r;
            panic[i] |= chamber_static_pressure_pa[i] <= 0.0_r;
        }
    }

    /*
     *             Q
     * Ts = Ts + -------
     *           AFR Cv
     *
     */

    void calc_combustion()
    {
        const real Q = energy_octane_j_per_kg;
        const real Cv = cv_j_per_kg_k;
        const real dTs = Q / (stoich_air_fuel_ratio * Cv);
        chamber_static_temperature_k[PY] += dTs;
    }

    fn void update()
    {
        /*
         * Calculate conservation of momentum.
         *
         */

        calc_chamber_static_pressures();
        calc_nozzle_machs();
        calc_nozzle_velocities();
        calc_nozzle_static_densities();
        calc_nozzle_mass_flow_rates();

        /*
         * Calculates transport.
         *
         */

        calc_nozzle_parcels();

        /*
         * Calculate conservation of energy.
         *
         */

        calc_compressions();
        calc_forward_energy_transfers();

        /*
         * Calculate conservation of mass.
         *
         */

        calc_mass_transfers();

        /*
         * Misc debug.
         *
         */

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
    real moment_of_inertia_kg_m2;

    /*
     * dw = dw/dt * dt
     *
     */

    fn void accelerate(const real angular_acceleration_r_per_s2)
    {
        angular_velocity_r_per_s += angular_acceleration_r_per_s2 * dt_s;
    }

    /*
     * dth = dth/dt * dt
     *
     */

    fn void turn()
    {
        last_theta_r = theta_r;
        theta_r += angular_velocity_r_per_s * dt_s;
    }

    fn bool otto_cycled()
    {
        const real t0 = modulos(last_theta_r, otto_cycle_r);
        const real t1 = modulos(theta_r, otto_cycle_r);
        return t0 > t1;
    }

    /*
     *      1         2
     * I = --- * m * r
     *      2
     *
     */

    fn void calc_moment_of_inertia()
    {
        moment_of_inertia_kg_m2 = 0.5_r * mass_kg * radius_m * radius_m;
    }

    fn bool update(const real angular_acceleration_r_per_s2)
    {
        angular_velocity_r_per_s = 300.0_r;
        accelerate(angular_acceleration_r_per_s2);
        turn();
        const bool cycled = otto_cycled();
        calc_moment_of_inertia();
        return cycled;
    }
};

template<size_t W>
struct sparkplugs
{
    static constexpr real fire_delay_theta_r = 1e-1_r;
    lane<W> engage_theta_r;
    mask<W> prev_fired;
    mask<W> fired;
    mask<W> rising_edge;

    real crankshaft_theta_r;

    fn void calc_fired()
    {
        for(size_t i = 0; i < W; i++)
        {
            prev_fired[i] = fired[i];
            const real mod_theta_r = modulos(crankshaft_theta_r, otto_cycle_r);
            const real mod_engage_theta_r = engage_theta_r[i] >= otto_cycle_r
                ? (engage_theta_r[i] - otto_cycle_r)
                : (engage_theta_r[i]);
            const real mod_engage_delayed_theta_r = mod_engage_theta_r + fire_delay_theta_r;
            fired[i] = mod_theta_r > mod_engage_delayed_theta_r;
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
            real mod_engage_theta_r = modulos(engage_theta_r[i], otto_cycle_r);
            if(mod_engage_theta_r < 0.0_r)
            {
                mod_engage_theta_r += otto_cycle_r;
            }
            real mod_theta_r = modulos(crankshaft_theta_r, otto_cycle_r);
            if(mod_theta_r < mod_engage_theta_r)
            {
                mod_theta_r += otto_cycle_r;
            }
            const real open_r = mod_theta_r - mod_engage_theta_r;
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
            open_ratio[i] = mod_theta_r < mod_engage_theta_r ? 0.0_r : R;
        }
    }

    void update()
    {
        calc_open_ratios();
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
    lane<W> chamber_static_pressure_pa;

    real crankshaft_angular_velocity_r_per_s;
    real crankshaft_theta_r;

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
            volumes_m3[i] = pi * radius * radius * h;
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
            head_mass_kg[i] = pi * r * r * h * p;
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
            const real Pg = chamber_static_pressure_pa[i];
            const real A = pi * diameter_m[i] * diameter_m[i];
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

    fn void calc_volumetrics()
    {
        calc_thetas();
        calc_sin_cos();
        calc_positions();
        calc_volumes();
    }

    fn void update()
    {
        calc_volumetrics();
        calc_masses();
        calc_moments_of_inertia();
        calc_gas_torques();
        calc_inertia_torques();
    }
};

#define FLUIDS(X)                     \
    X(chamber_volume_m3)              \
    X(chamber_nozzle_open_ratio)      \
    X(nozzle_mass_flow_rate_kg_per_s) \
    X(chamber_mass_kg)                \
    X(chamber_static_pressure_pa)     \
    X(chamber_static_temperature_k)

#define SPARKPLUGS(X) \
    X(fired)

#define PISTONS(X)    \
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

/*
 * Engines are modelled by W (Width) indepdent SIMD flow lanes
 * from source to sink of H (Height) volumes. The piston row
 * is specified by PY (Piston Y).
 *
 * For example, an inline-3 configuration has
 * W == 3, H == 9, and PY == 4.
 *
 * [ ]  [ ]  [ ] | <- Source
 * [ ]  [ ]  [ | | <- Air Filter
 * [ ]  [ ]  [ ] | <- Throttle Body
 * [ ]  [ ]  [ ] | <- Intake Runner
 * [ ]  [ ]  [ ] H <- Piston (PY)
 * [ ]  [ ]  [ ] | <- Exhaust Runner
 * [ ]  [ ]  [ ] | <- Collector
 * [ ]  [ ]  [ ] | <- Exhaust
 * [ ]  [ ]  [ ] | <- Sink
 * +---- W ----+ +
 *
 */

template<typename... T>
concept aggregate_of = (std::is_aggregate_v<T> and ...);

template<
    size_t W,
    size_t H,
    size_t PY,
    template<size_t> class P,
    template<size_t> class C,
    template<size_t> class S>
requires(H == 9 and aggregate_of<crankshaft, P<W>, flow<H, PY>, C<W>, S<W>>)
struct as_engine : engine
{
    crankshaft crankshaft = {};
    P<W> pistons = {};
    C<W> inlet_cam = {};
    C<W> outlet_cam = {};
    S<W> sparkplugs = {};
    std::array<flow<H, PY>, W> flows = {};
    grid front = grid(g_diags_size);
    grid back = grid(g_diags_size);

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

    fn void broadcast_states()
    {
        inlet_cam.crankshaft_theta_r = crankshaft.theta_r;
        outlet_cam.crankshaft_theta_r = crankshaft.theta_r;
        pistons.crankshaft_theta_r = crankshaft.theta_r;
        sparkplugs.crankshaft_theta_r = crankshaft.theta_r;
        pistons.crankshaft_angular_velocity_r_per_s = crankshaft.angular_velocity_r_per_s;
        for(size_t x = 0; x < W; x++)
        {
            flows[x].chamber_volume_m3[PY] = pistons.volumes_m3[x];
            flows[x].chamber_nozzle_open_ratio[PY - 1] = inlet_cam.open_ratio[x];
            flows[x].chamber_nozzle_open_ratio[PY + 0] = outlet_cam.open_ratio[x];
        }
        for(size_t x = 0; x < W; x++)
        {
            pistons.chamber_static_pressure_pa[x] = flows[x].chamber_static_pressure_pa[PY];
        }
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
        pistons.calc_volumetrics();
        broadcast_states();
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

    void update_crankshaft()
    {
        const bool otto_cycled = crankshaft.update(0.0_r);
        if(otto_cycled)
        {
            diags_swap();
        }
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

    void update_flows()
    {
        for(size_t x = 0; x < W; x++)
        {
            flows[x].update();
        }
    }

    fn void update_injection()
    {
        for(size_t x = 0; x < W; x++)
        {
            if(sparkplugs.rising_edge[x])
            {
                flows[x].calc_combustion();
            }
        }
    }

    grid new_grid(const size_t width, const size_t height) const
    {
        grid grid(height);
        for(auto& row : grid)
        {
            row.resize(width);
        }
        return grid;
    }

    void run(const size_t steps, const size_t x, const size_t y) override
    {
        for(size_t i = 0; i < steps; i++)
        {
            update_crankshaft();
            update_cams();
            update_sparkplugs();
            update_pistons();
            update_injection();
            update_flows();
            log_at(x, y);
            remember_volumes();
            broadcast_states();
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

    const line& get_static_temperature_signal_k() const override
    {
        return get_signal(g_chamber_static_temperature_k);
    }

    const line& get_static_pressure_signal_pa() const override
    {
        return get_signal(g_chamber_static_pressure_pa);
    }

    const line& get_volume_signal_m3() const override
    {
        return get_signal(g_chamber_volume_m3);
    }

    size_t get_w() const override { return W; }
    size_t get_h() const override { return H; }
    size_t get_y() const override { return PY; }
    size_t bytes() const override { return sizeof *this; }

    real get_port_open_ratio(const size_t x, const size_t y) const override
    {
        return flows[x].chamber_nozzle_open_ratio[y];
    }

    real get_panic_status(const size_t x, const size_t y) const override
    {
        return flows[x].panic[y];
    }
};

struct inline8 : as_engine<8, 9, 4, inline_pistons, simple_cam, sparkplugs>
{
    inline8()
    {
        setup_crankshaft();
        setup_pistons();
        setup_flow();
    }

    void setup_crankshaft()
    {
        this->crankshaft.mass_kg = 25.3_r;
        this->crankshaft.radius_m = 0.031_r;
    }

    void setup_pistons()
    {
        this->pistons.diameter_m.fill(0.085_r);
        this->pistons.crank_throw_length_m.fill(0.038_r);
        this->pistons.connecting_rod_length_m.fill(0.1_r);
        this->pistons.connecting_rod_mass_kg.fill(0.4_r);
        this->pistons.head_mass_density_kg_per_m3.fill(7800.0_r);
        this->pistons.head_compression_height_m.fill(0.025_r);
        this->pistons.head_clearance_height_m.fill(0.007_r);
        this->inlet_cam.ramp_theta_r.fill(pi * 0.9_r);
        this->outlet_cam.ramp_theta_r.fill(pi);
        real theta0_r = 0.0_r;
        const size_t width = get_w();
        for(size_t i = 0; i < width; i++)
        {
            const real shift = pi / 8.0_r;
            this->pistons.theta0_r[i] = theta0_r;
            this->inlet_cam.engage_theta_r[i] = theta0_r - shift + otto_intake_cycle_r;
            this->sparkplugs.engage_theta_r[i] = theta0_r + shift + otto_combustion_cycle_r;
            this->outlet_cam.engage_theta_r[i] = theta0_r - shift + otto_exhaust_cycle_r;
            theta0_r += otto_cycle_r / static_cast<real>(width);
        }
    }

    void setup_flow()
    {
        for(auto& flow : this->flows)
        {
            flow.chamber_nozzle_open_ratio.fill(1.0_r);
            flow.chamber_volume_m3 = {
                resevoir_volume_m3,
                0.3_r,
                0.2_r,
                0.1_r,
                0.0_r,
                0.1_r,
                0.2_r,
                0.3_r,
                resevoir_volume_m3,
            };
            flow.chamber_nozzle_flow_area_m2 = {
                6e-3_r,
                5e-3_r,
                4e-3_r,
                3e-3_r,
                2e-3_r,
                2e-3_r,
                2e-3_r,
                2e-3_r,
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
