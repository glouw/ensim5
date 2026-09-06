#include "ensim.hh"

#include <array>
#include <numbers>
#include <cmath>
#include <mutex>
#include <cassert>

#define fn __attribute__((used))

namespace ensim
{
    template<size_t N> using lane = std::array<real, N>;
    template<size_t N> using mask = std::array<bool, N>;

    static constexpr real g_dt_s = 1.0_r / g_sample_rate_hz;
    static constexpr real g_pi_r = std::numbers::pi_v<real>;
    static constexpr real g_otto_cycle_r = 4.0_r * g_pi_r;
    static constexpr real g_otto_intake_cycle_r = 0.0_r * g_pi_r;
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

    using std::sin;
    using std::cos;
    using std::fmax;
    using std::fmin;
    using std::log;
    using std::sqrt;
    using std::trunc;
    using std::exp;

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

    fn real frand()
    {
        const real random = 2.0_r * rand() / static_cast<real>(RAND_MAX);
        return random - 1.0_r;
    }

    template<size_t H, size_t PY>
    struct flow
    {
        static constexpr size_t N = H - 1;
        static_assert(N % 2 == 0);

        lane<H> chamber_prev_volume_m3 = {};
        lane<H> chamber_volume_m3 = {};
        lane<H> chamber_nozzle_flow_area_m2 = {};
        lane<H> chamber_nozzle_real_flow_area_m2 = {};
        lane<H> chamber_nozzle_open_ratio = {};
        lane<H> chamber_static_pressure_pa = {};
        lane<H> chamber_dynamic_pressure_pa = {};
        lane<H> chamber_total_pressure_pa = {};
        lane<H> chamber_static_temperature_k = {};
        lane<H> chamber_dynamic_temperature_k = {};
        lane<H> chamber_total_temperature_k = {};
        lane<H> chamber_mass_kg = {};
        lane<H> chamber_bulk_momentum_kg_m_per_s = {};
        lane<H> nozzle_mach = {};
        lane<H> nozzle_velocity_m_per_s = {};
        lane<H> nozzle_static_density_kg_per_m3 = {};
        lane<H> nozzle_static_temperature_k = {};
        lane<H> nozzle_mass_flow_rate_kg_per_s = {};
        lane<H> parcel_mass_kg = {};
        lane<H> parcel_static_temperature_k = {};
        mask<H> panic = {};
        real piston_injection_enabled = 0.0_r;
        real piston_chamber_flame_height_m = 0.0_r;
        real piston_chamber_mass_burned_m3 = 0.0_r;
        real piston_chamber_radius_m = 0.0_r;
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
                const real Ps = chamber_static_pressure_pa[i];
                const real V = chamber_volume_m3[i];
                const real Rs = g_specific_gas_constant_j_per_kg_k;
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
                const real Rs = g_specific_gas_constant_j_per_kg_k;
                const real Ts = chamber_static_temperature_k[i];
                const real V = chamber_volume_m3[i];
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
                const real X = (2.0_r / (g_gamma - 1.0_r));
                const real Pi = chamber_total_pressure_pa[i];
                const real Pj = chamber_total_pressure_pa[j];
                const real Pt = fmax(Pi, Pj);
                const real Ps = fmin(Pi, Pj);
                const real direction = Pi > Pj ? 1.0_r : -1.0_r;

                /*
                 *      y - 1                        3
                 *      ----- = 0.3333... where y = ---
                 *        y                          2
                 * Term
                 */

                static_assert(g_gamma == 3.0_r / 2.0_r);
                const real Y = cuberoot(Pt / Ps);
                const real M = direction * sqrt(X * (Y - 1.0_r));
                nozzle_mach[i] = clamper(M, -1.0_r, 1.0_r);
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
                const real Ps = chamber_static_pressure_pa[i];
                const real Pd = chamber_dynamic_pressure_pa[i];
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
                const real Ts = chamber_static_temperature_k[i];
                const real Td = chamber_dynamic_temperature_k[i];
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
                const real Rs = g_specific_gas_constant_j_per_kg_k;
                const real Tt = chamber_total_temperature_k[i];
                const real M = nozzle_mach[i];
                const real X = g_gamma * Rs * Tt;
                const real Y = 0.5_r * (g_gamma - 1.0_r) * M * M;
                const real u = M * sqrt(X / (1.0_r + Y));
                const real A = chamber_nozzle_real_flow_area_m2[i];
                const real mute = A == 0.0_r ? 0.0_r : 1.0_r;
                nozzle_velocity_m_per_s[i] = u * mute;
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
                const real Pt = chamber_total_pressure_pa[i];
                const real Rs = g_specific_gas_constant_j_per_kg_k;
                const real Tt = chamber_total_temperature_k[i];
                const real M = nozzle_mach[i];
                const real pt = Pt / (Rs * Tt);

                /*
                 *        1                  3
                 *      ----- = 2 where y = ---
                 *      y - 1                2
                 * Term
                 */

                static_assert(g_gamma == 3.0_r / 2.0_r);
                const real C = 1.0_r + 0.5_r * (g_gamma - 1.0_r) * M * M;
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
                const real Tt = chamber_total_temperature_k[i];
                const real M = nozzle_mach[i];
                const real Tns = Tt / (1.0_r + 0.5_r * (g_gamma - 1.0_r) * M * M);
                nozzle_static_temperature_k[i] = Tns;
            }
        }

        fn void calc_nozzle_real_flow_areas()
        {
            for(size_t i = 0; i < N; i++)
            {
                const real r = chamber_nozzle_open_ratio[i];
                const real A = chamber_nozzle_flow_area_m2[i];
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
         * Tsp = Ts,upstream
         *
         */

        fn void calc_nozzle_parcels()
        {
            for(size_t i = 0; i < N; i++)
            {
                const size_t j = i + 1;
                const real mdot = nozzle_mass_flow_rate_kg_per_s[i];
                const real dm = mdot * g_dt_s;
                parcel_mass_kg[i] = dm;
                const real Tsi = chamber_static_temperature_k[i];
                const real Tsj = chamber_static_temperature_k[j];
                const real Tsp = mdot > 0.0_r ? Tsi : Tsj;
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

                static_assert(g_gamma == 3.0_r / 2.0_r);
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
                const real Tsp = parcel_static_temperature_k[i];
                const real Ts0 = chamber_static_temperature_k[j];
                const real Ts1 = (Ts0 * m + Tsp * dm) / (m + dm);
                chamber_static_temperature_k[j] = dm > 0.0_r ? Ts1 : Ts0;
            }
        }

        fn void calc_reverse_energy_transfers()
        {
            for(size_t i = N; i > 0; i--)
            {
                const real dm = parcel_mass_kg[i];
                const real m = chamber_mass_kg[i];
                const real Tsp = parcel_static_temperature_k[i];
                const real Ts0 = chamber_static_temperature_k[i];
                const real Ts1 = (Ts0 * m - Tsp * dm) / (m - dm);
                chamber_static_temperature_k[i] = dm < 0.0_r ? Ts1 : Ts0;
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

                static_assert(g_gamma == 3.0_r / 2.0_r);
                chamber_static_temperature_k[i] *= sqrt(m0 / mi);
                chamber_static_temperature_k[j] *= sqrt(m1 / mj);
                chamber_mass_kg[i] = m0;
                chamber_mass_kg[j] = m1;
            }

            for(size_t i = 0; i < N; i++)
            {
                const size_t j = i + 1;
                const real dm = parcel_mass_kg[i];
                const real u = nozzle_velocity_m_per_s[i];
                const real p = dm * u;
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
                const real Rs = g_specific_gas_constant_j_per_kg_k;
                const real Ts = chamber_static_temperature_k[i];
                const real m = chamber_mass_kg[i];
                const real p = chamber_bulk_momentum_kg_m_per_s[i];
                const real pmax = m * sqrt(g_gamma * Rs * Ts);
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
                const real u = chamber_bulk_momentum_kg_m_per_s[i] / chamber_mass_kg[i];
                const real p = chamber_mass_kg[i] / chamber_volume_m3[i];
                const real q = 0.5_r * p * u * u;
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
                const real u = chamber_bulk_momentum_kg_m_per_s[i] / chamber_mass_kg[i];
                const real Cv = g_cv_j_per_kg_k;
                const real Cp = g_gamma * Cv;
                const real Td = 0.5_r * u * u / Cp;
                chamber_dynamic_temperature_k[i] = Td;
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
            if(piston_chamber_on_fire)
            {
                if(not piston_injection_enabled)
                {
                    piston_chamber_on_fire = false;
                    return;
                }
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

                const real Ts = chamber_static_temperature_k[PY];
                const real Ts0 = g_ambient_temperature_k;
                const real Ps = chamber_static_pressure_pa[PY];
                const real Ps0 = g_ambient_pressure_pa;
                const real Tr = Ts / Ts0;
                const real Pr = Ps/ Ps0;
                const real S = 0.4_r * Tr * Tr / sqrt(sqrt(sqrt(Pr)));

                /*
                 *                 2
                 * Vb = dh * pi * r
                 *
                 */

                const real randomness = 0.25_r;
                const real dh = S * g_dt_s;
                const real drh = dh * (1.0_r + randomness * frand());
                const real h1 = piston_chamber_flame_height_m;
                const real h2 = h1 + drh;
                const real r = piston_chamber_radius_m;
                const real Vb = (h2 - h1) * g_pi_r * r * r;

                /*
                 *      M
                 * p = ---
                 *      V
                 *
                 * Mburned = Vb * p
                 *
                 */

                const real p = M / V;
                const real Mb = Vb * p;

                /*
                 *             Q
                 * Ts = Ts + ------
                 *            M Cv
                 *
                 */

                const real TMb = piston_chamber_mass_burned_m3 + Mb;
                if(TMb / M < 1.0_r)
                {
                    const real MFb = Mb / (1.0_r + g_stoich_air_fuel_ratio);
                    const real Q = MFb * g_energy_octane_j_per_kg;
                    const real Cv = g_cv_j_per_kg_k;
                    const real dTs = Q / (M * Cv);
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
                panic[i] |= chamber_mass_kg[i] <= 0.0_r;
                panic[i] |= chamber_static_temperature_k[i] <= 0.0_r;
                panic[i] |= chamber_static_pressure_pa[i] <= 0.0_r;
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
    struct basic_sparkplugs
    {
        static constexpr real fire_delay_theta_r = 1e-1_r;
        lane<W> engage_theta_r = {};
        mask<W> prev_fired = {};
        mask<W> fired = {};
        mask<W> rising_edge = {};
        real crankshaft_theta_r = 0.0_r;

        fn void calc_fired()
        {
            for(size_t i = 0; i < W; i++)
            {
                prev_fired[i] = fired[i];
                const real theta0_r = modulos(crankshaft_theta_r, g_otto_cycle_r);
                const real theta1_r = engage_theta_r[i] >= g_otto_cycle_r ? (engage_theta_r[i] - g_otto_cycle_r) : engage_theta_r[i];
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
    struct basic_cams
    {
        lane<W> engage_theta_r = {};
        lane<W> ramp_theta_r = {};
        lane<W> open_ratio = {};
        real crankshaft_theta_r = 0.0_r;

        /*      4            1      2      3
         * r = t  [ 35 - 84 t + 70 t - 20 t ]
         *
         */

        fn void calc_open_ratios()
        {
            for(size_t i = 0; i < W; i++)
            {
                real theta0_r = modulos(engage_theta_r[i], g_otto_cycle_r);
                if(theta0_r < 0.0_r)
                {
                    theta0_r += g_otto_cycle_r;
                }
                real theta1_r = modulos(crankshaft_theta_r, g_otto_cycle_r);
                if(theta1_r < theta0_r)
                {
                    theta1_r += g_otto_cycle_r;
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
        real mass_kg = 0.0_r;
        real radius_m = 0.0_r;
        real moment_of_inertia_kg_m2 = 0.0_r;

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

    struct throttle
    {
        static constexpr size_t size = 4;
        std::array<real, size> table = {};

        real lookup(const real open_ratio)
        {
            const size_t last = size - 1;
            const real at = last * open_ratio;
            const size_t index = at;
            const real ratio = at - index;
            const size_t next = index + 1;
            const real delta = table[next] - table[index];
            return table[index] + delta * ratio;
        }
    };

    struct limiter
    {
        real max_angular_velocity_r_per_s = 600.0_r;
        real crankshaft_angular_velocity_r_per_s = 0.0_r;
        real limit_time_s = 0.1;
        real cycles = 0;
        bool limiting = false;

        void update()
        {
            if(not limiting)
            {
                if(crankshaft_angular_velocity_r_per_s > max_angular_velocity_r_per_s)
                {
                    limiting = true;
                }
            }
            if(limiting)
            {
                const real time_s = g_dt_s * cycles;
                if(time_s > limit_time_s)
                {
                    limiting = false;
                    cycles = 0;
                }
                cycles++;
            }
        }
    };

    struct crankshaft
    {
        real angular_velocity_r_per_s = 0.0_r;
        real angular_acceleration_r_per_s2 = 0.0_r;
        real mass_kg = 0.0_r;
        real radius_m = 0.0_r;
        real theta_r = 0.0_r;
        real last_theta_r = 0.0_r;
        real moment_of_inertia_kg_m2 = 0.0_r;

        /*
         * dw = a * dt
         *
         */

        fn void accelerate()
        {
            const real a = angular_acceleration_r_per_s2;
            angular_velocity_r_per_s += a * g_dt_s;
        }

        /*
         * dth = w * dt
         *
         */

        fn void turn()
        {
            last_theta_r = theta_r;
            const real w = angular_velocity_r_per_s;
            theta_r += w * g_dt_s;
        }

        fn bool otto_cycled()
        {
            const real t0 = modulos(last_theta_r, g_otto_cycle_r);
            const real t1 = modulos(theta_r, g_otto_cycle_r);
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

        lane<W> diameter_m = {};
        lane<W> crank_throw_length_m = {};
        lane<W> connecting_rod_length_m = {};
        lane<W> connecting_rod_mass_kg = {};
        lane<W> head_mass_density_kg_per_m3 = {};
        lane<W> head_compression_height_m = {};
        lane<W> head_clearance_height_m = {};
        lane<W> theta0_r = {};
        lane<W> theta_r = {};
        lane<W> sint = {};
        lane<W> cost = {};
        lane<W> pin_x_m = {};
        lane<W> pin_y_m = {};
        lane<W> bearing_x_m = {};
        lane<W> bearing_y_m = {};
        lane<W> volumes_m3 = {};
        lane<W> head_mass_kg = {};
        lane<W> moment_of_inertia_kg_m2 = {};
        lane<W> gas_torque_n_m = {};
        lane<W> inertia_torque_n_m = {};
        lane<W> friction_torque_n_m = {};
        lane<W> total_torque_n_m = {};
        lane<W> chamber_static_pressure_pa = {};
        lane<W> friction_n_m_s2_per_r2 = {};
        real crankshaft_angular_velocity_r_per_s = 0.0_r;
        real crankshaft_theta_r = 0.0_r;

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
                const real r = 0.5_r * diameter_m[i];
                const real h = 2.0_r * head_compression_height_m[i];
                const real p = head_mass_density_kg_per_m3[i];
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
                const real Pg = chamber_static_pressure_pa[i] - g_ambient_pressure_pa;
                const real ar = diameter_m[i] / 2.0_r;
                const real A = g_pi_r * ar * ar;
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
                const real K = friction_n_m_s2_per_r2[i];
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

    template<size_t W, size_t L, size_t S>
    struct pipe
    {
        static constexpr size_t M = L - 1;
        lane<W> piston_connect_m = {};
        real length_m = 0.0_r;
        real mic_position0_m = 0.0_r;
        real mic_position1_m = 0.0_r;
        line pipe_pressure_signal = {};

        pipe()
        {
            reset();
        }

        /*
         * Pipe Junction - mix mass W flow rates, velocities, static temperatures streams.
         *
         * --- +-------+
         *  |  |  0    | ---+
         *  |  +-------+    |
         *  |  +-------+    |
         *  W  |  1    | ---+ --> U0 ... UL-1
         *  |  +-------+    |
         *  |     ...       |
         *  |  +-------+    |
         *  |  |  W-1  | ---+
         * --- +-------+
         */

        lane<W> in_velocity_m_per_s = {};
        lane<W> in_static_density_kg_per_m3 = {};
        lane<W> in_static_temperature_k = {};

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

        lane<L> U_r = {};
        lane<L> U_ru = {};
        lane<L> U_rEs = {};
        lane<L> F_r = {};
        lane<L> F_ru = {};
        lane<L> F_rEs = {};
        lane<M> Ff_r = {};
        lane<M> Ff_ru = {};
        lane<M> Ff_rEs = {};

        lane<L> speed_of_sound_m_per_s = {};
        lane<L> local_speed_of_sound_m_per_s = {};
        lane<L> absolute_speed_of_sound_m_per_s = {};
        lane<L> static_pressure_pa = {};

        /*
         *                1   2    r Rs         1     2
         * rEs = Cv Ts + --- u  = ------- Ts + --- r u
         *                2        y - 1        2
         */

        real calc_specific_energy_density_from_static_temperature(const real r, const real u, const real Ts)
        {
            const real ru = r * u;
            const real Rs = g_specific_gas_constant_j_per_kg_k;
            const real rEs = r * Rs * Ts / (g_gamma - 1.0_r) + ru * ru / (2.0_r * r);
            return rEs;
        }

        /*
         *                     2
         *         Ps       r u
         * rEs = ------- + -----
         *        y - 1      2
         */

        real calc_specific_energy_density_from_static_pressure(const real r, const real u, const real Ps)
        {
            const real rEs = Ps / (g_gamma - 1.0_r) + 0.5_r * r * u * u;
            return rEs;
        }

        /*
         *                            1     2
         * Ps = (y - 1) * r * [ Es - --- * u ]
         *                            2
         */

        real calc_static_pressure_from_specific_energy(const real r, const real u, const real Es)
        {
            const real Ps = (g_gamma - 1.0_r) * r * (Es - 0.5_r * u * u);
            return Ps;
        }

        void as_cell(const size_t i, const real r, const real u, const real Ts)
        {
            U_r[i] = r;
            U_ru[i] = r * u;
            U_rEs[i] = calc_specific_energy_density_from_static_temperature(r, u, Ts);
        }

        void to_ambient(const size_t i)
        {
            const real r = g_ambient_density_kg_per_m3;
            const real u = 0.0_r;
            const real Ts = g_ambient_temperature_k;
            as_cell(i, r, u, Ts);
        }

        void reset()
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

        void calc_pipe_open_right()
        {
            const size_t Y = L - 2;
            const size_t Z = L - 1;
            const real r = U_r[Y];
            const real ru = U_ru[Y];
            const real u = ru / r;
            const real a = local_speed_of_sound_m_per_s[Y];
            if(u >= a)
            {
                /*
                 * Sonic or Super Sonic exit.
                 * Wave is so fast that next cell pressure value
                 * is overrided with this cell pressure value.
                 *
                 */

                U_r[Z] = r;
                U_ru[Z] = ru;
                U_rEs[Z] = U_rEs[Y];
            }
            else
            {
                /*
                 * Subsonic exit.
                 * Assume ambient conditions.
                 *
                 */

                const real Ps = g_ambient_pressure_pa;
                U_r[Z] = r;
                U_ru[Z] = ru;
                U_rEs[Z] = calc_specific_energy_density_from_static_pressure(r, u, Ps);
            }
        }

        /*
         * Sample two pipe positions to cancel noise,
         * kind of like a Stratocaster's single coil pickup selector
         * when in position 2 and 4.
         */

        real calc_audio_sample()
        {
            const size_t Z = L - 1;
            const size_t x = Z * mic_position0_m / length_m;
            const size_t y = Z * mic_position1_m / length_m;
            return static_pressure_pa[x] + static_pressure_pa[y];
        }

        void inject()
        {
            for(size_t i = 0; i < W; i++)
            {
                const real r = in_static_density_kg_per_m3[i];
                const real u = in_velocity_m_per_s[i];
                const real Ts = in_static_temperature_k[i];
                const real ratio = piston_connect_m[i] / length_m;
                as_cell(ratio * L, r, u, Ts);
            }
        }

        /*
         *  F = [ rr \ ruu + Ps \ u(rEs + Ps) ]
         */

        void calc_fluxes()
        {
            for(size_t i = 0; i < L; i++)
            {
                const real r = U_r[i];
                const real ru = U_ru[i];
                const real rEs = U_rEs[i];
                const real u = ru / r;
                const real Es = rEs / r;
                const real Ps = calc_static_pressure_from_specific_energy(r, u, Es);
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

        void calc_conserved()
        {
            const real dx_m = length_m / static_cast<real>(L);
            const real dt_s = g_dt_s / S;
            const real dt_dx = dt_s / dx_m;
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

        real calc_static_pressure(const size_t i)
        {
            const real r = U_r[i];
            const real ru = U_ru[i];
            const real rEs = U_rEs[i];
            const real u = ru / r;
            const real Es = rEs / r;
            const real Ps = calc_static_pressure_from_specific_energy(r, u, Es);
            return Ps;
        }

        void calc_static_pressures()
        {
            for(size_t i = 0; i < L; i++)
            {
                const real Ps = calc_static_pressure(i);
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
                const real Al = absolute_speed_of_sound_m_per_s[i];
                const real Ar = absolute_speed_of_sound_m_per_s[j];
                const real alpha = fmax(Al, Ar);
                Ff_r  [i] = 0.5_r * ((F_r  [i] + F_r  [j]) - alpha * (U_r  [j] - U_r  [i]));
                Ff_ru [i] = 0.5_r * ((F_ru [i] + F_ru [j]) - alpha * (U_ru [j] - U_ru [i]));
                Ff_rEs[i] = 0.5_r * ((F_rEs[i] + F_rEs[j]) - alpha * (U_rEs[j] - U_rEs[i]));
            }
        }

        void update()
        {
            inject();
            calc_speed_of_sounds();
            calc_local_speed_of_sounds();
            calc_absolute_speed_of_sounds();
            calc_static_pressures();
            for(size_t i = 0; i < S; i++)
            {
                calc_pipe_open_right();
                calc_fluxes();
                calc_flux_faces();
                calc_conserved();
            }
        }

        void gather_pipe_pressure_signal()
        {
            pipe_pressure_signal.clear();
            for(size_t i = 0; i < L; i++)
            {
                const real Ps = static_pressure_pa[i];
                pipe_pressure_signal.push_back(Ps);
            }
        }

        fn void calc_speed_of_sounds()
        {
            for(size_t i = 0; i < L; i++)
            {
                const real Ps = static_pressure_pa[i];
                const real C = sqrt(g_gamma * Ps / U_r[i]);
                speed_of_sound_m_per_s[i] = C;
            }
        }

        fn void calc_local_speed_of_sounds()
        {
            for(size_t i = 0; i < L; i++)
            {
                const real a = U_ru[i] / U_r[i];
                local_speed_of_sound_m_per_s[i] = a;
            }
        }

        fn void calc_absolute_speed_of_sounds()
        {
            for(size_t i = 0; i < L; i++)
            {
                const real a = local_speed_of_sound_m_per_s[i];
                const real c = speed_of_sound_m_per_s[i];
                absolute_speed_of_sound_m_per_s[i] = fabs(a) + c;
            }
        }
    };

    #define FLUIDS(X) \
        X(chamber_volume_m3) \
        X(chamber_nozzle_real_flow_area_m2) \
        X(chamber_static_pressure_pa) \
        X(chamber_static_temperature_k) \
        X(nozzle_static_temperature_k) \
        X(nozzle_static_density_kg_per_m3) \
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

    using grid = std::vector<line>;

    struct diags
    {
        grid front = grid(g_diags_size);
        grid back = grid(g_diags_size);
    };

    struct dc_filter
    {
        real x_prev = 0.0_r;
        real y_prev = 0.0_r;
        real alpha = 0.0_r;

        dc_filter()
        {
            set_cutoff_frequency(5.0_r);
        }

        void set_cutoff_frequency(const real cutoff_freq_hz)
        {
            const real rc = 1.0_r / (2.0_r * g_pi_r * cutoff_freq_hz);
            alpha = rc / (rc + g_dt_s);
        }

        real filter(const real x)
        {
            const real y = alpha * (y_prev + x - x_prev);
            x_prev = x;
            y_prev = y;
            return y;
        }
    };

    struct gain_filter
    {
        real ratio = 1.0_r;

        real filter(const real x)
        {
            return x * ratio;
        }
    };

    struct clamp_filter
    {
        real filter(const real x)
        {
            return clamper(x, -1.0_r, 1.0_r);
        }
    };

    const line g_impulse = {
        0.000028, -0.000030, 0.000031, -0.000033, 0.000034, -0.000036, 0.000037, -0.000039, 0.000040, -0.000042, 0.000044, -0.000046, 0.000048, -0.000050, 0.000052, -0.000054, 0.000056, -0.000058, 0.000060, -0.000063, 0.000065, -0.000068, 0.000070, -0.000073, 0.000076, -0.000079, 0.000081, -0.000084, 0.000087, -0.000090, 0.000093, -0.000097, 0.000100, -0.000103, 0.000107, -0.000110, 0.000114, -0.000118, 0.000121, -0.000125, 0.000129, -0.000133, 0.000137, -0.000142, 0.000146, -0.000151, 0.000155, -0.000160, 0.000164, -0.000169, 0.000174, -0.000179, 0.000184, -0.000190, 0.000195, -0.000201, 0.000206, -0.000212, 0.000218, -0.000224, 0.000230, -0.000237, 0.000243, -0.000250, 0.000257, -0.000264, 0.000271, -0.000279, 0.000286, -0.000294, 0.000302, -0.000311, 0.000319, -0.000328, 0.000337, -0.000346, 0.000355, -0.000365, 0.000375, -0.000386, 0.000396, -0.000407, 0.000418, -0.000430, 0.000442, -0.000455, 0.000468, -0.000481, 0.000495, -0.000509, 0.000524, -0.000539, 0.000555, -0.000572, 0.000589, -0.000607, 0.000626, -0.000646, 0.000667, -0.000688, 0.000711, -0.000735, 0.000761, -0.000788, 0.000816, -0.000847, 0.000880, -0.000915, 0.000953, -0.000995, 0.001040, -0.001090, 0.001146, -0.001210, 0.001282, -0.001367, 0.001469, -0.001596, 0.001760, -0.001988, 0.002339, -0.002983, 0.004715, -0.042354, -0.002361, 0.001148, 0.000734, -0.000083, 0.001500, -0.000687, 0.002051, -0.001261, 0.002699, -0.002077, 0.003822, -0.003869, 0.007416, -0.017900, -0.016933, 0.005648, -0.001045, 0.001214, 0.002679, -0.031037, -0.003070, 0.003451, 0.000024, 0.002112, 0.000893, -0.027258, 0.001503, 0.001761, 0.000655, -0.004511, 0.006395, -0.008128, -0.018708, 0.006919, -0.000971, 0.004604, -0.001652, 0.032017, 0.006589, -0.001354, 0.003117, -0.005089, -0.019754, 0.000085, 0.005485, -0.008882, 0.005540, -0.006535, 0.033626, 0.023457, -0.007225, 0.006971, -0.004153, -0.015464, 0.001900, 0.001695, -0.005079, 0.007834, 0.077609, -0.007658, 0.019889, 0.040200, -0.009383, 0.004776, -0.003616, 0.002796, 0.000172, -0.009220, 0.064130, -0.005378, 0.017512, 0.030225, -0.011337, 0.004208, -0.005386, 0.000937, 0.010359, 0.027093, -0.020046, 0.042271, -0.003561, -0.009258, -0.001443, -0.000920, -0.008932, -0.003845, 0.026366, -0.015544, 0.030439, 0.006740, -0.005165, -0.004578, 0.018473, 0.034563, 0.001170, -0.009770, -0.001490, -0.016098, -0.001502, -0.007095, 0.000313, 0.026248, -0.008325, -0.001226, -0.006186, -0.001429, -0.010045, -0.008117, 0.012246, 0.054544, -0.006552, 0.012291, 0.014317, -0.013106, 0.002000, -0.011357, 0.021270, 0.003192, 0.015892, -0.006928, -0.002986, -0.013948, -0.004157, -0.008491, 0.017516, -0.010183, 0.019686, 0.000857, -0.007068, -0.006245, -0.003001, -0.011510, 0.026341, 0.004819, 0.007867, 0.005143, -0.012173, 0.022885, -0.008110, -0.002458, -0.010060, -0.001940, -0.013048, -0.010728, -0.006681, -0.002614, 0.008533, -0.006177, 0.019049, -0.016950, 0.003738, 0.008432, -0.017701, -0.007873, -0.006520, -0.004633, -0.007691, -0.000939, 0.012690, -0.017515, -0.003993, -0.008311, 0.006906, -0.004823, -0.006944, -0.006136, -0.003983, -0.007931, 0.000571, 0.000225, 0.007747, 0.001041, -0.007989, -0.007824, 0.000282, 0.003374, -0.011752, -0.001102, -0.005296, -0.004511, -0.009493, 0.006481, 0.005848, -0.001000, -0.009805, -0.004208, -0.006704, -0.000827, -0.011124, -0.012471, -0.005354, -0.001618, -0.009595, 0.014062, 0.015609, 0.000767, 0.000044, -0.005640, -0.001751, -0.007264, 0.004812, -0.008666, 0.003176, -0.013457, 0.003746, 0.002956, 0.012041, -0.003489, 0.009181, -0.003203, 0.000144, 0.024732, 0.001496, -0.010815, 0.012227, 0.003479, 0.010618, 0.007356, -0.012405, 0.001724, -0.012615, 0.015259, 0.021565, 0.000751, 0.004400, 0.010718, 0.007126, 0.007662, -0.010993, -0.006237, -0.008336, -0.001683, -0.000668, 0.006280, 0.005320, -0.009571, -0.000256, -0.004064, -0.008022, -0.002748, -0.009062, 0.025445, 0.004798, -0.008833, -0.004323, -0.011360, 0.004305, 0.000291, -0.005710, -0.000090, -0.002090, 0.019034, -0.004474, -0.002704, -0.002225, -0.007578, 0.016369, -0.000132, -0.008842, -0.002406, -0.007930, 0.014106, -0.005687, 0.004570, -0.012783, 0.013067, 0.006546, -0.007353, -0.014336, -0.006662, 0.010711, -0.004512, 0.000894, -0.017394, -0.003741, -0.006117, -0.001286, -0.012649, -0.010578, -0.007275, -0.003622, -0.007935, -0.006393, 0.006080, -0.003476, -0.003305, 0.000295, 0.003925, -0.007214, -0.001899, -0.012098, 0.012259, 0.012686, 0.001892, -0.015533, -0.008929, 0.004046, 0.007281, -0.014267, -0.007037, 0.007636, 0.004458, -0.005500, -0.017130, -0.011762, 0.003731, -0.002683, 0.002465, -0.002635, 0.010446, 0.012180, -0.004249, -0.006332, 0.002608, 0.014221, -0.007337, -0.007619, -0.000967, 0.010992, 0.004491, -0.004385, 0.006785, 0.027559, 0.005183, -0.000944, -0.003772, 0.002375, -0.007481, -0.011809, 0.003082, 0.010495, 0.000595, 0.006033, 0.013446, 0.002164, -0.012010, 0.012015, 0.015480, -0.016132, 0.001995, 0.010585, 0.013128, 0.009227, 0.006187, 0.017670, 0.017940, 0.006794, -0.006721, -0.003508, -0.005932, -0.006328, -0.003950, -0.002625, -0.007048, 0.022957, 0.011088, 0.000936, -0.011364, 0.003926, -0.000834, -0.008998, 0.000029, -0.001838, 0.003735, 0.004813, -0.011387, -0.008588, 0.007559, -0.011257, -0.005499, 0.004117, -0.001786, -0.008627, 0.000321, 0.015385, -0.008028, 0.005479, 0.001247, -0.003476, -0.005572, -0.007028, 0.001932, 0.009971, -0.009676, -0.013539, 0.001157, -0.012627, 0.008205, 0.002315, 0.007738, 0.015599, -0.001137, 0.003014, 0.004595, 0.004709, -0.001101, 0.009723, 0.000550, -0.005890, 0.009097, 0.026656, -0.001995, 0.020619, 0.016080, -0.007802, -0.008138, -0.002976, 0.016189, 0.021692, 0.002129, 0.007870, -0.012752, -0.003675, -0.000874, 0.002367, 0.002676, 0.001830, -0.006161, 0.010657, 0.019326, 0.011089, 0.001951, 0.001021, -0.001462, -0.008302, 0.001716, 0.009993, 0.004689, -0.012809, 0.004198, 0.012816, -0.007857, 0.000273, 0.005943, 0.004763, 0.002727, -0.013666, 0.001755, 0.005166, -0.001333, -0.004083, -0.009399, -0.006290, -0.010574, -0.001481, 0.001727, -0.007428, 0.007894, 0.008605, -0.001526, 0.008424, 0.005009, -0.004754, -0.009622, 0.005821, 0.005437, -0.006389, -0.009236, 0.003698, -0.001899, -0.000695, 0.002327, -0.014142, -0.003515, -0.014891, 0.013596, 0.013225, -0.015671, 0.004976, -0.002388, 0.000993, 0.009753, 0.003775, -0.009225, 0.001139, 0.005715, -0.001035, 0.013704, -0.003493, -0.015305, -0.008611, -0.009169, -0.006208, -0.001168, -0.008042, -0.003358, -0.008876, 0.004285, 0.001899, -0.009948, 0.005889, -0.008215, 0.001398, 0.019472, -0.005671, 0.002683, -0.007034, -0.001902, 0.001078, -0.003905, -0.010639, 0.003795, -0.008638, -0.002902, 0.010712, -0.008659, -0.001424, -0.000144, 0.001393, -0.001118, -0.008715, 0.001355, -0.002856, 0.000570, -0.011979, -0.001916, 0.000831, -0.011840, 0.005797, 0.008447, -0.009099, -0.007433, -0.004530, -0.004309, -0.008694, -0.000067, -0.008419, 0.007291, -0.006405, 0.002545, -0.006404, 0.005718, -0.001300, 0.013341, 0.012583, -0.008146, -0.000108, -0.006187, -0.002298, 0.006656, -0.007448, -0.005107, 0.005730, 0.007275, 0.001850, -0.002865, -0.005645, 0.001764, 0.006782, -0.002215, -0.007274, -0.003299, -0.006687, 0.005860, -0.008756, -0.004494, -0.005297, 0.006909, 0.017144, -0.003487, -0.005140, -0.000364, -0.009167, -0.001891, -0.004166, -0.009338, -0.000059, -0.003003, -0.006797, 0.000240, -0.000011, -0.011692, 0.012666, -0.006014, -0.004438, -0.011118, -0.009013, 0.005562, -0.000812, -0.008446, 0.003696, 0.003774, 0.006935, -0.007500, 0.005827, 0.002108, -0.003650, -0.002037, 0.006334, -0.008364, 0.016108, 0.001676, -0.005525, 0.018551, 0.009582, -0.002962, 0.010244, -0.003249, -0.004792, -0.003174, -0.008605, -0.003636, 0.004748, 0.005890, 0.007994, 0.001957, 0.000805, -0.001771, 0.002413, 0.003451, 0.001584, 0.005515, -0.002981, 0.013681, -0.000508, 0.010366, 0.007760, 0.008924, 0.008849, -0.004575, -0.004285, -0.007625, -0.013232, -0.004544, 0.009821, 0.013647, 0.005163, 0.010273, -0.007456, 0.008278, 0.004717, -0.004808, -0.000603, -0.000049, 0.000160, 0.002857, -0.002875, 0.006816, 0.008237, 0.011990, -0.005719, 0.007491, 0.007628, -0.001277, 0.001540, -0.003474, 0.019355, 0.014981, 0.006842, -0.002513, -0.003795, 0.002570, -0.003646, -0.009801, -0.002846, 0.010389, -0.005025, 0.000685, -0.004954, 0.007485, 0.016762, -0.002610, -0.007114, -0.009003, 0.001042, -0.008555, 0.002861, 0.013629, 0.009105, 0.005973, -0.000329, 0.003793, -0.000469, -0.010915, -0.004875, -0.006968, -0.008699, -0.007380, 0.000796, -0.005672, 0.021732, -0.003480, -0.002589, -0.007001, -0.000170, 0.005000, -0.001027, 0.008756, 0.000872, 0.004062, 0.005132, 0.008560, 0.000265, -0.004206, -0.008084, -0.009419, 0.007280, -0.006490, -0.003531, -0.003278, 0.010964, 0.005905, -0.000782, -0.002341, -0.001657, -0.005804, -0.006557, -0.003278, 0.008041, 0.008133, -0.009808, 0.000547, 0.003384, -0.000500, -0.006111, 0.000464, -0.010957, -0.004618, -0.000100, -0.000389, 0.018872, -0.011684, -0.002546, 0.003824, -0.006857, -0.000960, -0.013718, 0.001649, -0.002181, 0.003584, 0.003807, -0.004214, 0.006284, -0.002678, -0.005112, -0.007607, -0.004347, -0.003473, 0.004039, -0.005875, 0.003083, 0.004838, 0.005444, 0.003299, -0.003266, -0.001775, -0.005123, -0.002473, 0.007584, -0.001069, 0.009626, 0.006708, -0.000700, 0.003408, -0.009925, 0.000068, 0.005402, -0.000763, -0.000931, 0.004407, 0.000484, 0.008775, 0.000713, -0.002148, -0.008497, -0.000602, -0.003338, -0.005218, 0.012799, -0.000081, 0.008283, -0.001640, -0.007877, -0.001693, -0.007095, -0.000561, 0.001443, -0.000846, 0.000456, 0.002948, 0.004224, 0.001219, 0.002343, 0.000780, -0.004726, 0.002618, 0.004762, 0.003060, -0.000436, 0.003437, -0.002970, -0.003983, -0.001902, -0.004594, -0.003705, -0.004485, 0.000457, -0.003812, -0.008437, -0.005239, -0.007058, -0.005846, -0.004513, -0.001485, 0.003381, 0.003841, 0.001932, 0.001882, -0.009028, 0.001453, 0.003522, -0.006369, -0.002707, -0.008987, -0.001262, -0.003135, -0.000354, -0.000161, 0.000687, -0.007104, 0.007536, 0.000511, 0.005322, 0.005379, -0.000808, -0.001048, -0.007882, 0.005948, -0.001071, -0.003045, -0.009942, -0.004985, 0.004284, -0.006173, 0.002945, 0.004817, -0.006915, -0.011947, -0.000536, 0.006772, -0.007153, -0.001947, -0.005916, -0.006169, 0.002785, 0.004848, -0.000615, -0.010406, -0.001483, 0.001662, -0.001078, -0.000763, -0.004468, -0.005091, -0.011554, 0.004065, -0.003695, 0.002615, -0.005821, 0.001475, -0.002898, 0.006145, 0.001395, -0.001932, -0.007915, 0.004505, -0.000501, 0.011306, 0.002712, -0.005316, -0.000413, 0.004963, 0.007920, 0.004059, -0.005346, -0.000610, -0.003154, 0.010990, 0.001415, 0.005220, -0.007521, 0.006702, 0.004887, -0.001635, 0.004616, -0.002906, -0.005021, 0.002469, 0.006896, 0.000765, -0.004993, -0.003440, -0.001791, 0.001993, 0.009013, 0.005992, -0.004943, 0.000638, 0.004751, -0.002706, 0.004622, 0.009585, -0.002796, 0.004107, 0.008902, 0.013084, -0.003251, -0.003643, 0.002809, 0.006900, 0.017828, 0.006467, -0.003326, -0.000888, 0.001690, 0.001113, 0.007078, -0.005468, 0.000859, -0.000589, -0.000352, -0.005251, -0.001302, -0.005137, 0.003568, 0.007628, 0.006126, 0.003825, -0.001229, -0.004009, -0.000939, 0.000340, 0.006169, 0.000850, -0.001195, 0.002038, 0.000018, -0.003623, 0.009340, 0.003990, 0.000022, 0.006188, -0.001638, -0.007370, -0.000872, 0.001918, 0.002743, 0.007154, -0.003769, 0.007022, -0.007463, 0.004649, -0.008075, -0.004564, -0.002233, -0.004946, -0.001526, 0.004024, -0.004408, 0.001343, 0.003890, 0.001929, 0.000332, -0.001527, -0.005337, -0.004202, -0.003300, 0.002647, 0.000710, -0.003735, -0.002194, 0.008999, -0.003155, 0.004098, -0.000401, 0.001246, -0.006768, 0.003672, 0.002099, -0.001225, -0.002354, -0.003797, 0.004293, 0.002271, -0.004114, -0.001825, 0.001816, -0.005224, -0.003434, 0.009547, -0.005731, -0.003293, 0.002694, 0.005269, 0.004385, -0.005526, -0.001797, -0.002096, 0.000838, -0.000205, -0.000330, -0.005445, 0.001414, -0.002917, 0.003097, -0.004458, 0.001948, 0.001792, -0.003040, 0.004139, 0.003422, -0.008233, 0.009180, -0.006827, 0.005709, -0.001700, -0.003866, 0.004720, 0.007054, -0.008376, 0.006992, 0.001758, -0.005291, -0.003283, 0.002263, 0.003268, -0.003710, 0.001865, -0.007430, 0.000351, -0.004260, -0.000760, 0.007190, -0.001973, -0.001574, -0.005340, -0.003178, 0.002234, -0.004614, 0.005684, 0.006128, -0.000080, -0.001981, 0.004522, -0.001962, -0.003947, -0.001769, -0.002310, 0.000850, -0.000664, 0.000755, 0.006598, 0.004347, -0.000831, 0.005646, 0.000596, -0.005116, 0.001871, -0.006592, -0.004944, 0.002574, -0.005874, 0.003209, -0.003640, -0.001465, -0.003556, -0.001386, -0.004404, 0.000697, -0.003157, 0.002678, 0.000423, -0.004376, 0.002956, 0.000160, -0.003149, 0.003601, 0.003124, -0.003927, 0.001553, 0.001331, 0.003659, -0.004849, 0.003385, 0.008925, -0.004068, -0.000781, 0.001521, 0.000314, -0.004763, -0.000104, -0.000624, -0.002277, 0.000461, 0.000033, 0.005734, -0.007491, -0.001233, 0.003517, -0.003491, 0.002166, -0.002954, -0.004071, -0.002714, -0.003335, -0.005182, 0.000938, -0.005685, -0.003632, -0.004651, -0.001200, -0.001380, -0.000430, -0.001371, 0.001853, -0.006408, 0.002217, -0.004181, 0.002138, 0.002109, -0.001975, -0.003660, -0.004048, 0.002079, 0.001651, -0.003307, -0.000603, -0.005880, -0.003595, -0.002189, -0.001791, -0.004620, -0.006983, 0.001717, -0.004827, -0.006534, 0.002373, 0.000868, -0.001635, -0.005562, -0.007534, -0.005132, 0.000286, -0.006468, -0.000294, -0.001935, 0.002704, 0.007276, -0.003285, 0.003865, 0.004129, -0.001479, -0.003637, -0.002152, -0.003527, 0.006384, -0.000419, -0.001905, 0.000356, 0.004432, 0.002230, -0.007377, 0.003018, 0.002195, -0.000849, -0.004433, 0.004405, 0.008600, -0.000784, 0.002244, 0.003413, 0.001565, 0.005959, 0.002269, 0.000345, -0.001464, 0.003051, 0.009309, 0.001797, -0.003390, 0.008670, 0.007661, -0.002750, 0.001126, 0.008079, 0.002343, -0.000856, 0.000163, 0.007256, -0.000343, 0.008646, -0.001365, -0.000770, 0.006435, 0.002558, 0.001325, -0.001782, -0.001496, 0.001169, -0.006298, 0.001338, 0.001325, 0.006268, 0.008084, -0.000978, 0.002487, 0.000606, -0.001855, -0.001803, 0.002697, 0.003134, 0.000626, -0.000985, 0.001919, 0.006352, -0.002269, 0.000807, -0.002226, -0.000930, -0.003021, 0.001089, -0.001006, 0.004952, -0.000622, 0.001862, 0.002628, 0.003628, 0.001614, -0.002264, -0.000484, -0.002710, -0.000460, 0.001877, 0.003398, 0.002368, 0.001511, 0.002515, -0.000936, 0.001260, 0.006288, 0.000088, -0.001192, -0.002195, -0.000899, 0.001400, -0.002186, 0.003264, -0.006819, 0.003807, -0.001342, -0.001433, -0.001670, 0.002133, -0.002073, 0.000807, 0.000613, 0.000703, 0.005178, -0.001569, 0.000899, 0.001724, -0.001224, 0.001321, -0.000470, 0.001631, 0.004887, -0.000743, 0.001850, -0.000139, 0.004881, -0.001500, 0.001688, 0.001427, -0.003767, -0.002179, -0.001610, -0.005879, 0.000978, 0.000885, 0.003224, -0.002615, 0.001217, -0.007824, 0.003341, -0.006325, -0.001257, -0.002658, -0.002997, 0.001599, -0.000747, 0.003948, -0.004759, -0.005593, -0.003647, -0.003751, -0.001053, -0.001865, 0.002831, -0.004731, -0.001043, 0.001928, -0.001419, -0.000643, -0.000320, -0.001942, -0.005545, 0.003291, -0.004703, 0.003490, -0.000043, -0.005034, -0.001518, -0.003970, -0.002179, -0.006363, 0.000621, -0.002599, -0.001645, 0.000889, -0.002816, 0.001957, -0.003145, 0.002031, -0.003342, -0.001213, -0.001955, -0.003311, 0.000895, -0.003139, -0.000695, -0.004227, -0.001171, 0.000797, -0.002253, 0.001104, 0.000535, -0.002357, 0.000477, -0.000426, 0.002687, -0.000815, 0.001515, 0.001208, -0.000456, 0.003938, 0.000561, 0.002303, -0.002019, 0.003943, 0.000070, 0.000831, 0.005919, 0.000097, 0.002662, -0.001255, -0.003167, 0.001549, -0.000846, -0.001726, 0.002844, 0.002905, 0.001481, -0.001842, 0.000883, 0.007083, -0.001367, 0.001208, -0.000879, 0.001847, -0.005688, 0.004857, -0.001543, 0.001455, -0.006834, -0.002292, 0.003316, -0.004939, 0.000323, -0.002436, 0.002228, -0.003709, 0.002876, -0.001397, -0.002688, -0.001812, -0.001784, -0.002800, -0.005007, 0.000610, -0.004325, -0.000532, -0.004665, -0.000294, -0.002979, -0.002282, -0.003109, 0.004470, -0.002725, -0.003679, 0.000887, 0.000710, -0.001003, -0.003082, 0.001062, -0.004754, 0.001186, -0.003502, -0.000056, 0.000934, -0.004934, 0.000733, 0.006257, -0.000988, -0.005551, -0.002315, 0.000034, -0.000480, -0.000428, 0.000241, -0.002646, -0.002267, -0.000022, -0.001060, -0.002774, -0.001437, -0.002961, 0.004788, 0.003490, 0.004627, 0.003710, 0.003419, -0.002239, -0.005269, 0.006444, 0.002722, -0.002245, 0.001917, -0.000816, -0.002699, 0.005199, -0.001839, -0.002204, -0.000840, -0.000016, 0.000979, -0.001088, -0.002114, -0.002498, 0.004001, 0.000877, 0.002384, 0.003149, 0.000829, 0.005661, 0.003714, 0.007295, 0.000409, 0.000599, 0.001703, 0.004112, -0.000348, 0.001029, -0.002441, 0.000381, 0.001514, 0.000954, 0.001165, -0.001512, 0.001568, -0.002065, 0.003040, 0.002838, 0.002496, 0.005364, -0.002563, 0.002161, 0.001782, 0.001292, 0.000133, -0.001201, 0.000738, 0.003138, -0.002321, 0.000875, -0.006310, 0.003550, 0.004216, 0.000392, -0.001824, -0.005148, -0.000118, 0.002088, 0.000748, 0.001400, -0.001264, -0.001023, 0.002417, 0.002380, 0.000013, -0.004789, 0.000686, 0.002343, 0.001272, 0.003569, -0.003026, 0.000289, -0.001195, 0.006555, 0.002622, -0.003652, -0.002263, -0.003133, -0.000431, 0.001323, -0.002832, 0.000477, 0.000019, 0.000924, 0.003987, -0.004366, -0.001458, -0.002665, 0.003212, 0.001485, 0.001454, 0.000317, 0.001822, 0.001409, 0.004940, 0.000362, -0.001985, -0.000031, -0.001590, -0.002275, -0.001643, -0.002825, -0.000587, 0.000659, -0.001943, 0.001854, -0.000213, -0.000248, -0.002240, -0.000031, 0.002157, -0.001704, 0.000231, 0.000602, 0.000718, 0.002275, 0.004326, -0.001738, 0.002944, -0.001802, -0.002526, -0.001972, -0.001475, 0.001654, 0.000466, -0.002778, 0.002064, -0.001113, 0.000227, 0.001737, -0.001606, 0.000023, -0.001258, 0.000955, -0.003851, -0.000994, 0.001530, 0.000325, 0.002076, 0.002712, -0.001099, 0.000156, 0.002579, 0.000887, -0.000630, 0.001444, -0.002754, -0.000546, -0.003047, 0.000077, -0.004066, -0.000813, 0.000097, -0.003356, -0.001022, -0.002792, -0.003070, 0.001148, -0.000160, 0.003923, -0.000785, 0.000169, 0.000266, 0.000532, 0.001350, -0.001271, -0.001383, 0.002430, 0.001062, -0.002678, 0.000654, 0.000144, 0.004143, -0.002713, -0.000851, -0.000937, 0.000834, 0.000009, 0.000368, -0.001103, 0.001740, -0.001520, -0.004279, 0.000607, 0.002565, 0.001822, 0.000909, -0.004223, -0.003670, -0.001719, -0.000392, -0.002047, 0.000473, -0.003278, -0.000779, -0.000000, -0.002592, -0.001575, -0.002276, 0.002112, -0.002087, 0.000229, -0.001624, -0.000379, 0.001001, 0.000325, 0.000579, -0.002234, -0.001549, -0.002968, 0.001163, -0.002942, -0.001162, -0.000582, -0.001660, 0.000419, -0.000249, -0.001369, -0.001515, 0.001266, -0.000003, -0.000293, -0.000719, 0.001073, -0.000341, -0.000474, 0.001133, -0.002919, -0.000805, -0.000400, -0.000373, -0.001870, 0.002041, -0.003684, -0.001337, 0.000377, -0.002291, -0.000042, 0.001559, -0.000586, -0.002546, 0.000779, -0.002939, 0.000546, -0.001175, -0.000106, -0.001277, 0.000275, -0.001816, -0.001779, 0.001964, -0.001323, 0.001896, -0.003010, -0.001604, -0.003333, -0.002017, -0.000726, -0.001022, -0.002804, -0.000796, -0.001144, 0.000657, -0.001944, -0.002255, 0.004951, -0.002044, 0.001448, -0.002522, -0.000728, 0.001244, 0.000605, 0.001477, -0.001692, 0.001708, 0.002538, -0.003131, -0.000508, 0.002523, 0.003566, 0.005046, -0.000476, -0.000213, -0.000105, 0.001689, -0.001502, -0.002221, 0.000164, -0.000609, -0.000986, -0.000105, -0.001453, -0.002010, 0.001793, 0.001330, 0.003503, 0.002832, 0.002424, -0.001017, -0.003701, 0.001437, 0.001816, 0.001490, 0.002494, 0.000601, 0.000199, -0.001388, -0.001473, -0.000537, -0.000079, 0.004673, -0.000025, 0.004096, 0.001548, 0.000300, 0.003678, -0.000032, 0.001289, 0.004907, 0.002014, -0.000396, -0.001960, 0.000675, 0.002383, -0.000127, 0.000465, -0.001928, 0.000434, 0.001493, -0.001676, -0.003801, 0.003279, 0.002445, 0.001476, 0.002716, 0.002861, 0.000625, 0.003614, -0.002944, -0.001905, -0.001032, 0.002978, -0.000325, 0.000586, 0.000895, 0.000678, 0.003505, -0.001843, -0.003512, -0.000963, 0.003621, 0.002131, -0.000797, 0.000934, 0.001138, 0.000239, -0.000811, 0.000424, 0.000622, 0.001822, 0.001044, -0.001556, 0.003268, -0.000318, 0.000407, -0.000325, 0.001895, 0.000017, -0.000735, 0.000935, -0.002407, -0.001919, 0.002802, 0.000397, -0.003034, -0.001451, -0.002653, -0.001774, 0.002551, -0.001601, 0.003027, -0.001262, -0.000033, 0.000120, -0.001955, -0.002323, -0.001005, -0.001653, -0.000860, -0.000531, 0.000766, -0.001775, 0.000114, 0.000643, 0.001455, 0.000286, 0.000828, 0.001293, -0.000967, 0.000061, 0.000561, -0.000142, -0.001509, 0.000603, 0.000488, 0.000919, 0.000588, 0.001702, -0.002230, -0.001705, 0.000575, -0.000028, -0.003011, -0.001721, 0.000296, -0.000054, -0.000472, 0.000291, 0.001432, 0.001352, 0.000779, 0.000744, -0.002029, 0.000936, -0.002919, 0.000397, 0.000425, -0.000856, -0.000285, -0.001732, -0.001336, -0.000062, 0.000942, -0.001234, -0.001789, -0.001840, -0.000004, 0.000947, 0.000376, 0.000898, 0.000175, -0.001524, -0.000361, -0.000576, -0.000854, -0.000938, -0.000270, 0.000632, -0.000130, -0.001109, 0.002542, 0.001614, 0.000630, 0.001218, -0.001060, -0.002644, -0.000021, -0.001244, 0.000123, -0.001701, -0.001431, -0.002126, 0.000426, -0.001136, -0.000106, -0.001108, -0.000557, 0.000816, 0.000040, -0.000357, -0.000667, 0.000245, 0.000242, -0.001206, -0.000203, -0.000202, -0.000633, -0.000279, -0.000374, -0.001235, -0.001536, -0.004277, 0.000402, -0.000274, 0.001264, -0.000034, -0.000603, 0.000285, -0.000569, -0.002066, 0.000209, -0.001147, -0.000221, -0.001981, -0.001415, -0.000088, -0.000943, 0.001159, -0.000476, 0.000379, -0.001988, -0.001336, -0.000096, -0.000636, 0.002485, -0.000120, -0.000135, -0.001751, 0.000983, -0.001840, -0.001416, -0.001546, -0.002486, 0.000534, -0.000337, -0.002546, -0.000350, -0.000443, -0.000037, -0.001545, -0.001531, 0.000679, -0.001636, -0.001676, -0.000118, 0.000171, -0.001154, -0.002750, -0.001927, -0.001026, 0.000167, -0.002853, 0.000365, -0.001717, -0.000626, -0.002624, -0.001132, 0.001335, -0.001836, 0.002590, -0.001620, 0.000526, -0.000748, 0.000312, -0.001733, -0.002365, -0.001571, 0.000767, -0.001767, 0.000521, -0.002301, -0.001562, 0.001914, -0.000390, 0.001466, 0.000107, -0.000361, 0.001463, 0.001485, 0.004321, -0.000303, 0.001796, 0.000797, 0.000286, 0.002185, -0.000690, -0.000437, -0.001733, -0.000667, 0.002220, 0.001072, -0.001330, 0.004855, -0.000271, 0.002775, -0.000151, 0.001913, 0.001513, -0.001375, 0.004159, 0.004960, 0.001749, -0.001393, 0.000042, 0.001933, -0.001878, 0.001902, 0.002162, 0.000205, -0.001841, -0.000131, 0.001716, 0.003703, -0.000554, 0.001422, 0.002636, 0.000279, 0.002313, -0.000241, 0.002286, -0.001127, 0.000125, 0.001903, 0.000267, -0.000736, -0.002807, 0.000195, 0.001522, 0.002499, 0.000116, -0.001078, 0.000993, 0.003802, 0.002846, 0.002009, 0.001429, 0.000069, -0.000642, 0.000648, -0.000028, -0.000064, -0.000728, -0.000779, -0.000428, -0.000487, 0.001329, 0.000501, 0.002361, 0.001088, -0.000698, 0.000405, 0.002052, -0.001293, -0.001243, 0.004841, -0.000383, -0.001148, -0.000909, -0.000733, -0.000595, -0.001170, -0.001817, -0.000328, -0.001765, 0.000432, -0.000168, 0.003785, -0.000816, -0.001830, 0.002111, 0.000792, -0.001147, -0.000744, -0.001316, -0.001756, -0.000753, -0.000933, -0.000909, -0.001391, -0.001338, -0.001197, 0.001563, -0.000468, -0.001678, -0.000662, -0.001490, 0.000303, 0.000453, 0.000726, -0.001726, 0.001212, -0.000758, -0.001091, -0.001402, 0.000959, 0.000633, 0.000381, -0.001176, 0.000822, -0.001558, -0.001108, -0.000031, 0.000479, 0.000364, -0.000473, -0.000752, -0.001868, -0.000937, -0.000244, -0.000075, -0.000841, -0.000202, 0.000005, -0.001193, 0.000366, -0.001923, -0.001501, 0.000199, 0.001580, -0.000826, -0.000504, 0.000895, 0.001701, -0.001011, -0.001735, -0.000915, -0.000741, 0.000699, 0.001681, 0.000561, 0.000780, -0.000112, -0.000284, 0.000165, -0.002178, -0.000388, -0.001104, -0.000168, 0.000485, -0.000021, 0.000263, 0.000076, -0.000884, -0.000616, 0.000094, 0.000632, -0.002247, -0.000655, -0.000508, 0.000532, -0.000119, -0.000428, -0.000061, -0.000177, 0.001702, 0.000607, -0.000929, 0.000278, -0.001550, 0.000901, 0.000482, 0.000128, -0.001296, -0.000840, 0.000535, -0.001074, -0.001900, -0.001165, -0.000194, -0.000431, 0.000127, 0.000715, -0.000621, 0.001117, 0.000528, 0.000440, 0.000687, 0.000674, 0.000951, -0.000423, 0.000785, -0.000532, -0.000308, 0.000184, -0.000591, 0.001856, -0.001147, 0.000455, -0.000930, -0.001027, 0.000475, 0.001247, -0.001864, -0.000505, 0.001021, 0.000335, -0.002368, -0.001135, 0.000313, -0.000725, -0.000064, 0.000046, -0.000808, -0.000233, -0.001801, 0.001117, -0.000074, -0.000017, 0.001685, -0.001164, -0.000956, -0.001481, -0.001300, -0.000843, -0.001761, -0.001603, -0.001498, -0.000919, -0.000491, -0.001988, 0.000807, -0.000605, -0.000821, -0.000171, -0.000790, -0.000560, -0.000625, 0.001073, -0.001280, -0.000512, -0.001711, -0.000580, -0.001444, -0.000838, -0.000766, -0.000650, -0.000355, -0.001280, -0.001013, -0.001053, -0.001857, -0.000630, 0.000032, 0.000304, -0.000600, -0.000780, -0.000381, -0.000177, -0.001162, -0.002172, -0.002708, -0.001812, 0.000205, 0.001706, 0.001502, -0.000835, 0.002914, 0.001675, 0.001661, -0.001715, 0.001932, -0.001521, -0.001056, -0.000508, -0.000968, 0.000462, -0.001552, 0.000150, -0.001076, 0.000461, 0.002104, -0.000660, -0.001273, -0.000394, 0.002948, 0.001166, 0.003562, 0.001018, -0.000105, 0.000479, 0.001949, 0.001220, 0.002541, 0.000683, -0.001553, 0.000877, 0.000929, 0.001499, 0.000691, -0.000422, 0.001850, 0.001761, 0.003208, 0.000637, -0.000576, 0.000188, 0.000491, -0.000616, 0.001336, -0.000074, -0.000094, 0.000935, 0.001768, 0.000961, 0.001555, 0.001789, -0.001271, 0.000154, 0.000032, 0.001574, -0.000740, -0.002183, 0.001505, 0.000526, 0.000360, -0.000908, -0.001016, 0.000090, 0.001063, 0.001277, -0.000168, -0.001852, -0.000390, 0.000243, 0.001648, 0.002210, 0.000549, -0.000494, -0.000979, 0.001478, 0.000861, 0.000311, -0.001468, -0.000750, 0.002504, 0.000497, -0.000671, 0.000048, 0.000894, -0.000089, 0.000477, 0.000207, -0.000949, -0.000499, -0.000821, 0.000280, 0.001488, 0.000361, 0.002116, -0.000260, -0.001443, 0.000682, -0.000678, 0.000213, -0.000348, 0.001916, 0.000434, 0.000624, -0.001294, -0.000156, -0.000249, -0.001001, -0.000503, -0.000278, 0.000793, -0.001050, 0.001298, 0.000410, -0.000589, -0.001200, -0.000938, 0.000573, 0.001563, 0.000732, -0.001068, -0.000873, -0.000301, 0.000645, -0.000321, -0.000146, -0.000500, 0.000732, -0.000969, -0.001779, 0.000181, 0.000516, -0.000411, -0.001037, 0.000303, 0.000132, 0.000120, -0.000415, -0.001293, -0.000501, -0.000481, -0.000002, -0.000449, -0.000534, -0.000727, -0.000941, 0.000175, 0.000594, 0.002355, -0.000718, 0.000154, 0.000213, -0.001315, 0.000122, -0.000998, -0.000025, -0.000961, 0.000121, 0.000564, 0.000577, -0.000160, -0.001340, -0.000224, -0.000684, 0.000160, 0.000230, -0.000958, -0.000716, 0.000305, 0.000139, -0.002319, 0.000496, 0.000082, -0.000966, -0.000721, -0.000570, -0.001259, 0.000237, -0.000414, -0.000939, 0.000484, 0.000960, -0.001586, 0.000144, -0.000813, -0.000996, -0.000145, -0.000412, -0.000845, -0.001899, 0.000413, -0.001808, -0.000536, -0.000937, 0.000054, 0.000628, 0.000284, 0.000930, -0.001006, -0.000013, -0.000232, -0.000338, -0.001757, -0.000672, -0.000091, -0.000596, -0.000416, -0.000844, 0.001635, 0.000026, 0.000082, -0.000102, -0.000133, 0.000233, -0.000438, -0.001260, -0.000084, -0.000363, 0.000100, -0.001473, 0.000158, -0.001020, 0.000212, -0.000109, 0.000564, -0.000498, 0.001306, -0.000340, 0.000578, -0.000289, -0.000078, 0.000477, 0.002033, -0.000841, -0.000840, -0.000845, -0.000303, -0.001651, -0.000591, -0.000197, -0.000476, -0.000080, -0.000753, -0.000911, -0.000418, -0.000803, -0.000847, -0.000503, -0.000559, -0.000788, -0.000807, 0.000226, 0.000540, -0.000302, 0.000005, -0.000530, -0.000701, 0.000238, -0.001297, -0.000469, -0.000229, 0.000409, -0.000387, -0.000439, -0.000423, -0.000697, -0.001821, 0.000619, -0.002651, -0.000911, -0.001162, -0.000017, 0.000098, 0.000040, -0.002669, -0.000067, -0.000542, -0.000150, -0.000082, 0.001787, -0.001332, -0.000263, 0.001421, 0.000856, -0.000257, 0.001572, -0.001451, 0.000025, -0.000258, -0.000512, -0.001528, 0.001213, 0.000531, 0.000386, 0.001717, 0.001149, -0.000574, -0.000041, 0.000207, 0.001104, 0.001789, 0.001397, 0.000047, 0.000687, 0.000254, 0.001793, -0.000202, 0.000820, -0.000317, -0.000234, 0.000501, -0.000344, 0.000313, 0.000117, 0.001953, 0.001599, 0.003341, -0.000235, 0.000435, 0.000402, 0.000161, 0.000217, 0.002963, 0.000214, -0.000778, 0.001522, -0.000259, 0.000961, -0.000176, 0.000066, 0.000321, 0.000049, 0.000851, -0.001382, 0.000929, 0.000068, -0.000361, 0.000237, 0.002253, -0.001549, 0.000167, 0.000751, 0.000001, 0.000226, 0.001725, 0.001772, -0.000682, 0.001251, -0.000915, 0.000089, 0.000760, 0.001175, 0.000474, -0.000025, 0.000245, -0.000231, 0.000159, 0.000172, -0.000220, 0.000702, 0.000726, -0.000470, -0.000769, -0.000975, 0.001127, 0.001666, -0.000043, -0.000549, -0.001056, 0.000374, -0.000089, 0.000465, 0.000034, 0.000319, 0.000527, -0.001154, -0.000376, 0.000473, 0.000042, -0.000317, 0.000013, 0.000159, -0.001925, 0.000687, -0.000665, 0.000386, 0.000492, 0.000593, -0.000879, 0.001305, -0.000410, -0.001061, -0.000347, -0.000336, -0.000452, -0.001465, -0.001483, -0.002013, -0.000124, -0.000250, -0.000115, -0.000253, 0.001744, 0.000433, -0.000962, 0.000273, 0.000102, -0.000106, -0.001235, -0.000737, -0.001289, 0.000773, -0.000369, -0.000303, -0.000715, 0.001472, -0.000029, -0.000919, -0.000826, -0.000053, -0.000007, -0.000091, 0.000492, -0.000538, 0.000452, -0.000343, 0.000166, 0.000101, 0.000064, -0.000865, -0.001064, -0.000854, -0.000213, -0.000570, 0.000812, -0.000237, 0.000338, -0.000167, -0.001123, -0.000261, 0.000013, 0.000748, -0.000338, 0.000384, 0.000904, -0.000028, -0.000050, 0.000108, 0.000317, -0.000368, -0.001076, -0.000135, -0.001399, 0.000396, -0.000937, -0.000284, 0.000170, 0.000148, -0.000058, -0.000189, -0.000128, -0.000212, -0.000547, -0.000593, -0.000372, 0.000323, -0.000210, -0.000202, 0.000449, -0.000157, -0.000305, -0.000676, 0.000275, -0.000625, -0.000488, -0.000886, -0.000597, -0.000488, -0.000417, -0.000905, -0.001312, -0.000650, -0.000267, 0.000021, 0.000457, 0.000095, -0.000148, 0.000127, 0.000672, 0.000114, -0.000410, 0.000130, -0.001211, 0.000264, -0.000563, -0.000689, 0.000072, -0.000088, -0.000519, 0.000049, -0.000383, -0.000725, 0.000418, -0.000903, -0.000688, -0.000155, 0.000466, 0.000833, 0.000497, 0.000549, 0.000277, -0.000895, 0.000499, -0.000744, -0.000293, -0.000991, 0.000242, -0.000960, -0.000126, 0.000004, -0.000559, -0.000477, -0.000583, -0.000416, -0.000103, -0.000936, 0.000482, -0.000586, 0.000185, -0.000560, -0.001098, -0.000970, 0.000111, -0.000669, -0.000159, -0.000307, -0.001890, 0.000078, -0.000086, -0.000946, -0.000774, -0.000968, -0.000881, -0.000684, -0.000192, -0.000706, 0.000495, -0.001450, -0.000192, 0.000280, -0.000799, -0.000133, 0.000375, -0.000947, 0.000180, 0.000028, 0.000074, -0.000177, 0.000609, -0.000138, -0.000402, -0.000348, 0.000320, -0.000231, 0.000878, 0.000132, 0.000898, 0.000544, 0.000496, 0.000171, 0.000476, -0.000189, 0.002449, 0.000583, 0.000523, -0.000407, 0.000196, 0.000576, 0.000335, -0.000653, 0.001437, 0.000864, 0.000959, 0.000698, 0.000796, 0.000697, 0.001191, -0.000469, 0.000236, 0.001178, 0.000483, 0.000103, 0.001200, -0.000236, 0.001965, 0.001267, -0.000844, 0.000262, 0.000964, 0.000056, 0.001730, -0.000039, -0.000728, 0.001576, 0.000741, -0.000209, -0.000343, -0.000146, -0.000177, 0.000321, -0.001103, 0.000714, 0.001285, -0.001165, -0.000302, -0.000844, -0.000609, 0.000256, -0.000135, 0.000532, 0.000210, 0.000144, -0.000454, 0.000446, 0.000077, 0.000807, -0.000458, 0.000025, -0.000223, 0.000242, 0.000092, -0.000227, -0.000503, 0.000753, -0.000136, 0.000513, 0.000434, -0.000734, 0.000350, -0.000079, -0.000546, -0.000005, -0.000043, 0.000343, -0.000914, -0.000491, 0.000468, 0.000459, 0.000398, 0.000217, 0.000154, 0.000019, -0.000941, 0.000877, -0.000684, 0.000048, 0.000643, 0.000403, 0.000113, 0.000464, -0.000416, -0.000684, -0.000050, 0.000207, 0.000060, -0.000402, -0.000575, -0.000405, -0.000335, -0.000883, 0.000517, -0.000261, 0.000520, 0.000487, -0.000036, -0.000038, 0.000201, 0.000080, -0.000456, 0.000424, 0.000175, 0.000868, -0.000870, -0.000069, 0.000084, 0.000360, -0.000752, 0.000725, -0.000120, -0.000314, 0.000016, -0.000885, 0.000311, -0.000110, -0.000558, -0.000307, -0.000283, -0.000260, -0.000127, -0.000803, 0.000157, 0.000860, -0.000177, 0.000063, -0.000229, -0.000305, -0.000124, -0.000024, -0.000003, 0.000261, -0.000407, 0.000110, -0.000388, -0.000174, 0.000601, -0.000294, -0.000631, -0.000559, -0.000059, -0.000710, -0.000492, -0.000021, 0.000188, -0.000187, -0.000156, 0.000081, 0.000206, 0.000037, -0.000974, -0.000294, 0.000181, -0.000151, -0.000570, -0.000565, 0.000292, 0.000240, -0.000406, -0.000157, -0.000068, 0.000161, -0.000122, -0.000113, 0.000172, -0.000058, 0.000090, 0.000012, -0.000092, -0.000594, -0.000813, -0.000054, -0.000404, 0.000307, 0.000334, -0.000528, -0.000676, -0.000039, 0.000176, -0.000142, -0.000358, 0.000168, 0.000544, -0.000071, 0.000134, -0.000434, -0.000551, 0.000120, -0.000336, -0.000130, -0.000399, -0.000450, -0.000204, -0.000151, -0.000428, 0.000074, -0.000060, -0.000089, -0.000620, 0.000173, 0.000159, -0.000559, -0.000333, -0.000635, -0.000300, -0.000149, -0.000342, -0.000289, 0.000308, -0.000320, 0.000452, 0.000391, 0.000081, -0.000741, -0.000501, -0.000488, -0.000605, -0.000516, -0.000626, -0.000397, -0.000057, -0.000901, 0.000238, -0.000906, -0.000424, -0.000443, -0.001067, -0.000139, -0.000420, -0.000866, -0.001008, -0.000732, 0.000062, 0.000003, -0.000118, -0.000207, -0.000250, -0.000207, -0.000596, -0.000809, -0.001359, -0.000432, 0.000279, -0.000239, -0.000352, -0.000178, -0.001078, -0.000011, -0.000463, -0.000225, -0.000689, -0.000897, -0.000283, -0.001073, -0.000357, -0.000143, -0.000016, -0.000473, -0.000664, -0.000113, 0.000851, 0.000367, -0.000427, 0.000851, 0.000562, 0.000258, -0.000371, -0.000124, 0.000818, -0.000461, 0.000655, 0.000863, 0.000587, 0.001020, -0.000410, 0.000896, -0.000378, 0.000312, 0.000777, -0.000580, 0.000255, 0.000298, 0.001228, 0.001243, -0.000119, 0.001623, 0.000833, 0.000182, 0.001150, -0.000678, 0.000212, 0.000399, -0.000308, 0.001147, 0.000226, 0.001643, 0.000859, -0.000054, 0.000720, 0.000115, 0.000926, -0.000003, -0.000316, 0.000780, 0.000979, 0.000614, -0.000158, 0.000533, -0.000286, 0.000260, -0.000190, -0.000348, -0.000900, 0.000382, 0.000304, 0.000216, 0.000245, 0.000170, 0.000520, 0.000499, 0.000714, 0.000911, 0.000052, 0.000474, 0.000117, 0.000127, 0.000628, 0.000012, -0.000483, 0.000666, -0.000109, 0.000641, 0.000365, 0.000166, -0.000361, -0.000516, 0.000355, 0.000007, 0.000118, 0.000169, -0.000125, -0.000555, 0.000164, 0.000535, 0.000436, 0.000699, 0.000365, -0.000307, -0.000343, 0.000680, -0.000764, 0.000068, -0.000194, 0.000347, 0.000381, 0.000433, 0.000546, 0.000447, -0.000484, -0.000064, -0.000406, -0.000368, -0.000161, -0.000228, 0.000287, -0.000156, -0.000444, 0.000551, -0.000430, -0.000768, -0.000460, 0.000325, 0.000076, -0.000207, 0.000015, -0.000113, -0.000742, 0.000047, 0.000084, -0.000183, 0.000091, -0.000570, -0.000710, -0.000366, -0.000210, -0.000374, -0.000287, -0.000692, -0.000049, -0.000135, -0.000357, 0.000121, 0.000085, 0.000372, -0.000839, -0.000680, -0.000200, -0.000527, -0.000359, -0.000086, -0.000426, 0.000597, 0.000499, -0.000640, 0.000328, -0.000128, 0.000107, -0.000379, -0.001162, 0.000061, -0.000357, -0.000403, -0.000524, 0.000019, -0.000081, -0.000413, 0.000216, -0.000250, 0.000059, -0.000131, 0.000082, -0.000455, -0.000276, -0.000004, -0.000120, 0.000468, -0.000002, -0.000142, -0.000122, -0.000563, -0.000125, -0.000175, -0.000159, -0.000042, -0.000283, 0.000562, -0.000229, 0.000430, -0.000011, -0.000247, -0.000463, -0.000865, 0.000130, -0.000071, -0.000067, -0.000427, -0.000256, 0.000032, -0.000145, -0.000399, 0.000381, 0.000086, -0.000321, -0.000122, -0.000397, -0.000152, -0.000074, 0.000315, 0.000093, 0.000199, -0.000205, -0.000091, -0.000162, 0.000115, -0.000168, -0.000032, -0.000016, -0.000246, 0.000182, 0.000190, -0.000324, -0.000134, 0.000512, -0.000178, -0.000424, -0.000028, -0.000085, -0.000147, 0.000197, -0.000834, 0.000030, -0.000530, -0.000258, -0.000505, 0.000061, -0.000222, -0.000487, -0.000309, -0.000388, 0.000232, 0.000085, -0.000025, 0.000409, 0.000037, 0.000024, -0.000663, -0.000107, -0.000619, -0.000708, -0.000213, -0.000624, -0.000049, 0.000239, -0.000044, 0.000003, -0.000198, 0.000154, -0.000248, -0.000249, -0.000588, -0.001104, -0.000401, -0.000593, 0.000129, 0.000136, -0.000093, -0.000126, -0.000373, -0.000359, -0.000267, -0.000401, -0.000360, -0.000283, -0.000097, -0.000031, -0.000688, 0.000186, -0.000578, 0.000036, -0.000596, -0.000221, -0.000672, -0.001273, -0.000191, -0.000500, 0.000238, -0.000796, -0.000691, 0.000024, 0.000215, 0.000847, 0.000125, 0.001043, 0.000075, 0.000085, 0.000715, 0.000065, 0.000371, -0.000519, -0.000063, 0.000177, 0.000358, 0.000051, -0.000651, 0.000138, -0.000180, 0.000637, -0.000281, -0.000376, 0.000380, 0.000537, 0.000218, 0.000735, 0.000697, -0.000098, 0.000416, 0.000142, 0.000709, 0.000191, 0.000341, 0.000105, 0.000067, 0.000278, 0.001575, 0.000926, 0.000596, 0.000272, 0.000405, 0.000964, 0.000357, 0.000330, 0.000415, -0.000082, 0.000366, -0.000020, 0.000787, -0.000281, -0.000360, 0.000563, 0.000434, 0.000527, -0.000412, 0.000033, 0.000479, -0.000768, 0.000787, 0.000912, 0.000505, 0.000089, 0.000144, 0.000870, -0.000662, 0.000115, 0.000262, 0.000097, 0.000070, 0.000118, 0.000339, -0.000126, -0.000129, 0.000136, 0.000564, 0.000566, -0.000274, 0.000202, -0.000295, -0.000050, 0.000090, 0.000559, 0.000066, -0.000200, -0.000170, -0.000164, -0.000138, -0.000332, -0.000197, -0.000403, -0.000042, 0.000333, 0.000258, -0.000564, -0.000785, -0.000299, -0.000082, -0.000354, -0.000056, 0.000008, 0.000291, 0.000100, -0.000074, -0.000110, -0.000463, 0.000180, 0.000172, 0.000290, -0.000335, -0.000065, -0.000535, -0.000094, -0.000085, 0.000295, -0.000204, -0.000244, -0.000151, -0.000368, 0.000153, -0.000100, 0.000507, -0.000058, -0.000180, 0.000084, -0.000080, -0.000273, -0.000345, 0.000531, 0.000169, 0.000021, -0.000148, 0.000038, -0.000141, 0.000197, -0.000021, 0.000175, 0.000139, -0.000493, 0.000075, -0.000736, -0.000113, -0.000241, -0.000381, -0.000390, -0.000083, 0.000087, -0.000059, -0.000080, -0.000318, -0.000238, -0.000150, 0.000276, 0.000315, -0.000005, -0.000239, 0.000105, -0.000309, 0.000068, -0.000252, 0.000036, -0.000005, -0.000532, -0.000201, -0.000104, 0.000049, 0.000019, 0.000150, 0.000255, -0.000171, -0.000313, -0.000192, 0.000210, -0.000260, 0.000202, -0.000168, -0.000124, -0.000234, -0.000346, -0.000190, 0.000169, 0.000057, -0.000356, 0.000151, -0.000042, 0.000170, -0.000499, 0.000010, -0.000090, -0.000198, -0.000048, -0.000364, -0.000029, -0.000463, -0.000067, -0.000141, -0.000298, -0.000217, -0.000109, -0.000315, -0.000099, 0.000055, -0.000111, 0.000021, -0.000206, -0.000319, -0.000209, -0.000037, -0.000524, -0.000117, -0.000303, 0.000054, -0.000232, -0.000058, -0.000460, -0.000194, -0.000260, -0.000107, -0.000266, 0.000074, -0.000296, -0.000002, 0.000145, -0.000211, -0.000160, 0.000163, 0.000451, -0.000373, -0.000004, -0.000032, 0.000057, 0.000304, -0.000061, 0.000058, -0.000154, -0.000382, -0.000271, 0.000002, -0.000048, -0.000048, -0.000044, -0.000134, -0.000412, -0.000038, -0.000501, -0.000026, -0.000409, 0.000225, -0.000182, 0.000259, -0.000105, -0.000395, -0.000488, -0.000119, 0.000176, -0.000255, -0.000199, -0.000132, -0.000048, -0.000131, -0.000305, -0.000283, 0.000191, -0.000519, -0.000381, -0.000343, 0.000187, -0.000422, -0.000436, 0.000186, -0.000276, -0.000384, -0.000336, -0.000208, -0.000485, -0.000812, -0.000018, -0.000534, 0.000162, -0.000851, 0.000041, -0.000366, -0.000256, -0.000061, 0.000339, 0.000041, 0.000082, 0.000208, -0.000236, -0.000178, 0.000273, 0.000076, 0.000072, 0.000210, 0.000451, 0.000042, -0.000307, 0.000283, 0.000080, 0.000495, -0.000437, 0.000411, 0.000129, 0.000045, -0.000617, 0.000309, 0.000511, 0.000079, 0.000888, 0.000309, 0.000247, -0.000108, 0.000374, -0.000236, 0.000390, 0.000322, 0.000470, 0.000844, 0.000013, 0.000855, 0.000259, 0.000859, -0.000019, 0.000238, 0.000373, 0.000269, -0.000080, 0.000027, 0.000198, 0.000241, 0.000517, -0.000299, 0.000523, -0.000083, 0.000013, -0.000028, 0.000538, 0.000265, -0.000001, 0.000437, 0.000292, 0.000248, 0.000121, 0.000462, 0.000321, 0.000212, -0.000131, 0.000020, 0.000033, 0.000092, -0.000118, 0.000306, 0.000187, -0.000116, 0.000094, 0.000264, 0.000194, 0.000133, 0.000141, -0.000201, 0.000088, -0.000363, -0.000073, -0.000137, -0.000173, -0.000042, 0.000031, -0.000154, 0.000023, -0.000112, 0.000077, 0.000255, 0.000368, 0.000104, -0.000051, -0.000022, -0.000133, -0.000118, 0.000350, 0.000387, -0.000114, -0.000131, 0.000094, -0.000103, 0.000004, -0.000129, -0.000064, 0.000095, 0.000151, -0.000128, -0.000020, -0.000037, -0.000096, 0.000078, 0.000062, -0.000170, -0.000201, 0.000136, -0.000162, -0.000129, 0.000168, 0.000594, 0.000205, -0.000047, -0.000131, 0.000319, -0.000077, 0.000008, 0.000165, -0.000082, -0.000439, 0.000025, -0.000120, -0.000347, -0.000388, -0.000250, -0.000240, -0.000340, 0.000103, -0.000125, 0.000111, -0.000155, 0.000048, -0.000506, -0.000131, -0.000297, -0.000298, -0.000021, -0.000339, -0.000062, -0.000139, 0.000206, -0.000329, 0.000098, -0.000035, -0.000117, 0.000057, -0.000224, -0.000132, -0.000225, 0.000252, -0.000037, -0.000073, -0.000286, -0.000233, -0.000433, -0.000130, -0.000066, -0.000136, 0.000059, -0.000590, 0.000041, -0.000148, -0.000107, -0.000215, 0.000085, 0.000276, 0.000110, -0.000087, -0.000248, -0.000331, -0.000306, -0.000041, -0.000134, -0.000138, 0.000043, 0.000215, -0.000056, -0.000221, 0.000099, 0.000011, -0.000227, -0.000264, -0.000140, -0.000191, 0.000108, -0.000182, -0.000000, -0.000001, -0.000397, -0.000206, -0.000230, 0.000008, -0.000291, -0.000124, -0.000127, -0.000219, 0.000055, 0.000182, -0.000166, -0.000049, 0.000210, -0.000125, -0.000007, 0.000044, 0.000140, -0.000168, -0.000138, -0.000120, -0.000405, -0.000150, -0.000015, -0.000082, -0.000174, 0.000038, -0.000152, -0.000209, -0.000107, -0.000130, -0.000454, -0.000084, 0.000398, 0.000116, -0.000072, 0.000122, 0.000042, -0.000051, 0.000092, -0.000279, -0.000021, -0.000289, -0.000398, -0.000296, -0.000035, -0.000175, -0.000077, -0.000387, -0.000216, -0.000330, -0.000223, -0.000206, -0.000480, -0.000223, -0.000063, -0.000067, -0.000204, -0.000174, 0.000244, -0.000529, -0.000007, -0.000144, -0.000045, -0.000335, -0.000357, -0.000777, -0.000006, -0.000484, -0.000224, -0.000383, -0.000082, -0.000360, -0.000209, -0.000299, -0.000537, -0.000440, -0.000090, -0.000219, -0.000259, -0.000197, 0.000133, -0.000019, -0.000367, 0.000077, 0.000643, -0.000194, 0.000328, 0.000225, 0.000332, 0.000187, 0.000177, 0.000036, 0.000010, -0.000200, 0.000361, 0.000405, -0.000139, 0.000118, 0.000450, 0.000227, 0.000079, 0.000527, 0.000274, 0.000192, 0.000114, 0.000144, 0.000531, 0.000765, -0.000040, 0.000621, 0.000284, 0.000389, 0.000019, 0.000541, -0.000119, 0.000613, 0.000334, 0.001127, 0.000264, 0.000530, 0.000131, 0.000494, -0.000266, 0.000707, 0.000290, -0.000298, -0.000425, 0.000797, 0.000017, 0.000336, 0.000096, 0.000044, -0.000163, 0.000369, -0.000284, 0.000399, 0.000109, 0.000056, -0.000091, 0.000176, 0.000172, 0.000199, -0.000074, 0.000127, 0.000160, 0.000311, 0.000181, -0.000385, -0.000385, 0.000313, 0.000208, 0.000059, 0.000298, 0.000110, -0.000326, 0.000200, -0.000114, 0.000119, -0.000132, -0.000343, -0.000179, 0.000072, -0.000216, -0.000035, -0.000051, -0.000010, 0.000014, 0.000228, -0.000077, -0.000185, -0.000183, -0.000150, -0.000103, 0.000204, 0.000229, 0.000109, 0.000010, -0.000088, -0.000087, 0.000011, -0.000010, -0.000408, -0.000268, 0.000323, -0.000372, -0.000463, -0.000073, 0.000042, -0.000285, 0.000433, 0.000018, -0.000160, -0.000123, 0.000005, -0.000194, 0.000152, -0.000076, 0.000032, 0.000005, 0.000045, 0.000004, -0.000033, -0.000082, 0.000000, -0.000026, -0.000225, -0.000098, -0.000238, -0.000100, 0.000079, 0.000292, 0.000100, -0.000053, 0.000011, -0.000161, -0.000261, 0.000043, -0.000027, 0.000067, -0.000157, 0.000038, -0.000070, -0.000015, 0.000059, -0.000016, -0.000080, -0.000070, -0.000169, -0.000092, -0.000211, -0.000241, -0.000075, 0.000189, -0.000064, -0.000171, 0.000026, -0.000327, 0.000014, 0.000182, -0.000092, -0.000098, -0.000014, -0.000145, 0.000116, -0.000125, -0.000002, 0.000002, -0.000027, -0.000060, -0.000264, -0.000091, -0.000185, -0.000067, -0.000107, -0.000087, 0.000160, 0.000021, -0.000097, -0.000088, -0.000072, 0.000011, -0.000242, -0.000099, -0.000107, -0.000082, -0.000020, -0.000015, 0.000037, -0.000008, -0.000344, -0.000066, -0.000207, -0.000219, -0.000149, 0.000113, -0.000302, 0.000026, -0.000110, -0.000142, 0.000099, 0.000246, -0.000273, 0.000183, -0.000136, -0.000143, -0.000079, -0.000205, -0.000010, 0.000070, 0.000110, -0.000015, -0.000050, -0.000091, -0.000025, -0.000047, -0.000048, -0.000321, -0.000020, 0.000005, 0.000005, -0.000164, 0.000274, 0.000004, -0.000160, -0.000088, -0.000210, -0.000163, -0.000231, -0.000226, 0.000088, -0.000005, -0.000036, 0.000143, 0.000128, -0.000048, 0.000044, -0.000148, -0.000152, -0.000025, -0.000172, -0.000097, -0.000013, -0.000513, -0.000186, -0.000183, -0.000088, 0.000047, -0.000114, -0.000259, -0.000411, -0.000125, -0.000378, -0.000387, -0.000030, -0.000021, -0.000285, 0.000016, -0.000170, -0.000377, -0.000073, -0.000019, -0.000299, 0.000192, -0.000115, -0.000061, 0.000068, -0.000094, -0.000077, -0.000348, -0.000489, -0.000151, -0.000183, -0.000157, -0.000356, -0.000194, -0.000170, -0.000256, -0.000217, -0.000063, -0.000163, -0.000189, -0.000076, -0.000030, -0.000035, 0.000386, -0.000034, 0.000320, 0.000321, 0.000015, -0.000076, 0.000004, 0.000239, 0.000029, 0.000240, 0.000169, -0.000181, -0.000086, 0.000265, -0.000080, 0.000054, 0.000241, -0.000030, 0.000004, 0.000208, -0.000342, -0.000005, 0.000101, 0.000204, 0.000288, 0.000515, 0.000147, 0.000507, 0.000122, 0.000277, 0.000005, 0.000529, 0.000400, 0.000219, 0.000125, 0.000203, 0.000061, 0.000148, -0.000014, 0.000297, 0.000147, 0.000253, -0.000087, 0.000151, -0.000178, 0.000184, 0.000142, 0.000307, 0.000162, 0.000159, 0.000075, 0.000312, -0.000011, 0.000271, 0.000167, 0.000091, 0.000010, 0.000136, 0.000212, 0.000265, -0.000021, 0.000050, 0.000176, 0.000182, 0.000096, 0.000020, -0.000263, 0.000261, 0.000311, 0.000081, 0.000017, 0.000015, 0.000152, -0.000162, -0.000102, 0.000028, 0.000026, 0.000107, 0.000003, 0.000148, 0.000208, -0.000041, -0.000067, 0.000122, 0.000121, 0.000049, 0.000182, -0.000198, -0.000074, 0.000011, 0.000063, 0.000023, -0.000193, 0.000031, 0.000168, -0.000198, -0.000118, -0.000217, -0.000008, 0.000122, -0.000082, 0.000015, 0.000082, -0.000147, -0.000085, 0.000068, -0.000069, 0.000062, -0.000048, -0.000199, -0.000018, -0.000144, 0.000187, -0.000149, -0.000025, -0.000028, -0.000170, -0.000169, 0.000080, -0.000087, -0.000123, 0.000014, -0.000143, -0.000037, 0.000045, 0.000070, -0.000014, -0.000166, 0.000111, -0.000084, -0.000021, -0.000041, -0.000105, 0.000184, -0.000026, -0.000223, -0.000046, -0.000074, -0.000069, -0.000088, -0.000113, -0.000215, -0.000211, -0.000211, 0.000109, -0.000075, -0.000017, -0.000004, -0.000218, -0.000052, 0.000090, -0.000131, -0.000082, 0.000043, 0.000086, -0.000092, -0.000108, -0.000217, -0.000153, -0.000021, 0.000056, -0.000171, -0.000019, -0.000095, -0.000110, -0.000093, -0.000017, -0.000084, 0.000047, -0.000161, 0.000052, -0.000128, -0.000101, -0.000069, -0.000117, -0.000083, -0.000004, 0.000010, -0.000179, 0.000035, -0.000012, -0.000012, 0.000109, -0.000025, -0.000088, -0.000048, -0.000030, -0.000160, -0.000095, -0.000160, 0.000071, 0.000145, -0.000004, -0.000141, -0.000076, -0.000150, -0.000123, 0.000027, -0.000010, -0.000084, -0.000117, -0.000200, -0.000109, 0.000003, 0.000036, -0.000010, -0.000085, -0.000062, -0.000103, -0.000066, -0.000078, 0.000147, 0.000064, -0.000083, -0.000037, -0.000092, -0.000088, -0.000069, -0.000080, -0.000027, -0.000094, -0.000112, -0.000003, -0.000079, 0.000075, 0.000057, -0.000050, -0.000072, -0.000008, -0.000018, -0.000086, -0.000058, -0.000177, -0.000012, -0.000062, -0.000027, 0.000016, -0.000034, -0.000303, -0.000274, -0.000282, -0.000022, -0.000050, -0.000235, -0.000076, -0.000019, -0.000196, -0.000356, -0.000108, -0.000167, -0.000057, -0.000253, -0.000127, -0.000259, -0.000131, 0.000098, -0.000123, -0.000149, -0.000118, -0.000193, -0.000105, -0.000016, -0.000126, -0.000153, -0.000261, -0.000242, -0.000121, -0.000139, -0.000340, -0.000366, -0.000116, -0.000102, -0.000300, -0.000182, 0.000087, 0.000016, -0.000145, -0.000055, -0.000170, 0.000057, 0.000320, 0.000020, 0.000262, 0.000254, 0.000022, -0.000020, 0.000128, 0.000019, 0.000177, 0.000262, -0.000091, 0.000145, 0.000040, -0.000036, -0.000103, 0.000102, 0.000067, 0.000022, -0.000065, 0.000092, 0.000103, -0.000131, 0.000297, 0.000129, 0.000256, 0.000438, 0.000241, 0.000298, 0.000382, 0.000267, 0.000177, 0.000345, 0.000053, 0.000104, 0.000303, 0.000077, -0.000097, 0.000296, 0.000329, 0.000014, 0.000260, -0.000140, 0.000228, 0.000080, -0.000009, 0.000123, 0.000023, 0.000391, 0.000240, 0.000363, 0.000279, 0.000284, 0.000068, 0.000128, 0.000186, 0.000062, 0.000150, 0.000158, -0.000050, 0.000170, 0.000262, -0.000151, 0.000107, -0.000093, 0.000039, 0.000140, -0.000172, 0.000206, 0.000128, -0.000088, -0.000148, -0.000037, 0.000140, 0.000025, 0.000006, -0.000005, 0.000056, -0.000136, 0.000060, 0.000011, 0.000039, -0.000016, 0.000107, -0.000135, 0.000047, 0.000050, -0.000076, 0.000099, -0.000191, 0.000118, 0.000087, -0.000242, -0.000042, -0.000002, -0.000103, 0.000017, -0.000131, -0.000028, -0.000075, -0.000176, 0.000163, 0.000044, 0.000149, -0.000059, -0.000030, 0.000044, -0.000100, -0.000090, 0.000064, -0.000030, -0.000212, -0.000006, -0.000129, -0.000220, 0.000054, 0.000010, -0.000149, -0.000023, -0.000030, -0.000006, -0.000122, -0.000046, 0.000086, 0.000044, 0.000089, 0.000159, -0.000052, -0.000039, -0.000139, -0.000066, -0.000052, -0.000198, -0.000138, -0.000070, -0.000165, 0.000006, -0.000126, -0.000106, -0.000097, 0.000042, -0.000114, -0.000094, 0.000011, 0.000016, -0.000076, -0.000145, 0.000071, 0.000092, 0.000015, -0.000048, 0.000098, 0.000043, -0.000139, -0.000099, -0.000144, -0.000150, -0.000076, 0.000022, -0.000142, 0.000021, -0.000075, -0.000209, -0.000128, -0.000071, 0.000081, -0.000006, -0.000066, -0.000120, -0.000083, -0.000178, -0.000022, -0.000063, 0.000114, 0.000055, 0.000012, -0.000011, -0.000125, -0.000120, -0.000072, 0.000047, 0.000026, 0.000067, -0.000079, -0.000160, -0.000102, -0.000041, -0.000045, -0.000126, -0.000073, -0.000105, -0.000077, -0.000130, -0.000051, -0.000027, 0.000036, -0.000003, 0.000023, -0.000025, -0.000085, -0.000000, -0.000074, -0.000036, 0.000045, 0.000077, -0.000024, -0.000109, -0.000163, 0.000033, -0.000191, -0.000040, -0.000042, -0.000051, -0.000117, 0.000012, -0.000015, -0.000117, -0.000049, -0.000131, -0.000112, -0.000080, 0.000064, -0.000058, -0.000123, 0.000076, -0.000089, 0.000015, 0.000111, -0.000020, 0.000157, -0.000006, -0.000099, -0.000046, -0.000073, -0.000107, -0.000012, 0.000102, -0.000016, -0.000148, -0.000239, -0.000214, -0.000172, -0.000099, -0.000060, -0.000060, -0.000016, 0.000010, -0.000150, -0.000018, -0.000070, -0.000077, 0.000004, -0.000071, -0.000175, -0.000047, -0.000122, -0.000119, 0.000026, -0.000085, -0.000057, -0.000124, -0.000218, -0.000186, -0.000210, -0.000218, -0.000078, -0.000122, -0.000128, 0.000015, -0.000197, -0.000059, -0.000145, -0.000137, -0.000025, -0.000000, -0.000160, -0.000056, -0.000101, -0.000107, -0.000125, 0.000121, 0.000314, 0.000136, 0.000192, 0.000061, -0.000001, 0.000001, -0.000058, 0.000079, 0.000356, 0.000009, 0.000068, 0.000258, -0.000037, 0.000013, 0.000195, -0.000033, 0.000166, -0.000013, -0.000045, 0.000009, -0.000016, -0.000015, 0.000128, 0.000279, 0.000405, 0.000204, 0.000105, 0.000046, 0.000134, 0.000076, 0.000257, 0.000196, 0.000150, 0.000032, 0.000052, 0.000143, 0.000230, 0.000060, -0.000013, 0.000049, 0.000191, 0.000079, -0.000051, 0.000002, -0.000087, 0.000048, 0.000154, 0.000196, 0.000099, 0.000067, 0.000019, 0.000176, 0.000053, 0.000017, 0.000142, 0.000073, -0.000037, 0.000057, 0.000092, 0.000119, -0.000056, -0.000034, 0.000034, 0.000040, 0.000061, -0.000085, -0.000086, -0.000103, 0.000126, 0.000142, 0.000118, -0.000074, -0.000041, 0.000185, -0.000013, 0.000018, -0.000073, -0.000049, 0.000196, 0.000112, 0.000046, 0.000111, 0.000057, 0.000009, -0.000045, 0.000041, -0.000041, 0.000018, -0.000131, 0.000042, -0.000008, 0.000022, 0.000020, -0.000078, 0.000013, 0.000019, 0.000075, 0.000031, -0.000177, -0.000093, -0.000025, 0.000034, 0.000033, -0.000086, 0.000054, -0.000026, -0.000032, 0.000065, 0.000037, 0.000018, 0.000086, -0.000006, -0.000033, 0.000033, -0.000061, -0.000065, 0.000098, 0.000132, 0.000118, 0.000024, -0.000185, -0.000177, -0.000116, -0.000083, -0.000013, -0.000079, -0.000111, 0.000101, -0.000009, -0.000159, -0.000151, -0.000097, -0.000014, -0.000022, 0.000009, 0.000035, -0.000019, -0.000137, -0.000062, 0.000055, 0.000097, -0.000043, -0.000004, -0.000056, -0.000081, -0.000043, -0.000042, -0.000102, -0.000114, 0.000018, -0.000011, -0.000159, -0.000109, -0.000103, 0.000023, 0.000133, 0.000024, 0.000039, -0.000083, -0.000233, -0.000049, -0.000048, -0.000003, 0.000074, 0.000049, -0.000048, -0.000052, -0.000032, -0.000204, -0.000023, 0.000034, -0.000075, -0.000090, -0.000077, -0.000125, -0.000054, -0.000053, -0.000004, 0.000023, -0.000145, -0.000091, -0.000039, -0.000053, -0.000024, 0.000124, 0.000008, -0.000064, 0.000012, -0.000005, -0.000108, -0.000054, 0.000018, 0.000018, 0.000034, -0.000097, -0.000077, -0.000100, -0.000067, -0.000036, 0.000026, -0.000033, 0.000029, -0.000010, -0.000062, -0.000080, -0.000037, -0.000064, -0.000041, -0.000068, -0.000003, -0.000001, -0.000089, -0.000041, -0.000006, -0.000027, 0.000031, -0.000044, -0.000015, -0.000115, -0.000022, 0.000001, 0.000052, 0.000090, 0.000003, 0.000032, -0.000164, -0.000105, -0.000060, -0.000097, 0.000017, 0.000014, -0.000057, 0.000035, -0.000098, -0.000114, -0.000082, 0.000025, -0.000095, -0.000201, -0.000073, -0.000091, -0.000061, -0.000055, -0.000056, -0.000005, -0.000223, -0.000141, -0.000087, -0.000151, -0.000072, 0.000023, -0.000153, -0.000138, -0.000050, -0.000012, -0.000093, -0.000090, 0.000017, -0.000052, -0.000105, -0.000190, -0.000215, -0.000124, -0.000203, -0.000087, -0.000142, -0.000013, -0.000088, -0.000242, -0.000079, -0.000080, -0.000149, -0.000129, -0.000124, -0.000101, -0.000076, -0.000192, 0.000049, -0.000126, 0.000127, 0.000134, 0.000084, 0.000039, -0.000011, 0.000004, -0.000002, -0.000034, 0.000266, 0.000253, 0.000070, 0.000041, -0.000025, 0.000072, 0.000180, -0.000011, 0.000160, 0.000028, 0.000156, -0.000050, -0.000019, 0.000190, 0.000063, 0.000386, 0.000263, 0.000167, 0.000233, 0.000094, 0.000144, -0.000021, 0.000142, 0.000347, -0.000035, -0.000006, 0.000221, -0.000026, 0.000201, 0.000005, 0.000114, 0.000165, 0.000011, 0.000177, -0.000035, 0.000101, 0.000102, 0.000189, 0.000215, 0.000087, 0.000107, 0.000166, 0.000054, 0.000023, -0.000020, 0.000228, 0.000135, -0.000098, 0.000108, 0.000012, 0.000036, 0.000007, 0.000073, 0.000277, -0.000027, 0.000084, -0.000064, 0.000070, 0.000031, 0.000072, 0.000061, -0.000064, -0.000028, 0.000036, -0.000050, 0.000019, -0.000000, 0.000117, 0.000119, -0.000015, -0.000053, 0.000029, -0.000105, -0.000006, -0.000064, 0.000144, -0.000032, 0.000010, -0.000042, -0.000009, 0.000026, -0.000009, 0.000014, 0.000002, -0.000070, -0.000000, -0.000094, -0.000085, -0.000058, 0.000033, 0.000045, -0.000055, 0.000024, 0.000004, -0.000197, -0.000047, -0.000054, 0.000038, 0.000041, -0.000087, -0.000078, -0.000070, -0.000062, 0.000091, 0.000114, 0.000039, -0.000046, -0.000078, -0.000098, -0.000033, -0.000110, -0.000012, -0.000106, -0.000063, 0.000025, -0.000083, -0.000092, -0.000123, -0.000063, 0.000070, -0.000016, -0.000051, -0.000044, 0.000008, -0.000012, 0.000095, 0.000071, 0.000095, -0.000024, -0.000081, -0.000052, -0.000088, -0.000094, -0.000041, -0.000105, 0.000016, 0.000081, -0.000106, -0.000026, -0.000030, -0.000010, 0.000043, -0.000024, -0.000060, -0.000106, -0.000021, 0.000046, 0.000131, -0.000004, 0.000085, -0.000017, -0.000052, -0.000058, 0.000008, 0.000009, -0.000028, -0.000050, -0.000060, -0.000000, -0.000067, -0.000008, 0.000046, -0.000047, -0.000094, -0.000057, -0.000145, -0.000093, 0.000001, 0.000019, 0.000014, -0.000088, -0.000040, -0.000100, -0.000106, -0.000062, -0.000024, -0.000027, -0.000063, -0.000042, -0.000115, 0.000000, -0.000042, 0.000028, 0.000035, -0.000046, -0.000078, -0.000075, -0.000076, -0.000054, -0.000003, -0.000015, -0.000066, -0.000098, -0.000021, -0.000080, 0.000030, -0.000055, 0.000037, -0.000075, -0.000020, -0.000061, -0.000053, -0.000044, 0.000068, 0.000057, -0.000056, -0.000013, -0.000119, -0.000017, -0.000044, 0.000002, 0.000031, -0.000028, -0.000046, 0.000040, -0.000075, -0.000019, -0.000018, -0.000011, 0.000005, -0.000020, -0.000028, -0.000110, -0.000035, -0.000037, 0.000056, 0.000006, 0.000001, -0.000057, -0.000012, -0.000043, -0.000146, -0.000043, -0.000082, -0.000036, -0.000029, -0.000039, -0.000129, -0.000168, -0.000003, -0.000065, -0.000133, 0.000011, -0.000118, -0.000069, -0.000039, -0.000005, -0.000018, -0.000126, -0.000025, -0.000087, -0.000066, -0.000101, -0.000094, -0.000044, -0.000118, -0.000047, -0.000030, -0.000075, -0.000089, -0.000133, -0.000037, -0.000055, -0.000100, -0.000029, -0.000113, -0.000078, -0.000061, -0.000064, -0.000016, -0.000026, 0.000014, 0.000060, 0.000082, 0.000108, 0.000100, -0.000042, 0.000018, 0.000002, 0.000077, 0.000182, -0.000051, 0.000011, 0.000122, -0.000034, 0.000109, 0.000007, 0.000045, 0.000049, 0.000009, 0.000090, 0.000056, 0.000017, 0.000034, 0.000153, 0.000213, 0.000149, 0.000123, 0.000078, 0.000070, 0.000027, 0.000151, 0.000189, -0.000036, 0.000107, 0.000190, 0.000083, 0.000082, 0.000107, 0.000113, 0.000050, 0.000064, 0.000095, -0.000005, 0.000113, 0.000091, 0.000129, 0.000043, 0.000129, 0.000077, -0.000003, 0.000035, 0.000057, 0.000110, 0.000057, -0.000113, 0.000003, 0.000112, -0.000030, 0.000055, 0.000028, 0.000034, -0.000029, 0.000036, 0.000103, -0.000010, 0.000028, -0.000025, 0.000030, 0.000059, 0.000087, 0.000047, -0.000046, 0.000130, 0.000100, 0.000047, -0.000020, -0.000104, 0.000032, 0.000024, 0.000011, 0.000089, -0.000013, -0.000080, -0.000021, 0.000017, 0.000028, 0.000024, -0.000039, -0.000059, -0.000067, 0.000026, -0.000027, -0.000067, 0.000041, 0.000067, 0.000031, -0.000000, -0.000084, -0.000010, -0.000048, -0.000043, 0.000046, 0.000022, -0.000059, 0.000034, -0.000044, 0.000008, 0.000006, 0.000014, 0.000019, -0.000052, -0.000016, -0.000061, -0.000035, -0.000038, 0.000007, -0.000026, -0.000011, -0.000040, -0.000049, -0.000026, -0.000020, 0.000025, -0.000007, -0.000023, -0.000049, -0.000033, -0.000020, -0.000026, 0.000067, 0.000012, -0.000086, -0.000063, -0.000030, 0.000036, -0.000056, 0.000021, -0.000082, -0.000007, 0.000012, -0.000094, -0.000004, -0.000049, -0.000034, -0.000018, -0.000049, -0.000003, -0.000001, -0.000001, -0.000033, -0.000013, -0.000033, -0.000064, -0.000101, 0.000006, 0.000024, -0.000049, -0.000049, -0.000109, -0.000021, -0.000052, -0.000077, -0.000063, -0.000052, -0.000025, -0.000021, -0.000076, 0.000000, -0.000020, -0.000004, -0.000017, 0.000050, 0.000009, -0.000005, -0.000040, -0.000051, 0.000024, -0.000055, -0.000020, -0.000058, 0.000008, 0.000019, -0.000043, 0.000015, -0.000010, -0.000085, 0.000007, -0.000042, -0.000019, 0.000016, -0.000020, -0.000031, -0.000035, -0.000032, -0.000018, -0.000035, -0.000023, -0.000022, -0.000037, -0.000039, -0.000016, -0.000066, 0.000015, 0.000025, -0.000029, -0.000016, 0.000001, -0.000019, -0.000029, -0.000016, 0.000044, 0.000036, 0.000029, 0.000016, -0.000029, 0.000018, 0.000029, 0.000012, -0.000016, 0.000001, -0.000019, -0.000033, -0.000032, 0.000009, -0.000011, -0.000023, -0.000035, -0.000045, -0.000045, -0.000036, -0.000012, -0.000012, 0.000005, -0.000003, -0.000007, -0.000095, 0.000027, -0.000048, -0.000036, -0.000092, -0.000068, -0.000074, -0.000121, -0.000067, -0.000068, -0.000036, -0.000149, -0.000061, -0.000032, -0.000051, -0.000016, -0.000079, -0.000079, -0.000021, -0.000050, -0.000039, -0.000086, -0.000013, -0.000082, -0.000076, -0.000074, -0.000090, -0.000165, -0.000104, -0.000083, -0.000010, -0.000079, -0.000056, -0.000108, -0.000098, -0.000053, -0.000073, -0.000062, -0.000016, -0.000128, -0.000035, -0.000075, -0.000002, -0.000108, 0.000033, -0.000103, 0.000073, 0.000025, 0.000019, 0.000061, 0.000046, 0.000056, 0.000002, 0.000138, 0.000073, -0.000019, 0.000137, -0.000081, 0.000122, 0.000020, 0.000053, 0.000000, 0.000046, -0.000008, 0.000009, 0.000091, 0.000034, 0.000042, 0.000079, 0.000212, 0.000069, 0.000188, 0.000023, 0.000023, 0.000131, 0.000025, 0.000224, 0.000091, 0.000152, 0.000033, 0.000057, 0.000145, 0.000024, 0.000120, 0.000038, 0.000092, -0.000026, 0.000095, 0.000040, 0.000157, 0.000083, 0.000060, 0.000062, 0.000061, 0.000105, -0.000008, 0.000100, 0.000015, 0.000040, 0.000091, -0.000038, 0.000100, -0.000066, 0.000068, -0.000010, 0.000043, -0.000075, -0.000009, 0.000021, 0.000114, 0.000013, -0.000003, 0.000049, 0.000018, 0.000021, 0.000013, 0.000062, 0.000044, -0.000048, 0.000037, -0.000037, 0.000043, -0.000014, 0.000055, -0.000011, -0.000018, -0.000037, 0.000070, -0.000016, 0.000033, -0.000054, -0.000032, 0.000023, 0.000009, 0.000011, -0.000014, 0.000052, -0.000030, -0.000042, 0.000019, -0.000034, 0.000045, -0.000030, 0.000050, -0.000006, -0.000022, -0.000009, 0.000073, 0.000012, 0.000015, -0.000007, -0.000007, -0.000006, -0.000049, -0.000043, -0.000065, 0.000039, -0.000090, 0.000012, -0.000008, -0.000021, 0.000016, 0.000010, 0.000001, -0.000024, -0.000044, -0.000025, -0.000030, 0.000058, 0.000043, -0.000019, -0.000020, -0.000056, -0.000029, -0.000050, -0.000013, -0.000015, -0.000075, -0.000007, 0.000011, -0.000021, -0.000029, -0.000024, -0.000033, -0.000031, -0.000005, -0.000024, -0.000053, -0.000008, -0.000058, -0.000009, -0.000045, 0.000027, -0.000024, -0.000034, -0.000045, -0.000020, -0.000026, -0.000051, -0.000013, -0.000061, -0.000005, -0.000011, -0.000079, -0.000058, -0.000039, -0.000055, -0.000002, 0.000012, -0.000034, -0.000006, -0.000032, -0.000023, -0.000022, 0.000046, -0.000008, -0.000029, -0.000026, 0.000002, -0.000011, 0.000026, -0.000019, -0.000050, -0.000032, -0.000039, 0.000003, -0.000030, -0.000053, -0.000054, -0.000025, -0.000018, -0.000010, -0.000003, 0.000001, -0.000007, -0.000015, -0.000039, -0.000018, 0.000028, 0.000003, -0.000036, -0.000017, -0.000003, -0.000033, -0.000006, -0.000039, -0.000001, -0.000042, -0.000049, -0.000038, -0.000045, -0.000044, -0.000005, 0.000005, -0.000049, -0.000033, -0.000029, -0.000041, -0.000019, 0.000017, -0.000032, -0.000030, -0.000060, -0.000033, -0.000051, 0.000002, 0.000017, -0.000012, -0.000007, -0.000056, 0.000003, -0.000053, -0.000035, 0.000000, -0.000022, 0.000006, -0.000016, -0.000022, -0.000011, 0.000003, -0.000011, -0.000031, 0.000010, -0.000012, -0.000089, -0.000048, -0.000034, -0.000027, 0.000000, -0.000048, -0.000055, -0.000105, -0.000022, -0.000063, -0.000040, -0.000023, -0.000070, 0.000024, -0.000054, -0.000027, -0.000039, -0.000033, -0.000030, -0.000039, -0.000037, -0.000036, -0.000087, -0.000047, -0.000082, -0.000038, -0.000041, -0.000034, -0.000008, -0.000140, -0.000066, -0.000034, -0.000014, -0.000042, -0.000076, -0.000003, -0.000044, -0.000038, -0.000063, -0.000063, -0.000037, -0.000023, 0.000062, 0.000037, 0.000071, 0.000037, 0.000038, 0.000068, 0.000076, 0.000110, 0.000019, 0.000092, -0.000044, 0.000014, 0.000106, 0.000042, 0.000011, 0.000030, 0.000016, 0.000016, 0.000047, -0.000009, 0.000077, 0.000014, 0.000105, 0.000092, 0.000130, 0.000072, 0.000100, 0.000060, -0.000020, 0.000203, 0.000073, 0.000055, 0.000073, 0.000052, 0.000104, 0.000084, 0.000004, 0.000011, 0.000058, 0.000067, 0.000051, 0.000031, 0.000074, 0.000053, 0.000077, 0.000064, 0.000093, -0.000001, 0.000045, 0.000052, -0.000048, 0.000020, 0.000056, 0.000006, -0.000027, -0.000020, 0.000042, 0.000011, -0.000019, -0.000041, 0.000044, 0.000087, 0.000023, 0.000016, 0.000025, 0.000046, -0.000011, 0.000051, 0.000022, -0.000042, -0.000011, 0.000010, 0.000048, -0.000050, 0.000057, -0.000016, -0.000034, -0.000018, 0.000058, 0.000045, -0.000001, -0.000051, -0.000022, 0.000052, -0.000002, -0.000004, 0.000007, 0.000018, -0.000010, 0.000012, -0.000027, 0.000006, -0.000012, -0.000001, 0.000028, -0.000025, 0.000042, 0.000056, -0.000016, 0.000003, 0.000014, 0.000005, -0.000007, -0.000042, -0.000035, -0.000000, 0.000001, -0.000039, -0.000016, -0.000020, 0.000001, -0.000008, -0.000006, -0.000063, -0.000021, -0.000035, 0.000012, -0.000009, 0.000040, 0.000016, -0.000070, -0.000031, -0.000042, -0.000014, -0.000026, -0.000052, -0.000054, 0.000031, -0.000041, -0.000038, -0.000057, -0.000009, 0.000010, 0.000027, 0.000007, -0.000022, 0.000012, -0.000036, -0.000003, 0.000026, -0.000000, -0.000018, -0.000068, -0.000044, -0.000028, -0.000012, -0.000009, -0.000043, -0.000013, -0.000005, -0.000024, -0.000026, -0.000039, -0.000047, 0.000042, 0.000027, -0.000004, -0.000009, -0.000033, -0.000033, -0.000052, 0.000008, -0.000010, -0.000015, -0.000033, -0.000007, 0.000012, -0.000017, -0.000005, -0.000056, -0.000001, -0.000014, -0.000023, -0.000032, -0.000019, -0.000013, 0.000012, -0.000003, -0.000025, -0.000040, -0.000010, -0.000042, -0.000018, -0.000017, -0.000029, -0.000022, -0.000018, 0.000012, 0.000004, -0.000019, -0.000029, -0.000023, -0.000006, -0.000019, -0.000006, -0.000013, -0.000051, 0.000018, -0.000020, 0.000000, -0.000059, -0.000009, -0.000025, -0.000014, -0.000013, 0.000017, -0.000055, -0.000036, -0.000013, -0.000021, 0.000017, 0.000017, -0.000014, -0.000009, -0.000005, -0.000052, -0.000010, 0.000008, -0.000013, 0.000024, 0.000005, -0.000024, -0.000024, -0.000008, -0.000018, 0.000009, 0.000014, -0.000017, -0.000029, -0.000007, -0.000014, -0.000018, -0.000001, 0.000003, -0.000023, -0.000012, -0.000018, -0.000030, -0.000073, -0.000047, -0.000015, 0.000014, -0.000031, -0.000030, -0.000080, -0.000048, -0.000027, 0.000010, -0.000017, -0.000045, -0.000023, -0.000023, -0.000021, -0.000051, -0.000016, 0.000006, -0.000038, -0.000083, -0.000014, -0.000037, -0.000092, -0.000067, -0.000048, -0.000087, -0.000022, -0.000056, -0.000066, -0.000106, -0.000012, -0.000064, -0.000053, -0.000045, -0.000028, -0.000006, -0.000045, -0.000068, -0.000041, -0.000007, -0.000062, -0.000022, 0.000069, 0.000039, 0.000047, 0.000024, 0.000021, -0.000007, 0.000107, 0.000078, -0.000024, 0.000048, -0.000005, 0.000057, -0.000037, 0.000046, 0.000016, 0.000045, -0.000046, -0.000002, 0.000016, 0.000074, 0.000005, 0.000049, 0.000097, 0.000119, 0.000051, 0.000138, -0.000018, 0.000051, 0.000108, 0.000138, -0.000032, 0.000135, 0.000024, 0.000055, 0.000014, 0.000039, 0.000067, 0.000117, 0.000009, 0.000035, 0.000074, 0.000070, 0.000052, 0.000074, 0.000071, 0.000028, 0.000022, 0.000055, -0.000012, -0.000026, 0.000080, -0.000047, 0.000001, 0.000054, 0.000054, -0.000002, -0.000009, 0.000002, 0.000063, 0.000040, 0.000044, 0.000030, 0.000067, -0.000019, 0.000038, 0.000005, 0.000077, -0.000014, 0.000032, -0.000008, -0.000033, -0.000003, 0.000032, -0.000043, 0.000021, 0.000029, 0.000011, -0.000015, -0.000005, -0.000026, 0.000052, -0.000036, 0.000042, 0.000011, 0.000027, -0.000049, -0.000000, -0.000031, -0.000008, 0.000003, 0.000014, 0.000013, 0.000008, 0.000004, -0.000026, -0.000012, 0.000025, -0.000002, -0.000028, -0.000049, -0.000005, -0.000036, 0.000012, -0.000048, 0.000045, -0.000012, -0.000028, -0.000028, -0.000036, -0.000001, -0.000029, 0.000010, -0.000015, -0.000008, 0.000019, -0.000021, -0.000048, -0.000027, 0.000009, -0.000023, -0.000043, -0.000006, 0.000020, -0.000027, -0.000038, -0.000024, 0.000009, 0.000034, 0.000010, -0.000001, -0.000019, -0.000022, 0.000002, 0.000030, 0.000006, -0.000027, -0.000007, -0.000043, -0.000026, -0.000016, -0.000016, -0.000025, -0.000021, -0.000016, -0.000005, -0.000023, -0.000051, 0.000000, 0.000023, 0.000030, 0.000018, -0.000005, -0.000013, -0.000037, -0.000009, -0.000011, 0.000020, 0.000008, -0.000010, -0.000018, -0.000049, -0.000030, -0.000021, -0.000002, -0.000020, -0.000023, 0.000000, -0.000022, -0.000007, 0.000004, -0.000001, -0.000013, 0.000008, -0.000037, -0.000023, -0.000003, -0.000031, -0.000030, 0.000008, -0.000006, -0.000004, -0.000026, -0.000031, -0.000008, -0.000027, -0.000005, -0.000013, -0.000006, -0.000026, 0.000008, -0.000012, -0.000027, -0.000005, 0.000000, -0.000006, -0.000019, -0.000016, -0.000012, -0.000043, -0.000012, -0.000033, 0.000002, -0.000018, -0.000014, -0.000012, -0.000015, -0.000031, -0.000008, -0.000005, -0.000001, 0.000002, -0.000032, -0.000023, -0.000039, -0.000015, -0.000021, -0.000003, -0.000022, 0.000002, -0.000013, -0.000024, -0.000023, -0.000016, 0.000002, -0.000001, -0.000028, -0.000017, 0.000006, -0.000042, 0.000004, -0.000005, 0.000011, -0.000012, -0.000018, -0.000034, -0.000016, 0.000005, -0.000015, -0.000021, -0.000019, -0.000032, -0.000024, -0.000045, -0.000037, -0.000068, -0.000002, -0.000059, -0.000015, -0.000044, -0.000039, -0.000045, -0.000019, -0.000037, -0.000005, -0.000022, -0.000036, -0.000037, -0.000032, -0.000039, -0.000015, -0.000017, -0.000022, -0.000073, -0.000021, -0.000042, -0.000051, -0.000049, -0.000040, -0.000043, 0.000005, -0.000032, -0.000017, -0.000025, -0.000046, -0.000017, 0.000000, -0.000026, -0.000025, -0.000027, -0.000036, 0.000044, 0.000055, 0.000043, 0.000014, 0.000016, 0.000018, 0.000097, -0.000002, 0.000061, -0.000011, 0.000033, 0.000027, 0.000051, 0.000021, 0.000008, -0.000012, 0.000025, 0.000018, 0.000046, 0.000021, 0.000038, 0.000021, 0.000076, 0.000053, 0.000101, 0.000016, 0.000061, 0.000035, 0.000084, 0.000039, 0.000092, 0.000036, 0.000022, 0.000066, 0.000041, 0.000063, 0.000015, 0.000022, 0.000031, 0.000028, 0.000053, 0.000064, 0.000060, 0.000023, 0.000019, 0.000013, 0.000031, 0.000019, 0.000020, 0.000001, 0.000010, 0.000029, 0.000052, 0.000023, -0.000010, 0.000035, 0.000013, 0.000036, 0.000041, 0.000017, 0.000039, 0.000002, 0.000010, 0.000019, 0.000025, 0.000005, -0.000006, 0.000008, 0.000020, -0.000011, 0.000025, -0.000004, -0.000008, 0.000014, 0.000002, -0.000014, -0.000014, 0.000003, 0.000039, -0.000005, -0.000020, 0.000013, 0.000003, -0.000004, 0.000007, -0.000020, -0.000026, -0.000007, 0.000038, -0.000001, -0.000008, -0.000012, 0.000014, 0.000014, -0.000005, -0.000012, 0.000021, -0.000043, 0.000012, -0.000012, -0.000001, 0.000013, -0.000023, -0.000026, -0.000014, 0.000007, 0.000004, 0.000017, -0.000033, -0.000005, 0.000022, -0.000003, -0.000012, 0.000011, -0.000011, -0.000020, -0.000023, -0.000000, 0.000001, -0.000039, -0.000034, -0.000005, 0.000024, 0.000008, -0.000010, -0.000031, -0.000039, 0.000004, 0.000018, 0.000000, -0.000031, -0.000006, -0.000026, -0.000001, -0.000003, -0.000016, -0.000022, -0.000045, -0.000003, -0.000012, -0.000024, -0.000032, -0.000029, 0.000011, 0.000003, -0.000014, -0.000024, -0.000020, -0.000036, -0.000002, 0.000007, -0.000001, -0.000025, -0.000032, -0.000021, -0.000029, -0.000023, -0.000019, -0.000019, -0.000028, -0.000013, -0.000015, -0.000015, -0.000001, 0.000006, -0.000023, -0.000023, -0.000003, -0.000011, -0.000025, -0.000022, -0.000008, 0.000001, 0.000001, -0.000015, -0.000001, -0.000014, -0.000026, 0.000024, -0.000016, -0.000007, -0.000016, -0.000012, 0.000013, 0.000005, -0.000018, -0.000015, -0.000009, -0.000016, -0.000007, -0.000008, -0.000020, -0.000011, -0.000010, -0.000011, -0.000017, -0.000010, -0.000000, 0.000000, 0.000001, 0.000011, -0.000008, -0.000013, -0.000009, 0.000005, -0.000017, -0.000001, 0.000002, -0.000001, -0.000023, 0.000001, -0.000019, -0.000009, -0.000016, 0.000006, -0.000003, -0.000016, -0.000009, -0.000002, -0.000021, 0.000009, -0.000015, -0.000004, -0.000005, 0.000004, -0.000004, -0.000000, -0.000008, 0.000007, -0.000001, -0.000011, -0.000001, -0.000000, -0.000015, -0.000019, -0.000006, -0.000007, -0.000020, -0.000052, -0.000046, -0.000030, -0.000018, -0.000007, -0.000042, -0.000036, -0.000024, -0.000020, -0.000017, -0.000019, -0.000018, -0.000023, -0.000035, -0.000000, -0.000003, -0.000020, -0.000036, -0.000031, -0.000016, -0.000013, -0.000004, -0.000051, -0.000051, -0.000026, -0.000028, -0.000037, -0.000034, -0.000062, -0.000038, -0.000040, -0.000011, -0.000036, -0.000030, -0.000045, -0.000026, -0.000021, -0.000002, -0.000027, -0.000016, -0.000057, 0.000016, 0.000049, 0.000041, 0.000036, 0.000008, 0.000022, 0.000066, 0.000039, 0.000002, 0.000026, -0.000013, 0.000051, -0.000003, 0.000022, 0.000016, -0.000020, -0.000009, 0.000019, 0.000052, 0.000046, 0.000014, 0.000008, 0.000000, 0.000082, 0.000067, 0.000026, 0.000056, 0.000011, 0.000032, 0.000102, 0.000036, 0.000043, 0.000003, 0.000040, 0.000060, 0.000052, 0.000029, 0.000049, 0.000004, 0.000024, 0.000069, 0.000059, 0.000029, 0.000024, 0.000000, -0.000001, 0.000057, -0.000006, -0.000018, 0.000030, 0.000031, 0.000022, 0.000028, 0.000013, 0.000043, 0.000005, 0.000012, 0.000049, 0.000022, 0.000012, 0.000006, -0.000015, 0.000005, 0.000032, 0.000014, -0.000010, 0.000004, 0.000003, 0.000030, 0.000025, -0.000006, 0.000003, 0.000016, 0.000002, 0.000005, -0.000027, 0.000028, -0.000019, -0.000008, -0.000015, 0.000030, -0.000005, -0.000010, -0.000030, 0.000003, 0.000025, 0.000033, -0.000024, 0.000012, 0.000008, -0.000004, -0.000008, 0.000005, -0.000001, -0.000016, -0.000024, 0.000003, 0.000004, -0.000005, -0.000019, -0.000024, -0.000017, 0.000007, 0.000018, -0.000020, -0.000015, -0.000006, 0.000012, 0.000004, -0.000008, -0.000018, -0.000009, 0.000004, -0.000018, -0.000019, -0.000001, -0.000013, 0.000021, 0.000008, 0.000012, 0.000001, -0.000019, -0.000028, 0.000012, 0.000007, 0.000003, -0.000029, -0.000015, -0.000015, -0.000008, -0.000010, -0.000002, -0.000006, -0.000018, -0.000001, -0.000018, -0.000004, -0.000026, -0.000015, 0.000002, -0.000004, -0.000005, -0.000016, -0.000026, 0.000004, 0.000020, -0.000006, -0.000021, -0.000028, -0.000006, -0.000016, 0.000001, -0.000017, -0.000014, -0.000020, 0.000007, -0.000003, -0.000011, -0.000010, -0.000008, -0.000018, 0.000016, -0.000002, -0.000017, -0.000016, -0.000007, 0.000004, 0.000000, 0.000000, 0.000001, -0.000025, -0.000014, -0.000013, -0.000008, -0.000016, -0.000032, -0.000019, -0.000002, 0.000002, 0.000001, -0.000014, -0.000011, -0.000020, -0.000000, -0.000002, -0.000012, 0.000002, -0.000011, -0.000013, -0.000000, -0.000009, -0.000001, -0.000007, -0.000008, 0.000006, -0.000013, -0.000015, -0.000019, -0.000018, -0.000006, 0.000010, -0.000011, -0.000011, -0.000013, -0.000019, -0.000020, -0.000005, -0.000013, -0.000015, -0.000020, -0.000015, -0.000001, -0.000006, 0.000003, -0.000006, -0.000004, -0.000001, -0.000007, -0.000003, -0.000016, -0.000003, 0.000005, -0.000004, -0.000009, -0.000014, -0.000021, 0.000000, -0.000013, 0.000004, -0.000012, -0.000023, -0.000013, -0.000012, 0.000001, 0.000003, -0.000006, -0.000017, -0.000003, -0.000023, -0.000028, -0.000020, -0.000016, -0.000012, -0.000022, -0.000040, -0.000025, -0.000023, -0.000016, -0.000014, -0.000017, -0.000007, -0.000021, -0.000019, -0.000000, -0.000017, -0.000026, -0.000032, -0.000023, -0.000016, -0.000004, -0.000015, -0.000038, -0.000004, -0.000043, -0.000024, -0.000009, -0.000031, -0.000042, -0.000016, -0.000019, 0.000007, -0.000049, -0.000019, -0.000023, -0.000014, -0.000007, -0.000004, -0.000048, -0.000019, -0.000029, 0.000019, 0.000047, 0.000033, -0.000003, 0.000027, 0.000000, 0.000044, 0.000017, 0.000021, 0.000004, 0.000010, 0.000026, 0.000014, -0.000004, -0.000011, 0.000008, 0.000016, 0.000041, 0.000028, 0.000002, 0.000022, 0.000016, 0.000030, 0.000076, 0.000008, 0.000037, 0.000036, 0.000038, 0.000076, 0.000053, 0.000014, 0.000042, -0.000003, 0.000072, 0.000029, 0.000029, 0.000014, 0.000020, 0.000026, 0.000064, 0.000016, 0.000011, 0.000032, 0.000019, 0.000019, 0.000031, -0.000029, 0.000024, 0.000025, 0.000010, 0.000016, 0.000025, 0.000007, 0.000028, 0.000010, 0.000050, 0.000013, -0.000001, 0.000011, 0.000013, 0.000014, 0.000040, -0.000014, -0.000007, 0.000014, 0.000009, 0.000014, 0.000014, -0.000020, 0.000000, 0.000006, 0.000018, -0.000013, 0.000000, 0.000004, 0.000013, -0.000001, 0.000015, -0.000009, -0.000000, -0.000012, 0.000006, 0.000018, 0.000014, -0.000008, -0.000001, 0.000002, -0.000005, 0.000011, 0.000013, -0.000019, 0.000001, -0.000007, 0.000011, -0.000003, -0.000013, -0.000014, -0.000007, -0.000004, 0.000004, -0.000016, 0.000005, 0.000003, 0.000000, 0.000016, 0.000006, -0.000022, -0.000001, 0.000004, -0.000015, -0.000023, -0.000005, -0.000006, 0.000012, -0.000023, 0.000004, -0.000003, -0.000021, 0.000006, -0.000002, -0.000008, -0.000002, -0.000015, -0.000002, -0.000002, -0.000014, -0.000003, -0.000015, 0.000004, -0.000011, -0.000009, -0.000003, -0.000007, -0.000010, -0.000005, -0.000001, -0.000012, -0.000012, -0.000012, -0.000004, 0.000006, 0.000003, -0.000019, -0.000030, -0.000008, 0.000007, -0.000008, -0.000015, -0.000011, -0.000008, -0.000003, -0.000016, -0.000010, -0.000008, -0.000014, -0.000004, -0.000007, -0.000019, -0.000015, -0.000014, -0.000006, 0.000006, 0.000004, -0.000004, -0.000017, -0.000000, -0.000012, -0.000013, -0.000017, -0.000016, -0.000010, 0.000003, -0.000003, -0.000009, -0.000017, -0.000014, -0.000016, 0.000002, -0.000001, -0.000012, -0.000011, 0.000002, -0.000001, -0.000007, -0.000001, -0.000010, -0.000000, -0.000002, -0.000013, 0.000001, -0.000022, -0.000010, -0.000001, -0.000010, -0.000009, -0.000012, -0.000015, -0.000014, -0.000020, -0.000008, -0.000010, -0.000003, -0.000009, -0.000007, 0.000004, 0.000004, -0.000006, -0.000001, 0.000003, -0.000009, 0.000000, -0.000008, -0.000007, 0.000007, 0.000009, -0.000001, -0.000005, -0.000001, -0.000007, -0.000010, -0.000001, -0.000017, -0.000015, -0.000003, -0.000020, -0.000006, -0.000002, -0.000011, -0.000007, -0.000009, -0.000004, 0.000004, 0.000003, -0.000006, 0.000004, -0.000001, -0.000005, -0.000001, -0.000017, -0.000011, -0.000018, -0.000008, -0.000005, -0.000032, -0.000038, -0.000023, -0.000025, -0.000009, -0.000033, -0.000017, -0.000017, -0.000023, -0.000005, -0.000006, -0.000014, -0.000018, -0.000014, -0.000014, -0.000011, -0.000001, -0.000011, -0.000014, -0.000028, -0.000018, -0.000018, -0.000011, -0.000033, -0.000026, -0.000016, -0.000018, -0.000030, -0.000016, -0.000022, -0.000015, -0.000008, -0.000016, -0.000020, -0.000025, -0.000010, -0.000015, 0.000017, 0.000031, 0.000022, 0.000026, 0.000003, 0.000022, 0.000054, 0.000011, 0.000003, 0.000009, 0.000036, 0.000012, 0.000010, 0.000002, -0.000002, 0.000009, 0.000024, 0.000018, 0.000013, 0.000011, 0.000017, 0.000001, 0.000021, 0.000008, 0.000066, 0.000048, 0.000012, 0.000027, 0.000057, 0.000025, 0.000020, 0.000025, 0.000053, 0.000017, 0.000028, 0.000010, 0.000015, 0.000026, 0.000032, 0.000027, 0.000028, 0.000026, 0.000020, 0.000016, -0.000000, -0.000010, 0.000036, 0.000030, -0.000010, 0.000016, 0.000029, 0.000010, 0.000012, 0.000017, 0.000019, 0.000004, 0.000006, 0.000017, 0.000005, 0.000007, 0.000006, 0.000006, 0.000009, 0.000010, 0.000005, 0.000021, -0.000025, -0.000013, 0.000025, -0.000003, -0.000002, -0.000007, 0.000012, 0.000001, 0.000003, 0.000003, 0.000000, 0.000001, 0.000005, 0.000020, 0.000003, 0.000007, -0.000011, 0.000019, -0.000012, 0.000001, 0.000012, 0.000000, -0.000005, -0.000008, 0.000003, -0.000005, -0.000013, -0.000003, -0.000016, -0.000000, 0.000006, -0.000004, -0.000012, -0.000006, -0.000005, 0.000022, -0.000020, -0.000011, -0.000003, -0.000007, -0.000014, -0.000012, 0.000004, 0.000006, -0.000003, 0.000003, -0.000006, -0.000003, -0.000003, -0.000007, -0.000011, 0.000006, -0.000002, -0.000001, -0.000016, -0.000011, -0.000006, -0.000001, -0.000016, -0.000016, 0.000005, -0.000006, -0.000008, -0.000012, -0.000005, 0.000005, -0.000002, -0.000002, -0.000008, -0.000004, 0.000002, -0.000017, -0.000019, -0.000005, -0.000011, -0.000006, -0.000014, -0.000007, -0.000005, -0.000004, -0.000003, -0.000008, -0.000001, 0.000004, -0.000002, -0.000011, -0.000008, -0.000010, -0.000000, -0.000006, 0.000002, 0.000007, -0.000008, -0.000009, -0.000016, -0.000001, -0.000005, -0.000015, -0.000002, -0.000004, -0.000003, -0.000002, -0.000015, -0.000014, -0.000012, 0.000000, 0.000009, -0.000010, -0.000008, -0.000004, -0.000006, -0.000003, 0.000001, 0.000003, -0.000005, -0.000006, -0.000007, -0.000007, -0.000001, -0.000003, -0.000002, -0.000003, -0.000008, -0.000010, -0.000009, -0.000008, -0.000006, -0.000009, 0.000007, -0.000009, 0.000001, 0.000001, -0.000012, -0.000004, -0.000010, 0.000006, -0.000000, -0.000009, -0.000004, -0.000006, 0.000003, -0.000001, -0.000009, -0.000002, -0.000008, -0.000009, -0.000009, -0.000013, -0.000010, -0.000014, -0.000002, -0.000003, -0.000002, -0.000003, -0.000005, -0.000007, 0.000004, 0.000003, -0.000008, -0.000012, -0.000001, -0.000002, -0.000002, -0.000007, -0.000001, -0.000001, -0.000002, -0.000006, -0.000002, -0.000004, -0.000010, -0.000011, -0.000015, -0.000029, -0.000020, -0.000017, -0.000008, -0.000026, -0.000019, -0.000014, -0.000017, -0.000015, -0.000024, -0.000010, -0.000014, -0.000007, -0.000007, -0.000010, -0.000016, -0.000011, -0.000006, -0.000010, -0.000008, -0.000013, -0.000011, -0.000027, -0.000021, -0.000014, -0.000016, -0.000015, -0.000023, -0.000007, -0.000019, -0.000021, -0.000019, -0.000014, -0.000010, -0.000012, -0.000014, -0.000010, -0.000012, -0.000013, -0.000009, -0.000014, 0.000029, 0.000027, 0.000015, 0.000010, -0.000003, 0.000043, 0.000014, 0.000010, 0.000014, 0.000017, 0.000016, 0.000002, 0.000011, -0.000003, 0.000008, 0.000013, 0.000021, 0.000005, 0.000018, -0.000007, 0.000011, 0.000010, 0.000010, 0.000057, 0.000024, 0.000022, 0.000025, 0.000031, 0.000034, 0.000016, 0.000037, 0.000019, 0.000022, 0.000022, -0.000001, 0.000015, 0.000020, 0.000015, 0.000028, 0.000026, 0.000017, 0.000017, -0.000007, 0.000014, 0.000010, 0.000017, 0.000013, 0.000008, 0.000026, 0.000002, 0.000011, 0.000013, 0.000010, 0.000020, 0.000007, 0.000014, -0.000001, 0.000007, 0.000021, 0.000018, 0.000007, 0.000025, 0.000006, -0.000006, -0.000008, -0.000009, 0.000009, 0.000003, -0.000002, -0.000004, 0.000013, 0.000007, 0.000004, 0.000003, 0.000008, 0.000010, 0.000018, -0.000004, 0.000002, 0.000007, -0.000004, -0.000006, 0.000002, 0.000017, -0.000001, -0.000008, -0.000003, -0.000002, -0.000008, -0.000005, -0.000003, -0.000001, 0.000008, 0.000001, -0.000006, -0.000005, 0.000003, 0.000002, 0.000007, -0.000012, -0.000002, 0.000002, -0.000012, -0.000013, 0.000008, 0.000010, -0.000002, -0.000005, -0.000003, -0.000008, -0.000002, -0.000013, -0.000010, 0.000004, 0.000006, -0.000002, -0.000011, -0.000006, -0.000005, -0.000003, -0.000010, -0.000005, -0.000001, -0.000005, -0.000011, -0.000017, 0.000001, -0.000000, 0.000001, 0.000002, -0.000003, -0.000006, -0.000003, -0.000006, -0.000008, 0.000005, -0.000006, -0.000012, -0.000003, -0.000003, -0.000010, -0.000002, -0.000009, -0.000007, 0.000004, -0.000005, -0.000008, -0.000015, -0.000007, -0.000010, -0.000001, 0.000005, -0.000004, -0.000001, -0.000019, -0.000008, -0.000005, -0.000007, -0.000012, -0.000009, -0.000003, -0.000001, -0.000007, -0.000009, -0.000007, -0.000005, 0.000009, -0.000006, -0.000014, -0.000007, -0.000007, -0.000002, -0.000002, 0.000001, -0.000010, -0.000003, -0.000009, -0.000004, 0.000001, -0.000006, -0.000004, -0.000002, -0.000009, -0.000012, -0.000005, -0.000009, -0.000006, -0.000003, -0.000003, -0.000003, -0.000004, -0.000009, -0.000008, -0.000004, 0.000002, 0.000003, -0.000009, 0.000000, -0.000007, 0.000003, 0.000000, -0.000007, -0.000006, -0.000008, -0.000003, -0.000004, -0.000009, -0.000006, -0.000008, -0.000000, -0.000007, -0.000004, -0.000002, -0.000008, -0.000000, -0.000005, -0.000005, -0.000002, -0.000005, 0.000003, -0.000003, 0.000002, 0.000001, 0.000002, -0.000002, -0.000007, 0.000003, -0.000009, -0.000005, -0.000009, -0.000008, -0.000000, -0.000004, -0.000006, -0.000001, -0.000004, 0.000003, -0.000001, -0.000023, -0.000016, -0.000015, -0.000004, -0.000009, -0.000019, -0.000002, -0.000010, -0.000002, -0.000013, -0.000005, -0.000009, -0.000013, -0.000007, -0.000010, -0.000015, -0.000013, -0.000012, 0.000001, -0.000010, -0.000018, -0.000004, -0.000019, -0.000022, -0.000017, -0.000017, -0.000010, -0.000018, -0.000012, -0.000018, -0.000012, -0.000010, -0.000008, -0.000006, -0.000018, -0.000011, -0.000008, -0.000015, -0.000012, -0.000006, -0.000006, -0.000000, 0.000012, 0.000020, 0.000013, 0.000005, 0.000007, 0.000028, 0.000009, 0.000008, -0.000001, 0.000012, 0.000007, 0.000016, 0.000004, 0.000003, 0.000009, 0.000011, 0.000015, 0.000006, 0.000012, 0.000008, 0.000018, -0.000007, 0.000028, 0.000024, 0.000018, 0.000038, 0.000008, 0.000030, 0.000023, 0.000028, 0.000012, 0.000015, 0.000019, 0.000014, 0.000023, 0.000013, 0.000014, 0.000028, 0.000014, 0.000016, 0.000017, -0.000006, 0.000022, 0.000010, 0.000013, 0.000013, 0.000010, 0.000013, 0.000015, 0.000014, 0.000010, 0.000008, 0.000010, 0.000000, 0.000013, 0.000000, 0.000007, 0.000020, 0.000003, 0.000015, -0.000004, 0.000001, 0.000008, -0.000007, -0.000000, 0.000015, 0.000004, -0.000004, 0.000009, 0.000002, 0.000001, 0.000012, 0.000009, 0.000015, -0.000003, -0.000006, 0.000003, -0.000002, -0.000011, 0.000010, 0.000003, -0.000002, -0.000002, 0.000001, -0.000008, -0.000005, -0.000006, -0.000003, 0.000003, -0.000004, -0.000004, -0.000003, -0.000004, -0.000002, -0.000000, 0.000007, 0.000001, -0.000005, -0.000010, -0.000006, -0.000005, 0.000007, 0.000001, -0.000004, -0.000007, -0.000003, 0.000004, -0.000011, -0.000006, -0.000001, 0.000004, 0.000000, -0.000010, -0.000008, -0.000006, 0.000000, -0.000001, -0.000004, -0.000004, -0.000004, -0.000009, -0.000008, -0.000004, 0.000000, 0.000003, 0.000001, -0.000007, -0.000008, 0.000001, -0.000002, -0.000008, -0.000008, -0.000001, -0.000007, -0.000003, -0.000003, -0.000004, -0.000001, -0.000003, 0.000003, -0.000010, -0.000010, -0.000008, -0.000003, -0.000005, -0.000002, -0.000001, -0.000005, 0.000002, -0.000003, -0.000007, -0.000005, -0.000001, -0.000005, -0.000006, -0.000004, -0.000002, -0.000005, -0.000003, -0.000002, 0.000000, -0.000004, -0.000003, -0.000005, -0.000010, -0.000000, -0.000003, 0.000005, -0.000008, -0.000005, -0.000007, -0.000003, 0.000003, -0.000009, -0.000006, -0.000008, -0.000004, -0.000008, -0.000006, -0.000007, -0.000003, -0.000005, 0.000001, 0.000001, -0.000003, -0.000007, -0.000006, -0.000004, 0.000000, -0.000002, -0.000005, -0.000004, -0.000008, -0.000002, -0.000002, -0.000003, -0.000006, -0.000006, -0.000009, -0.000005, -0.000002, -0.000006, 0.000002, -0.000008, 0.000002,
    };

    struct convolution_filter
    {
        static constexpr size_t size = 8192;
        static constexpr size_t buffer_size = 2 * size;
        static constexpr size_t mask = size - 1;

        std::array<real, buffer_size> history = {};
        size_t head = 0;

        const line* impulse = nullptr;

        /*
         *         N - 1
         *         _____
         *         \
         *  y[n] = /____ h[i] * x[n - i]
         *         i = 0
         */

        real filter(const real x)
        {
            history[head] = x;
            history[head + size] = x;
            real y = 0;
            for(size_t i = 0; i < size; i++)
            {
                y += impulse->operator[](i) * history[i + head];
            }
            head = (head - 1) & mask;
            return y;
        }

        void set_impulse(const line& impulse)
        {
            this->impulse = &impulse;
        }
    };

    template<size_t W, size_t H>
    struct alignas(std::hardware_destructive_interference_size) mailbox
    {
        /*
         * Recieve
         */

        std::atomic<real> throttle_open_ratio = 0.0_r;
        std::atomic<size_t> log_x = -1;
        std::atomic<size_t> log_y = -1;
        std::atomic<bool> injection_enabled = true;

        /*
         * Send
         */

        std::atomic<size_t> swap_drops = 0;
        std::atomic<real> engine_angular_velocity_r_per_s = 0.0_r;
        std::atomic<real> engine_load_torque_n_m = 0.0_r;
        std::array<std::array<std::atomic<real>, W>, H> port_open_ratios = {};
        std::array<std::array<std::atomic<bool>, W>, H> panics = {};
    };

    /*
     * [ ]  ...  [ ] | <- Source
     * [ ]  ...  [ | | <- Intake
     * [ ]  ...  [ ] | <- Intake Manifold
     * [ ]  ...  [ ] | <- Intake Runner
     * [ ]  ...  [ ] H <- Piston (PISTON_Y)
     * [ ]  ...  [ ] | <- Exhaust Runner
     * [ ]  ...  [ ] | <- Exhaust Manifold
     * [ ]  ...  [ ] | <- Exhaust
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
        template<size_t> class PISTONS,
        template<size_t> class CAMS,
        template<size_t> class SPARKPLUGS>
    struct as_engine : engine
    {
        real lumped_drag_torque_n_m = {};
        struct PISTONS<W> pistons = {};
        struct CAMS<W> inlet_cam = {};
        struct CAMS<W> outlet_cam = {};
        struct SPARKPLUGS<W> sparkplugs = {};
        std::array<struct flow<H, PISTON_Y>, W> flows = {};
        struct limiter limiter = {};
        struct throttle throttle = {};
        struct crankshaft crankshaft = {};
        struct flywheel flywheel = {};
        struct dc_filter dc = {};
        struct gain_filter gain = {};
        struct clamp_filter clamp = {};
        struct convolution_filter convolution = {};
        struct diags diags = {};
        struct pipe<W, PIPE_CELLS, PIPE_SUBSTEPS> pipe = {};
        line audio_signal = {};
        struct mailbox<W, H> mailbox = {};
        mutable std::vector<float> audio_data = {};
        line pipe_pressure = {};
        std::mutex swap_mutex = {};

        void log_step(const size_t x, const size_t y)
        {
            if(x < W and y < H)
            {
                #define X(name) diags.back[g_##name].push_back(flows[x].name[y]);
                FLUIDS(X)
                #undef X
                if(y == PISTON_Y)
                {
                    #define X(name) diags.back[g_##name].push_back(pistons.name[x]);
                    PISTONS(X)
                    #undef X
                }
            }
        }

        real calc_system_acceleration(const real load_torque_n_m)
        {
            /*
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
            t -= lumped_drag_torque_n_m;
            t -= load_torque_n_m;
            return t / I;
        }

        fn void broadcast(const real throttle_open_ratio, const bool injection_enabled, const real load_torque_n_m)
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
            limiter.crankshaft_angular_velocity_r_per_s = crankshaft.angular_velocity_r_per_s;

            /*
             * Cam open ratios -> chamber open ratios.
             *
             */

            for(size_t x = 0; x < W; x++)
            {
                flows[x].chamber_nozzle_open_ratio[PISTON_Y - 1] = inlet_cam.open_ratio[x];
                flows[x].chamber_nozzle_open_ratio[PISTON_Y + 0] = outlet_cam.open_ratio[x];
            }

            /*
             * Throttle open ratios -> chamber open ratios.
             *
             */

            for(size_t x = 0; x < W; x++)
            {
                const real open_ratio = throttle.lookup(throttle_open_ratio);
                flows[x].chamber_nozzle_open_ratio[THROTTLE_Y] = open_ratio;
            }

            /*
             * Piston shapes <-> flow shapes.
             *
             */

            for(size_t x = 0; x < W; x++)
            {
                flows[x].piston_injection_enabled = injection_enabled;
                flows[x].piston_chamber_radius_m = pistons.diameter_m[x] / 2.0_r;
                flows[x].chamber_volume_m3[PISTON_Y] = pistons.volumes_m3[x];
                pistons.chamber_static_pressure_pa[x] = flows[x].chamber_static_pressure_pa[PISTON_Y];
            }

            crankshaft.angular_acceleration_r_per_s2 = calc_system_acceleration(load_torque_n_m);
            crankshaft.angular_velocity_r_per_s = fmax(crankshaft.angular_velocity_r_per_s, 0.0_r);
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
            flywheel.update();
            crankshaft.update();
            pistons.calc_volumetrics();
            broadcast(0.0_r, false, 0.0_r);
            remember_volumes();
            reset_chambers();
        }

        bool diags_swap()
        {
            if(swap_mutex.try_lock())
            {
                for(size_t i = 0; i < diags.front.size(); i++)
                {
                    std::swap(diags.front[i], diags.back[i]);
                }
                for(auto& line : diags.back)
                {
                    line.clear();
                }
                pipe.gather_pipe_pressure_signal();
                std::swap(pipe_pressure, pipe.pipe_pressure_signal);
                swap_mutex.unlock();
                return true;
            }
            else
            {
                /*
                 * Discard back if front in use by renderer.
                 * Addng more samples to back will distort diags oscilloscope trigger.
                 */

                for(auto& line : diags.back)
                {
                    line.clear();
                }
                return false;
            }
        }

        void sample_audio()
        {
            real x0 = pipe.calc_audio_sample();
            x0 = dc.filter(x0);
            x0 = convolution.filter(x0);
            x0 = gain.filter(x0);
            x0 = clamp.filter(x0);
            audio_signal.push_back(x0);
        }

        void update_pipe()
        {
            for(size_t x = 0; x < W; x++)
            {
                const real u = flows[x].nozzle_velocity_m_per_s[AUDIO_Y];
                const real Ts = flows[x].nozzle_static_temperature_k[AUDIO_Y];
                const real r = flows[x].nozzle_static_density_kg_per_m3[AUDIO_Y];
                pipe.in_velocity_m_per_s[x] = u;
                pipe.in_static_temperature_k[x] = Ts;
                pipe.in_static_density_kg_per_m3[x] = r;
            }
            pipe.update();
        }

        void update_limiter()
        {
            limiter.update();
        }

        void update_flywheel()
        {
            flywheel.update();
        }

        bool update_crankshaft()
        {
            const bool otto_cycled = crankshaft.update();
            if(otto_cycled)
            {
                return diags_swap();
            }
            return false;
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

        void post_mailbox(const size_t swap_drops)
        {
            mailbox.engine_angular_velocity_r_per_s = crankshaft.angular_velocity_r_per_s;
            for(size_t y = 0; y < H; y++)
            for(size_t x = 0; x < W; x++)
            {
                mailbox.port_open_ratios[y][x] = flows[x].chamber_nozzle_open_ratio[y];
                mailbox.panics[y][x] = flows[x].panic[y];
            }
            mailbox.swap_drops += swap_drops;
        }

        void run(const size_t steps) override
        {
            const real throttle_open_ratio = mailbox.throttle_open_ratio;
            const real load_torque_n_m = mailbox.engine_load_torque_n_m;
            const size_t log_x = mailbox.log_x;
            const size_t log_y = mailbox.log_y;
            const bool injection_enabled = mailbox.injection_enabled;
            audio_signal.clear();
            audio_signal.reserve(steps);
            size_t swap_drops = 0;
            for(size_t step = 0; step < steps; step++)
            {
                update_limiter();
                update_flywheel();
                if(update_crankshaft())
                {
                    swap_drops++;
                }
                update_cams();
                update_sparkplugs();
                update_pistons();
                update_ignition();
                update_flows();
                update_pipe();
                log_step(log_x, log_y);
                remember_volumes();
                const bool injection_overrided = injection_enabled && not limiter.limiting;
                broadcast(throttle_open_ratio, injection_overrided, load_torque_n_m);
                sample_audio();
            }
            post_mailbox(swap_drops);
        }

        size_t get_width() const override
        {
            return W;
        }

        size_t get_height() const override
        {
            return H;
        }

        size_t get_piston_y() const override
        {
            return PISTON_Y;
        }

        size_t get_audio_y() const override
        {
            return AUDIO_Y;
        }

        size_t get_throttle_y() const override
        {
            return THROTTLE_Y;
        }

        size_t get_bytes() const override
        {
            return sizeof *this;
        }

        std::string_view get_signal_name(const size_t index) const override
        {
            return g_signal_names[index];
        }

        const std::atomic<real>& get_angular_velocity_r_per_s() const override
        {
            return mailbox.engine_angular_velocity_r_per_s;
        }

        const std::atomic<real>& get_port_open_ratio(const size_t x, const size_t y) const override
        {
            return mailbox.port_open_ratios[y][x];
        }

        const std::atomic<bool>& get_panic(const size_t x, const size_t y) const override
        {
            return mailbox.panics[y][x];
        }

        size_t get_swap_drops() const override
        {
            return mailbox.swap_drops;
        }

        const line& get_signal(const size_t index) const override
        {
            return diags.front[index];
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

        const line& get_audio_signal() const override
        {
            return audio_signal;
        }

        const line& get_impulse_signal() const override
        {
            return g_impulse;
        }

        const std::vector<float>& get_audio_data() const override
        {
            audio_data.clear();
            const line& audio_signal = get_audio_signal();
            for(const real& x : audio_signal)
            {
                audio_data.push_back(static_cast<float>(x));
            }
            return audio_data;
        }

        const line& get_pipe_pressure_signal() const override
        {
            return pipe_pressure;
        }

        void set_throttle_open_ratio(const real open_ratio) override
        {
            mailbox.throttle_open_ratio = open_ratio;
        }

        void set_injection_on() override
        {
            mailbox.injection_enabled = true;
        }

        void set_injection_off() override
        {
            mailbox.injection_enabled = false;
        }

        void set_logger(const size_t x, const size_t y) override
        {
            mailbox.log_x = x;
            mailbox.log_y = y;
        }

        void set_swap_lock_on() override
        {
            swap_mutex.lock();
        }

        void set_swap_lock_off() override
        {
            swap_mutex.unlock();
        }

        void set_load_torque_n_m(const real load_torque_n_m) override
        {
            mailbox.engine_load_torque_n_m = load_torque_n_m;
        }
    };

    struct inline4 : as_engine<
        /* W             */ 4,
        /* H             */ 9,
        /* THROTTLE_Y    */ 2,
        /* PISTON_Y      */ 4,
        /* AUDIO_Y       */ 7,
        /* PIPE_CELLS    */ 256,
        /* PIPE_SUBSTEPS */ 10,
        inline_pistons,
        basic_cams,
        basic_sparkplugs>
    {
        inline4()
        {
            this->convolution.set_impulse(g_impulse);
            this->lumped_drag_torque_n_m = 30.2_r;
            this->limiter.max_angular_velocity_r_per_s = 700.0_r;
            this->limiter.limit_time_s = 0.07_r;
            this->flywheel.mass_kg = 15.0_r;
            this->flywheel.radius_m = 0.18_r;
            this->crankshaft.mass_kg = 12.5_r;
            this->crankshaft.radius_m = 0.045_r;
            this->crankshaft.angular_velocity_r_per_s = 150.0_r;
            this->pistons.diameter_m.fill(0.086_r);
            this->pistons.crank_throw_length_m.fill(0.043_r);
            this->pistons.connecting_rod_length_m.fill(0.145_r);
            this->pistons.connecting_rod_mass_kg.fill(0.45_r);
            this->pistons.head_mass_density_kg_per_m3.fill(2700.0_r);
            this->pistons.head_compression_height_m.fill(0.030_r);
            this->pistons.head_clearance_height_m.fill(0.007_r);
            this->pistons.friction_n_m_s2_per_r2.fill(0.00005_r);
            this->inlet_cam.ramp_theta_r.fill(g_pi_r * 1.0_r);
            this->outlet_cam.ramp_theta_r.fill(g_pi_r * 0.6_r);
            real theta0_r = 0.0_r;
            for(size_t i = 0; i < get_width(); i++)
            {
                this->pistons.theta0_r[i] = theta0_r;
                this->inlet_cam.engage_theta_r[i]  = theta0_r + g_otto_intake_cycle_r - 0.4_r;
                this->sparkplugs.engage_theta_r[i] = theta0_r + g_otto_combustion_cycle_r - 0.4_r;
                this->outlet_cam.engage_theta_r[i] = theta0_r + g_otto_exhaust_cycle_r + 0.9_r;
                theta0_r += g_otto_cycle_r / static_cast<real>(get_width());
            }
            for(auto& flow : this->flows)
            {
                flow.chamber_nozzle_open_ratio.fill(1.0_r);
                flow.chamber_nozzle_flow_area_m2 = {
                    0.00250_r, /* Atmospheric Source -> Intake           */
                    0.00120_r, /* Intake             -> Throttle         */
                    0.00085_r, /* Throttle           -> Runner           */
                    0.00090_r, /* Runner             -> Piston           */
                    0.00120_r, /* Piston             -> Runner           */
                    0.00150_r, /* Runner             -> Chamber1         */
                    0.00175_r, /* Chamber1           -> Chamber2         */
                    0.00200_r, /* Chamber2           -> Atmospheric Sink */
                };
            }
            /*                             Atmospheric Source    Intake   Throttle  Runner    Piston     Runner    Chamber1  Chamber2   Atmospheric Sink    */
            flows[0].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.0005_r, 0.0005_r, 0.0005_r, g_resevoir_volume_m3 };
            flows[1].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.0005_r, 0.0005_r, 0.0005_r, g_resevoir_volume_m3 };
            flows[2].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.0005_r, 0.0005_r, 0.0005_r, g_resevoir_volume_m3 };
            flows[3].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.0005_r, 0.0005_r, 0.0005_r, g_resevoir_volume_m3 };
            this->throttle.table = {
                0.00100_r,
                0.02500_r,
                0.25000_r,
                1.00000_r,
            };
            this->pipe.piston_connect_m = { 0.0_r, 0.38_r, 0.17_r, 0.63_r };
            this->pipe.mic_position0_m = 0.7_r;
            this->pipe.mic_position1_m = 1.2_r;
            this->pipe.length_m = 1.50_r;
            this->dc.set_cutoff_frequency(2.0_r);
            this->gain.ratio = 0.00003_r;
        }
    };

    std::unique_ptr<engine> new_engine(const type type)
    {
        std::unique_ptr<engine> engine;
        switch(type)
        {
        default:
        case type::inline4: engine = std::make_unique<ensim::inline4>(); break;
        }
        engine->reset();
        return engine;
    }
}
