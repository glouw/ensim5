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
        lane<W> piston_connect_ratio = {};
        real length_m = 0.0_r;
        real mic_position0_ratio = 0.5_r;
        real mic_position1_ratio = 0.5_r;
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
            const size_t x = Z * mic_position0_ratio;
            const size_t y = Z * mic_position1_ratio;
            return static_pressure_pa[x] + static_pressure_pa[y];
        }

        void inject()
        {
            for(size_t i = 0; i < W; i++)
            {
                const real r = in_static_density_kg_per_m3[i];
                const real u = in_velocity_m_per_s[i];
                const real Ts = in_static_temperature_k[i];
                const real ratio = piston_connect_ratio[i];
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

    static const line g_impulse = {
        -0.000013, 0.000013, -0.000014, 0.000014, -0.000015, 0.000016, -0.000016, 0.000017, -0.000018, 0.000018, -0.000019, 0.000020, -0.000020, 0.000021, -0.000022, 0.000023, -0.000024, 0.000025, -0.000025, 0.000026, -0.000027, 0.000028, -0.000029, 0.000030, -0.000031, 0.000032, -0.000033, 0.000034, -0.000035, 0.000036, -0.000038, 0.000039, -0.000040, 0.000041, -0.000043, 0.000044, -0.000045, 0.000047, -0.000048, 0.000050, -0.000051, 0.000053, -0.000055, 0.000057, -0.000058, 0.000060, -0.000062, 0.000064, -0.000066, 0.000068, -0.000071, 0.000073, -0.000075, 0.000078, -0.000080, 0.000083, -0.000085, 0.000088, -0.000090, 0.000093, -0.000096, 0.000099, -0.000102, 0.000105, -0.000108, 0.000112, -0.000115, 0.000119, -0.000122, 0.000126, -0.000130, 0.000134, -0.000138, 0.000142, -0.000146, 0.000150, -0.000155, 0.000159, -0.000164, 0.000169, -0.000174, 0.000179, -0.000184, 0.000190, -0.000195, 0.000201, -0.000207, 0.000213, -0.000220, 0.000226, -0.000233, 0.000240, -0.000247, 0.000255, -0.000263, 0.000271, -0.000279, 0.000288, -0.000297, 0.000307, -0.000316, 0.000327, -0.000337, 0.000348, -0.000360, 0.000372, -0.000385, 0.000399, -0.000413, 0.000428, -0.000444, 0.000460, -0.000478, 0.000497, -0.000518, 0.000539, -0.000563, 0.000588, -0.000616, 0.000646, -0.000679, 0.000716, -0.000757, 0.000803, -0.000855, 0.000915, -0.000986, 0.001070, -0.001173, 0.001303, -0.001472, 0.001706, -0.002055, 0.002642, -0.003868, 0.008250, 0.023161, -0.004363, 0.001524, -0.001355, 0.000166, -0.000503, -0.000487, 0.000094, -0.001123, 0.000872, 0.000899, 0.017692, 0.005497, -0.003132, 0.000822, -0.001665, 0.000065, -0.001179, -0.000267, -0.000947, -0.000286, 0.003208, -0.001375, 0.000068, -0.001718, 0.000938, -0.003856, 0.016271, 0.006191, -0.004664, 0.002443, -0.005633, 0.047883, 0.002839, -0.003666, -0.000070, 0.003364, 0.015304, -0.004623, 0.001542, -0.007249, 0.016966, 0.037648, 0.007252, -0.004204, 0.001451, -0.001009, -0.003773, -0.001236, -0.004770, 0.003958, 0.006402, -0.005924, -0.000496, -0.004395, 0.034034, -0.005106, -0.001146, -0.005082, -0.000831, -0.005631, 0.000399, -0.008616, 0.019346, 0.067162, -0.006790, 0.000182, -0.011654, 0.038287, 0.022475, -0.018464, 0.012850, 0.026170, -0.002134, -0.014779, -0.001872, 0.002316, -0.005282, -0.006644, -0.002298, -0.008944, 0.000666, -0.020176, 0.015713, 0.012605, -0.006821, -0.009414, 0.016805, 0.025503, -0.014244, -0.000386, 0.036903, 0.005740, -0.016572, 0.002294, 0.024388, -0.000846, -0.009732, -0.003715, -0.008140, -0.008245, -0.008746, -0.003905, -0.007818, -0.004891, -0.008253, -0.003287, -0.009037, 0.001166, 0.013570, -0.010291, 0.000925, 0.015359, -0.003453, -0.010487, 0.001864, 0.015345, -0.005987, -0.002768, -0.013055, -0.006844, -0.005590, -0.006323, -0.004984, -0.006207, -0.003953, -0.007925, 0.000654, -0.002170, -0.012452, 0.003036, 0.014376, 0.005842, -0.003449, -0.006066, -0.009597, -0.005310, -0.004576, -0.007702, -0.005395, -0.004828, -0.004688, -0.004604, 0.004618, -0.010552, 0.006073, 0.001183, -0.004520, -0.003473, -0.000864, -0.012605, 0.026577, 0.018652, -0.014559, -0.000494, -0.008262, -0.002863, 0.023678, -0.006177, -0.004825, 0.006476, -0.005978, -0.002854, -0.002395, -0.010478, 0.020542, 0.013572, -0.012496, 0.000033, -0.010260, 0.003594, 0.010925, -0.012558, 0.007118, 0.022337, -0.011010, -0.002058, -0.007488, -0.003585, 0.005552, -0.003600, -0.006399, -0.003966, -0.006028, -0.003821, -0.006499, -0.002504, 0.013573, -0.013572, 0.009466, 0.005398, -0.005849, 0.005987, 0.004860, 0.018221, -0.005314, -0.005978, -0.000108, 0.000630, 0.011442, -0.008505, -0.003292, 0.009640, -0.005508, -0.006236, 0.002412, 0.010720, 0.000835, 0.004427, -0.004496, 0.001681, 0.005715, -0.011325, -0.000441, 0.011750, -0.008082, -0.004389, -0.004745, 0.006163, -0.000862, 0.002937, -0.009502, -0.001569, -0.012205, 0.003853, 0.010016, 0.000264, 0.005059, 0.006359, -0.000690, 0.008628, -0.002520, -0.005545, -0.004931, -0.006641, 0.001000, 0.002637, -0.009330, 0.005936, 0.002804, 0.010465, 0.005858, -0.001098, -0.006591, -0.004576, -0.005833, -0.004813, -0.005369, -0.005006, -0.004919, -0.005192, 0.006957, -0.005088, 0.003550, 0.003240, -0.002321, -0.001045, -0.007619, -0.002536, -0.006972, -0.001705, -0.005234, -0.007405, -0.001886, 0.003751, -0.004534, 0.001522, -0.007363, -0.004802, -0.006119, -0.002241, -0.003177, -0.002905, -0.005460, 0.005252, -0.006159, -0.000600, -0.000977, -0.007654, 0.007341, -0.000191, -0.005729, -0.001008, -0.004217, 0.004093, -0.006704, 0.000934, -0.003471, 0.002033, 0.009512, -0.005712, 0.000488, 0.001593, 0.009475, -0.001673, 0.001208, 0.003667, 0.005791, 0.005397, -0.003960, -0.000504, -0.007107, 0.009616, 0.009500, -0.000793, -0.009651, 0.016735, 0.015491, -0.001327, 0.009455, 0.003399, 0.008210, -0.007293, -0.003487, -0.004370, 0.002208, 0.007227, -0.000818, 0.016863, 0.012633, 0.018765, -0.007153, 0.005532, 0.000542, 0.001261, -0.008980, 0.001210, -0.008037, 0.002520, 0.005848, 0.001245, 0.009354, -0.007503, -0.003722, 0.007763, 0.007207, -0.003296, -0.004789, 0.003340, 0.006097, -0.000182, 0.004341, 0.011304, 0.005267, 0.006714, 0.002172, -0.001319, 0.007203, 0.006548, 0.002380, 0.004358, -0.003443, -0.003222, 0.000314, -0.004703, -0.003501, 0.000609, 0.008980, -0.005565, -0.003097, -0.007135, 0.015216, -0.000244, -0.002150, 0.001111, 0.007497, 0.001388, 0.002480, -0.002741, -0.006644, -0.004690, 0.002167, -0.006099, -0.003487, 0.007960, 0.009128, -0.011047, 0.001326, -0.000709, -0.004854, -0.007101, -0.000183, 0.003899, 0.002615, -0.005617, 0.007837, 0.010403, 0.006682, 0.000787, 0.003909, 0.010112, -0.007078, -0.008068, -0.003060, -0.006635, 0.001483, -0.000813, -0.003722, 0.003218, -0.006213, 0.006929, -0.010721, 0.001395, 0.002114, -0.006898, -0.001989, -0.004012, 0.005136, 0.000841, 0.008850, 0.009855, -0.002898, -0.003621, -0.002071, -0.004769, -0.007498, -0.002646, 0.007701, -0.004624, -0.003371, -0.004072, -0.007267, -0.001232, -0.003708, -0.005624, 0.007631, 0.001900, -0.000618, -0.001286, -0.004651, -0.001733, 0.000658, 0.001675, 0.005003, -0.010049, -0.003604, -0.008865, -0.006780, -0.006059, 0.002943, -0.005462, 0.004842, -0.001115, -0.005001, 0.007314, 0.000517, -0.000473, -0.007277, 0.010063, 0.001123, -0.006187, 0.006346, 0.000218, 0.005868, -0.003976, 0.005513, -0.007525, 0.005045, -0.007452, -0.000367, 0.000805, 0.005510, -0.002250, -0.002999, 0.002683, -0.001716, 0.005014, -0.000061, 0.002206, -0.000911, -0.003159, -0.004080, -0.004160, -0.004923, -0.000913, 0.002450, -0.000457, 0.000342, 0.007212, -0.004150, -0.004695, -0.003782, -0.004433, -0.006524, 0.003451, 0.005940, 0.004657, 0.000368, 0.008385, -0.006012, 0.010463, 0.002239, -0.005526, -0.004939, 0.005988, -0.001769, 0.001042, 0.000704, -0.000115, 0.000843, -0.005619, -0.006630, 0.000015, -0.003902, 0.000411, 0.003425, 0.000648, -0.000269, 0.002331, 0.000899, 0.004158, -0.001609, -0.005446, -0.002127, -0.002787, 0.002072, 0.002266, -0.004892, -0.003231, -0.002172, -0.008250, 0.007821, -0.000641, -0.005650, -0.002904, -0.000927, -0.000875, 0.001721, -0.002037, -0.004669, -0.001125, -0.006085, -0.001376, -0.002748, 0.001197, -0.002138, -0.001807, 0.001887, 0.001511, -0.000114, 0.000048, 0.000049, -0.003093, 0.006079, 0.003816, -0.004272, -0.003390, -0.002591, -0.004323, -0.000592, -0.003605, -0.004773, -0.000892, -0.005300, -0.004091, -0.000289, 0.007558, 0.002810, -0.005233, -0.002318, -0.003707, 0.005090, 0.000747, -0.003627, -0.003037, -0.006108, -0.001925, 0.003111, -0.006860, 0.000225, -0.000276, -0.003687, 0.004747, 0.004975, 0.000056, -0.000494, -0.004024, -0.003071, 0.001162, 0.005378, -0.001830, -0.001702, -0.000137, -0.002421, -0.001984, 0.006600, -0.000208, 0.006899, 0.002170, 0.000053, 0.001602, -0.003958, 0.010540, 0.001917, 0.004543, 0.001241, -0.001075, 0.003046, -0.002140, 0.005605, 0.008097, 0.006370, 0.002481, 0.001351, -0.003659, -0.000311, -0.000174, 0.003121, 0.009021, 0.000125, 0.002754, -0.002199, -0.001574, 0.003158, 0.001794, 0.000117, 0.007019, 0.003550, -0.004003, -0.001631, 0.002011, 0.003314, 0.005592, -0.005454, 0.006758, -0.003425, 0.005229, 0.003058, 0.001463, 0.002266, 0.001729, 0.005294, -0.002345, -0.000812, -0.005712, -0.001484, 0.003869, 0.002814, -0.000124, 0.005952, -0.001957, 0.004023, -0.003034, 0.004303, 0.001013, 0.003006, 0.003227, 0.005063, -0.002279, -0.003596, -0.001145, -0.002245, 0.004733, -0.000286, 0.004405, 0.000878, 0.000302, -0.000843, -0.002301, -0.001168, 0.000213, -0.001272, -0.003852, -0.001845, -0.006168, 0.000997, 0.009563, -0.006202, -0.001330, -0.002944, 0.001170, 0.006512, -0.001212, 0.000685, -0.006713, -0.001728, 0.003531, -0.002341, -0.005151, -0.003846, 0.001065, 0.004019, -0.000253, -0.003661, -0.002448, 0.006182, 0.003700, -0.001714, -0.003315, -0.005095, -0.000326, 0.000437, -0.000101, -0.004561, 0.000213, -0.000484, -0.005445, 0.001885, -0.000489, -0.000242, 0.001613, -0.000773, -0.003229, -0.001568, -0.003852, -0.001603, -0.001546, -0.002347, 0.001385, -0.001887, 0.003701, -0.000930, 0.001569, -0.004836, 0.001536, -0.000273, -0.002496, -0.000544, -0.000242, 0.000772, -0.003806, 0.000566, -0.001983, 0.001907, -0.002195, 0.003570, -0.000128, -0.000385, -0.003236, -0.005902, -0.001787, 0.003077, 0.000512, -0.000855, -0.003351, -0.003330, 0.003457, 0.001418, 0.003189, -0.000445, 0.000229, 0.001404, -0.001196, -0.001593, 0.001744, -0.001674, 0.004620, -0.005608, 0.002396, -0.002422, 0.001851, 0.001809, 0.000664, 0.000453, 0.000885, -0.001176, -0.002222, 0.003551, 0.000528, -0.003795, -0.000625, -0.002576, 0.004555, 0.002935, -0.002772, -0.001419, -0.005199, 0.004027, 0.002520, -0.003647, -0.003881, -0.002543, 0.002870, -0.003367, -0.002017, -0.002792, 0.006103, 0.003311, -0.003766, -0.001984, 0.001634, 0.004000, -0.002542, -0.004275, -0.002821, 0.000735, -0.001307, -0.003050, -0.001635, -0.002273, 0.003321, -0.000750, -0.003368, -0.000099, -0.003044, 0.000759, -0.000606, 0.001515, -0.000072, -0.003523, -0.001095, 0.002140, 0.005003, -0.000309, -0.004011, -0.003447, -0.002734, 0.001407, -0.000968, -0.001019, 0.002119, -0.002933, -0.000205, 0.002858, -0.000383, 0.000258, -0.004003, 0.000793, 0.001063, -0.003481, -0.001428, -0.002140, 0.000547, -0.000051, -0.002465, -0.001476, 0.003169, -0.001801, -0.003244, -0.003035, 0.001905, -0.000791, -0.000889, -0.002039, 0.001540, -0.002294, -0.002267, -0.000613, -0.000097, 0.001004, -0.002529, 0.000696, -0.000299, 0.000801, -0.002233, 0.005026, 0.002425, 0.001715, -0.004790, 0.000999, 0.002540, -0.003386, -0.001708, -0.002696, 0.002120, -0.001581, -0.001691, -0.000194, 0.002419, -0.001152, -0.002714, -0.002015, 0.001594, -0.001734, -0.002925, -0.002422, 0.000303, 0.000220, -0.002865, -0.001735, 0.000927, 0.002277, -0.002626, -0.002316, 0.000518, 0.003954, -0.002162, 0.000363, 0.000698, -0.000202, 0.000731, 0.000169, 0.000022, -0.001188, -0.000460, 0.000946, 0.000997, 0.001750, -0.001112, 0.003291, 0.004244, 0.004110, 0.000327, -0.000422, -0.000492, 0.000271, 0.002392, 0.001015, -0.000426, -0.001229, 0.002701, 0.003301, 0.002368, -0.001875, -0.001330, 0.000038, 0.001643, 0.001996, 0.004269, 0.003721, 0.000475, 0.003938, 0.002748, 0.002888, 0.000237, 0.000482, 0.001846, 0.001521, 0.002162, 0.001162, 0.003146, -0.000980, -0.000037, 0.003148, 0.000943, -0.001952, 0.000121, 0.000351, 0.001475, 0.000318, -0.000005, -0.000520, -0.002473, 0.002050, 0.003032, 0.000590, -0.002556, 0.001689, 0.003148, 0.000943, -0.003206, -0.000239, 0.002092, 0.001831, -0.000322, 0.002589, -0.002057, 0.000422, 0.001257, 0.003171, 0.000342, -0.001285, -0.000747, 0.004199, -0.002081, -0.002377, -0.000480, -0.000269, -0.001496, -0.002992, 0.000193, -0.001166, -0.000563, -0.001162, 0.002269, -0.003339, 0.000235, -0.000257, -0.000265, -0.001917, 0.000864, -0.002263, -0.002358, -0.001466, 0.000558, -0.001736, -0.001252, 0.000164, 0.000420, -0.000131, -0.001604, -0.001588, 0.000120, -0.000662, 0.003146, -0.000134, -0.001921, -0.000053, 0.000194, 0.001654, -0.000432, -0.001923, 0.000014, -0.001035, 0.002879, 0.001656, -0.002448, -0.001378, -0.000424, 0.000215, -0.000237, -0.002110, 0.000584, 0.001257, 0.000342, -0.000494, -0.001759, -0.001762, 0.000838, -0.002536, -0.003260, -0.000792, 0.002950, 0.000470, -0.003795, 0.000624, -0.000138, -0.000692, 0.000843, 0.000420, -0.002142, -0.002244, 0.002212, 0.003492, -0.003252, -0.000684, -0.000122, 0.000140, -0.002194, 0.001578, 0.001735, -0.000599, -0.001147, 0.002644, -0.000809, -0.002563, -0.000461, -0.002131, 0.000017, -0.001854, -0.000888, -0.002660, -0.000783, -0.000385, -0.002868, -0.002053, -0.001463, 0.000817, 0.000459, 0.005471, 0.002235, -0.002518, -0.000339, 0.000198, 0.000853, -0.000919, -0.001474, -0.001458, 0.000507, -0.000625, -0.001144, -0.001725, 0.001941, -0.003203, -0.001955, 0.001719, 0.000509, 0.000146, -0.002906, -0.002760, -0.001115, -0.000691, 0.000879, -0.000089, -0.001013, -0.000750, -0.001367, 0.000167, 0.000042, -0.001428, -0.002904, 0.000525, 0.000056, -0.000578, -0.002589, 0.001049, -0.000012, -0.001926, 0.000055, 0.000178, -0.002196, 0.001356, -0.000960, -0.002805, -0.000059, 0.001337, -0.000691, -0.001343, 0.001692, -0.000104, -0.000349, 0.000131, -0.001352, -0.001791, -0.001669, -0.000507, 0.000739, -0.001197, 0.000593, -0.000533, -0.001641, 0.002150, -0.001821, -0.001431, 0.000939, -0.000296, -0.000897, 0.000392, -0.000844, -0.001613, -0.000975, 0.001177, 0.002830, 0.000932, -0.000623, -0.000620, 0.001988, -0.001687, 0.000228, 0.000631, -0.000482, -0.000128, 0.000579, -0.000155, 0.000470, -0.000075, -0.001851, -0.000323, 0.000889, -0.001389, -0.001177, -0.000293, 0.000077, 0.000212, 0.000107, -0.000078, -0.000375, -0.000812, -0.000459, -0.001381, -0.000121, 0.003058, -0.000997, -0.001971, 0.001645, -0.000992, 0.002288, 0.002511, -0.001479, -0.000942, 0.000995, 0.001018, -0.002053, -0.000293, -0.000196, -0.001305, -0.000481, 0.000129, -0.001461, -0.002273, 0.000050, 0.000065, 0.001782, 0.000036, 0.000499, 0.001852, 0.000629, -0.000053, 0.001679, 0.001034, -0.001332, 0.000200, 0.002586, 0.000180, -0.000733, 0.000794, 0.001594, 0.000523, -0.001719, -0.000771, 0.001720, -0.001564, -0.000232, -0.000213, -0.001492, 0.001235, -0.000054, -0.000859, 0.001902, 0.000217, -0.000379, 0.001980, 0.001760, -0.000518, 0.000694, 0.002334, 0.002263, 0.001903, 0.000123, 0.000904, 0.003597, 0.000894, -0.001883, 0.002557, -0.000273, 0.000420, 0.001389, 0.000084, 0.000465, 0.000653, -0.000215, -0.001432, 0.001986, 0.000646, -0.000600, 0.001643, 0.000510, -0.001200, -0.000696, -0.000806, -0.000231, 0.000657, -0.001561, 0.001047, 0.000108, 0.001155, -0.001121, -0.000622, 0.000455, -0.000190, 0.002238, -0.000345, 0.000247, 0.002252, -0.000737, -0.001249, 0.000298, -0.001576, 0.000064, -0.000178, 0.001765, -0.001371, -0.000710, 0.000016, -0.000060, -0.000666, 0.000153, 0.001564, 0.000714, -0.001295, 0.000576, -0.000821, -0.000283, -0.000117, -0.000078, -0.000347, 0.000549, 0.003441, 0.000228, -0.002258, 0.002326, -0.000155, -0.000731, 0.000020, -0.001218, 0.000297, 0.000816, 0.000160, -0.002812, -0.000732, -0.001935, 0.000891, 0.000030, -0.002519, -0.000305, 0.000924, -0.000165, -0.000900, 0.000136, 0.001252, 0.000413, -0.000635, -0.000180, -0.000015, -0.001442, -0.001258, -0.001100, -0.000442, -0.001345, 0.001005, -0.000872, -0.000207, -0.000728, 0.000127, -0.001312, 0.000789, -0.000856, 0.000838, 0.001326, -0.000505, -0.000926, -0.001822, -0.000099, -0.000930, -0.000379, 0.000388, 0.000713, -0.000298, -0.001482, 0.000572, -0.000464, 0.000127, -0.000070, -0.001240, 0.000503, -0.001164, -0.000785, -0.001182, -0.000759, -0.000856, 0.000330, -0.001019, -0.000766, 0.001121, -0.000277, -0.000956, 0.000433, -0.001880, 0.001007, -0.000406, 0.000170, -0.000459, -0.001126, -0.000676, -0.000084, -0.001430, -0.000745, 0.000301, -0.001996, -0.000727, -0.001565, -0.000469, -0.000854, -0.001252, -0.000298, 0.001331, 0.000537, -0.000035, 0.000206, 0.001431, 0.000357, -0.001584, -0.000349, 0.001991, -0.000929, 0.000118, -0.001065, -0.001187, -0.001709, -0.000422, -0.000884, -0.001353, 0.000133, -0.000462, 0.000264, -0.000231, -0.000401, 0.000830, -0.000224, 0.000039, -0.000214, -0.001128, -0.000944, -0.000270, -0.000802, 0.000867, -0.000759, -0.000492, 0.000327, -0.001069, -0.000881, -0.000602, 0.000383, -0.000273, -0.000526, -0.000386, 0.000667, -0.000978, -0.000127, -0.000363, -0.000550, -0.001167, -0.000379, 0.000152, -0.000351, 0.000001, -0.000293, 0.000114, -0.000439, 0.000730, 0.000906, -0.001401, -0.000278, -0.001276, 0.000370, -0.000966, 0.000176, 0.001111, -0.000093, -0.000305, 0.000139, -0.000731, 0.000370, 0.000084, -0.000942, -0.000048, -0.001433, -0.000213, 0.000499, -0.001379, -0.001018, 0.000117, -0.000252, 0.000022, -0.000642, -0.000018, 0.000228, -0.000333, 0.000192, -0.000504, -0.000373, 0.000500, -0.000761, 0.000300, 0.000075, -0.000564, 0.001319, 0.000053, -0.000810, 0.000311, 0.000082, -0.000001, -0.000228, 0.000048, 0.000423, 0.000073, -0.000616, 0.000033, -0.000716, 0.000471, -0.000258, -0.000050, 0.000304, 0.000524, 0.000050, -0.000169, 0.000977, 0.000929, 0.000216, 0.000461, 0.000052, 0.000127, 0.000844, 0.000039, 0.001891, 0.001151, -0.000432, 0.000946, -0.000606, -0.000916, 0.000775, -0.000024, -0.000462, 0.000815, -0.000207, 0.000675, 0.000185, 0.000410, -0.000200, -0.000960, 0.000925, -0.000101, 0.000469, 0.000497, 0.000079, 0.001757, 0.001156, 0.000617, 0.001077, 0.000502, -0.000544, 0.000009, 0.002356, -0.000176, -0.001090, 0.000912, -0.000975, -0.000122, 0.002431, 0.000306, -0.000004, -0.000074, 0.000020, 0.000521, 0.000711, 0.002095, -0.000552, -0.000218, 0.000858, 0.000110, -0.000352, 0.000278, 0.000775, -0.000668, -0.000678, 0.000464, 0.000056, -0.000529, -0.000589, 0.001012, 0.000179, 0.000181, 0.000116, -0.000418, -0.000110, 0.000847, -0.000241, 0.000243, 0.000868, -0.000859, 0.000337, 0.000232, 0.000374, 0.000032, -0.000001, 0.000607, -0.000970, -0.000593, -0.000427, 0.000299, -0.000914, 0.000339, 0.000616, -0.001286, 0.000496, 0.000507, 0.000688, 0.001345, 0.000572, 0.000064, -0.000048, -0.000258, 0.000867, -0.000529, -0.000207, 0.000428, -0.001056, 0.000376, -0.000463, -0.000141, -0.000395, 0.000043, -0.000248, -0.000760, -0.000400, 0.000254, -0.000577, -0.000186, 0.000540, -0.000265, -0.000517, -0.000791, -0.000177, 0.000189, 0.000145, -0.000789, -0.001066, -0.001061, 0.000153, -0.000011, -0.000718, 0.000671, -0.000420, -0.000488, -0.000577, -0.000661, -0.000220, -0.000804, -0.000391, 0.000168, -0.000404, -0.000770, -0.000091, -0.001394, -0.000061, 0.000695, -0.000521, -0.000257, -0.000074, -0.000582, -0.000083, 0.000302, -0.001013, -0.000787, -0.000250, -0.001010, -0.000503, -0.000303, -0.000752, -0.000406, -0.000341, -0.000682, -0.000167, -0.000064, -0.000962, 0.001255, -0.000798, -0.000723, 0.000140, -0.000590, -0.000312, 0.000124, -0.000516, -0.000090, -0.000831, -0.000047, 0.000156, -0.000838, -0.000578, -0.000595, -0.000046, -0.000523, -0.000377, -0.000132, -0.000221, 0.000322, -0.000159, -0.000540, -0.000850, -0.000623, -0.000644, 0.000695, -0.000547, -0.000337, -0.000085, -0.000435, -0.000172, -0.000388, -0.000467, -0.001143, -0.000742, -0.000073, 0.000365, 0.000050, 0.000496, 0.000413, -0.000233, 0.000217, -0.000100, 0.000000, 0.000272, 0.000112, 0.000104, -0.000869, -0.000183, -0.000425, 0.000104, -0.000169, 0.000422, 0.000368, -0.000439, 0.000190, -0.000595, 0.000198, -0.000252, -0.000652, -0.000447, -0.000062, -0.000216, 0.000036, -0.000331, -0.000311, -0.000600, -0.000126, 0.000108, 0.000117, 0.000219, -0.000025, -0.000903, 0.000305, -0.000970, 0.000067, -0.000606, 0.000227, -0.000363, -0.000709, -0.000275, 0.000209, -0.000269, 0.000076, 0.000195, -0.000280, 0.000466, 0.000607, 0.000425, 0.000268, 0.000262, -0.000412, 0.000846, -0.000169, -0.000177, 0.000182, -0.000465, -0.000255, -0.000525, -0.000166, -0.000202, -0.000294, 0.000387, -0.000479, -0.000214, -0.000244, -0.000301, -0.000296, -0.000004, 0.000270, -0.000360, 0.000131, -0.000516, -0.000484, -0.000095, -0.000470, -0.000244, -0.000697, -0.000202, -0.000674, 0.000213, -0.000665, -0.000083, -0.000762, 0.000541, -0.000285, -0.000273, -0.000323, 0.000114, -0.000532, -0.000147, -0.000220, -0.000121, -0.000383, 0.000138, 0.000990, 0.000304, 0.000215, 0.001139, 0.000490, 0.000616, 0.000493, 0.000272, 0.000325, -0.000326, 0.000396, 0.000557, 0.000074, 0.000430, 0.000070, 0.000295, -0.000002, 0.000819, -0.000320, -0.000276, -0.000045, -0.000450, 0.000104, 0.000266, 0.000214, -0.000582, -0.000246, 0.000296, 0.001104, 0.000943, 0.000414, 0.001069, 0.000120, 0.000268, 0.001788, 0.000153, 0.000871, -0.000112, 0.000123, 0.001120, 0.000128, -0.000031, -0.000092, 0.000165, 0.000604, 0.000192, 0.001424, 0.000169, 0.000102, -0.000154, 0.000261, 0.000671, 0.000036, 0.000357, 0.001214, -0.000366, 0.000505, 0.001119, -0.000354, 0.000084, 0.000425, -0.000751, 0.000739, -0.000541, -0.000030, -0.000134, -0.000264, 0.000273, -0.000536, 0.000677, -0.000283, -0.000067, -0.000107, -0.000318, 0.000238, -0.000431, 0.000023, -0.000474, -0.000963, -0.000586, 0.000186, -0.000006, 0.000286, 0.000428, 0.000221, 0.000531, -0.000103, -0.000065, -0.000174, -0.000175, 0.000250, -0.000300, 0.000306, -0.000383, 0.000017, 0.000163, 0.000101, -0.000169, -0.000490, 0.000343, -0.000062, 0.000206, 0.000515, -0.000009, -0.000424, -0.000083, -0.000469, -0.000001, -0.000009, -0.000273, -0.000344, -0.000174, -0.000228, 0.000120, -0.000233, -0.000349, -0.000793, -0.000044, -0.000183, -0.000059, -0.000825, -0.000027, -0.000481, -0.000270, 0.000108, -0.000122, -0.000311, -0.000237, -0.000490, 0.000174, -0.000064, -0.000325, -0.000344, -0.000101, -0.000271, -0.000197, 0.000348, -0.000566, -0.000425, 0.000096, -0.000096, -0.000333, -0.000231, -0.000220, -0.000369, -0.000217, -0.000083, -0.000043, -0.000037, 0.000074, -0.000507, -0.000238, -0.000118, -0.000107, -0.000684, -0.000750, -0.000193, -0.000220, -0.000302, -0.000202, -0.000628, -0.000170, -0.000240, -0.000397, -0.000360, -0.000730, -0.000203, -0.000470, -0.000339, -0.000436, -0.000246, 0.000002, -0.000337, 0.000103, 0.000002, -0.000417, 0.000004, -0.000101, -0.000015, -0.000014, -0.000043, -0.000499, -0.000129, -0.000263, 0.000037, -0.000307, 0.000154, -0.000203, -0.000385, 0.000019, -0.000062, -0.000364, -0.000052, -0.000481, -0.000212, 0.000120, -0.000470, -0.000308, -0.000516, -0.000137, -0.000066, -0.000098, -0.000639, -0.000039, -0.000529, 0.000097, -0.000388, -0.000236, 0.000161, -0.000213, 0.000467, 0.000297, -0.000416, -0.000104, -0.000259, 0.000181, -0.000073, -0.000329, 0.000171, -0.000240, 0.000516, 0.000192, 0.000079, -0.000036, -0.000138, -0.000230, -0.000079, -0.000064, -0.000175, 0.000024, -0.000290, -0.000165, -0.000007, -0.000184, 0.000146, 0.000002, -0.000135, -0.000153, 0.000016, -0.000039, -0.000407, 0.000086, -0.000390, -0.000165, 0.000308, -0.000341, -0.000391, -0.000201, -0.000061, -0.000517, 0.000061, -0.000141, 0.000052, -0.000064, -0.000285, -0.000019, 0.000039, -0.000243, 0.000250, -0.000068, 0.000308, -0.000043, -0.000304, -0.000046, 0.000037, 0.000144, 0.000274, -0.000208, 0.000132, -0.000160, -0.000190, 0.000116, -0.000247, -0.000080, -0.000371, -0.000379, 0.000071, -0.000283, 0.000042, -0.000437, -0.000138, -0.000008, -0.000136, -0.000121, -0.000231, -0.000238, 0.000100, -0.000229, 0.000097, -0.000392, -0.000052, 0.000066, -0.000164, 0.000291, 0.000599, -0.000047, 0.000624, -0.000178, 0.000222, 0.001033, -0.000048, 0.000351, -0.000122, 0.000570, 0.000476, -0.000122, 0.000095, 0.000021, 0.000824, 0.000511, -0.000132, 0.000386, -0.000014, 0.000177, 0.000358, -0.000116, 0.000213, 0.000088, -0.000319, 0.000835, -0.000602, 0.000596, 0.000304, 0.000510, 0.000849, 0.000119, 0.000102, 0.000627, 0.000389, 0.000235, -0.000213, 0.000700, 0.000180, 0.000257, 0.000330, -0.000107, 0.000706, 0.000141, -0.000458, 0.000276, 0.000371, 0.000500, 0.000159, 0.000151, 0.000290, 0.000057, 0.000198, 0.000055, 0.000069, 0.000124, -0.000348, 0.000451, -0.000203, 0.000248, 0.000326, -0.000061, 0.000367, -0.000285, -0.000147, 0.000059, -0.000056, 0.000214, -0.000054, -0.000126, 0.000444, -0.000220, -0.000047, 0.000323, 0.000087, -0.000002, -0.000332, -0.000083, -0.000050, 0.000046, 0.000029, -0.000135, 0.000207, 0.000042, -0.000015, -0.000049, -0.000248, 0.000499, -0.000166, -0.000049, -0.000187, -0.000577, 0.000085, -0.000107, -0.000441, -0.000344, -0.000274, -0.000119, 0.000010, -0.000281, -0.000022, -0.000093, 0.000108, -0.000197, -0.000110, -0.000307, -0.000158, -0.000030, -0.000465, 0.000237, -0.000050, -0.000209, -0.000350, -0.000115, -0.000002, -0.000090, -0.000120, -0.000178, 0.000129, -0.000145, -0.000262, -0.000126, 0.000008, -0.000173, -0.000163, -0.000330, 0.000153, -0.000341, 0.000036, -0.000059, -0.000244, -0.000378, -0.000567, -0.000260, -0.000200, -0.000313, -0.000090, -0.000059, 0.000005, -0.000416, 0.000112, -0.000231, -0.000163, -0.000331, -0.000362, 0.000104, -0.000381, -0.000171, -0.000209, -0.000263, -0.000361, -0.000028, -0.000015, 0.000181, -0.000040, 0.000205, -0.000013, 0.000061, -0.000153, 0.000007, -0.000048, -0.000167, -0.000163, -0.000125, -0.000004, -0.000163, 0.000073, -0.000265, -0.000124, -0.000317, -0.000118, -0.000027, -0.000346, -0.000014, -0.000223, -0.000124, -0.000027, -0.000292, -0.000130, -0.000097, -0.000372, -0.000133, -0.000325, -0.000038, -0.000061, 0.000044, -0.000081, -0.000062, -0.000078, -0.000011, -0.000413, -0.000011, 0.000191, -0.000134, -0.000070, -0.000064, -0.000001, 0.000075, -0.000217, -0.000135, -0.000252, -0.000274, -0.000272, -0.000239, -0.000031, -0.000017, -0.000046, -0.000053, -0.000154, -0.000226, 0.000104, -0.000107, -0.000255, -0.000273, -0.000049, -0.000215, -0.000277, -0.000128, 0.000122, -0.000108, -0.000273, 0.000007, 0.000023, -0.000022, -0.000066, -0.000320, 0.000093, -0.000098, -0.000238, -0.000027, -0.000021, -0.000179, 0.000105, 0.000114, 0.000016, -0.000170, -0.000119, -0.000234, -0.000178, 0.000012, 0.000051, -0.000073, -0.000003, -0.000204, -0.000152, 0.000051, -0.000339, 0.000146, -0.000107, 0.000061, -0.000160, -0.000117, 0.000059, -0.000140, -0.000109, -0.000066, 0.000029, 0.000010, 0.000019, -0.000032, -0.000026, -0.000017, -0.000270, -0.000100, -0.000215, -0.000038, -0.000287, -0.000114, -0.000012, 0.000067, 0.000110, 0.000013, 0.000052, 0.000152, 0.000069, -0.000126, -0.000089, -0.000042, 0.000032, -0.000229, -0.000122, 0.000008, -0.000154, -0.000057, -0.000046, -0.000011, -0.000123, -0.000058, -0.000207, -0.000115, -0.000257, -0.000160, -0.000189, 0.000268, 0.000228, 0.000164, 0.000350, 0.000167, 0.000730, 0.000229, 0.000078, 0.000132, 0.000567, 0.000114, -0.000035, 0.000101, 0.000226, 0.000408, -0.000253, 0.000065, -0.000017, 0.000438, -0.000010, 0.000127, 0.000038, 0.000106, 0.000475, -0.000120, 0.000159, -0.000099, 0.000197, 0.000157, 0.000185, 0.000281, 0.000611, 0.000169, 0.000589, 0.000061, 0.000158, 0.000443, 0.000376, 0.000204, 0.000110, -0.000021, 0.000741, 0.000315, 0.000081, 0.000325, 0.000139, 0.000395, 0.000011, 0.000173, 0.000310, 0.000138, 0.000338, -0.000164, 0.000077, 0.000134, 0.000407, -0.000021, -0.000211, -0.000027, 0.000178, -0.000232, 0.000102, -0.000046, 0.000281, -0.000047, 0.000158, 0.000126, -0.000164, -0.000162, 0.000019, -0.000315, -0.000117, -0.000018, -0.000043, 0.000237, -0.000259, 0.000161, 0.000050, 0.000145, 0.000013, -0.000120, 0.000001, 0.000048, -0.000000, -0.000061, -0.000181, -0.000032, -0.000106, -0.000305, -0.000068, -0.000043, -0.000085, 0.000104, -0.000080, -0.000179, -0.000172, -0.000101, -0.000137, -0.000133, -0.000119, -0.000296, 0.000254, -0.000200, -0.000060, -0.000039, 0.000074, -0.000080, -0.000273, -0.000311, -0.000110, -0.000221, -0.000278, -0.000029, -0.000156, 0.000027, -0.000038, -0.000103, -0.000149, -0.000067, 0.000101, -0.000280, -0.000133, -0.000202, -0.000077, -0.000182, -0.000033, -0.000020, -0.000068, -0.000098, 0.000013, 0.000062, -0.000044, -0.000054, -0.000146, -0.000179, -0.000119, 0.000008, -0.000166, -0.000165, -0.000155, -0.000038, -0.000053, -0.000034, -0.000185, -0.000194, -0.000145, -0.000056, -0.000210, 0.000058, -0.000134, -0.000066, -0.000259, 0.000009, -0.000116, -0.000059, -0.000187, -0.000137, -0.000224, -0.000059, -0.000178, -0.000196, -0.000130, -0.000099, 0.000006, -0.000033, -0.000022, -0.000046, -0.000108, -0.000180, -0.000021, -0.000172, -0.000017, -0.000098, -0.000195, -0.000069, -0.000015, -0.000116, -0.000162, -0.000135, -0.000040, -0.000230, -0.000054, -0.000077, -0.000024, -0.000115, 0.000175, -0.000027, -0.000081, -0.000064, -0.000076, -0.000014, 0.000003, -0.000097, -0.000142, -0.000194, -0.000050, -0.000104, -0.000092, 0.000007, -0.000096, -0.000092, -0.000191, -0.000059, -0.000051, -0.000151, -0.000063, -0.000183, 0.000022, -0.000157, -0.000011, -0.000121, -0.000006, -0.000128, -0.000069, -0.000075, -0.000100, -0.000169, -0.000020, -0.000180, -0.000007, -0.000074, 0.000108, -0.000068, -0.000003, 0.000055, -0.000110, -0.000125, -0.000121, 0.000012, -0.000054, -0.000114, -0.000044, -0.000055, -0.000065, -0.000038, -0.000005, 0.000075, -0.000098, -0.000010, -0.000057, 0.000084, -0.000180, -0.000063, -0.000143, -0.000067, -0.000011, -0.000078, 0.000033, -0.000067, -0.000075, -0.000138, -0.000146, 0.000027, -0.000112, 0.000045, -0.000116, 0.000050, 0.000074, -0.000004, -0.000070, -0.000071, -0.000083, -0.000041, -0.000100, -0.000059, -0.000081, -0.000028, -0.000093, 0.000090, -0.000138, 0.000049, -0.000061, -0.000080, -0.000015, 0.000024, -0.000058, -0.000066, -0.000091, 0.000059, -0.000054, 0.000072, 0.000029, 0.000013, -0.000075, -0.000030, 0.000050, 0.000045, -0.000054, 0.000020, -0.000045, 0.000124, -0.000026, 0.000050, 0.000026, 0.000193, 0.000139, 0.000290, 0.000001, 0.000340, 0.000176, 0.000051, 0.000109, 0.000201, 0.000143, 0.000019, 0.000099, 0.000286, -0.000118, -0.000045, 0.000038, 0.000228, 0.000013, 0.000190, 0.000161, 0.000151, 0.000152, 0.000070, 0.000031, 0.000088, 0.000009, -0.000019, 0.000142, 0.000336, 0.000240, 0.000373, 0.000166, 0.000428, 0.000302, 0.000129, 0.000075, 0.000029, 0.000281, -0.000032, 0.000171, 0.000269, -0.000017, 0.000008, 0.000212, 0.000244, 0.000016, 0.000089, 0.000070, 0.000222, -0.000026, 0.000147, -0.000092, 0.000127, -0.000034, -0.000142, 0.000148, 0.000094, -0.000022, -0.000068, 0.000151, 0.000062, -0.000147, -0.000002, -0.000102, 0.000034, 0.000015, 0.000124, 0.000127, -0.000098, -0.000008, -0.000015, 0.000017, -0.000054, 0.000004, 0.000043, 0.000012, 0.000030, 0.000006, -0.000025, -0.000059, -0.000041, -0.000079, -0.000060, -0.000033, -0.000074, -0.000016, -0.000045, 0.000046, -0.000053, -0.000071, -0.000041, 0.000028, -0.000141, -0.000091, 0.000000, -0.000032, -0.000078, -0.000000, 0.000075, -0.000175, -0.000191, -0.000033, -0.000159, 0.000021, -0.000084, -0.000029, -0.000028, -0.000032, -0.000140, -0.000008, 0.000082, -0.000109, -0.000026, -0.000056, -0.000002, -0.000065, 0.000037, -0.000032, -0.000099, -0.000084, -0.000150, -0.000098, 0.000006, -0.000094, -0.000170, -0.000118, -0.000100, -0.000139, -0.000047, -0.000074, -0.000027, -0.000040, -0.000080, -0.000057, -0.000039, -0.000152, -0.000200, 0.000011, -0.000064, -0.000097, -0.000091, -0.000012, -0.000003, 0.000004, -0.000070, -0.000114, -0.000156, -0.000039, -0.000031, -0.000102, -0.000101, -0.000137, -0.000036, -0.000099, -0.000059, -0.000042, -0.000167, -0.000040, -0.000153, -0.000079, -0.000065, -0.000132, -0.000044, -0.000037, 0.000062, -0.000127, -0.000075, -0.000046, -0.000058, -0.000076, -0.000035, -0.000092, -0.000028, -0.000059, 0.000010, 0.000001, -0.000096, -0.000155, -0.000006, -0.000120, -0.000071, -0.000056, -0.000025, -0.000110, -0.000044, -0.000098, -0.000029, -0.000083, -0.000066, -0.000109, -0.000088, 0.000002, -0.000066, -0.000019, -0.000101, -0.000054, -0.000077, -0.000056, -0.000055, 0.000012, -0.000015, -0.000049, -0.000030, -0.000064, -0.000053, -0.000041, -0.000042, -0.000050, -0.000017, -0.000044, -0.000067, -0.000118, -0.000048, -0.000046, -0.000057, -0.000034, -0.000106, -0.000068, -0.000118, -0.000072, -0.000134, -0.000047, -0.000096, -0.000035, -0.000008, -0.000037, -0.000024, 0.000008, 0.000020, -0.000039, -0.000040, -0.000094, -0.000038, -0.000081, -0.000037, -0.000065, -0.000065, -0.000012, -0.000053, 0.000001, -0.000079, -0.000091, -0.000090, -0.000040, -0.000035, -0.000053, -0.000054, -0.000028, -0.000062, -0.000035, -0.000023, -0.000037, -0.000012, 0.000013, 0.000001, -0.000033, 0.000050, 0.000019, 0.000009, 0.000001, 0.000016, -0.000090, -0.000071, -0.000045, 0.000010, -0.000028, -0.000022, 0.000017, -0.000021, 0.000008, -0.000074, -0.000029, -0.000054, -0.000080, -0.000036, -0.000050, -0.000059, 0.000004, -0.000021, 0.000022, -0.000043, -0.000016, 0.000014, 0.000012, -0.000015, -0.000015, 0.000060, -0.000006, 0.000041, 0.000037, 0.000037, -0.000012, 0.000057, 0.000032, 0.000201, 0.000152, 0.000162, 0.000136, 0.000289, 0.000053, 0.000109, 0.000163, 0.000128, 0.000041, 0.000179, 0.000142, 0.000027, 0.000046, 0.000149, 0.000117, 0.000037, 0.000069, 0.000091, 0.000042, -0.000055, 0.000057, 0.000083, -0.000025, 0.000074, 0.000141, -0.000012, 0.000048, 0.000172, 0.000271, 0.000105, 0.000213, 0.000158, 0.000180, 0.000129, 0.000198, 0.000185, 0.000026, 0.000220, 0.000139, 0.000001, 0.000003, 0.000125, 0.000147, 0.000115, 0.000124, 0.000103, 0.000030, 0.000005, -0.000029, 0.000047, -0.000009, 0.000111, 0.000129, 0.000038, -0.000027, 0.000057, 0.000065, 0.000015, -0.000039, 0.000031, 0.000048, 0.000004, 0.000016, 0.000010, -0.000092, -0.000035, -0.000018, -0.000079, -0.000020, -0.000046, 0.000045, 0.000003, -0.000027, -0.000046, -0.000016, -0.000142, -0.000078, 0.000034, 0.000022, 0.000015, 0.000040, -0.000044, -0.000045, -0.000048, 0.000042, 0.000010, 0.000031, -0.000026, -0.000083, -0.000081, 0.000000, -0.000038, -0.000063, -0.000022, -0.000074, -0.000045, 0.000001, -0.000107, -0.000072, -0.000031, -0.000035, -0.000008, -0.000044, -0.000110, -0.000036, -0.000044, -0.000031, -0.000100, -0.000051, -0.000039, -0.000043, -0.000069, -0.000056, -0.000042, -0.000063, -0.000035, -0.000127, -0.000032, -0.000077, -0.000002, -0.000015, -0.000091, -0.000039, 0.000001, -0.000035, -0.000063, 0.000021, -0.000063, -0.000025, -0.000057, -0.000111, -0.000059, -0.000033, -0.000035, -0.000060, -0.000091, -0.000053, -0.000038, -0.000048, -0.000042, -0.000040, -0.000043, -0.000038, -0.000044, -0.000020, -0.000107, -0.000037, -0.000065, -0.000070, -0.000079, -0.000002, -0.000077, -0.000074, -0.000057, -0.000057, -0.000058, -0.000027, -0.000067, -0.000064, -0.000076, -0.000006, -0.000084, -0.000079, -0.000066, -0.000031, -0.000054, 0.000045, -0.000057, -0.000060, -0.000053, 0.000014, -0.000050, -0.000075, -0.000088, -0.000032, -0.000053, -0.000050, -0.000008, -0.000049, -0.000019, 0.000016, -0.000010, -0.000027, -0.000031, -0.000073, -0.000057, -0.000081, -0.000057, -0.000042, -0.000067, -0.000070, -0.000042, -0.000031, -0.000008, -0.000039, -0.000063, -0.000037, -0.000056, -0.000020, -0.000042, -0.000050, -0.000042, -0.000055, -0.000018, -0.000042, -0.000057, -0.000057, -0.000010, -0.000020, -0.000031, -0.000010, -0.000069, -0.000037, -0.000034, -0.000018, -0.000033, -0.000029, -0.000039, -0.000002, -0.000032, -0.000046, -0.000037, -0.000050, -0.000033, -0.000043, -0.000066, -0.000045, -0.000030, -0.000006, -0.000006, -0.000026, -0.000067, -0.000007, -0.000022, 0.000032, -0.000014, -0.000012, -0.000027, -0.000026, 0.000017, -0.000019, -0.000016, 0.000005, 0.000005, -0.000036, -0.000019, -0.000049, 0.000003, -0.000061, -0.000042, -0.000070, -0.000032, -0.000026, -0.000054, -0.000035, -0.000016, 0.000024, 0.000008, -0.000016, -0.000027, 0.000022, -0.000021, -0.000009, 0.000013, -0.000004, -0.000046, 0.000034, -0.000014, -0.000022, 0.000017, 0.000002, -0.000038, 0.000005, -0.000010, -0.000014, -0.000002, -0.000017, -0.000019, -0.000021, -0.000046, -0.000017, -0.000021, -0.000007, -0.000008, -0.000033, -0.000027, 0.000008, -0.000033, -0.000003, -0.000028, -0.000012, -0.000048, 0.000017, 0.000105, 0.000111, 0.000096, 0.000178, 0.000142, 0.000027, 0.000136, 0.000100, 0.000003, 0.000138, 0.000116, -0.000030, 0.000061, 0.000106, 0.000074, 0.000033, 0.000124, 0.000080, -0.000011, 0.000036, 0.000049, -0.000027, 0.000096, 0.000062, 0.000018, 0.000036, 0.000081, 0.000036, 0.000154, 0.000075, 0.000152, 0.000200, 0.000148, 0.000070, 0.000148, 0.000122, 0.000025, 0.000147, 0.000034, -0.000006, 0.000094, 0.000142, 0.000043, 0.000071, 0.000058, 0.000059, 0.000082, 0.000064, 0.000100, 0.000022, 0.000071, 0.000013, 0.000019, -0.000016, 0.000110, 0.000081, 0.000007, 0.000004, 0.000082, -0.000005, -0.000001, 0.000034, -0.000019, 0.000005, -0.000007, 0.000010, -0.000021, 0.000016, 0.000059, -0.000019, -0.000042, 0.000025, -0.000056, -0.000006, 0.000026, 0.000013, 0.000018, -0.000049, -0.000020, -0.000045, -0.000024, -0.000045, -0.000022, -0.000051, -0.000030, 0.000025, 0.000006, -0.000050, -0.000009, -0.000056, -0.000076, 0.000017, -0.000040, -0.000014, -0.000030, -0.000019, -0.000018, -0.000050, -0.000050, -0.000048, -0.000033, -0.000057, -0.000030, -0.000045, -0.000060, -0.000006, -0.000017, -0.000073, -0.000020, -0.000098, -0.000028, -0.000027, -0.000050, -0.000005, -0.000015, -0.000014, -0.000013, -0.000060, -0.000014, -0.000044, -0.000051, -0.000037, -0.000052, -0.000006, -0.000042, -0.000063, -0.000052, -0.000065, 0.000011, -0.000005, -0.000069, -0.000009, -0.000021, -0.000079, -0.000062, -0.000051, -0.000051, -0.000027, -0.000029, -0.000032, -0.000036, 0.000005, -0.000058, -0.000049, -0.000052, -0.000037, -0.000047, -0.000029, -0.000056, -0.000021, -0.000044, -0.000031, -0.000037, -0.000024, -0.000026, -0.000001, -0.000030, -0.000037, -0.000031, -0.000054, -0.000048, -0.000063, -0.000064, -0.000034, -0.000034, -0.000018, -0.000021, -0.000042, -0.000039, -0.000045, -0.000033, -0.000053, -0.000039, -0.000019, -0.000047, -0.000039, -0.000024, -0.000048, -0.000039, -0.000021, -0.000038, -0.000046, -0.000020, -0.000018, -0.000037, -0.000016, -0.000054, -0.000023, -0.000074, -0.000062, -0.000063, -0.000008, -0.000016, 0.000013, -0.000024, -0.000027, -0.000015, -0.000040, -0.000045, -0.000007, -0.000017, -0.000030, -0.000053, -0.000018, -0.000030, -0.000042, -0.000013, -0.000062, -0.000029, -0.000008, -0.000037, -0.000016, -0.000009, -0.000039, 0.000001, -0.000004, -0.000003, -0.000029, -0.000005, -0.000018, -0.000015, -0.000016, -0.000032, -0.000028, 0.000010, -0.000010, -0.000015, -0.000001, -0.000016, -0.000034, -0.000029, -0.000041, -0.000009, -0.000053, -0.000030, -0.000018, 0.000004, -0.000005, -0.000019, -0.000002, -0.000016, -0.000014, 0.000017, -0.000013, -0.000022, 0.000002, -0.000015, -0.000022, 0.000002, -0.000010, 0.000010, -0.000017, -0.000010, -0.000027, -0.000019, -0.000005, -0.000013, -0.000016, -0.000012, -0.000015, -0.000031, -0.000008, 0.000002, -0.000016, -0.000012, -0.000009, -0.000012, -0.000008, -0.000020, -0.000013, -0.000024, 0.000010, 0.000014, -0.000030, 0.000022, -0.000028, -0.000028, 0.000014, 0.000009, 0.000011, -0.000012, -0.000005, -0.000009, -0.000008, -0.000018, 0.000006, -0.000014, -0.000014, -0.000016, -0.000018, 0.000001, -0.000005, -0.000026, 0.000016, 0.000102, 0.000071, 0.000072, 0.000156, 0.000042, 0.000088, 0.000117, -0.000002, 0.000106, 0.000097, -0.000016, 0.000076, 0.000057, 0.000048, 0.000038, 0.000076, 0.000025, 0.000009, 0.000063, -0.000010, 0.000037, 0.000071, -0.000004, 0.000047, 0.000047, 0.000041, -0.000009, 0.000018, 0.000145, 0.000086, 0.000128, 0.000146, 0.000049, 0.000111, 0.000105, 0.000002, 0.000111, 0.000099, 0.000002, 0.000055, 0.000070, 0.000006, 0.000071, 0.000067, 0.000075, 0.000048, 0.000032, -0.000023, 0.000049, 0.000015, 0.000055, 0.000036, 0.000034, 0.000051, -0.000031, 0.000000, 0.000047, -0.000012, 0.000048, -0.000005, 0.000015, 0.000009, 0.000017, -0.000008, 0.000001, 0.000015, -0.000016, 0.000011, -0.000035, -0.000013, 0.000015, -0.000020, 0.000060, -0.000032, -0.000034, -0.000018, -0.000026, 0.000003, 0.000018, -0.000012, -0.000006, 0.000014, -0.000036, -0.000048, 0.000006, -0.000030, -0.000011, -0.000022, -0.000025, -0.000027, -0.000021, -0.000058, -0.000044, -0.000023, -0.000005, -0.000047, 0.000006, -0.000009, -0.000052, 0.000008, -0.000012, -0.000036, 0.000001, -0.000065, 0.000004, 0.000004, -0.000019, -0.000040, -0.000006, -0.000016, -0.000058, -0.000030, -0.000016, -0.000007, -0.000024, -0.000068, -0.000025, -0.000032, -0.000025, -0.000041, -0.000027, -0.000018, -0.000048, -0.000029, -0.000025, -0.000044, -0.000025, -0.000021, -0.000027, -0.000021, -0.000039, -0.000050, -0.000006, -0.000039, -0.000049, -0.000045, -0.000026, -0.000013, -0.000042, -0.000019, -0.000008, -0.000031, -0.000035, -0.000027, -0.000016, -0.000004, -0.000031, -0.000040, -0.000018, -0.000039, -0.000043, -0.000039, -0.000025, -0.000009, 0.000006, -0.000022, -0.000035, -0.000036, -0.000040, -0.000016, -0.000026, -0.000025, -0.000041, -0.000040, -0.000011, -0.000053, -0.000027, -0.000018, -0.000020, -0.000029, -0.000020, -0.000028, -0.000009, -0.000017, -0.000036, -0.000022, -0.000038, -0.000026, 0.000002, -0.000003, -0.000013, -0.000005, -0.000025, -0.000024, -0.000014, -0.000021, -0.000019, -0.000018, -0.000025, -0.000042, -0.000011, -0.000022, -0.000030, -0.000015, -0.000028, 0.000002, -0.000017, -0.000035, -0.000029, -0.000015, -0.000015, -0.000028, -0.000003, -0.000012, -0.000002, -0.000025, -0.000014, 0.000003, -0.000017, -0.000024, -0.000011, -0.000019, -0.000019, -0.000019, -0.000005, -0.000020, -0.000019, -0.000035, -0.000021, -0.000008, -0.000026, -0.000015, -0.000018, -0.000017, -0.000016, -0.000018, -0.000007, 0.000005, -0.000019, -0.000013, -0.000020, -0.000016, -0.000010, -0.000001, -0.000010, -0.000003, -0.000003, -0.000014, -0.000028, -0.000020, -0.000019, -0.000004, -0.000012, -0.000020, -0.000010, -0.000003, -0.000027, -0.000019, -0.000013, -0.000011, -0.000015, -0.000006, -0.000020, 0.000013, -0.000004, -0.000009, 0.000008, -0.000025, 0.000006, -0.000008, 0.000001, -0.000009, -0.000005, -0.000010, -0.000010, -0.000016, -0.000008, -0.000005, -0.000005, -0.000014, -0.000009, -0.000012, 0.000000, -0.000008, -0.000007, -0.000007, 0.000009, -0.000007, -0.000014, -0.000006, -0.000012, -0.000007, 0.000009, -0.000012, 0.000000, -0.000011, -0.000009, -0.000015, -0.000009, 0.000001, 0.000011, 0.000006, -0.000015, 0.000012, 0.000089, 0.000059, 0.000086, 0.000089, 0.000022, 0.000080, 0.000023, 0.000044, 0.000088, -0.000002, 0.000033, 0.000059, 0.000039, 0.000043, 0.000063, 0.000012, 0.000029, 0.000020, -0.000007, 0.000061, 0.000005, 0.000040, 0.000036, 0.000027, 0.000005, 0.000031, 0.000030, 0.000029, 0.000115, 0.000052, 0.000109, 0.000086, 0.000021, 0.000088, 0.000062, 0.000063, 0.000066, 0.000004, 0.000017, 0.000075, 0.000044, 0.000069, 0.000048, -0.000017, 0.000053, 0.000030, 0.000045, 0.000058, 0.000002, 0.000048, 0.000024, 0.000036, -0.000021, 0.000012, 0.000016, 0.000016, 0.000032, 0.000006, 0.000029, -0.000004, -0.000003, 0.000012, 0.000030, -0.000021, 0.000016, -0.000035, -0.000015, 0.000036, -0.000007, 0.000020, -0.000029, -0.000000, -0.000002, 0.000010, 0.000015, 0.000005, -0.000009, -0.000007, 0.000005, -0.000025, -0.000017, -0.000017, -0.000005, -0.000011, -0.000030, -0.000017, -0.000016, -0.000029, -0.000044, 0.000009, -0.000025, -0.000027, -0.000001, -0.000039, 0.000016, -0.000034, -0.000024, -0.000010, -0.000030, -0.000012, -0.000011, -0.000011, -0.000027, -0.000021, -0.000036, 0.000006, -0.000027, -0.000020, -0.000017, -0.000011, -0.000018, -0.000048, -0.000030, -0.000034, -0.000018, -0.000043, -0.000007, -0.000009, -0.000021, -0.000015, -0.000046, -0.000015, -0.000013, -0.000031, -0.000008, -0.000026, -0.000030, -0.000022, -0.000038, -0.000025, -0.000028, -0.000020, -0.000008, -0.000029, -0.000013, -0.000020, -0.000022, -0.000007, -0.000039, -0.000025, -0.000028, -0.000022, -0.000016, -0.000015, -0.000015, -0.000011, -0.000028, -0.000009, -0.000020, -0.000022, -0.000023, -0.000023, -0.000008, -0.000034, -0.000016, -0.000031, -0.000019, -0.000020, -0.000021, -0.000026, -0.000011, -0.000013, -0.000026, -0.000023, -0.000016, -0.000023, -0.000034, -0.000021, -0.000014, -0.000006, 0.000008, -0.000013, -0.000009, -0.000021, -0.000021, -0.000005, -0.000013, -0.000032, -0.000012, -0.000015, -0.000025, -0.000003, -0.000021, -0.000013, -0.000023, -0.000008, -0.000018, -0.000012, -0.000026, -0.000009, -0.000015, -0.000011, 0.000007, -0.000024, -0.000012, -0.000015, -0.000019, -0.000010, -0.000015, -0.000004, -0.000016, -0.000015, -0.000019, -0.000008, -0.000023, -0.000011, -0.000008, 0.000001, 0.000007, -0.000010, -0.000012, -0.000009, -0.000009, -0.000006, -0.000016, -0.000010, -0.000006, -0.000004, 0.000004, 0.000006, -0.000011, -0.000014, -0.000012, -0.000013, -0.000018, -0.000011, -0.000009, -0.000013, -0.000016, 0.000001, -0.000008, -0.000009, -0.000008, -0.000010, -0.000014, -0.000001, -0.000010, -0.000003, 0.000001, -0.000014, -0.000003, -0.000017, -0.000012, -0.000009, -0.000006, -0.000003, -0.000002, -0.000023, -0.000011, -0.000012, 0.000000, -0.000001, -0.000008, -0.000009, -0.000011, -0.000003, -0.000011, -0.000000, -0.000015, -0.000006, -0.000004, -0.000004, -0.000003, -0.000012, -0.000008, -0.000007, -0.000001, -0.000011, -0.000008, -0.000002, -0.000003, 0.000003, -0.000003, -0.000005, -0.000011, -0.000006, 0.000004, -0.000008, -0.000006, -0.000000, -0.000012, 0.000003, -0.000006, 0.000001, -0.000005, -0.000012, -0.000006, -0.000002, -0.000006, -0.000010, -0.000003, 0.000001, -0.000008, 0.000014, 0.000066, 0.000040, 0.000089, 0.000047, 0.000041, 0.000051, 0.000016, 0.000078, 0.000006, 0.000029, 0.000045, 0.000038, 0.000034, 0.000049, 0.000004, 0.000033, 0.000005, 0.000023, 0.000026, 0.000002, 0.000026, 0.000029, 0.000000, 0.000013, 0.000034, 0.000009, 0.000017, 0.000026, 0.000079, 0.000041, 0.000076, 0.000050, 0.000069, 0.000042, 0.000034, 0.000041, 0.000019, 0.000038, 0.000059, 0.000028, 0.000020, 0.000041, 0.000007, 0.000053, -0.000011, 0.000043, 0.000016, 0.000007, 0.000023, 0.000021, 0.000004, 0.000014, 0.000017, 0.000015, 0.000013, -0.000008, 0.000008, 0.000016, 0.000017, 0.000023, 0.000008, -0.000010, 0.000007, -0.000025, 0.000010, -0.000014, 0.000017, -0.000009, 0.000011, -0.000009, 0.000020, -0.000009, -0.000008, 0.000006, -0.000013, -0.000004, -0.000008, -0.000002, -0.000010, -0.000003, -0.000009, -0.000001, -0.000022, -0.000005, -0.000016, -0.000002, -0.000018, -0.000004, -0.000013, -0.000013, -0.000001, -0.000017, -0.000024, -0.000010, -0.000001, -0.000006, -0.000019, -0.000001, -0.000009, -0.000031, 0.000005, -0.000013, 0.000004, -0.000027, -0.000008, -0.000028, -0.000022, -0.000030, -0.000011, -0.000024, -0.000007, -0.000020, -0.000022, -0.000013, -0.000022, -0.000002, -0.000021, -0.000021, 0.000002, -0.000023, -0.000009, -0.000006, -0.000021, -0.000016, -0.000013, -0.000003, -0.000023, -0.000017, -0.000022, -0.000011, -0.000019, -0.000015, -0.000013, -0.000009, -0.000018, -0.000015, -0.000019, -0.000011, -0.000014, -0.000020, -0.000024, -0.000017, -0.000020, -0.000024, -0.000019, -0.000017, -0.000020, -0.000019, -0.000009, -0.000017, -0.000016, -0.000018, -0.000010, -0.000021, -0.000010, -0.000019, -0.000012, -0.000022, -0.000006, -0.000005, -0.000005, -0.000010, -0.000011, -0.000022, -0.000011, -0.000018, -0.000008, -0.000024, -0.000012, -0.000008, -0.000018, -0.000012, -0.000010, -0.000012, -0.000010, -0.000013, -0.000006, -0.000012, -0.000018, -0.000002, -0.000004, -0.000007, -0.000014, -0.000009, -0.000017, -0.000003, -0.000006, -0.000009, -0.000016, -0.000008, -0.000012, -0.000015, -0.000016, -0.000017, -0.000013, -0.000011, 0.000000, -0.000007, -0.000010, -0.000005, -0.000015, -0.000009, -0.000013, -0.000004, -0.000007, -0.000001, -0.000001, -0.000007, -0.000009, -0.000005, -0.000006, -0.000017, -0.000005, -0.000005, -0.000010, -0.000006, -0.000008, -0.000009, -0.000011, -0.000004, -0.000004, -0.000008, -0.000004, -0.000002, -0.000006, -0.000009, -0.000008, -0.000011, -0.000008, -0.000003, -0.000004, -0.000007, -0.000007, -0.000010, -0.000004, -0.000011, -0.000003, -0.000010, -0.000008, -0.000004, 0.000002, -0.000008, -0.000006, -0.000008, -0.000003, -0.000008, -0.000006, -0.000004, -0.000007, -0.000002, -0.000004, 0.000001, -0.000004, -0.000002, -0.000005, -0.000006, -0.000003, -0.000009, -0.000001, -0.000008, 0.000003, -0.000004, -0.000000, 0.000002, -0.000007, -0.000002, -0.000008, -0.000002, -0.000004, -0.000005, -0.000004, -0.000007, -0.000002, -0.000003, -0.000004, -0.000008, -0.000008, 0.000003, 0.000002, -0.000002, 0.000001, 0.000000, -0.000000, 0.000001, -0.000002, -0.000003, -0.000007, 0.000001, -0.000006, 0.000000, -0.000004, -0.000001, -0.000007, 0.000020, 0.000050, 0.000033, 0.000064, 0.000017, 0.000050, 0.000020, 0.000042, 0.000027, 0.000005, 0.000044, 0.000026, 0.000034, 0.000029, 0.000009, 0.000025, 0.000004, 0.000024, 0.000013, 0.000010, 0.000025, 0.000010, 0.000003, 0.000016, 0.000011, 0.000007, 0.000010, 0.000002, 0.000033, 0.000043, 0.000057, 0.000063, 0.000016, 0.000042, 0.000012, 0.000053, 0.000035, 0.000020, 0.000033, 0.000028, 0.000028, 0.000043, 0.000006, 0.000033, 0.000007, 0.000023, 0.000000, 0.000011, 0.000024, 0.000006, 0.000007, 0.000017, 0.000009, 0.000002, 0.000016, 0.000014, 0.000010, -0.000002, 0.000020, 0.000004, 0.000001, -0.000014, -0.000001, 0.000010, 0.000015, 0.000003, -0.000008, -0.000007, -0.000002, 0.000018, -0.000007, 0.000009, -0.000009, -0.000005, -0.000016, -0.000008, -0.000002, -0.000010, -0.000003, 0.000002, -0.000019, -0.000004, -0.000008, 0.000003, -0.000012, -0.000008, -0.000003, -0.000013, -0.000010, -0.000015, -0.000012, -0.000005, 0.000001, -0.000007, -0.000011, -0.000005, 0.000003, -0.000005, -0.000012, -0.000012, -0.000014, -0.000020, -0.000019, -0.000007, -0.000008, -0.000018, -0.000009, -0.000018, -0.000012, -0.000018, -0.000013, -0.000012, -0.000019, -0.000004, -0.000012, -0.000011, -0.000011, -0.000011, -0.000014, -0.000007, -0.000012, -0.000016, -0.000019, -0.000005, -0.000005, -0.000012, -0.000009, -0.000015, -0.000014, -0.000010, -0.000005, -0.000003, -0.000014, -0.000007, -0.000014, -0.000015, -0.000013, -0.000012, -0.000013, -0.000011, -0.000017, -0.000007, -0.000014, -0.000009, -0.000010, -0.000003, -0.000013, -0.000013, -0.000014, -0.000014, -0.000008, -0.000006, -0.000004, -0.000005, -0.000010, -0.000015, -0.000014, -0.000008, -0.000011, -0.000009, -0.000011, -0.000010, -0.000009, -0.000008, -0.000010, -0.000010, -0.000004, -0.000011, -0.000004, -0.000007, -0.000011, -0.000006, -0.000009, -0.000011, -0.000012, -0.000008, -0.000009, -0.000009, -0.000007, -0.000006, -0.000011, -0.000008, -0.000007, -0.000008, -0.000010, -0.000005, -0.000004, -0.000008, -0.000013, -0.000008, -0.000009, -0.000003, -0.000006, -0.000002, -0.000012, -0.000003, -0.000001, -0.000007, -0.000006, -0.000010, -0.000008, -0.000010, -0.000006, -0.000009, -0.000009, -0.000005, -0.000005, -0.000003, -0.000007, -0.000005, -0.000004, -0.000008, -0.000006, -0.000009, -0.000008, -0.000008, -0.000007, -0.000007, -0.000004, -0.000003, 0.000003, -0.000008, -0.000001, -0.000004, -0.000002, -0.000005, -0.000006, -0.000005, -0.000007, -0.000004, -0.000007, -0.000005, -0.000002, -0.000005, -0.000003, -0.000008, 0.000002, 0.000001, 0.000001, -0.000005, -0.000004, -0.000006, -0.000003, -0.000003, -0.000004, -0.000005, -0.000004, -0.000002, -0.000000, -0.000004, -0.000002, -0.000005, 0.000000, -0.000001, -0.000004, -0.000005, 0.000000, -0.000001, -0.000004, -0.000004, -0.000004, -0.000003, -0.000006, -0.000003, -0.000000, -0.000001, -0.000003, -0.000007, 0.000000, -0.000004, 0.000000, -0.000002, 0.000001, -0.000004, 0.000002, -0.000005, -0.000002, -0.000004, 0.000002, -0.000004, 0.000003, -0.000006, -0.000001, -0.000003, -0.000002, -0.000006, -0.000002, -0.000003, 0.000005, -0.000007, 0.000002, -0.000004, 0.000005, -0.000009, 0.000026, 0.000028, 0.000038, 0.000033, 0.000019, 0.000029, 0.000010, 0.000036, 0.000001, 0.000028, 0.000021, 0.000021, 0.000026, 0.000006, 0.000019, 0.000005, 0.000025, 0.000005, 0.000015, 0.000019, 0.000005, 0.000010, 0.000013, 0.000010, 0.000011, 0.000002, 0.000011, -0.000009, 0.000040, 0.000030, 0.000037, 0.000026, 0.000020, 0.000036, 0.000027, 0.000035, 0.000000, 0.000032, 0.000026, 0.000028, 0.000022, 0.000018, 0.000008, 0.000007, 0.000013, 0.000008, 0.000010, 0.000024, 0.000006, 0.000007, 0.000009, 0.000012, 0.000006, 0.000008, 0.000002, 0.000007, 0.000011, 0.000000, -0.000002, -0.000002, 0.000009, 0.000008, 0.000011, -0.000010, -0.000000, 0.000003, 0.000007, 0.000002, -0.000001, 0.000002, -0.000013, 0.000002, -0.000008, 0.000002, 0.000003, -0.000004, -0.000001, -0.000014, 0.000003, -0.000007, 0.000004, -0.000006, -0.000009, -0.000000, -0.000011, -0.000012, -0.000005, -0.000000, -0.000008, -0.000003, -0.000001, -0.000008, -0.000007, -0.000004, -0.000014, -0.000007, -0.000009, -0.000008, -0.000009, -0.000006, -0.000016, -0.000007, -0.000010, -0.000009, -0.000013, -0.000006, -0.000008, -0.000003, -0.000002, -0.000010, -0.000005, -0.000006, -0.000012, -0.000010, -0.000007, -0.000014, -0.000008, -0.000002, -0.000001, -0.000008, -0.000010, -0.000013, -0.000011, -0.000010, -0.000009, -0.000005, -0.000009, -0.000014, -0.000008, -0.000009, -0.000010, -0.000011, -0.000009, -0.000014, -0.000007, -0.000006, -0.000006, -0.000007, -0.000006, -0.000012, -0.000010, -0.000006, -0.000012, -0.000005, -0.000006, -0.000004, -0.000009, -0.000010, -0.000011, -0.000006, -0.000004, -0.000010, -0.000007, -0.000006, -0.000010, -0.000007, -0.000005, -0.000008, -0.000005, -0.000004, -0.000009, -0.000009, -0.000007, -0.000008, -0.000008, -0.000005, -0.000011, -0.000008, -0.000003, -0.000006, -0.000008, -0.000009, -0.000008, -0.000009, -0.000006, -0.000008, -0.000008, -0.000002, -0.000010, -0.000003, -0.000003, -0.000002, -0.000007, -0.000006, -0.000007, -0.000005, -0.000002, -0.000009, -0.000007, -0.000006, -0.000007, -0.000008, -0.000003, -0.000005, -0.000003, -0.000002, -0.000007, -0.000003, -0.000004, -0.000002, -0.000007, -0.000002, -0.000005, -0.000002, -0.000006, -0.000008, -0.000005, -0.000001, -0.000002, -0.000003, -0.000007, -0.000003, -0.000003, -0.000003, -0.000008, -0.000004, -0.000004, -0.000004, -0.000004, -0.000003, -0.000003, -0.000006, -0.000004, 0.000000, -0.000000, -0.000001, -0.000004, -0.000006, -0.000003, -0.000003, -0.000002, -0.000004, -0.000003, -0.000002, -0.000001, -0.000002, -0.000002, -0.000001, -0.000001, -0.000002, -0.000002, -0.000003, -0.000001, -0.000002, -0.000003, -0.000005, -0.000001, -0.000004, -0.000004, -0.000001, 0.000003, -0.000001, -0.000000, -0.000003, -0.000001, -0.000002, 0.000001, -0.000002, -0.000001, -0.000001, 0.000000, -0.000003, -0.000000, -0.000002, 0.000000, -0.000003, -0.000001, -0.000003, -0.000001, -0.000004, -0.000001, -0.000003, -0.000001, -0.000001, -0.000001, -0.000005, -0.000001, -0.000004, -0.000001, -0.000000, -0.000002, 0.000001, 0.000001, -0.000002, -0.000001, -0.000002, 0.000001, -0.000003, 0.000003, -0.000005, 0.000004, -0.000004, 0.000005, -0.000006, 0.000023, 0.000018, 0.000031, 0.000021, 0.000023, 0.000014, 0.000023, 0.000011, 0.000008, 0.000023, 0.000011, 0.000022, 0.000007, 0.000015, 0.000004, 0.000016, 0.000003, 0.000013, 0.000013, -0.000000, 0.000010, 0.000006, 0.000009, 0.000004, 0.000012, -0.000001, 0.000009, 0.000003, 0.000024, 0.000019, 0.000030, 0.000024, 0.000029, 0.000011, 0.000019, 0.000013, 0.000017, 0.000022, 0.000015, 0.000022, 0.000001, 0.000013, 0.000003, 0.000018, 0.000003, 0.000019, 0.000005, 0.000004, 0.000012, 0.000008, 0.000006, 0.000004, 0.000008, -0.000000, 0.000007, 0.000004, 0.000005, -0.000003, 0.000004, 0.000007, 0.000003, 0.000001, -0.000004, 0.000002, 0.000003, -0.000001, -0.000002, 0.000002, -0.000009, -0.000001, -0.000004, 0.000001, -0.000003, 0.000003, -0.000009, 0.000001, 0.000001, -0.000004, -0.000004, -0.000004, -0.000005, -0.000003, -0.000006, -0.000006, 0.000000, 0.000002, -0.000000, -0.000008, -0.000001, -0.000006, -0.000005, -0.000004, -0.000003, -0.000010, -0.000005, -0.000006, -0.000008, -0.000005, -0.000008, -0.000008, -0.000009, -0.000000, -0.000006, -0.000001, -0.000005, -0.000007, -0.000005, -0.000009, -0.000010, -0.000002, -0.000004, -0.000003, -0.000001, -0.000003, -0.000009, -0.000004, -0.000009, -0.000011, -0.000004, -0.000005, -0.000007, -0.000007, -0.000006, -0.000006, -0.000007, -0.000009, -0.000008, -0.000008, -0.000004, -0.000005, -0.000008, -0.000006, -0.000007, -0.000006, -0.000007, -0.000007, -0.000009, -0.000004, -0.000004, -0.000005, -0.000005, -0.000004, -0.000005, -0.000005, -0.000007, -0.000006, -0.000006, -0.000007, -0.000006, -0.000005, -0.000004, -0.000006, -0.000006, -0.000006, -0.000009, -0.000004, -0.000005, -0.000007, -0.000005, -0.000007, -0.000005, -0.000006, -0.000004, -0.000005, -0.000004, -0.000006, -0.000007, -0.000006, -0.000005, -0.000003, -0.000005, -0.000002, -0.000006, -0.000005, -0.000004, -0.000004, -0.000004, -0.000002, -0.000005, -0.000006, -0.000005, -0.000007, -0.000004, -0.000004, -0.000005, -0.000003, -0.000003, -0.000003, -0.000004, -0.000005, -0.000004, -0.000002, -0.000002, -0.000005, -0.000005, -0.000005, -0.000002, -0.000003, -0.000001, -0.000005, -0.000003, -0.000003, -0.000003, -0.000005, -0.000003, -0.000004, -0.000004, -0.000003, -0.000004, -0.000005, -0.000001, -0.000000, -0.000001, -0.000002, -0.000003, -0.000003, -0.000003, -0.000004, -0.000003, -0.000002, -0.000002, -0.000003, -0.000001, -0.000003, -0.000002, -0.000002, -0.000002, -0.000001, 0.000000, -0.000003, -0.000003, -0.000002, -0.000001, -0.000003, -0.000004, -0.000003, -0.000003, -0.000000, -0.000001, -0.000001, -0.000002, -0.000003, -0.000002, -0.000002, -0.000003, -0.000002, -0.000002, 0.000000, -0.000002, -0.000001, -0.000002, -0.000000, -0.000000, -0.000002, -0.000002, 0.000000, 0.000000, -0.000001, -0.000002, -0.000001, -0.000002, -0.000002, -0.000002, -0.000000, -0.000001, -0.000001, -0.000000, -0.000001, 0.000000, -0.000000, 0.000000, -0.000001, 0.000000, -0.000002, 0.000001, -0.000001, 0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, 0.000000, 0.000001, -0.000001, -0.000002, -0.000001, -0.000002, -0.000001, -0.000003, 0.000002, -0.000002, 0.000002, -0.000002, 0.000018, 0.000016, 0.000023, 0.000012, 0.000018, 0.000010, 0.000018, 0.000003, 0.000017, 0.000010, 0.000018, 0.000008, 0.000010, 0.000005, 0.000011, 0.000005, 0.000010, 0.000007, 0.000001, 0.000007, 0.000007, 0.000005, 0.000006, 0.000004, 0.000002, 0.000010, 0.000002, 0.000001, 0.000015, 0.000022, 0.000027, 0.000013, 0.000014, 0.000013, 0.000021, 0.000004, 0.000020, 0.000010, 0.000013, 0.000006, 0.000011, 0.000006, 0.000012, 0.000003, 0.000007, 0.000009, 0.000002, 0.000008, 0.000003, 0.000007, 0.000002, 0.000001, 0.000001, 0.000010, 0.000002, 0.000002, -0.000003, 0.000008, 0.000005, -0.000000, -0.000001, 0.000003, -0.000001, 0.000003, 0.000000, -0.000001, -0.000001, -0.000002, -0.000001, 0.000001, -0.000000, -0.000001, -0.000001, -0.000001, 0.000001, -0.000004, -0.000003, 0.000001, -0.000005, -0.000005, 0.000002, 0.000000, -0.000003, -0.000002, -0.000007, -0.000002, -0.000002, -0.000002, -0.000005, -0.000005, -0.000007, -0.000002, -0.000007, -0.000002, -0.000006, -0.000004, -0.000004, -0.000003, -0.000004, -0.000001, -0.000006, -0.000004, -0.000006, -0.000007, -0.000002, -0.000001, -0.000004, -0.000003, -0.000002, -0.000004, -0.000006, -0.000007, -0.000006, -0.000006, -0.000002, -0.000004, -0.000006, -0.000005, -0.000005, -0.000006, -0.000005, -0.000005, -0.000006, -0.000003, -0.000006, -0.000004, -0.000003, -0.000004, -0.000005, -0.000006, -0.000007, -0.000004, -0.000003, -0.000002, -0.000004, -0.000006, -0.000004, -0.000004, -0.000007, -0.000003, -0.000005, -0.000004, -0.000002, -0.000005, -0.000005, -0.000005, -0.000006, -0.000004, -0.000005, -0.000004, -0.000004, -0.000006, -0.000006, -0.000003, -0.000004, -0.000004, -0.000004, -0.000005, -0.000005, -0.000005, -0.000003, -0.000001, -0.000002, -0.000004, -0.000002, -0.000004, -0.000003, -0.000004, -0.000003, -0.000003, -0.000003, -0.000004, -0.000003, -0.000005, -0.000002, -0.000003, -0.000003, -0.000002, -0.000005, -0.000004, -0.000004, -0.000002, -0.000004, -0.000004, -0.000004, -0.000002, -0.000004, -0.000002, -0.000002, -0.000002, -0.000003, -0.000003, -0.000003, -0.000001, -0.000004, -0.000003, -0.000003, -0.000004, -0.000004, -0.000001, -0.000001, -0.000001, -0.000002, -0.000002, -0.000000, -0.000003, -0.000003, -0.000002, -0.000001, -0.000001, -0.000003, -0.000002, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000003, -0.000003, -0.000002, -0.000002, -0.000003, -0.000002, -0.000002, -0.000002, -0.000000, -0.000000, -0.000001, -0.000002, -0.000002, -0.000002, -0.000001, -0.000002, -0.000001, -0.000001, -0.000000, -0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, 0.000000, -0.000001, -0.000000, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, 0.000000, -0.000001, -0.000001, -0.000001, -0.000002, -0.000000, 0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, 0.000001, -0.000001, 0.000000, -0.000002, -0.000000, -0.000001, 0.000000, -0.000001, 0.000001, -0.000002, 0.000000, -0.000001, 0.000002, -0.000002, 0.000001, -0.000002, 0.000003, 0.000011, 0.000017, 0.000013, 0.000012, 0.000008, 0.000012, 0.000006, 0.000008, 0.000011, 0.000010, 0.000009, 0.000007, 0.000005, 0.000007, 0.000004, 0.000007, 0.000007, 0.000002, 0.000007, 0.000004, 0.000004, 0.000005, 0.000003, 0.000005, 0.000005, 0.000001, -0.000001, 0.000003, 0.000017, 0.000015, 0.000013, 0.000012, 0.000014, 0.000012, 0.000010, 0.000006, 0.000009, 0.000008, 0.000011, 0.000005, 0.000008, 0.000004, 0.000006, 0.000008, 0.000007, 0.000001, 0.000008, 0.000004, 0.000002, 0.000003, 0.000004, 0.000004, 0.000006, -0.000002, 0.000001, 0.000002, 0.000004, 0.000000, -0.000001, -0.000001, 0.000004, -0.000001, 0.000004, -0.000002, -0.000002, -0.000001, 0.000001, -0.000001, 0.000003, -0.000005, 0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000003, -0.000000, -0.000002, -0.000001, -0.000003, -0.000002, -0.000002, 0.000001, -0.000002, -0.000003, -0.000003, -0.000004, -0.000003, -0.000002, -0.000005, -0.000002, -0.000002, -0.000004, -0.000001, -0.000002, -0.000002, -0.000003, -0.000005, -0.000001, -0.000003, -0.000001, -0.000003, -0.000005, -0.000001, -0.000002, -0.000006, -0.000002, -0.000005, -0.000004, -0.000003, -0.000002, -0.000004, -0.000003, -0.000004, -0.000004, -0.000004, -0.000004, -0.000005, -0.000004, -0.000003, -0.000002, -0.000003, -0.000004, -0.000004, -0.000005, -0.000005, -0.000002, -0.000001, -0.000003, -0.000003, -0.000003, -0.000003, -0.000004, -0.000004, -0.000002, -0.000003, -0.000003, -0.000003, -0.000004, -0.000003, -0.000003, -0.000004, -0.000002, -0.000003, -0.000004, -0.000004, -0.000003, -0.000002, -0.000003, -0.000004, -0.000004, -0.000005, -0.000003, -0.000002, -0.000001, -0.000003, -0.000003, -0.000002, -0.000003, -0.000004, -0.000003, -0.000003, -0.000002, -0.000003, -0.000003, -0.000003, -0.000003, -0.000002, -0.000003, -0.000003, -0.000002, -0.000003, -0.000002, -0.000002, -0.000003, -0.000001, -0.000002, -0.000003, -0.000003, -0.000002, -0.000002, -0.000002, -0.000003, -0.000002, -0.000002, -0.000002, -0.000003, -0.000003, -0.000003, -0.000001, -0.000001, -0.000002, -0.000002, -0.000002, -0.000001, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000001, -0.000002, -0.000001, -0.000001, -0.000002, -0.000001, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000001, -0.000001, -0.000001, -0.000002, -0.000000, -0.000000, -0.000001, -0.000002, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000002, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000000, -0.000001, -0.000001, -0.000002, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, 0.000001, 0.000001, 0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, 0.000000, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000000, 0.000001, 0.000000, -0.000001, 0.000002, 0.000010, 0.000011, 0.000011, 0.000008, 0.000008, 0.000008, 0.000004, 0.000009, 0.000009, 0.000007, 0.000006, 0.000003, 0.000006, 0.000004, 0.000005, 0.000006, 0.000001, 0.000005, 0.000004, 0.000003, 0.000003, 0.000002, 0.000004, 0.000004, 0.000000, 0.000001, 0.000003, 0.000005, 0.000010, 0.000011, 0.000013, 0.000009, 0.000010, 0.000009, 0.000004, 0.000007, 0.000008, 0.000007, 0.000006, 0.000002, 0.000006, 0.000005, 0.000005, 0.000004, 0.000002, 0.000004, 0.000001, 0.000004, 0.000004, 0.000002, 0.000003, 0.000002, 0.000001, -0.000000, 0.000004, 0.000002, -0.000001, 0.000001, 0.000002, 0.000002, 0.000001, -0.000001, -0.000001, -0.000000, 0.000001, 0.000001, -0.000000, -0.000002, -0.000001, 0.000001, -0.000001, -0.000001, 0.000002, -0.000002, -0.000002, -0.000001, -0.000002, -0.000000, -0.000001, -0.000001, -0.000001, -0.000002, -0.000002, -0.000004, -0.000002, -0.000003, -0.000002, -0.000001, -0.000002, -0.000001, -0.000002, -0.000003, -0.000002, -0.000003, -0.000000, -0.000002, -0.000003, -0.000002, -0.000003, -0.000002, -0.000002, -0.000002, -0.000003, -0.000003, -0.000003, -0.000002, -0.000002, -0.000003, -0.000004, -0.000002, -0.000003, -0.000003, -0.000004, -0.000002, -0.000001, -0.000002, -0.000002, -0.000003, -0.000003, -0.000004, -0.000003, -0.000002, -0.000002, -0.000001, -0.000003, -0.000003, -0.000002, -0.000003, -0.000001, -0.000003, -0.000003, -0.000002, -0.000003, -0.000002, -0.000003, -0.000002, -0.000003, -0.000003, -0.000003, -0.000003, -0.000003, -0.000002, -0.000003, -0.000003, -0.000003, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000001, -0.000002, -0.000001, -0.000002, -0.000002, -0.000001, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000001, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000001, -0.000001, -0.000001, -0.000002, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000002, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000002, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000002, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000001, -0.000001, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, 0.000000, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000000, -0.000001, 0.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000000, 0.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, 0.000000, -0.000000, 0.000000, -0.000001, 0.000002, 0.000008, 0.000009, 0.000008, 0.000006, 0.000007, 0.000005, 0.000004, 0.000007, 0.000005, 0.000005, 0.000003, 0.000004, 0.000004, 0.000003, 0.000005, 0.000001, 0.000004, 0.000003, 0.000003, 0.000003, 0.000002, 0.000004, 0.000002, 0.000000, 0.000001, 0.000003, 0.000001, 0.000002, 0.000009, 0.000010, 0.000008, 0.000007, 0.000006, 0.000004, 0.000005, 0.000007, 0.000006, 0.000004, 0.000003, 0.000005, 0.000003, 0.000003, 0.000005, 0.000001, 0.000003, 0.000003, 0.000003, 0.000003, 0.000001, 0.000004, 0.000002, 0.000000, 0.000001, 0.000002, 0.000001, -0.000001, 0.000002, 0.000001, 0.000000, 0.000001, -0.000001, 0.000001, -0.000000, 0.000000, -0.000000, -0.000002, 0.000001, -0.000001, -0.000000, 0.000001, -0.000001, -0.000002, -0.000002, -0.000000, 0.000000, -0.000001, -0.000002, -0.000000, -0.000002, -0.000001, -0.000002, -0.000002, -0.000000, -0.000002, -0.000000, -0.000001, -0.000002, -0.000002, -0.000002, -0.000000, -0.000001, -0.000001, -0.000003, -0.000001, -0.000001, -0.000003, -0.000000, -0.000002, -0.000002, -0.000002, -0.000003, -0.000001, -0.000002, -0.000002, -0.000002, -0.000002, -0.000003, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000001, -0.000003, -0.000002, -0.000002, -0.000002, -0.000001, -0.000001, -0.000002, -0.000002, -0.000001, -0.000002, -0.000002, -0.000002, -0.000003, -0.000002, -0.000002, -0.000002, -0.000002, -0.000003, -0.000002, -0.000001, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000001, -0.000002, -0.000002, -0.000001, -0.000002, -0.000002, -0.000001, -0.000002, -0.000002, -0.000002, -0.000001, -0.000001, -0.000002, -0.000001, -0.000002, -0.000002, -0.000002, -0.000001, -0.000002, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000002, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000002, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000001, -0.000001, 0.000002, 0.000007, 0.000007, 0.000006, 0.000004, 0.000005, 0.000003, 0.000004, 0.000005, 0.000004, 0.000003, 0.000002, 0.000003, 0.000002, 0.000004, 0.000001, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, 0.000003, 0.000001, -0.000000, 0.000002, 0.000003, 0.000001, -0.000000, 0.000004, 0.000007, 0.000007, 0.000006, 0.000004, 0.000005, 0.000003, 0.000005, 0.000004, 0.000004, 0.000003, 0.000002, 0.000003, 0.000003, 0.000003, 0.000001, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, -0.000001, 0.000001, 0.000002, 0.000000, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, -0.000001, 0.000000, -0.000000, 0.000000, -0.000000, -0.000001, 0.000001, 0.000000, -0.000001, -0.000000, -0.000002, -0.000000, 0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, 0.000000, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000002, -0.000000, -0.000001, -0.000001, -0.000002, -0.000002, -0.000001, -0.000001, -0.000001, -0.000002, -0.000002, -0.000002, -0.000002, -0.000002, -0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000002, -0.000002, -0.000001, -0.000001, -0.000001, -0.000002, -0.000002, -0.000001, -0.000001, -0.000001, -0.000002, -0.000002, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, 0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000002, 0.000005, 0.000005, 0.000004, 0.000003, 0.000004, 0.000002, 0.000004, 0.000003, 0.000003, 0.000001, 0.000003, 0.000002, 0.000003, 0.000002, 0.000002, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000000, 0.000001, 0.000002, -0.000000, 0.000001, 0.000001, 0.000003, 0.000006, 0.000004, 0.000004, 0.000004, 0.000004, 0.000002, 0.000003, 0.000004, 0.000003, 0.000001, 0.000003, 0.000001, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002, 0.000000, 0.000000, 0.000001, 0.000002, -0.000000, 0.000001, 0.000001, 0.000001, 0.000001, -0.000001, 0.000000, -0.000000, 0.000001, -0.000001, 0.000000, 0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, 0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000002, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000002, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000002, 0.000004, 0.000004, 0.000003, 0.000003, 0.000002, 0.000002, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, 0.000002, 0.000001, 0.000001, 0.000001, 0.000002, 0.000001, -0.000000, 0.000001, 0.000001, 0.000000, 0.000001, 0.000001, 0.000001, 0.000003, 0.000003, 0.000004, 0.000003, 0.000003, 0.000002, 0.000002, 0.000003, 0.000002, 0.000002, 0.000001, 0.000001, 0.000002, 0.000002, 0.000001, 0.000002, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, -0.000000, 0.000001, 0.000001, 0.000000, 0.000001, 0.000001, 0.000001, 0.000000, -0.000001, 0.000000, 0.000000, 0.000000, 0.000000, -0.000000, 0.000000, -0.000001, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000002, 0.000003, 0.000003, 0.000002, 0.000002, 0.000001, 0.000002, 0.000002, 0.000002, 0.000001, 0.000001, 0.000001, 0.000002, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, 0.000001, 0.000001, 0.000000, 0.000001, 0.000001, 0.000001, 0.000000, 0.000002, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, 0.000001, 0.000001, 0.000002, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, 0.000001, 0.000000, 0.000000, 0.000001, 0.000001, 0.000000, 0.000001, 0.000001, 0.000000, -0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000000, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000001, -0.000001, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, 0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, 0.000002, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, 0.000001, 0.000001, 0.000000, 0.000001, 0.000001, 0.000001, 0.000000, 0.000000, 0.000002, 0.000002, 0.000002, 0.000002, 0.000002, 0.000001, 0.000002, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, 0.000001, 0.000001, 0.000001, -0.000000, 0.000001, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -0.000000, 0.000001, 0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000001, -0.000001, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000001, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, 0.000001, 0.000002, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, 0.000000, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000002, 0.000002, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000000, 0.000001, 0.000000, 0.000000, 0.000000, 0.000001, -0.000000, 0.000000, 0.000001, -0.000000, 0.000000, 0.000000, 0.000000, -0.000000, 0.000000, 0.000000, 0.000000, -0.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000,
    };

    struct convolution_filter
    {
        static constexpr size_t size = 8192;
        static constexpr size_t buffer_size = 2 * size;
        static constexpr std::size_t mask = size - 1;

        std::array<real, buffer_size> history = {};
        size_t head = 0;

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
                y += g_impulse[i] * history[i + head];
            }
            head = (head - 1) & mask;
            return y;
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
        /* AUDIO_Y       */ 6,
        /* PIPE_CELLS    */ 256,
        /* PIPE_SUBSTEPS */ 10,
        inline_pistons,
        basic_cams,
        basic_sparkplugs>
    {
        inline4()
        {
            this->lumped_drag_torque_n_m = 29.2_r;
            this->limiter.max_angular_velocity_r_per_s = 900.0_r;
            this->limiter.limit_time_s = 0.033_r;
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
                this->inlet_cam.engage_theta_r[i]  = theta0_r + g_otto_intake_cycle_r - g_pi_r / 8.0_r;
                this->sparkplugs.engage_theta_r[i] = theta0_r + g_otto_combustion_cycle_r - g_pi_r / 6.2_r;
                this->outlet_cam.engage_theta_r[i] = theta0_r + g_otto_exhaust_cycle_r + 0.6_r;
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
                    0.00110_r, /* Runner             -> Chamber1         */
                    0.00200_r, /* Chamber1           -> Chamber2         */
                    0.00250_r, /* Chamber2           -> Atmospheric Sink */
                };
            }
            /*                             Atmospheric Source    Intake   Throttle  Runner    Piston     Runner    Chamber1  Chamber2   Atmospheric Sink    */
            flows[0].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.0003_r, 0.0003_r, 0.0003_r, g_resevoir_volume_m3 };
            flows[1].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.0003_r, 0.0003_r, 0.0003_r, g_resevoir_volume_m3 };
            flows[2].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.0003_r, 0.0003_r, 0.0003_r, g_resevoir_volume_m3 };
            flows[3].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.0003_r, 0.0003_r, 0.0003_r, g_resevoir_volume_m3 };
            this->throttle.table = {
                0.00100_r,
                0.02500_r,
                0.25000_r,
                1.00000_r,
            };
            this->pipe.piston_connect_ratio = { 0.0_r, 0.38_r, 0.17_r, 0.63_r };
            this->pipe.mic_position0_ratio = 0.8_r;
            this->pipe.mic_position1_ratio = 0.9_r;
            this->pipe.length_m = 1.0_r;
            this->dc.set_cutoff_frequency(5.0_r);
            this->gain.ratio = 0.00004_r;
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
