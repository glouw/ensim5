#include "ensim.hh"

#include <array>
#include <numbers>
#include <cmath>
#include <mutex>
#include <cassert>

#define fn __attribute__((used))

namespace ensim
{
    static constexpr double g_dt_s = 1.0 / g_sample_rate_hz;
    static constexpr double g_pi_r = std::numbers::pi_v<double>;
    static constexpr double g_otto_cycle_r = 4.0 * g_pi_r;
    static constexpr double g_otto_intake_cycle_r = 0.0 * g_pi_r;
    static constexpr double g_otto_combustion_cycle_r = 2.0 * g_pi_r;
    static constexpr double g_otto_exhaust_cycle_r = 3.0 * g_pi_r;
    static constexpr double g_resevoir_volume_m3 = 1e9;
    static constexpr double g_ambient_temperature_k = 300.0;
    static constexpr double g_ambient_pressure_pa = 132'800.0;
    static constexpr double g_ambient_density_kg_per_m3 = 1.225;
    static constexpr double g_gamma = 3.0 / 2.0;
    static constexpr double g_universal_gas_constant_j_per_mol_k = 8.3144598;
    static constexpr double g_cv_j_per_mol_k = g_universal_gas_constant_j_per_mol_k / (g_gamma - 1.0);
    static constexpr double g_molar_mass_kg_per_mol = 0.023;
    static constexpr double g_cv_j_per_kg_k = g_cv_j_per_mol_k / g_molar_mass_kg_per_mol;
    static constexpr double g_specific_gas_constant_j_per_kg_k = g_universal_gas_constant_j_per_mol_k / g_molar_mass_kg_per_mol;
    static constexpr double g_energy_octane_j_per_kg = 47.9e6;
    static constexpr double g_stoich_air_fuel_ratio = 14.7;

    static constexpr float g_pi_r_f = g_pi_r;
    static constexpr float g_gamma_f = g_gamma;
    static constexpr float g_dt_s_f = g_dt_s;

    using std::sin;
    using std::cos;
    using std::fmax;
    using std::fmin;
    using std::log;
    using std::sqrt;
    using std::trunc;
    using std::exp;
    using std::fabs;

    template <typename T>
    fn constexpr T clamper(const T value, const T lower, const T upper)
    {
        return fmax(fmin(value, upper), lower);
    }

    template <typename T>
    fn constexpr T modulos(const T value, const T by)
    {
        return value - trunc(value / by) * by;
    }

    template <typename T>
    fn T cuberoot(const T x)
    {
        return exp(log(x) / 3.0);
    }

    fn double frand()
    {
        const double random = 2.0 * rand() / RAND_MAX;
        return random - 1.0;
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
        double piston_injection_enabled = 0.0;
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

                static_assert(g_gamma == 3.0 / 2.0);
                const double Y = cuberoot(Pt / Ps);
                const double M = direction * sqrt(X * (Y - 1.0));
                nozzle_mach[i] = clamper(M, -1.0, 1.0);
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
                const double A = chamber_nozzle_real_flow_area_m2[i];
                const double mute = A == 0.0 ? 0.0 : 1.0;
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

                static_assert(g_gamma == 3.0 / 2.0);
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

                static_assert(g_gamma == 3.0 / 2.0);
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

                static_assert(g_gamma == 3.0 / 2.0);
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
                if(not piston_injection_enabled)
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

                const double randomness = 0.25;
                const double dh = S * g_dt_s;
                const double drh = dh * (1.0 + randomness * frand());
                const double h1 = piston_chamber_flame_height_m;
                const double h2 = h1 + drh;
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
    struct basic_sparkplugs
    {
        static constexpr double fire_delay_theta_r = 1e-1;
        std::array<double, W> engage_theta_r = {};
        std::array<bool, W> prev_fired = {};
        std::array<bool, W> fired = {};
        std::array<bool, W> rising_edge = {};
        double crankshaft_theta_r = 0.0;

        fn void calc_fired()
        {
            for(size_t i = 0; i < W; i++)
            {
                prev_fired[i] = fired[i];
                const double theta0_r = modulos(crankshaft_theta_r, g_otto_cycle_r);
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
    struct basic_cams
    {
        std::array<double, W> engage_theta_r = {};
        std::array<double, W> ramp_theta_r = {};
        std::array<double, W> temp_open_ratio = {};
        std::array<double, W> open_ratio = {};
        double crankshaft_theta_r = 0.0;
        double crankshaft_angular_velocity_r_per_s = 0.0;

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
                double theta1_r = modulos(crankshaft_theta_r, g_otto_cycle_r);
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
                const double R = clamper(A - B + C - D, 0.0, 1.0);
                temp_open_ratio[i] = theta1_r < theta0_r ? 0.0 : R;
            }
        }

        fn virtual void set_open_ratios()
        {
            for(size_t i = 0; i < W; i++)
            {
                open_ratio[i] = temp_open_ratio[i];
            }
        }

        void update()
        {
            calc_open_ratios();
            set_open_ratios();
        }
    };

    template<size_t W>
    struct vtec_cams : basic_cams<W>
    {
        double vtec_engage_r_per_s = 400.0;
        double vtec_boost = 10.0;

        fn void set_open_ratios() override
        {
            for(size_t i = 0; i < W; i++)
            {
                if(this->crankshaft_angular_velocity_r_per_s > vtec_engage_r_per_s)
                {
                    const double boost = clamper(vtec_boost * this->temp_open_ratio[i], 0.0, 1.0);
                    this->open_ratio[i] = boost;
                }
                else
                {
                    this->open_ratio[i] = this->temp_open_ratio[i];
                }
            }
        }
    };

    struct flywheel
    {
        double mass_kg = 0.0;
        double radius_m = 0.0;
        double moment_of_inertia_kg_m2 = 0.0;

        /*
         *      1     2
         * I = --- m r
         *      2
         *
         */

        fn void calc_moment_of_inertia()
        {
            moment_of_inertia_kg_m2 = 0.5 * mass_kg * radius_m * radius_m;
        }

        void update()
        {
            calc_moment_of_inertia();
        }
    };

    struct throttle
    {
        static constexpr size_t size = 4;
        std::array<double, size> table = {};

        double lookup(const double open_ratio)
        {
            const size_t last = size - 1;
            const double at = last * open_ratio;
            const size_t index = at;
            const double ratio = at - index;
            const size_t next = index + 1;
            const double delta = table[next] - table[index];
            return table[index] + delta * ratio;
        }
    };

    struct limiter
    {
        double max_angular_velocity_r_per_s = 600.0;
        double crankshaft_angular_velocity_r_per_s = 0.0;
        double limit_time_s = 0.1;
        double cycles = 0;
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
                const double time_s = g_dt_s * cycles;
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
        double angular_velocity_r_per_s = 0.0;
        double angular_acceleration_r_per_s2 = 0.0;
        double mass_kg = 0.0;
        double radius_m = 0.0;
        double theta_r = 0.0;
        double last_theta_r = 0.0;
        double moment_of_inertia_kg_m2 = 0.0;

        /*
         * dw = a * dt
         *
         */

        fn void accelerate()
        {
            const double a = angular_acceleration_r_per_s2;
            angular_velocity_r_per_s += a * g_dt_s;
        }

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

        /*
         *      1     2
         * I = --- m r
         *      2
         *
         */

        fn void calc_moment_of_inertia()
        {
            moment_of_inertia_kg_m2 = 0.5 * mass_kg * radius_m * radius_m;
        }

        fn bool update()
        {
            calc_moment_of_inertia();
            accelerate();
            turn();
            return otto_cycled();
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
        double crankshaft_angular_velocity_r_per_s = 0.0;
        double crankshaft_theta_r = 0.0;

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
                const double t = theta_r[i];
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
                const double r = crank_throw_length_m[i];
                const double l = connecting_rod_length_m[i];
                const double x = r * sint[i];
                const double y = r * cost[i];
                bearing_x_m[i] = x;
                bearing_y_m[i] = y;
                pin_x_m[i] = 0.0;
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
                const double w = crankshaft_angular_velocity_r_per_s;
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
         *           2
         * Tf = - K w
         *
         */

        fn void calc_friction_torque()
        {
            for(size_t i = 0; i < W; i++)
            {
                const double K = friction_n_m_s2_per_r2[i];
                const double w = crankshaft_angular_velocity_r_per_s;
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
        std::array<float, W> piston_connect_m = {};
        float length_m = 0.0f;
        float mic_position0_m = 0.0f;
        float mic_position1_m = 0.0f;

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

        std::array<float, W> in_velocity_m_per_s = {};
        std::array<float, W> in_static_density_kg_per_m3 = {};
        std::array<float, W> in_static_temperature_k = {};

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
        std::array<float, L> local_speed_of_sound_m_per_s = {};
        std::array<float, L> absolute_speed_of_sound_m_per_s = {};
        std::array<float, L> static_pressure_pa = {};

        /*
         *                1   2    r Rs         1     2
         * rEs = Cv Ts + --- u  = ------- Ts + --- r u
         *                2        y - 1        2
         */

        fn float calc_specific_energy_density_from_static_temperature(const float r, const float u, const float Ts)
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

        fn float calc_specific_energy_density_from_static_pressure(const float r, const float u, const float Ps)
        {
            const float rEs = Ps / (g_gamma_f - 1.0f) + 0.5f * r * u * u;
            return rEs;
        }

        /*
         *                            1     2
         * Ps = (y - 1) * r * [ Es - --- * u ]
         *                            2
         */

        fn float calc_static_pressure_from_specific_energy(const float r, const float u, const float Es)
        {
            const float Ps = (g_gamma_f - 1.0f) * r * (Es - 0.5f * u * u);
            return Ps;
        }

        fn void as_cell(const size_t i, const float r, const float u, const float Ts)
        {
            U_r[i] = r;
            U_ru[i] = r * u;
            U_rEs[i] = calc_specific_energy_density_from_static_temperature(r, u, Ts);
        }

        fn void to_ambient(const size_t i)
        {
            const float r = g_ambient_density_kg_per_m3;
            const float u = 0.0f;
            const float Ts = g_ambient_temperature_k;
            as_cell(i, r, u, Ts);
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
            const float r = U_r[Y];
            const float ru = U_ru[Y];
            const float u = ru / r;
            const float a = local_speed_of_sound_m_per_s[Y];
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

                const float Ps = g_ambient_pressure_pa;
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

        fn float calc_audio_sample()
        {
            const size_t Z = L - 1;
            const size_t x = Z * mic_position0_m / length_m;
            const size_t y = Z * mic_position1_m / length_m;
            return static_pressure_pa[x] + static_pressure_pa[y];
        }

        fn void inject()
        {
            for(size_t i = 0; i < W; i++)
            {
                const float r = in_static_density_kg_per_m3[i];
                const float u = in_velocity_m_per_s[i];
                const float Ts = in_static_temperature_k[i];
                const float ratio = piston_connect_m[i] / length_m;
                as_cell(ratio * L, r, u, Ts);
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
                const float Ps = calc_static_pressure_from_specific_energy(r, u, Es);
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
            const float Ps = calc_static_pressure_from_specific_energy(r, u, Es);
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

        fn void calc_speed_of_sounds()
        {
            for(size_t i = 0; i < L; i++)
            {
                const float Ps = static_pressure_pa[i];
                const float C = sqrt(g_gamma_f * Ps / U_r[i]);
                speed_of_sound_m_per_s[i] = C;
            }
        }

        fn void calc_local_speed_of_sounds()
        {
            for(size_t i = 0; i < L; i++)
            {
                const float a = U_ru[i] / U_r[i];
                local_speed_of_sound_m_per_s[i] = a;
            }
        }

        fn void calc_absolute_speed_of_sounds()
        {
            for(size_t i = 0; i < L; i++)
            {
                const float a = local_speed_of_sound_m_per_s[i];
                const float c = speed_of_sound_m_per_s[i];
                absolute_speed_of_sound_m_per_s[i] = fabs(a) + c;
            }
        }
    };

    #define FLUIDS(X) \
        X(chamber_volume_m3) \
        X(chamber_nozzle_real_flow_area_m2) \
        X(chamber_static_pressure_pa) \
        X(chamber_static_temperature_k) \
        X(chamber_mass_kg) \
        X(nozzle_static_temperature_k) \
        X(nozzle_static_density_kg_per_m3) \
        X(nozzle_velocity_m_per_s)

    #define PISTONS(X) \
        X(total_torque_n_m)

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
        float ratio = 1.0f;

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

    struct convolution_filter
    {
        static constexpr size_t size = 4096;
        static constexpr size_t buffer_size = 2 * size;
        static constexpr size_t mask = size - 1;

        std::array<float, buffer_size> history = {};
        size_t head = 0;

        const std::vector<float>* impulse = nullptr;

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
                y += impulse->operator[](i) * history[i + head];
            }
            head = (head - 1) & mask;
            return y;
        }

        fn void set_impulse(const std::vector<float>& impulse)
        {
            assert(impulse.size() == size);
            this->impulse = &impulse;
        }
    };

    template<size_t W, size_t H>
    struct alignas(std::hardware_destructive_interference_size) mailbox
    {
        /*
         * Recieve
         */

        std::atomic<double> throttle_open_ratio = 0.0;
        std::atomic<size_t> log_x = -1;
        std::atomic<size_t> log_y = -1;
        std::atomic<bool> injection_enabled = true;

        /*
         * Send
         */

        std::atomic<size_t> swap_drops = 0;
        std::atomic<double> engine_angular_velocity_r_per_s = 0.0;
        std::atomic<double> engine_load_torque_n_m = 0.0;
        std::array<std::array<std::atomic<double>, W>, H> port_open_ratios = {};
        std::array<std::array<std::atomic<bool>, W>, H> panics = {};
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
        template<size_t> class PISTONS,
        template<size_t> class CAMS,
        template<size_t> class SPARKPLUGS>
    struct as_engine : engine
    {
        double lumped_drag_torque_n_m = {};
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
        std::vector<float> audio_signal = {};
        struct mailbox<W, H> mailbox = {};
        std::vector<float> pipe_static_pressure_pa = {};
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

        double calc_system_acceleration(const double load_torque_n_m)
        {
            /*
             *      t
             * a = ---
             *      I
             */

            double I = 0.0;
            for(size_t x = 0; x < W; x++)
            {
                I += pistons.moment_of_inertia_kg_m2[x];
            }
            I += flywheel.moment_of_inertia_kg_m2;
            I += crankshaft.moment_of_inertia_kg_m2;
            double t = 0.0;
            for(size_t x = 0; x < W; x++)
            {
                t += pistons.total_torque_n_m[x];
            }
            t -= lumped_drag_torque_n_m;
            t -= load_torque_n_m;
            return t / I;
        }

        fn void broadcast(const double throttle_open_ratio, const bool injection_enabled, const double load_torque_n_m)
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
            inlet_cam.crankshaft_angular_velocity_r_per_s = crankshaft.angular_velocity_r_per_s;
            outlet_cam.crankshaft_angular_velocity_r_per_s = crankshaft.angular_velocity_r_per_s;

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
                const double open_ratio = throttle.lookup(throttle_open_ratio);
                flows[x].chamber_nozzle_open_ratio[THROTTLE_Y] = open_ratio;
            }

            /*
             * Piston shapes <-> flow shapes.
             *
             */

            for(size_t x = 0; x < W; x++)
            {
                flows[x].piston_injection_enabled = injection_enabled;
                flows[x].piston_chamber_radius_m = pistons.diameter_m[x] / 2.0;
                flows[x].chamber_volume_m3[PISTON_Y] = pistons.volumes_m3[x];
                pistons.chamber_static_pressure_pa[x] = flows[x].chamber_static_pressure_pa[PISTON_Y];
            }

            crankshaft.angular_acceleration_r_per_s2 = calc_system_acceleration(load_torque_n_m);
            crankshaft.angular_velocity_r_per_s = fmax(crankshaft.angular_velocity_r_per_s, 0.0);
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
            broadcast(0.0, false, 0.0);
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
                pipe_static_pressure_pa.assign(pipe.static_pressure_pa.begin(), pipe.static_pressure_pa.end());
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

        void update_pipe()
        {
            for(size_t x = 0; x < W; x++)
            {
                const double u = flows[x].nozzle_velocity_m_per_s[AUDIO_Y];
                const double Ts = flows[x].nozzle_static_temperature_k[AUDIO_Y];
                const double r = flows[x].nozzle_static_density_kg_per_m3[AUDIO_Y];
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
            const double throttle_open_ratio = mailbox.throttle_open_ratio;
            const double load_torque_n_m = mailbox.engine_load_torque_n_m;
            const size_t log_x = mailbox.log_x;
            const size_t log_y = mailbox.log_y;
            const bool injection_enabled = mailbox.injection_enabled;
            audio_signal.clear();
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
                float x = pipe.calc_audio_sample();
                x = dc.filter(x);
                x = convolution.filter(x);
                x = gain.filter(x);
                x = clamp.filter(x);
                audio_signal.push_back(x);
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

        const std::atomic<double>& get_angular_velocity_r_per_s() const override
        {
            return mailbox.engine_angular_velocity_r_per_s;
        }

        const std::atomic<double>& get_port_open_ratio(const size_t x, const size_t y) const override
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

        const std::vector<double>& get_signal(const size_t index) const override
        {
            return diags.front[index];
        }

        const std::vector<double>& get_static_temperature_signal_k() const override
        {
            return get_signal(g_chamber_static_temperature_k);
        }

        const std::vector<double>& get_static_pressure_signal_pa() const override
        {
            return get_signal(g_chamber_static_pressure_pa);
        }

        const std::vector<double>& get_volume_signal_m3() const override
        {
            return get_signal(g_chamber_volume_m3);
        }

        const std::vector<float>& get_audio_signal() const override
        {
            return audio_signal;
        }

        const std::vector<float>& get_impulse_signal() const override
        {
            return *convolution.impulse;
        }

        const std::vector<float>& get_pipe_pressure_signal() const override
        {
            return pipe_static_pressure_pa;
        }

        void set_throttle_open_ratio(const double open_ratio) override
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

        void set_load_torque_n_m(const double load_torque_n_m) override
        {
            mailbox.engine_load_torque_n_m = load_torque_n_m;
        }
    };

    const std::vector<float> g_impulse = { -0.001799f, 0.001870f, -0.001897f, 0.001972f, -0.002002f, 0.002081f, -0.002116f, 0.002199f, -0.002238f, 0.002326f, -0.002370f, 0.002463f, -0.002513f, 0.002611f, -0.002668f, 0.002773f, -0.002836f, 0.002949f, -0.003020f, 0.003141f, -0.003223f, 0.003353f, -0.003446f, 0.003588f, -0.003693f, 0.003848f, -0.003969f, 0.004140f, -0.004278f, 0.004468f, -0.004628f, 0.004842f, -0.005028f, 0.005270f, -0.005490f, 0.005769f, -0.006030f, 0.006356f, -0.006672f, 0.007061f, -0.007450f, 0.007923f, -0.008414f, 0.009007f, -0.009644f, 0.010418f, -0.011280f, 0.012338f, -0.013572f, 0.015125f, -0.017043f, 0.019576f, -0.022974f, 0.027908f, -0.035615f, 0.049615f, -0.083066f, 0.278194f, 0.186277f, -0.082104f, 0.245069f, 0.951810f, 0.079744f, -0.015058f, 0.771661f, 0.232401f, -0.053057f, -0.062631f, 0.060914f, 0.138221f, 0.264735f, 0.285731f, 0.495165f, 0.586556f, 0.347162f, 0.139050f, -0.097358f, -0.031730f, -0.123942f, 0.028056f, -0.001174f, 0.423864f, 0.054578f, 0.285723f, 0.172777f, 0.281736f, 0.301603f, 0.141729f, 0.361272f, -0.306274f, 0.093129f, -0.117172f, 0.176220f, 0.125096f, 0.357301f, -0.132178f, 0.130336f, -0.264050f, 0.340199f, -0.132581f, 0.155613f, 0.046111f, 0.083742f, 0.311026f, 0.089839f, 0.331934f, -0.161811f, 0.114417f, -0.028843f, -0.035006f, -0.034659f, -0.079416f, -0.066057f, -0.103561f, 0.041329f, 0.252538f, -0.096950f, 0.398147f, -0.280634f, 0.209116f, -0.014877f, 0.309493f, 0.081499f, 0.231043f, 0.145140f, 0.244074f, 0.234016f, -0.134282f, 0.208181f, 0.444186f, 0.454864f, 0.220934f, -0.043356f, -0.030735f, 0.256108f, 0.055381f, 0.154575f, 0.210247f, 0.250979f, 0.214179f, -0.168535f, 0.134936f, 0.241044f, 0.226504f, -0.003273f, 0.019007f, 0.028781f, 0.025961f, 0.013170f, 0.295959f, 0.098273f, 0.045917f, -0.023161f, -0.007438f, -0.041941f, -0.066265f, 0.005479f, 0.058636f, 0.191558f, 0.190543f, -0.027828f, 0.053463f, 0.113749f, 0.071644f, -0.373187f, -0.028434f, -0.100035f, 0.049823f, 0.056068f, 0.214854f, 0.015167f, 0.058959f, -0.110835f, -0.147336f, -0.148769f, -0.070717f, -0.013135f, 0.074031f, -0.040389f, -0.008819f, 0.105544f, -0.090028f, -0.127972f, -0.068276f, -0.034432f, -0.037090f, 0.025421f, -0.086912f, 0.027017f, -0.041168f, -0.164510f, -0.092424f, 0.028436f, -0.179778f, -0.116844f, 0.033783f, -0.146903f, -0.117936f, 0.087371f, -0.043601f, -0.055109f, 0.022780f, -0.075432f, -0.132115f, -0.108146f, -0.182503f, -0.023194f, 0.038661f, -0.185314f, -0.074016f, 0.027115f, -0.195624f, 0.022350f, -0.093809f, -0.054003f, -0.048363f, -0.142136f, -0.033996f, -0.254705f, -0.130017f, -0.024248f, -0.160407f, 0.004045f, -0.024259f, -0.027768f, 0.056496f, -0.039432f, -0.253253f, -0.035510f, -0.112600f, -0.066573f, -0.075646f, -0.145380f, -0.016447f, -0.212485f, -0.120460f, -0.115753f, 0.117897f, 0.371646f, -0.007276f, 0.184891f, 0.328140f, -0.165394f, -0.072129f, 0.171259f, 0.312006f, -0.053843f, 0.071113f, 0.128257f, -0.045530f, 0.077139f, 0.073022f, 0.128061f, -0.040720f, -0.049558f, -0.114468f, -0.103073f, 0.162806f, -0.048967f, 0.044881f, 0.037149f, -0.049956f, 0.037838f, 0.083957f, 0.111023f, 0.007305f, -0.017689f, -0.045348f, 0.011146f, -0.071549f, 0.015242f, -0.046898f, -0.110946f, -0.037606f, -0.085484f, -0.052864f, -0.125851f, -0.035586f, -0.098230f, 0.214256f, -0.071263f, -0.011016f, 0.019949f, -0.191426f, 0.098025f, 0.095547f, -0.032444f, -0.085893f, -0.180388f, -0.059316f, -0.076177f, -0.073012f, -0.050855f, -0.115151f, -0.069305f, -0.017748f, -0.063978f, -0.092078f, -0.128437f, 0.010558f, -0.116286f, 0.184267f, 0.306916f, 0.035686f, 0.301601f, -0.016923f, -0.005224f, 0.138499f, 0.119391f, 0.005878f, 0.069971f, -0.035157f, 0.051438f, 0.139140f, -0.029178f, -0.075146f, -0.141807f, 0.109912f, 0.060580f, 0.112749f, 0.148563f, 0.015495f, 0.062012f, 0.025774f, -0.008427f, -0.050441f, -0.002154f, -0.052169f, -0.033071f, -0.045143f, -0.061633f, 0.021944f, -0.082090f, 0.010580f, -0.113866f, 0.057723f, -0.040965f, -0.022256f, -0.033617f, -0.065203f, 0.170966f, -0.051701f, -0.059784f, -0.138424f, -0.040055f, -0.033097f, -0.066497f, -0.018655f, -0.084484f, 0.026938f, -0.054901f, -0.049185f, -0.106279f, -0.079742f, -0.096190f, -0.052624f, -0.012356f, -0.121906f, -0.031627f, -0.046556f, -0.104577f, 0.024599f, -0.025226f, -0.052790f, 0.036456f, -0.129811f, 0.006905f, -0.180176f, -0.074949f, -0.056450f, -0.048025f, -0.078997f, -0.069409f, -0.154249f, -0.041052f, -0.098004f, -0.087419f, -0.063104f, -0.053944f, 0.025238f, -0.024094f, 0.028460f, -0.119135f, -0.074987f, -0.043959f, -0.076860f, -0.077455f, -0.096753f, -0.092641f, -0.132838f, -0.115200f, -0.093733f, -0.032898f, -0.069305f, 0.013105f, -0.035825f, 0.019520f, -0.086144f, -0.017441f, -0.107220f, -0.075743f, -0.121200f, -0.079415f, -0.072979f, -0.077916f, -0.118618f, -0.001369f, -0.042012f, -0.008791f, -0.096152f, -0.003497f, -0.053085f, -0.068068f, -0.116905f, -0.034894f, -0.107631f, 0.100961f, 0.110643f, 0.163708f, -0.017788f, 0.173203f, 0.059587f, 0.083201f, 0.008253f, 0.105562f, 0.065965f, -0.015885f, -0.000265f, 0.112866f, 0.057698f, -0.028554f, -0.019836f, 0.045435f, 0.004447f, -0.023734f, 0.057474f, 0.095409f, -0.048049f, 0.006126f, -0.018996f, 0.029955f, -0.013183f, 0.004036f, -0.067787f, -0.027826f, 0.052216f, -0.025644f, -0.042664f, -0.000504f, 0.014166f, -0.015359f, 0.009294f, 0.010726f, -0.057498f, -0.068513f, -0.053842f, -0.000938f, 0.015332f, -0.045396f, -0.036732f, -0.009084f, -0.036910f, -0.081864f, -0.019490f, -0.022143f, -0.042810f, -0.049783f, 0.035754f, 0.008924f, -0.039359f, -0.034464f, -0.047786f, -0.074193f, -0.047534f, -0.091062f, -0.055740f, -0.002583f, -0.001429f, -0.040092f, -0.027285f, 0.056631f, 0.100438f, 0.136638f, 0.022592f, 0.098183f, 0.063429f, 0.054274f, 0.081574f, 0.098657f, 0.000860f, 0.056887f, 0.089501f, 0.004431f, -0.008686f, 0.025457f, 0.021082f, -0.031621f, 0.063403f, 0.031204f, 0.008351f, 0.016055f, 0.025403f, -0.009167f, 0.009595f, -0.037834f, -0.001617f, 0.049250f, -0.074994f, -0.026439f, 0.013793f, 0.004526f, 0.005848f, 0.009373f, -0.069093f, -0.035057f, -0.007151f, -0.005710f, -0.056450f, -0.020186f, -0.023126f, -0.041376f, -0.016612f, 0.015209f, -0.049768f, -0.046633f, -0.026051f, -0.014327f, -0.029718f, -0.047383f, -0.033838f, -0.061792f, -0.055610f, -0.051546f, -0.012590f, 0.007740f, -0.038688f, -0.038396f, -0.054453f, -0.034078f, -0.037464f, -0.064193f, -0.056282f, -0.026198f, -0.058646f, -0.026371f, -0.019380f, -0.035276f, -0.054687f, -0.041510f, -0.069640f, -0.041515f, -0.037090f, -0.042761f, -0.006467f, -0.049704f, -0.059538f, -0.049886f, -0.075691f, -0.032679f, -0.053022f, -0.025205f, -0.025477f, -0.057751f, -0.019783f, -0.017751f, -0.062500f, -0.053738f, -0.051207f, -0.063332f, -0.036576f, -0.033063f, -0.066367f, -0.051128f, -0.047697f, -0.021727f, -0.027127f, -0.023329f, -0.026211f, -0.054969f, -0.033997f, -0.049775f, -0.057335f, -0.035776f, -0.055481f, -0.037573f, -0.040004f, -0.026725f, -0.046362f, -0.032557f, -0.046399f, -0.028837f, -0.032704f, -0.032061f, -0.046842f, -0.050731f, 0.061860f, 0.120198f, 0.064361f, 0.104736f, 0.071804f, 0.045334f, 0.105398f, -0.003035f, 0.044805f, 0.055236f, -0.019191f, 0.057259f, 0.019481f, 0.033449f, 0.055377f, -0.000761f, 0.007682f, 0.020781f, 0.043147f, -0.017136f, 0.021617f, 0.017426f, -0.029675f, 0.006988f, -0.005484f, 0.010020f, 0.019858f, -0.034530f, -0.010112f, 0.026666f, -0.003558f, -0.009766f, 0.001402f, -0.037272f, 0.000251f, -0.040137f, -0.028243f, 0.020460f, -0.008433f, -0.021355f, -0.016569f, -0.029686f, -0.037329f, 0.006188f, -0.003823f, -0.026652f, -0.023249f, -0.015997f, -0.014252f, -0.035486f, -0.026564f, -0.018124f, -0.032205f, -0.011239f, -0.015078f, -0.034137f, -0.050901f, -0.030889f, -0.032785f, -0.009180f, -0.002942f, -0.014080f, -0.017990f, -0.060111f, 0.025081f, 0.116552f, 0.059694f, 0.093554f, 0.052119f, 0.061837f, 0.089927f, -0.014216f, 0.093741f, 0.005454f, 0.025266f, 0.017837f, 0.007350f, 0.049290f, 0.026965f, 0.036592f, 0.028551f, 0.018086f, -0.023031f, 0.033531f, -0.036356f, 0.007319f, 0.011025f, 0.007296f, 0.015628f, -0.035449f, 0.005687f, 0.004208f, -0.008375f, 0.010483f, -0.032113f, -0.002580f, -0.041037f, -0.024509f, 0.005961f, 0.003071f, -0.013033f, -0.019207f, -0.042755f, -0.013558f, 0.003915f, -0.034743f, -0.019735f, -0.028923f, -0.007738f, -0.037812f, -0.022823f, -0.010541f, -0.024606f, -0.009779f, -0.038955f, -0.039550f, -0.042867f, -0.034121f, -0.026397f, 0.001182f, -0.028218f, -0.012426f, -0.037212f, -0.022451f, -0.040321f, -0.031196f, -0.034340f, -0.023552f, -0.022789f, -0.036474f, -0.028925f, -0.024455f, -0.035206f, -0.037884f, -0.039152f, -0.021257f, -0.016186f, -0.025265f, -0.030688f, -0.037463f, -0.031422f, -0.047306f, -0.023126f, -0.022111f, -0.032552f, -0.021836f, -0.028060f, -0.032736f, -0.022162f, -0.029822f, -0.034184f, -0.036620f, -0.029992f, -0.030885f, -0.027843f, -0.015338f, -0.038512f, -0.025558f, -0.031855f, -0.021915f, -0.045252f, -0.018953f, -0.032100f, -0.012693f, -0.029502f, -0.015120f, -0.030574f, -0.019145f, -0.043226f, -0.027101f, -0.037971f, -0.013225f, -0.035435f, -0.014216f, -0.042633f, -0.008658f, -0.029596f, -0.002976f, -0.042951f, 0.057427f, 0.079876f, 0.076685f, 0.052328f, 0.053981f, 0.057366f, 0.032857f, 0.063170f, 0.017810f, 0.042581f, 0.014990f, 0.041838f, -0.000735f, 0.037180f, 0.019571f, 0.003604f, 0.026482f, 0.003039f, 0.001743f, 0.023425f, 0.026375f, -0.010365f, 0.007745f, 0.013099f, -0.007578f, 0.009110f, -0.009529f, -0.007603f, -0.015759f, 0.011038f, -0.005631f, -0.001573f, -0.021553f, -0.017179f, 0.022066f, 0.000296f, -0.013912f, -0.006305f, -0.016266f, -0.018930f, -0.010254f, -0.012620f, -0.002841f, -0.012317f, -0.026184f, -0.032865f, -0.013463f, 0.006624f, 0.001800f, -0.018780f, -0.025364f, -0.020086f, -0.022307f, -0.016571f, -0.005244f, -0.010993f, -0.021025f, -0.017189f, -0.021884f, -0.025180f, -0.028159f, -0.016538f, 0.010156f, -0.024290f, 0.008708f, 0.084360f, 0.057376f, 0.053053f, 0.049688f, 0.056057f, 0.043604f, 0.050254f, 0.017638f, 0.038759f, 0.018276f, 0.022022f, -0.001910f, 0.041661f, 0.002010f, 0.042157f, -0.002806f, 0.013701f, 0.007712f, 0.008907f, -0.019653f, 0.022474f, -0.001176f, 0.008345f, -0.007988f, 0.004772f, -0.029964f, 0.015031f, -0.002386f, -0.006245f, -0.027411f, -0.010640f, 0.012650f, -0.003098f, -0.009762f, -0.011774f, -0.026702f, -0.008491f, -0.017501f, -0.003462f, -0.013652f, -0.024973f, -0.023965f, -0.018094f, 0.003752f, -0.010606f, -0.022513f, -0.022779f, -0.016584f, -0.023664f, -0.016755f, -0.009450f, -0.026862f, -0.017865f, -0.017989f, -0.024995f, -0.029521f, -0.014458f, -0.002857f, -0.020695f, -0.020298f, -0.023905f, -0.031549f, -0.015249f, -0.019390f, -0.015716f, -0.020870f, -0.026750f, -0.019046f, -0.021736f, -0.025748f, -0.028169f, -0.025607f, -0.014724f, -0.009103f, -0.015496f, -0.020012f, -0.023883f, -0.025953f, -0.019234f, -0.017849f, -0.021083f, -0.024974f, -0.018566f, -0.025042f, -0.022077f, -0.018308f, -0.019515f, -0.016733f, -0.024097f, -0.018851f, -0.014590f, -0.019152f, -0.020147f, -0.022684f, -0.018279f, -0.021268f, -0.021523f, -0.012603f, -0.017320f, -0.013227f, -0.017180f, -0.021428f, -0.018934f, -0.019099f, -0.016081f, -0.023493f, -0.017097f, -0.021455f, -0.013994f, -0.014920f, -0.011336f, -0.019194f, -0.014545f, -0.020752f, 0.034230f, 0.078463f, 0.050348f, 0.054549f, 0.046567f, 0.037659f, 0.039584f, 0.027459f, 0.030024f, 0.031522f, 0.014714f, 0.030080f, 0.009142f, 0.020026f, 0.009339f, 0.009799f, 0.023990f, -0.002206f, 0.013216f, 0.017994f, 0.006275f, 0.007920f, -0.004017f, -0.002436f, 0.015985f, -0.001307f, -0.008218f, -0.009445f, 0.017372f, -0.008805f, -0.002215f, -0.011389f, -0.005147f, -0.001910f, 0.001631f, -0.002017f, -0.014827f, -0.014248f, -0.005191f, 0.008987f, -0.010256f, -0.013724f, -0.012406f, -0.012166f, -0.008852f, -0.001084f, -0.011339f, -0.009061f, -0.012732f, -0.019521f, -0.011492f, 0.003432f, -0.013204f, -0.014783f, -0.019451f, -0.015777f, -0.011144f, -0.010634f, -0.004672f, -0.012807f, -0.007579f, -0.016483f, -0.013143f, -0.019296f, 0.004753f, 0.073424f, 0.055120f, 0.044931f, 0.036889f, 0.027617f, 0.042590f, 0.026763f, 0.034019f, 0.023542f, 0.019377f, 0.017718f, 0.010068f, 0.018949f, 0.002422f, 0.016458f, 0.008339f, 0.005309f, 0.012764f, 0.002244f, 0.008453f, -0.006289f, -0.000013f, 0.013726f, -0.001165f, -0.014535f, -0.002047f, 0.003604f, -0.009454f, -0.003652f, -0.012842f, -0.001498f, -0.004566f, 0.001242f, -0.012866f, -0.016068f, -0.009410f, 0.004625f, -0.011232f, -0.017110f, -0.014307f, -0.016697f, -0.002647f, -0.006455f, -0.013913f, -0.011157f, -0.019950f, -0.016043f, -0.001945f, -0.011263f, -0.018242f, -0.016175f, -0.017664f, -0.014876f, -0.008102f, -0.007770f, -0.017986f, -0.013404f, -0.019445f, -0.017091f, -0.016358f, -0.009445f, -0.012104f, -0.013782f, -0.018192f, -0.019230f, -0.014186f, -0.011177f, -0.013918f, -0.014981f, -0.016354f, -0.016739f, -0.013637f, -0.013832f, -0.015234f, -0.017776f, -0.011303f, -0.015021f, -0.018138f, -0.017998f, -0.013347f, -0.015970f, -0.010933f, -0.010302f, -0.012873f, -0.013103f, -0.016470f, -0.014716f, -0.014307f, -0.014394f, -0.016347f, -0.017075f, -0.011377f, -0.010431f, -0.013476f, -0.014561f, -0.012599f, -0.013648f, -0.012756f, -0.016199f, -0.010307f, -0.010315f, -0.013373f, -0.010395f, -0.010882f, -0.009786f, -0.012323f, -0.015300f, -0.015013f, -0.010735f, -0.011968f, -0.011615f, -0.010325f, -0.008162f, -0.014319f, 0.017017f, 0.065784f, 0.041376f, 0.048228f, 0.030220f, 0.040276f, 0.022114f, 0.030701f, 0.017116f, 0.023197f, 0.017485f, 0.021176f, 0.009618f, 0.018200f, 0.013847f, 0.003634f, 0.016586f, 0.004671f, 0.008221f, -0.004260f, 0.010576f, 0.009668f, -0.006010f, 0.001857f, 0.013347f, -0.001849f, 0.000494f, -0.005781f, 0.000816f, 0.004007f, -0.001915f, -0.013338f, -0.004401f, 0.009149f, -0.002643f, -0.010126f, -0.005430f, -0.007564f, -0.000016f, -0.006578f, -0.006780f, -0.010536f, -0.013266f, 0.004875f, -0.000809f, -0.010765f, -0.009631f, -0.009678f, -0.007007f, -0.003709f, -0.010065f, -0.009513f, -0.010320f, -0.011212f, -0.007890f, -0.004272f, -0.005412f, -0.008959f, -0.012472f, -0.011414f, -0.006543f, -0.008916f, -0.007987f, -0.010998f, -0.006026f, 0.054260f, 0.046802f, 0.037857f, 0.029998f, 0.036691f, 0.021340f, 0.024605f, 0.017917f, 0.018261f, 0.016858f, 0.017022f, 0.009007f, 0.016340f, 0.002572f, 0.008647f, 0.009084f, 0.005951f, -0.001674f, 0.000923f, 0.015421f, -0.006698f, -0.002433f, 0.007939f, -0.000269f, -0.004383f, -0.005267f, 0.001412f, 0.001997f, -0.008668f, -0.012234f, -0.001028f, 0.003698f, -0.013450f, -0.007973f, -0.007283f, -0.003091f, -0.002835f, -0.006924f, -0.009225f, -0.017439f, -0.003168f, -0.003411f, -0.010896f, -0.012510f, -0.010783f, -0.006978f, -0.005847f, -0.010844f, -0.009835f, -0.012450f, -0.013849f, -0.007684f, -0.005367f, -0.008705f, -0.012526f, -0.015063f, -0.010172f, -0.010382f, -0.010811f, -0.009775f, -0.012086f, -0.009720f, -0.009851f, -0.013031f, -0.010197f, -0.007936f, -0.012880f, -0.012133f, -0.013362f, -0.010536f, -0.008057f, -0.010027f, -0.010295f, -0.012463f, -0.010596f, -0.011144f, -0.011112f, -0.012293f, -0.011094f, -0.006666f, -0.010720f, -0.010782f, -0.010730f, -0.010610f, -0.012038f, -0.010415f, -0.008590f, -0.011063f, -0.009319f, -0.007877f, -0.007923f, -0.011223f, -0.011265f, -0.009506f, -0.009681f, -0.010097f, -0.008086f, -0.009905f, -0.008306f, -0.009841f, -0.007509f, -0.008873f, -0.008190f, -0.007923f, -0.008901f, -0.009630f, -0.008430f, -0.008864f, -0.008166f, -0.008403f, -0.005230f, -0.009276f, -0.003911f, -0.012277f, 0.009796f, 0.050451f, 0.036713f, 0.034378f, 0.029207f, 0.028120f, 0.022787f, 0.022176f, 0.015443f, 0.018596f, 0.012183f, 0.012194f, 0.013492f, 0.007577f, 0.013172f, 0.008526f, 0.006597f, 0.004158f, 0.012990f, -0.003147f, 0.000585f, 0.009255f, 0.001916f, 0.001521f, -0.001127f, 0.006570f, -0.000050f, -0.007654f, -0.000707f, 0.007186f, -0.006736f, -0.002782f, -0.003798f, 0.002499f, -0.002593f, -0.002999f, -0.009352f, -0.003398f, 0.002879f, -0.006735f, -0.006882f, -0.006093f, -0.001760f, -0.002602f, -0.007072f, -0.005026f, -0.008493f, -0.006706f, -0.003381f, -0.003230f, -0.007578f, -0.008191f, -0.004955f, -0.004835f, -0.005946f, -0.006163f, -0.008866f, -0.005943f, -0.007634f, -0.007488f, -0.003774f, -0.007237f, -0.007081f, -0.005841f, -0.007223f, 0.040488f, 0.039124f, 0.029483f, 0.026358f, 0.024085f, 0.019645f, 0.018513f, 0.013148f, 0.018394f, 0.009401f, 0.012010f, 0.010048f, 0.005779f, 0.008926f, 0.007741f, 0.001892f, 0.005070f, 0.007195f, -0.005124f, 0.008069f, -0.001335f, -0.000076f, -0.002734f, 0.003188f, 0.002931f, -0.008619f, -0.003798f, 0.003676f, -0.005939f, -0.004934f, -0.004422f, -0.000058f, -0.004851f, -0.004865f, -0.010964f, -0.001236f, -0.002396f, -0.007535f, -0.007800f, -0.005825f, -0.003231f, -0.007271f, -0.007884f, -0.009312f, -0.008400f, -0.003467f, -0.004639f, -0.008362f, -0.010981f, -0.007326f, -0.006237f, -0.007439f, -0.007991f, -0.008333f, -0.006377f, -0.009358f, -0.007373f, -0.006173f, -0.010119f, -0.009001f, -0.009354f, -0.006460f, -0.006298f, -0.007855f, -0.009216f, -0.008479f, -0.008102f, -0.008551f, -0.009816f, -0.007111f, -0.006515f, -0.008694f, -0.007303f, -0.008333f, -0.008783f, -0.006827f, -0.007781f, -0.009316f, -0.006532f, -0.006520f, -0.008144f, -0.009394f, -0.006737f, -0.008517f, -0.006742f, -0.007041f, -0.007619f, -0.007919f, -0.006280f, -0.006913f, -0.006600f, -0.006853f, -0.006996f, -0.007689f, -0.007165f, -0.007767f, -0.006920f, -0.006016f, -0.004830f, -0.005495f, -0.006625f, -0.006939f, -0.007586f, -0.006666f, -0.004966f, -0.005180f, -0.006067f, -0.005737f, -0.005864f, -0.006504f, -0.005456f, -0.007217f, -0.002520f, -0.007288f, 0.003725f, 0.040927f, 0.031065f, 0.027141f, 0.022947f, 0.020147f, 0.020147f, 0.015549f, 0.016650f, 0.011235f, 0.012225f, 0.009898f, 0.007151f, 0.010980f, 0.006037f, 0.004935f, 0.010502f, -0.000612f, 0.004783f, 0.007158f, -0.000707f, 0.000949f, 0.003961f, 0.005406f, -0.005891f, 0.004734f, 0.002445f, -0.003085f, -0.003155f, 0.000998f, -0.000019f, -0.001583f, -0.006035f, 0.000379f, 0.001638f, -0.005100f, -0.004582f, -0.002616f, -0.000495f, -0.004184f, -0.003120f, -0.005997f, -0.002051f, -0.001818f, -0.004226f, -0.007257f, -0.003725f, -0.002869f, -0.003357f, -0.005169f, -0.005152f, -0.003425f, -0.006745f, -0.002595f, -0.005424f, -0.005504f, -0.005033f, -0.004866f, -0.003902f, -0.004601f, -0.005904f, -0.004162f, -0.005373f, -0.004657f, -0.007449f, 0.028417f, 0.033288f, 0.022736f, 0.023532f, 0.016277f, 0.018721f, 0.013099f, 0.013260f, 0.009343f, 0.011772f, 0.006489f, 0.006018f, 0.008974f, 0.002673f, 0.006124f, 0.006923f, -0.003391f, 0.006594f, 0.001388f, 0.000001f, 0.000481f, 0.004420f, -0.004880f, -0.002379f, 0.003100f, -0.003027f, -0.003094f, -0.000448f, -0.000647f, -0.003023f, -0.007557f, -0.002288f, 0.000346f, -0.005705f, -0.005391f, -0.003131f, -0.002854f, -0.005071f, -0.005282f, -0.006893f, -0.001444f, -0.003976f, -0.006376f, -0.007633f, -0.004321f, -0.005935f, -0.005061f, -0.006570f, -0.003643f, -0.006793f, -0.004634f, -0.005490f, -0.007215f, -0.006577f, -0.006524f, -0.004309f, -0.005742f, -0.006855f, -0.005850f, -0.006113f, -0.006652f, -0.006928f, -0.004761f, -0.006129f, -0.005379f, -0.006712f, -0.007078f, -0.005218f, -0.006421f, -0.006467f, -0.004762f, -0.006110f, -0.006595f, -0.005641f, -0.005929f, -0.005297f, -0.005983f, -0.007019f, -0.005788f, -0.005554f, -0.004874f, -0.004879f, -0.005675f, -0.005868f, -0.006076f, -0.006399f, -0.004930f, -0.004909f, -0.004806f, -0.004947f, -0.005659f, -0.005784f, -0.005453f, -0.004148f, -0.004172f, -0.005076f, -0.005348f, -0.005114f, -0.004689f, -0.005483f, -0.003546f, -0.004480f, -0.004238f, -0.004514f, -0.004887f, -0.004921f, -0.004726f, -0.003266f, -0.003356f, -0.004217f, -0.004199f, -0.004610f, -0.003490f, -0.004903f, -0.000403f, 0.032204f, 0.026660f, 0.022187f, 0.018626f, 0.016688f, 0.013813f, 0.013126f, 0.011368f, 0.010822f, 0.009755f, 0.006027f, 0.010525f, 0.002082f, 0.008646f, 0.003853f, 0.002768f, 0.007489f, 0.000756f, 0.001307f, 0.004882f, 0.000700f, -0.002358f, 0.006477f, -0.001211f, -0.000716f, 0.001111f, 0.001469f, -0.002098f, -0.003901f, 0.002365f, -0.000498f, -0.004335f, -0.000819f, 0.000022f, -0.002133f, -0.003060f, -0.003505f, -0.001077f, -0.000762f, -0.004952f, -0.003213f, -0.001807f, -0.001684f, -0.003872f, -0.003425f, -0.003368f, -0.003784f, -0.001912f, -0.004264f, -0.003911f, -0.003624f, -0.002514f, -0.003254f, -0.004372f, -0.003057f, -0.003751f, -0.004167f, -0.003859f, -0.002932f, -0.004915f, -0.002579f, -0.005684f, -0.001491f, -0.006182f, 0.018081f, 0.028921f, 0.018325f, 0.017379f, 0.014739f, 0.013243f, 0.011431f, 0.009915f, 0.007229f, 0.008583f, 0.004821f, 0.008832f, 0.000929f, 0.007754f, 0.000656f, 0.002578f, 0.004213f, 0.000323f, 0.001444f, 0.003821f, -0.003908f, 0.000380f, 0.002174f, -0.002776f, -0.000280f, -0.000111f, -0.001174f, -0.005002f, -0.000999f, 0.000273f, -0.004407f, -0.002879f, -0.001218f, -0.003213f, -0.004150f, -0.004443f, -0.001572f, -0.001581f, -0.005719f, -0.004181f, -0.003467f, -0.003274f, -0.005333f, -0.003062f, -0.004171f, -0.003758f, -0.003741f, -0.005536f, -0.004792f, -0.004753f, -0.003174f, -0.004737f, -0.004546f, -0.004472f, -0.004699f, -0.005089f, -0.003643f, -0.005088f, -0.003845f, -0.005479f, -0.004496f, -0.004335f, -0.005471f, -0.004041f, -0.004567f, -0.005238f, -0.003860f, -0.004371f, -0.004087f, -0.005286f, -0.005295f, -0.004229f, -0.004326f, -0.004028f, -0.004071f, -0.004513f, -0.004777f, -0.004780f, -0.003726f, -0.003811f, -0.004317f, -0.004568f, -0.004531f, -0.004437f, -0.003460f, -0.003697f, -0.003969f, -0.004361f, -0.003984f, -0.004493f, -0.003319f, -0.003159f, -0.003947f, -0.003656f, -0.004039f, -0.003872f, -0.003847f, -0.002791f, -0.003407f, -0.003306f, -0.003884f, -0.003370f, -0.003642f, -0.002994f, -0.002900f, -0.002925f, -0.003560f, -0.002960f, -0.003092f, -0.003210f, -0.003030f, -0.002160f, -0.003023f, -0.002796f, -0.002649f, 0.023061f, 0.022697f, 0.017622f, 0.016113f, 0.013209f, 0.011675f, 0.010313f, 0.008383f, 0.008398f, 0.005950f, 0.008168f, 0.003491f, 0.007708f, 0.002000f, 0.004401f, 0.004261f, 0.000529f, 0.004285f, 0.002158f, -0.001039f, 0.004871f, -0.000888f, 0.000076f, 0.001663f, 0.000968f, -0.002817f, 0.001177f, 0.001187f, -0.003076f, -0.000668f, 0.000522f, -0.001382f, -0.002011f, -0.002227f, 0.000532f, -0.002330f, -0.003009f, -0.001454f, -0.001008f, -0.002606f, -0.002037f, -0.002714f, -0.001827f, -0.001735f, -0.003106f, -0.002756f, -0.002194f, -0.002218f, -0.003298f, -0.002316f, -0.002597f, -0.003073f, -0.001918f, -0.003150f, -0.002579f, -0.003345f, -0.002688f, -0.002657f, -0.003256f, -0.002479f, -0.002687f, -0.003994f, -0.001228f, -0.005033f, 0.012159f, 0.022861f, 0.014575f, 0.015039f, 0.011376f, 0.011026f, 0.008376f, 0.007759f, 0.006136f, 0.005609f, 0.006374f, 0.002850f, 0.005075f, 0.000988f, 0.003360f, 0.003683f, -0.000513f, 0.004108f, -0.001508f, 0.000046f, 0.001996f, -0.001576f, -0.000247f, 0.001182f, -0.001851f, -0.002830f, 0.001162f, -0.001724f, -0.002727f, -0.000471f, -0.002376f, -0.002550f, -0.003441f, -0.000415f, -0.002037f, -0.003835f, -0.002448f, -0.002001f, -0.003650f, -0.003077f, -0.002859f, -0.002435f, -0.002707f, -0.003844f, -0.003583f, -0.002947f, -0.003067f, -0.003656f, -0.002948f, -0.003738f, -0.003540f, -0.003163f, -0.003596f, -0.003299f, -0.003935f, -0.002736f, -0.003927f, -0.003627f, -0.003433f, -0.004095f, -0.003339f, -0.003221f, -0.002769f, -0.004118f, -0.003859f, -0.003450f, -0.003101f, -0.003569f, -0.003533f, -0.003520f, -0.003754f, -0.003111f, -0.002947f, -0.003276f, -0.003627f, -0.003831f, -0.003128f, -0.002643f, -0.003240f, -0.003547f, -0.003187f, -0.003554f, -0.002938f, -0.002660f, -0.003081f, -0.002908f, -0.003396f, -0.003230f, -0.002634f, -0.002443f, -0.003018f, -0.003120f, -0.002985f, -0.002699f, -0.002760f, -0.002151f, -0.002745f, -0.002647f, -0.002857f, -0.002704f, -0.002806f, -0.001863f, -0.002409f, -0.002339f, -0.002543f, -0.002458f, -0.002557f, -0.001996f, -0.002202f, -0.001912f, -0.002449f, -0.001816f, -0.002757f, -0.001302f, -0.003173f, 0.016794f, 0.019518f, 0.013525f, 0.013006f, 0.010141f, 0.009960f, 0.007419f, 0.008087f, 0.005207f, 0.006588f, 0.003420f, 0.005978f, 0.002422f, 0.004262f, 0.003143f, 0.001077f, 0.004589f, -0.001069f, 0.003159f, 0.001239f, -0.000456f, 0.002035f, 0.001061f, -0.001705f, 0.001619f, 0.000238f, -0.001965f, 0.000713f, -0.000510f, -0.000944f, -0.001481f, 0.000427f, -0.001224f, -0.002170f, -0.000583f, -0.001126f, -0.001905f, -0.001374f, -0.001674f, -0.000845f, -0.002003f, -0.002096f, -0.001562f, -0.001649f, -0.002370f, -0.001590f, -0.002124f, -0.002031f, -0.001455f, -0.002370f, -0.002006f, -0.002592f, -0.001484f, -0.002685f, -0.001434f, -0.002702f, -0.002022f, -0.002293f, -0.001554f, -0.002651f, -0.002253f, -0.002791f, -0.000923f, -0.003868f, 0.006713f, 0.019352f, 0.011595f, 0.012628f, 0.008797f, 0.008913f, 0.005968f, 0.006631f, 0.004528f, 0.005995f, 0.002385f, 0.004533f, 0.001100f, 0.003316f, 0.002147f, 0.001000f, 0.002980f, -0.001764f, 0.002325f, -0.000140f, -0.000318f, 0.001107f, -0.001143f, -0.001727f, 0.001353f, -0.001996f, -0.000935f, -0.000468f, -0.001548f, -0.002378f, -0.001368f, -0.000600f, -0.002613f, -0.001667f, -0.001425f, -0.002590f, -0.002233f, -0.002206f, -0.001505f, -0.002449f, -0.002782f, -0.002242f, -0.002174f, -0.002782f, -0.002239f, -0.002545f, -0.002651f, -0.002434f, -0.002866f, -0.002534f, -0.002911f, -0.001989f, -0.003123f, -0.002204f, -0.003458f, -0.002347f, -0.002818f, -0.002225f, -0.003135f, -0.002731f, -0.002755f, -0.002274f, -0.002889f, -0.002647f, -0.003163f, -0.002369f, -0.002369f, -0.002537f, -0.002939f, -0.002803f, -0.002592f, -0.002078f, -0.002801f, -0.002475f, -0.002908f, -0.002472f, -0.002169f, -0.002330f, -0.002625f, -0.002569f, -0.002645f, -0.001964f, -0.002191f, -0.002284f, -0.002659f, -0.002320f, -0.002564f, -0.001757f, -0.002198f, -0.002043f, -0.002340f, -0.002117f, -0.002407f, -0.001583f, -0.002233f, -0.001819f, -0.002316f, -0.001832f, -0.002117f, -0.001507f, -0.002131f, -0.001656f, -0.002188f, -0.001584f, -0.001924f, -0.001362f, -0.001948f, -0.001323f, -0.002267f, -0.001300f, -0.002099f, -0.000833f, -0.002095f, -0.000439f, -0.003463f, 0.011685f, 0.016215f, 0.011046f, 0.010294f, 0.008520f, 0.007464f, 0.006594f, 0.005382f, 0.005388f, 0.003794f, 0.004370f, 0.002867f, 0.003264f, 0.002684f, 0.001596f, 0.003277f, -0.000448f, 0.003407f, -0.000240f, 0.001162f, 0.001579f, -0.000671f, 0.000821f, 0.000917f, -0.001573f, 0.001094f, -0.000402f, -0.000719f, -0.000709f, 0.000607f, -0.001780f, -0.000522f, -0.000419f, -0.001101f, -0.001188f, -0.001031f, -0.000789f, -0.001420f, -0.001574f, -0.000812f, -0.001312f, -0.001460f, -0.001195f, -0.001780f, -0.001134f, -0.001654f, -0.001559f, -0.001812f, -0.001105f, -0.001823f, -0.001129f, -0.002215f, -0.001454f, -0.001682f, -0.001324f, -0.002192f, -0.001403f, -0.001762f, -0.001218f, -0.002146f, -0.001373f, -0.002569f, -0.000487f, -0.002969f, 0.003666f, 0.015055f, 0.010348f, 0.009824f, 0.007287f, 0.006397f, 0.005365f, 0.004520f, 0.004832f, 0.003239f, 0.003320f, 0.001863f, 0.002272f, 0.002186f, 0.001065f, 0.002061f, -0.000759f, 0.002412f, -0.001026f, 0.001183f, 0.000250f, -0.001351f, 0.000593f, -0.000740f, -0.001081f, 0.000455f, -0.001420f, -0.001391f, -0.000824f, -0.000686f, -0.002183f, -0.000407f, -0.001722f, -0.001479f, -0.001904f, -0.001036f, -0.001736f, -0.001840f, -0.001761f, -0.001653f, -0.002140f, -0.001596f, -0.002129f, -0.001449f, -0.002190f, -0.001907f, -0.002238f, -0.001725f, -0.002310f, -0.001590f, -0.002560f, -0.001910f, -0.001966f, -0.001853f, -0.002555f, -0.001896f, -0.001959f, -0.001991f, -0.002215f, -0.002183f, -0.002279f, -0.001641f, -0.002114f, -0.002101f, -0.002412f, -0.001639f, -0.001953f, -0.002132f, -0.002200f, -0.002042f, -0.001848f, -0.001712f, -0.002056f, -0.002014f, -0.002285f, -0.001576f, -0.001766f, -0.001887f, -0.002050f, -0.001845f, -0.001874f, -0.001536f, -0.001955f, -0.001657f, -0.001908f, -0.001691f, -0.001566f, -0.001574f, -0.001897f, -0.001675f, -0.001801f, -0.001321f, -0.001510f, -0.001506f, -0.001678f, -0.001539f, -0.001642f, -0.001384f, -0.001382f, -0.001212f, -0.001694f, -0.001285f, -0.001499f, -0.001231f, -0.001375f, -0.001181f, -0.001445f, -0.001016f, -0.001488f, -0.000979f, -0.001465f, -0.000794f, -0.001748f, -0.000309f, -0.002533f, 0.007326f, 0.013953f, 0.008692f, 0.008686f, 0.006318f, 0.006444f, 0.004668f, 0.004876f, 0.003613f, 0.003738f, 0.002737f, 0.002731f, 0.002346f, 0.001539f, 0.002509f, 0.000461f, 0.002599f, -0.000272f, 0.001690f, 0.000503f, -0.000257f, 0.001514f, -0.000664f, 0.000362f, 0.000502f, -0.000708f, -0.000252f, 0.000286f, -0.001135f, -0.000144f, -0.000294f, -0.000876f, -0.000735f, -0.000613f, -0.000647f, -0.001311f, -0.000491f, -0.000947f, -0.001001f, -0.001032f, -0.001119f, -0.000963f, -0.001112f, -0.001248f, -0.000964f, -0.001213f, -0.001090f, -0.001296f, -0.001450f, -0.001033f, -0.001076f, -0.001439f, -0.001370f, -0.001081f, -0.001221f, -0.001391f, -0.001467f, -0.001375f, -0.000847f, -0.001385f, -0.001364f, -0.001688f, -0.000641f, -0.001965f, 0.001192f, 0.012163f, 0.008536f, 0.007967f, 0.005765f, 0.005381f, 0.003876f, 0.004083f, 0.003287f, 0.002826f, 0.002153f, 0.001893f, 0.001877f, 0.001211f, 0.001661f, -0.000111f, 0.001724f, -0.000766f, 0.001531f, -0.000370f, -0.000334f, 0.000524f, -0.001187f, 0.000264f, -0.000376f, -0.001163f, -0.000375f, -0.000610f, -0.001467f, -0.000302f, -0.001142f, -0.001223f, -0.001228f, -0.000779f, -0.001533f, -0.001129f, -0.001189f, -0.001520f, -0.001294f, -0.001535f, -0.001367f, -0.001257f, -0.001700f, -0.001437f, -0.001490f, -0.001532f, -0.001555f, -0.001701f, -0.001453f, -0.001459f, -0.001736f, -0.001737f, -0.001441f, -0.001530f, -0.001622f, -0.001785f, -0.001566f, -0.001390f, -0.001729f, -0.001781f, -0.001633f, -0.001258f, -0.001648f, -0.001621f, -0.001811f, -0.001461f, -0.001392f, -0.001548f, -0.001668f, -0.001677f, -0.001294f, -0.001455f, -0.001612f, -0.001519f, -0.001508f, -0.001250f, -0.001463f, -0.001511f, -0.001543f, -0.001485f, -0.001155f, -0.001246f, -0.001463f, -0.001414f, -0.001423f, -0.001199f, -0.001312f, -0.001191f, -0.001347f, -0.001259f, -0.001229f, -0.001147f, -0.001194f, -0.001255f, -0.001174f, -0.001079f, -0.001170f, -0.000986f, -0.001253f, -0.001028f, -0.001131f, -0.000957f, -0.001045f, -0.000893f, -0.001143f, -0.000908f, -0.001101f, -0.000832f, -0.001024f, -0.000641f, -0.001189f, -0.000548f, -0.001285f, -0.000282f, -0.001997f, 0.004463f, 0.011485f, 0.007179f, 0.007053f, 0.005036f, 0.005096f, 0.003767f, 0.003873f, 0.002811f, 0.003042f, 0.001944f, 0.002508f, 0.001192f, 0.002299f, 0.000437f, 0.002178f, -0.000175f, 0.001698f, 0.000013f, 0.000684f, 0.000737f, -0.000379f, 0.000700f, -0.000115f, -0.000424f, 0.000524f, -0.000626f, -0.000173f, -0.000122f, -0.000605f, -0.000523f, -0.000338f, -0.000566f, -0.000718f, -0.000358f, -0.000845f, -0.000656f, -0.000855f, -0.000714f, -0.000692f, -0.000951f, -0.000760f, -0.000860f, -0.000827f, -0.000997f, -0.001040f, -0.000697f, -0.000891f, -0.001170f, -0.000919f, -0.000870f, -0.000989f, -0.001104f, -0.001059f, -0.000724f, -0.000977f, -0.001186f, -0.001130f, -0.000741f, -0.001066f, -0.001028f, -0.001187f, -0.000751f, -0.001106f, -0.000046f, 0.009254f, 0.007263f, 0.006271f, 0.004777f, 0.004190f, 0.003274f, 0.003188f, 0.002519f, 0.002272f, 0.001541f, 0.001750f, 0.000977f, 0.001724f, 0.000121f, 0.001374f, -0.000393f, 0.001200f, -0.000328f, 0.000251f, 0.000109f, -0.000700f, 0.000421f, -0.000703f, -0.000396f, -0.000217f, -0.001024f, -0.000403f, -0.000592f, -0.000898f, -0.000843f, -0.000608f, -0.001133f, -0.000857f, -0.000830f, -0.001126f, -0.000984f, -0.001149f, -0.000936f, -0.001216f, -0.001080f, -0.001146f, -0.001118f, -0.001159f, -0.001394f, -0.001040f, -0.001096f, -0.001394f, -0.001226f, -0.001110f, -0.001287f, -0.001287f, -0.001352f, -0.001042f, -0.001201f, -0.001403f, -0.001384f, -0.000993f, -0.001271f, -0.001257f, -0.001375f, -0.001102f, -0.001202f, -0.001239f, -0.001363f, -0.001125f, -0.001018f, -0.001283f, -0.001289f, -0.001259f, -0.001045f, -0.001034f, -0.001203f, -0.001192f, -0.001232f, -0.000954f, -0.001167f, -0.001056f, -0.001185f, -0.001013f, -0.000997f, -0.001053f, -0.001115f, -0.001032f, -0.000989f, -0.000935f, -0.000969f, -0.001047f, -0.001024f, -0.000911f, -0.000916f, -0.000789f, -0.001048f, -0.000818f, -0.000990f, -0.000819f, -0.000898f, -0.000756f, -0.000898f, -0.000723f, -0.000922f, -0.000680f, -0.000872f, -0.000700f, -0.000835f, -0.000587f, -0.000906f, -0.000535f, -0.000869f, -0.000487f, -0.000900f, -0.000395f, -0.001028f, -0.000083f, -0.001491f, 0.002504f, 0.009148f, 0.006070f, 0.005473f, 0.004333f, 0.003910f, 0.003237f, 0.002808f, 0.002464f, 0.001996f, 0.002046f, 0.001417f, 0.001683f, 0.000838f, 0.001396f, 0.000436f, 0.001068f, 0.000407f, 0.000491f, 0.000642f, -0.000274f, 0.000853f, -0.000579f, 0.000475f, -0.000080f, -0.000320f, 0.000029f, -0.000276f, -0.000459f, -0.000085f, -0.000411f, -0.000491f, -0.000295f, -0.000577f, -0.000510f, -0.000608f, -0.000424f, -0.000586f, -0.000678f, -0.000568f, -0.000680f, -0.000631f, -0.000846f, -0.000506f, -0.000668f, -0.000868f, -0.000698f, -0.000642f, -0.000779f, -0.000866f, -0.000735f, -0.000579f, -0.000845f, -0.000929f, -0.000678f, -0.000666f, -0.000863f, -0.000861f, -0.000755f, -0.000630f, -0.000844f, -0.000865f, -0.000819f, -0.000551f, -0.000714f, 0.006893f, 0.006332f, 0.004748f, 0.004110f, 0.003164f, 0.002842f, 0.002249f, 0.002309f, 0.001470f, 0.001583f, 0.000918f, 0.001345f, 0.000527f, 0.000966f, 0.000116f, 0.000678f, 0.000056f, 0.000240f, 0.000168f, -0.000499f, 0.000434f, -0.000808f, 0.000253f, -0.000631f, -0.000413f, -0.000343f, -0.000589f, -0.000628f, -0.000385f, -0.000857f, -0.000591f, -0.000694f, -0.000721f, -0.000793f, -0.000803f, -0.000720f, -0.000965f, -0.000835f, -0.000821f, -0.000845f, -0.000955f, -0.000938f, -0.000788f, -0.001064f, -0.000944f, -0.000880f, -0.000899f, -0.001083f, -0.001002f, -0.000813f, -0.000928f, -0.001144f, -0.000901f, -0.000911f, -0.001025f, -0.001037f, -0.000913f, -0.000843f, -0.000997f, -0.001085f, -0.000980f, -0.000770f, -0.000989f, -0.000980f, -0.000995f, -0.000806f, -0.000969f, -0.000906f, -0.000961f, -0.000852f, -0.000814f, -0.000952f, -0.000938f, -0.000862f, -0.000744f, -0.000867f, -0.000863f, -0.000896f, -0.000815f, -0.000794f, -0.000727f, -0.000855f, -0.000770f, -0.000802f, -0.000715f, -0.000798f, -0.000750f, -0.000700f, -0.000750f, -0.000678f, -0.000710f, -0.000717f, -0.000687f, -0.000678f, -0.000624f, -0.000651f, -0.000669f, -0.000613f, -0.000639f, -0.000618f, -0.000570f, -0.000583f, -0.000573f, -0.000597f, -0.000582f, -0.000538f, -0.000512f, -0.000546f, -0.000458f, -0.000589f, -0.000443f, -0.000594f, -0.000273f, -0.000835f, 0.001001f, 0.007335f, 0.005023f, 0.004543f, 0.003403f, 0.003225f, 0.002397f, 0.002419f, 0.001751f, 0.001896f, 0.001281f, 0.001453f, 0.000858f, 0.001120f, 0.000584f, 0.000876f, 0.000458f, 0.000427f, 0.000514f, -0.000055f, 0.000619f, -0.000327f, 0.000483f, -0.000239f, -0.000010f, -0.000008f, -0.000301f, -0.000113f, -0.000129f, -0.000434f, -0.000175f, -0.000388f, -0.000345f, -0.000415f, -0.000302f, -0.000475f, -0.000482f, -0.000437f, -0.000484f, -0.000518f, -0.000560f, -0.000332f, -0.000688f, -0.000553f, -0.000519f, -0.000513f, -0.000694f, -0.000536f, -0.000487f, -0.000594f, -0.000782f, -0.000429f, -0.000615f, -0.000619f, -0.000725f, -0.000425f, -0.000641f, -0.000592f, -0.000802f, -0.000357f, -0.000723f, -0.000517f, -0.000835f, -0.000265f, -0.000935f, 0.004990f, 0.005381f, 0.003737f, 0.003416f, 0.002547f, 0.002215f, 0.001890f, 0.001657f, 0.001344f, 0.001069f, 0.000963f, 0.000729f, 0.000731f, 0.000434f, 0.000367f, 0.000305f, 0.000094f, 0.000321f, -0.000343f, 0.000372f, -0.000614f, 0.000280f, -0.000557f, -0.000130f, -0.000356f, -0.000453f, -0.000357f, -0.000371f, -0.000655f, -0.000320f, -0.000687f, -0.000496f, -0.000657f, -0.000461f, -0.000760f, -0.000558f, -0.000735f, -0.000561f, -0.000841f, -0.000553f, -0.000730f, -0.000797f, -0.000708f, -0.000618f, -0.000811f, -0.000769f, -0.000716f, -0.000684f, -0.000917f, -0.000638f, -0.000703f, -0.000778f, -0.000882f, -0.000682f, -0.000696f, -0.000751f, -0.000885f, -0.000622f, -0.000753f, -0.000783f, -0.000772f, -0.000657f, -0.000721f, -0.000754f, -0.000815f, -0.000655f, -0.000629f, -0.000741f, -0.000754f, -0.000685f, -0.000656f, -0.000688f, -0.000678f, -0.000675f, -0.000649f, -0.000625f, -0.000694f, -0.000652f, -0.000612f, -0.000592f, -0.000623f, -0.000642f, -0.000585f, -0.000626f, -0.000537f, -0.000582f, -0.000578f, -0.000576f, -0.000567f, -0.000521f, -0.000541f, -0.000522f, -0.000510f, -0.000552f, -0.000489f, -0.000502f, -0.000475f, -0.000495f, -0.000475f, -0.000459f, -0.000441f, -0.000467f, -0.000418f, -0.000477f, -0.000405f, -0.000420f, -0.000375f, -0.000428f, -0.000380f, -0.000436f, -0.000323f, -0.000411f, -0.000268f, -0.000537f, 0.000252f, 0.005587f, 0.004406f, 0.003503f, 0.002895f, 0.002430f, 0.002063f, 0.001801f, 0.001571f, 0.001311f, 0.001151f, 0.000966f, 0.000843f, 0.000761f, 0.000586f, 0.000582f, 0.000292f, 0.000514f, 0.000008f, 0.000532f, -0.000182f, 0.000418f, -0.000228f, 0.000128f, -0.000082f, -0.000154f, -0.000013f, -0.000210f, -0.000214f, -0.000156f, -0.000311f, -0.000239f, -0.000272f, -0.000279f, -0.000394f, -0.000279f, -0.000380f, -0.000375f, -0.000424f, -0.000255f, -0.000569f, -0.000342f, -0.000404f, -0.000416f, -0.000552f, -0.000313f, -0.000457f, -0.000544f, -0.000457f, -0.000336f, -0.000568f, -0.000466f, -0.000457f, -0.000370f, -0.000613f, -0.000448f, -0.000438f, -0.000414f, -0.000611f, -0.000372f, -0.000518f, -0.000338f, -0.000700f, -0.000224f, -0.000861f, 0.003401f, 0.004657f, 0.002907f, 0.002851f, 0.001989f, 0.001877f, 0.001392f, 0.001365f, 0.001028f, 0.000928f, 0.000661f, 0.000648f, 0.000491f, 0.000370f, 0.000301f, 0.000152f, 0.000231f, -0.000067f, 0.000182f, -0.000293f, 0.000109f, -0.000333f, -0.000069f, -0.000266f, -0.000325f, -0.000206f, -0.000436f, -0.000304f, -0.000387f, -0.000413f, -0.000467f, -0.000389f, -0.000497f, -0.000479f, -0.000473f, -0.000518f, -0.000588f, -0.000441f, -0.000526f, -0.000629f, -0.000508f, -0.000562f, -0.000603f, -0.000575f, -0.000469f, -0.000671f, -0.000631f, -0.000499f, -0.000577f, -0.000661f, -0.000575f, -0.000533f, -0.000633f, -0.000629f, -0.000504f, -0.000572f, -0.000656f, -0.000615f, -0.000511f, -0.000551f, -0.000613f, -0.000600f, -0.000513f, -0.000573f, -0.000596f, -0.000542f, -0.000517f, -0.000550f, -0.000574f, -0.000546f, -0.000507f, -0.000488f, -0.000564f, -0.000523f, -0.000516f, -0.000474f, -0.000512f, -0.000481f, -0.000500f, -0.000482f, -0.000467f, -0.000479f, -0.000458f, -0.000463f, -0.000433f, -0.000453f, -0.000444f, -0.000422f, -0.000450f, -0.000391f, -0.000417f, -0.000400f, -0.000420f, -0.000397f, -0.000378f, -0.000372f, -0.000369f, -0.000374f, -0.000382f, -0.000339f, -0.000348f, -0.000330f, -0.000353f, -0.000340f, -0.000322f, -0.000291f, -0.000325f, -0.000305f, -0.000320f, -0.000280f, -0.000286f, -0.000255f, -0.000308f, -0.000260f, -0.000214f, 0.004212f, 0.003787f, 0.002745f, 0.002410f, 0.001895f, 0.001723f, 0.001397f, 0.001261f, 0.001016f, 0.000936f, 0.000717f, 0.000742f, 0.000502f, 0.000550f, 0.000285f, 0.000458f, 0.000109f, 0.000390f, -0.000038f, 0.000278f, -0.000089f, 0.000106f, -0.000029f, -0.000083f, 0.000014f, -0.000206f, -0.000076f, -0.000182f, -0.000176f, -0.000190f, -0.000174f, -0.000286f, -0.000220f, -0.000256f, -0.000291f, -0.000309f, -0.000191f, -0.000413f, -0.000270f, -0.000299f, -0.000345f, -0.000396f, -0.000226f, -0.000402f, -0.000390f, -0.000294f, -0.000336f, -0.000452f, -0.000313f, -0.000332f, -0.000376f, -0.000452f, -0.000255f, -0.000410f, -0.000381f, -0.000414f, -0.000262f, -0.000443f, -0.000331f, -0.000452f, -0.000192f, -0.000560f, -0.000147f, -0.000788f, 0.002227f, 0.003961f, 0.002324f, 0.002309f, 0.001612f, 0.001549f, 0.001075f, 0.001105f, 0.000771f, 0.000780f, 0.000472f, 0.000580f, 0.000272f, 0.000426f, 0.000072f, 0.000302f, -0.000080f, 0.000229f, -0.000197f, 0.000132f, -0.000244f, -0.000036f, -0.000207f, -0.000195f, -0.000175f, -0.000313f, -0.000227f, -0.000280f, -0.000338f, -0.000294f, -0.000367f, -0.000370f, -0.000349f, -0.000345f, -0.000458f, -0.000375f, -0.000392f, -0.000467f, -0.000378f, -0.000422f, -0.000500f, -0.000426f, -0.000371f, -0.000528f, -0.000447f, -0.000405f, -0.000520f, -0.000448f, -0.000418f, -0.000442f, -0.000553f, -0.000426f, -0.000411f, -0.000493f, -0.000488f, -0.000406f, -0.000468f, -0.000478f, -0.000459f, -0.000388f, -0.000480f, -0.000472f, -0.000439f, -0.000384f, -0.000460f, -0.000436f, -0.000448f, -0.000386f, -0.000446f, -0.000408f, -0.000420f, -0.000375f, -0.000434f, -0.000391f, -0.000407f, -0.000371f, -0.000387f, -0.000384f, -0.000389f, -0.000363f, -0.000365f, -0.000354f, -0.000362f, -0.000350f, -0.000368f, -0.000316f, -0.000354f, -0.000309f, -0.000363f, -0.000292f, -0.000336f, -0.000286f, -0.000340f, -0.000287f, -0.000312f, -0.000260f, -0.000308f, -0.000270f, -0.000315f, -0.000228f, -0.000291f, -0.000235f, -0.000299f, -0.000222f, -0.000282f, -0.000188f, -0.000290f, -0.000185f, -0.000301f, -0.000141f, -0.000289f, -0.000136f, -0.000323f, -0.000067f, -0.000408f, 0.003050f, 0.003229f, 0.002196f, 0.001959f, 0.001542f, 0.001386f, 0.001126f, 0.000969f, 0.000839f, 0.000692f, 0.000643f, 0.000498f, 0.000480f, 0.000315f, 0.000363f, 0.000184f, 0.000280f, 0.000076f, 0.000182f, 0.000024f, 0.000063f, 0.000023f, -0.000061f, 0.000027f, -0.000152f, -0.000020f, -0.000148f, -0.000118f, -0.000133f, -0.000157f, -0.000203f, -0.000153f, -0.000206f, -0.000234f, -0.000157f, -0.000271f, -0.000247f, -0.000194f, -0.000281f, -0.000283f, -0.000179f, -0.000310f, -0.000297f, -0.000211f, -0.000297f, -0.000327f, -0.000228f, -0.000268f, -0.000331f, -0.000286f, -0.000219f, -0.000364f, -0.000272f, -0.000283f, -0.000242f, -0.000370f, -0.000237f, -0.000306f, -0.000243f, -0.000392f, -0.000154f, -0.000409f, -0.000096f, -0.000659f, 0.001342f, 0.003314f, 0.001923f, 0.001853f, 0.001308f, 0.001239f, 0.000886f, 0.000850f, 0.000629f, 0.000597f, 0.000424f, 0.000384f, 0.000279f, 0.000236f, 0.000172f, 0.000110f, 0.000081f, -0.000001f, 0.000012f, -0.000051f, -0.000100f, -0.000080f, -0.000183f, -0.000069f, -0.000271f, -0.000120f, -0.000292f, -0.000183f, -0.000264f, -0.000230f, -0.000335f, -0.000218f, -0.000339f, -0.000279f, -0.000303f, -0.000350f, -0.000330f, -0.000263f, -0.000409f, -0.000310f, -0.000334f, -0.000385f, -0.000336f, -0.000302f, -0.000425f, -0.000351f, -0.000324f, -0.000356f, -0.000419f, -0.000289f, -0.000397f, -0.000377f, -0.000344f, -0.000319f, -0.000406f, -0.000364f, -0.000334f, -0.000329f, -0.000400f, -0.000324f, -0.000347f, -0.000347f, -0.000365f, -0.000313f, -0.000337f, -0.000339f, -0.000359f, -0.000301f, -0.000329f, -0.000320f, -0.000333f, -0.000309f, -0.000315f, -0.000302f, -0.000315f, -0.000296f, -0.000303f, -0.000288f, -0.000300f, -0.000279f, -0.000296f, -0.000260f, -0.000290f, -0.000260f, -0.000295f, -0.000237f, -0.000274f, -0.000239f, -0.000279f, -0.000232f, -0.000256f, -0.000216f, -0.000261f, -0.000225f, -0.000239f, -0.000198f, -0.000246f, -0.000202f, -0.000239f, -0.000178f, -0.000220f, -0.000183f, -0.000233f, -0.000176f, -0.000200f, -0.000153f, -0.000221f, -0.000146f, -0.000225f, -0.000110f, -0.000230f, -0.000085f, -0.000273f, -0.000007f, -0.000420f, 0.002082f, 0.002789f, 0.001726f, 0.001664f, 0.001190f, 0.001155f, 0.000841f, 0.000831f, 0.000608f, 0.000619f, 0.000444f, 0.000452f, 0.000295f, 0.000325f, 0.000203f, 0.000234f, 0.000125f, 0.000135f, 0.000075f, 0.000046f, 0.000048f, -0.000034f, 0.000039f, -0.000100f, -0.000004f, -0.000115f, -0.000072f, -0.000104f, -0.000116f, -0.000144f, -0.000102f, -0.000199f, -0.000123f, -0.000164f, -0.000207f, -0.000146f, -0.000186f, -0.000237f, -0.000138f, -0.000225f, -0.000235f, -0.000159f, -0.000229f, -0.000250f, -0.000172f, -0.000220f, -0.000262f, -0.000193f, -0.000202f, -0.000280f, -0.000196f, -0.000213f, -0.000229f, -0.000266f, -0.000173f, -0.000253f, -0.000228f, -0.000240f, -0.000168f, -0.000283f, -0.000174f, -0.000296f, -0.000084f, -0.000467f, 0.000699f, 0.002711f, 0.001613f, 0.001503f, 0.001042f, 0.001006f, 0.000717f, 0.000687f, 0.000482f, 0.000483f, 0.000319f, 0.000327f, 0.000196f, 0.000203f, 0.000107f, 0.000114f, 0.000023f, 0.000022f, -0.000011f, -0.000044f, -0.000048f, -0.000128f, -0.000071f, -0.000165f, -0.000101f, -0.000195f, -0.000168f, -0.000173f, -0.000201f, -0.000214f, -0.000204f, -0.000265f, -0.000187f, -0.000247f, -0.000273f, -0.000222f, -0.000282f, -0.000261f, -0.000227f, -0.000310f, -0.000280f, -0.000225f, -0.000315f, -0.000276f, -0.000253f, -0.000309f, -0.000298f, -0.000223f, -0.000320f, -0.000296f, -0.000266f, -0.000265f, -0.000315f, -0.000266f, -0.000265f, -0.000303f, -0.000284f, -0.000244f, -0.000288f, -0.000283f, -0.000278f, -0.000244f, -0.000281f, -0.000268f, -0.000261f, -0.000250f, -0.000275f, -0.000247f, -0.000258f, -0.000239f, -0.000263f, -0.000240f, -0.000250f, -0.000224f, -0.000251f, -0.000224f, -0.000248f, -0.000213f, -0.000236f, -0.000208f, -0.000242f, -0.000194f, -0.000228f, -0.000196f, -0.000227f, -0.000190f, -0.000209f, -0.000180f, -0.000219f, -0.000183f, -0.000195f, -0.000165f, -0.000202f, -0.000168f, -0.000197f, -0.000146f, -0.000191f, -0.000149f, -0.000195f, -0.000127f, -0.000180f, -0.000135f, -0.000183f, -0.000121f, -0.000169f, -0.000109f, -0.000182f, -0.000102f, -0.000182f, -0.000065f, -0.000195f, -0.000051f, -0.000233f, 0.000031f, -0.000385f, 0.001375f, 0.002357f, 0.001401f, 0.001358f, 0.000960f, 0.000926f, 0.000677f, 0.000662f, 0.000492f, 0.000481f, 0.000355f, 0.000344f, 0.000244f, 0.000253f, 0.000161f, 0.000178f, 0.000085f, 0.000121f, 0.000024f, 0.000083f, -0.000028f, 0.000051f, -0.000071f, 0.000010f, -0.000084f, -0.000043f, -0.000074f, -0.000087f, -0.000099f, -0.000087f, -0.000149f, -0.000075f, -0.000159f, -0.000125f, -0.000117f, -0.000177f, -0.000137f, -0.000128f, -0.000204f, -0.000117f, -0.000169f, -0.000194f, -0.000142f, -0.000161f, -0.000210f, -0.000139f, -0.000172f, -0.000208f, -0.000155f, -0.000166f, -0.000195f, -0.000183f, -0.000146f, -0.000209f, -0.000174f, -0.000169f, -0.000163f, -0.000206f, -0.000152f, -0.000189f, -0.000152f, -0.000225f, -0.000088f, -0.000302f, 0.000294f, 0.002134f, 0.001392f, 0.001199f, 0.000873f, 0.000785f, 0.000592f, 0.000536f, 0.000401f, 0.000364f, 0.000265f, 0.000253f, 0.000158f, 0.000152f, 0.000079f, 0.000092f, 0.000011f, 0.000037f, -0.000054f, 0.000006f, -0.000094f, -0.000035f, -0.000137f, -0.000058f, -0.000143f, -0.000115f, -0.000150f, -0.000151f, -0.000150f, -0.000156f, -0.000204f, -0.000146f, -0.000214f, -0.000170f, -0.000194f, -0.000229f, -0.000179f, -0.000196f, -0.000247f, -0.000167f, -0.000253f, -0.000213f, -0.000191f, -0.000230f, -0.000249f, -0.000177f, -0.000250f, -0.000220f, -0.000211f, -0.000219f, -0.000254f, -0.000185f, -0.000225f, -0.000236f, -0.000219f, -0.000201f, -0.000235f, -0.000211f, -0.000213f, -0.000205f, -0.000234f, -0.000194f, -0.000210f, -0.000206f, -0.000217f, -0.000193f, -0.000207f, -0.000191f, -0.000211f, -0.000184f, -0.000201f, -0.000185f, -0.000197f, -0.000179f, -0.000189f, -0.000171f, -0.000196f, -0.000168f, -0.000179f, -0.000164f, -0.000181f, -0.000164f, -0.000171f, -0.000150f, -0.000174f, -0.000154f, -0.000158f, -0.000142f, -0.000165f, -0.000142f, -0.000156f, -0.000124f, -0.000156f, -0.000132f, -0.000155f, -0.000107f, -0.000148f, -0.000115f, -0.000150f, -0.000104f, -0.000137f, -0.000100f, -0.000141f, -0.000096f, -0.000131f, -0.000083f, -0.000143f, -0.000073f, -0.000141f, -0.000053f, -0.000153f, -0.000035f, -0.000182f, 0.000031f, -0.000300f, 0.000843f, 0.001958f, 0.001156f, 0.001106f, 0.000778f, 0.000746f, 0.000543f, 0.000529f, 0.000390f, 0.000381f, 0.000284f, 0.000268f, 0.000200f, 0.000187f, 0.000137f, 0.000119f, 0.000088f, 0.000069f, 0.000048f, 0.000023f, 0.000019f, -0.000009f, -0.000017f, -0.000025f, -0.000054f, -0.000029f, -0.000087f, -0.000046f, -0.000091f, -0.000085f, -0.000075f, -0.000115f, -0.000089f, -0.000100f, -0.000142f, -0.000073f, -0.000151f, -0.000111f, -0.000107f, -0.000152f, -0.000122f, -0.000110f, -0.000168f, -0.000110f, -0.000129f, -0.000161f, -0.000125f, -0.000127f, -0.000159f, -0.000132f, -0.000124f, -0.000161f, -0.000136f, -0.000129f, -0.000144f, -0.000148f, -0.000127f, -0.000140f, -0.000145f, -0.000138f, -0.000121f, -0.000155f, -0.000113f, -0.000181f, 0.000055f, 0.001645f, 0.001209f, 0.000946f, 0.000736f, 0.000614f, 0.000487f, 0.000416f, 0.000340f, 0.000276f, 0.000220f, 0.000182f, 0.000142f, 0.000106f, 0.000074f, 0.000046f, 0.000037f, -0.000003f, -0.000015f, -0.000037f, -0.000033f, -0.000067f, -0.000075f, -0.000083f, -0.000098f, -0.000087f, -0.000137f, -0.000104f, -0.000132f, -0.000130f, -0.000130f, -0.000165f, -0.000130f, -0.000144f };

    struct inline4 : as_engine<
        /* W             */ 4,
        /* H             */ 9,
        /* THROTTLE_Y    */ 2,
        /* PISTON_Y      */ 4,
        /* AUDIO_Y       */ 5,
        /* PIPE_CELLS    */ 256,
        /* PIPE_SUBSTEPS */ 10,
        inline_pistons,
        vtec_cams,
        basic_sparkplugs>
    {
        inline4()
        {
            this->convolution.set_impulse(g_impulse);
            this->lumped_drag_torque_n_m = 35.2;
            this->limiter.max_angular_velocity_r_per_s = 700.0;
            this->limiter.limit_time_s = 0.04;
            this->flywheel.mass_kg = 18.0;
            this->flywheel.radius_m = 0.20;
            this->crankshaft.mass_kg = 12.5;
            this->crankshaft.radius_m = 0.045;
            this->crankshaft.angular_velocity_r_per_s = 150.0;
            this->pistons.diameter_m.fill(0.086);
            this->pistons.crank_throw_length_m.fill(0.043);
            this->pistons.connecting_rod_length_m.fill(0.145);
            this->pistons.connecting_rod_mass_kg.fill(0.45);
            this->pistons.head_mass_density_kg_per_m3.fill(2700.0);
            this->pistons.head_compression_height_m.fill(0.030);
            this->pistons.head_clearance_height_m.fill(0.007);
            this->pistons.friction_n_m_s2_per_r2.fill(0.00005);
            this->inlet_cam.ramp_theta_r.fill(g_pi_r * 0.85);
            this->inlet_cam.vtec_engage_r_per_s = 450.0;
            this->outlet_cam.ramp_theta_r.fill(g_pi_r * 0.5);
            this->outlet_cam.vtec_engage_r_per_s = 450.0;
            double theta0_r = 0.0;
            for(size_t i = 0; i < get_width(); i++)
            {
                this->pistons.theta0_r[i] = theta0_r;
                this->inlet_cam.engage_theta_r[i]  = theta0_r + g_otto_intake_cycle_r - 1.0;
                this->sparkplugs.engage_theta_r[i] = theta0_r + g_otto_combustion_cycle_r - 0.4;
                this->outlet_cam.engage_theta_r[i] = theta0_r + g_otto_exhaust_cycle_r + 0.9;
                theta0_r += g_otto_cycle_r / get_width();
            }
            for(auto& flow : this->flows)
            {
                flow.chamber_nozzle_open_ratio.fill(1.0);
                flow.chamber_nozzle_flow_area_m2 = {
                    0.00250, /* Source   -> Intake    */
                    0.00120, /* Intake   -> Throttle  */
                    0.00085, /* Throttle -> Runner    */
                    0.00090, /* Runner   -> Piston    */
                    0.00120, /* Piston   -> Chamber0  */
                    0.00170, /* Chamber0 -> Chamber1  */
                    0.00185, /* Chamber1 -> Chamber2  */
                    0.00250, /* Chamber2 -> Sink      */
                };
                flow.chamber_volume_m3 = {
                    g_resevoir_volume_m3, /* Source   */
                    0.0030,               /* Intake   */
                    0.0008,               /* Throttle */
                    0.0003,               /* Runner   */
                    0.0000,               /* Piston   */
                    0.0003,               /* Chamber0 */
                    0.0004,               /* Chamber1 */
                    0.0005,               /* Chamber2 */
                    g_resevoir_volume_m3  /* Sink     */
                };
            }
            this->throttle.table = {
                0.001,
                0.025,
                0.250,
                1.000,
            };
            this->pipe.piston_connect_m = { 0.00, 0.38, 0.16, 0.26 };
            this->pipe.mic_position0_m = 0.61;
            this->pipe.mic_position1_m = 0.74;
            this->pipe.length_m = 1.0;
            this->dc.set_cutoff_frequency(10.0);
            this->gain.ratio = 0.000002;
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
