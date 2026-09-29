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

    using std::sin, std::cos, std::fmax, std::fmin, std::log, std::sqrt, std::trunc, std::exp, std::fabs;

    fn auto clamper(const auto value, const auto lower, const auto upper)
    {
        return fmax(fmin(value, upper), lower);
    }

    fn auto modulos(const auto value, const auto by)
    {
        return value - trunc(value / by) * by;
    }

    fn auto cuberoot(const auto value)
    {
        return exp(log(value) / 3.0);
    }

    fn double frand()
    {
        return 2.0 * rand() / RAND_MAX - 1.0;
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

                const double dh = S * g_dt_s;
                const double drh = dh / 2.0 + dh / 2.0 * frand();
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

    template<size_t W, size_t S = 0>
    struct cams
    {
        std::array<double, W> engage_theta_r = {};
        std::array<double, W> temp_ramp_theta_r = {};
        std::array<double, W> temp_open_ratio = {};
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
                const double t = open_r / temp_ramp_theta_r[i];
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

        fn virtual void set_ramp_thetas()
        {
            for(size_t i = 0; i < W; i++)
            {
                temp_ramp_theta_r[i] = ramp_theta_r[i];
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
            set_ramp_thetas();
            calc_open_ratios();
            set_open_ratios();
        }
    };

    template<size_t W, size_t S>
    requires(S > 0)
    struct vtec_cams : cams<W>
    {
        std::array<double, S> vtec_engage_r_per_s = {};
        std::array<double, S> vtec_open_boost = {};
        std::array<double, S> vtec_ramp_boost = {};

        fn void set_ramp_thetas() override
        {
            for(size_t i = 0; i < W; i++)
            for(size_t j = 0; j < S; j++)
            {
                if(this->flywheel_angular_velocity_r_per_s > vtec_engage_r_per_s[j])
                {
                    this->temp_ramp_theta_r[i] = vtec_ramp_boost[j] * this->ramp_theta_r[i];
                }
            }
        }

        fn void set_open_ratios() override
        {
            for(size_t i = 0; i < W; i++)
            for(size_t j = 0; j < S; j++)
            {
                if(this->flywheel_angular_velocity_r_per_s > vtec_engage_r_per_s[j])
                {
                    this->open_ratio[i] = clamper(vtec_open_boost[j] * this->temp_open_ratio[i], 0.0, 1.0);
                }
            }
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
        double max_angular_velocity_r_per_s = 0.0;
        double flywheel_angular_velocity_r_per_s = 0.0;
        double limit_time_s = 0.1;
        double cycles = 0;
        bool limiting = false;

        void update()
        {
            if(not limiting)
            {
                if(flywheel_angular_velocity_r_per_s > max_angular_velocity_r_per_s)
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

    struct disk
    {
        double angular_velocity_r_per_s = 0.0;
        double mass_kg = 0.0;
        double radius_m = 0.0;
        double moment_of_inertia_kg_m2 = 0.0;

        /*
         * dw = a * dt
         *
         */

        fn void accelerate(const double angular_acceleration_r_per_s_s)
        {
            const double a = angular_acceleration_r_per_s_s;
            angular_velocity_r_per_s += a * g_dt_s;
            angular_velocity_r_per_s = fmax(angular_velocity_r_per_s, 0.0);
        }

        fn void calc_moment_of_inertia()
        {
            moment_of_inertia_kg_m2 = 0.5 * mass_kg * radius_m * radius_m;
        }
    };

    struct load : disk
    {
        double friction_n_m_s2_per_r2 = 0.0;
        double friction_torque_n_m = 0.0;
        double total_torque_n_m = 0.0;

        /*
         * Tf = -K w
         *
         */

        fn void calc_friction_torque()
        {
            const double K = friction_n_m_s2_per_r2;
            const double w = angular_velocity_r_per_s;
            friction_torque_n_m = -K * w;
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

    struct flywheel : disk
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

        double calc_torque(const bool clutch_engaged, double wl, double wr)
        {
            const double Tk = clutch_engaged ? damping_coefficient_n_m_s : 0.0;
            const double Tc = Tk * (wl - wr);
            return Tc;
        }
    };

    template<size_t N>
    struct gearbox
    {
        static constexpr size_t last_gear = N - 1;
        std::array<double, N> table = {};
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
        float mic_position_m = 0.0f;

        pipe()
        {
            reset();
        }

        /*
         * --- +-------+
         *  |  |  0    | ---+
         *  |  +-------+    |
         *  |  +-------+    |
         *  WW |  1    | ---+ --> U0 ... UL-1
         *  |  +-------+    |
         *  |     ...       |
         *  |  +-------+    |
         *  |  |  WW-1 | ---+
         * --- +-------+
         */

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
            const size_t x = Z * mic_position_m / length_m;
            return static_pressure_pa[x];
        }

        fn void inject()
        {
            for(size_t i = 0; i < WW; i++)
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
            for(size_t i = 0; i < S; i++)
            {
                calc_static_pressures();
                calc_speed_of_sounds();
                calc_local_speed_of_sounds();
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

    #define FLUIDS(X)                       \
        X(chamber_volume_m3)                \
        X(chamber_nozzle_real_flow_area_m2) \
        X(chamber_static_pressure_pa)       \
        X(chamber_static_temperature_k)     \
        X(chamber_mass_kg)                  \
        X(nozzle_static_temperature_k)      \
        X(nozzle_static_density_kg_per_m3)  \
        X(nozzle_mach)

    #define PISTONS(X) X(total_torque_n_m)

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
        static constexpr size_t size = 16384;
        static constexpr size_t mask = size - 1;

        std::vector<float> history = {};
        size_t head = 0;

        convolution_filter()
        {
            history.resize(2 * size);
        }

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
        std::atomic<bool> clutch_engaged = true;
        std::atomic<size_t> gear = 0;

        /*
         * Send
         */

        std::atomic<size_t> swap_drops = 0;
        std::atomic<double> load_angular_velocity_r_per_s = 0.0;
        std::atomic<double> engine_angular_velocity_r_per_s = 0.0;
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
        size_t PIPE_COUNT,
        size_t VTEC_STEPS,
        size_t GEARS,
        template<size_t> typename PISTONS,
        template<size_t, size_t> typename CAMS,
        template<size_t> typename SPARKPLUGS>
    struct as_engine : engine
    {
        double lumped_parasitic_torque_n_m = {};
        struct PISTONS<W> pistons = {};
        struct CAMS<W, VTEC_STEPS> inlet_cam = {};
        struct CAMS<W, VTEC_STEPS> outlet_cam = {};
        struct SPARKPLUGS<W> sparkplugs = {};
        std::array<struct flow<H, PISTON_Y>, W> flows = {};
        struct limiter limiter = {};
        struct throttle throttle = {};
        struct flywheel flywheel = {};
        struct dc_filter dc = {};
        struct gain_filter gain = {};
        struct clamp_filter clamp = {};
        struct convolution_filter convolution = {};
        struct gearbox<GEARS> gearbox = {};
        struct clutch clutch = {};
        struct load load = {};
        struct diags diags = {};
        static constexpr size_t WW = W / PIPE_COUNT;
        std::array<struct pipe<WW, PIPE_CELLS, PIPE_SUBSTEPS>, PIPE_COUNT> pipes = {};
        std::vector<float> audio_signal = {};
        struct mailbox<W, H> mailbox = {};
        std::array<std::vector<float>, PIPE_COUNT> pipe_static_pressures_pa = {};
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

        fn void broadcast(const double throttle_open_ratio, const bool injection_enabled)
        {
            /*
             * flywheel theta -> inlet/outlet cams + pistons + sparkplugs thetas.
             *
             */

            inlet_cam.flywheel_theta_r = flywheel.theta_r;
            outlet_cam.flywheel_theta_r = flywheel.theta_r;
            pistons.flywheel_theta_r = flywheel.theta_r;
            sparkplugs.flywheel_theta_r = flywheel.theta_r;
            pistons.flywheel_angular_velocity_r_per_s = flywheel.angular_velocity_r_per_s;
            limiter.flywheel_angular_velocity_r_per_s = flywheel.angular_velocity_r_per_s;
            inlet_cam.flywheel_angular_velocity_r_per_s = flywheel.angular_velocity_r_per_s;
            outlet_cam.flywheel_angular_velocity_r_per_s = flywheel.angular_velocity_r_per_s;

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
            pistons.calc_volumetrics();
            broadcast(0.0, false);
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
                size_t i = 0;
                for(auto& pipe : pipes)
                {
                    auto& pressure = pipe.static_pressure_pa;
                    pipe_static_pressures_pa[i++].assign(pressure.begin(), pressure.end());
                }
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
                const size_t i = x / WW;
                const size_t j = x % WW;
                pipes[i].in_velocity_m_per_s[j] = flows[x].nozzle_velocity_m_per_s[AUDIO_Y];
                pipes[i].in_static_temperature_k[j] = flows[x].nozzle_static_temperature_k[AUDIO_Y];
                pipes[i].in_static_density_kg_per_m3[j] = flows[x].nozzle_static_density_kg_per_m3[AUDIO_Y];
            }
            for(auto& pipe : pipes)
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

        void update_flywheel_and_load(const size_t gear, const bool clutch_engaged)
        {
            flywheel.update();
            load.update();
            const double N = gearbox.table[gear];
            const double Il = calc_engine_moment_of_inertia();
            const double Tl = calc_engine_torque();
            const double wl = flywheel.angular_velocity_r_per_s;
            const double Ir = load.moment_of_inertia_kg_m2;
            const double Tr = load.total_torque_n_m;
            const double wr = load.angular_velocity_r_per_s;
            const double Tc = clutch.calc_torque(clutch_engaged, wl, N * wr);
            const double al = (Tl - Tc) / Il;
            const double ar = (Tr + N * Tc) / Ir;
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

        void post_mailbox(const size_t swap_drops)
        {
            mailbox.load_angular_velocity_r_per_s = load.angular_velocity_r_per_s;
            mailbox.engine_angular_velocity_r_per_s = flywheel.angular_velocity_r_per_s;
            for(size_t y = 0; y < H; y++)
            for(size_t x = 0; x < W; x++)
            {
                mailbox.port_open_ratios[y][x] = flows[x].chamber_nozzle_open_ratio[y];
                mailbox.panics[y][x] = flows[x].panic[y];
            }
            mailbox.swap_drops += swap_drops;
        }

        float calc_audio_sample()
        {
            float x = 0.0f;
            for(auto& pipe : pipes)
            {
                x += pipe.calc_audio_sample();
            }
            x = dc.filter(x);
            x = convolution.filter(x);
            x = gain.filter(x);
            x = clamp.filter(x);
            return x;
        }

        void run(const size_t steps) override
        {
            const double throttle_open_ratio = mailbox.throttle_open_ratio;
            const size_t log_x = mailbox.log_x;
            const size_t log_y = mailbox.log_y;
            const bool injection_enabled = mailbox.injection_enabled;
            const bool clutch_engaged = mailbox.clutch_engaged;
            const size_t gear = mailbox.gear;
            audio_signal.clear();
            size_t swap_drops = 0;
            for(size_t step = 0; step < steps; step++)
            {
                update_limiter();
                update_flywheel_and_load(gear, clutch_engaged);
                if(flywheel.otto_cycled())
                {
                    diags_swap();
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
                broadcast(throttle_open_ratio, injection_overrided);
                const float sample = calc_audio_sample();
                audio_signal.push_back(sample);
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

        const std::atomic<double>& get_load_angular_velocity_r_per_s() const override
        {
            return mailbox.load_angular_velocity_r_per_s;
        }

        const std::atomic<double>& get_engine_angular_velocity_r_per_s() const override
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

        const std::span<const std::vector<float>> get_pipe_pressure_signals() const override
        {
            return pipe_static_pressures_pa;
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

        void engage_clutch() override
        {
            mailbox.clutch_engaged = true;
        }

        void disengage_clutch() override
        {
            mailbox.clutch_engaged = false;
        }

        void increment_gear() override
        {
            if(mailbox.gear < gearbox.last_gear)
            {
                mailbox.gear++;
            }
        }

        void decrement_gear() override
        {
            if(mailbox.gear > 0)
            {
                mailbox.gear--;
            }
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
    };

    const std::vector<float> g_impulse =
;

    struct inline8 : as_engine<
          8, /* W             */
          9, /* H             */
          2, /* THROTTLE_Y    */
          4, /* PISTON_Y      */
          7, /* AUDIO_Y       */
        192, /* PIPE_CELLS    */
          8, /* PIPE_SUBSTEPS */
          2, /* PIPE_COUNT    */
          5, /* VTEC_STEPS    */
          6, /* GEARS         */
        inline_pistons,
        vtec_cams,
        sparkplugs>
    {
        inline8()
        {
            this->gearbox.table = { 7.0, 6.0, 5.0, 4.25, 3.7, 3.0 };
            this->clutch.damping_coefficient_n_m_s = 7.0;
            this->convolution.set_impulse(g_impulse);
            this->lumped_parasitic_torque_n_m = 15.0;
            this->pistons.friction_n_m_s2_per_r2.fill(0.00003);
            this->limiter.max_angular_velocity_r_per_s = 1350.0;
            this->limiter.limit_time_s = 0.02;
            this->load.mass_kg = 1000.0;
            this->load.radius_m = 0.225;
            this->load.friction_n_m_s2_per_r2 = 1.0;
            this->flywheel.mass_kg = 18.5;
            this->flywheel.radius_m = 0.19;
            this->flywheel.angular_velocity_r_per_s = this->load.angular_velocity_r_per_s = 100.0;
            this->pistons.diameter_m.fill(0.086);
            this->pistons.crank_throw_length_m.fill(0.043);
            this->pistons.connecting_rod_length_m.fill(0.145);
            this->pistons.connecting_rod_mass_kg.fill(0.45);
            this->pistons.head_mass_density_kg_per_m3.fill(2700.0);
            this->pistons.head_compression_height_m.fill(0.030);
            this->pistons.head_clearance_height_m.fill(0.007);
            this->inlet_cam.ramp_theta_r.fill(g_pi_r * 0.85);
            this->outlet_cam.ramp_theta_r.fill(g_pi_r * 0.5);
            this->inlet_cam.vtec_engage_r_per_s  = { 0.0, 200.0, 400.0, 600.0, 800.0 };
            this->outlet_cam.vtec_engage_r_per_s = { 0.0, 200.0, 400.0, 600.0, 800.0 };
            this->inlet_cam.vtec_open_boost  = { 1.0, 1.2, 1.4, 1.6, 1.8 };
            this->outlet_cam.vtec_open_boost = { 1.0, 1.2, 1.4, 1.6, 1.8 };
            this->inlet_cam.vtec_ramp_boost  = { 1.00, 1.10, 1.20, 1.30, 1.40 };
            this->outlet_cam.vtec_ramp_boost = { 1.00, 1.05, 1.10, 1.15, 1.20 };
            double theta0_r = 0.0;
            for(size_t i = 0; i < get_width(); i++)
            {
                this->pistons.theta0_r[i] = theta0_r;
                this->inlet_cam.engage_theta_r[i]  = theta0_r + g_otto_intake_cycle_r - 1.0;
                this->sparkplugs.engage_theta_r[i] = theta0_r + g_otto_combustion_cycle_r - 0.8;
                this->outlet_cam.engage_theta_r[i] = theta0_r + g_otto_exhaust_cycle_r + 0.7;
                theta0_r += g_otto_cycle_r / get_width();
            }
            for(auto& flow : this->flows)
            {
                flow.chamber_nozzle_open_ratio.fill(1.0);
                flow.chamber_nozzle_flow_area_m2 = {
                    0.00250, /* Source   -> Intake    */
                    0.00120, /* Intake   -> Throttle  */
                    0.00135, /* Throttle -> Runner    */
                    0.00090, /* Runner   -> Piston    */
                    0.00120, /* Piston   -> Chamber0  */
                    0.00220, /* Chamber0 -> Chamber1  */
                    0.00320, /* Chamber1 -> Chamber2  */
                    0.00320, /* Chamber2 -> Sink      */
                };
                flow.chamber_volume_m3 = {
                    g_resevoir_volume_m3, /* Source   */
                    0.0030,               /* Intake   */
                    0.0008,               /* Throttle */
                    0.0003,               /* Runner   */
                    0.0000,               /* Piston   */
                    0.0002,               /* Chamber0 */
                    0.0003,               /* Chamber1 */
                    0.0003,               /* Chamber2 */
                    g_resevoir_volume_m3  /* Sink     */
                };
            }
            this->throttle.table = {
                0.001,
                0.050,
                0.100,
                1.000,
            };
            pipes[0].piston_connect_m = { 0.00, 0.12, 0.31, 0.19 };
            pipes[1].piston_connect_m = { 0.31, 0.00, 0.19, 0.12 };
            for(auto& pipe : this->pipes)
            {
                pipe.mic_position_m = pipe.length_m = 0.8;
            }
            this->dc.set_cutoff_frequency(5.0);
            this->gain.ratio = 0.00001;
        }
    };

    std::unique_ptr<engine> new_engine(const type type)
    {
        std::unique_ptr<engine> engine;
        switch(type)
        {
        default:
        case type::inline8: engine = std::make_unique<ensim::inline8>(); break;
        }
        engine->reset();
        return engine;
    }
}
