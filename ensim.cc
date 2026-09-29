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

    struct linear_table
    {
        static constexpr size_t size = 4;
        std::array<double, size> table = {};

        double at(const double ratio)
        {
            const double x = ratio * (size - 1);
            const size_t i = x;
            if(i >= size - 1)
            {
                return table[size - 1];
            }
            const double t = x - i;
            return table[i] + (table[i + 1] - table[i]) * t;
        }
    };

    struct throttle
    {
        linear_table mapping;
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
            U_r[Z] = r;
            U_ru[Z] = ru;
            if(u >= a)
            {
                /*
                 * Sonic or Super Sonic exit.
                 * Wave is so fast that next cell pressure value
                 * is overrided with this cell pressure value.
                 *
                 */

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
                U_rEs[Z] = calc_specific_energy_density_from_static_pressure(r, u, Ps);
            }
        }

        fn float calc_audio_sample()
        {
            const size_t Z = L - 1;
            const size_t x = Z * mic_position_m / length_m;
            const float pressure = static_pressure_pa[x] + static_pressure_pa[x - 1];
            return 0.5f * pressure;
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
                calc_local_speed_of_sounds();
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
        X(nozzle_velocity_m_per_s)

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
                const double open_ratio = throttle.mapping.at(throttle_open_ratio);
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

        const std::atomic<size_t>& get_gear() const override
        {
            return mailbox.gear;
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
        { 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, -0.000000f, 0.000000f, -0.000000f, 0.000000f, -0.000000f, 0.000000f, -0.000000f, 0.000000f, -0.000000f, 0.000000f, -0.000000f, 0.000000f, -0.000000f, 0.000000f, -0.000000f, 0.000001f, -0.000001f, 0.000001f, -0.000001f, 0.000001f, -0.000001f, 0.000001f, -0.000001f, 0.000001f, -0.000001f, 0.000001f, -0.000002f, 0.000002f, -0.000002f, 0.000002f, -0.000002f, 0.000002f, -0.000002f, 0.000002f, -0.000003f, 0.000003f, -0.000003f, 0.000003f, -0.000003f, 0.000003f, -0.000003f, 0.000003f, -0.000004f, 0.000004f, -0.000004f, 0.000004f, -0.000004f, 0.000004f, -0.000005f, 0.000005f, -0.000005f, 0.000005f, -0.000005f, 0.000005f, -0.000006f, 0.000006f, -0.000006f, 0.000006f, -0.000006f, 0.000007f, -0.000007f, 0.000007f, -0.000007f, 0.000008f, -0.000008f, 0.000008f, -0.000008f, 0.000008f, -0.000009f, 0.000009f, -0.000009f, 0.000009f, -0.000010f, 0.000010f, -0.000010f, 0.000010f, -0.000011f, 0.000011f, -0.000011f, 0.000011f, -0.000012f, 0.000012f, -0.000012f, 0.000013f, -0.000013f, 0.000013f, -0.000014f, 0.000014f, -0.000014f, 0.000015f, -0.000015f, 0.000015f, -0.000016f, 0.000016f, -0.000016f, 0.000017f, -0.000017f, 0.000018f, -0.000018f, 0.000018f, -0.000019f, 0.000019f, -0.000020f, 0.000020f, -0.000021f, 0.000021f, -0.000022f, 0.000022f, -0.000023f, 0.000023f, -0.000024f, 0.000024f, -0.000025f, 0.000025f, -0.000026f, 0.000026f, -0.000027f, 0.000028f, -0.000028f, 0.000029f, -0.000030f, 0.000030f, -0.000031f, 0.000032f, -0.000033f, 0.000033f, -0.000034f, 0.000035f, -0.000036f, 0.000037f, -0.000038f, 0.000039f, -0.000040f, 0.000041f, -0.000042f, 0.000043f, -0.000044f, 0.000045f, -0.000047f, 0.000048f, -0.000049f, 0.000051f, -0.000052f, 0.000054f, -0.000056f, 0.000057f, -0.000059f, 0.000061f, -0.000063f, 0.000065f, -0.000067f, 0.000070f, -0.000072f, 0.000075f, -0.000077f, 0.000080f, -0.000083f, 0.000087f, -0.000090f, 0.000094f, -0.000098f, 0.000103f, -0.000107f, 0.000112f, -0.000118f, 0.000124f, -0.000130f, 0.000138f, -0.000145f, 0.000154f, -0.000163f, 0.000174f, -0.000186f, 0.000199f, -0.000214f, 0.000230f, -0.000250f, 0.000272f, -0.000298f, 0.000328f, -0.000364f, 0.000407f, -0.000461f, 0.000528f, -0.000614f, 0.000729f, -0.000889f, 0.001126f, -0.001511f, 0.002236f, -0.004105f, 0.019367f, 0.007532f, -0.003686f, 0.001935f, -0.001989f, 0.001104f, -0.001498f, 0.000788f, -0.001282f, 0.000640f, -0.001177f, 0.000574f, -0.001138f, 0.000566f, -0.001151f, 0.000613f, -0.001229f, 0.000741f, -0.001418f, 0.001040f, -0.001899f, 0.001936f, -0.004047f, 0.018160f, 0.005027f, -0.002813f, 0.000855f, -0.001498f, 0.000206f, -0.001103f, -0.000048f, -0.000914f, -0.000179f, -0.000803f, -0.000256f, -0.000731f, -0.000305f, -0.000679f, -0.000337f, -0.000640f, -0.000358f, -0.000609f, -0.000372f, -0.000583f, -0.000381f, -0.000561f, -0.000386f, -0.000541f, -0.000389f, -0.000524f, -0.000390f, -0.000508f, -0.000388f, -0.000493f, -0.000386f, -0.000479f, -0.000382f, -0.000466f, -0.000378f, -0.000453f, -0.000372f, -0.000441f, -0.000367f, -0.000429f, -0.000360f, -0.000417f, -0.000354f, -0.000406f, -0.000347f, -0.000395f, -0.000339f, -0.000384f, -0.000332f, -0.000373f, -0.000324f, -0.000362f, -0.000316f, -0.000352f, -0.000308f, -0.000342f, -0.000300f, -0.000332f, -0.000292f, -0.000322f, -0.000284f, -0.000312f, -0.000276f, -0.000302f, -0.000267f, -0.000293f, -0.000259f, -0.000283f, -0.000251f, -0.000274f, -0.000243f, -0.000265f, -0.000235f, -0.000256f, -0.000227f, -0.000247f, -0.000219f, -0.000238f, -0.000212f, -0.000229f, -0.000204f, -0.000221f, -0.000196f, -0.000212f, -0.000189f, -0.000204f, -0.000182f, -0.000196f, -0.000175f, -0.000188f, -0.000168f, -0.000180f, -0.000161f, -0.000173f, -0.000154f, -0.000165f, -0.000148f, -0.000158f, -0.000141f, -0.000150f, -0.000135f, -0.000143f, -0.000129f, -0.000136f, -0.000123f, -0.000129f, -0.000118f, -0.000123f, -0.000112f, -0.000116f, -0.000107f, -0.000110f, -0.000102f, -0.000104f, -0.000096f, -0.000097f, -0.000092f, -0.000092f, -0.000087f, -0.000086f, -0.000082f, -0.000080f, -0.000078f, -0.000074f, -0.000074f, -0.000069f, -0.000070f, -0.000063f, -0.000066f, -0.000058f, -0.000063f, -0.000053f, -0.000059f, -0.000048f, -0.000056f, -0.000043f, -0.000053f, -0.000038f, -0.000050f, -0.000034f, -0.000047f, -0.000029f, -0.000045f, -0.000025f, -0.000042f, -0.000020f, -0.000040f, -0.000016f, -0.000038f, -0.000012f, -0.000036f, -0.000008f, -0.000034f, -0.000004f, -0.000032f, 0.000000f, -0.000031f, 0.000004f, -0.000030f, 0.000008f, -0.000029f, 0.000011f, -0.000028f, 0.000015f, -0.000027f, 0.000019f, -0.000026f, 0.000022f, -0.000026f, 0.000026f, -0.000025f, 0.000029f, -0.000025f, 0.000033f, -0.000025f, 0.000036f, -0.000025f, 0.000039f, -0.000026f, 0.000043f, -0.000026f, 0.000046f, -0.000027f, 0.000049f, -0.000027f, 0.000053f, -0.000028f, 0.000056f, -0.000029f, 0.000059f, -0.000031f, 0.000062f, -0.000032f, 0.000066f, -0.000034f, 0.000069f, -0.000035f, 0.000073f, -0.000037f, 0.000076f, -0.000040f, 0.000080f, -0.000042f, 0.000084f, -0.000045f, 0.000088f, -0.000048f, 0.000092f, -0.000051f, 0.000096f, -0.000055f, 0.000101f, -0.000059f, 0.000106f, -0.000063f, 0.000111f, -0.000068f, 0.000116f, -0.000073f, 0.000122f, -0.000079f, 0.000129f, -0.000085f, 0.000135f, -0.000092f, 0.000143f, -0.000100f, 0.000152f, -0.000109f, 0.000161f, -0.000119f, 0.000171f, -0.000130f, 0.000184f, -0.000143f, 0.000197f, -0.000158f, 0.000214f, -0.000176f, 0.000233f, -0.000197f, 0.000256f, -0.000223f, 0.000284f, -0.000256f, 0.000320f, -0.000297f, 0.000368f, -0.000353f, 0.000433f, -0.000432f, 0.000529f, -0.000552f, 0.000683f, -0.000760f, 0.000977f, -0.001209f, 0.001752f, -0.002913f, 0.009819f, 0.007063f, -0.002777f, 0.001438f, -0.001296f, 0.000730f, -0.000877f, 0.000452f, -0.000677f, 0.000302f, -0.000560f, 0.000209f, -0.000482f, 0.000145f, -0.000425f, 0.000098f, -0.000382f, 0.000061f, -0.000347f, 0.000032f, -0.000318f, 0.000008f, -0.000293f, -0.000012f, -0.000271f, -0.000031f, -0.000250f, -0.000047f, -0.000230f, -0.000064f, -0.000210f, -0.000080f, -0.000189f, -0.000097f, -0.000166f, -0.000118f, -0.000139f, -0.000143f, -0.000105f, -0.000180f, -0.000054f, -0.000241f, 0.000037f, -0.000376f, 0.000288f, -0.000995f, 0.011524f, 0.000941f, -0.000895f, 0.000155f, -0.000610f, 0.000008f, -0.000513f, -0.000054f, -0.000462f, -0.000088f, -0.000428f, -0.000109f, -0.000404f, -0.000123f, -0.000385f, -0.000133f, -0.000369f, -0.000140f, -0.000355f, -0.000145f, -0.000342f, -0.000148f, -0.000330f, -0.000151f, -0.000320f, -0.000152f, -0.000310f, -0.000153f, -0.000300f, -0.000153f, -0.000291f, -0.000152f, -0.000282f, -0.000152f, -0.000274f, -0.000150f, -0.000266f, -0.000149f, -0.000258f, -0.000147f, -0.000250f, -0.000145f, -0.000242f, -0.000143f, -0.000235f, -0.000141f, -0.000227f, -0.000138f, -0.000220f, -0.000136f, -0.000213f, -0.000133f, -0.000206f, -0.000130f, -0.000199f, -0.000127f, -0.000193f, -0.000124f, -0.000186f, -0.000121f, -0.000180f, -0.000118f, -0.000173f, -0.000115f, -0.000167f, -0.000112f, -0.000161f, -0.000109f, -0.000155f, -0.000105f, -0.000150f, -0.000102f, -0.000144f, -0.000099f, -0.000138f, -0.000096f, -0.000133f, -0.000093f, -0.000127f, -0.000090f, -0.000122f, -0.000087f, -0.000117f, -0.000084f, -0.000112f, -0.000081f, -0.000107f, -0.000078f, -0.000103f, -0.000075f, -0.000098f, -0.000072f, -0.000093f, -0.000069f, -0.000089f, -0.000066f, -0.000085f, -0.000063f, -0.000081f, -0.000060f, -0.000077f, -0.000058f, -0.000073f, -0.000055f, -0.000069f, -0.000052f, -0.000065f, -0.000050f, -0.000062f, -0.000047f, -0.000058f, -0.000045f, -0.000055f, -0.000043f, -0.000051f, -0.000040f, -0.000048f, -0.000038f, -0.000045f, -0.000036f, -0.000042f, -0.000034f, -0.000039f, -0.000031f, -0.000037f, -0.000029f, -0.000034f, -0.000027f, -0.000032f, -0.000025f, -0.000029f, -0.000024f, -0.000027f, -0.000022f, -0.000025f, -0.000020f, -0.000022f, -0.000018f, -0.000020f, -0.000016f, -0.000019f, -0.000015f, -0.000017f, -0.000013f, -0.000015f, -0.000012f, -0.000013f, -0.000010f, -0.000012f, -0.000008f, -0.000010f, -0.000007f, -0.000009f, -0.000006f, -0.000008f, -0.000004f, -0.000006f, -0.000003f, -0.000005f, -0.000002f, -0.000004f, -0.000000f, -0.000003f, 0.000001f, -0.000002f, 0.000002f, -0.000001f, 0.000003f, -0.000000f, 0.000004f, 0.000000f, 0.000005f, 0.000001f, 0.000006f, 0.000002f, 0.000007f, 0.000002f, 0.000008f, 0.000003f, 0.000009f, 0.000003f, 0.000010f, 0.000003f, 0.000011f, 0.000004f, 0.000012f, 0.000004f, 0.000012f, 0.000004f, 0.000013f, 0.000004f, 0.000014f, 0.000005f, 0.000015f, 0.000005f, 0.000016f, 0.000005f, 0.000016f, 0.000005f, 0.000017f, 0.000005f, 0.000018f, 0.000004f, 0.000018f, 0.000004f, 0.000019f, 0.000004f, 0.000019f, 0.000004f, 0.000020f, 0.000004f, 0.000021f, 0.000003f, 0.000021f, 0.000003f, 0.000022f, 0.000003f, 0.000022f, 0.000002f, 0.000023f, 0.000002f, 0.000024f, 0.000001f, 0.000024f, 0.000001f, 0.000025f, 0.000000f, 0.000025f, -0.000001f, 0.000026f, -0.000001f, 0.000026f, -0.000002f, 0.000027f, -0.000003f, 0.000028f, -0.000003f, 0.000028f, -0.000004f, 0.000029f, -0.000005f, 0.000030f, -0.000006f, 0.000030f, -0.000007f, 0.000031f, -0.000008f, 0.000032f, -0.000009f, 0.000033f, -0.000010f, 0.000033f, -0.000011f, 0.000034f, -0.000013f, 0.000035f, -0.000014f, 0.000036f, -0.000015f, 0.000037f, -0.000017f, 0.000039f, -0.000018f, 0.000040f, -0.000020f, 0.000041f, -0.000022f, 0.000043f, -0.000024f, 0.000044f, -0.000026f, 0.000046f, -0.000028f, 0.000048f, -0.000031f, 0.000050f, -0.000033f, 0.000053f, -0.000036f, 0.000055f, -0.000039f, 0.000058f, -0.000043f, 0.000061f, -0.000047f, 0.000065f, -0.000051f, 0.000069f, -0.000056f, 0.000074f, -0.000061f, 0.000079f, -0.000068f, 0.000085f, -0.000075f, 0.000093f, -0.000083f, 0.000102f, -0.000093f, 0.000112f, -0.000105f, 0.000125f, -0.000120f, 0.000141f, -0.000139f, 0.000162f, -0.000163f, 0.000190f, -0.000196f, 0.000229f, -0.000244f, 0.000286f, -0.000317f, 0.000381f, -0.000444f, 0.000561f, -0.000723f, 0.001048f, -0.001812f, 0.006717f, 0.003987f, -0.001697f, 0.000885f, -0.000835f, 0.000467f, -0.000587f, 0.000303f, -0.000470f, 0.000216f, -0.000401f, 0.000163f, -0.000357f, 0.000127f, -0.000326f, 0.000102f, -0.000304f, 0.000085f, -0.000287f, 0.000072f, -0.000274f, 0.000063f, -0.000264f, 0.000056f, -0.000257f, 0.000052f, -0.000252f, 0.000050f, -0.000248f, 0.000051f, -0.000247f, 0.000053f, -0.000247f, 0.000058f, -0.000251f, 0.000065f, -0.000257f, 0.000076f, -0.000267f, 0.000092f, -0.000283f, 0.000115f, -0.000308f, 0.000150f, -0.000349f, 0.000208f, -0.000422f, 0.000315f, -0.000579f, 0.000581f, -0.001105f, 0.002266f, 0.006815f, -0.001558f, 0.000576f, -0.000696f, 0.000215f, -0.000492f, 0.000088f, -0.000400f, 0.000025f, -0.000348f, -0.000011f, -0.000315f, -0.000034f, -0.000291f, -0.000049f, -0.000274f, -0.000059f, -0.000260f, -0.000065f, -0.000249f, -0.000070f, -0.000240f, -0.000072f, -0.000233f, -0.000073f, -0.000227f, -0.000072f, -0.000221f, -0.000071f, -0.000217f, -0.000069f, -0.000213f, -0.000067f, -0.000209f, -0.000064f, -0.000207f, -0.000060f, -0.000204f, -0.000055f, -0.000203f, -0.000050f, -0.000202f, -0.000045f, -0.000201f, -0.000039f, -0.000201f, -0.000032f, -0.000202f, -0.000024f, -0.000204f, -0.000015f, -0.000207f, -0.000005f, -0.000211f, 0.000006f, -0.000217f, 0.000020f, -0.000226f, 0.000037f, -0.000239f, 0.000059f, -0.000259f, 0.000090f, -0.000291f, 0.000140f, -0.000353f, 0.000242f, -0.000530f, 0.000720f, 0.007939f, -0.000604f, -0.000069f, -0.000180f, -0.000274f, -0.000009f, -0.000460f, 0.000280f, -0.001101f, 0.007199f, 0.001097f, -0.000920f, 0.000156f, -0.000594f, -0.000006f, -0.000490f, -0.000067f, -0.000444f, -0.000088f, -0.000432f, -0.000070f, -0.000474f, 0.000051f, -0.000820f, 0.007730f, 0.000405f, -0.000767f, -0.000017f, -0.000656f, -0.000014f, -0.000728f, 0.000187f, -0.001260f, 0.007055f, 0.000681f, -0.000921f, -0.000114f, -0.000641f, -0.000249f, -0.000548f, -0.000302f, -0.000499f, -0.000329f, -0.000465f, -0.000344f, -0.000440f, -0.000353f, -0.000419f, -0.000359f, -0.000400f, -0.000363f, -0.000383f, -0.000364f, -0.000367f, -0.000365f, -0.000352f, -0.000365f, -0.000337f, -0.000364f, -0.000322f, -0.000364f, -0.000308f, -0.000362f, -0.000294f, -0.000361f, -0.000279f, -0.000360f, -0.000265f, -0.000358f, -0.000251f, -0.000357f, -0.000236f, -0.000356f, -0.000221f, -0.000355f, -0.000206f, -0.000355f, -0.000191f, -0.000355f, -0.000175f, -0.000356f, -0.000158f, -0.000357f, -0.000141f, -0.000359f, -0.000123f, -0.000363f, -0.000104f, -0.000368f, -0.000083f, -0.000375f, -0.000060f, -0.000384f, -0.000034f, -0.000398f, -0.000003f, -0.000417f, 0.000035f, -0.000446f, 0.000085f, -0.000492f, 0.000158f, -0.000572f, 0.000286f, -0.000746f, 0.000594f, -0.001351f, 0.002648f, 0.006339f, -0.001774f, 0.000627f, -0.000906f, 0.000324f, -0.000823f, 0.000465f, -0.001411f, 0.006347f, 0.001517f, -0.001204f, 0.000317f, -0.000850f, 0.000179f, -0.000782f, 0.000171f, -0.000792f, 0.000220f, -0.000853f, 0.000329f, -0.001001f, 0.000594f, -0.001567f, 0.007623f, -0.000070f, -0.000021f, -0.001112f, 0.001132f, 0.005811f, -0.000359f, -0.001548f, 0.003929f, 0.003962f, -0.002081f, 0.000403f, -0.001144f, -0.000058f, -0.000838f, -0.000266f, -0.000658f, -0.000403f, -0.000516f, -0.000523f, -0.000374f, -0.000661f, -0.000190f, -0.000879f, 0.000155f, -0.001481f, 0.001998f, 0.005411f, -0.001605f, -0.000258f, 0.000131f, 0.006104f, -0.001880f, 0.000210f, -0.001191f, -0.000101f, -0.000979f, -0.000236f, -0.000854f, -0.000320f, -0.000759f, -0.000385f, -0.000676f, -0.000445f, -0.000591f, -0.000513f, -0.000488f, -0.000612f, -0.000332f, -0.000809f, 0.000031f, -0.001580f, 0.004854f, 0.001916f, -0.001667f, 0.000192f, -0.001122f, -0.000067f, -0.000942f, -0.000171f, -0.000844f, -0.000226f, -0.000777f, -0.000259f, -0.000726f, -0.000281f, -0.000683f, -0.000295f, -0.000646f, -0.000305f, -0.000613f, -0.000312f, -0.000583f, -0.000317f, -0.000554f, -0.000320f, -0.000527f, -0.000322f, -0.000501f, -0.000323f, -0.000476f, -0.000324f, -0.000451f, -0.000325f, -0.000426f, -0.000325f, -0.000402f, -0.000326f, -0.000377f, -0.000328f, -0.000351f, -0.000332f, -0.000323f, -0.000338f, -0.000293f, -0.000348f, -0.000257f, -0.000367f, -0.000208f, -0.000408f, -0.000119f, -0.000537f, 0.000267f, 0.006200f, -0.000793f, -0.000129f, -0.000471f, -0.000248f, -0.000384f, -0.000292f, -0.000330f, -0.000319f, -0.000286f, -0.000343f, -0.000241f, -0.000370f, -0.000189f, -0.000410f, -0.000118f, -0.000480f, 0.000005f, -0.000647f, 0.000332f, -0.001378f, 0.004273f, 0.003505f, -0.002549f, 0.005969f, 0.000736f, -0.000730f, -0.000220f, -0.000489f, -0.000292f, -0.000452f, -0.000287f, -0.000452f, -0.000260f, -0.000466f, -0.000220f, -0.000495f, -0.000165f, -0.000542f, -0.000083f, -0.000627f, 0.000058f, -0.000818f, 0.000413f, -0.001575f, 0.004453f, 0.002377f, -0.001435f, 0.000284f, -0.000827f, -0.000004f, -0.000634f, -0.000123f, -0.000528f, -0.000193f, -0.000450f, -0.000247f, -0.000379f, -0.000302f, -0.000299f, -0.000377f, -0.000180f, -0.000527f, 0.000100f, -0.001144f, 0.005019f, 0.001326f, -0.001166f, 0.000203f, -0.000788f, 0.000023f, -0.000662f, -0.000049f, -0.000594f, -0.000085f, -0.000548f, -0.000106f, -0.000515f, -0.000118f, -0.000488f, -0.000124f, -0.000466f, -0.000127f, -0.000446f, -0.000128f, -0.000429f, -0.000127f, -0.000414f, -0.000124f, -0.000400f, -0.000120f, -0.000388f, -0.000114f, -0.000377f, -0.000108f, -0.000367f, -0.000101f, -0.000359f, -0.000092f, -0.000351f, -0.000083f, -0.000346f, -0.000072f, -0.000342f, -0.000059f, -0.000341f, -0.000043f, -0.000344f, -0.000022f, -0.000353f, 0.000006f, -0.000372f, 0.000048f, -0.000412f, 0.000123f, -0.000509f, 0.000305f, -0.000857f, 0.001450f, 0.004672f, -0.001171f, 0.000326f, -0.000587f, 0.000100f, -0.000463f, 0.000046f, -0.000432f, 0.000056f, -0.000457f, 0.000130f, -0.000574f, 0.000365f, -0.001064f, 0.002150f, 0.004531f, -0.001525f, 0.000552f, -0.000831f, 0.000269f, -0.000672f, 0.000187f, -0.000615f, 0.000167f, -0.000602f, 0.000184f, -0.000623f, 0.000236f, -0.000690f, 0.000354f, -0.000863f, 0.000677f, -0.001554f, 0.003822f, 0.002877f, -0.001261f, 0.000298f, -0.000602f, -0.000008f, -0.000409f, -0.000129f, -0.000309f, -0.000199f, -0.000238f, -0.000251f, -0.000174f, -0.000306f, -0.000097f, -0.000394f, 0.000059f, -0.000731f, 0.004925f, 0.000309f, -0.000513f, -0.000137f, -0.000334f, -0.000233f, -0.000255f, -0.000287f, -0.000194f, -0.000335f, -0.000131f, -0.000391f, -0.000051f, -0.000478f, 0.000079f, -0.000653f, 0.000387f, -0.001281f, 0.003020f, 0.003057f, -0.001410f, 0.000403f, -0.000767f, 0.000108f, -0.000582f, -0.000006f, -0.000487f, -0.000067f, -0.000426f, -0.000107f, -0.000380f, -0.000136f, -0.000341f, -0.000161f, -0.000304f, -0.000184f, -0.000268f, -0.000209f, -0.000227f, -0.000241f, -0.000177f, -0.000287f, -0.000106f, -0.000365f, 0.000016f, -0.000531f, 0.000312f, -0.001127f, 0.002621f, 0.003785f, -0.001578f, 0.000580f, -0.000888f, 0.000276f, -0.000704f, 0.000168f, -0.000620f, 0.000120f, -0.000576f, 0.000100f, -0.000554f, 0.000102f, -0.000557f, 0.000136f, -0.000611f, 0.000279f, -0.001000f, 0.003337f, 0.006927f, 0.000668f, -0.000730f, -0.000059f, -0.000579f, -0.000109f, -0.000510f, 0.005245f, -0.000779f, -0.000021f, -0.000772f, 0.000049f, -0.000887f, 0.000287f, -0.001433f, 0.004239f, 0.000939f, -0.000915f, -0.000173f, -0.000552f, -0.000350f, -0.000427f, -0.000427f, -0.000350f, -0.000480f, -0.000285f, -0.000531f, -0.000212f, -0.000600f, -0.000104f, -0.000732f, 0.000131f, -0.001174f, 0.001662f, 0.003678f, -0.001432f, 0.000160f, -0.000806f, -0.000104f, -0.000641f, -0.000198f, -0.000559f, -0.000244f, -0.000508f, -0.000269f, -0.000470f, -0.000284f, -0.000441f, -0.000293f, -0.000417f, -0.000298f, -0.000395f, -0.000300f, -0.000376f, -0.000300f, -0.000359f, -0.000298f, -0.000343f, -0.000296f, -0.000328f, -0.000292f, -0.000315f, -0.000288f, -0.000302f, -0.000282f, -0.000291f, -0.000275f, -0.000281f, -0.000267f, -0.000273f, -0.000255f, -0.000270f, -0.000239f, -0.000274f, -0.000211f, -0.000296f, -0.000154f, -0.000373f, 0.000025f, -0.000821f, 0.003720f, 0.006942f, -0.000359f, -0.000157f, -0.000534f, -0.000250f, 0.000162f, 0.006980f, 0.002141f, -0.000823f, -0.000395f, -0.000402f, -0.000560f, -0.000287f, -0.000635f, -0.000198f, -0.000714f, -0.000081f, -0.000855f, 0.000166f, -0.001314f, 0.001752f, 0.003645f, -0.001427f, 0.000047f, -0.000775f, -0.000235f, -0.000593f, -0.000344f, -0.000492f, -0.000410f, -0.000416f, -0.000463f, -0.000344f, -0.000519f, -0.000260f, -0.000599f, -0.000133f, -0.000755f, 0.000148f, -0.001302f, 0.002328f, 0.002833f, -0.001474f, 0.000203f, -0.000866f, -0.000067f, -0.000690f, -0.000169f, -0.000598f, -0.000222f, -0.000538f, -0.000253f, -0.000492f, -0.000275f, -0.000454f, -0.000290f, -0.000420f, -0.000302f, -0.000389f, -0.000313f, -0.000359f, -0.000323f, -0.000329f, -0.000334f, -0.000297f, -0.000347f, -0.000262f, -0.000366f, -0.000220f, -0.000394f, -0.000165f, -0.000441f, -0.000079f, -0.000540f, 0.000102f, -0.000861f, 0.001112f, 0.004489f, -0.001176f, 0.000130f, -0.000586f, -0.000152f, -0.000289f, 0.004716f, -0.000713f, -0.000122f, -0.000524f, -0.000204f, -0.000451f, -0.000248f, -0.000392f, -0.000292f, -0.000319f, -0.000371f, -0.000158f, -0.000725f, 0.001880f, 0.006845f, 0.000245f, -0.000462f, -0.000456f, -0.000342f, -0.000475f, -0.000327f, -0.000464f, -0.000322f, -0.000449f, -0.000317f, -0.000438f, -0.000307f, -0.000435f, -0.000283f, -0.000459f, -0.000194f, -0.000680f, 0.001122f, 0.006396f, 0.000726f, -0.000678f, -0.000367f, -0.000455f, -0.000448f, -0.000391f, -0.000481f, -0.000341f, -0.000514f, -0.000281f, -0.000568f, -0.000183f, -0.000690f, 0.000055f, -0.001206f, 0.003495f, 0.001124f, -0.001095f, 0.000006f, -0.000737f, -0.000159f, -0.000619f, -0.000224f, -0.000553f, -0.000258f, -0.000505f, -0.000283f, -0.000460f, -0.000314f, -0.000389f, -0.000443f, 0.003909f, -0.000036f, -0.000720f, 0.004148f, -0.000446f, -0.000406f, -0.000498f, -0.000383f, -0.000488f, -0.000378f, -0.000470f, -0.000376f, -0.000450f, -0.000377f, -0.000426f, -0.000384f, -0.000393f, -0.000406f, -0.000334f, -0.000475f, -0.000170f, -0.000867f, 0.003775f, 0.000914f, -0.001243f, 0.000401f, -0.001476f, 0.002396f, 0.002043f, -0.001237f, 0.000041f, -0.000783f, -0.000128f, -0.000693f, -0.000142f, -0.000699f, -0.000075f, -0.000801f, 0.000151f, -0.001261f, 0.001929f, 0.002716f, -0.001117f, 0.004232f, -0.001244f, 0.000089f, -0.000925f, -0.000042f, -0.000838f, -0.000071f, -0.000814f, -0.000051f, -0.000836f, 0.000018f, -0.000885f, 0.004884f, 0.002587f, 0.002127f, 0.001835f, -0.001353f, -0.000102f, -0.000845f, -0.000328f, -0.000692f, -0.000412f, -0.000610f, -0.000459f, -0.000544f, -0.000505f, -0.000459f, -0.000614f, -0.000127f, 0.003799f, -0.001006f, -0.000230f, -0.000761f, -0.000297f, -0.000721f, -0.000284f, -0.000736f, -0.000214f, -0.000829f, -0.000003f, -0.001268f, 0.001905f, 0.002224f, -0.001441f, 0.000098f, -0.000960f, -0.000087f, -0.000852f, -0.000116f, -0.000838f, -0.000069f, -0.000908f, 0.000110f, -0.001266f, 0.001482f, 0.002585f, -0.001220f, -0.000090f, -0.000704f, -0.000304f, -0.000558f, -0.000385f, -0.000473f, -0.000436f, -0.000398f, -0.000494f, -0.000299f, -0.000604f, -0.000079f, -0.001064f, 0.002867f, 0.001134f, -0.001160f, 0.000087f, -0.000848f, -0.000029f, -0.000785f, -0.000014f, -0.000835f, 0.000155f, -0.001247f, 0.004009f, 0.004438f, -0.001081f, -0.000140f, -0.000691f, -0.000293f, -0.000588f, -0.000338f, -0.000538f, -0.000353f, -0.000508f, -0.000355f, -0.000486f, -0.000350f, -0.000470f, -0.000341f, -0.000457f, -0.000329f, -0.000447f, -0.000315f, -0.000440f, -0.000298f, -0.000436f, -0.000278f, -0.000435f, -0.000254f, -0.000433f, -0.000258f, -0.000168f, 0.005608f, 0.001487f, -0.000868f, -0.000181f, -0.000544f, -0.000300f, -0.000462f, -0.000334f, -0.000420f, -0.000346f, -0.000391f, -0.000350f, -0.000368f, -0.000350f, -0.000348f, -0.000347f, -0.000330f, -0.000344f, -0.000312f, -0.000340f, -0.000295f, -0.000336f, -0.000278f, -0.000333f, -0.000261f, -0.000330f, -0.000244f, -0.000328f, -0.000225f, -0.000328f, -0.000204f, -0.000331f, -0.000179f, -0.000341f, -0.000143f, -0.000368f, -0.000078f, -0.000462f, 0.000217f, 0.003698f, -0.000574f, -0.000133f, -0.000339f, -0.000216f, -0.000277f, -0.000243f, -0.000241f, -0.000257f, -0.000214f, -0.000265f, -0.000190f, -0.000272f, -0.000165f, -0.000280f, -0.000140f, -0.000291f, -0.000111f, -0.000307f, -0.000075f, -0.000332f, -0.000025f, -0.000379f, 0.000057f, -0.000482f, 0.000248f, -0.000849f, 0.001582f, 0.002755f, -0.001001f, 0.000248f, -0.000514f, 0.000040f, -0.000383f, -0.000035f, -0.000317f, -0.000073f, -0.000275f, -0.000095f, -0.000244f, -0.000109f, -0.000220f, -0.000118f, -0.000198f, -0.000125f, -0.000179f, -0.000131f, -0.000162f, -0.000136f, -0.000145f, -0.000141f, -0.000128f, -0.000146f, -0.000110f, -0.000153f, -0.000091f, -0.000161f, -0.000070f, -0.000173f, -0.000044f, -0.000191f, -0.000011f, -0.000222f, 0.000045f, 0.004156f, 0.000041f, -0.000428f, 0.000263f, -0.000907f, 0.003510f, 0.001186f, -0.000787f, 0.000188f, -0.000471f, 0.000036f, -0.000370f, -0.000024f, -0.000318f, -0.000055f, -0.000284f, -0.000074f, -0.000260f, -0.000086f, -0.000240f, -0.000094f, -0.000224f, -0.000099f, -0.000209f, -0.000103f, -0.000197f, -0.000106f, -0.000185f, -0.000108f, -0.000174f, -0.000109f, -0.000163f, -0.000111f, -0.000152f, -0.000112f, -0.000142f, -0.000113f, -0.000132f, -0.000115f, -0.000121f, -0.000116f, -0.000111f, -0.000119f, -0.000099f, -0.000122f, -0.000087f, -0.000127f, -0.000073f, -0.000133f, -0.000057f, -0.000142f, -0.000038f, -0.000156f, -0.000014f, -0.000176f, 0.000018f, -0.000208f, 0.000068f, -0.000265f, 0.000156f, -0.000388f, 0.000377f, -0.000861f, 0.003030f, 0.001135f, -0.000275f, 0.003772f, -0.000824f, 0.000276f, -0.000484f, 0.000109f, -0.000361f, 0.000005f, -0.000225f, -0.000253f, 0.003825f, 0.000482f, -0.000621f, 0.000044f, 0.001195f, 0.009092f, 0.000655f, -0.000557f, -0.000304f, -0.000295f, -0.000416f, -0.000222f, -0.000459f, -0.000179f, -0.000491f, -0.000136f, -0.000530f, -0.000076f, -0.000602f, 0.000054f, -0.000869f, 0.003623f, -0.000231f, -0.000128f, -0.000932f, 0.002388f, 0.000954f, -0.000931f, -0.000040f, -0.000617f, -0.000192f, -0.000510f, -0.000260f, -0.000443f, -0.000310f, -0.000372f, -0.000399f, -0.000133f, 0.002989f, -0.000749f, -0.000151f, -0.000568f, -0.000210f, -0.000529f, -0.000217f, -0.000518f, -0.000206f, -0.000517f, -0.000210f, -0.000272f, 0.005396f, 0.003460f, 0.002616f, 0.000644f, -0.001011f, -0.000186f, -0.000743f, -0.000306f, -0.000658f, -0.000348f, -0.000613f, -0.000368f, -0.000578f, -0.000387f, -0.000530f, -0.000470f, 0.002849f, -0.000206f, -0.000715f, -0.000287f, -0.000686f, -0.000272f, -0.000699f, -0.000221f, -0.000758f, -0.000090f, -0.000997f, 0.000717f, 0.003024f, -0.001126f, -0.000154f, -0.000702f, -0.000311f, -0.000606f, -0.000353f, -0.000449f, 0.004701f, 0.001014f, -0.001044f, -0.000158f, -0.000865f, -0.000108f, -0.001139f, 0.002991f, 0.000213f, -0.000873f, -0.000283f, -0.000840f, 0.002643f, -0.000209f, -0.000757f, -0.000370f, -0.000715f, -0.000348f, -0.000737f, -0.000279f, -0.000816f, -0.000111f, -0.001137f, 0.001411f, 0.003885f, 0.001606f, 0.001719f, -0.001248f, -0.000228f, -0.000851f, -0.000393f, -0.000733f, -0.000452f, -0.000667f, -0.000480f, -0.000621f, -0.000496f, -0.000583f, -0.000505f, -0.000550f, -0.000511f, -0.000519f, -0.000515f, -0.000489f, -0.000520f, -0.000458f, -0.000525f, -0.000425f, -0.000534f, -0.000388f, -0.000548f, -0.000343f, -0.000573f, -0.000284f, -0.000619f, -0.000191f, -0.000725f, 0.000019f, -0.001124f, 0.001811f, 0.001357f, -0.000998f, -0.000134f, -0.000555f, -0.000394f, -0.000264f, -0.000865f, 0.002469f, 0.001046f, -0.001152f, 0.000111f, -0.000860f, -0.000003f, -0.000792f, 0.000004f, -0.000828f, 0.000150f, -0.001172f, 0.001956f, 0.001257f, -0.000960f, -0.000067f, -0.000625f, -0.000193f, -0.000547f, -0.000210f, -0.000533f, -0.000177f, -0.000573f, -0.000063f, -0.000782f, 0.000676f, 0.002616f, -0.000975f, -0.000026f, -0.000595f, -0.000172f, -0.000496f, -0.000217f, -0.000446f, -0.000235f, -0.000412f, -0.000243f, -0.000386f, -0.000245f, -0.000365f, -0.000244f, -0.000346f, -0.000242f, -0.000328f, -0.000239f, -0.000311f, -0.000235f, -0.000296f, -0.000230f, -0.000280f, -0.000226f, -0.000265f, -0.000222f, -0.000250f, -0.000218f, -0.000235f, -0.000215f, -0.000218f, -0.000214f, -0.000200f, -0.000216f, -0.000178f, -0.000224f, -0.000148f, -0.000245f, -0.000095f, -0.000308f, 0.000043f, -0.000598f, 0.001285f, 0.002838f, -0.001650f, 0.003630f, 0.000431f, -0.000350f, -0.000312f, 0.000112f, 0.003195f, -0.000662f, -0.000038f, -0.000434f, -0.000113f, -0.000387f, -0.000123f, -0.000375f, -0.000110f, -0.000384f, -0.000076f, -0.000417f, -0.000006f, -0.000511f, 0.000182f, -0.000910f, 0.002212f, 0.001393f, -0.000873f, 0.000143f, -0.000528f, -0.000012f, -0.000426f, -0.000066f, -0.000376f, -0.000091f, -0.000345f, -0.000102f, -0.000325f, -0.000106f, -0.000310f, -0.000104f, -0.000300f, -0.000099f, -0.000294f, -0.000089f, -0.000292f, -0.000076f, -0.000296f, -0.000056f, -0.000307f, -0.000026f, -0.000331f, 0.000021f, -0.000381f, 0.000105f, -0.000491f, 0.000301f, -0.000848f, 0.001375f, 0.002484f, 0.000438f, 0.003181f, -0.000281f, 0.005689f, -0.000984f, 0.003143f, 0.000218f, -0.000384f, -0.000436f, -0.000176f, -0.000539f, -0.000080f, -0.000634f, 0.000082f, -0.000978f, 0.002915f, 0.000357f, -0.000596f, -0.000249f, -0.000381f, -0.000357f, -0.000294f, -0.000416f, -0.000225f, -0.000475f, -0.000137f, -0.000581f, 0.000056f, -0.000977f, 0.002221f, 0.001040f, -0.000887f, -0.000011f, -0.000578f, -0.000150f, -0.000484f, -0.000197f, -0.000442f, -0.000208f, -0.000431f, -0.000185f, -0.000470f, -0.000067f, -0.000760f, 0.001908f, 0.000897f, -0.000941f, 0.000166f, -0.000779f, 0.000167f, -0.000859f, 0.000383f, -0.001335f, 0.002681f, 0.001051f, -0.000456f, -0.000847f, 0.003744f, 0.006000f, -0.001868f, 0.000475f, -0.001148f, 0.000186f, -0.000992f, 0.000133f, -0.000996f, 0.000252f, -0.001335f, 0.003032f, -0.000147f, 0.002512f, 0.000162f, -0.000758f, -0.000303f, -0.000567f, -0.000395f, -0.000491f, -0.000437f, -0.000435f, -0.000480f, -0.000350f, -0.000637f, 0.002572f, -0.000135f, 0.002442f, -0.000379f, -0.000579f, -0.000405f, -0.000547f, -0.000408f, -0.000530f, -0.000396f, -0.000537f, -0.000333f, -0.000696f, 0.000632f, 0.002932f, 0.000904f, 0.003112f, 0.001961f, -0.001309f, -0.000100f, -0.000892f, -0.000281f, -0.000755f, -0.000383f, -0.000596f, -0.000685f, 0.001543f, 0.000959f, -0.001678f, 0.002686f, -0.000400f, -0.000400f, -0.000825f, -0.000234f, -0.000919f, -0.000115f, -0.001044f, 0.000109f, -0.001453f, 0.001778f, 0.001246f, -0.000911f, -0.000514f, -0.000344f, -0.001039f, 0.001922f, 0.000450f, -0.000909f, -0.000370f, -0.000602f, -0.000544f, -0.000422f, -0.000774f, 0.002486f, -0.000285f, -0.000607f, -0.000535f, -0.000471f, -0.000609f, -0.000386f, -0.000662f, -0.000314f, -0.000688f, -0.000362f, 0.000781f, 0.002576f, 0.000496f, 0.001328f, 0.000728f, 0.000985f, 0.000888f, -0.001228f, -0.000191f, -0.000782f, -0.000517f, -0.000049f, 0.003906f, 0.000358f, -0.001286f, 0.001668f, 0.000625f, -0.001095f, -0.000357f, -0.000748f, -0.000568f, -0.000506f, -0.000899f, 0.000552f, 0.001765f, -0.001341f, -0.000130f, -0.000941f, -0.000301f, -0.000811f, -0.000370f, -0.000731f, -0.000411f, -0.000668f, -0.000442f, -0.000610f, -0.000473f, -0.000547f, -0.000516f, -0.000460f, -0.000603f, -0.000280f, -0.000947f, 0.001655f, 0.000836f, -0.001098f, -0.000072f, -0.000830f, -0.000166f, -0.000783f, -0.000105f, -0.001017f, 0.004995f, 0.000240f, 0.002169f, -0.000570f, -0.000564f, -0.000476f, -0.000594f, -0.000416f, -0.000643f, -0.000306f, -0.000807f, 0.000226f, 0.001857f, -0.000896f, -0.000328f, -0.000566f, -0.000485f, -0.000414f, -0.000618f, -0.000209f, -0.000923f, 0.000653f, 0.001548f, 0.000449f, 0.001845f, -0.001453f, 0.000170f, -0.001013f, -0.000027f, -0.000873f, -0.000098f, -0.000801f, -0.000127f, -0.000757f, -0.000134f, -0.000732f, -0.000126f, -0.000723f, -0.000099f, -0.000736f, -0.000041f, -0.000799f, 0.000109f, -0.001086f, 0.001294f, 0.001404f, -0.000832f, -0.000208f, -0.000454f, -0.000395f, -0.000226f, 0.002522f, -0.000573f, -0.000313f, -0.000429f, -0.000361f, -0.000380f, -0.000375f, -0.000348f, -0.000380f, -0.000322f, -0.000382f, -0.000297f, -0.000384f, -0.000271f, -0.000389f, -0.000239f, -0.000404f, -0.000192f, -0.000448f, -0.000083f, -0.000676f, 0.002292f, 0.000218f, -0.000545f, -0.000153f, -0.000408f, -0.000207f, -0.000359f, -0.000223f, -0.000330f, -0.000225f, -0.000311f, -0.000222f, -0.000296f, -0.000216f, -0.000284f, -0.000206f, -0.000276f, -0.000194f, -0.000271f, -0.000178f, -0.000271f, -0.000155f, -0.000281f, -0.000119f, -0.000314f, -0.000039f, -0.000449f, 0.000417f, 0.002161f, -0.000650f, 0.000044f, -0.000389f, -0.000052f, -0.000325f, -0.000078f, -0.000295f, -0.000085f, -0.000277f, -0.000084f, -0.000265f, -0.000078f, -0.000258f, -0.000069f, -0.000254f, -0.000056f, -0.000253f, -0.000040f, -0.000258f, -0.000019f, -0.000269f, 0.000011f, -0.000291f, 0.000055f, -0.000335f, 0.000132f, -0.000435f, 0.000314f, -0.000778f, 0.001500f, 0.002462f, -0.000660f, -0.000110f, 0.000677f, 0.002923f, -0.001246f, 0.000661f, -0.001027f, 0.001120f, 0.002538f, -0.000346f, -0.000410f, 0.000479f, 0.002767f, -0.000706f, -0.000007f, -0.000361f, -0.000168f, -0.000219f, -0.000373f, 0.001958f, 0.002222f, 0.002605f, 0.000128f, -0.000538f, -0.000138f, -0.000450f, -0.000172f, -0.000420f, -0.000181f, -0.000404f, -0.000181f, -0.000395f, -0.000176f, -0.000389f, -0.000168f, -0.000386f, -0.000156f, -0.000387f, -0.000141f, -0.000391f, -0.000122f, -0.000401f, -0.000094f, -0.000422f, -0.000050f, -0.000472f, 0.000052f, -0.000685f, 0.002233f, 0.000060f, -0.000320f, -0.000281f, -0.000174f, -0.000369f, -0.000083f, -0.000457f, 0.000050f, -0.000674f, 0.002572f, -0.000740f, 0.001833f, 0.001562f, -0.000998f, 0.000185f, -0.000733f, 0.000646f, 0.003617f, -0.000045f, 0.003147f, 0.005693f, 0.000794f, -0.001054f, 0.000350f, 0.002145f, -0.000796f, -0.000412f, -0.000425f, -0.000749f, 0.000636f, 0.001705f, -0.001502f, 0.000819f, 0.001818f, -0.000919f, -0.000333f, -0.000653f, -0.000402f, -0.000624f, -0.000384f, -0.000644f, -0.000325f, -0.000711f, -0.000194f, -0.000922f, 0.000395f, 0.001277f, 0.000650f, 0.001220f, -0.001233f, -0.000039f, -0.000882f, -0.000196f, -0.000772f, -0.000250f, -0.000719f, -0.000266f, -0.000695f, -0.000256f, -0.000699f, -0.000209f, -0.000762f, -0.000042f, -0.001149f, 0.003222f, 0.004624f, 0.000179f, 0.001214f, 0.000715f, -0.000955f, -0.000350f, -0.000705f, -0.000439f, -0.000641f, -0.000459f, -0.000609f, -0.000461f, -0.000587f, -0.000455f, -0.000571f, -0.000446f, -0.000558f, -0.000433f, -0.000548f, -0.000417f, -0.000543f, -0.000394f, -0.000547f, -0.000359f, -0.000570f, -0.000288f, -0.000671f, 0.000041f, 0.001745f, -0.000523f, 0.002073f, -0.000836f, -0.000224f, -0.000652f, -0.000283f, -0.000604f, -0.000289f, -0.000587f, -0.000270f, -0.000600f, 0.002102f, -0.000596f, -0.000242f, -0.000671f, -0.000099f, -0.000944f, 0.001329f, 0.000621f, -0.000795f, -0.000219f, -0.000548f, -0.000324f, -0.000464f, -0.000364f, -0.000411f, -0.000388f, -0.000366f, -0.000412f, -0.000313f, -0.000455f, -0.000218f, -0.000600f, 0.000263f, 0.002338f, -0.000875f, -0.000059f, -0.000612f, -0.000122f, -0.000608f, -0.000017f, -0.000894f, 0.001662f, 0.001244f, -0.001088f, 0.000242f, -0.001049f, 0.002095f, -0.000049f, -0.000335f, -0.000523f, -0.000071f, 0.001269f, 0.002084f, 0.000308f, -0.000527f, -0.000548f, 0.000393f, 0.001636f, -0.000888f, -0.000238f, -0.000145f, 0.004578f, -0.001167f, 0.000043f, -0.000834f, -0.000086f, -0.000745f, -0.000124f, -0.000704f, -0.000131f, -0.000686f, -0.000118f, -0.000688f, -0.000083f, -0.000716f, -0.000012f, -0.000801f, 0.000159f, -0.001112f, 0.001256f, 0.001669f, -0.000626f, -0.000657f, 0.001576f, 0.001081f, -0.001041f, 0.000027f, -0.000717f, -0.000120f, -0.000614f, -0.000166f, -0.000571f, -0.000167f, -0.000574f, -0.000099f, -0.000736f, 0.000715f, 0.001437f, 0.001992f, -0.000103f, -0.000852f, 0.000731f, 0.001919f, -0.001018f, -0.000046f, -0.000658f, -0.000221f, 0.003781f, 0.000619f, -0.001363f, 0.001819f, 0.002846f, -0.000755f, -0.000410f, -0.000489f, -0.000518f, -0.000403f, -0.000566f, -0.000340f, -0.000610f, -0.000265f, -0.000682f, -0.000132f, -0.000885f, 0.000461f, 0.002048f, -0.000615f, -0.000746f, 0.001146f, 0.001389f, -0.001201f, 0.000022f, -0.000845f, -0.000124f, -0.000756f, -0.000131f, -0.000828f, 0.002435f, -0.000384f, -0.000410f, -0.000552f, -0.000280f, -0.000783f, 0.001522f, 0.001970f, 0.000473f, -0.000505f, 0.001957f, 0.000626f, 0.001586f, -0.001293f, 0.000042f, -0.001053f, 0.000105f, 0.000970f, 0.003474f, 0.001164f, -0.001137f, -0.000238f, -0.000839f, -0.000316f, -0.000831f, -0.000208f, -0.001105f, 0.002050f, -0.000364f, 0.000591f, 0.001430f, -0.001437f, 0.000055f, -0.001121f, -0.000038f, -0.001098f, 0.000051f, -0.001309f, 0.000930f, 0.000960f, -0.000678f, -0.000702f, -0.000195f, -0.001184f, 0.001504f, 0.000788f, -0.000878f, -0.000509f, -0.000356f, -0.001023f, 0.000903f, 0.003102f, 0.000159f, 0.000918f, -0.001586f, 0.000565f, 0.001353f, -0.000662f, -0.000650f, -0.000496f, -0.000688f, 0.001378f, 0.000564f, 0.001699f, -0.001316f, -0.000119f, -0.001003f, -0.000223f, -0.000934f, -0.000231f, -0.000921f, -0.000193f, -0.000961f, -0.000068f, 0.003429f, 0.001924f, -0.000075f, -0.000673f, -0.000697f, -0.000417f, -0.000866f, -0.000188f, -0.001232f, 0.001430f, 0.000526f, -0.000801f, -0.000753f, 0.000915f, 0.000453f, -0.001012f, -0.000379f, -0.000421f, 0.001760f, -0.001526f, 0.001303f, 0.000222f, -0.000815f, -0.000466f, -0.000673f, -0.000443f, -0.000826f, 0.001131f, 0.001036f, 0.000978f, -0.000755f, -0.000554f, -0.000547f, -0.000681f, -0.000036f, 0.003302f, -0.000005f, -0.000741f, -0.000531f, -0.000628f, -0.000551f, -0.000597f, -0.000544f, -0.000579f, -0.000528f, -0.000568f, -0.000506f, -0.000565f, -0.000472f, -0.000582f, -0.000397f, -0.000714f, 0.002153f, 0.002041f, -0.001076f, 0.000612f, 0.001721f, 0.002586f, -0.000822f, -0.000480f, -0.000543f, 0.003766f, -0.000905f, -0.000396f, -0.000791f, -0.000418f, -0.000763f, -0.000396f, -0.000813f, 0.001259f, 0.001686f, -0.000003f, -0.000987f, -0.000235f, -0.000966f, -0.000127f, -0.001199f, 0.000879f, 0.001100f, 0.001204f, -0.000481f, -0.000630f, -0.000603f, -0.000524f, -0.000657f, -0.000450f, -0.000697f, -0.000379f, -0.000743f, -0.000294f, -0.000813f, -0.000167f, -0.000959f, 0.000110f, -0.001493f, 0.003155f, 0.001575f, -0.001041f, -0.000480f, 0.000121f, 0.001246f, -0.000998f, -0.000399f, -0.000011f, 0.003059f, -0.001046f, -0.000192f, -0.000819f, -0.000309f, -0.000704f, -0.000386f, -0.000583f, -0.000562f, 0.001555f, -0.000054f, -0.000804f, -0.000272f, -0.000699f, -0.000311f, -0.000648f, -0.000325f, -0.000611f, -0.000329f, -0.000582f, -0.000329f, -0.000556f, -0.000325f, -0.000533f, -0.000319f, -0.000512f, -0.000312f, -0.000493f, -0.000304f, -0.000475f, -0.000294f, -0.000459f, -0.000283f, -0.000445f, -0.000270f, -0.000434f, -0.000253f, -0.000428f, -0.000228f, -0.000437f, -0.000173f, -0.000539f, 0.001574f, -0.000136f, -0.000384f, -0.000280f, -0.000320f, -0.000294f, -0.000293f, -0.000292f, -0.000276f, -0.000285f, -0.000262f, -0.000275f, -0.000250f, -0.000265f, -0.000239f, -0.000254f, -0.000229f, -0.000243f, -0.000220f, -0.000231f, -0.000211f, -0.000220f, -0.000203f, -0.000209f, -0.000195f, -0.000197f, -0.000187f, -0.000186f, -0.000180f, -0.000175f, -0.000173f, -0.000164f, -0.000167f, -0.000152f, -0.000161f, -0.000141f, -0.000156f, -0.000130f, -0.000151f, -0.000118f, -0.000148f, -0.000105f, -0.000147f, 0.001871f, -0.000166f, -0.000107f, -0.000167f, -0.000091f, -0.000170f, -0.000073f, -0.000177f, -0.000051f, -0.000190f, -0.000017f, -0.000222f, 0.000049f, 0.001548f, 0.002093f, 0.000146f, -0.000245f, -0.000108f, -0.000141f, -0.000160f, -0.000101f, -0.000129f, 0.002028f, -0.000640f, 0.000696f, 0.001662f, -0.000402f, -0.000170f, -0.000013f, -0.000525f, 0.000785f, 0.001372f, 0.000980f, 0.001194f, -0.000842f, 0.000183f, -0.000411f, -0.000262f, 0.000933f, 0.005385f, 0.000359f, 0.001721f, 0.000272f, -0.000198f, -0.000550f, -0.000024f, -0.000617f, 0.000028f, -0.000652f, 0.000079f, -0.000702f, 0.000164f, -0.000824f, 0.000452f, 0.000971f, 0.004291f, 0.001642f, -0.000362f, -0.000601f, -0.000015f, -0.000883f, 0.000372f, -0.001608f, 0.004437f, 0.001873f, 0.000309f, 0.000598f, -0.001002f, -0.000097f, -0.000651f, -0.000373f, -0.000182f, 0.001819f, -0.000984f, 0.001654f, -0.000631f, -0.000333f, -0.000584f, 0.000958f, 0.002118f, 0.002432f, 0.000188f, -0.000537f, -0.000603f, -0.000457f, -0.000598f, -0.000459f, -0.000572f, -0.000468f, -0.000542f, -0.000480f, -0.000508f, -0.000498f, -0.000462f, -0.000541f, -0.000337f, 0.001531f, -0.000587f, -0.000475f, -0.000454f, -0.000551f, -0.000349f, -0.000669f, -0.000129f, -0.001097f, 0.001677f, 0.002121f, -0.001444f, 0.000158f, -0.001000f, -0.000006f, -0.000923f, 0.000020f, -0.001046f, 0.000631f, 0.000689f, -0.000521f, -0.000698f, 0.001114f, 0.000094f, -0.000695f, -0.000325f, -0.000543f, -0.000385f, -0.000490f, -0.000401f, -0.000460f, -0.000403f, -0.000439f, -0.000398f, -0.000424f, -0.000389f, -0.000411f, -0.000378f, -0.000401f, -0.000365f, -0.000393f, -0.000349f, -0.000388f, -0.000331f, -0.000386f, -0.000309f, -0.000389f, -0.000280f, -0.000404f, -0.000232f, -0.000449f, -0.000124f, -0.000654f, 0.000847f, 0.000800f, -0.000759f, 0.000038f, 0.001385f, -0.000282f, -0.000392f, -0.000246f, -0.000379f, -0.000244f, -0.000362f, -0.000241f, -0.000345f, -0.000237f, -0.000331f, -0.000231f, -0.000318f, -0.000224f, -0.000307f, -0.000214f, -0.000301f, -0.000196f, -0.000311f, -0.000139f, -0.000455f, 0.001923f, 0.000413f, -0.000887f, 0.001539f, 0.000824f, -0.000680f, -0.000005f, -0.000460f, -0.000091f, -0.000400f, -0.000112f, -0.000374f, -0.000114f, -0.000361f, -0.000105f, -0.000359f, -0.000085f, -0.000371f, -0.000045f, -0.000414f, 0.000055f, -0.000621f, 0.001320f, 0.000433f, -0.000442f, -0.000092f, -0.000277f, -0.000165f, -0.000218f, -0.000195f, -0.000180f, -0.000216f, -0.000144f, -0.000241f, -0.000096f, -0.000291f, 0.000007f, -0.000498f, 0.001381f, 0.000371f, -0.000337f, 0.001967f, -0.000456f, -0.000012f, -0.000350f, -0.000046f, -0.000331f, -0.000028f, -0.000374f, 0.000150f, 0.001925f, -0.000303f, -0.000271f, 0.001638f, 0.002476f, -0.000721f, 0.000641f, 0.003310f, -0.000365f, -0.000086f, 0.001845f, -0.000643f, 0.000008f, -0.000731f, 0.000785f, 0.000947f, -0.000653f, -0.000164f, -0.000332f, -0.000464f, 0.000660f, 0.001069f, -0.000949f, 0.000181f, -0.000780f, 0.000190f, -0.000920f, 0.003169f, -0.000159f, 0.000560f, 0.001118f, -0.000759f, -0.000161f, -0.000448f, -0.000318f, -0.000321f, 0.001306f, 0.000219f, 0.001032f, -0.000672f, -0.000236f, -0.000412f, -0.000409f, -0.000153f, 0.001103f, 0.000569f, 0.001086f, -0.001022f, 0.000071f, -0.000774f, -0.000018f, -0.000735f, 0.000031f, 0.001500f, -0.000143f, -0.000719f, 0.000039f, -0.000984f, 0.000930f, 0.002530f, -0.000854f, -0.000299f, -0.000176f, 0.001728f, -0.000967f, 0.000319f, 0.002639f, -0.000358f, -0.000491f, -0.000417f, -0.000479f, -0.000379f, -0.000594f, 0.000961f, 0.001009f, 0.000827f, 0.003932f, -0.000244f, -0.000590f, -0.000494f, -0.000505f, -0.000525f, -0.000466f, -0.000539f, -0.000434f, -0.000550f, -0.000400f, -0.000566f, -0.000360f, -0.000593f, -0.000298f, -0.000657f, -0.000165f, -0.000901f, 0.000895f, 0.000816f, -0.000852f, -0.000227f, -0.000582f, -0.000343f, -0.000498f, -0.000284f, 0.002953f, 0.000598f, -0.001112f, 0.001038f, 0.001080f, -0.001246f, 0.000447f, 0.000586f, 0.000770f, 0.000469f, -0.000379f, 0.001180f, 0.002053f, 0.000227f, -0.001236f, 0.000612f, 0.000367f, -0.000684f, -0.000519f, -0.000515f, -0.000563f, -0.000479f, -0.000562f, -0.000464f, -0.000550f, -0.000453f, -0.000535f, -0.000445f, -0.000519f, -0.000437f, -0.000502f, -0.000430f, -0.000484f, -0.000425f, -0.000464f, -0.000421f, -0.000443f, -0.000420f, -0.000416f, -0.000428f, -0.000377f, -0.000457f, -0.000290f, -0.000630f, 0.001240f, 0.000686f, 0.001805f, -0.001156f, 0.003717f, 0.002464f, -0.001230f, 0.000087f, 0.000419f, 0.001122f, -0.000137f, -0.000427f, -0.000756f, 0.000192f, 0.001380f, -0.000842f, 0.001567f, -0.000660f, -0.000420f, -0.000567f, -0.000456f, -0.000522f, -0.000468f, -0.000490f, -0.000472f, -0.000463f, -0.000474f, -0.000437f, -0.000477f, -0.000408f, -0.000483f, -0.000374f, -0.000499f, -0.000323f, -0.000544f, -0.000214f, -0.000746f, 0.000968f, 0.000283f, -0.000665f, -0.000258f, -0.000466f, -0.000375f, -0.000311f, -0.000630f, 0.002380f, 0.000516f, -0.000866f, -0.000026f, 0.001571f, -0.000188f, -0.000638f, -0.000166f, -0.000641f, -0.000124f, -0.000687f, -0.000003f, -0.000927f, 0.001135f, 0.000808f, -0.000737f, -0.000196f, -0.000484f, -0.000290f, -0.000423f, -0.000289f, -0.000458f, -0.000040f, 0.002096f, -0.000900f, 0.000232f, 0.001594f, -0.000431f, -0.000493f, 0.001300f, -0.000319f, -0.000350f, -0.000447f, -0.000255f, -0.000521f, -0.000128f, -0.000738f, 0.000918f, 0.000437f, -0.000616f, -0.000229f, -0.000401f, -0.000346f, -0.000264f, -0.000547f, 0.000868f, 0.000371f, -0.000712f, 0.000006f, -0.000777f, 0.002194f, 0.001727f, 0.001628f, 0.002070f, -0.000413f, -0.000317f, -0.000524f, -0.000274f, -0.000540f, -0.000232f, -0.000573f, -0.000163f, -0.000649f, -0.000012f, -0.000923f, 0.001126f, 0.000525f, 0.000849f, 0.001933f, 0.001818f, -0.000016f, 0.001399f, -0.000012f, 0.000083f, 0.001583f, 0.000113f, -0.000281f, -0.000894f, 0.000598f, 0.001110f, -0.001116f, 0.000479f, 0.000560f, -0.000868f, -0.000177f, -0.000901f, 0.000587f, 0.000456f, -0.000426f, 0.001657f, 0.000280f, -0.000300f, -0.000841f, 0.000054f, 0.003049f, 0.001527f, -0.000852f, 0.001398f, 0.001195f, -0.000474f, -0.001020f, 0.002176f, 0.000481f, -0.000298f, 0.001049f, -0.001306f, -0.000103f, -0.001140f, 0.000527f, 0.001883f, -0.000489f, -0.000751f, -0.000505f, -0.000780f, 0.001635f, -0.000615f, -0.000608f, -0.000659f, -0.000570f, -0.000668f, -0.000517f, -0.000766f, 0.000322f, 0.003676f, -0.000289f, -0.000741f, -0.000569f, -0.000687f, -0.000524f, -0.000815f, 0.000327f, 0.003231f, -0.001297f, 0.000902f, 0.002374f, -0.000782f, -0.000777f, 0.000791f, -0.000355f, -0.000520f, 0.001871f, 0.000248f, -0.000946f, -0.000522f, -0.000829f, -0.000507f, -0.000861f, -0.000394f, -0.001044f, 0.000259f, 0.000797f, -0.000839f, -0.000780f, 0.000374f, 0.000445f, -0.001017f, -0.000498f, -0.000343f, 0.001522f, -0.001369f, 0.001474f, -0.000860f, -0.000379f, -0.000885f, -0.000355f, -0.000884f, 0.001843f, -0.000953f, -0.000159f, 0.001085f, -0.000402f, -0.001126f, 0.001483f, 0.001377f, -0.001250f, 0.000035f, 0.002678f, -0.000245f, -0.000832f, -0.000534f, -0.000734f, -0.000559f, -0.000692f, -0.000571f, -0.000549f, 0.002664f, -0.000082f, -0.000820f, -0.000502f, -0.000746f, -0.000446f, -0.000900f, 0.000531f, 0.002191f, -0.001062f, -0.000288f, -0.000911f, -0.000292f, -0.000925f, 0.001504f, -0.001513f, 0.002236f, 0.001853f, -0.001663f, 0.001502f, -0.000130f, -0.001058f, 0.000967f, 0.001453f, -0.001480f, 0.000275f, 0.000910f, -0.000723f, -0.000661f, 0.001341f, -0.001115f, 0.000920f, 0.000178f, -0.001045f, -0.000149f, 0.001317f, 0.001255f, -0.000915f, 0.001470f, 0.000792f, -0.000829f, -0.000810f, 0.000372f, 0.002135f, -0.001454f, 0.000794f, 0.000000f, -0.000910f, -0.000522f, -0.000929f, 0.000766f, 0.000645f, 0.002749f, -0.000743f, 0.000227f, 0.002825f, -0.001629f, 0.002064f, 0.000369f, -0.001207f, -0.000154f, 0.001083f, -0.001192f, -0.000458f, -0.001018f, -0.000471f, -0.001023f, -0.000380f, -0.001174f, 0.000195f, 0.000600f, 0.001240f, 0.003141f, 0.000997f, -0.000895f, -0.000744f, -0.000668f, 0.003659f, 0.001358f, -0.000935f, 0.001134f, -0.000707f, -0.000938f, -0.000697f, -0.000928f, -0.000690f, -0.000705f, 0.002588f, -0.001092f, 0.002326f, 0.000073f, 0.000022f, 0.001945f, 0.000460f, -0.000644f, -0.000530f, 0.002740f, -0.001011f, -0.000718f, -0.001030f, -0.000726f, -0.000998f, -0.000680f, -0.001106f, 0.000140f, 0.002366f, -0.001478f, 0.000589f, -0.000622f, -0.000668f, -0.001346f, 0.002252f, 0.001928f, -0.001426f, -0.000516f, -0.001088f, -0.000669f, -0.000926f, -0.000844f, 0.000284f, 0.000134f, 0.000173f, 0.000458f, -0.000803f, -0.000930f, -0.000311f, 0.000723f, -0.000263f, 0.000369f, -0.001796f, 0.001201f, -0.000412f, 0.000409f, 0.000017f, 0.000347f, -0.000469f, -0.001055f, -0.000521f, -0.001142f, 0.001046f, -0.000453f, -0.000852f, -0.000768f, -0.000642f, -0.001012f, 0.000282f, 0.000539f, -0.000928f, 0.001138f, -0.001415f, -0.000145f, -0.001279f, -0.000047f, -0.001593f, 0.001891f, 0.002301f, -0.001340f, -0.000385f, -0.000926f, -0.000406f, 0.001015f, -0.000971f, -0.000405f, -0.001085f, 0.000682f, 0.000315f, -0.001050f, -0.000357f, 0.000690f, -0.000561f, -0.000854f, 0.000184f, 0.000423f, -0.001290f, 0.000596f, 0.000365f, -0.000852f, -0.000573f, -0.000380f, 0.002429f, 0.003377f, 0.001131f, 0.000423f, -0.000562f, -0.000854f, 0.000778f, -0.000617f, -0.000707f, -0.000648f, -0.000687f, -0.000614f, -0.000695f, 0.000792f, -0.000820f, 0.002098f, 0.000372f, -0.001084f, -0.000292f, -0.001082f, 0.000533f, 0.000259f, -0.000948f, -0.000419f, -0.000788f, -0.000447f, 0.000421f, 0.000811f, 0.000025f, -0.000842f, -0.000463f, -0.000689f, -0.000527f, -0.000438f, 0.002827f, 0.001262f, -0.000988f, -0.000337f, -0.000940f, 0.000388f, 0.001786f, -0.001299f, 0.001025f, -0.000095f, -0.000767f, -0.000525f, -0.000619f, -0.000598f, -0.000454f, 0.001150f, -0.000783f, -0.000435f, -0.000688f, -0.000255f, 0.003288f, 0.000103f, -0.000933f, -0.000297f, -0.000851f, -0.000248f, -0.000983f, 0.000469f, 0.000213f, -0.000706f, -0.000522f, -0.000470f, -0.000679f, 0.001726f, 0.000774f, -0.000830f, -0.000388f, -0.000618f, -0.000463f, -0.000545f, -0.000491f, -0.000492f, -0.000512f, -0.000433f, -0.000573f, 0.000942f, -0.000415f, -0.000466f, -0.000496f, -0.000390f, -0.000542f, -0.000301f, -0.000584f, 0.003011f, 0.000688f, -0.000847f, -0.000240f, -0.000596f, -0.000348f, -0.000468f, 0.001205f, 0.000789f, -0.000463f, -0.000455f, -0.000425f, -0.000488f, 0.000839f, 0.000883f, 0.000527f, -0.000701f, -0.000328f, -0.000473f, -0.000464f, -0.000229f, 0.000826f, -0.000524f, -0.000531f, 0.000434f, 0.000591f, -0.000871f, -0.000070f, 0.000666f, -0.000323f, -0.000135f, 0.001874f, 0.000315f, 0.001164f, -0.000891f, -0.000011f, -0.000941f, 0.000367f, 0.000214f, 0.000931f, 0.001724f, -0.001167f, 0.000010f, -0.000829f, -0.000043f, -0.000972f, 0.001332f, 0.001753f, -0.000956f, -0.000115f, -0.000755f, 0.001288f, 0.000608f, 0.000647f, 0.000033f, 0.000794f, 0.000153f, -0.000785f, -0.000263f, -0.000639f, -0.000301f, 0.000696f, -0.000347f, -0.000605f, -0.000270f, -0.000686f, 0.003508f, 0.000362f, -0.000799f, -0.000309f, -0.000563f, -0.000512f, 0.000436f, 0.000346f, -0.001013f, 0.000093f, -0.001200f, 0.002259f, 0.000893f, 0.001504f, 0.005690f, -0.000694f, 0.000189f, 0.000365f, -0.000730f, -0.000531f, -0.000581f, -0.000567f, -0.000546f, -0.000567f, -0.000523f, -0.000563f, -0.000500f, -0.000564f, -0.000467f, -0.000582f, -0.000400f, -0.000682f, -0.000008f, 0.000897f, -0.000838f, -0.000283f, -0.000642f, -0.000368f, 0.000814f, -0.000281f, -0.000664f, -0.000282f, -0.000680f, -0.000193f, -0.000853f, 0.000671f, 0.001918f, -0.000722f, -0.000393f, -0.000250f, 0.000931f, -0.001175f, 0.001473f, 0.001659f, -0.000828f, -0.000307f, -0.000592f, -0.000404f, -0.000485f, -0.000540f, 0.000129f, 0.002361f, -0.000667f, -0.000254f, 0.000599f, -0.000308f, -0.000670f, -0.000241f, -0.000787f, 0.000701f, 0.002385f, 0.000264f, -0.000801f, -0.000231f, -0.000747f, -0.000167f, -0.000898f, 0.000583f, 0.000437f, -0.000729f, -0.000280f, 0.000718f, 0.000659f, -0.000270f, -0.000635f, -0.000310f, -0.000617f, -0.000290f, -0.000619f, -0.000260f, -0.000627f, -0.000223f, -0.000643f, -0.000175f, -0.000680f, -0.000061f, 0.002796f, 0.000122f, -0.001037f, 0.001124f, -0.000293f, -0.000099f, 0.000750f, -0.000378f, -0.000670f, -0.000046f, -0.001053f, 0.001491f, 0.001357f, 0.002190f, 0.002164f, -0.000425f, 0.002601f, 0.000504f, -0.000835f, 0.000082f, 0.000327f, -0.000748f, -0.000450f, -0.000543f, -0.000550f, -0.000446f, -0.000611f, -0.000363f, -0.000663f, 0.001071f, -0.000768f, 0.002523f, 0.001591f, -0.000216f, -0.001205f, 0.002945f, 0.000989f, -0.000421f, 0.000725f, -0.000865f, -0.000609f, -0.000038f, 0.000553f, -0.000763f, -0.000704f, 0.000409f, 0.001830f, -0.000871f, 0.002066f, 0.001434f, -0.001477f, 0.002244f, 0.001323f, 0.000152f, 0.000157f, -0.000111f, 0.000647f, -0.000108f, -0.000808f, -0.000237f, 0.000838f, 0.002486f, 0.001189f, -0.000788f, -0.000582f, -0.000860f, -0.000602f, -0.000835f, 0.000128f, 0.001955f, 0.001667f, -0.001178f, -0.000518f, -0.000840f, -0.000757f, -0.000365f, 0.000202f, 0.000067f, 0.000232f, -0.001259f, -0.000289f, -0.001096f, -0.000234f, 0.000216f, -0.000695f, -0.000720f, -0.000687f, -0.000501f, 0.001962f, -0.001061f, -0.000435f, -0.000864f, -0.000507f, -0.000780f, -0.000557f, -0.000678f, -0.000688f, -0.000188f, 0.000685f, -0.001102f, -0.000200f, 0.001626f, -0.000343f, -0.000839f, -0.000434f, -0.000750f, -0.000540f, 0.000627f, -0.000163f, -0.000917f, -0.000293f, -0.000859f, -0.000286f, -0.000854f, -0.000234f, -0.000919f, 0.000042f, 0.000346f, -0.000500f, -0.000684f, -0.000285f, -0.000867f, 0.000822f, -0.000424f, -0.000385f, -0.000715f, -0.000165f, -0.000972f, 0.000523f, 0.000201f, 0.002391f, 0.001505f, -0.001202f, -0.000096f, -0.000773f, -0.000385f, 0.000473f, 0.000170f, -0.000965f, -0.000033f, -0.000990f, 0.000350f, 0.000383f, 0.000584f, 0.000176f, -0.000727f, -0.000334f, -0.000549f, -0.000403f, -0.000483f, -0.000425f, -0.000441f, -0.000435f, -0.000403f, -0.000446f, -0.000358f, -0.000479f, -0.000264f, -0.000649f, 0.001676f, 0.002269f, 0.001719f, -0.000010f, -0.000784f, 0.000108f, 0.000527f, -0.000610f, -0.000360f, -0.000470f, -0.000395f, -0.000434f, -0.000395f, -0.000412f, -0.000389f, -0.000394f, -0.000381f, -0.000376f, -0.000343f, 0.001858f, 0.001457f, -0.000670f, -0.000156f, -0.000708f, 0.000912f, -0.000181f, -0.000414f, -0.000435f, -0.000265f, -0.000570f, 0.001812f, 0.001391f, 0.000435f, 0.000358f, -0.000789f, -0.000138f, -0.000614f, -0.000211f, -0.000550f, -0.000236f, -0.000512f, -0.000244f, -0.000485f, -0.000244f, -0.000466f, -0.000236f, -0.000459f, -0.000205f, -0.000526f, 0.000537f, 0.000647f, 0.000657f, -0.000569f, -0.000216f, -0.000399f, -0.000326f, 0.001730f, 0.000178f, -0.000601f, -0.000162f, -0.000413f, 0.001099f, -0.000534f, 0.001279f, -0.000712f, -0.000032f, -0.000621f, -0.000082f, 0.000321f, 0.001108f, -0.000294f, -0.000133f, -0.000648f, 0.000165f, 0.000540f, -0.000232f, -0.000564f, 0.000008f, -0.000818f, 0.000747f, 0.001704f, -0.000667f, -0.000216f, -0.000405f, 0.000841f, -0.000117f, -0.000150f, 0.001897f, -0.000487f, 0.000128f, -0.001119f, 0.001692f, 0.000647f, 0.001548f, 0.002367f, -0.001017f, -0.000087f, -0.000581f, -0.000320f, -0.000391f, -0.000472f, -0.000209f, -0.000691f, 0.002630f, -0.000681f, 0.001708f, 0.000783f, -0.000808f, 0.000785f, 0.000272f, 0.001661f, 0.000185f, 0.000275f, -0.000780f, -0.000430f, 0.000400f, 0.000961f, -0.001093f, 0.001810f, -0.001118f, 0.000542f, 0.000425f, 0.001728f, -0.000587f, 0.000993f, -0.000731f, -0.000335f, -0.000659f, -0.000365f, -0.000609f, -0.000400f, -0.000533f, -0.000520f, 0.000399f, -0.000031f, -0.000747f, -0.000250f, -0.000661f, -0.000274f, -0.000640f, 0.000787f, 0.002303f, -0.000207f, -0.000773f, 0.000024f, 0.001082f, -0.000858f, 0.000888f, 0.000303f, 0.000712f, -0.000605f, -0.000538f, 0.001474f, 0.002459f, 0.000310f, -0.000737f, -0.000469f, -0.000563f, -0.000552f, -0.000450f, -0.000724f, 0.000946f, -0.000125f, -0.000700f, -0.000422f, -0.000544f, -0.000555f, -0.000200f, 0.000624f, -0.000834f, 0.000711f, -0.000432f, -0.000631f, -0.000083f, 0.000479f, 0.001006f, 0.000717f, -0.000878f, 0.001116f, -0.000665f, 0.001027f, -0.001399f, 0.000726f, 0.000323f, -0.000577f, -0.000431f, 0.000776f, -0.001061f, 0.000019f, -0.001165f, 0.001332f, -0.000719f, 0.001178f, -0.000320f, -0.000033f, 0.000886f, -0.000930f, 0.000438f, 0.000506f, 0.001432f, 0.001501f, -0.000795f, -0.000389f, -0.000723f, -0.000367f, -0.000738f, -0.000316f, -0.000768f, -0.000273f, -0.000474f, 0.002200f, 0.000100f, -0.000428f, -0.000975f, 0.001317f, 0.000733f, 0.001052f, -0.000615f, 0.001602f, 0.000191f, -0.000738f, -0.000659f, 0.000451f, -0.000402f, -0.000446f, -0.000989f, 0.001019f, 0.001520f, 0.002289f, 0.003094f, 0.000927f, -0.001057f, 0.000572f, -0.000675f, 0.000654f, -0.000042f, -0.000690f, 0.002223f, -0.000632f, 0.001206f, 0.000242f, 0.000338f, 0.000027f, 0.000053f, -0.001155f, -0.000420f, -0.000481f, 0.000844f, -0.000631f, -0.000520f, -0.000847f, 0.000833f, 0.000597f, 0.001123f, -0.000582f, 0.001860f, -0.000639f, -0.000131f, 0.000971f, -0.000774f, -0.000586f, 0.002012f, 0.000920f, 0.001497f, -0.000569f, 0.002018f, -0.001558f, 0.000558f, -0.000423f, -0.000707f, -0.001150f, 0.000424f, -0.000078f, -0.000973f, -0.000824f, -0.000237f, -0.000053f, -0.001111f, -0.000449f, 0.000608f, 0.001107f, -0.000693f, -0.000819f, -0.000728f, -0.000814f, -0.000659f, -0.000928f, 0.000206f, -0.000130f, -0.001125f, 0.000528f, -0.000632f, 0.001987f, -0.000815f, 0.000395f, 0.000067f, -0.001226f, -0.000251f, -0.001310f, 0.000508f, -0.000071f, 0.000005f, 0.000201f, -0.000288f, 0.000254f, -0.001181f, -0.000355f, 0.000517f, -0.000522f, 0.000118f, -0.000472f, -0.000842f, -0.000590f, -0.000484f, 0.000871f, -0.001373f, 0.000496f, 0.001069f, -0.000846f, 0.000330f, -0.000292f, -0.000954f, -0.000331f, -0.001021f, 0.000009f, 0.000562f, -0.000635f, -0.000823f, 0.000020f, 0.001074f, 0.001064f, -0.000303f, -0.000728f, -0.000581f, -0.000662f, 0.000295f, 0.001312f, 0.000630f, -0.000566f, -0.000793f, -0.000336f, -0.001096f, 0.001350f, 0.001025f, -0.001245f, 0.000181f, 0.000156f, -0.000160f, 0.003938f, 0.001540f, -0.000901f, -0.000643f, -0.000559f, -0.000888f, 0.001548f, 0.000332f, 0.000177f, 0.000340f, -0.001296f, 0.000656f, -0.000195f, -0.000831f, -0.000496f, -0.000847f, 0.000647f, -0.000638f, 0.000526f, -0.000655f, -0.000454f, -0.001040f, 0.001673f, 0.000507f, -0.001167f, -0.000197f, -0.001079f, 0.001031f, 0.000419f, 0.000127f, -0.000201f, -0.001020f, -0.000079f, -0.001471f, 0.002288f, 0.001692f, -0.001282f, 0.001546f, -0.001581f, 0.001257f, 0.000463f, -0.000870f, 0.000566f, 0.000087f, -0.001174f, 0.000666f, -0.000656f, -0.000411f, -0.000951f, 0.000238f, 0.000376f, 0.000418f, 0.000260f, 0.000062f, -0.001012f, 0.000204f, -0.000047f, -0.000755f, -0.000518f, -0.000628f, -0.000549f, -0.000587f, -0.000546f, 0.000779f, -0.000548f, -0.000577f, -0.000450f, 0.001193f, 0.000390f, -0.001031f, 0.000464f, 0.000413f, -0.000867f, -0.000343f, -0.000690f, -0.000383f, -0.000680f, -0.000197f, 0.000990f, 0.000151f, 0.000353f, -0.000160f, -0.000683f, -0.000442f, -0.000076f, 0.001749f, 0.000619f, 0.000338f, -0.000535f, -0.000671f, 0.000582f, 0.000126f, -0.000903f, -0.000172f, -0.000898f, 0.000261f, -0.000003f, -0.000549f, -0.000529f, -0.000394f, -0.000618f, -0.000164f, 0.000796f, 0.000479f, -0.000742f, -0.000201f, -0.000832f, 0.001388f, 0.000671f, -0.000539f, -0.000645f, 0.000197f, 0.001794f, 0.000066f, 0.000013f, -0.000820f, -0.000197f, -0.000718f, -0.000217f, -0.000708f, -0.000157f, -0.000825f, 0.000429f, 0.000261f, 0.002204f, 0.000742f, -0.000850f, -0.000280f, -0.000478f, 0.000635f, -0.000707f, -0.000287f, -0.000548f, -0.000435f, 0.000455f, 0.001373f, 0.000212f, 0.000654f, -0.000162f, -0.000111f, 0.001602f, -0.000342f, -0.000607f, -0.000286f, -0.000732f, 0.000138f, 0.000401f, -0.000602f, -0.000414f, -0.000477f, 0.000601f, 0.002317f, 0.000685f, -0.000695f, -0.000500f, -0.000128f, 0.000331f, 0.000342f, 0.001035f, -0.000655f, -0.000374f, -0.000499f, 0.000573f, 0.001776f, -0.000075f, -0.000759f, -0.000248f, -0.000787f, 0.001546f, -0.000210f, -0.000603f, -0.000409f, -0.000547f, -0.000396f, -0.000561f, -0.000330f, -0.000652f, -0.000052f, 0.000269f, 0.001075f, 0.000716f, -0.001051f, 0.000046f, -0.001068f, 0.001306f, 0.001953f, 0.000921f, -0.000680f, 0.000090f, 0.000171f, -0.000775f, -0.000254f, -0.000717f, -0.000116f, 0.000707f, 0.000980f, 0.000846f, -0.000835f, -0.000320f, -0.000586f, -0.000423f, -0.000491f, -0.000489f, -0.000379f, -0.000653f, 0.000283f, 0.000108f, 0.001174f, 0.000373f, -0.000962f, -0.000050f, -0.000388f, 0.002040f, 0.003055f, -0.000837f, -0.000346f, 0.000427f, -0.000594f, -0.000385f, -0.000692f, -0.000059f, 0.000226f, 0.000972f, 0.000929f, -0.000214f, -0.000331f, -0.000725f, -0.000293f, -0.000108f, 0.001704f, 0.000207f, -0.000481f, -0.000976f, 0.002233f, 0.001647f, -0.001496f, 0.000203f, -0.001302f, 0.000687f, 0.000872f, 0.000904f, -0.000002f, -0.000458f, 0.001293f, -0.000542f, -0.000687f, 0.000120f, 0.000327f, -0.001010f, 0.000010f, 0.001267f, 0.000593f, -0.000450f, -0.000537f, -0.000690f, -0.000178f, 0.000321f, -0.000819f, -0.000362f, -0.000740f, 0.000679f, 0.001296f, -0.000677f, -0.000579f, 0.000356f, -0.000189f, 0.000408f, 0.000399f, 0.001375f, -0.000063f, -0.000887f, 0.001973f, -0.000617f, 0.001810f, -0.000309f, -0.000717f, -0.000546f, -0.000626f, -0.000572f, -0.000585f, -0.000578f, -0.000553f, -0.000581f, -0.000519f, -0.000589f, -0.000473f, -0.000636f, 0.000405f, -0.000513f, -0.000472f, -0.000606f, -0.000353f, -0.000770f, 0.000253f, 0.001057f, 0.000522f, 0.001113f, -0.000347f, 0.001433f, -0.000263f, 0.000456f, -0.000149f, -0.000849f, 0.000497f, -0.000445f, 0.000144f, 0.000536f, -0.001178f, 0.000413f, -0.000002f, -0.000479f, -0.000767f, 0.000442f, 0.000135f, -0.000759f, -0.000362f, -0.000617f, -0.000393f, -0.000609f, -0.000245f, 0.000886f, 0.000452f, -0.000779f, -0.000053f, 0.001049f, 0.000976f, 0.000265f, -0.000967f, -0.000088f, -0.000994f, 0.001058f, -0.000380f, 0.000214f, 0.000635f, 0.000116f, -0.000288f, 0.000447f, -0.000718f, -0.000351f, -0.000618f, -0.000360f, -0.000604f, -0.000331f, -0.000626f, -0.000227f, 0.001562f, 0.000661f, -0.000486f, -0.000397f, -0.000601f, -0.000179f, 0.001140f, -0.000506f, -0.000490f, -0.000393f, -0.000546f, -0.000305f, -0.000663f, 0.000809f, -0.000341f, -0.000391f, -0.000542f, -0.000223f, -0.000767f, 0.000893f, -0.000179f, 0.000160f, 0.000719f, -0.000620f, 0.001325f, -0.001321f, 0.002679f, 0.000737f, -0.000787f, -0.000249f, -0.000596f, -0.000292f, -0.000499f, 0.002170f, 0.000065f, -0.000499f, -0.000562f, 0.000115f, 0.002036f, 0.000005f, 0.000695f, -0.000844f, 0.000146f, 0.000490f, 0.000362f, -0.000668f, -0.000347f, -0.000553f, -0.000378f, -0.000519f, -0.000371f, -0.000528f, -0.000221f, 0.000584f, -0.000831f, 0.000641f, 0.000221f, -0.000693f, -0.000247f, -0.000556f, 0.000886f, -0.000728f, 0.000819f, -0.000155f, -0.000536f, 0.000771f, -0.000327f, -0.000570f, 0.000131f, 0.000409f, -0.000847f, 0.000505f, 0.000005f, 0.000515f, -0.000562f, -0.000277f, -0.000559f, -0.000197f, -0.000708f, 0.001374f, 0.000006f, 0.000040f, 0.000525f, -0.000945f, 0.000411f, 0.001484f, -0.000830f, -0.000111f, -0.000747f, 0.000295f, 0.001425f, 0.000183f, -0.000393f, -0.000412f, -0.000464f, -0.000231f, 0.000739f, -0.000545f, -0.000347f, -0.000421f, -0.000330f, 0.001524f, 0.000515f, 0.000997f, 0.000791f, 0.000201f, 0.001943f, 0.000364f, -0.000207f, -0.000692f, -0.000257f, -0.000650f, -0.000259f, -0.000637f, -0.000240f, -0.000639f, -0.000225f, 0.000771f, 0.000544f, -0.001559f, 0.002828f, 0.002141f, -0.000653f, 0.001253f, -0.001256f, 0.000084f, -0.000950f, -0.000071f, -0.000652f, 0.000963f, -0.001620f, 0.001915f, 0.002371f, -0.000369f, 0.001529f, 0.000766f, -0.000505f, 0.002349f, -0.000039f, 0.000786f, -0.000357f, 0.000921f, 0.000433f, 0.000129f, -0.000488f, -0.000646f, 0.000201f, 0.000133f, 0.000111f, -0.000782f, -0.000522f, -0.000554f, -0.000743f, 0.000122f, 0.002439f, -0.001055f, -0.000224f, 0.000015f, 0.002466f, -0.000865f, 0.000546f, -0.000072f, -0.000569f, 0.000500f, 0.000052f, -0.000596f, -0.000683f, -0.000529f, -0.000724f, 0.000254f, 0.000943f, -0.000478f, -0.000796f, 0.000684f, 0.000619f, -0.000109f, -0.000370f, 0.000735f, -0.000652f, 0.000726f, 0.001021f, -0.000625f, 0.000572f, -0.000338f, 0.000883f, -0.000564f, -0.000710f, -0.000515f, 0.000508f, 0.001687f, -0.000125f, -0.000845f, -0.000004f, -0.000012f, -0.001180f, 0.000690f, 0.000202f, -0.000113f, 0.001215f, 0.002380f, 0.001211f, -0.001147f, -0.000296f, 0.000822f, 0.000112f, -0.000144f, -0.000971f, -0.000517f, -0.000807f, 0.001426f, 0.000273f, -0.000567f, -0.000648f, 0.001165f, -0.000763f, 0.001345f, -0.000660f, 0.000250f, -0.000321f, 0.002187f, 0.000114f, 0.000229f, 0.000914f, -0.000796f, 0.001481f, -0.000925f, -0.000330f, -0.001464f, 0.001447f, 0.000933f, -0.001277f, 0.000301f, -0.000373f, 0.000133f, 0.000719f, 0.000608f, -0.001247f, -0.000301f, -0.001449f, 0.001382f, 0.000778f, -0.000189f, -0.000500f, -0.001044f, -0.000123f, 0.001099f, -0.001222f, -0.000328f, -0.001346f, 0.000964f, -0.000919f, 0.000780f, 0.001845f, -0.001337f, -0.000490f, -0.000907f, -0.000724f, -0.000173f, 0.000175f, 0.000321f, 0.001164f, 0.000926f, 0.002162f, -0.001364f, -0.000463f, -0.001003f, -0.000615f, -0.000862f, -0.000730f, -0.000056f, -0.000122f, -0.000026f, -0.000938f, -0.000537f, -0.000980f, 0.001203f, 0.000004f, 0.001889f, 0.000271f, -0.001090f, -0.000627f, 0.001419f, 0.000619f, -0.000264f, -0.000710f, 0.001683f, 0.000530f, 0.000502f, 0.001000f, -0.000534f, 0.000356f, -0.001379f, 0.000760f, -0.000676f, -0.000697f, -0.000908f, -0.000572f, -0.001008f, -0.000245f, 0.000231f, -0.000165f, 0.000858f, 0.000120f, -0.000370f, -0.001112f, 0.000046f, -0.000395f, 0.000211f, -0.000043f, -0.001003f, 0.001697f, 0.000037f, 0.000772f, -0.000001f, 0.000594f, 0.000866f, 0.000970f, 0.001481f, 0.000018f, -0.000076f, -0.000089f, -0.001164f, -0.000513f, -0.001016f, 0.001235f, 0.000340f, -0.000484f, -0.000674f, 0.000328f, -0.001242f, 0.000539f, -0.001083f, -0.000388f, -0.001217f, 0.000454f, 0.001388f, 0.000058f, -0.000265f, -0.000846f, -0.000750f, -0.000630f, -0.000957f, 0.000257f, -0.000128f, 0.000201f, -0.000269f, -0.001103f, 0.000299f, -0.000529f, 0.001307f, 0.001145f, 0.000028f, 0.000028f, -0.001257f, 0.000434f, 0.001703f, -0.000416f, 0.000148f, -0.000990f, -0.000519f, -0.000869f, -0.000533f, -0.000850f, -0.000481f, -0.000932f, 0.000737f, 0.000232f, -0.000703f, -0.000759f, 0.000034f, 0.000655f, -0.000236f, -0.000684f, -0.000682f, -0.000358f, 0.000108f, -0.000584f, 0.001086f, -0.000958f, 0.001531f, 0.002039f, -0.000950f, 0.000563f, -0.000936f, -0.000460f, -0.000744f, -0.000563f, -0.000544f, 0.000338f, -0.000884f, 0.000796f, 0.001592f, 0.000053f, -0.000882f, -0.000422f, -0.000766f, -0.000441f, -0.000737f, -0.000417f, -0.000745f, -0.000344f, -0.000861f, 0.000273f, 0.000764f, 0.001135f, -0.000640f, 0.000331f, 0.000736f, -0.001088f, 0.000420f, -0.000324f, -0.000593f, -0.000517f, -0.000599f, 0.000984f, 0.003110f, 0.000059f, -0.000573f, -0.000725f, 0.000082f, -0.000112f, -0.000711f, -0.000452f, -0.000615f, 0.000365f, -0.000548f, 0.000319f, 0.001230f, -0.000302f, -0.000094f, 0.000926f, 0.001029f, -0.000297f, -0.000693f, 0.000316f, 0.000141f, -0.000106f, -0.000883f, 0.000446f, -0.000472f, -0.000034f, 0.000057f, -0.000907f, 0.000058f, 0.000564f, 0.001875f, 0.000948f, -0.000802f, -0.000453f, -0.000584f, -0.000559f, -0.000439f, -0.000775f, 0.000540f, 0.002785f, 0.000133f, 0.000341f, -0.001178f, 0.000633f, 0.000513f, 0.000252f, -0.000032f, -0.000981f, 0.000316f, -0.000091f, 0.000115f, 0.000842f, -0.001221f, 0.001640f, 0.000920f, 0.000068f, -0.000841f, 0.000116f, -0.000257f, -0.000635f, -0.000566f, -0.000536f, -0.000465f, 0.001090f, 0.000868f, 0.000019f, -0.000231f, -0.000766f, 0.000295f, 0.000473f, -0.000224f, -0.000638f, -0.000557f, -0.000468f, 0.000092f, 0.000775f, 0.000291f, -0.000032f, 0.001276f, 0.000187f, -0.000533f, -0.000702f, 0.000048f, 0.000098f, -0.000957f, 0.000022f, 0.000067f, -0.000040f, 0.001203f, -0.001237f, 0.001191f, -0.000478f, 0.000118f, -0.000116f, 0.000133f, -0.000466f, -0.000673f, -0.000212f, 0.000292f, -0.000871f, -0.000206f, -0.000995f, 0.000981f, 0.000331f, 0.000588f, 0.000512f, -0.000255f, -0.000361f, -0.000730f, -0.000333f, -0.000719f, -0.000299f, -0.000764f, -0.000022f, 0.000273f, 0.000222f, -0.000720f, -0.000272f, -0.000744f, -0.000136f, -0.001026f, 0.001661f, 0.000414f, 0.000258f, -0.000435f, -0.000065f, 0.000257f, -0.000958f, 0.000393f, 0.000675f, -0.000783f, 0.000354f, -0.000327f, -0.000568f, -0.000329f, -0.000723f, 0.001500f, -0.000117f, 0.000411f, -0.000022f, -0.000426f, 0.000952f, -0.000377f, -0.000190f, -0.000887f, 0.000350f, 0.000746f, -0.000451f, -0.000707f, 0.000156f, 0.000350f, 0.000448f, 0.001878f, -0.000514f, -0.000528f, -0.000122f, 0.000069f, -0.000002f, -0.000140f, -0.000527f, 0.000267f, 0.000180f, 0.000210f, -0.000762f, -0.000262f, -0.000591f, -0.000331f, 0.000750f, -0.000481f, 0.000775f, 0.001121f, -0.000697f, -0.000447f, 0.000533f, -0.000223f, 0.000266f, 0.002323f, -0.000327f, 0.000048f, -0.000033f, 0.001138f, 0.000221f, -0.000951f, 0.000547f, 0.000331f, -0.000302f, 0.000843f, 0.000051f, -0.000483f, -0.000569f, -0.000430f, -0.000542f, -0.000449f, 0.001058f, 0.000761f, -0.000446f, -0.000619f, 0.000159f, 0.002102f, -0.000199f, 0.000494f, -0.000819f, -0.000210f, 0.000206f, 0.000535f, -0.000852f, -0.000203f, -0.000870f, 0.000392f, 0.000750f, -0.000221f, -0.000867f, 0.000413f, 0.000665f, -0.000909f, -0.000069f, 0.000332f, -0.000622f, -0.000478f, 0.000147f, -0.000170f, 0.000498f, 0.001040f, -0.001117f, 0.000531f, -0.000647f, 0.001214f, 0.000510f, 0.000402f, -0.000181f, -0.000447f, 0.000501f, -0.000682f, -0.000427f, -0.000570f, 0.000154f, 0.000457f, -0.000045f, -0.000788f, 0.000062f, 0.000440f, -0.000919f, 0.000423f, 0.000374f, -0.000214f, 0.000910f, -0.000573f, -0.000414f, 0.001471f, 0.000221f, 0.000392f, -0.000977f, 0.001050f, 0.000026f, -0.000795f, 0.000319f, -0.000425f, -0.000581f, -0.000418f, -0.000656f, 0.001268f, 0.000314f, -0.000744f, -0.000369f, -0.000594f, -0.000426f, -0.000524f, -0.000496f, -0.000040f, 0.000470f, 0.000589f, 0.000053f, 0.000383f, -0.000395f, 0.000246f, -0.000574f, -0.000411f, -0.000574f, -0.000059f, 0.001368f, -0.000225f, -0.000559f, -0.000454f, -0.000170f, 0.002621f, 0.000407f, -0.000700f, -0.000404f, -0.000539f, -0.000454f, -0.000488f, -0.000469f, -0.000451f, -0.000482f, -0.000394f, 0.000233f, -0.000512f, -0.000175f, 0.000513f, 0.000543f, -0.000348f, -0.000704f, 0.000854f, 0.000828f, -0.000652f, -0.000480f, 0.000968f, 0.000833f, -0.000706f, -0.000308f, -0.000296f, 0.000655f, 0.000050f, -0.000394f, -0.000560f, -0.000310f, -0.000578f, -0.000272f, -0.000596f, -0.000213f, -0.000665f, 0.001487f, -0.000118f, -0.000356f, -0.000670f, 0.000932f, 0.000971f, -0.000772f, -0.000181f, -0.000644f, 0.000003f, 0.000928f, 0.000000f, 0.000173f, 0.001107f, 0.002464f, 0.000037f, 0.001196f, -0.000384f, -0.000433f, 0.001072f, -0.000486f, -0.000599f, 0.000221f, 0.001020f, -0.000018f, -0.000228f, -0.000788f, 0.000744f, 0.000136f, -0.000395f, -0.000631f, -0.000278f, -0.000769f, 0.000471f, 0.001574f, -0.000024f, -0.000015f, 0.000245f, -0.000250f, 0.000529f, 0.000212f, -0.000085f, -0.000042f, -0.000399f, 0.000538f, -0.000452f, 0.000387f, 0.000467f, -0.000010f, 0.000742f, 0.000975f, 0.000073f, 0.001133f, -0.000555f, 0.000502f, -0.000922f, -0.000228f, -0.000061f, -0.000272f, -0.000683f, 0.001011f, -0.000533f, -0.000267f, -0.000229f, 0.000449f, 0.001370f, 0.000590f, 0.000606f, -0.000716f, 0.000556f, -0.000145f, -0.000569f, -0.000581f, -0.000524f, -0.000618f, 0.001475f, -0.000248f, 0.000035f, 0.000025f, -0.000187f, -0.000524f, -0.000663f, -0.000376f, -0.000801f, 0.000037f, 0.000275f, -0.000242f, 0.000485f, -0.000958f, -0.000166f, -0.000964f, 0.000453f, 0.001140f, 0.000718f, 0.001557f, -0.000641f, -0.000196f, -0.001076f, 0.000728f, 0.000807f, -0.000250f, 0.001752f, -0.000397f, -0.000167f, 0.000338f, -0.000363f, -0.000519f, -0.000730f, -0.000279f, -0.001045f, 0.000797f, 0.001468f, -0.000133f, 0.000215f, 0.000063f, 0.000536f, -0.000129f, 0.000840f, -0.001013f, -0.000304f, 0.000115f, 0.000938f, -0.000833f, -0.000331f, -0.000110f, -0.000425f, -0.000724f, 0.000721f, -0.000817f, -0.000143f, 0.001336f, -0.000209f, -0.000561f, -0.000759f, 0.001033f, 0.000312f, -0.000015f, -0.000060f, -0.000371f, 0.000428f, 0.000868f, -0.000596f, -0.000717f, -0.000359f, -0.000861f, -0.000065f, 0.000168f, 0.000117f, 0.000252f, -0.000161f, -0.001007f, 0.000722f, 0.000407f, -0.000717f, -0.000586f, -0.000233f, -0.000063f, -0.000606f, -0.000597f, -0.000207f, 0.000004f, 0.000363f, -0.000596f, 0.000164f, 0.000217f, -0.000091f, 0.000065f, 0.000100f, -0.000872f, 0.001429f, 0.000362f, 0.000978f, 0.000195f, -0.000220f, -0.000052f, -0.000651f, 0.000409f, -0.000390f, 0.000651f, -0.001082f, 0.000329f, 0.000570f, -0.000087f, 0.000789f, -0.001006f, -0.000154f, -0.000978f, 0.000547f, -0.000393f, -0.000031f, 0.000017f, -0.000384f, 0.001104f, 0.000087f, -0.000385f, -0.000698f, 0.000153f, 0.000994f, 0.000205f, 0.000763f, -0.000527f, 0.000398f, 0.000345f, -0.000870f, -0.000267f, -0.000831f, 0.001164f, 0.000714f, -0.000038f, -0.000411f, -0.000614f, -0.000457f, 0.000265f, -0.000764f, 0.000334f, 0.000591f, -0.001057f, 0.000399f, -0.000671f, -0.000012f, 0.000154f, 0.000184f, 0.000317f, -0.000589f, 0.000809f, -0.000315f, -0.000612f, -0.000384f, -0.000755f, 0.000977f, 0.000136f, -0.000059f, -0.000625f, 0.000415f, 0.000749f, -0.000529f, -0.000477f, -0.000455f, 0.000313f, -0.000716f, 0.000231f, -0.000399f, 0.000882f, -0.000030f, 0.000139f, -0.000751f, 0.000383f, -0.000463f, -0.000408f, -0.000673f, 0.000300f, 0.001906f, 0.000807f, 0.000126f, -0.000916f, -0.000019f, 0.000019f, -0.000042f, -0.000140f, -0.000456f, -0.000678f, 0.000085f, 0.000113f, -0.000528f, -0.000631f, 0.000337f, 0.001232f, -0.000666f, 0.001287f, -0.000316f, 0.000391f, -0.000524f, -0.000578f, 0.000096f, 0.000490f, 0.000122f, -0.000281f, -0.000586f, -0.000505f, 0.000720f, 0.001081f, -0.000485f, 0.002938f, 0.000565f, 0.000070f, 0.000925f, -0.000123f, 0.000245f, -0.000670f, -0.000500f, -0.000539f, 0.000290f, -0.000491f, -0.000570f, -0.000451f, -0.000624f, -0.000022f, 0.000405f, -0.000180f, -0.000662f, 0.000592f, 0.000602f, 0.000423f, -0.000259f, 0.000247f, 0.000389f, 0.000482f, 0.000706f, -0.000394f, 0.000268f, 0.001668f, -0.000395f, -0.000452f, -0.000685f, 0.000245f, -0.000603f, -0.000462f, -0.000608f, -0.000437f, 0.000954f, 0.000416f, -0.000993f, 0.000335f, 0.001454f, 0.000727f, -0.000204f, -0.000340f, 0.000559f, -0.000651f, 0.000170f, -0.000630f, 0.000530f, 0.000008f, -0.000748f, 0.001037f, -0.000369f, 0.000189f, 0.000272f, -0.000850f, -0.000383f, -0.000433f, 0.000504f, 0.000335f, -0.000552f, 0.000561f, -0.000718f, -0.000487f, -0.000161f, 0.000210f, 0.000650f, 0.000104f, 0.000373f, -0.000408f, 0.000095f, -0.000420f, 0.000684f, 0.000151f, 0.000197f, -0.000016f, 0.000115f, -0.000701f, -0.000501f, 0.002255f, 0.000589f, 0.002196f, -0.000532f, -0.000640f, -0.000586f, -0.000548f, 0.000053f, -0.000543f, 0.000243f, 0.001168f, -0.000182f, 0.000129f, -0.000690f, -0.000574f, -0.000433f, 0.000201f, -0.000026f, 0.001363f, 0.000343f, -0.000332f, -0.000945f, 0.000544f, -0.000373f, 0.000927f, -0.000202f, 0.000056f, -0.000474f, 0.001918f, 0.001494f, -0.000957f, 0.000229f, -0.000480f, 0.000074f, 0.000984f, 0.000392f, -0.000140f, 0.000167f, -0.001212f, 0.000539f, -0.000880f, -0.000257f, -0.001172f, 0.000855f, 0.000035f, 0.000026f, -0.000599f, -0.000438f, 0.000719f, -0.000235f, 0.001252f, -0.000234f, -0.000844f, 0.002021f, 0.000750f, -0.000683f, 0.000274f, -0.000282f, -0.000196f, 0.000021f, -0.000142f, 0.000893f, -0.001139f, 0.000203f, -0.000685f, -0.000448f, 0.000208f, 0.000329f, 0.000587f, -0.000756f, 0.000256f, -0.001017f, 0.000527f, -0.000271f, 0.001201f, -0.000034f, 0.001789f, -0.000377f, -0.001034f, 0.000221f, 0.000875f, -0.001210f, 0.000298f, 0.001102f, -0.001333f, 0.000915f, -0.000672f, -0.000493f, 0.000697f, 0.000748f, 0.000366f, 0.000009f, 0.001098f, 0.000193f, 0.000272f, -0.000690f, -0.000838f, -0.000429f, -0.000191f, -0.000339f, 0.001676f, -0.000169f, -0.000304f, -0.000765f, -0.000397f, 0.000546f, -0.000156f, 0.000119f, -0.001043f, -0.000132f, -0.000050f, 0.000151f, 0.000638f, 0.000844f, 0.001474f, -0.000577f, 0.000305f, -0.000843f, 0.000026f, -0.000827f, -0.000784f, 0.000094f, -0.000249f, -0.000864f, -0.000607f, -0.000701f, -0.000650f, 0.000185f, -0.000976f, -0.000156f, -0.000252f, 0.000798f, 0.000015f, 0.000538f, 0.000846f, -0.000573f, -0.000293f, 0.000236f, 0.000305f, -0.000742f, 0.000525f, -0.000503f, -0.000824f, 0.000773f, 0.001115f, -0.000035f, 0.000024f, -0.000042f, 0.000965f, -0.000889f, -0.000502f, -0.000785f, -0.000578f, -0.000607f, 0.000518f, -0.000514f, 0.000209f, 0.000599f, -0.000856f, 0.000898f, 0.000048f, 0.000687f, -0.000802f, 0.000575f, -0.000979f, 0.001474f, 0.000162f, -0.000819f, -0.000801f, 0.000834f, 0.002209f, -0.001044f, -0.000284f, -0.000351f, -0.000193f, -0.000131f, 0.000507f, -0.000620f, 0.001038f, 0.000827f, -0.000180f, -0.001231f, 0.000777f, 0.001910f, -0.001146f, -0.000225f, 0.000208f, 0.000647f, -0.000801f, 0.000145f, -0.001127f, 0.001464f, -0.000515f, 0.000090f, -0.001113f, -0.000161f, -0.001356f, 0.001016f, -0.000281f, 0.001013f, 0.002995f, -0.000051f, 0.000330f, 0.000179f, 0.001853f, -0.000364f, -0.000257f, -0.001125f, -0.000314f, -0.000066f, 0.000517f, -0.000178f, 0.001608f, 0.000655f, 0.000595f, -0.001103f, -0.000498f, -0.000955f, -0.000441f, 0.000369f, -0.000798f, 0.000005f, 0.000410f, -0.001084f, -0.000370f, 0.000743f, -0.000310f, 0.000069f, -0.000577f, -0.000854f, 0.000846f, -0.000620f, -0.000054f, -0.000118f, 0.000753f, 0.000163f, -0.001063f, 0.000136f, 0.000654f, 0.001591f, 0.000782f, -0.000815f, 0.000590f, -0.000501f, -0.000211f, -0.000206f, 0.000652f, 0.000489f, -0.001089f, -0.000496f, -0.000912f, 0.000028f, 0.000147f, 0.000110f, -0.000396f, 0.000488f, 0.001988f, 0.001202f, -0.001052f, 0.000688f, -0.000894f, -0.000416f, -0.001180f, 0.000116f, 0.000378f, 0.001166f, 0.000135f, -0.000489f, 0.001861f, -0.000442f, -0.000370f, -0.001028f, -0.000243f, -0.000344f, -0.001065f, 0.000644f, 0.001091f, 0.003816f, 0.001374f, -0.000324f, 0.001546f, -0.000949f, 0.000392f, -0.000709f, -0.000036f, -0.001519f, 0.000581f, -0.000190f, -0.000212f, 0.001453f, -0.000323f, 0.000280f, 0.000588f, 0.000006f, 0.001184f, -0.001759f, 0.000971f, -0.000712f, 0.001310f, 0.000181f, -0.000686f, -0.000090f, 0.000225f, 0.000834f, -0.000155f, -0.000649f, -0.000084f, -0.000180f, 0.000224f, -0.000342f, 0.000842f, 0.000213f, -0.001082f, -0.000811f, 0.000267f, -0.000068f, -0.000319f, -0.000918f, 0.000652f, 0.000613f, 0.000877f, 0.000695f, -0.001341f, 0.000103f, -0.001085f, -0.000664f, -0.000465f, 0.000216f, -0.000959f, -0.000274f, -0.000634f, -0.000130f, -0.000289f, -0.000768f, 0.000106f, -0.001038f, -0.000607f, -0.000907f, -0.000674f, 0.000368f, 0.000636f, 0.000651f, 0.000199f, -0.000194f, 0.000925f, 0.000640f, 0.000075f, -0.000442f, 0.000109f, -0.000675f, 0.000831f, 0.000458f, 0.000289f, -0.000365f, 0.000071f, -0.000153f, 0.000657f, -0.000403f, -0.000813f, -0.000936f, 0.000008f, 0.000193f, 0.001748f, 0.000540f, 0.000570f, 0.000222f, -0.000463f, -0.000292f, -0.000349f, -0.000877f, -0.000746f, -0.000822f, -0.000716f, -0.000640f, 0.001448f, 0.000656f, -0.000023f, 0.000797f, -0.000682f, -0.000744f, -0.000277f, 0.001080f, 0.000277f, 0.001859f, 0.000136f, -0.000878f, -0.000714f, -0.000162f, -0.000840f, 0.001501f, -0.000956f, -0.000233f, -0.000561f, -0.000817f, -0.000744f, -0.000743f, -0.000763f, 0.000075f, -0.000146f, -0.000180f, -0.000992f, 0.000212f, -0.000318f, 0.000212f, 0.000338f, -0.000169f, -0.000118f, 0.000253f, 0.001629f, -0.001278f, -0.000316f, -0.001072f, -0.000097f, 0.000354f, -0.000505f, -0.000886f, -0.000431f, -0.000987f, 0.000071f, 0.000990f, -0.000182f, -0.000928f, -0.000301f, -0.000375f, 0.000433f, -0.000578f, 0.001145f, 0.000739f, -0.001075f, -0.000397f, -0.000781f, -0.000529f, -0.000645f, -0.000641f, -0.000262f, 0.000334f, 0.000706f, 0.000320f, 0.000829f, 0.002260f, 0.000002f, -0.000655f, -0.000074f, -0.000127f, -0.000431f, -0.000819f, 0.000786f, -0.000425f, -0.000614f, -0.000697f, -0.000218f, 0.001272f, -0.000160f, -0.000533f, -0.000637f, 0.000086f, 0.000226f, -0.000343f, -0.000750f, 0.001263f, -0.000706f, -0.000357f, 0.000270f, 0.000485f, 0.000099f, -0.000338f, -0.000134f, 0.000841f, 0.001627f, -0.000130f, 0.000004f, -0.000621f, 0.000074f, -0.000275f, -0.000677f, -0.000018f, -0.000479f, -0.000699f, -0.000341f, -0.000920f, 0.001129f, 0.000994f, 0.000033f, -0.000022f, 0.000784f, -0.000116f, -0.000220f, 0.000080f, -0.000826f, -0.000457f, -0.000164f, -0.000370f, 0.000817f, 0.000559f, -0.001043f, 0.000405f, 0.000095f, 0.000568f, -0.001249f, 0.000155f, 0.000336f, 0.002157f, 0.000465f, 0.000008f, 0.000185f, -0.000438f, 0.000140f, -0.000594f, 0.000388f, 0.000424f, 0.000898f, 0.000850f, -0.000599f, 0.001038f, -0.000157f, -0.000184f, -0.000914f, -0.000058f, 0.000547f, 0.000146f, 0.000252f, 0.001300f, 0.001216f, -0.000737f, 0.000149f, -0.000813f, 0.001353f, -0.000802f, -0.000445f, -0.000826f, -0.000004f, 0.000010f, -0.000267f, -0.000233f, 0.000956f, -0.000290f, -0.000652f, 0.000376f, -0.000063f, 0.000431f, -0.000022f, 0.001318f, 0.001436f, 0.001229f, -0.000025f, -0.000101f, 0.001365f, -0.000287f, -0.000911f, -0.000370f, -0.000208f, -0.000537f, -0.000926f, 0.000414f, 0.000085f, -0.000166f, -0.000828f, -0.000142f, 0.000098f, -0.000069f, -0.001053f, 0.001464f, 0.000355f, -0.000068f, -0.000177f, 0.000476f, 0.001085f, -0.000131f, -0.000868f, -0.000389f, -0.000233f, -0.000688f, -0.000677f, -0.000382f, 0.000847f, 0.001224f, -0.001219f, 0.000385f, -0.000421f, 0.000598f, -0.000116f, -0.000125f, 0.000842f, 0.000508f, 0.001130f, 0.000341f, -0.000094f, -0.000425f, -0.000217f, -0.000902f, 0.000468f, -0.000562f, -0.000761f, 0.000228f, -0.000663f, 0.000286f, 0.000151f, 0.000377f, 0.000106f, -0.000325f, -0.000716f, -0.000607f, -0.000729f, 0.000682f, 0.000706f, -0.000154f, 0.000374f, -0.000216f, -0.000580f, 0.000557f, 0.000612f, 0.000781f, 0.000777f, -0.000948f, -0.000605f, -0.000599f, -0.000830f, -0.000029f, 0.000341f, -0.000438f, 0.000349f, -0.000773f, -0.000692f, 0.000040f, 0.000651f, 0.000227f, 0.000213f, 0.000050f, 0.000049f, -0.000897f, -0.000172f, -0.000225f, -0.000203f, -0.000672f, -0.000275f, 0.001233f, 0.000293f, -0.000143f, -0.000514f, -0.000714f, 0.001058f, 0.001283f, -0.001012f, 0.000203f, -0.000450f, -0.000590f, -0.000705f, -0.000432f, -0.000868f, 0.000150f, 0.000060f, -0.000044f, -0.000010f, 0.000140f, -0.000394f, -0.000094f, -0.000245f, -0.000353f, 0.000187f, -0.000079f, -0.000145f, -0.000215f, 0.001404f, -0.000186f, -0.000126f, 0.000472f, 0.000487f, -0.000004f, -0.000622f, -0.000603f, 0.000483f, 0.000175f, -0.000806f, -0.000316f, -0.000891f, 0.000158f, -0.000454f, 0.000536f, 0.000836f, 0.000483f, 0.000748f, 0.000585f, -0.000277f, -0.000704f, -0.000637f, 0.001463f, 0.001657f, -0.000888f, -0.000241f, -0.000963f, 0.000199f, 0.000652f, -0.000558f, -0.000448f, 0.000196f, -0.000077f, 0.000488f, 0.000099f, 0.001605f, -0.001379f, 0.001014f, 0.000128f, -0.000693f, -0.000569f, -0.000582f, 0.001715f, 0.000432f, -0.000843f, -0.000347f, -0.000858f, 0.000457f, -0.000594f, -0.000472f, -0.000683f, -0.000413f, -0.000700f, -0.000367f, -0.000701f, 0.000914f, -0.000533f, -0.000206f, -0.000253f, -0.000117f, -0.000190f, -0.000622f, -0.000389f, 0.000554f, -0.000571f, 0.000570f, 0.000374f, -0.000450f, 0.000035f, -0.000611f, 0.000135f, 0.000397f, -0.000648f, 0.000047f, 0.000881f, 0.000204f, -0.000435f, -0.000462f, -0.000558f, -0.000366f, -0.000643f, 0.000424f, 0.000535f, 0.000687f, 0.000972f, -0.000296f, 0.000603f, -0.000146f, -0.000272f, 0.000346f, -0.000700f, 0.000894f, -0.000035f, -0.000120f, 0.000780f, 0.000803f, -0.000048f, -0.000766f, 0.000483f, 0.000082f, -0.000029f, 0.000186f, 0.002720f, 0.000088f, 0.000334f, -0.000865f, 0.000413f, 0.000346f, -0.000149f, -0.000255f, -0.000085f, -0.000299f, -0.000639f, -0.000492f, 0.001069f, 0.000763f, -0.000802f, 0.000815f, 0.000490f, -0.000288f, -0.000800f, 0.000490f, -0.000058f, -0.000321f, -0.000259f, -0.000847f, 0.000077f, 0.000465f, -0.000471f, 0.000283f, -0.000237f, 0.000502f, -0.000725f, 0.000674f, 0.000322f, -0.000731f, -0.000341f, 0.001328f, 0.000530f, -0.000167f, 0.000232f, 0.000299f, -0.000010f, -0.000904f, -0.000253f, -0.000871f, 0.000060f, -0.000292f, -0.000573f, -0.000215f, 0.000017f, -0.000327f, -0.000398f, -0.000730f, 0.000068f, 0.000029f, 0.001128f, -0.000066f, -0.000561f, -0.000625f, 0.000176f, 0.000704f, -0.000458f, 0.000456f, -0.000886f, 0.000115f, 0.000905f, -0.000120f, 0.000516f, 0.000607f, 0.000631f, -0.000452f, -0.000603f, -0.000415f, 0.000161f, -0.000653f, 0.000564f, -0.000285f, -0.000650f, -0.000374f, 0.000174f, 0.000079f, 0.001061f, -0.000386f, -0.000402f, 0.000674f, -0.000587f, 0.000482f, 0.000411f, -0.000415f, -0.000684f, 0.000517f, 0.000648f, 0.000033f, -0.000015f, -0.000726f, -0.000047f, -0.000274f, -0.000148f, -0.000149f, 0.000019f, 0.000170f, -0.000342f, 0.000474f, 0.000030f, 0.001688f, -0.000407f, 0.000279f, -0.000186f, 0.001037f, 0.000457f, 0.000649f, 0.000427f, -0.000482f, -0.000303f, -0.000722f, 0.001099f, 0.000072f, -0.000815f, 0.000961f, 0.000315f, -0.000863f, 0.000679f, -0.001232f, 0.000300f, -0.000607f, 0.000136f, -0.000519f, 0.000067f, 0.000480f, -0.000340f, -0.000093f, -0.000197f, 0.000042f, 0.000825f, 0.000702f, -0.000838f, -0.000366f, -0.000707f, -0.000033f, -0.000302f, -0.000636f, -0.000414f, -0.000639f, 0.000069f, -0.000168f, -0.000099f, -0.000657f, -0.000297f, 0.000216f, -0.000170f, -0.000696f, -0.000203f, 0.000070f, 0.001991f, -0.000482f, -0.000416f, -0.000589f, -0.000350f, 0.000442f, 0.000017f, -0.000322f, -0.000309f, -0.000030f, 0.000267f, 0.001039f, 0.000279f, -0.000111f, 0.000927f, 0.000464f, -0.000087f, -0.000141f, -0.000459f, 0.000384f, 0.000079f, -0.001059f, 0.000738f, 0.000484f, 0.000625f, -0.000044f, 0.000284f, 0.000190f, -0.000071f, 0.000667f, 0.000582f, 0.000522f, -0.000838f, -0.000384f, -0.000454f, 0.000777f, -0.000677f, 0.000503f, -0.000212f, 0.000868f, -0.000427f, -0.000427f, -0.000635f, -0.000402f, -0.000115f, 0.000144f, 0.001726f, 0.001488f, 0.000538f, 0.001059f, 0.000854f, -0.000035f, 0.000036f, -0.000309f, 0.000340f, 0.000296f, -0.000093f, -0.000518f, -0.000118f, 0.000727f, -0.000599f, 0.000216f, 0.000351f, 0.000554f, 0.000249f, 0.001003f, 0.000236f, 0.000501f, -0.000205f, -0.000345f, 0.000448f, 0.000790f, -0.001168f, 0.001672f, 0.000145f, 0.000095f, 0.000693f, 0.000707f, 0.000100f, -0.000092f, 0.000110f, 0.000530f, -0.000231f, -0.000004f, -0.000361f, 0.000061f, 0.000283f, -0.000527f, 0.000470f, -0.001375f, 0.000674f, 0.000053f, -0.000105f, -0.000119f, -0.000756f, 0.002224f, 0.000080f, -0.000266f, -0.000864f, 0.001004f, 0.000305f, -0.000595f, 0.000804f, 0.000459f, -0.000316f, -0.000680f, 0.000832f, -0.000466f, 0.000976f, -0.000686f, 0.000525f, 0.001071f, -0.000657f, -0.000925f, -0.000116f, -0.000893f, -0.000108f, -0.000832f, 0.001098f, 0.001194f, 0.001554f, 0.000569f, 0.000224f, -0.000361f, -0.001124f, -0.000429f, 0.000127f, 0.000329f, -0.001291f, -0.000409f, -0.000148f, 0.001505f, -0.001271f, 0.001403f, -0.000917f, -0.000398f, -0.000563f, 0.000292f, -0.000350f, -0.000338f, -0.000793f, 0.001156f, -0.000687f, 0.001708f, -0.000962f, -0.000172f, 0.000312f, -0.000556f, 0.001433f, 0.000557f, 0.000320f, -0.001196f, -0.000186f, -0.000323f, -0.000881f, -0.000691f, -0.000900f, 0.000529f, 0.000793f, -0.000040f, 0.001453f, 0.000332f, -0.000776f, -0.000851f, -0.000224f, -0.000871f, -0.000663f, -0.000502f, 0.001023f, 0.000841f, 0.000052f, -0.000041f, -0.000783f, -0.000796f, -0.000662f, -0.000914f, 0.000009f, 0.000360f, 0.001499f, -0.000455f, -0.000793f, -0.000795f, 0.000867f, 0.000300f, -0.000024f, -0.000552f, -0.000297f, 0.001651f, -0.000244f, 0.000438f, 0.001303f, -0.000504f, -0.000044f, -0.000280f, -0.000475f, 0.000849f, -0.000715f, 0.000479f, 0.001814f, -0.000586f, -0.000164f, -0.000272f, -0.001013f, -0.000521f, -0.001084f, 0.001190f, 0.000182f, -0.001037f, -0.000527f, 0.000328f, -0.000429f, -0.000252f, 0.000612f, 0.000142f, 0.001162f, -0.000686f, 0.000850f, -0.000365f, 0.000027f, -0.000296f, 0.001116f, -0.000338f, -0.000263f, -0.001103f, 0.000861f, 0.000487f, 0.000155f, -0.000992f, 0.000919f, 0.000535f, -0.000680f, -0.000237f, -0.000993f, -0.000583f, -0.000873f, -0.000216f, 0.000521f, -0.000324f, -0.000510f, 0.001673f, -0.000487f, -0.000225f, -0.000132f, 0.000836f, -0.000411f, 0.000143f, -0.000453f, 0.001280f, -0.000782f, -0.000721f, -0.000713f, -0.000765f, -0.000469f, 0.000277f, -0.000860f, 0.000240f, -0.000964f, -0.000234f, 0.001355f, -0.000330f, -0.000485f, -0.000420f, -0.000385f, -0.000671f, -0.000840f, 0.001151f, 0.000996f, -0.001253f, -0.000125f, -0.000103f, -0.000093f, -0.000088f, -0.000715f, -0.000459f, 0.000027f, -0.001165f, 0.002855f, 0.001087f, -0.001083f, -0.000380f, -0.000118f, -0.000614f, -0.000636f, -0.000650f, 0.000044f, 0.001066f, 0.000096f, 0.000762f, 0.000144f, 0.000789f, 0.000087f, -0.000008f, -0.000488f, 0.000075f, -0.000750f, 0.001374f, 0.000122f, -0.000588f, -0.000696f, -0.000557f, -0.000722f, -0.000243f, 0.000509f, 0.001737f, -0.000083f, -0.000094f, 0.000077f, 0.001291f, -0.000476f, -0.000415f, -0.000192f, 0.000081f, -0.000471f, -0.000191f, 0.001508f, 0.000189f, 0.000299f, -0.000251f, 0.000263f, 0.000154f, -0.000530f, 0.000298f, 0.001456f, -0.000803f, 0.001016f, 0.000687f, -0.000199f, -0.000134f, 0.000623f, -0.000548f, 0.000090f, -0.001266f, 0.001086f, 0.001417f, -0.001030f, -0.000614f, -0.000305f, -0.000089f, 0.000041f, 0.000206f, -0.000561f, 0.000822f, -0.000300f, 0.000388f, 0.000071f, -0.000256f, 0.000007f, -0.000594f, -0.000011f, -0.000826f, 0.000924f, 0.001069f, -0.000805f, -0.000627f, -0.000578f, -0.000168f, -0.000779f, 0.000181f, 0.000856f, 0.000761f, -0.000046f, 0.000318f, 0.000275f, -0.000190f, -0.000877f, 0.000199f, -0.000434f, 0.000667f, -0.000134f, -0.000370f, 0.000242f, -0.000686f, -0.000186f, 0.000954f, -0.001082f, -0.000388f, -0.000787f, 0.000089f, -0.001006f, 0.000347f, -0.000254f, -0.000103f, -0.000851f, -0.000282f, 0.000305f, -0.000354f, -0.000144f, 0.000216f, 0.000455f, 0.000833f, -0.000237f, 0.000711f, -0.000528f, -0.000156f, -0.000296f, -0.000472f, -0.000027f, 0.002225f, -0.000318f, 0.000858f, -0.000278f, -0.000732f, -0.000222f, -0.000476f, -0.000769f, -0.000338f, -0.000206f, 0.000510f, 0.000201f, -0.000783f, -0.000343f, 0.000781f, -0.000789f, 0.000071f, -0.000255f, 0.000955f, 0.000298f, 0.000199f, -0.000074f, 0.000159f, 0.000077f, 0.000366f, -0.000372f, -0.000261f, -0.000916f, 0.000996f, 0.000276f, -0.000248f, -0.000027f, -0.000084f, 0.000897f, 0.000609f, -0.000337f, -0.000230f, -0.000288f, 0.000354f, -0.000548f, 0.000318f, -0.000117f, -0.000112f, -0.000453f, -0.000227f, 0.000462f, -0.000072f, 0.000836f, -0.000351f, 0.000686f, -0.000335f, -0.000040f, 0.000614f, -0.000397f, -0.000080f, 0.000887f, -0.000813f, 0.001042f, -0.000255f, 0.001086f, -0.000337f, -0.000671f, -0.000654f, 0.000524f, 0.000409f, -0.000951f, -0.000021f, 0.000199f, 0.000221f, -0.000193f, -0.000757f, -0.000515f, -0.000634f, 0.000241f, -0.000054f, -0.000002f, 0.000237f, -0.000560f, -0.000088f, -0.000560f, -0.000522f, -0.000002f, -0.000656f, -0.000503f, 0.000280f, 0.000883f, -0.000513f, 0.000349f, -0.000405f, -0.000083f, -0.000740f, -0.000372f, 0.000312f, 0.000271f, 0.000191f, 0.000402f, 0.000185f, 0.000006f, -0.000553f, -0.000585f, -0.000468f, -0.000702f, 0.001001f, 0.000863f, -0.000503f, 0.000386f, 0.000726f, -0.000586f, -0.000016f, -0.000689f, 0.000610f, -0.000281f, -0.000043f, 0.001290f, 0.000705f, -0.000824f, 0.000088f, -0.000401f, -0.000044f, -0.000110f, -0.000404f, -0.000102f, -0.000484f, 0.001075f, -0.000632f, 0.000109f, -0.000385f, -0.000421f, -0.000110f, -0.000697f, -0.000177f, -0.000098f, 0.001649f, 0.000279f, -0.000851f, -0.000276f, -0.000747f, -0.000256f, -0.000832f, 0.000703f, 0.000313f, -0.000271f, 0.000662f, 0.000950f, -0.000538f, 0.000199f, 0.000529f, -0.000129f, -0.000173f, -0.000107f, 0.001125f, -0.000240f, 0.000537f, -0.000943f, 0.000420f, -0.000041f, 0.000068f, 0.000362f, -0.000610f, -0.000404f, 0.000155f, -0.000018f, -0.000562f, 0.000442f, 0.000505f, -0.000864f, -0.000191f, 0.000462f, 0.000792f, -0.000954f, 0.000023f, 0.000056f, 0.001211f, 0.000197f, 0.000282f, 0.000057f, -0.000163f, -0.000043f, -0.000405f, 0.000206f, -0.000736f, 0.000416f, 0.000476f, 0.000041f, -0.000112f, 0.000197f, -0.000203f, -0.000030f, -0.000118f, -0.000145f, -0.000440f, -0.000605f, -0.000119f, 0.000124f, -0.000216f, -0.000702f, 0.000669f, 0.000641f, -0.000868f, 0.000195f, 0.000462f, -0.000186f, -0.000049f, 0.000490f, 0.000918f, -0.000432f, 0.000155f, 0.000452f, -0.000185f, -0.000086f, 0.000777f, 0.000605f, -0.000580f, -0.000181f, -0.000449f, -0.000183f, -0.000326f, 0.000923f, 0.000768f, -0.000044f, -0.000352f, -0.000318f, 0.001380f, -0.000183f, 0.000009f, -0.000941f, 0.000481f, -0.000037f, -0.000641f, -0.000203f, 0.000299f, -0.000575f, 0.000577f, 0.000133f, 0.000459f, 0.001469f, -0.000255f, 0.000309f, -0.000488f, -0.000404f, 0.000536f, 0.001064f, -0.000060f, -0.000585f, -0.000225f, 0.000639f, 0.001334f, -0.000000f, -0.000901f, 0.000005f, -0.000178f, 0.000266f, -0.000049f, 0.000540f, -0.000185f, -0.000807f, -0.000358f, 0.000550f, 0.000548f, -0.000444f, -0.000513f, -0.000684f, -0.000031f, -0.000091f, -0.000509f, -0.000656f, -0.000452f, 0.000074f, 0.001038f, 0.000687f, -0.000071f, -0.000717f, -0.000025f, -0.000353f, 0.000946f, 0.000021f, 0.000553f, 0.000146f, 0.000394f, 0.000502f, -0.000011f, 0.000380f, -0.000642f, -0.000714f, 0.000480f, 0.000420f, -0.000504f, 0.000127f, -0.000073f, -0.000068f, -0.000497f, 0.000242f, 0.000355f, -0.000611f, -0.000728f, 0.000613f, 0.000822f, 0.000476f, 0.000107f, 0.000234f, 0.000024f, 0.000596f, -0.000544f, 0.000156f, -0.000990f, -0.000189f, -0.000385f, 0.000342f, -0.000469f, -0.000089f, -0.000714f, -0.000400f, 0.000059f, 0.000317f, -0.000871f, 0.000453f, 0.000326f, -0.000432f, 0.001149f, 0.000921f, 0.000250f, 0.000198f, 0.000902f, 0.000594f, 0.000864f, -0.000115f, 0.001103f, -0.000873f, -0.000033f, -0.000812f, -0.000286f, -0.000030f, -0.000311f, -0.000123f, -0.000770f, -0.000217f, -0.000215f, 0.000398f, -0.000318f, -0.000579f, 0.000382f, 0.001818f, 0.000873f, -0.000847f, -0.000259f, 0.000298f, 0.000227f, -0.000237f, 0.000751f, -0.001048f, -0.000420f, -0.000417f, -0.000048f, -0.000592f, 0.000093f, 0.001336f, -0.000830f, -0.000412f, 0.000932f, 0.001248f, -0.000338f, -0.000222f, 0.001287f, 0.000840f, -0.000166f, 0.000722f, 0.000253f, -0.000549f, 0.000876f, 0.000393f, 0.000246f, -0.000437f, 0.000493f, 0.000460f, -0.000535f, -0.000414f, 0.000082f, -0.000418f, -0.000256f, 0.000778f, 0.001697f, 0.000161f, -0.000308f, -0.000715f, -0.000505f, -0.000095f, -0.000091f, -0.000462f, -0.000612f, 0.000953f, -0.000313f, -0.000550f, -0.000841f, -0.000306f, -0.000517f, -0.000818f, -0.000371f, 0.000872f, -0.000256f, 0.000647f, 0.000172f, 0.001079f, -0.000310f, 0.000321f, -0.000874f, -0.000491f, 0.000515f, 0.000133f, -0.000933f, -0.000414f, -0.000924f, 0.000461f, -0.000477f, 0.000983f, 0.000893f, -0.000619f, -0.000137f, 0.000034f, -0.000171f, 0.000129f, -0.000195f, -0.000025f, -0.000958f, 0.000025f, 0.000680f, -0.000156f, -0.000490f, -0.000779f, 0.000382f, 0.000191f, -0.000731f, 0.001608f, 0.001078f, -0.000136f, -0.000582f, 0.000631f, 0.000421f, -0.000374f, -0.000188f, 0.000128f, -0.000674f, -0.000693f, 0.000985f, 0.000640f, -0.000493f, 0.000814f, 0.000551f, -0.000132f, -0.000068f, 0.000542f, -0.000516f, -0.000569f, 0.001136f, 0.000006f, -0.000644f, -0.000437f, 0.000345f, -0.000189f, -0.000211f, -0.000055f, -0.000503f, -0.000405f, 0.000109f, 0.000421f, 0.000344f, 0.000527f, -0.000288f, -0.000375f, -0.000013f, 0.000966f, 0.000815f, -0.000217f, -0.000885f, 0.000143f, 0.000272f, 0.000338f, 0.000082f, -0.000063f, -0.000592f, 0.000014f, -0.000238f, -0.000060f, -0.000492f, -0.000118f, -0.000011f, -0.000758f, 0.000309f, -0.000039f, -0.000074f, 0.000195f, -0.000260f, -0.000375f, -0.000862f, 0.000373f, -0.000322f, -0.000269f, -0.000662f, -0.000198f, -0.000871f, -0.000216f, 0.000299f, 0.000332f, 0.000030f, -0.000059f, 0.000802f, -0.000625f, 0.000584f, -0.000388f, -0.000461f, -0.000509f, -0.000352f, -0.000767f, -0.000565f, -0.000336f, 0.000989f, -0.000410f, -0.000171f, -0.000133f, 0.000603f, 0.000542f, -0.000234f, -0.000087f, 0.000709f, -0.000723f, -0.000299f, 0.000424f, -0.000817f, -0.000356f, -0.000278f, 0.000149f, -0.000567f, 0.000392f, 0.001134f, 0.000286f, 0.000938f, -0.000726f, 0.000472f, -0.000976f, 0.000899f, 0.000163f, 0.000211f, -0.000747f, 0.000608f, -0.000790f, 0.000291f, 0.000031f, -0.000400f, -0.000489f, -0.000130f, 0.000216f, -0.000519f, -0.000166f, 0.001085f, -0.000145f, 0.000167f, -0.000640f, 0.000018f, 0.000087f, -0.000220f, 0.000223f, -0.000298f, 0.000217f, 0.001183f, -0.000200f, -0.000481f, 0.001279f, 0.000252f, -0.000011f, -0.000397f, -0.000329f, -0.000362f, -0.000371f, -0.000046f, 0.000534f, 0.000153f, 0.000888f, -0.000209f, -0.000686f, -0.000652f, 0.000053f, 0.000133f, -0.000050f, 0.000144f, 0.000671f, -0.000075f, -0.000227f, -0.000503f, -0.000064f, -0.000243f, 0.000151f, 0.000870f, -0.000984f, 0.001397f, 0.000470f, 0.000849f, -0.000795f, 0.000032f, -0.000502f, -0.000211f, 0.001065f, 0.001375f, 0.000112f, -0.000614f, -0.000297f, 0.000573f, -0.000016f, -0.000924f, -0.000415f, -0.000191f, 0.000556f, -0.000529f, -0.000094f, 0.000317f, 0.000687f, -0.000051f, -0.000876f, -0.000214f, 0.000289f, -0.000205f, -0.000499f, 0.000781f, -0.000283f, -0.000356f, 0.001778f, 0.000331f, 0.001538f, -0.000136f, -0.000471f, -0.000244f, -0.000672f, -0.000337f, 0.000012f, 0.000362f, 0.000861f, 0.000880f, -0.000339f, 0.000255f, 0.000828f, -0.000694f, -0.000644f, -0.000597f, -0.000014f, 0.000249f, 0.000830f, 0.000416f, -0.000197f, 0.000088f, 0.000187f, -0.000012f, -0.000024f, 0.000569f, 0.001210f, -0.000433f, -0.000126f, 0.000646f, -0.001250f, 0.000015f, -0.000686f, 0.000225f, -0.000236f, -0.000035f, -0.000850f, -0.000127f, -0.000542f, 0.000275f, 0.000425f, -0.000774f, -0.000228f, -0.000730f, 0.000132f, 0.000603f, -0.000268f, -0.000869f, -0.000066f, 0.000044f, 0.000032f, -0.000214f, 0.000065f, 0.001142f, 0.000191f, -0.001137f, 0.000183f, -0.000504f, -0.000247f, 0.000355f, 0.000587f, -0.000508f, 0.000096f, 0.000199f, -0.000590f, -0.000322f, -0.000529f, -0.000421f, 0.001061f, -0.000038f, 0.000555f, -0.000544f, 0.000232f, 0.000224f, -0.000369f, -0.000286f, -0.000265f, 0.000107f, -0.000159f, 0.000197f, 0.000203f, 0.000087f, 0.000041f, -0.000030f, 0.000120f, 0.000794f, -0.000456f, -0.000306f, -0.000631f, 0.001813f, 0.000206f, -0.000726f, 0.000199f, -0.000059f, 0.000314f, -0.000281f, 0.001466f, -0.000790f, -0.000536f, -0.000331f, 0.000228f, -0.000717f, -0.000145f, -0.000248f, 0.000465f, 0.000322f, 0.000545f, -0.000894f, 0.000489f, 0.000155f, 0.000216f, -0.000456f, -0.000736f, 0.000023f, 0.000187f, -0.000261f, 0.000107f, -0.000786f, -0.000286f, -0.000488f, -0.000269f, -0.000534f, 0.000686f, -0.000356f, 0.000426f, 0.000887f, -0.000285f, 0.000951f, -0.000877f, 0.000991f, 0.000775f, -0.000785f, -0.000397f, -0.000041f, -0.000683f, -0.000263f, -0.000541f, 0.001055f, 0.000396f, -0.000028f, 0.000786f, -0.000457f, -0.000557f, -0.000780f, -0.000062f, -0.000850f, 0.001103f, 0.000167f, -0.000149f, -0.000363f, -0.000225f, -0.000167f, 0.000031f, 0.000253f, -0.000657f, 0.000220f, 0.000071f, -0.000243f, 0.000094f, 0.001066f, -0.000423f, -0.000270f, -0.000467f, 0.000113f, -0.000265f, 0.000054f, 0.000489f, -0.000714f, -0.000502f, 0.000294f, 0.000264f, 0.000059f, -0.000530f, 0.001110f, 0.000373f, -0.001109f, 0.000468f, -0.000487f, -0.000519f, -0.000253f, 0.000463f, -0.000279f, 0.000484f, 0.000344f, -0.000496f, 0.000114f, 0.000880f, 0.000024f, 0.000057f, 0.001596f, 0.000645f, -0.000303f, 0.000073f, -0.001081f, 0.000659f, 0.000379f, -0.000365f, 0.000290f, -0.000312f, 0.000002f, -0.000813f, -0.000276f, 0.000779f, -0.000777f, -0.000460f, -0.000264f, 0.000532f, 0.000742f, -0.000380f, -0.000143f, -0.000152f, -0.000488f, 0.000389f, -0.000792f, 0.000876f, 0.000064f, 0.000712f, 0.000232f, -0.000309f, 0.000041f, -0.000618f, 0.000024f, 0.000472f, -0.000541f, -0.000041f, -0.000433f, 0.000133f, -0.000226f, -0.000464f, -0.000163f, -0.000172f, 0.001549f, 0.001404f, 0.001130f, -0.000058f, 0.000299f, 0.000298f, -0.000711f, 0.000474f, 0.000914f, -0.000824f, -0.000457f, 0.000289f, 0.000305f, -0.000687f, 0.000460f, -0.000912f, -0.000437f, -0.000392f, -0.000497f, -0.000539f, 0.000601f, -0.000248f, 0.000052f, 0.000354f, -0.000266f, 0.000600f, 0.000445f, -0.000429f, 0.000375f, -0.000866f, -0.000309f, -0.000535f, -0.000235f, 0.000684f, -0.000305f, 0.000195f, -0.000375f, -0.000019f, -0.000149f, -0.000715f, -0.000388f, 0.000515f, -0.000664f, 0.000323f, 0.000628f, -0.000470f, 0.000133f, -0.000071f, -0.000111f, 0.000145f, -0.000260f, -0.000104f, -0.000541f, -0.000409f, 0.000764f, -0.000538f, -0.000320f, -0.000015f, -0.000325f, -0.000490f, 0.000292f, -0.000647f, -0.000505f, 0.000702f, 0.000668f, -0.000891f, 0.000799f, -0.000056f, -0.000481f, 0.000568f, 0.000925f, -0.000512f, -0.000204f, -0.000261f, 0.000029f, -0.000455f, 0.000732f, -0.000659f, -0.000372f, -0.000786f, 0.000616f, 0.000720f, -0.000089f, 0.000115f, 0.000069f, 0.000870f, -0.000756f, 0.000046f, -0.000669f, 0.000113f, 0.000553f, -0.000075f, -0.000148f, 0.000453f, -0.000263f, -0.000455f, -0.000220f, 0.000174f, -0.000120f, -0.000051f, -0.000226f, -0.000374f, -0.000552f, -0.000109f, -0.000257f, 0.000480f, 0.000154f, -0.000559f, -0.000507f, 0.000347f, -0.000378f, -0.000494f, 0.000470f, -0.000664f, -0.000353f, -0.000585f, 0.000192f, -0.000410f, 0.000802f, 0.000668f, -0.000258f, 0.000001f, 0.000447f, -0.000244f, -0.000693f, 0.000821f, -0.000073f, 0.000047f, 0.000014f, -0.000282f, 0.000315f, -0.000538f, -0.000394f, -0.000586f, 0.000037f, -0.000023f, -0.000232f, -0.000306f, -0.000562f, -0.000320f, -0.000579f, 0.000603f, 0.000525f, -0.000144f, -0.000164f, -0.000505f, 0.000048f, -0.000329f, 0.000627f, -0.000129f, 0.000066f, 0.001501f, 0.000642f, -0.000164f, 0.000256f, 0.000298f, -0.001009f, 0.001207f, 0.000830f, -0.000119f, -0.000460f, 0.000086f, 0.000819f, -0.000109f, -0.000549f, -0.000439f, -0.000371f, 0.000507f, 0.000109f, -0.000005f, -0.000193f, 0.000871f, 0.000646f, -0.000566f, 0.000238f, -0.000171f, 0.000980f, -0.000310f, 0.000869f, 0.000473f, -0.000659f, 0.000246f, -0.000448f, -0.000018f, -0.000034f, -0.000405f, 0.000118f, -0.000164f, 0.000102f, -0.000056f, -0.000390f, -0.000389f, -0.000642f, 0.000089f, 0.000402f, -0.000292f, 0.000360f, 0.000102f, -0.000538f, 0.000339f, 0.000508f, -0.000443f, 0.000047f, 0.000567f, 0.000435f, -0.000242f, -0.000257f, 0.000009f, -0.000493f, -0.000062f, -0.000367f, 0.000344f, 0.001391f, -0.000523f, -0.000189f, 0.000198f, 0.000680f, 0.000058f, -0.000450f, 0.000454f, -0.000182f, -0.000599f, -0.000218f, -0.000272f, 0.000120f, -0.000092f, -0.000487f, -0.000365f, 0.000905f, 0.000297f, 0.000645f, -0.000428f, -0.000055f, -0.000818f, 0.000173f, 0.000876f, 0.000386f, -0.000348f, -0.000290f, -0.000267f, -0.000396f, -0.000593f, -0.000299f, -0.000680f, 0.000004f, 0.000018f, -0.000116f, -0.000049f, -0.000020f, 0.000133f, 0.000353f, 0.000534f, 0.000234f, 0.000007f, -0.000156f, 0.000569f, -0.000498f, 0.000553f, -0.000185f, -0.000587f, -0.000316f, -0.000666f, 0.000205f, 0.000638f, 0.000646f, -0.000302f, -0.000177f, -0.000209f, 0.000031f, -0.000245f, 0.001668f, 0.001287f, -0.000226f, 0.000413f, 0.000251f, 0.000036f, -0.000103f, 0.000108f, -0.000420f, -0.000375f, -0.000609f, -0.000345f, -0.000636f, 0.001014f, -0.000305f, 0.000324f, 0.000889f, 0.000490f, 0.000264f, -0.000241f, -0.000031f, -0.000078f, 0.000331f, -0.000673f, -0.000054f, -0.000033f, -0.000211f, -0.000311f, -0.000079f, 0.000537f, -0.000212f, -0.000274f, -0.000750f, 0.000043f, 0.000072f, 0.000488f, -0.000361f, -0.000377f, 0.000194f, -0.000450f, 0.000678f, -0.000552f, -0.000416f, -0.000537f, -0.000297f, -0.000064f, 0.000599f, 0.000892f, 0.000551f, 0.001017f, -0.000013f, 0.000832f, -0.000209f, 0.000987f, 0.000376f, 0.000060f, -0.000777f, -0.000106f, 0.000658f, -0.000241f, 0.001722f, 0.000609f, -0.000287f, -0.000188f, 0.000818f, -0.001002f, 0.000390f, 0.001138f, -0.000182f, -0.000115f, -0.000431f, 0.000017f, 0.000448f, 0.000809f, 0.000735f, -0.000439f, -0.000097f, -0.000069f, -0.000114f, -0.000518f, -0.000270f, -0.000337f, 0.000305f, 0.000360f, 0.000073f, -0.000221f, 0.000622f, -0.000188f, 0.000450f, 0.000592f, -0.000118f, -0.000219f, 0.000753f, 0.001242f, -0.000482f, 0.001137f, 0.000639f, 0.000362f, 0.000408f, 0.000514f, -0.000657f, -0.000683f, -0.000040f, -0.000448f, 0.001117f, 0.000594f, -0.000060f, 0.000113f, 0.001277f, 0.001141f, -0.000306f, -0.000774f, 0.000724f, 0.000263f, 0.001244f, -0.000619f, 0.000863f, -0.000222f, -0.000846f, -0.000246f, -0.000743f, -0.000353f, -0.000885f, -0.000438f, 0.000426f, 0.000237f, 0.001131f, 0.000360f, -0.000006f, 0.000830f, 0.000090f, -0.000185f, -0.000357f, 0.000930f, -0.000524f, -0.000819f, -0.000433f, 0.000948f, -0.000200f, 0.000209f, -0.000937f, -0.000491f, -0.000659f, 0.000330f, -0.000768f, 0.000268f, -0.000360f, -0.000589f, -0.000212f, 0.000336f, 0.000137f, -0.000223f, 0.001391f, 0.000105f, 0.000142f, 0.000434f, -0.000714f, -0.000547f, 0.000045f, 0.000236f, -0.000869f, 0.000204f, -0.000352f, -0.000205f, 0.000243f, -0.000971f, -0.000301f, -0.000379f, 0.000896f, -0.000124f, 0.000516f, 0.000048f, -0.000770f, -0.000832f, 0.000051f, 0.000490f, -0.000493f, 0.000275f, 0.000604f, -0.000021f, -0.000163f, -0.000104f, -0.000900f, 0.002016f, 0.002103f, 0.001344f, -0.000800f, -0.000289f, -0.000132f, -0.000016f, -0.000916f, -0.000422f, 0.000223f, 0.000862f, -0.000777f, 0.000777f, 0.000816f, -0.000697f, -0.000004f, -0.000173f, 0.000255f, 0.000188f, -0.000683f, 0.001721f, 0.001355f, 0.000465f, -0.000707f, 0.000833f, 0.002282f, -0.000975f, -0.000572f, 0.000023f, -0.000194f, -0.000978f, 0.000256f, -0.000205f, 0.000225f, -0.000576f, -0.000522f, -0.000317f, 0.000516f, -0.000012f, -0.001058f, 0.000693f, 0.001355f, -0.001280f, 0.001020f, 0.000719f, -0.000191f, -0.000728f, -0.000509f, 0.000329f, -0.000554f, -0.000036f, -0.000066f, -0.000712f, 0.002145f, -0.000438f, -0.000341f, -0.000995f, -0.000254f, -0.000928f, 0.001885f, 0.001911f, -0.000340f, -0.001246f, -0.000037f, -0.000566f, -0.000270f, 0.001115f, -0.000466f, -0.000603f, -0.000621f, -0.000501f, -0.000522f, -0.000900f, -0.000141f, -0.000067f, 0.001813f, -0.000430f, -0.000020f, -0.000820f, 0.000524f, 0.000086f, 0.000184f, -0.000555f, -0.000905f, -0.000206f, -0.000515f, 0.000134f, -0.000919f, -0.000568f, -0.000103f, 0.000265f, 0.000428f, -0.000768f, 0.000277f, -0.000460f, 0.000264f, -0.000235f, 0.000304f, -0.000571f, 0.000013f, -0.000041f, -0.000784f, -0.000270f, 0.000118f, -0.000653f, -0.000139f, 0.000977f, 0.000088f, -0.000948f, 0.000214f, 0.000385f, 0.000331f, 0.000464f, 0.000213f, -0.000613f, 0.001668f, -0.000430f, -0.000721f, 0.001500f, 0.000281f, -0.000590f, 0.001580f, 0.002169f, -0.001362f, -0.000170f, -0.000362f, -0.000370f, -0.000901f, 0.000218f, -0.000693f, 0.001009f, 0.001669f, -0.000736f, -0.000047f, 0.000281f, 0.000097f, -0.001318f, 0.000663f, -0.000052f, -0.000512f, -0.000599f, -0.000561f, 0.000469f, -0.000989f, -0.000018f, -0.000135f, 0.000820f, 0.000183f, -0.000071f, -0.000311f, -0.000747f, 0.000723f, -0.000254f, 0.000743f, 0.000691f, -0.000577f, -0.000603f, -0.000321f, -0.000710f, -0.000083f, 0.000071f, -0.000493f, -0.000652f, -0.000649f, -0.000899f, -0.000384f, -0.000440f, 0.000292f, 0.001056f, 0.001268f, -0.000185f, 0.000377f, -0.000580f, -0.000603f, -0.000268f, -0.000593f, -0.000121f, -0.000931f, 0.000112f, -0.000489f, -0.000512f, -0.000406f, -0.000448f, 0.000444f, 0.000899f, -0.000324f, -0.000440f, -0.000089f, -0.000514f, -0.000396f, -0.000611f, 0.000647f, 0.000111f, -0.000321f, -0.000631f, 0.000924f, 0.000502f, -0.000454f, 0.000327f, 0.002048f, 0.000755f, -0.000826f, -0.000491f, -0.000777f, 0.000014f, -0.000659f, 0.000413f, 0.000708f, 0.000179f, -0.000060f, 0.000950f, 0.001374f, -0.000473f, -0.000251f, -0.000006f, 0.000253f, -0.000059f, -0.000156f, 0.000147f, 0.000044f, -0.000307f, -0.000896f, 0.000871f, -0.000131f, 0.000039f, 0.001194f, 0.001604f, -0.000551f, -0.000311f, -0.000274f, -0.000186f, 0.001171f, 0.000372f, -0.000263f, -0.000110f, 0.000432f, -0.000315f, -0.000477f, -0.000582f, -0.000546f, 0.001651f, 0.000454f, 0.000564f, -0.000829f, -0.000624f, -0.000804f, -0.000343f, -0.000672f, -0.000516f, -0.000333f, -0.000504f, 0.000268f, -0.000067f, -0.000659f, -0.000074f, -0.000442f, 0.000571f, -0.000463f, -0.000397f, -0.000040f, -0.000711f, -0.000254f, 0.002055f, 0.000953f, -0.000842f, -0.000079f, 0.000499f, -0.000900f, -0.000329f, 0.000219f, -0.000497f, -0.000634f, -0.000256f, -0.000652f, 0.001244f, -0.000138f, 0.000653f, 0.000213f, 0.000024f, -0.000272f, -0.000429f, 0.000986f, -0.000629f, -0.000086f, 0.000540f, -0.000606f, 0.001001f, 0.000868f, 0.000802f, 0.000668f, -0.000062f, -0.000382f, -0.000460f, -0.000140f, 0.000581f, -0.000296f, 0.000623f, -0.000816f, 0.000259f, -0.000237f, -0.000243f, -0.000414f, 0.000954f, 0.000260f, -0.000011f, 0.000551f, -0.000690f, -0.000496f, 0.000657f, 0.000176f, -0.000597f, 0.000154f, -0.000118f, -0.000322f, -0.000442f, -0.000877f, 0.000382f, 0.002025f, 0.000123f, 0.001128f, 0.000039f, -0.000456f, 0.000230f, 0.000883f, -0.000454f, 0.000298f, -0.000146f, 0.000386f, -0.000616f, 0.000821f, -0.000256f, -0.000380f, 0.001004f, -0.000492f, 0.000454f, -0.000816f, -0.000249f, -0.000416f, 0.000284f, -0.000191f, -0.000398f, -0.000057f, 0.001466f, -0.000215f, -0.000167f, -0.000390f, 0.000273f, -0.000171f, -0.000595f, 0.000191f, 0.000198f, -0.000719f, -0.000109f, 0.000750f, 0.000485f, -0.000932f, 0.000558f, -0.000540f, 0.000829f, 0.000447f, -0.000033f, -0.000189f, -0.000232f, 0.000030f, -0.000427f, -0.000169f, -0.000182f, -0.000021f, 0.000285f, -0.000342f, -0.000215f, 0.000216f, -0.000445f, -0.000528f, -0.000266f, 0.000534f, -0.000538f, 0.000129f, 0.000363f, 0.001206f, -0.000586f, 0.000471f, 0.000354f, -0.000817f, -0.000634f, -0.000676f, 0.000625f, -0.000092f, 0.000847f, -0.000100f, -0.000296f, 0.000106f, 0.000113f, 0.000836f, 0.000152f, 0.000256f, -0.000250f, -0.000590f, 0.000421f, 0.000074f, -0.000457f, 0.000728f, 0.000434f, 0.000459f, -0.000030f, -0.000050f, -0.000260f, -0.000281f, -0.000465f, 0.000223f, 0.000522f, -0.000254f, -0.000548f, 0.000356f, -0.000360f, 0.000094f, 0.000604f, 0.001097f, -0.000026f, -0.000307f, -0.000124f, -0.000068f, -0.000256f, -0.000261f, -0.000068f, 0.000743f, -0.000117f, -0.000018f, -0.000025f, -0.000522f, -0.000709f, 0.000900f, -0.000279f, 0.000531f, 0.000947f, -0.001011f, -0.000424f, -0.000110f, -0.000641f, -0.000798f, -0.000357f, -0.000280f, -0.000844f, -0.000316f, -0.000459f, -0.000806f, -0.000268f, -0.000667f, -0.000436f, 0.000313f, -0.000305f, 0.000012f, -0.000608f, -0.000155f, 0.000026f, 0.000529f, -0.000203f, 0.000159f, 0.000133f, -0.000524f, -0.000043f, -0.000635f, -0.000500f, -0.000335f, -0.000313f, 0.000374f, 0.001341f, 0.000569f, 0.000017f, 0.000049f, -0.000215f, -0.000881f, 0.000232f, 0.000013f, -0.000647f, 0.000124f, -0.000043f, -0.000493f, 0.000240f, -0.000192f, 0.000205f, -0.000267f, 0.000760f, -0.000845f, 0.001898f, 0.000273f, -0.000572f, -0.000652f, 0.000252f, 0.000120f, -0.000967f, -0.000028f, -0.000547f, 0.000729f, 0.000439f, -0.000466f, 0.000317f, -0.000361f, -0.000285f, -0.000437f, -0.000193f, -0.000282f, -0.000350f, 0.000982f, -0.000625f, -0.000264f, -0.000514f, -0.000389f, 0.001303f, 0.000416f, -0.000431f, -0.000006f, -0.000501f, -0.000549f, -0.000222f, -0.000415f, -0.000347f, -0.000374f, -0.000512f, 0.000334f, -0.000622f, -0.000103f, -0.000146f, -0.000233f, -0.000170f, -0.000254f, -0.000601f, -0.000302f, 0.000118f, -0.000464f, -0.000140f, 0.000186f, -0.000107f, -0.000414f, -0.000523f, 0.000007f, -0.000299f, -0.000283f, -0.000108f, -0.000445f, 0.001047f, 0.000079f, -0.000178f, 0.000102f, -0.000335f, -0.000487f, -0.000306f, -0.000062f, -0.000165f, -0.000243f, 0.000003f, -0.000147f, 0.000369f, 0.000385f, -0.000326f, 0.000569f, 0.000715f, -0.000410f, 0.001133f, 0.001092f, -0.000535f, -0.000298f, 0.000500f, -0.000239f, -0.000154f, -0.000208f, -0.000159f, 0.000090f, -0.000141f, 0.000267f, 0.000206f, 0.000160f, -0.000619f, 0.000575f, 0.000051f, -0.000271f, -0.000560f, 0.001740f, 0.000401f, -0.000670f, -0.000095f, 0.000441f, 0.001425f, 0.000509f, 0.000138f, 0.000367f, -0.000110f, 0.000284f, 0.000447f, 0.000791f, -0.000133f, -0.000422f, 0.000725f, 0.000114f, -0.000241f, -0.000550f, 0.000255f, 0.000147f, -0.000412f, -0.000079f, 0.000078f, -0.000652f, -0.000383f, -0.000335f, -0.000008f, -0.000230f, 0.000890f, 0.000068f, -0.000146f, 0.000238f, 0.000067f, -0.000315f, -0.000115f, -0.000467f, 0.000016f, -0.000310f, -0.000591f, -0.000200f, 0.000292f, -0.000309f, -0.000015f, -0.000833f, 0.000749f, 0.001618f, -0.000159f, -0.000125f, 0.000730f, 0.000303f, -0.000564f, 0.000146f, 0.000270f, -0.000329f, 0.000857f, 0.000741f, 0.000766f, 0.000014f, 0.000126f, -0.000261f, -0.000554f, -0.000294f, -0.000217f, 0.001258f, 0.000176f, 0.000375f, 0.000244f, 0.000298f, -0.000612f, 0.000169f, 0.000082f, 0.000635f, 0.000251f, 0.001070f, -0.000538f, 0.000777f, 0.000849f, 0.001636f, 0.000513f, -0.000081f, -0.000039f, -0.000054f, -0.000286f, 0.000525f, 0.000635f, 0.000293f, -0.000405f, -0.000139f, 0.000155f, 0.000582f, 0.000360f, -0.000500f, 0.000042f, -0.000473f, 0.000427f, -0.000221f, 0.000095f, -0.000273f, -0.000300f, 0.000434f, -0.000499f, 0.000114f, 0.000728f, 0.000616f, -0.000516f, 0.000354f, 0.001287f, -0.000154f, -0.000187f, -0.000063f, -0.000044f, -0.000379f, 0.000498f, 0.000407f, 0.000511f, 0.001296f, -0.000742f, -0.000310f, 0.000634f, -0.000256f, 0.000058f, -0.000131f, 0.000083f, -0.000353f, -0.000443f, -0.000160f, 0.000515f, 0.000434f, 0.000569f, 0.001169f, 0.000287f, 0.000302f, 0.001386f, -0.000593f, -0.000297f, -0.000076f, 0.000174f, -0.000069f, 0.001057f, -0.000039f, -0.000494f, -0.000581f, -0.000256f, 0.000040f, 0.000507f, 0.000094f, -0.000166f, -0.000733f, -0.000303f, 0.000429f, 0.001681f, 0.000121f, -0.000901f, 0.000087f, -0.000506f, -0.000928f, 0.000182f, 0.000051f, 0.000135f, 0.000006f, -0.000291f, -0.000481f, -0.000261f, -0.000822f, -0.000438f, 0.000358f, -0.000317f, 0.001202f, -0.000358f, -0.000593f, 0.000371f, -0.000496f, -0.000329f, 0.000882f, -0.000627f, 0.000183f, 0.000367f, -0.000972f, -0.000084f, 0.001710f, -0.000020f, 0.000295f, 0.000347f, -0.000449f, -0.000484f, 0.000343f, 0.000178f, 0.000290f, 0.000846f, -0.000450f, -0.000199f, 0.000306f, -0.000539f, -0.000071f, -0.000219f, -0.000792f, 0.000467f, -0.000689f, 0.000318f, 0.000191f, -0.000700f, -0.000704f, -0.000119f, -0.000806f, -0.000055f, 0.000629f, -0.000715f, 0.000379f, -0.000291f, -0.000545f, -0.000383f, 0.000494f, -0.000242f, -0.000537f, -0.000881f, -0.000645f, -0.000398f, 0.000666f, -0.000803f, 0.000585f, -0.000420f, -0.000816f, -0.000139f, -0.000468f, -0.000319f, -0.000315f, -0.000366f, 0.000025f, 0.000258f, 0.000863f, -0.000194f, 0.000608f, -0.000139f, 0.001113f, 0.000699f, -0.000819f, -0.000290f, -0.000492f, -0.000132f, 0.000226f, 0.001232f, -0.000251f, 0.000290f, 0.001001f, -0.000602f, 0.000226f, 0.000421f, -0.000431f, 0.000203f, -0.000576f, -0.000878f, -0.000071f, 0.000978f, -0.000271f, -0.000025f, -0.000351f, -0.000906f, 0.000768f, -0.000205f, -0.000245f, -0.000081f, -0.001161f, 0.000130f, -0.000546f, 0.000019f, 0.000518f, -0.000089f, -0.000641f, 0.000178f, 0.000097f, -0.000587f, 0.000774f, -0.000733f, 0.000271f, -0.000251f, -0.000305f, -0.000536f, -0.000082f, -0.000676f, -0.000412f, -0.000301f, -0.000760f, -0.000197f, 0.000084f, 0.000523f, 0.000390f, -0.000539f, -0.000606f, -0.000546f, 0.000115f, -0.000427f, 0.000290f, 0.000238f, -0.000222f, -0.000119f, -0.000894f, -0.000295f, 0.000773f, 0.001598f, -0.000071f, 0.000888f, 0.000273f, -0.000289f, -0.000288f, 0.000508f, 0.000851f, 0.000909f, 0.000211f, -0.000310f, -0.000197f, -0.000595f, -0.000369f, -0.000252f, -0.000029f, 0.000326f, 0.000971f, 0.000303f, -0.000712f, -0.000356f, 0.000506f, -0.000251f, -0.000688f, 0.000484f, -0.000266f, -0.000367f, -0.000427f, -0.000584f, 0.000277f, 0.000261f, 0.001416f, -0.000251f, -0.000023f, -0.000945f, -0.000260f, 0.000139f, -0.000052f, 0.000049f, 0.000949f, -0.000931f, -0.000215f, -0.000109f, 0.001179f, 0.000942f, -0.000025f, -0.000779f, -0.000255f, -0.000148f, 0.000201f, 0.000338f, -0.000253f, 0.000982f, 0.000106f, -0.000207f, 0.000096f, 0.000247f, -0.001028f, 0.000431f, 0.000105f, -0.000514f, 0.000045f, -0.000608f, -0.000219f, -0.000801f, -0.000074f, -0.000605f, 0.000583f, 0.000456f, 0.000404f, -0.000747f, 0.000143f, 0.000654f, -0.000540f, 0.000131f, -0.000137f, -0.000131f, -0.000530f, -0.000599f, -0.000181f, 0.000037f, -0.000427f, -0.000603f, 0.000152f, -0.000575f, -0.000210f, 0.000949f, 0.000193f, -0.000394f, 0.000838f, 0.000357f, -0.000566f, -0.000322f, -0.000385f, -0.000316f, 0.000993f, 0.000013f, -0.000406f, -0.000067f, -0.000447f, 0.000098f, 0.000080f, -0.000004f, 0.000907f, 0.000971f, -0.000790f, 0.000094f, 0.000381f, -0.000203f, 0.001738f, -0.000121f, -0.000181f, -0.000047f, -0.000437f, 0.000192f, 0.001051f, -0.000539f, 0.000405f, 0.000552f, -0.000417f, -0.000152f, 0.000486f, -0.000291f, 0.000056f, -0.000281f, -0.000457f, -0.000129f, 0.000628f, -0.000170f, 0.000450f, -0.000494f, 0.000210f, 0.000660f, -0.000347f, -0.000671f, 0.000010f, 0.000092f, -0.000034f, -0.000217f, 0.000498f, 0.000623f, -0.000499f, -0.000130f, 0.000086f, -0.000689f, -0.000266f, -0.000505f, 0.000134f, -0.000860f, 0.000160f, 0.000099f, 0.000283f, -0.000074f, 0.000833f, 0.000200f, -0.000576f, -0.000532f, 0.000251f, 0.000334f, -0.000733f, 0.000265f, -0.000204f, -0.000071f, 0.000161f, 0.000084f, 0.001338f, 0.000739f, 0.000558f, 0.000088f, 0.000301f, -0.000478f, 0.000461f, -0.000047f, -0.000037f, 0.001229f, -0.000622f, -0.000437f, 0.000404f, 0.000338f, -0.000636f, 0.000075f, 0.000169f, 0.000154f, 0.001554f, -0.000330f, 0.000502f, 0.000710f, -0.000097f, 0.001656f, -0.000351f, -0.000207f, 0.000652f, -0.000316f, -0.000665f, 0.000740f, -0.000265f, -0.000050f, -0.000054f, -0.000394f, -0.000395f, -0.000239f, -0.001037f, 0.000743f, -0.000256f, -0.000124f, 0.000508f, 0.000485f, -0.000393f, -0.000274f, -0.000307f, 0.000138f, -0.000059f, -0.000658f, 0.000884f, -0.000269f, -0.000716f, 0.000010f, -0.000448f, 0.001618f, 0.000001f, -0.000014f, -0.000021f, -0.000382f, 0.000744f, 0.000510f, -0.000444f, -0.000443f, 0.000074f, 0.000406f, 0.000070f, -0.000667f, -0.000306f, 0.000230f, -0.000423f, -0.000615f, -0.000320f, 0.000284f, 0.000033f, -0.000130f, -0.000441f, -0.000670f, -0.000397f, -0.000655f, -0.000643f, -0.000402f, 0.000327f, 0.000044f, -0.000391f, 0.000219f, -0.000701f, -0.000107f, -0.000384f, 0.000570f, 0.000191f, -0.000744f, 0.000783f, -0.000027f, 0.000112f, 0.000823f, 0.000408f, 0.000198f, 0.000971f, -0.000974f, 0.000673f, 0.000045f, 0.000280f, -0.000356f, 0.000150f, -0.000147f, -0.000058f, -0.000370f, 0.000198f, -0.000347f, -0.000101f, 0.001217f, 0.000093f, -0.000595f, -0.000389f, -0.000500f, -0.000484f, -0.000237f, 0.000062f, 0.001093f, 0.000671f, -0.000101f, 0.000954f, -0.000082f, 0.000055f, 0.000173f, 0.000113f, 0.000335f, -0.000214f, -0.000662f, 0.000122f, -0.000173f, -0.000402f, 0.000424f, 0.000449f, -0.000043f, 0.000225f, 0.000467f, -0.000126f, -0.000331f, 0.000103f, 0.000751f, 0.000516f, -0.000547f, 0.000149f, -0.000729f, -0.000340f, -0.000866f, -0.000464f, -0.000321f, 0.001129f, -0.000330f, 0.000350f, 0.000157f, -0.000208f, 0.000454f, 0.000044f, -0.000610f, 0.001172f, 0.001700f, -0.000859f, -0.000284f, 0.000997f, 0.000177f, -0.000357f, -0.000571f, 0.000044f, -0.000674f, -0.000254f, -0.000670f, -0.000355f, -0.000701f, 0.000489f, -0.000244f, -0.000200f, 0.000350f, 0.000750f, -0.000222f, -0.000393f, -0.000144f, 0.000958f, 0.001035f, -0.000694f, -0.000034f, -0.000324f, 0.000136f, -0.000355f, -0.000476f, -0.000888f, 0.000779f, 0.000047f, 0.000764f, 0.001139f, 0.000472f, -0.000065f, 0.000319f, -0.000484f, -0.000092f, -0.000714f, -0.000013f, -0.000347f, -0.000612f, -0.000621f, 0.001288f, 0.000103f, -0.000468f, -0.000103f, 0.000063f, -0.000108f, -0.000346f, 0.000070f, -0.000045f, -0.000162f, -0.000123f, -0.000208f, -0.000247f, 0.000421f, -0.000037f, -0.000264f, 0.000072f, 0.000330f, -0.000092f, -0.000356f, -0.000789f, 0.000278f, 0.000360f, -0.000482f, -0.000164f, -0.000400f, -0.000457f, 0.000075f, -0.000143f, -0.000044f, 0.000948f, -0.000689f, -0.000080f, 0.000029f, -0.000656f, 0.000646f, -0.000021f, -0.000194f, -0.000557f, -0.000251f, -0.000052f, 0.000255f, 0.000270f, -0.000033f, 0.000233f, 0.000765f, -0.000157f, -0.000251f, -0.000239f, -0.000366f, 0.000432f, -0.000339f, 0.000007f, -0.000115f, 0.000471f, -0.000411f, 0.000432f, 0.001564f, -0.000156f, -0.000701f, 0.000218f, -0.000916f, 0.000363f, -0.000274f, -0.000226f, -0.000420f, -0.000471f, -0.000499f, 0.000091f, -0.000184f, 0.000834f, 0.001133f, -0.000639f, -0.000653f, 0.000253f, 0.000119f, -0.000256f, 0.000473f, -0.000841f, -0.000265f, -0.000022f, 0.001556f, -0.000411f, -0.000263f, 0.000214f, 0.000253f, -0.000204f, -0.000652f, -0.000180f, -0.000267f, -0.000092f, 0.000061f, -0.000602f, 0.000364f, -0.000143f, 0.000061f, 0.000694f, -0.000792f, -0.000080f, -0.000527f, 0.000210f, -0.000013f, -0.000245f, 0.000053f, 0.001484f, -0.000019f, 0.001096f, 0.000226f, 0.000676f, 0.000772f, 0.000600f, -0.000057f, -0.000590f, -0.000057f, 0.000040f, 0.000045f, -0.000637f, -0.000002f, -0.000584f, -0.000575f, -0.000635f, -0.000019f, 0.000745f, -0.000281f, -0.000075f, 0.000571f, 0.000440f, 0.000128f, 0.000174f, 0.000394f, -0.000129f, 0.000292f, -0.000546f, -0.000541f, 0.000487f, 0.000214f, -0.000080f, -0.000417f, -0.000483f, -0.000484f, -0.000374f, 0.000048f, 0.000273f, -0.000117f, 0.000274f, -0.000251f, 0.000422f, -0.000275f, -0.000387f, -0.000335f, 0.000209f, -0.000000f, 0.000385f, -0.000188f, -0.000472f, -0.000145f, -0.000244f, -0.000197f, -0.000035f, 0.000141f, 0.000446f, 0.000599f, -0.000114f, -0.000712f, 0.000540f, -0.000009f, 0.000220f, -0.000625f, -0.000348f, -0.000098f, 0.000648f, -0.000048f, 0.000687f, -0.000580f, 0.000229f, 0.001037f, -0.000504f, 0.000748f, -0.000323f, 0.000275f, -0.000354f, -0.000484f, -0.000374f, -0.000038f, -0.000262f, 0.000742f, 0.000413f, -0.000363f, 0.000088f, -0.000112f, -0.000064f, -0.000039f, -0.000419f, -0.000554f, -0.000306f, -0.000523f, -0.000104f, -0.000080f, 0.000522f, 0.000195f, -0.000034f, 0.000301f, 0.000464f, 0.000146f, 0.001241f, 0.000470f, 0.000459f, -0.000292f, -0.000584f, 0.000177f, 0.000376f, -0.000439f, 0.000294f, -0.000002f, -0.000285f, -0.000094f, -0.000469f, 0.000007f, 0.000793f, -0.000242f, 0.000751f, 0.000503f, -0.000069f, 0.000073f, -0.000067f, 0.000077f, -0.000729f, -0.000090f, -0.000189f, -0.000575f, -0.000077f, -0.000224f, -0.000562f, -0.000039f, -0.000037f, 0.000187f, 0.000462f, -0.000554f, 0.000964f, -0.000198f, -0.000172f, 0.000039f, -0.000423f, -0.000423f, -0.000254f, 0.000405f, 0.000030f, 0.000250f, -0.000260f, -0.000285f, -0.000486f, 0.000054f, 0.000686f, 0.000244f, 0.000113f, 0.000131f, 0.000342f, -0.000652f, 0.000448f, -0.000624f, 0.000224f, -0.000776f, -0.000279f, -0.000003f, 0.000044f, -0.000421f, 0.000130f, 0.001021f, 0.000364f, 0.001364f, -0.000207f, -0.000225f, -0.000122f, 0.000229f, -0.000013f, -0.000500f, -0.000528f, -0.000301f, -0.000297f, -0.000382f, 0.000238f, 0.000620f, 0.000524f, -0.000389f, 0.000968f, 0.000520f, -0.000824f, 0.000569f, 0.000598f, -0.000390f, 0.000175f, 0.000487f, -0.000028f, -0.000210f, -0.000020f, 0.000039f, 0.000391f, 0.000096f, 0.000384f, -0.000670f, -0.000087f, 0.000134f, 0.000606f, 0.000262f, -0.000263f, 0.000248f, 0.000175f, -0.000744f, 0.000505f, 0.000566f, 0.000136f, -0.000162f, 0.000060f, 0.000630f, 0.000485f, 0.000212f, -0.000264f, 0.000248f, -0.000311f, -0.000191f, 0.000803f, -0.000553f, -0.000730f, -0.000111f, -0.000547f, 0.000182f, -0.000192f, 0.000163f, -0.000268f, -0.000808f, 0.000303f, -0.000128f, -0.000538f, 0.000123f, 0.000536f, 0.000210f, 0.000584f, 0.000062f, -0.000141f, -0.000421f, -0.000091f, -0.000461f, -0.000674f, 0.000400f, 0.000514f, -0.000498f, 0.000467f, 0.000011f, -0.000298f, 0.000235f, -0.000075f, 0.000060f, -0.000163f, -0.000383f, -0.000182f, -0.000509f, 0.000235f, 0.000532f, 0.000371f, -0.000774f, -0.000176f, -0.000560f, 0.000005f, -0.000083f, 0.000539f, 0.000282f, -0.000356f, 0.000095f, 0.000398f, -0.000144f, 0.000136f, -0.000471f, 0.000709f, 0.000195f, -0.000578f, 0.000811f, 0.000309f, -0.000299f, -0.000270f, -0.000696f, 0.000382f, 0.000641f, 0.001407f, 0.001041f, 0.000101f, -0.000033f, 0.000181f, -0.000474f, -0.000635f, -0.000464f, -0.000592f, 0.000194f, 0.000249f, -0.000240f, 0.000243f, 0.000277f, -0.000223f, 0.000235f, 0.000533f, -0.000527f, -0.000395f, -0.000017f, -0.000003f, -0.000536f, -0.000209f, 0.000111f, -0.000321f, -0.000149f, 0.000015f, -0.000321f, 0.000121f, 0.000620f, -0.000328f, -0.000345f, -0.000130f, 0.000516f, -0.000329f, 0.001104f, -0.000372f, -0.000317f, 0.000045f, 0.000195f, -0.000324f, 0.000578f, 0.000055f, 0.000163f, -0.000271f, 0.000064f, 0.000421f, -0.000154f, 0.000216f, 0.001059f, -0.000661f, -0.000234f, 0.000249f, -0.000707f, 0.000370f, 0.000267f, 0.001359f, 0.000373f, -0.000139f, 0.000277f, -0.000022f, -0.000284f, 0.000851f, -0.000426f, -0.000411f, -0.000129f, -0.000426f, 0.000115f, -0.000280f, 0.000094f, -0.000148f, -0.000572f, -0.000560f, -0.000484f, -0.000773f, 0.000270f, -0.000570f, -0.000143f, -0.000269f, 0.000251f, -0.000545f, 0.000302f, -0.000279f, -0.000257f, 0.000368f, -0.000600f, 0.000059f, -0.000184f, -0.000634f, 0.000212f, -0.000343f, -0.000109f, 0.000924f, 0.000017f, 0.000160f, 0.000369f, -0.000013f, -0.000360f, -0.000627f, 0.000406f, -0.000364f, 0.001370f, 0.001169f, -0.000548f, -0.000337f, 0.000794f, 0.000189f, -0.000862f, -0.000583f, 0.000521f, 0.000923f, -0.000566f, 0.000826f, -0.001104f, 0.000156f, 0.000303f, -0.000381f, -0.000454f, 0.000119f, -0.001145f, 0.000597f, -0.000372f, 0.000437f, 0.000994f, 0.000198f, 0.000728f, -0.000355f, -0.000516f, -0.000652f, -0.000003f, -0.000534f, -0.000272f, -0.000682f, -0.000058f, -0.000312f, -0.000460f, 0.000442f, -0.000203f, -0.000285f, -0.000164f, -0.000307f, 0.000235f, -0.000205f, 0.000189f, 0.000093f, -0.000258f, -0.000143f, -0.000297f, -0.000233f, -0.000100f, 0.000093f, -0.000433f, -0.000083f, -0.000333f, -0.000030f, 0.000471f, 0.000164f, -0.000852f, 0.000384f, 0.000126f, -0.000532f, -0.000613f, -0.000071f, -0.000626f, 0.000033f, -0.000255f, 0.000148f, 0.000128f, 0.000761f, 0.000434f, 0.000486f, 0.000669f, -0.000083f, -0.000315f, 0.000243f, 0.000670f, -0.000392f, -0.000297f, 0.000110f, 0.000421f, -0.000432f, -0.000096f, -0.000445f, 0.000868f, -0.000400f, -0.000138f, -0.000057f, 0.000021f, -0.000121f, 0.000311f, 0.000569f, 0.000716f, -0.000055f, -0.000255f, -0.000041f, -0.000558f, -0.000487f, -0.000044f, 0.000848f, -0.000630f, -0.000394f, 0.000062f, 0.001336f, 0.000881f, -0.000303f, -0.000306f, 0.000266f, -0.000094f, 0.000157f, 0.000343f, 0.000140f, 0.000245f, -0.000454f, -0.000522f, 0.000148f, 0.000327f, 0.000420f, 0.000139f, -0.000310f, 0.000121f, -0.000413f, -0.000132f, 0.000101f, -0.000101f, 0.001047f, 0.000225f, -0.000603f, -0.000571f, 0.000189f, -0.000598f, 0.000652f, -0.000015f, 0.000423f, -0.000383f, 0.000262f, -0.000111f, -0.000609f, -0.000144f, 0.000736f, -0.000368f, 0.001013f, -0.000168f, -0.000644f, -0.000415f, 0.000219f, 0.000173f, -0.000028f, 0.000751f, -0.000096f, -0.000541f, -0.000192f, -0.000021f, -0.000639f, 0.000815f, -0.000751f, -0.000203f, 0.001215f, -0.000002f, 0.000456f, 0.000131f, -0.000225f, -0.000066f, -0.000433f, 0.000662f, 0.000534f, -0.001024f, 0.000768f, -0.000429f, 0.000343f, 0.000472f, -0.000296f, 0.000384f, -0.000042f, -0.000502f, 0.000045f, -0.000428f, -0.000204f, -0.000168f, 0.000759f, 0.000998f, -0.000408f, -0.000239f, -0.000284f, -0.000135f, 0.000047f, -0.000769f, -0.000013f, -0.000044f, 0.000784f, 0.000477f, 0.000172f, -0.000536f, 0.000140f, 0.000417f, -0.000684f, 0.000763f, -0.000255f, 0.000050f, -0.000001f, -0.000322f, -0.000153f, -0.000152f, 0.000076f, 0.000938f, -0.000350f, 0.000501f, 0.000034f, 0.001292f, -0.000431f, -0.000307f, -0.000324f, -0.000512f, 0.000078f, -0.000414f, -0.000775f, 0.000080f, 0.000033f, -0.000302f, 0.000455f, -0.000335f, -0.000710f, 0.000195f, -0.000069f, -0.000172f, 0.000634f, 0.001175f, -0.000340f, -0.000051f, -0.000043f, -0.000407f, 0.000076f, 0.000486f, 0.000449f, -0.000226f, -0.000272f, -0.000046f, 0.000075f, 0.000308f, 0.000166f, -0.000615f, 0.000149f, -0.000260f, 0.000431f, 0.000268f, 0.000734f, -0.000062f, 0.000975f, 0.000904f, -0.000489f, -0.000173f, -0.000030f, -0.000677f, -0.000540f, 0.000108f, 0.000832f, 0.001211f, -0.000685f, 0.000112f, -0.000155f, 0.000197f, -0.000382f, -0.000079f, -0.000036f, 0.000870f, -0.000315f, 0.000205f, 0.000589f, 0.000368f, -0.000291f, 0.000398f, 0.000254f, 0.000959f, -0.000287f, -0.000694f, -0.000214f, -0.000763f, -0.000314f, -0.000270f, -0.000016f, 0.000305f, -0.000697f, 0.000027f, 0.000296f, -0.000108f, 0.000390f, -0.000491f, 0.000280f, 0.000045f, -0.000800f, 0.000086f, 0.000529f, -0.000463f, 0.000770f, -0.000448f, -0.000247f, 0.000047f, 0.000642f, -0.000049f, -0.000476f, 0.000108f, -0.000350f, -0.000700f, 0.000342f, 0.000332f, 0.000436f, 0.000479f, -0.000658f, -0.000108f, -0.000127f, -0.000130f, -0.000472f, -0.000067f, 0.000302f, -0.000112f, -0.000403f, 0.000531f, -0.000734f, 0.000121f, -0.000057f, 0.000119f, -0.000066f, -0.000465f, 0.000102f, 0.001198f, -0.000064f, 0.000493f, -0.000364f, -0.000291f, -0.000108f, -0.000489f, 0.000448f, 0.000415f, -0.000295f, 0.000282f, -0.000212f, 0.001115f, -0.000308f, 0.000163f, 0.000541f, -0.000790f, 0.000153f, -0.000130f, -0.000285f, 0.000326f, -0.000470f, -0.000517f, -0.000104f, -0.000731f, 0.000411f, -0.000418f, 0.000081f, 0.000081f, -0.000529f, 0.000220f, 0.000184f, -0.000533f, 0.000013f, -0.000284f, -0.000736f, -0.000568f, -0.000580f, 0.000181f, -0.000598f, -0.000018f, -0.000442f, -0.000552f, -0.000382f, 0.000064f, -0.000255f, 0.000063f, -0.000010f, -0.000097f, -0.000350f, -0.000492f, -0.000461f, -0.000213f, 0.001591f, 0.000800f, -0.000362f, -0.000261f, -0.000352f, -0.000454f, 0.001113f, -0.000490f, -0.000032f, -0.000493f, -0.000285f, -0.000502f, 0.000599f, -0.000062f, -0.000121f, 0.000299f, 0.000463f, -0.000143f, 0.000468f, 0.000641f, -0.000285f, -0.000231f, -0.000091f, 0.000034f, -0.000288f, 0.000080f, 0.000215f, 0.000142f, 0.000267f, 0.000392f, -0.000740f, -0.000001f, -0.000541f, 0.000505f, 0.000714f, -0.000561f, -0.000083f, -0.000396f, 0.000251f, -0.000292f, -0.000220f, -0.000230f, -0.000426f, 0.000582f, 0.000467f, -0.000149f, 0.000236f, 0.000294f, 0.000131f, -0.000329f, 0.000052f, 0.000210f, -0.000161f, -0.000315f, -0.000174f, 0.000095f, -0.000841f, 0.000288f, 0.000311f, -0.000670f, -0.000443f, 0.000177f, -0.000222f, -0.000371f, 0.000436f, 0.001169f, -0.000370f, -0.000095f, 0.000084f, 0.000345f, 0.000310f, 0.000009f, -0.000597f, 0.000078f, -0.000246f, -0.000330f, -0.000015f, 0.000004f, -0.000117f, 0.000327f, 0.000107f, -0.000608f, 0.000451f, -0.000034f, -0.000313f, -0.000588f, -0.000056f, -0.000130f, 0.000487f, 0.000691f, -0.000551f, 0.000725f, 0.000276f, -0.000169f, -0.000276f, -0.000265f, 0.001258f, 0.000943f, -0.000282f, -0.000030f, -0.000636f, 0.000241f, 0.000289f, 0.000101f, 0.000305f, -0.000594f, 0.000353f, 0.000257f, -0.000239f, 0.000054f, 0.000278f, 0.000785f, 0.000260f, -0.000423f, -0.000143f, -0.000259f, 0.000690f, 0.000002f, 0.000216f, 0.000246f, 0.000321f, 0.000505f, -0.000459f, -0.000537f, -0.000071f, 0.000028f, -0.000166f, -0.000584f, 0.000338f, 0.000964f, -0.000568f, 0.000233f, -0.000511f, -0.000131f, -0.000706f, -0.000214f, -0.000096f, 0.000411f, -0.000684f, 0.000276f, 0.000168f, 0.000615f, -0.000672f, 0.001345f, 0.000269f, -0.000390f, -0.000390f, 0.000128f, 0.000361f, 0.001130f, -0.000664f, -0.000379f, -0.000506f, 0.000207f, 0.000283f, 0.001590f, 0.000443f, 0.000086f, -0.000382f, 0.000509f, -0.000327f, -0.000283f, -0.000420f, -0.000206f, 0.000012f, -0.000321f, -0.000268f, 0.001103f, -0.000279f, -0.000064f, -0.000720f, -0.000108f, -0.000809f, 0.000027f, -0.000445f, 0.000322f, 0.000352f, 0.000466f, -0.000823f, 0.000433f, 0.000058f, 0.000219f, 0.000409f, -0.000556f, 0.000360f, -0.000173f, 0.000604f, 0.000513f, 0.000067f, 0.000012f, -0.000827f, 0.000013f, 0.000221f, -0.000165f, 0.000293f, -0.000489f, -0.000169f, -0.000195f, -0.000309f, 0.000674f, -0.000158f, 0.000216f, 0.000275f, -0.000383f, 0.000224f, -0.000401f, 0.000885f, -0.000042f, -0.000252f, -0.000239f, -0.000397f, 0.000438f, -0.000161f, -0.000162f, -0.000341f, 0.000270f, -0.000374f, 0.000086f, 0.000461f, 0.000075f, 0.000466f, -0.000073f, -0.000087f, 0.000018f, 0.000349f, 0.000259f, -0.000252f, -0.000205f, 0.000334f, -0.000592f, -0.000294f, 0.000387f, -0.000102f, -0.000691f, -0.000499f, -0.000398f, 0.000213f, 0.000154f, 0.000020f, -0.000168f, 0.000522f, -0.000469f, -0.000042f, 0.000452f, 0.000618f, -0.000338f, -0.000185f, -0.000343f, -0.000324f, -0.000211f, -0.000076f, 0.000287f, 0.000032f, -0.000394f, -0.000345f, -0.000158f, 0.000827f, -0.000111f, 0.000021f, -0.000094f, -0.000032f, 0.000336f, -0.000376f, 0.000169f, -0.000244f, -0.000242f, 0.000593f, -0.000387f, -0.000403f, 0.000074f, -0.000150f, 0.000650f, -0.000549f, 0.000205f, -0.000011f, -0.000352f, 0.000969f, -0.000030f, 0.000048f, -0.000300f, -0.000119f, -0.000020f, -0.000254f, -0.000487f, -0.000296f, -0.000651f, -0.000413f, 0.000157f, -0.000155f, -0.000505f, -0.000083f, 0.000282f, -0.000290f, 0.000236f, 0.000458f, -0.000182f, -0.000360f, -0.000505f, -0.000275f, -0.000256f, 0.000212f, -0.000023f, 0.000933f, -0.000144f, -0.000468f, -0.000384f, -0.000044f, 0.000078f, -0.000158f, 0.000264f, -0.000605f, -0.000317f, -0.000325f, 0.000101f, 0.000929f, 0.000899f, -0.000352f, 0.000350f, 0.000253f, -0.000446f, 0.000134f, 0.000552f, -0.000334f, -0.000010f, -0.000511f, -0.000116f, -0.000490f, -0.000403f, 0.000170f, -0.000277f, 0.000479f, -0.000271f, 0.000216f, 0.000953f, -0.000640f, -0.000181f, 0.000232f, 0.000206f, -0.000191f, -0.000366f, 0.000215f, -0.000444f, 0.000101f, -0.000139f, -0.000281f, 0.000601f, 0.000620f, -0.000827f, 0.000153f, 0.000277f, -0.000400f, 0.000289f, 0.000036f, -0.000039f, -0.000355f, -0.000452f, 0.000139f, 0.000208f, 0.000046f, -0.000185f, -0.000052f, -0.000130f, -0.000327f, 0.000165f, 0.000733f, 0.000360f, -0.000259f, 0.000779f, 0.000625f, -0.000117f, -0.000110f, -0.000184f, -0.000259f, -0.000148f, 0.000680f, -0.000379f, 0.000789f, 0.000700f, -0.000277f, -0.000149f, 0.000517f, 0.000200f, 0.000785f, 0.001274f, -0.000646f, -0.000100f, 0.000352f, -0.000162f, -0.000313f, 0.000630f, -0.000163f, -0.000614f, 0.000033f, -0.000110f, 0.000578f, 0.000315f, 0.000033f, 0.000251f, 0.000410f, 0.000700f, 0.000622f, 0.000036f, 0.000605f, 0.000185f, -0.000461f, -0.000187f, -0.000342f, -0.000046f, -0.000403f, -0.000371f, 0.000149f, 0.000044f, -0.000384f, 0.000909f, 0.000079f, 0.000062f, -0.000382f, -0.000154f, 0.000637f, 0.000468f, 0.000090f, -0.000363f, -0.000444f, -0.000041f, -0.000392f, -0.000076f, 0.000304f, 0.000203f, 0.000170f, 0.001053f, 0.000070f, -0.000429f, -0.000243f, 0.000844f, -0.000122f, -0.000149f, 0.000468f, -0.000649f, -0.000231f, 0.000637f, 0.000597f, 0.001193f, 0.000790f, 0.000398f, -0.000063f, 0.000193f, -0.000232f, 0.000343f, 0.000511f, -0.000255f, -0.000216f, -0.000005f, -0.000242f, 0.000448f, 0.000378f, -0.000122f, 0.000235f, -0.000168f, -0.000408f, -0.000124f, 0.000796f, 0.000971f, -0.000320f, -0.000500f, -0.000699f, 0.000167f, 0.000455f, 0.000097f, -0.000626f, 0.000345f, -0.000157f, -0.000785f, 0.000254f, -0.000153f, 0.000594f, 0.000364f, -0.000505f, 0.000128f, -0.000438f, -0.000119f, 0.000214f, 0.000275f, 0.000487f, -0.000873f, 0.000632f, 0.000136f, 0.000229f, 0.001120f, -0.000808f, 0.000425f, -0.000081f, 0.000377f, 0.001051f, -0.000088f, 0.000615f, -0.000570f, 0.000262f, 0.001038f, -0.000110f, -0.000236f, 0.000051f, 0.000130f, -0.000417f, -0.000156f, 0.000949f, 0.000090f, -0.000180f, 0.000263f, -0.000864f, -0.000156f, 0.000122f, -0.000272f, 0.000116f, -0.000127f, -0.000252f, -0.000779f, 0.000196f, -0.000439f, -0.000507f, -0.000141f, -0.000038f, -0.000531f, -0.000196f, -0.000076f, -0.000151f, -0.000205f, -0.000365f, -0.000336f, 0.000137f, -0.000275f, -0.000014f, 0.000022f, -0.000154f, 0.000745f, -0.000955f, 0.000019f, 0.000352f, -0.000358f, -0.000104f, 0.000383f, -0.000088f, -0.000328f, 0.000141f, 0.000907f, -0.000351f, 0.000094f, -0.000963f, -0.000488f, 0.000189f, -0.000750f, 0.000216f, -0.000160f, -0.000196f, 0.000375f, -0.000517f, 0.000796f, -0.000367f, -0.000151f, 0.000174f, -0.000136f, 0.000259f, -0.000509f, -0.000104f, -0.000351f, -0.000192f, 0.001182f, -0.000874f, 0.000638f, 0.000771f, -0.000370f, 0.000086f, 0.000666f, -0.000083f, -0.000211f, -0.000559f, -0.000333f, -0.000057f, -0.000198f, -0.000891f, 0.000362f, 0.000599f, -0.000479f, 0.000176f, 0.000632f, 0.000768f, 0.000086f, 0.000547f, 0.001021f, -0.000066f, -0.000359f, -0.000253f, 0.000334f, -0.000733f, -0.000349f, -0.000495f, -0.000306f, 0.000084f, 0.000485f, 0.000100f, 0.000621f, -0.000249f, -0.000730f, 0.000413f, 0.000721f, -0.000445f, -0.000201f, 0.000443f, -0.000347f, -0.000486f, 0.000340f, 0.000705f, -0.000088f, -0.000734f, -0.000458f, -0.000532f, 0.000026f, -0.000469f, 0.000265f, 0.000889f, -0.000827f, -0.000054f, 0.000097f, 0.000191f, 0.000251f, -0.000407f, -0.000161f, -0.000284f, -0.000208f, 0.000117f, 0.000948f, -0.000014f, -0.000710f, 0.000139f, -0.000063f, 0.000438f, -0.000136f, -0.000608f, 0.000270f, -0.000335f, -0.000213f, -0.000080f, -0.000287f, -0.000595f, -0.000540f, -0.000073f, -0.000101f, 0.001368f, 0.000085f, 0.000433f, 0.000847f, -0.000621f, -0.000114f, -0.000164f, -0.000112f, -0.000315f, -0.000068f, -0.000264f, -0.000300f, 0.000465f, -0.000079f, 0.000129f, -0.000071f, -0.000340f, 0.000022f, -0.000315f, 0.000681f, -0.000161f, 0.000114f, 0.000171f, -0.000488f, 0.000212f, -0.000333f, 0.000079f, 0.000143f, -0.000345f, -0.000314f, -0.000667f, 0.000184f, 0.001227f, 0.000079f, -0.000071f, -0.000827f, 0.000766f, 0.000587f, -0.000179f, 0.000049f, -0.000785f, -0.000236f, -0.000246f, 0.000739f, 0.000341f, -0.000521f, -0.000535f, -0.000575f, -0.000132f, 0.000582f, -0.000472f, 0.000371f, -0.000270f, -0.000339f, -0.000446f, -0.000568f, 0.000134f, 0.000126f, -0.000065f, -0.000486f, -0.000439f, 0.000341f, -0.000160f, 0.000394f, 0.000174f, 0.000040f, -0.000039f, 0.000997f, 0.000495f, 0.000148f, -0.000141f, -0.000591f, 0.001110f, 0.000163f, -0.000401f, 0.000453f, -0.000318f, -0.000297f, 0.000511f, 0.000598f, 0.000162f, -0.000266f, 0.000806f, -0.000230f, 0.000036f, 0.000052f, -0.000177f, 0.001021f, -0.000080f, -0.000013f, 0.000415f, 0.000310f, -0.000206f, 0.000649f, -0.000287f, -0.000444f, 0.000739f, -0.000152f, -0.000084f, -0.000347f, -0.000540f, 0.000449f, -0.000294f, -0.000348f, 0.000592f, -0.000124f, -0.000443f, -0.000261f, -0.000123f, -0.000190f, -0.000234f, 0.000269f, -0.000326f, -0.000078f, -0.000619f, -0.000638f, 0.000196f, 0.000641f, 0.000141f, -0.000454f, -0.000061f, -0.000593f, 0.000929f, 0.001183f, -0.000142f, -0.000035f, 0.000032f, -0.000164f, -0.000058f, 0.000359f, -0.000608f, -0.000358f, 0.000124f, -0.000520f, 0.000132f, -0.000008f, 0.000142f, 0.000155f, 0.000524f, -0.000257f, 0.000301f, 0.000477f, -0.000742f, -0.000288f, -0.000342f, -0.000402f, -0.000422f, -0.000139f, 0.000616f, 0.000066f, -0.000260f, -0.000442f, -0.000392f, 0.000540f, 0.000048f, -0.000346f, -0.000300f, -0.000095f, 0.000811f, -0.000158f, 0.000253f, -0.000194f, -0.000066f, -0.000540f, 0.000255f, 0.000649f, 0.000179f, -0.000080f, -0.000378f, -0.000040f, 0.000163f, -0.000054f, 0.001694f, 0.000737f, -0.000727f, 0.000509f, -0.000393f, 0.000213f, -0.000022f, 0.000446f, -0.000107f, -0.000393f, -0.000527f, 0.000120f, 0.000069f, -0.000263f, -0.000272f, -0.000585f, 0.000267f, 0.000032f, -0.000335f, 0.000003f, -0.000177f, -0.000333f, 0.000106f, 0.000365f, -0.000394f, 0.000332f, 0.000195f, -0.000445f, -0.000295f, -0.000304f, 0.000424f, 0.000596f, -0.000315f, -0.000327f, 0.000009f, -0.000238f, -0.000136f, -0.000336f, 0.000262f, -0.000334f, -0.000221f, 0.000302f, -0.000237f, 0.000498f, 0.000125f, -0.000353f, -0.000098f, 0.000379f, -0.000358f, 0.000376f, 0.000087f, 0.000185f, 0.000597f, -0.000666f, -0.000482f, -0.000604f, -0.000064f, -0.000086f, -0.000448f, 0.000301f, 0.000084f, 0.000229f, -0.000310f, 0.000252f, 0.000660f, -0.000227f, 0.000494f, -0.000196f, -0.000265f, 0.000799f, 0.000141f, -0.000154f, 0.000214f, -0.000051f, -0.000140f, 0.000661f, 0.000152f, -0.000230f, 0.000980f, 0.000322f, 0.000217f, -0.000223f, -0.000321f, -0.000409f, 0.000287f, -0.000270f, -0.000640f, 0.001023f, 0.000046f, -0.000699f, 0.000452f, -0.000065f, -0.000222f, -0.000398f, 0.000157f, 0.000144f, -0.000065f, -0.000473f, -0.000516f, -0.000465f, -0.000082f, -0.000321f, 0.001055f, 0.001009f, -0.000401f, -0.000183f, 0.000325f, -0.000436f, 0.000288f, -0.000029f, -0.000332f, 0.000079f, -0.000400f, 0.000246f, -0.000072f, 0.000920f, -0.000561f, 0.000371f, 0.000060f, -0.000259f, -0.000172f, 0.000040f, -0.000253f, -0.000810f, -0.000098f, -0.000234f, -0.000521f, 0.000553f, -0.000082f, 0.001021f, -0.000312f, -0.000075f, 0.000515f, -0.000499f, -0.000399f, -0.000563f, -0.000269f, 0.000065f, 0.000533f, 0.000539f, 0.000182f, -0.000428f, -0.000149f, 0.000494f, 0.000688f, -0.000302f, 0.000070f, -0.000445f, -0.000115f, 0.000143f, 0.000416f, 0.000388f, -0.000560f, -0.000357f, -0.000529f, 0.000127f, -0.000037f, 0.000107f, 0.000808f, -0.000178f, 0.000085f, -0.000404f, -0.000097f, 0.000373f, -0.000031f, -0.000006f, -0.000435f, -0.000175f, -0.000163f, 0.000129f, -0.000213f, 0.000043f, 0.000055f, 0.000024f, -0.000503f, 0.000337f, 0.000752f, -0.000284f, -0.000036f, -0.000375f, -0.000467f, 0.000236f, -0.000298f, 0.000359f, -0.000290f, 0.000530f, 0.000110f, -0.000424f, -0.000207f, -0.000261f, 0.000260f, -0.000223f, 0.000097f, 0.000525f, 0.000899f, 0.000304f, -0.000686f, -0.000101f, -0.000131f, 0.000137f, -0.000078f, 0.000796f, 0.000451f, -0.001054f, -0.000003f, -0.000147f, 0.000622f, 0.000158f, 0.000608f, -0.000506f, 0.000361f, 0.001218f, 0.000421f, -0.000666f, -0.000541f, -0.000099f, 0.000300f, -0.000172f, -0.000485f, -0.000303f, -0.000403f, 0.000591f, -0.000115f, -0.000158f, -0.000068f, 0.000145f, 0.000248f, 0.000023f, -0.000021f, 0.000012f, -0.000617f, -0.000285f, -0.000083f, -0.000287f, 0.000084f, -0.000550f, -0.000002f, 0.000671f, -0.000509f, 0.000057f, -0.000558f, 0.000223f, 0.000177f, -0.000542f, -0.000243f, -0.000491f, 0.000207f, -0.000719f, 0.000504f, -0.000196f, -0.000268f, 0.000083f, -0.000313f, 0.000118f, -0.000515f, 0.001045f, -0.000121f, -0.000324f, -0.000340f, -0.000340f, 0.000384f, -0.000474f, -0.000462f, -0.000353f, 0.000876f, 0.000191f, -0.000061f, -0.000089f, -0.000491f, 0.000065f, 0.000470f, -0.000465f, -0.000269f, -0.000417f, 0.000114f, -0.000378f, -0.000269f, -0.000347f, -0.000139f, 0.000091f, -0.000284f, 0.000134f, 0.000305f, 0.000499f, 0.000383f, -0.000004f, -0.000438f, -0.000271f, 0.000267f, -0.000274f, 0.000934f, 0.000211f, -0.000777f, -0.000080f, 0.000010f, -0.000441f, -0.000165f, -0.000058f, -0.000033f, 0.000294f, 0.000348f, -0.000107f, -0.000049f, -0.000331f, 0.000083f, -0.000060f, 0.000124f, 0.000688f, 0.000382f, -0.000073f, -0.000212f, -0.000245f, -0.000011f, 0.000458f, -0.000199f, 0.000137f, 0.000463f, 0.000190f, 0.000243f, 0.000223f, 0.000307f, -0.000822f, 0.000150f, -0.000337f, -0.000224f, 0.001088f, 0.000069f, -0.000217f, -0.000687f, 0.000285f, -0.000092f, 0.000038f, 0.000198f, -0.000407f, 0.000060f, -0.000645f, 0.000826f, -0.000055f, -0.000387f, -0.000121f, 0.000234f, 0.000506f, 0.000618f, 0.000986f, 0.000077f, 0.000274f, 0.000325f, 0.001306f, 0.000833f, -0.000563f, 0.000039f, 0.000070f, -0.000389f, 0.000141f, 0.000446f, -0.000288f, 0.000658f, -0.000160f, -0.000415f, -0.000789f, -0.000097f, 0.000677f, 0.000045f, -0.000083f, -0.000045f, -0.000123f, 0.000311f, 0.000329f, 0.000208f, -0.000464f, -0.000155f, 0.000291f, -0.000451f, -0.000352f, -0.000186f, 0.000156f, -0.000425f, -0.000087f, -0.000571f, -0.000299f, 0.000038f, 0.000017f, 0.000150f, -0.000335f, 0.000241f, 0.000510f, -0.000391f, -0.000345f, 0.000251f, 0.000090f, -0.000229f, 0.000426f, -0.000435f, -0.000343f, -0.000095f, 0.000552f, 0.000693f, -0.000088f, 0.000701f, -0.000750f, 0.000093f, 0.000248f, -0.000209f, 0.000714f, -0.000638f, 0.000255f, -0.000300f, -0.000488f, -0.000550f, -0.000252f, -0.000251f, -0.000240f, 0.000544f, -0.000224f, 0.000895f, -0.000475f, 0.000636f, 0.000027f, 0.000429f, 0.001345f, 0.000028f, -0.000272f, -0.000243f, -0.000124f, -0.000597f, -0.000044f, -0.000552f, -0.000423f, 0.000062f, 0.000623f, 0.001184f, -0.000292f, -0.000008f, 0.000361f, -0.000202f, 0.000125f, 0.000142f, -0.000143f, -0.000247f, -0.000365f, -0.000182f, 0.000220f, 0.000983f, -0.000150f, 0.000166f, -0.000311f, 0.000472f, -0.000387f, 0.000618f, -0.000357f, -0.000631f, -0.000646f, -0.000105f, 0.000720f, -0.000554f, -0.000200f, -0.000144f, -0.000330f, -0.000554f, 0.000868f, 0.000654f, -0.000154f, -0.000452f, -0.000437f, 0.001156f, -0.000634f, 0.000333f, 0.000347f, -0.000673f, -0.000314f, -0.000361f, 0.000742f, 0.000162f, 0.000012f, 0.000447f, -0.000343f, 0.000239f, -0.000377f, 0.000022f, -0.000238f, 0.000380f, -0.000229f, -0.000504f, -0.000067f, 0.000248f, 0.000686f, 0.000036f, -0.000347f, 0.000536f, -0.000151f, -0.000328f, 0.000429f, 0.000334f, 0.000467f, 0.000472f, -0.000268f, -0.000063f, 0.000329f, -0.000446f, -0.000175f, 0.000673f, -0.000342f, 0.000380f, 0.000210f, -0.000227f, -0.000199f, 0.000225f, 0.000061f, -0.000399f, -0.000144f, -0.000480f, 0.000447f, 0.000559f, -0.000141f, 0.000171f, -0.000435f, 0.000054f, 0.000161f, -0.000181f, 0.000215f, -0.000439f, 0.000769f, -0.000149f, 0.000408f, 0.000282f, 0.000024f, -0.000222f, 0.000142f, -0.000260f, -0.000111f, -0.000328f, -0.000607f, 0.000043f, -0.000204f, -0.000396f, 0.000281f, -0.000268f, -0.000327f, -0.000119f, 0.000515f, -0.000707f, -0.000086f, -0.000250f, 0.000267f, -0.000073f, -0.000749f, 0.000146f, -0.000456f, 0.000125f, 0.000267f, -0.000679f, -0.000093f, -0.000115f, -0.000484f, 0.000330f, -0.000356f, -0.000448f, -0.000416f, -0.000087f, -0.000365f, -0.000285f, -0.000279f, 0.000275f, 0.000018f, 0.000817f, -0.000762f, -0.000262f, -0.000310f, 0.000346f, 0.000007f, -0.000069f, 0.000052f, 0.000102f, 0.000454f, -0.000396f, 0.000163f, -0.000027f, -0.000330f, 0.000004f, 0.000062f, 0.000212f, -0.000449f, 0.000445f, -0.000012f, 0.000156f, -0.000213f, 0.000083f, 0.000602f, 0.000014f, 0.000172f, -0.000687f, 0.000026f, -0.000082f, 0.000095f, 0.000101f, 0.000093f, 0.000386f, -0.000357f, -0.000311f, 0.000087f, 0.000700f, -0.000283f, -0.000039f, -0.000099f, -0.000071f, 0.000672f, -0.000299f, -0.000122f, -0.000384f, 0.001148f, -0.000411f, 0.000206f, -0.000463f, 0.000095f, 0.000775f, -0.000290f, 0.000441f, 0.000167f, 0.000237f, 0.000355f, -0.000314f, 0.000697f, 0.000038f, -0.000151f, 0.000104f, 0.000602f, -0.000223f, 0.000010f, 0.000669f, -0.000317f, 0.000247f, -0.000706f, 0.000607f, 0.000064f, 0.000217f, -0.000134f, 0.000170f, -0.000514f, 0.000098f, 0.000513f, -0.000212f, 0.000841f, 0.000622f, 0.000279f, -0.000134f, 0.000053f, 0.000019f, -0.000215f, -0.000319f, -0.000129f, 0.000293f, 0.000027f, -0.000011f, 0.000542f, -0.000177f, -0.000501f, -0.000167f, -0.000192f, -0.000469f, -0.000344f, -0.000083f, 0.000042f, -0.000628f, 0.000686f, 0.000191f, 0.000634f, 0.000266f, 0.000909f, 0.000422f, -0.000141f, 0.000165f, 0.000444f, -0.000066f, -0.000714f, 0.000198f, 0.000310f, -0.000139f, -0.000211f, -0.000016f, 0.000526f, -0.000511f, -0.000053f, -0.000044f, -0.000585f, -0.000637f, 0.000107f, 0.000437f, 0.000033f, -0.000105f, -0.000170f, -0.000181f, -0.000148f, -0.000391f, 0.000739f, 0.000633f, 0.000047f, -0.000639f, 0.000007f, -0.000480f, 0.000094f, 0.000272f, -0.000677f, -0.000597f, -0.000146f, 0.000613f, 0.001333f, 0.000160f, -0.000262f, -0.000434f, -0.000568f, 0.000352f, 0.000384f, 0.000347f, -0.000230f, -0.000143f, -0.000469f, 0.000034f, 0.000098f, 0.000783f, 0.000181f, -0.000258f, 0.000268f, -0.000698f, -0.000082f, -0.000022f, 0.000129f, -0.000129f, -0.000242f, 0.001264f, -0.000424f, 0.000335f, -0.000378f, -0.000154f, -0.000321f, -0.000655f, -0.000217f, 0.000066f, -0.000584f, -0.000221f, -0.000307f, -0.000259f, -0.000261f, -0.000062f, 0.000487f, 0.000976f, 0.000489f, 0.000852f, 0.000181f, 0.000097f, -0.000157f, 0.000357f, 0.000021f, 0.000318f, -0.000324f, 0.000084f, -0.000377f, -0.000426f, -0.000124f, 0.000481f, -0.000587f, -0.000075f, 0.000356f, -0.000095f, -0.000349f, 0.000211f, -0.000006f, -0.000278f, -0.000874f, -0.000013f, -0.000486f, 0.000564f, -0.000451f, 0.000480f, 0.000216f, -0.000159f, -0.000174f, -0.000241f, 0.000223f, -0.000137f, 0.000207f, -0.000069f, -0.000492f, -0.000239f, 0.000353f, -0.000201f, 0.000664f, 0.000260f, -0.000407f, -0.000580f, 0.000184f, 0.000479f, 0.000058f, -0.000260f, -0.000175f, 0.000072f, -0.000276f, 0.000765f, 0.001094f, -0.000386f, -0.000092f, 0.000131f, -0.000127f, -0.000161f, -0.000156f, -0.000102f, -0.000179f, 0.000013f, 0.000287f, -0.000193f, -0.000461f, -0.000350f, -0.000015f, 0.000202f, -0.000465f, 0.000595f, -0.000029f, 0.000003f, -0.000512f, 0.000172f, 0.000291f, -0.000554f, 0.000159f, -0.000046f, -0.000098f, 0.000091f, 0.000328f, -0.000256f, -0.000122f, -0.000890f, 0.000849f, 0.000814f, 0.000154f, -0.000406f, -0.000354f, -0.000581f, 0.000150f, 0.000545f, 0.000381f, -0.000461f, 0.000512f, -0.000167f, 0.000039f, -0.000363f, -0.000065f, -0.000241f, -0.000377f, 0.000212f, 0.000325f, 0.000195f, -0.000212f, -0.000768f, -0.000162f, 0.000144f, -0.000366f, -0.000264f, -0.000131f, -0.000257f, -0.000583f, -0.000001f, -0.000021f, -0.000233f, -0.000167f, -0.000581f, 0.000022f, -0.000109f, 0.000214f, 0.000035f, 0.000207f, -0.000061f, -0.000189f, -0.000094f, 0.000762f, 0.000333f, -0.000085f, -0.000447f, 0.000270f, -0.000034f, 0.000013f, 0.000396f, -0.000521f, -0.000143f, 0.000143f, 0.000405f, -0.000138f, -0.000357f, 0.000078f, 0.000014f, -0.000596f, -0.000147f, -0.000234f, -0.000077f, 0.000022f, 0.000607f, 0.000248f, 0.000373f, -0.000225f, 0.000442f, 0.000014f, -0.000348f, 0.000778f, -0.000030f, 0.000013f, -0.000066f, 0.000211f, -0.000025f, -0.000139f, -0.000731f, 0.000427f, -0.000110f, -0.000445f, -0.000103f, -0.000340f, 0.000337f, 0.000232f, 0.000862f, 0.000479f, 0.000069f, -0.000480f, -0.000228f, -0.000146f, -0.000160f, 0.000704f, 0.000410f, 0.000511f, 0.000595f, 0.000176f, 0.000216f, 0.000426f, -0.000410f, -0.000245f, 0.000494f, 0.001497f, -0.000333f, 0.000678f, -0.000169f, -0.000496f, 0.000158f, -0.000393f, -0.000010f, 0.000107f, 0.000162f, 0.000279f, 0.000075f, -0.000252f, -0.000187f, 0.000811f, -0.000263f, -0.000266f, -0.000408f, -0.000154f, -0.000154f, -0.000579f, -0.000388f, -0.000154f, 0.000047f, 0.000760f, -0.000097f, -0.000538f, -0.000288f, -0.000041f, 0.000297f, -0.000867f, -0.000329f, -0.000292f, 0.000069f, 0.000727f, 0.000466f, 0.000351f, -0.000065f, -0.000338f, 0.000435f, 0.000072f, -0.000135f, -0.000111f, -0.000380f, -0.000313f, -0.000367f, -0.000523f, 0.000059f, -0.000531f, -0.000275f, -0.000724f, 0.000077f, -0.000638f, 0.000408f, 0.000074f, 0.000280f, 0.000080f, -0.000143f, -0.000066f, -0.000197f, 0.000854f, -0.000825f, -0.000011f, -0.000545f, 0.000448f, 0.000152f, 0.000067f, -0.000074f, 0.000571f, -0.000324f, 0.000222f, -0.000551f, -0.000255f, 0.000584f, 0.000356f, -0.000029f, 0.000408f, 0.000254f, 0.000340f, 0.000796f, 0.000429f, -0.000419f, -0.000156f, -0.000538f, 0.000102f, -0.000294f, 0.000278f, -0.000050f, -0.000086f, -0.000372f, 0.000773f, -0.000294f, -0.000306f, 0.000201f, 0.000130f, -0.000465f, -0.000322f, -0.000220f, 0.000328f, 0.000198f, -0.000142f, 0.000409f, 0.000137f, 0.000126f, -0.000120f, -0.000484f, 0.000293f, 0.000391f, -0.000077f, 0.000168f, -0.000566f, 0.000125f, 0.000334f, -0.000178f, -0.000292f, -0.000038f, -0.000130f, -0.000039f, -0.000448f, 0.000279f, -0.000079f, -0.000458f, -0.000599f, 0.000056f, -0.000136f, -0.000335f, 0.000007f, 0.000043f, 0.000641f, -0.000470f, -0.000174f, 0.000217f, 0.000026f, -0.000046f, 0.000462f, -0.000241f, -0.000465f, 0.000175f, 0.000185f, -0.000022f, -0.000030f, 0.000153f, 0.000300f, -0.000604f, 0.000105f, 0.000075f, -0.000026f, -0.000406f, -0.000087f, -0.000276f, -0.000508f, 0.000180f, 0.000025f, 0.000474f, -0.000468f, -0.000546f, 0.000209f, -0.000210f, -0.000668f, 0.000395f, -0.000321f, -0.000098f, 0.000294f, 0.000071f, 0.000200f, 0.000157f, -0.000092f, -0.000370f, -0.000380f, -0.000499f, -0.000235f, -0.000135f, 0.000171f, -0.000042f, -0.000091f, -0.000277f, -0.000290f, 0.000101f, 0.000089f, -0.000041f, -0.000105f, -0.000488f, 0.000397f, -0.000210f, 0.000843f, 0.000236f, 0.000302f, -0.000283f, 0.000204f, -0.000373f, -0.000227f, -0.000115f, 0.000226f, 0.000030f, -0.000162f, -0.000370f, 0.000281f, -0.000117f, -0.000089f, 0.000260f, -0.000040f, 0.000822f, -0.000084f, 0.000290f, 0.000009f, -0.000177f, 0.000130f, 0.000078f, -0.000589f, -0.000452f, 0.000017f, 0.000492f, -0.000401f, -0.000116f, 0.000163f, 0.000353f, 0.000408f, 0.000698f, 0.000203f, -0.000165f, -0.000038f, 0.000935f, -0.000658f, -0.000225f, -0.000228f, 0.000329f, -0.000147f, 0.000351f, -0.000049f, 0.000836f, -0.000176f, -0.000106f, -0.000119f, -0.000352f, 0.000478f, -0.000313f, 0.000457f, 0.000401f, 0.000223f, 0.000310f, -0.000141f, -0.000130f, 0.000075f, 0.000068f, -0.000246f, 0.000064f, -0.000416f, 0.000390f, 0.000308f, 0.000448f, 0.000553f, 0.000042f, 0.000265f, 0.000307f, -0.000048f, 0.000541f, 0.000971f, -0.000057f, 0.000442f, 0.000131f, 0.000315f, 0.000126f, 0.000024f, 0.000002f, 0.001230f, -0.000149f, 0.000502f, -0.000391f, -0.000103f, 0.000105f, 0.000214f, -0.000059f, 0.000257f, -0.000297f, -0.000221f, 0.000280f, -0.000066f, -0.000097f, -0.000346f, 0.000082f, 0.000104f, -0.000307f, -0.000095f, 0.000188f, 0.000603f, -0.000440f, 0.000173f, 0.000065f, 0.000187f, -0.000245f, -0.000093f, 0.000154f, 0.000011f, -0.000375f, -0.000022f, -0.000156f, 0.000126f, 0.000568f, -0.000057f, 0.000138f, 0.000127f, -0.000010f, 0.000232f, 0.000195f, 0.000542f, -0.000127f, -0.000368f, -0.000609f, 0.000312f, 0.000207f, -0.000085f, -0.000131f, 0.000051f, 0.000078f, -0.000550f, -0.000451f, -0.000550f, -0.000367f, 0.000081f, 0.001005f, -0.000549f, -0.000342f, -0.000253f, 0.000513f, 0.000423f, -0.000077f, 0.000074f, -0.000588f, 0.000350f, 0.000401f, -0.000313f, 0.000198f, 0.000235f, -0.000465f, 0.000205f, 0.000195f, 0.000865f, 0.000092f, -0.000444f, -0.000472f, 0.000553f, 0.000153f, -0.000523f, -0.000081f, -0.000442f, -0.000310f, 0.000259f, 0.000624f, 0.000020f, 0.000106f, -0.000626f, 0.000414f, -0.000336f, 0.000343f, 0.000227f, -0.000567f, 0.000435f, -0.000192f, -0.000057f, 0.000215f, 0.000563f, -0.000507f, 0.000781f, 0.000078f, -0.000248f, -0.000442f, -0.000197f, 0.000239f, 0.000166f, -0.000272f, -0.000185f, -0.000170f, -0.000360f, -0.000076f, -0.000193f, -0.000018f, 0.000541f, 0.000883f, 0.000030f, 0.000210f, 0.001284f, 0.000326f, -0.000777f, -0.000146f, -0.000027f, 0.000480f, 0.000279f, 0.000528f, 0.000041f, -0.000561f, 0.000681f, 0.000389f, 0.000190f, 0.000128f, 0.000910f, 0.000917f, -0.000096f, -0.000340f, -0.000386f, -0.000662f, 0.000673f, 0.001113f, -0.000750f, -0.000521f, -0.000249f, 0.000120f, -0.000882f, 0.000218f, 0.000161f, -0.000574f, 0.000241f, -0.000286f, -0.000542f, -0.000240f, -0.000207f, -0.000322f, 0.000213f, -0.000092f, 0.000279f, -0.000312f, -0.000040f, 0.000357f, -0.000078f, 0.000097f, 0.000246f, 0.000377f, -0.000724f, -0.000510f, -0.000452f, 0.000386f, -0.000173f, -0.000237f, -0.000337f, -0.000323f, 0.000214f, 0.000506f, -0.000219f, -0.000514f, 0.000179f, -0.000070f, -0.000312f, 0.000603f, -0.000458f, -0.000367f, -0.000393f, 0.000626f, -0.000458f, -0.000430f, -0.000038f, 0.000150f, -0.000194f, -0.000302f, -0.000135f, -0.000250f, -0.000346f, 0.000064f, -0.000098f, -0.000113f, -0.000156f, 0.000614f, -0.000054f, 0.000184f, -0.000123f, 0.000334f, -0.000442f, 0.000322f, 0.000817f, -0.000260f, -0.000573f, -0.000376f, -0.000145f, -0.000482f, 0.000438f, -0.000268f, 0.000148f, -0.000189f, 0.000836f, 0.001221f, -0.000255f, -0.000403f, -0.000344f, -0.000365f, -0.000076f, 0.001105f, -0.000562f, -0.000007f, -0.000772f, -0.000461f, 0.000186f, 0.000674f, 0.000843f, -0.000182f, 0.000275f, 0.000320f, -0.000207f, 0.000431f, -0.000283f, -0.000398f, 0.000284f, -0.000330f, 0.000139f, -0.000325f, -0.000560f, -0.000047f, -0.000279f, 0.000162f, 0.000006f, -0.000319f, -0.000232f, -0.000225f, -0.000101f, -0.000070f, 0.000031f, -0.000208f, 0.000029f, -0.000086f, 0.000097f, 0.000763f, -0.000462f, -0.000226f, 0.000485f, -0.000158f, -0.000414f, -0.000042f, 0.000247f, -0.000052f, -0.000138f, 0.000162f, -0.000524f, 0.000364f, 0.000719f, -0.000651f, 0.000082f, -0.000778f, -0.000040f, -0.000133f, 0.000027f, -0.000128f, 0.000859f, -0.000060f, 0.000479f, 0.000880f, -0.000233f, 0.000148f, -0.000330f, 0.000310f, -0.000259f, -0.000333f, -0.000471f, 0.000311f, -0.000280f, -0.000048f, -0.000303f, 0.000155f, 0.000112f, -0.000363f, 0.000335f, -0.000438f, 0.000125f, -0.000183f, -0.000461f, -0.000333f, -0.000315f, 0.000312f, 0.000332f, 0.000155f, -0.000386f, -0.000568f, 0.000420f, -0.000248f, -0.000537f, -0.000387f, 0.000306f, 0.000308f, -0.000531f, -0.000531f, 0.000176f, 0.000113f, -0.000596f, -0.000343f, -0.000017f, -0.000411f, -0.000435f, 0.000288f, -0.000195f, 0.000854f, -0.000177f, -0.000087f, -0.000105f, -0.000268f, -0.000039f, -0.000267f, -0.000610f, 0.000450f, 0.000444f, 0.000101f, -0.000041f, -0.000308f, -0.000381f, -0.000275f, -0.000457f, 0.000317f, -0.000143f, -0.000259f, 0.000562f, -0.000375f, 0.000484f, 0.000542f, -0.000294f, 0.000296f, -0.000158f, 0.000401f, 0.000898f, 0.000246f, -0.000415f, -0.000460f, -0.000124f, 0.000368f, 0.000159f, 0.000074f, -0.000263f, -0.000414f, -0.000151f, -0.000073f, 0.000232f, -0.000103f, -0.000203f, -0.000003f, 0.000344f, 0.000056f, 0.001280f, -0.000042f, -0.000178f, 0.000218f, 0.000320f, 0.000245f, -0.000131f, 0.000284f, 0.000393f, -0.000233f, 0.000588f, -0.000167f, -0.000554f, 0.000126f, -0.000170f, -0.000179f, 0.000380f, -0.000296f, -0.000336f, -0.000003f, -0.000755f, 0.000492f, -0.000280f, 0.000208f, 0.000037f, -0.000227f, 0.000124f, 0.000004f, -0.000393f, -0.000177f, -0.000088f, -0.000057f, 0.000149f, 0.000301f, 0.000857f, -0.000316f, -0.000173f, 0.000489f, -0.000178f, -0.000418f, 0.000087f, -0.000347f, -0.000054f, -0.000272f, 0.000116f, 0.000127f, -0.000012f, 0.000051f, -0.000239f, -0.000019f, 0.000219f, 0.000393f, 0.000092f, 0.000648f, -0.000149f, -0.000210f, -0.000193f, -0.000204f, -0.000290f, -0.000238f, 0.000004f, -0.000317f, 0.000207f, 0.000327f, -0.000075f, -0.000502f, -0.000179f, -0.000419f, 0.000090f, -0.000296f, -0.000128f, 0.000530f, -0.000033f, -0.000291f, 0.000780f, -0.000141f, -0.000149f, -0.000409f, -0.000089f, -0.000352f, 0.000211f, -0.000195f, -0.000038f, -0.000362f, 0.000379f, 0.000046f, -0.000570f, 0.000232f, -0.000052f, -0.000184f, -0.000184f, -0.000329f, -0.000005f, 0.000599f, 0.000049f, 0.000076f, -0.000666f, 0.000391f, 0.000262f, -0.000217f, 0.000130f, 0.000138f, -0.000237f, 0.000527f, 0.000291f, 0.000028f, 0.000682f, -0.000295f, 0.000063f, 0.000605f, 0.000333f, -0.000032f, -0.000177f, -0.000073f, 0.000401f, -0.000247f, 0.000266f, 0.000477f, -0.000300f, -0.000187f, -0.000156f, 0.000185f, -0.000198f, 0.000243f, 0.000211f, 0.000185f, 0.000626f, -0.000195f, -0.000279f, -0.000002f, -0.000307f, 0.000060f, 0.000519f, 0.000245f, -0.000298f, -0.000465f, -0.000486f, 0.000189f, -0.000057f, 0.000064f, -0.000180f, 0.000079f, 0.000020f, 0.000517f, 0.000449f, -0.000192f, -0.000273f, -0.000345f, -0.000023f, 0.000295f, 0.000706f, -0.000475f, -0.000029f, -0.000081f, -0.000212f, -0.000462f, -0.000191f, -0.000480f, -0.000077f, -0.000202f, 0.000140f, 0.000404f, 0.000295f, 0.000272f, 0.000273f, 0.000473f, 0.000229f, -0.000522f, 0.000318f, -0.000122f, 0.000341f, 0.000025f, -0.000746f, 0.000026f, 0.000240f, 0.000019f, 0.000389f, 0.000081f, 0.000262f, -0.000037f, -0.000093f, 0.000238f, 0.000056f, -0.000029f, -0.000239f, 0.000147f, -0.000484f, -0.000386f, -0.000109f, 0.000227f, 0.000097f, -0.000185f, 0.000161f, -0.000555f, -0.000137f, -0.000076f, 0.000372f, -0.000088f, -0.000075f, -0.000585f, 0.000003f, -0.000366f, -0.000165f, -0.000071f, 0.001200f, 0.000167f, 0.000060f, -0.000089f, 0.000417f, 0.000484f, -0.000221f, -0.000328f, 0.000505f, -0.000315f, 0.000333f, 0.000133f, -0.000670f, 0.000297f, -0.000070f, 0.000281f, 0.000096f, -0.000482f, -0.000005f, -0.000048f, 0.000218f, 0.000108f, -0.000238f, -0.000181f, -0.000416f, 0.000145f, -0.000538f, 0.000233f, -0.000145f, 0.000295f, -0.000114f, 0.000334f, -0.000361f, -0.000146f, 0.000109f, 0.000146f, 0.000030f, -0.000258f, 0.000059f, 0.000616f, -0.000547f, 0.000294f, 0.000341f, 0.000113f, -0.000081f, -0.000118f, -0.000129f, 0.000073f, 0.000013f, -0.000200f, -0.000218f, 0.000193f, 0.000036f, -0.000198f, 0.000378f, -0.000278f, -0.000213f, 0.000099f, 0.000405f, -0.000316f, 0.000323f, -0.000132f, -0.000252f, -0.000445f, -0.000270f, -0.000110f, 0.000061f, -0.000345f, -0.000282f, -0.000312f, 0.000072f, -0.000029f, 0.000143f, -0.000272f, 0.000074f, -0.000507f, -0.000308f, -0.000100f, 0.000224f, 0.000503f, 0.000038f, 0.000220f, 0.000262f, 0.000393f, -0.000634f, 0.000156f, 0.000086f, -0.000228f, -0.000278f, 0.000068f, 0.000456f, 0.000320f, 0.000243f, 0.000676f, -0.000375f, -0.000147f, -0.000243f, 0.000106f, 0.000228f, 0.000045f, -0.000279f, 0.000338f, -0.000025f, -0.000181f, 0.000433f, 0.000077f, -0.000392f, 0.000273f, 0.000165f, 0.000084f, -0.000263f, -0.000340f, 0.000102f, -0.000506f, -0.000345f, -0.000049f, -0.000082f, -0.000269f, -0.000296f, -0.000296f, 0.000719f, 0.000466f, 0.000109f, -0.000256f, -0.000276f, 0.000296f, 0.000009f, -0.000343f, -0.000156f, 0.000066f, 0.000687f, 0.000572f, 0.000124f, -0.000548f, -0.000019f, 0.000498f, -0.000275f, 0.000228f, -0.000372f, 0.000077f, 0.000537f, -0.000297f, 0.000128f, -0.000185f, -0.000025f, 0.000468f, 0.000085f, 0.000181f, 0.000289f, -0.000494f, 0.000209f, -0.000278f, 0.000092f, -0.000313f, 0.000042f, 0.000065f, -0.000343f, -0.000405f, -0.000237f, 0.000343f, 0.000217f, -0.000382f, 0.000270f, -0.000114f, -0.000262f, 0.000579f, -0.000041f, 0.000116f, 0.000350f, 0.000118f, 0.000522f, -0.000321f, -0.000066f, -0.000199f, -0.000318f, -0.000541f, 0.000756f, 0.000417f, 0.000451f, -0.000535f, 0.000490f, -0.000009f, -0.000168f, -0.000112f, -0.000407f, -0.000215f, -0.000204f, -0.000379f, -0.000058f, -0.000067f, -0.000361f, -0.000056f, -0.000201f, 0.000193f, 0.000027f, -0.000170f, 0.000069f, 0.000848f, -0.000421f, -0.000283f, -0.000051f, 0.000766f, -0.000465f, 0.000111f, -0.000081f, -0.000500f, -0.000578f, 0.000244f, -0.000007f, -0.000546f, 0.000474f, -0.000008f, 0.000072f, 0.000985f, 0.000337f, 0.000112f, -0.000106f, -0.000028f, 0.000260f, 0.000328f, -0.000265f, -0.000289f, -0.000300f, 0.000055f, -0.000542f, -0.000210f, -0.000200f, -0.000010f, 0.000134f, 0.000054f, 0.000088f, 0.000120f, -0.000409f, 0.000409f, 0.000078f, -0.000676f, 0.000259f, 0.000206f, -0.000415f, 0.000002f, 0.000615f, -0.000392f, -0.000237f, 0.000058f, 0.000230f, -0.000570f, 0.000252f, -0.000137f, 0.000112f, -0.000123f, 0.000642f, 0.000096f, 0.000360f, -0.000488f, -0.000028f, 0.000329f, 0.000116f, 0.000212f, 0.000022f, 0.000609f, 0.000261f, 0.000295f, 0.000064f, -0.000183f, -0.000403f, -0.000234f, 0.000087f, -0.000173f, -0.000113f, 0.000372f, -0.000191f, -0.000201f, -0.000273f, 0.000667f, -0.000242f, -0.000387f, -0.000350f, -0.000276f, -0.000028f, 0.001168f, 0.000363f, 0.000008f, 0.000502f, -0.000747f, 0.000717f, 0.000420f, -0.000272f, -0.000147f, 0.000289f, 0.000440f, 0.000071f, 0.000247f, 0.001197f, -0.000254f, 0.000285f, 0.000200f, -0.000220f, 0.000830f, -0.000148f, -0.000480f, -0.000088f, -0.000095f, 0.000043f, -0.000040f, 0.000225f, 0.000009f, 0.000044f, 0.000508f, 0.000279f, -0.000185f, -0.000061f, 0.000331f, -0.000334f, -0.000316f, 0.000522f, -0.000147f, -0.000531f, 0.000413f, -0.000193f, 0.000076f, -0.000034f, 0.000038f, 0.000219f, -0.000379f, 0.000412f, 0.000233f, 0.000139f, -0.000068f, 0.000116f, 0.000142f, 0.000421f, -0.000005f, -0.000347f, 0.000267f, 0.000014f, -0.000547f, -0.000325f, -0.000058f, -0.000127f, -0.000441f, -0.000168f, -0.000240f, 0.000258f, 0.000407f, -0.000247f, 0.000403f, 0.000199f, 0.000490f, 0.000095f, 0.000011f, -0.000214f, 0.000447f, -0.000323f, 0.000335f, -0.000385f, -0.000364f, 0.000098f, 0.000025f, 0.000242f, 0.000173f, -0.000492f, -0.000275f, 0.000039f, 0.000494f, -0.000126f, -0.000032f, -0.000560f, 0.000276f, 0.000490f, 0.000107f, 0.000056f, -0.000531f, 0.000138f, 0.000475f, 0.000011f, -0.000118f, 0.000968f, 0.000279f, -0.000623f, -0.000071f, 0.000004f, 0.000092f, -0.000515f, 0.000213f, 0.001030f, -0.000005f, 0.000467f, -0.000071f, 0.000219f, -0.000021f, -0.000277f, 0.000490f, -0.000294f, 0.000430f, -0.000767f, 0.000069f, -0.000142f, 0.000258f, 0.000084f, 0.000213f, -0.000396f, 0.000141f, 0.000277f, -0.000057f, -0.000499f, -0.000318f, 0.000156f, -0.000462f, 0.000335f, -0.000334f, 0.000440f, 0.000063f, -0.000047f, 0.000012f, 0.000505f, 0.000111f, 0.000326f, -0.000129f, 0.000245f, 0.000198f, 0.000026f, 0.000874f, 0.000536f, -0.000696f, -0.000334f, 0.000042f, -0.000228f, -0.000637f, 0.000060f, -0.000250f, -0.000240f, 0.000301f, 0.000202f, 0.000177f, -0.000052f, -0.000040f, -0.000019f, -0.000430f, -0.000569f, 0.000475f, 0.000628f, -0.000319f, -0.000292f, -0.000092f, 0.000189f, -0.000658f, -0.000178f, -0.000305f, -0.000077f, -0.000260f, -0.000206f, -0.000045f, 0.000120f, 0.000036f, 0.000121f, -0.000230f, -0.000150f, 0.000020f, 0.000107f, 0.000341f, -0.000606f, 0.000355f, 0.000080f, 0.000053f, 0.000285f, -0.000428f, 0.000175f, 0.000006f, 0.000358f, 0.000234f, -0.000040f, -0.000473f, -0.000265f, -0.000154f, -0.000401f, 0.000240f, 0.000508f, -0.000201f, -0.000213f, -0.000044f, -0.000459f, 0.000083f, 0.000081f, -0.000385f, -0.000274f, 0.000294f, 0.000103f, -0.000326f, 0.000087f, -0.000258f, -0.000324f, 0.000052f, 0.000119f, -0.000145f, -0.000090f, -0.000234f, -0.000600f, 0.000202f, 0.000098f, -0.000068f, 0.000110f, -0.000010f, 0.000132f, -0.000168f, 0.000100f, -0.000565f, -0.000334f, -0.000126f, 0.000798f, 0.001403f, -0.000426f, -0.000283f, 0.000025f, -0.000443f, -0.000041f, -0.000160f, -0.000500f, -0.000067f, -0.000386f, 0.000390f, -0.000187f, -0.000016f, 0.000300f, 0.000160f, 0.000185f, 0.000472f, 0.000378f, 0.000017f, -0.000356f, -0.000674f, 0.000131f, -0.000268f, -0.000276f, 0.000707f, -0.000177f, 0.000107f, -0.000536f, -0.000157f, 0.000201f, 0.000353f, 0.000104f, -0.000160f, 0.000355f, -0.000009f, -0.000672f, 0.000462f, -0.000176f, 0.000031f, -0.000204f, -0.000157f, -0.000260f, 0.000424f, -0.000012f, 0.000262f, 0.000261f, 0.000008f, 0.000093f, 0.000175f, -0.000184f, 0.000288f, -0.000146f, -0.000608f, 0.000461f, -0.000087f, 0.000050f, 0.000120f, -0.000021f, -0.000076f, 0.000272f, -0.000243f, 0.000101f, -0.000165f, -0.000210f, 0.000083f, 0.000017f, 0.000370f, 0.000331f, -0.000251f, -0.000225f, 0.000212f, 0.000651f, -0.000220f, -0.000227f, 0.000123f, 0.000025f, 0.000404f, 0.000492f, 0.000041f, -0.000451f, -0.000186f, -0.000067f, -0.000028f, 0.000195f, 0.000112f, -0.000209f, 0.000079f, -0.000348f, 0.000615f, -0.000081f, -0.000179f, -0.000443f, 0.000335f, 0.000361f, -0.000362f, -0.000276f, -0.000483f, -0.000548f, -0.000333f, -0.000075f, 0.000200f, -0.000388f, -0.000309f, 0.000138f, 0.000149f, 0.000426f, 0.000549f, -0.000069f, -0.000268f, -0.000384f, 0.000879f, -0.000027f, -0.000393f, 0.000117f, -0.000600f, -0.000194f, 0.000465f, 0.000141f, -0.000226f, -0.000517f, -0.000065f, -0.000115f, -0.000234f, -0.000307f, -0.000194f, 0.000324f, 0.000595f, 0.000575f, -0.000489f, -0.000429f, -0.000201f, -0.000100f, -0.000474f, -0.000164f, 0.000575f, 0.000485f, -0.000178f, -0.000029f, -0.000286f, 0.000228f, 0.000698f, -0.000029f, 0.000507f, -0.000185f, -0.000530f, 0.000070f, 0.000304f, -0.000310f, -0.000056f, -0.000290f, 0.000123f, -0.000190f, 0.000248f, 0.000303f, 0.000298f, 0.000551f, -0.000125f, 0.000129f, -0.000044f, 0.000006f, -0.000481f, 0.000285f, -0.000216f, -0.000112f, 0.000084f, -0.000178f, -0.000226f, 0.000160f, 0.000182f, 0.000596f, 0.000521f, -0.000526f, 0.000008f, -0.000216f, 0.000004f, -0.000626f, -0.000080f, 0.000356f, -0.000143f, -0.000285f, -0.000512f, -0.000082f, 0.000201f, -0.000099f, 0.000007f, -0.000608f, 0.000316f, -0.000081f, -0.000124f, 0.000017f, -0.000141f, -0.000168f, 0.000687f, 0.000544f, 0.000123f, 0.000401f, 0.000002f, -0.000257f, -0.000318f, -0.000155f, -0.000190f, 0.000041f, 0.000391f, -0.000001f, -0.000024f, -0.000051f, -0.000412f, -0.000146f, -0.000488f, 0.000019f, -0.000546f, -0.000217f, -0.000384f, -0.000206f, -0.000114f, -0.000078f, 0.000466f, 0.000291f, -0.000214f, 0.000053f, 0.000039f, -0.000096f, -0.000043f, -0.000131f, -0.000273f, 0.000050f, 0.000341f, 0.000237f, -0.000148f, -0.000193f, 0.000722f, 0.000358f, -0.000012f, 0.000205f, -0.000332f, 0.000194f, 0.000640f, -0.000186f, -0.000576f, -0.000214f, -0.000532f, 0.000799f, -0.000152f, 0.000178f, -0.000302f, 0.000599f, 0.000209f, -0.000282f, -0.000106f, 0.000120f, -0.000350f, 0.000116f, 0.000346f, 0.000000f, -0.000324f, 0.000393f, 0.000056f, 0.000205f, 0.000483f, -0.000391f, -0.000369f, 0.000189f, -0.000008f, 0.000169f, 0.000081f, 0.000183f, -0.000239f, 0.000124f, -0.000398f, -0.000221f, -0.000435f, 0.000315f, -0.000276f, -0.000137f, 0.000021f, 0.000626f, 0.000227f, -0.000395f, 0.000168f, 0.000016f, 0.000177f, -0.000113f, 0.000458f, -0.000110f, -0.000156f, 0.000077f, 0.000346f, 0.000208f, -0.000349f, -0.000473f, 0.000232f, 0.000497f, 0.000407f, -0.000112f, 0.000316f, 0.000076f, -0.000141f, 0.000762f, -0.000155f, -0.000631f, 0.000201f, 0.000055f, 0.000588f, -0.000132f, -0.000143f, -0.000187f, -0.000060f, 0.000313f, -0.000010f, 0.000070f, 0.000300f, 0.000187f, 0.000362f, -0.000343f, -0.000127f, -0.000055f, -0.000313f, -0.000394f, -0.000025f, 0.000191f, -0.000132f, -0.000142f, -0.000254f, 0.000259f, 0.000619f, 0.000005f, -0.000422f, 0.000107f, -0.000338f, -0.000473f, 0.000159f, 0.000348f, -0.000276f, -0.000447f, 0.000203f, -0.000007f, -0.000111f, 0.000041f, -0.000116f, 0.000465f, 0.000333f, 0.000378f, -0.000291f, 0.000210f, 0.000299f, 0.000201f, 0.000117f, -0.000141f, 0.000103f, -0.000041f, -0.000100f, -0.000188f, -0.000033f, 0.001006f, 0.000497f, -0.000232f, 0.000013f, -0.000269f, 0.000412f, 0.000247f, 0.000300f, 0.000723f, -0.000339f, 0.000628f, -0.000049f, 0.000355f, 0.000323f, -0.000626f, 0.000629f, -0.000062f, 0.000131f, -0.000286f, -0.000021f, 0.000204f, -0.000010f, 0.000087f, 0.000132f, -0.000056f, 0.000051f, 0.000002f, 0.000337f, -0.000340f, -0.000219f, 0.000834f, 0.000616f, -0.000624f, 0.000147f, -0.000087f, 0.000300f, 0.000489f, -0.000011f, -0.000624f, -0.000099f, 0.000127f, -0.000044f, 0.000004f, 0.000291f, 0.000055f, 0.000298f, 0.000049f, 0.000135f, -0.000743f, 0.000150f, -0.000050f, -0.000347f, 0.000434f, 0.000666f, -0.000243f, -0.000579f, -0.000345f, -0.000251f, -0.000459f, 0.000107f, -0.000351f, -0.000176f, -0.000293f, 0.000179f, -0.000167f, 0.000417f, 0.000094f, 0.000281f, -0.000307f, -0.000182f, -0.000034f, 0.000035f, -0.000559f, 0.000026f, 0.000173f, 0.000035f, -0.000096f, -0.000400f, -0.000447f, 0.000124f, -0.000085f, -0.000223f, 0.000327f, 0.000146f, -0.000439f, -0.000157f, -0.000342f, -0.000507f, -0.000073f, -0.000019f, -0.000353f, 0.000015f, -0.000243f, -0.000108f, 0.000011f, 0.000263f, 0.000171f, -0.000389f, 0.000258f, 0.000044f, 0.000099f, 0.000466f, -0.000298f, 0.000586f, 0.000272f, 0.000220f, -0.000306f, -0.000378f, 0.000032f, -0.000057f, 0.000277f, 0.000017f, 0.000024f, 0.000277f, -0.000328f, -0.000219f, -0.000122f, -0.000093f, -0.000112f, 0.000240f, 0.000088f, 0.000433f, -0.000265f, 0.000249f, -0.000267f, 0.000089f, -0.000225f, 0.000113f, -0.000391f, 0.000567f, -0.000556f, 0.000414f, 0.000322f, -0.000480f, 0.000432f, 0.000052f, -0.000282f, 0.000466f, 0.000784f, 0.000430f, 0.000697f, -0.000298f, -0.000297f, -0.000359f, 0.000085f, 0.000017f, -0.000666f, 0.000063f, -0.000266f, 0.000365f, 0.000194f, -0.000591f, 0.000057f, -0.000072f, -0.000300f, 0.000184f, 0.000133f, 0.000097f, -0.000095f, -0.000255f, 0.000024f, 0.000048f, 0.000301f, -0.000282f, -0.000024f, -0.000111f, -0.000257f, -0.000080f, -0.000103f, 0.000468f, -0.000062f, -0.000203f, -0.000304f, 0.000148f, 0.000282f, -0.000252f, -0.000050f, 0.000046f, 0.000055f, -0.000014f, -0.000433f, 0.000076f, -0.000492f, 0.000062f, -0.000444f, 0.000088f, -0.000212f, -0.000328f, 0.000112f, 0.000091f, -0.000020f, -0.000709f, 0.000326f, 0.000039f, -0.000704f, -0.000386f, 0.000053f, -0.000366f, 0.000454f, -0.000107f, -0.000401f, -0.000409f, 0.000074f, -0.000084f, -0.000051f, 0.000475f, -0.000658f, 0.000020f, -0.000215f, 0.000300f, -0.000397f, -0.000478f, 0.000353f, -0.000102f, -0.000201f, -0.000323f, 0.000242f, 0.000097f, 0.000310f, -0.000233f, -0.000254f, 0.000105f, -0.000282f, 0.000445f, -0.000260f, -0.000370f, -0.000028f, 0.000266f, -0.000236f, 0.000052f, -0.000301f, 0.000082f, 0.000159f, 0.000116f, -0.000244f, -0.000322f, 0.000235f, -0.000240f, -0.000143f, 0.000437f, -0.000022f, -0.000150f, -0.000143f, 0.000073f, -0.000019f, 0.000109f, -0.000138f, 0.000112f, 0.000112f, -0.000093f, -0.000144f, -0.000545f, 0.000429f, -0.000090f, -0.000326f, -0.000147f, 0.000304f, 0.000179f, -0.000432f, -0.000050f, 0.000178f, -0.000044f, 0.000495f, 0.000005f, -0.000321f, 0.000054f, -0.000322f, -0.000076f, -0.000042f, 0.000064f, -0.000025f, 0.000065f, 0.000177f, -0.000111f, 0.000009f, 0.000498f, -0.000256f, 0.000059f, -0.000238f, -0.000084f, 0.000041f, -0.000501f, 0.000104f, 0.000357f, -0.000153f, -0.000277f, 0.000243f, -0.000243f, -0.000005f, 0.000269f, -0.000100f, 0.000073f, 0.000449f, 0.000078f, 0.000117f, 0.000601f, 0.000022f, -0.000207f, -0.000113f, 0.000378f, -0.000159f, 0.000079f, 0.000105f, -0.000233f, -0.000008f, -0.000131f, -0.000535f, 0.000050f, -0.000116f, -0.000032f, -0.000206f, 0.000489f, -0.000077f, -0.000264f, -0.000084f, -0.000102f, -0.000062f, -0.000351f, -0.000108f, -0.000222f, -0.000043f, -0.000169f, -0.000003f, 0.000262f, 0.000035f, 0.000103f, -0.000315f, -0.000024f, -0.000096f, -0.000111f, -0.000156f, -0.000112f, 0.000338f, -0.000045f, -0.000045f, -0.000008f, 0.000693f, -0.000087f, -0.000355f, -0.000383f, 0.000004f, -0.000117f, -0.000264f, -0.000377f, 0.000078f, 0.000091f, -0.000106f, 0.000238f, -0.000124f, 0.000302f, 0.000503f, -0.000388f, 0.000120f, 0.000358f, 0.000158f, 0.000012f, 0.000240f, 0.000137f, 0.000453f, 0.000507f, -0.000044f, -0.000016f, 0.000250f, -0.000122f, 0.000152f, 0.000278f, -0.000143f, 0.000234f, 0.000119f, 0.000412f, -0.000388f, 0.000258f, 0.000012f, 0.000122f, 0.000056f, 0.000430f, -0.000305f, 0.000120f, -0.000130f, 0.000161f, 0.000247f, 0.000033f, 0.000629f, 0.000269f, 0.000259f, -0.000127f, -0.000319f, 0.000169f, -0.000003f, 0.000527f, -0.000232f, -0.000006f, 0.000227f, 0.000215f, 0.000775f, -0.000059f, 0.000125f, 0.000228f, -0.000310f, -0.000155f, -0.000534f, 0.000270f, -0.000706f, -0.000101f, 0.000718f, 0.000385f, 0.000039f, 0.000225f, -0.000187f, 0.000119f, -0.000371f, -0.000114f, -0.000489f, -0.000069f, 0.000395f, -0.000208f, -0.000276f, 0.000358f, 0.000210f, -0.000380f, -0.000066f, 0.000470f, 0.000791f, -0.000476f, -0.000141f, -0.000319f, -0.000213f, 0.000108f, 0.000389f, -0.000446f, -0.000428f, 0.000257f, 0.000500f, -0.000667f, 0.000366f, -0.000409f, 0.000280f, -0.000381f, -0.000162f, 0.000025f, -0.000041f, 0.000429f, 0.000112f, -0.000490f, -0.000425f, -0.000438f, -0.000254f, -0.000173f, 0.000584f, -0.000015f, -0.000134f, -0.000179f, 0.000287f, -0.000079f, -0.000499f, 0.000173f, 0.000331f, -0.000458f, 0.000123f, -0.000065f, 0.000011f, 0.000024f, 0.000162f, 0.000175f, -0.000011f, 0.000261f, 0.000040f, -0.000193f, 0.000396f, 0.000465f, -0.000316f, -0.000014f, -0.000142f, -0.000434f, -0.000317f, -0.000250f, 0.000299f, 0.000675f, 0.000212f, 0.000039f, 0.000423f, -0.000351f, -0.000011f, -0.000094f, 0.000074f, -0.000501f, 0.000442f, 0.000331f, -0.000354f, 0.000263f, -0.000537f, -0.000067f, -0.000115f, -0.000116f, -0.000438f, 0.000223f, 0.000156f, -0.000158f, 0.000050f, -0.000298f, -0.000157f, -0.000045f, -0.000392f, 0.000692f, 0.000365f, -0.000028f, 0.000196f, -0.000232f, -0.000459f, -0.000229f, 0.000388f, 0.000349f, -0.000595f, 0.000146f, 0.000037f, 0.000036f, -0.000538f, 0.000329f, -0.000318f, 0.000469f, 0.000855f, 0.000200f, -0.000205f, 0.000175f, -0.000258f, 0.000078f, 0.000181f, -0.000303f, -0.000127f, 0.000126f, -0.000198f, -0.000296f, -0.000229f, 0.000561f, 0.000263f, -0.000109f, -0.000220f, -0.000076f, 0.000096f, -0.000202f, 0.000060f, -0.000371f, 0.000127f, 0.000568f, -0.000237f, -0.000336f, 0.000249f, -0.000038f, -0.000385f, -0.000071f, -0.000003f, 0.000357f, -0.000361f, -0.000273f, -0.000255f, -0.000313f, -0.000167f, -0.000037f, 0.000386f, 0.000224f, 0.000105f, 0.000232f, 0.000031f, -0.000005f, 0.000393f, -0.000322f, 0.000103f, 0.000101f, 0.000192f, 0.000185f, 0.000195f, -0.000032f, 0.000058f, 0.000559f, -0.000221f, 0.000079f, 0.000441f, -0.000250f, 0.000003f, 0.000208f, -0.000380f, -0.000009f, 0.000311f, -0.000382f, -0.000366f, 0.000475f, 0.000116f, -0.000264f, -0.000003f, -0.000498f, 0.000098f, 0.000391f, 0.000587f, -0.000090f, 0.000281f, -0.000323f, -0.000003f, 0.000092f, -0.000213f, 0.000004f, 0.000413f, -0.000021f, 0.000011f, 0.000196f, -0.000314f, 0.000347f, 0.000093f, -0.000169f, 0.000755f, -0.000098f, -0.000111f, -0.000346f, -0.000084f, -0.000169f, 0.000030f, 0.000014f, -0.000135f, -0.000466f, -0.000170f, 0.000251f, 0.000057f, -0.000243f, 0.000004f, -0.000339f, 0.000252f, 0.000566f, 0.000317f, -0.000266f, -0.000113f, 0.000137f, 0.000364f, -0.000112f, -0.000233f, 0.000016f, -0.000399f, 0.000133f, 0.000052f, -0.000064f, 0.000220f, 0.000285f, -0.000131f, -0.000023f, -0.000331f, 0.000513f, 0.000202f, -0.000128f, -0.000283f, 0.000303f, 0.000582f, -0.000192f, 0.000601f, -0.000252f, -0.000299f, -0.000211f, -0.000215f, 0.000228f, -0.000658f, 0.000034f, 0.000120f, 0.000408f, 0.000008f, -0.000281f, -0.000003f, 0.000080f, -0.000081f, -0.000189f, 0.000041f, -0.000523f, -0.000001f, 0.000088f, -0.000279f, -0.000111f, 0.000276f, 0.000169f, -0.000580f, -0.000258f, -0.000310f, -0.000311f, 0.000638f, -0.000074f, -0.000001f, -0.000549f, 0.000480f, 0.000014f, 0.000184f, 0.000802f, 0.000110f, 0.000365f, -0.000416f, -0.000070f, 0.000063f, -0.000186f, 0.000071f, -0.000297f, -0.000052f, 0.000294f, 0.000364f, -0.000054f, -0.000017f, -0.000176f, 0.000016f, 0.000084f, -0.000405f, -0.000180f, -0.000527f, -0.000143f, 0.000155f, -0.000340f, -0.000128f, -0.000214f, 0.000013f, 0.000053f, 0.000151f, -0.000027f, 0.000055f, 0.000596f, -0.000177f, -0.000131f, -0.000249f, 0.000594f, -0.000030f, -0.000253f, -0.000072f, 0.000058f, -0.000346f, -0.000226f, 0.000495f, -0.000012f, 0.000239f, 0.000001f, -0.000286f, -0.000437f, -0.000141f, -0.000162f, -0.000041f, 0.000016f, -0.000138f, 0.000107f, -0.000033f, -0.000017f, 0.000031f, 0.000028f, 0.000043f, 0.000068f, 0.000250f, -0.000355f, -0.000010f, -0.000054f, -0.000195f, -0.000087f, 0.000300f, 0.000064f, -0.000039f, 0.000098f, 0.000629f, -0.000179f, 0.000154f, -0.000145f, -0.000606f, 0.000160f, -0.000340f, 0.000172f, 0.000354f, -0.000410f, -0.000131f, 0.000244f, -0.000216f, -0.000024f, -0.000213f, -0.000129f, -0.000324f, -0.000325f, 0.000316f, 0.000071f, -0.000006f, -0.000042f, 0.000097f, -0.000007f, -0.000035f, -0.000178f, -0.000464f, 0.000343f, 0.000139f, 0.000101f, -0.000194f, -0.000239f, -0.000098f, 0.000236f, 0.000074f, 0.000599f, 0.000500f, -0.000193f, -0.000335f, 0.000383f, 0.000437f, -0.000054f, -0.000355f, -0.000462f, -0.000055f, 0.000220f, -0.000417f, 0.000134f, -0.000216f, 0.000078f, -0.000050f, -0.000443f, 0.000189f, -0.000085f, -0.000383f, 0.000307f, 0.000368f, 0.000243f, -0.000137f, -0.000081f, 0.000161f, 0.000238f, 0.000679f, 0.000082f, 0.000248f, 0.000090f, -0.000112f, 0.000566f, -0.000221f, -0.000282f, 0.000053f, 0.000038f, 0.000282f, 0.000404f, 0.000441f, 0.000283f, 0.000787f, -0.000019f, -0.000257f, 0.000114f, 0.000038f, -0.000447f, -0.000166f, 0.000219f, 0.000152f, 0.000173f, 0.000286f, -0.000041f, -0.000225f, 0.000107f, 0.000570f, -0.000163f, -0.000159f, 0.000065f, 0.000304f, 0.000063f, -0.000025f, 0.000003f, 0.000698f, -0.000305f, 0.000864f, -0.000169f, 0.000378f, -0.000144f, -0.000072f, -0.000259f, 0.000337f, 0.000209f, -0.000097f, -0.000241f, -0.000040f, -0.000209f, -0.000004f, 0.000136f, 0.000054f, -0.000229f, -0.000322f, -0.000528f, 0.000452f, -0.000179f, -0.000056f, -0.000129f, 0.000372f, 0.000102f, 0.000199f, 0.000054f, -0.000266f, -0.000022f, -0.000317f, -0.000306f, -0.000100f, -0.000266f, 0.000223f, 0.000274f, -0.000093f, -0.000191f, -0.000223f, 0.000103f, 0.000247f, -0.000114f, 0.000248f, 0.000249f, -0.000130f, 0.000009f, -0.000135f, -0.000301f, 0.000641f, 0.000161f, 0.000139f, -0.000158f, -0.000335f, 0.000091f, -0.000395f, 0.000131f, -0.000461f, 0.000099f, -0.000342f, -0.000094f, -0.000217f, -0.000106f, 0.000372f, -0.000083f, -0.000213f, -0.000151f, 0.000047f, -0.000301f, -0.000235f, -0.000153f, -0.000177f, -0.000187f, -0.000338f, 0.000042f, -0.000150f, 0.000028f, -0.000258f, -0.000163f, 0.000279f, -0.000081f, 0.000157f, -0.000009f, 0.000872f, -0.000325f, -0.000323f, -0.000034f, 0.000158f, 0.000040f, -0.000486f, 0.000428f, 0.000260f, -0.000574f, -0.000339f, -0.000136f, -0.000038f, 0.000406f, -0.000026f, -0.000331f, -0.000182f, 0.000312f, -0.000260f, 0.000240f, -0.000047f, 0.000207f, 0.000017f, 0.000095f, 0.000142f, -0.000221f, -0.000058f, 0.000537f, 0.000232f, -0.000249f, -0.000238f, -0.000299f, -0.000594f, -0.000059f, 0.000169f, 0.000368f, 0.000224f, 0.000218f, 0.000035f, -0.000441f, -0.000113f, -0.000074f, 0.000251f, -0.000045f, -0.000029f, -0.000327f, -0.000038f, 0.000493f, -0.000201f, -0.000133f, -0.000428f, -0.000232f, -0.000284f, -0.000339f, 0.000169f, 0.000276f, 0.000135f, -0.000115f, 0.000390f, -0.000052f, -0.000206f, 0.000354f, -0.000157f, -0.000071f, -0.000255f, 0.000305f, -0.000325f, -0.000326f, 0.000708f, 0.000887f, 0.000238f, -0.000635f, -0.000144f, -0.000443f, -0.000059f, 0.000166f, -0.000209f, -0.000098f, 0.000032f, 0.000018f, 0.000190f, -0.000085f, 0.000305f, 0.000295f, -0.000505f, 0.000087f, 0.000075f, 0.000184f, 0.000139f, -0.000243f, 0.000238f, -0.000210f, 0.000025f, -0.000017f, 0.000057f, -0.000578f, 0.000406f, -0.000070f, 0.000515f, -0.000103f, -0.000229f, 0.000651f, 0.000410f, 0.000135f, -0.000230f, 0.000044f, -0.000148f, 0.000291f, -0.000259f, -0.000273f, -0.000061f, 0.000114f, -0.000037f, 0.000163f, 0.000424f, 0.000239f, -0.000263f, -0.000065f, 0.000100f, -0.000235f, -0.000296f, 0.000275f, -0.000227f, -0.000148f, 0.000489f, 0.000402f, -0.000275f, -0.000133f, -0.000400f, 0.000004f, 0.000403f, 0.000239f, -0.000525f, -0.000020f, 0.000336f, 0.000163f, -0.000052f, -0.000253f, 0.000231f, 0.000681f, 0.000029f, -0.000084f, -0.000480f, -0.000116f, 0.000397f, -0.000089f, 0.000113f, 0.000126f, -0.000317f, -0.000146f, -0.000204f, -0.000170f, -0.000184f, 0.000019f, -0.000009f, 0.000573f, 0.000325f, -0.000171f, 0.000517f, -0.000091f, -0.000027f, 0.000156f, 0.000052f, 0.000031f, 0.000383f, -0.000012f, -0.000036f, 0.000141f, -0.000176f, -0.000132f, 0.000563f, 0.000513f, -0.000150f, -0.000104f, -0.000469f, 0.000031f, -0.000185f, 0.000122f, -0.000253f, 0.000295f, 0.000073f, 0.000066f, 0.000067f, -0.000181f, 0.000169f, -0.000180f, -0.000424f, -0.000251f, 0.000414f, -0.000154f, -0.000058f, -0.000319f, 0.000066f, -0.000000f, -0.000346f, 0.000150f, -0.000434f, -0.000021f, -0.000016f, 0.000391f, 0.000130f, -0.000168f, -0.000169f, -0.000293f, 0.000359f, -0.000342f, -0.000420f, -0.000265f, 0.000352f, 0.000290f, 0.000358f, -0.000451f, -0.000207f, -0.000108f, 0.000107f, 0.000001f, 0.000435f, -0.000047f, -0.000358f, 0.000079f, -0.000154f, -0.000098f, 0.000193f, 0.000362f, 0.000019f, -0.000124f, -0.000340f, 0.000038f, 0.000521f, -0.000243f, 0.000103f, -0.000076f, -0.000251f, 0.000033f, -0.000181f, -0.000235f, 0.000017f, -0.000225f, -0.000353f, 0.000334f, 0.000062f, -0.000110f, -0.000328f, 0.000032f, -0.000010f, -0.000144f, 0.000264f, 0.000088f, -0.000354f, 0.000304f, 0.000064f, -0.000199f, -0.000284f, -0.000186f, 0.000199f, -0.000138f, -0.000011f, -0.000565f, -0.000117f, -0.000133f, -0.000234f, 0.000320f, 0.000476f, -0.000329f, 0.000076f, -0.000023f, -0.000210f, -0.000068f, 0.000565f, -0.000110f, -0.000147f, 0.000733f, 0.000063f, -0.000005f, -0.000337f, -0.000243f, -0.000199f, 0.000116f, 0.000262f, 0.000628f, -0.000082f, 0.000316f, 0.000106f, -0.000343f, -0.000167f, 0.000342f, 0.000218f, -0.000082f, -0.000046f, -0.000126f, -0.000120f, -0.000165f, -0.000263f, -0.000064f, 0.000061f, 0.000092f, -0.000276f, -0.000048f, 0.000232f, 0.000253f, -0.000011f, 0.000249f, -0.000205f, -0.000142f, 0.000419f, 0.000194f, -0.000427f, -0.000220f, -0.000115f, 0.000011f, -0.000160f, -0.000513f, -0.000247f, 0.000179f, 0.000556f, -0.000248f, -0.000265f, -0.000107f, -0.000316f, -0.000295f, -0.000049f, 0.000282f, 0.000365f, 0.000537f, 0.000123f, -0.000267f, -0.000418f, -0.000303f, -0.000107f, -0.000048f, 0.000066f, 0.000464f, -0.000146f, 0.000255f, 0.000365f, -0.000325f, 0.000202f, 0.000004f, -0.000249f, 0.000390f, 0.000208f, -0.000284f, 0.000525f, -0.000151f, -0.000151f, 0.000856f, 0.000614f, -0.000310f, 0.000715f, 0.000681f, -0.000346f, -0.000177f, -0.000297f, -0.000225f, 0.000149f, 0.000182f, 0.000167f, 0.000036f, -0.000235f, 0.000201f, 0.000397f, 0.000333f, 0.000413f, -0.000456f, 0.000328f, 0.000083f, 0.000187f, 0.000145f, 0.000108f, -0.000234f, -0.000342f, 0.000158f, -0.000255f, 0.000238f, -0.000526f, 0.000194f, 0.000258f, 0.000237f, 0.000143f, 0.000197f, 0.000447f, -0.000054f, 0.000678f, -0.000017f, -0.000364f, 0.000032f, 0.000281f, -0.000002f, -0.000312f, 0.000092f, 0.000172f, -0.000252f, 0.000273f, 0.000407f, 0.000023f, -0.000070f, -0.000304f, 0.000249f, 0.000295f, 0.000161f, 0.000091f, 0.000002f, -0.000036f, 0.000385f, -0.000409f, 0.000052f, -0.000202f, -0.000052f, -0.000311f, 0.000087f, 0.000234f, 0.000204f, 0.000152f, -0.000192f, 0.000490f, -0.000196f, -0.000035f, 0.000403f, 0.000044f, 0.000191f, 0.000015f, -0.000682f, -0.000145f, 0.000273f, -0.000041f, 0.000253f, 0.000007f, -0.000448f, -0.000148f, -0.000244f, -0.000056f, -0.000237f, -0.000377f, 0.000078f, 0.000065f, -0.000254f, 0.000209f, 0.000150f, -0.000153f, -0.000025f, -0.000161f, -0.000088f, 0.000025f, 0.000149f, 0.000496f, 0.000700f, -0.000331f, -0.000548f, 0.000450f, -0.000470f, -0.000238f, -0.000055f, -0.000448f, -0.000148f, 0.000144f, -0.000148f, 0.000016f, -0.000375f, -0.000306f, 0.000397f, 0.000275f, -0.000248f, -0.000152f, 0.000328f, -0.000110f, -0.000055f, -0.000371f, -0.000188f, -0.000352f, -0.000093f, 0.000199f, -0.000044f, -0.000328f, -0.000372f, -0.000008f, 0.000100f, -0.000104f, -0.000139f, 0.000142f, -0.000413f, 0.000319f, 0.000105f, -0.000047f, -0.000051f, -0.000423f, 0.000146f, -0.000079f, 0.000060f, 0.000009f, 0.000078f, -0.000189f, -0.000233f, -0.000153f, -0.000208f, -0.000298f, -0.000186f, -0.000205f, 0.000525f, -0.000140f, 0.000127f, 0.000461f, 0.000021f, -0.000144f, 0.000308f, 0.000169f, -0.000211f, -0.000072f, 0.000113f, -0.000209f, -0.000211f, -0.000326f, -0.000387f, -0.000289f, -0.000131f, 0.000129f, 0.000001f, 0.000558f, -0.000378f, -0.000420f, 0.000517f, 0.000718f, 0.000110f, -0.000228f, 0.000012f, -0.000233f, 0.000066f, -0.000320f, 0.000210f, -0.000215f, 0.000157f, -0.000003f, 0.000266f, -0.000220f, -0.000118f, -0.000449f, -0.000236f, -0.000015f, -0.000029f, -0.000343f, -0.000069f, 0.000115f, -0.000146f, -0.000332f, 0.000159f, 0.000067f, -0.000258f, 0.000081f, -0.000139f, -0.000479f, -0.000055f, -0.000080f, -0.000150f, 0.000025f, 0.000203f, 0.000379f, 0.000001f, 0.000222f, 0.000603f, -0.000010f, -0.000157f, -0.000123f, -0.000120f, -0.000239f, -0.000128f, 0.000083f, 0.000343f, 0.000652f, -0.000233f, 0.000264f, -0.000004f, 0.000131f, 0.000173f, -0.000404f, -0.000197f, 0.000166f, -0.000158f, -0.000153f, -0.000029f, -0.000129f, -0.000328f, -0.000183f, 0.000243f, 0.000240f, -0.000326f, -0.000010f, -0.000095f, 0.000070f, 0.000195f, 0.000022f, -0.000282f, 0.000069f, -0.000180f, 0.000075f, 0.000104f, 0.000143f, -0.000293f, 0.000117f, -0.000042f, 0.000305f, 0.000537f, -0.000420f, 0.000184f, 0.000118f, -0.000240f, -0.000321f, 0.000276f, -0.000329f, 0.000179f, 0.000084f, -0.000267f, 0.000153f, -0.000043f, 0.000076f, -0.000331f, -0.000285f, 0.000317f, -0.000034f, 0.000558f, 0.000302f, 0.000085f, -0.000141f, -0.000134f, -0.000178f, -0.000030f, -0.000112f, -0.000037f, 0.000202f, -0.000024f, -0.000004f, -0.000024f, -0.000131f, -0.000214f, -0.000307f, -0.000130f, -0.000030f, -0.000451f, -0.000134f, -0.000034f, 0.000243f, 0.000092f, 0.000220f, -0.000328f, -0.000436f, -0.000007f, -0.000185f, -0.000295f, -0.000168f, -0.000081f, -0.000282f, -0.000267f, -0.000114f, 0.000072f, 0.000008f, -0.000425f, 0.000031f, 0.000080f, -0.000024f, 0.000050f, -0.000183f, -0.000176f, -0.000100f, 0.000022f, -0.000129f, -0.000012f, -0.000289f, 0.000087f, -0.000021f, 0.000353f, -0.000092f, -0.000046f, 0.000095f, -0.000215f, 0.000248f, 0.000171f, 0.000134f, -0.000232f, 0.000464f, -0.000124f, 0.000601f, 0.000348f, 0.000068f, 0.000376f, -0.000143f, 0.000111f, 0.000225f, 0.000286f, 0.000234f, -0.000022f, 0.000050f, 0.000423f, 0.000051f, -0.000106f, -0.000155f, 0.000033f, 0.000304f, 0.000460f, 0.000534f, -0.000204f, -0.000197f, 0.000380f, 0.000489f, 0.000023f, -0.000193f, -0.000060f, 0.000173f, 0.000020f, -0.000300f, 0.000349f, -0.000209f, 0.000262f, 0.000072f, -0.000180f, -0.000165f, -0.000262f, -0.000054f, 0.000257f, 0.000064f, -0.000280f, 0.000279f, -0.000251f, -0.000232f, -0.000101f, 0.000070f, -0.000145f, -0.000020f, -0.000355f, 0.000371f, 0.000039f, 0.000015f, -0.000031f, 0.000222f, -0.000064f, -0.000085f, -0.000060f, 0.000001f, 0.000209f, 0.000360f, 0.000466f, -0.000369f, -0.000083f, -0.000079f, 0.000120f, -0.000283f, -0.000372f, -0.000384f, -0.000106f, -0.000060f, -0.000261f, -0.000015f, 0.000584f, -0.000069f, 0.000416f, -0.000254f, 0.000187f, -0.000181f, -0.000055f, -0.000169f, -0.000239f, 0.000093f, 0.000068f, 0.000213f, 0.000199f, 0.000305f, -0.000241f, 0.000009f, 0.000202f, 0.000048f, 0.000196f, -0.000297f, -0.000137f, -0.000117f, -0.000104f, -0.000078f, 0.000192f, 0.000336f, 0.000849f, 0.000266f, 0.000197f, 0.000346f, 0.000196f, -0.000438f, -0.000056f, -0.000149f, 0.000122f, -0.000043f, -0.000097f, 0.000030f, 0.000230f, -0.000154f, 0.000459f, 0.000165f, -0.000360f, 0.000374f, -0.000035f, -0.000334f, 0.000030f, -0.000061f, 0.000480f, 0.000035f, 0.000022f, -0.000010f, 0.000129f, 0.000220f, 0.000240f, -0.000063f, 0.000100f, -0.000028f, 0.000064f, 0.000196f, 0.000303f, -0.000327f, 0.000193f, 0.000068f, 0.000212f, -0.000323f, -0.000002f, 0.000524f, 0.000100f, 0.000060f, -0.000288f, -0.000180f, -0.000256f, 0.000013f, -0.000152f, -0.000381f, -0.000248f, 0.000089f, -0.000279f, -0.000218f, 0.000015f, -0.000225f, 0.000495f, -0.000171f, -0.000049f, 0.000019f, 0.000026f, -0.000115f, 0.000148f, -0.000144f, -0.000300f, 0.000069f, 0.000493f, -0.000345f, 0.000409f, 0.000023f, -0.000177f, 0.000189f, -0.000493f, 0.000156f, }
    ;

    struct inline8 : as_engine<
          8, /* W             */
          9, /* H             */
          2, /* THROTTLE_Y    */
          4, /* PISTON_Y      */
          5, /* AUDIO_Y       */
        200, /* PIPE_CELLS    */
         10, /* PIPE_SUBSTEPS */
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
            this->clutch.damping_coefficient_n_m_s = 10.0;
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
            this->throttle.mapping.table = {
                0.001,
                0.100,
                0.300,
                1.000,
            };
            // 1-8-4-3-6-5-7-2
            pipes[0].piston_connect_m = { 0.00, 0.34, 0.08, 0.23 };
            pipes[1].piston_connect_m = { 0.34, 0.23, 0.00, 0.08 };
            for(auto& pipe : this->pipes)
            {
                pipe.mic_position_m = pipe.length_m = 1.0;
            }
            this->dc.set_cutoff_frequency(5.0);
            this->gain.ratio = 0.0001;
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
