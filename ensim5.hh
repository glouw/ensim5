#pragma once

#include <cmath>

#define fn //__attribute__((used))

namespace ensim5
{
    using std::sin, std::cos, std::sqrt, std::fmax, std::fmin;

    template<typename T> requires(std::same_as<T, float> || std::same_as<T, double>)
    class dbl_c
    {
        T value;
    public:
        constexpr dbl_c() = default;
        constexpr dbl_c(const double value): value(static_cast<T>(value)) {}
        constexpr dbl_c& operator+=(const dbl_c other) { value += other.value; return *this; }
        constexpr operator T() const { return value; }
    };

    using val_t = dbl_c<float>;

    constexpr size_t g_sample_rate_hz = 48000;
    constexpr val_t g_dt_s = val_t(1.0) / g_sample_rate_hz;
    constexpr val_t g_sink_m3 = 1e6;
    constexpr val_t g_ambient_static_temperature_k = 293.0;
    constexpr val_t g_ambient_static_pressure_pa = 101325.0;

    template<size_t H> requires(H % 2 == 1)
    struct flow_s
    {
        /* Assume chamber static pressure == chamber total pressure. */
        /* Assume chamber static temperature == chamber total temperature. */
        /* Assume chamber static density == chamber total density. */
        /* Assume nozzle flow fixes gamma at 1.333 to generate SIMD friendly instructions. */
        /* Assume nozzle static density is calculated from nozzle mach, nozzle gamma, and upstream static chamber conditions. */
        /* Assume chambers do not fix gamma and track gas mixture at constant volume (cv) to mix static temperature. */
        /* Assume chamber static pressure mechanics are derived from PV=nRT. */
        /* Assume chamber constant volume (cv) constants are averaged across Otto Cycle static temperatures and their specific cycle lifetimes,
         * and align roughly with the nozzle gamma value of 1.333. */

        static constexpr val_t nozzle_gamma = 1.333;
        static constexpr val_t chamber_cv_n2_j_per_mol_k = 23.0;
        static constexpr val_t chamber_cv_o2_j_per_mol_k = 22.5;
        static constexpr val_t chamber_cv_ar_j_per_mol_k = 12.5;
        static constexpr val_t chamber_cv_c8h18_j_per_mol_k = 250.0;
        static constexpr val_t chamber_cv_co2_j_per_mol_k = 33.0;
        static constexpr val_t chamber_cv_h2o_j_per_mol_k = 30.0;
        static constexpr val_t chamber_molar_mass_kg_per_mol_c8h18 = 0.1142285200;
        static constexpr val_t chamber_molar_mass_kg_per_mol_o2 = 0.0319988000;
        static constexpr val_t chamber_molar_mass_kg_per_mol_n2 = 0.0280134000;
        static constexpr val_t chamber_molar_mass_kg_per_mol_ar = 0.0399480000;
        static constexpr val_t chamber_molar_mass_kg_per_mol_co2 = 0.0440095000;
        static constexpr val_t chamber_molar_mass_kg_per_mol_h2o = 0.0180152800;
        static constexpr val_t chamber_universal_gas_constant_j_per_mol_k = 8.3144598;

        using lane_t = val_t[H];
        using mask_t = bool[H];

        lane_t chamber_volume_m3 = {};
        lane_t chamber_nozzle_flow_area_m2 = {};
        lane_t chamber_molar_mass_kg_per_mol = {};
        lane_t chamber_cv_j_per_mol_k = {};
        lane_t chamber_specific_gas_constant_j_per_kg_k = {};
        lane_t chamber_static_pressure_pa = {};
        lane_t chamber_moles = {};
        lane_t chamber_total_cv_j_per_k = {};
        lane_t chamber_mol_ratio_n2 = {};
        lane_t chamber_mol_ratio_o2 = {};
        lane_t chamber_mol_ratio_ar = {};
        lane_t chamber_mol_ratio_c8h18 = {};
        lane_t chamber_mol_ratio_co2 = {};
        lane_t chamber_mol_ratio_h2o = {};
        lane_t chamber_static_temperature_k = {};
        lane_t chamber_mass_kg = {};

        lane_t nozzle_mach = {};
        lane_t nozzle_velocity_m_per_s = {};
        lane_t nozzle_static_density_kg_per_m3 = {};
        lane_t nozzle_mass_flow_rate_kg_per_s = {};

        /* Parcels are tiny chamber masses determined by nozzle mass flow rate. */
        lane_t parcel_molar_mass_kg_per_mol = {};
        lane_t parcel_cv_j_per_mol_k = {};
        lane_t parcel_moles = {};
        lane_t parcel_total_cv_j_per_k = {};
        lane_t parcel_mol_ratio_n2 = {};
        lane_t parcel_mol_ratio_o2 = {};
        lane_t parcel_mol_ratio_ar = {};
        lane_t parcel_mol_ratio_c8h18 = {};
        lane_t parcel_mol_ratio_co2 = {};
        lane_t parcel_mol_ratio_h2o = {};
        lane_t parcel_static_temperature_k = {};
        lane_t parcel_mass_kg = {};
        mask_t parcel_flow_right = {};

        fn val_t clamp(const val_t val, const val_t lower, const val_t upper)
        {
            return fmax(fmin(val, val_t(upper)), val_t(lower));
        }

        fn val_t calc_mix(const val_t v1, const val_t w1, const val_t v2, const val_t w2)
        {
            return (v1 * w1 + v2 * w2) / (w1 + w2);
        }

        /*
         *     Ps * V
         * m = -------
         *     Rs * Ts
         */

        fn void calc_chamber_ambients()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                chamber_static_temperature_k[i] = g_ambient_static_temperature_k;
                chamber_static_pressure_pa[i] = g_ambient_static_pressure_pa;
            }
            for(size_t i = 0; i < H - 1; i++)
            {
                const val_t Ps = chamber_static_pressure_pa[i];
                const val_t V = chamber_volume_m3[i];
                const val_t Rs = chamber_specific_gas_constant_j_per_kg_k[i];
                const val_t Ts = chamber_static_temperature_k[i];
                chamber_mass_kg[i] = Ps * V / (Rs * Ts);
            }
        }

        /*
         *      ____
         *      \
         * cv = /___ mr[i] * M
         *         i
         */

        fn void calc_chamber_molar_masses()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                chamber_molar_mass_kg_per_mol [i]
                    = chamber_mol_ratio_n2    [i] * chamber_molar_mass_kg_per_mol_n2
                    + chamber_mol_ratio_o2    [i] * chamber_molar_mass_kg_per_mol_o2
                    + chamber_mol_ratio_ar    [i] * chamber_molar_mass_kg_per_mol_ar
                    + chamber_mol_ratio_c8h18 [i] * chamber_molar_mass_kg_per_mol_c8h18
                    + chamber_mol_ratio_co2   [i] * chamber_molar_mass_kg_per_mol_co2
                    + chamber_mol_ratio_h2o   [i] * chamber_molar_mass_kg_per_mol_h2o;
            }
        }

        /*
         *      ____
         *      \
         * cv = /___ mr[i] * cv
         *         i
         */

        fn void calc_chamber_molar_cvs()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                chamber_cv_j_per_mol_k[i]
                    = chamber_mol_ratio_n2    [i] * chamber_cv_n2_j_per_mol_k
                    + chamber_mol_ratio_o2    [i] * chamber_cv_o2_j_per_mol_k
                    + chamber_mol_ratio_ar    [i] * chamber_cv_ar_j_per_mol_k
                    + chamber_mol_ratio_c8h18 [i] * chamber_cv_c8h18_j_per_mol_k
                    + chamber_mol_ratio_co2   [i] * chamber_cv_co2_j_per_mol_k
                    + chamber_mol_ratio_h2o   [i] * chamber_cv_h2o_j_per_mol_k;
            }
        }

        /*
         *       R
         * Rs = ---
         *       M
         */

        fn void calc_chamber_specific_gas_constants()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                const val_t R = chamber_universal_gas_constant_j_per_mol_k;
                const val_t M = chamber_molar_mass_kg_per_mol[i];
                chamber_specific_gas_constant_j_per_kg_k[i] = R / M;
            }
        }

        /*
         *      m * Rs * Ts
         * Ps = -----------
         *           V
         */

        fn void calc_chamber_static_pressures()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                const val_t m = chamber_mass_kg[i];
                const val_t Rs = chamber_specific_gas_constant_j_per_kg_k[i];
                const val_t Ts = chamber_static_temperature_k[i];
                const val_t V = chamber_volume_m3[i];
                chamber_static_pressure_pa[i] = m * Rs * Ts / V;
            }
        }

        /*
         *      m
         * n = ---
         *      M
         */

        fn void calc_chamber_moles()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                const val_t m = chamber_mass_kg[i];
                const val_t M = chamber_molar_mass_kg_per_mol[i];
                chamber_moles[i] = m / M;
            }
        }

        /*
         * cv(total) = n * cv
         */

        fn void calc_chamber_total_cvs()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                chamber_total_cv_j_per_k[i] = chamber_moles[i] * chamber_cv_j_per_mol_k[i];
            }
        }

        /*
         *              ____________________
         *             /
         *            /            y - 1
         *           /             -----
         *          /                y              ____________
         *         /   2         Pt                /
         * M = _  /  ----- * [ (----) - 1 ]  = _  / A * (B - 1)
         *      \/   y - 1       Ps             \/
         */

        fn void calc_nozzle_machs()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                const size_t j = i + 1;
                const val_t X = (val_t(2.0) / (nozzle_gamma - val_t(1.0)));
                const val_t Px = chamber_static_pressure_pa[i];
                const val_t Py = chamber_static_pressure_pa[j];
                const val_t Pt = fmax(Px, Py);
                const val_t Ps = fmin(Px, Py);
                const val_t direction = Px > Py ? val_t(1.0) : val_t(-1.0);
                parcel_flow_right[i] = direction > val_t(0.0);

                /*
                 *      y - 1    1.333 - 1                        1     1          __________
                 *      ----- = ---------- = 0.25                --- * ---        /    ______
                 *        y       1.333                0.25       2     2        /    /
                 *  term                         = term    = term          = _  / _  /  term
                 *                                                            \/   \/
                 */

                static_assert(nozzle_gamma == val_t(1.333), "nozzle_gamma (y) must be 1.333 to use manual pow optimization");
                const val_t Y = sqrt(sqrt(Pt / Ps));
                const val_t M = direction * sqrt(X * (Y - val_t(1.0)));
                nozzle_mach[i] = clamp(M, val_t(-1.0), val_t(1.0));
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
            for(size_t i = 0; i < H - 1; i++)
            {
                const val_t Rs = chamber_specific_gas_constant_j_per_kg_k[i];
                const val_t Tt = chamber_static_temperature_k[i];
                const val_t M = nozzle_mach[i];
                const val_t X = nozzle_gamma * Rs * Tt;
                const val_t Y = val_t(0.5) * (nozzle_gamma - val_t(1.0)) * M * M;
                const val_t u = M * sqrt(X / (val_t(1.0) + Y));
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
            for(size_t i = 0; i < H - 1; i++)
            {
                const val_t Pt = chamber_static_pressure_pa[i];
                const val_t Rs = chamber_specific_gas_constant_j_per_kg_k[i];
                const val_t Tt = chamber_static_temperature_k[i];
                const val_t M = nozzle_mach[i];
                const val_t X = Pt / (Rs * Tt);
                const val_t C = val_t(1.0) + val_t(0.5) * (nozzle_gamma - val_t(1.0)) * M * M;

                /*
                 *         1         1
                 *       ----- = --------- = 3.00
                 *       y - 1   1.333 - 1
                 *  term                         = term * term * term
                 */

                static_assert(nozzle_gamma == val_t(1.333), "nozzle_gamma (y) must be 1.333 to use manual pow optimization");
                nozzle_static_density_kg_per_m3[i] = X / (C * C * C);
            }
        }

        fn void calc_parcel_property(lane_t& __restrict nozzle, const lane_t& __restrict chamber)
        {
            for(size_t i = 0; i < H - 1; i++) nozzle[i] = val_t(0.0);
            for(size_t i = 0; i < H - 1; i++) nozzle[i] += parcel_flow_right[i] ? chamber[i + 0] : val_t(0.0);
            for(size_t i = 0; i < H - 1; i++) nozzle[i] += parcel_flow_right[i] ? val_t(0.0) : chamber[i + 1];
        }

        fn void calc_parcel_properties()
        {
            calc_parcel_property(parcel_molar_mass_kg_per_mol, chamber_molar_mass_kg_per_mol);
            calc_parcel_property(parcel_cv_j_per_mol_k,        chamber_cv_j_per_mol_k);
            calc_parcel_property(parcel_mol_ratio_n2,          chamber_mol_ratio_n2);
            calc_parcel_property(parcel_mol_ratio_o2,          chamber_mol_ratio_o2);
            calc_parcel_property(parcel_mol_ratio_ar,          chamber_mol_ratio_ar);
            calc_parcel_property(parcel_mol_ratio_c8h18,       chamber_mol_ratio_c8h18);
            calc_parcel_property(parcel_mol_ratio_co2,         chamber_mol_ratio_co2);
            calc_parcel_property(parcel_mol_ratio_h2o,         chamber_mol_ratio_h2o);
            calc_parcel_property(parcel_static_temperature_k,  chamber_static_temperature_k);
        }

        /* .
         * m = ps A u
         *
         */

        fn void calc_nozzle_mass_flow_rates()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                const val_t ps = nozzle_static_density_kg_per_m3[i];
                const val_t A = chamber_nozzle_flow_area_m2[i];
                const val_t u = nozzle_velocity_m_per_s[i];
                const val_t mdot = ps * A * u;
                nozzle_mass_flow_rate_kg_per_s[i] = mdot;
                const val_t dm = mdot * g_dt_s;
                parcel_mass_kg[i] = dm;
                parcel_moles[i] = dm / parcel_molar_mass_kg_per_mol[i];
            }
        }

        /*
         *        nr1 n1 + nr2 n2
         * nr1 = --------------------
         *            n1 + n2
         */

        fn void calc_chamber_mol_ratios()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                chamber_mol_ratio_n2   [i] = calc_mix(chamber_mol_ratio_n2   [i], chamber_moles[i], parcel_mol_ratio_n2   [i], parcel_moles[i]);
                chamber_mol_ratio_o2   [i] = calc_mix(chamber_mol_ratio_o2   [i], chamber_moles[i], parcel_mol_ratio_o2   [i], parcel_moles[i]);
                chamber_mol_ratio_ar   [i] = calc_mix(chamber_mol_ratio_ar   [i], chamber_moles[i], parcel_mol_ratio_ar   [i], parcel_moles[i]);
                chamber_mol_ratio_c8h18[i] = calc_mix(chamber_mol_ratio_c8h18[i], chamber_moles[i], parcel_mol_ratio_c8h18[i], parcel_moles[i]);
                chamber_mol_ratio_co2  [i] = calc_mix(chamber_mol_ratio_co2  [i], chamber_moles[i], parcel_mol_ratio_co2  [i], parcel_moles[i]);
                chamber_mol_ratio_h2o  [i] = calc_mix(chamber_mol_ratio_h2o  [i], chamber_moles[i], parcel_mol_ratio_h2o  [i], parcel_moles[i]);
            }
        }

        /*
         * cv(total) = n * cv
         */

        fn void calc_parcel_total_cvs()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                parcel_total_cv_j_per_k[i] = parcel_moles[i] * parcel_cv_j_per_mol_k[i];
            }
        }

        /*
         *        n1 cv1 Ts1 + n2 cv2 Ts2
         * Ts1 = -------------------------
         *            n1 cv1 + n2 cv2
         */

        fn void calc_chamber_static_temperatures()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                chamber_static_temperature_k[i] = calc_mix(
                    chamber_static_temperature_k[i],
                    chamber_total_cv_j_per_k[i],
                    parcel_static_temperature_k[i],
                    parcel_total_cv_j_per_k[i]
                );
            }
        }

        fn void calc_chamber_masses()
        {
            for(size_t i = 0; i < H - 1; i++)
            {
                chamber_mass_kg[i] += parcel_mass_kg[i];
            }
        }

        fn void calc_chamber_setup()
        {
            calc_chamber_molar_masses();
            calc_chamber_molar_cvs();
            calc_chamber_specific_gas_constants();
        }

        fn void calc_chamber_ambient_reset()
        {
            calc_chamber_setup();
            calc_chamber_ambients();
        }

        fn void calc_flow()
        {
            calc_chamber_setup();
            calc_chamber_static_pressures();
            calc_chamber_moles();
            calc_chamber_total_cvs();
            calc_nozzle_machs();
            calc_nozzle_velocities();
            calc_nozzle_static_densities();
            calc_parcel_properties();
            calc_nozzle_mass_flow_rates();
            calc_chamber_mol_ratios();
            calc_parcel_total_cvs();
            calc_chamber_static_temperatures();
            calc_chamber_masses();
        }
    };

    struct crankshaft_s
    {
        val_t mass_kg = {};
        val_t radius_m = {};
        val_t theta_r = {};
        val_t angular_velocity_r_per_s = {};
        val_t moment_of_inertia_kg_m2 = {};

        /*
         * dw = dw/dt * dt
         *
         */

        fn void accelerate(const val_t angular_acceleration_r_per_s2)
        {
            angular_velocity_r_per_s += angular_acceleration_r_per_s2 * g_dt_s;
        }

        /*
         * dth = dth/dt * dt
         *
         */

        fn void turn()
        {
            theta_r += angular_velocity_r_per_s * g_dt_s;
        }

        /*
         *      1         2
         * I = --- * m * r
         *      2
         *
         */

        fn void calc_moment_of_inertia()
        {
            moment_of_inertia_kg_m2 = val_t(0.5) * mass_kg * radius_m * radius_m;
        }

        fn void calc_crankshaft()
        {
            calc_moment_of_inertia();
        }
    };

    template<size_t W, size_t PY>
    struct pistons_s
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

        using lane_t = val_t[W];

        lane_t diameter_m = {};
        lane_t crank_throw_length_m = {};
        lane_t connecting_rod_length_m = {};
        lane_t connecting_rod_mass_kg = {};
        lane_t head_mass_density_kg_per_m3 = {};
        lane_t head_compression_height_m = {};
        lane_t head_clearance_height_m = {};
        lane_t theta0_r = {};
        lane_t theta_r = {};
        lane_t sint = {};
        lane_t cost = {};
        lane_t pin_x_m = {};
        lane_t pin_y_m = {};
        lane_t bearing_x_m = {};
        lane_t bearing_y_m = {};
        lane_t volumes_m3 = {};
        lane_t head_mass_kg = {};
        lane_t moment_of_inertia_kg_m2 = {};
        lane_t gas_torque_n_m = {};
        lane_t inertia_torque_n_m = {};

        /*
         * t = t0 + t1
         */

        fn void calc_thetas(const crankshaft_s& __restrict crankshaft)
        {
            for(size_t i = 0; i < W; i++)
            {
                theta_r[i] = theta0_r[i] + crankshaft.theta_r;
            }
        }

        fn void calc_sin_cos()
        {
            for(size_t i = 0; i < W; i++)
            {
                sint[i] = sin(theta_r[i]);
                cost[i] = cos(theta_r[i]);
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
                const val_t r = crank_throw_length_m[i];
                const val_t l = connecting_rod_length_m[i];
                const val_t t = theta_r[i];
                const val_t x = r * sint[i];
                const val_t y = r * cost[i];
                bearing_x_m[i] = x;
                bearing_y_m[i] = y;
                pin_x_m[i] = val_t(0.0);
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
                const val_t r = crank_throw_length_m[i];
                const val_t l = connecting_rod_length_m[i];
                const val_t cm = head_compression_height_m[i];
                const val_t cl = head_clearance_height_m[i];
                const val_t block_deck_surface_m = r + l + cm + cl;
                const val_t y = pin_y_m[i] + cm;
                const val_t radius = diameter_m[i] / val_t(2.0);
                const val_t h = block_deck_surface_m - y;
                volumes_m3[i] = val_t(M_PI) * radius * radius * h;
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
                const val_t r = val_t(0.5) * diameter_m[i];
                const val_t h = val_t(2.0) * head_compression_height_m[i];
                const val_t p = head_mass_density_kg_per_m3[i];
                head_mass_kg[i] = val_t(M_PI) * r * r * h * p;
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
                const val_t r = crank_throw_length_m[i];
                const val_t mp = head_mass_kg[i];
                const val_t mr = connecting_rod_mass_kg[i];
                moment_of_inertia_kg_m2[i] = (mp + (val_t(1.0 / 3.0) * mr)) * r * r;
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

        fn void calc_gas_torques(const auto& __restrict flow)
        {
            for(size_t i = 0; i < W; i++)
            {
                const val_t Pg = flow.chamber_static_pressure_pa[PY];
                const val_t A = val_t(M_PI) * diameter_m[i] * diameter_m[i];
                const val_t r = crank_throw_length_m[i];
                const val_t l = connecting_rod_length_m[i];
                const val_t t = theta_r[i];
                const val_t X = Pg * A * r * sint[i];
                const val_t Y = val_t(1.0) + (r / l) * cost[i];
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

        fn void calc_inertia_torques(const crankshaft_s& __restrict crankshaft)
        {
            for(size_t i = 0; i < W; i++)
            {
                const val_t r = crank_throw_length_m[i];
                const val_t l = connecting_rod_length_m[i];
                const val_t I = moment_of_inertia_kg_m2[i];
                const val_t w = crankshaft.angular_velocity_r_per_s;
                const val_t t = theta_r[i];
                const val_t rl = r / l;
                const val_t s = sint[i];
                const val_t c = cost[i];
                const val_t X = val_t(0.25) * rl * s;
                const val_t Y = s * c;
                const val_t Z = val_t(0.75) * rl * (val_t(3.0) * s - val_t(4.0) * s * s * s);
                inertia_torque_n_m[i] = I * w * w * (X - Y - Z);
            }
        }

        fn void calc_volumetrics(const crankshaft_s& __restrict crankshaft)
        {
            calc_thetas(crankshaft);
            calc_sin_cos();
            calc_positions();
            calc_volumes();
        }

        fn void calc_pistons(const auto& __restrict flow, const crankshaft_s& __restrict crankshaft)
        {
            calc_volumetrics(crankshaft);
            calc_masses();
            calc_moments_of_inertia();
            calc_gas_torques(flow);
            calc_inertia_torques(crankshaft);
        }
    };

    template <size_t W, size_t H, size_t PY>
    struct engine_s
    {
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

        crankshaft_s crankshaft;
        pistons_s<W, PY> pistons;
        flow_s<H> flow[W];

        fn void rig_pistons()
        {
            pistons.calc_volumetrics(crankshaft);
        }

        fn void rig_flow_chambers()
        {
            for(size_t x = 0; x < W; x++)
            {
                flow[x].chamber_volume_m3[PY] = pistons.volumes_m3[x];
            }
        }

        fn void reset_flow_chambers()
        {
            for(size_t x = 0; x < W; x++)
            {
                flow[x].calc_chamber_ambient_reset();
            }
        }

        fn void setup_engine()
        {
            rig_pistons();
            rig_flow_chambers();
            reset_flow_chambers();
        }

        fn void step_crankshaft()
        {
            crankshaft.calc_crankshaft();
        }

        fn void step_pistons()
        {
            for(size_t i = 0; i < W; i++)
            {
                pistons.calc_pistons(flow[i], crankshaft);
            }
        }

        fn void step_flow()
        {
            for(size_t i = 0; i < W; i++)
            {
                flow[i].calc_flow();
            }
        }

        fn void step_engine()
        {
            step_crankshaft();
            step_pistons();
            step_flow();
        }
    };

    auto init_engine_8_9_5()
    {
        engine_s<8, 9, 5> engine = {
            .crankshaft = {
                .mass_kg = 0.1,
                .radius_m = 0.25,
            },
            .pistons = {
                .diameter_m                  = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 },
                .crank_throw_length_m        = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 },
                .connecting_rod_length_m     = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 },
                .connecting_rod_mass_kg      = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 },
                .head_mass_density_kg_per_m3 = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 },
                .head_compression_height_m   = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 },
                .head_clearance_height_m     = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 },
                .theta0_r                    = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 },
            },
            .flow = {{
                .chamber_volume_m3           = { g_sink_m3,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1, g_sink_m3 },
                .chamber_nozzle_flow_area_m2 = {      0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01      },
            },{
                .chamber_volume_m3           = { g_sink_m3,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1, g_sink_m3 },
                .chamber_nozzle_flow_area_m2 = {      0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01      },
            },{
                .chamber_volume_m3           = { g_sink_m3,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1, g_sink_m3 },
                .chamber_nozzle_flow_area_m2 = {      0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01      },
            },{
                .chamber_volume_m3           = { g_sink_m3,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1, g_sink_m3 },
                .chamber_nozzle_flow_area_m2 = {      0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01      },
            },{
                .chamber_volume_m3           = { g_sink_m3,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1, g_sink_m3 },
                .chamber_nozzle_flow_area_m2 = {      0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01      },
            },{
                .chamber_volume_m3           = { g_sink_m3,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1, g_sink_m3 },
                .chamber_nozzle_flow_area_m2 = {      0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01      },
            },{
                .chamber_volume_m3           = { g_sink_m3,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1, g_sink_m3 },
                .chamber_nozzle_flow_area_m2 = {      0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01      },
            },{
                .chamber_volume_m3           = { g_sink_m3,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1, g_sink_m3 },
                .chamber_nozzle_flow_area_m2 = {      0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01      },
            }}
        };
        engine.setup_engine();
        return engine;
    }
}
