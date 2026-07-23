#include "ensim.hh"

#include <cmath>
#include <array>
#include <numbers>

#define fn __attribute__((used))

namespace ensim
{
    using std::sin, std::cos, std::sqrt, std::fmax, std::fmin, std::cbrt;

    static constexpr double pi = std::numbers::pi_v<double>;
    static constexpr int sample_rate_hz = 48000;
    static constexpr double dt_s = 1.0 / sample_rate_hz;
    static constexpr double otto_cycle_r = 4.0 * pi;
    static constexpr double otto_intake_cycle_r = 0.0 * pi;
    static constexpr double otto_compression_cycle_r = 1.0 * pi;
    static constexpr double otto_combustion_cycle_r = 2.0 * pi;
    static constexpr double otto_exhaust_cycle_r = 3.0 * pi;
    static constexpr double resevoir_volume_m3 = 1e9;

    consteval bool is_odd(const int value) { return value % 2 == 1; }
    consteval bool is_evn(const int value) { return value % 2 == 0; }

    constexpr double clamp(const double value, const double lower, const double upper)
    {
        return fmax(fmin(value, upper), lower);
    }

    constexpr double modulos(const double value, const double by)
    {
        return value - trunc(value / by) * by;
    }

    template<int H>
    requires(is_odd(H))
    struct flow
    {
        static constexpr int N = H - 1;
        static_assert(is_evn(N));

        static constexpr double nozzle_gamma = 1.33333;
        static constexpr double universal_gas_constant_j_per_mol_k = 8.3144598;
        static constexpr double cv_j_per_mol_k = universal_gas_constant_j_per_mol_k / (nozzle_gamma - 1.0);
        static constexpr double molar_mass_kg_per_mol = 0.023;
        static constexpr double specific_gas_constant_j_per_kg_k = universal_gas_constant_j_per_mol_k / molar_mass_kg_per_mol;
        static constexpr double ambient_static_temperature_k = 300.0;
        static constexpr double ambient_static_pressure_pa = 101325.0;

        static_assert(nozzle_gamma == 1.33333, "nozzle_gamma (y) must be 1.33333 to use manual pow optimization");

        using lane = std::array<double, H>;
        using mask = std::array<int, H>;

        lane chamber_prev_volume_m3;
        lane chamber_volume_m3;
        lane chamber_nozzle_flow_area_m2;
        lane chamber_nozzle_open_ratio;
        lane chamber_static_pressure_pa;
        lane chamber_static_temperature_k;
        lane chamber_mass_kg;
        lane nozzle_mach;
        lane nozzle_velocity_m_per_s;
        lane nozzle_static_density_kg_per_m3;
        lane nozzle_mass_flow_rate_kg_per_s;
        lane parcel_mass_kg;
        lane parcel_static_temperature_k;
        mask panic;

        /*
         *     Ps * V
         * m = -------
         *     Rs * Ts
         */

        fn void calc_chamber_ambients()
        {
            for(int i = 0; i < H; i++)
            {
                chamber_static_temperature_k[i] = ambient_static_temperature_k;
                chamber_static_pressure_pa[i] = ambient_static_pressure_pa;
            }
            for(int i = 0; i < H; i++)
            {
                const double Ps = chamber_static_pressure_pa[i];
                const double V = chamber_volume_m3[i];
                const double Rs = specific_gas_constant_j_per_kg_k;
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
            for(int i = 0; i < N; i++)
            {
                const double m = chamber_mass_kg[i];
                const double Rs = specific_gas_constant_j_per_kg_k;
                const double Ts = chamber_static_temperature_k[i];
                const double V = chamber_volume_m3[i];
                chamber_static_pressure_pa[i] = m * Rs * Ts / V;
            }
        }

        /*
         * Assume upstream chamber total pressure (Pt) is equal to upstream chamber static pressure.
         * Downstream chamber static pressure (Ps) remains as-is.
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
            for(int i = 0; i < N; i++)
            {
                const int j = i + 1;
                const double X = (2.0 / (nozzle_gamma - 1.0));
                const double Pi = chamber_static_pressure_pa[i];
                const double Pj = chamber_static_pressure_pa[j];
                const double Pt = fmax(Pi, Pj);
                const double Ps = fmin(Pi, Pj);
                const double direction = Pi > Pj ? 1.0 : -1.0;

                /*
                 *     y - 1    1.333 - 1                        1     1          __________
                 *     ----- = ---------- = 0.25                --- * ---        /    ______
                 *       y       1.333                0.25       2     2        /    /
                 * Term                         = Term    = Term          = _  / _  /  Term
                 *                                                           \/   \/
                 */

                const double Y = sqrt(sqrt(Pt / Ps));
                const double M = direction * sqrt(X * (Y - 1.0));
                nozzle_mach[i] = clamp(M, -1.0, 1.0);
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
            for(int i = 0; i < N; i++)
            {
                const double Rs = specific_gas_constant_j_per_kg_k;
                const double Tt = chamber_static_temperature_k[i];
                const double M = nozzle_mach[i];
                const double X = nozzle_gamma * Rs * Tt;
                const double Y = 0.5 * (nozzle_gamma - 1.0) * M * M;
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
            for(int i = 0; i < N; i++)
            {
                const double Pt = chamber_static_pressure_pa[i];
                const double Rs = specific_gas_constant_j_per_kg_k;
                const double Tt = chamber_static_temperature_k[i];
                const double M = nozzle_mach[i];
                const double X = Pt / (Rs * Tt);
                const double C = 1.0 + 0.5 * (nozzle_gamma - 1.0) * M * M;

                /*
                 *        1         1
                 *      ----- = --------- = 3
                 *      y - 1   1.333 - 1
                 * Term                     = Term * Term * Term
                 *
                 */

                nozzle_static_density_kg_per_m3[i] = X / (C * C * C);
            }
        }

        /* .
         * m = ps A u
         *
         */

        fn void calc_nozzle_mass_flow_rates()
        {
            for(int i = 0; i < N; i++)
            {
                const double ps = nozzle_static_density_kg_per_m3[i];
                const double A = chamber_nozzle_open_ratio[i] * chamber_nozzle_flow_area_m2[i];
                const double u = nozzle_velocity_m_per_s[i];
                const double mdot = ps * A * u;
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
            for(int i = 0; i < N; i++)
            {
                const int j = i + 1;
                const double mdot = nozzle_mass_flow_rate_kg_per_s[i];
                const double dm = mdot * dt_s;
                parcel_mass_kg[i] = dm;
                const double Ti = chamber_static_temperature_k[i];
                const double Tj = chamber_static_temperature_k[j];
                parcel_static_temperature_k[i] = mdot > 0.0 ? Ti : Tj;
            }
        }

        /*                  y - 1
         *               V1
         * Ts2 = Ts1 * (----)
         *               V2
         */

        fn void calc_compressions()
        {
            for(int i = 0; i < N; i++)
            {
                const double Ts1 = chamber_static_temperature_k[i];
                const double V1 = chamber_prev_volume_m3[i];
                const double V2 = chamber_volume_m3[i];
                const double dv = V1 / V2;

                /*                                     _____
                 *                         1          /
                 *      y = 1.333 - 1.0 = --- = _ 3  / Term
                 * Term                    3     \  /
                 *                                \/
                 */

                chamber_static_temperature_k[i] = Ts1 * cbrt(dv);
            }
        }

        /*
         *       Ts m + Tp dm
         * Ts = --------------
         *          m + dm
         */

        fn void calc_forward_energy_transfers()
        {
            for(int i = 0; i < N; i++)
            {
                const int j = i + 1;
                const double dm = parcel_mass_kg[i];
                const double m = chamber_mass_kg[j];
                const double Tp = parcel_static_temperature_k[i];
                const double Ts0 = chamber_static_temperature_k[j];
                const double Ts1 = (Ts0 * m + Tp * dm) / (m + dm);
                chamber_static_temperature_k[j] = dm > 0.0 ? Ts1 : Ts0;
            }
        }

        /*                  y - 1
         *               m1
         * Ts2 = Ts1 * (----)
         *               m2
         */

        fn void calc_mass_transfers()
        {
            for(int i = 0; i < N; i++)
            {
                const int j = i + 1;
                const double mi = chamber_mass_kg[i];
                const double mj = chamber_mass_kg[j];
                const double dm = parcel_mass_kg[i];
                const double m0 = mi - dm;
                const double m1 = mj + dm;

                /*                                     _____
                 *                         1          /
                 *      y = 1.333 - 1.0 = --- = _ 3  / Term
                 * Term                    3     \  /
                 *                                \/
                 */

                chamber_static_temperature_k[i] *= cbrt(m0 / mi);
                chamber_static_temperature_k[j] *= cbrt(m1 / mj);
                chamber_mass_kg[i] = m0;
                chamber_mass_kg[j] = m1;
            }
        }

        fn void calc_panics()
        {
            for(int i = 0; i < N; i++)
            {
                panic[i] |= chamber_mass_kg[i] <= 0.0;
                panic[i] |= chamber_static_temperature_k[i] <= 0.0;
                panic[i] |= chamber_static_pressure_pa[i] <= 0.0;
            }
        }

        fn void calc()
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
        double mass_kg;
        double radius_m;
        double theta_r;
        double last_theta_r;
        double angular_velocity_r_per_s;
        double moment_of_inertia_kg_m2;

        /*
         * dw = dw/dt * dt
         *
         */

        fn void accelerate(const double angular_acceleration_r_per_s2)
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
            const double t0 = modulos(last_theta_r, otto_cycle_r);
            const double t1 = modulos(theta_r, otto_cycle_r);
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
            moment_of_inertia_kg_m2 = 0.5 * mass_kg * radius_m * radius_m;
        }

        fn bool calc(const double angular_acceleration_r_per_s2)
        {
            angular_velocity_r_per_s = 100.0;
            accelerate(angular_acceleration_r_per_s2);
            turn();
            const bool cycled = otto_cycled();
            calc_moment_of_inertia();
            return cycled;
        }
    };

    template<int W>
    struct sparkplugs
    {
        using lane = std::array<double, W>;
        using mask = std::array<int, W>;

        static constexpr double fire_delay_theta_r = 1e-1;
        lane engage_theta_r;
        mask prev_fired;
        mask fired;
        mask rising_edge;

        double crankshaft_theta_r;

        fn void calc_fired()
        {
            for(int i = 0; i < W; i++)
            {
                prev_fired[i] = fired[i];
                const double mod_theta_r = modulos(crankshaft_theta_r, otto_cycle_r);
                const double mod_engage_theta_r = engage_theta_r[i] >= otto_cycle_r
                    ? (engage_theta_r[i] - otto_cycle_r)
                    : (engage_theta_r[i]);
                const double mod_engage_delayed_theta_r = mod_engage_theta_r + fire_delay_theta_r;
                fired[i] = mod_theta_r > mod_engage_delayed_theta_r;
            }
        }

        fn void calc_rising_edge()
        {
            for(int i = 0; i < W; i++)
            {
                rising_edge[i] = fired[i] & !prev_fired[i];
            }
        }

        fn void calc()
        {
            calc_fired();
            calc_rising_edge();
        }
    };

    template<int W>
    struct simple_cam
    {
        using lane = std::array<double, W>;

        lane engage_theta_r;
        lane ramp_theta_r;
        lane open_ratio;

        double crankshaft_theta_r;

        /*         4       5       6       7
         * r = 35 t  - 84 t  + 70 t  - 20 t
         *
         * With friendly powers:
         *
         *      4            1      2      3
         * r = t  [ 35 - 84 t + 70 t - 20 t ]
         *
         */

        fn void calc_open_ratios()
        {
            for(int i = 0; i < W; i++)
            {
                double mod_engage_theta_r = modulos(engage_theta_r[i], otto_cycle_r);
                if(mod_engage_theta_r < 0.0)
                {
                    mod_engage_theta_r += otto_cycle_r;
                }
                double mod_theta_r = modulos(crankshaft_theta_r, otto_cycle_r);
                if(mod_theta_r < mod_engage_theta_r)
                {
                    mod_theta_r += otto_cycle_r;
                }
                const double open_r = mod_theta_r - mod_engage_theta_r;
                const double t = open_r / ramp_theta_r[i];
                const double a = t * t * t * t;
                const double b = t * a;
                const double c = t * b;
                const double d = t * c;
                const double A = 35.0 * a;
                const double B = 84.0 * b;
                const double C = 70.0 * c;
                const double D = 20.0 * d;
                const double R = clamp(A - B + C - D, 0.0, 1.0);
                open_ratio[i] = mod_theta_r < mod_engage_theta_r ? 0.0 : R;
            }
        }

        void calc()
        {
            calc_open_ratios();
        }
    };

    template<int W>
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

        using lane = std::array<double, W>;

        lane diameter_m;
        lane crank_throw_length_m;
        lane connecting_rod_length_m;
        lane connecting_rod_mass_kg;
        lane head_mass_density_kg_per_m3;
        lane head_compression_height_m;
        lane head_clearance_height_m;
        lane theta0_r;
        lane theta_r;
        lane sint;
        lane cost;
        lane pin_x_m;
        lane pin_y_m;
        lane bearing_x_m;
        lane bearing_y_m;
        lane volumes_m3;
        lane head_mass_kg;
        lane moment_of_inertia_kg_m2;
        lane gas_torque_n_m;
        lane inertia_torque_n_m;
        lane chamber_static_pressure_pa;

        double crankshaft_angular_velocity_r_per_s;
        double crankshaft_theta_r;

        /*
         * t = t0 + t1
         */

        fn void calc_thetas()
        {
            for(int i = 0; i < W; i++)
            {
                theta_r[i] = crankshaft_theta_r - theta0_r[i];
            }
        }

        fn void calc_sin_cos()
        {
            for(int i = 0; i < W; i++)
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
            for(int i = 0; i < W; i++)
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
            for(int i = 0; i < W; i++)
            {
                const double r = crank_throw_length_m[i];
                const double l = connecting_rod_length_m[i];
                const double cm = head_compression_height_m[i];
                const double cl = head_clearance_height_m[i];
                const double block_deck_surface_m = r + l + cm + cl;
                const double y = pin_y_m[i] + cm;
                const double radius = diameter_m[i] / 2.0;
                const double h = block_deck_surface_m - y;
                volumes_m3[i] = pi * radius * radius * h;
            }
        }

        /*           2
         * M = pi * r  * h * p
         *
         */

        fn void calc_masses()
        {
            for(int i = 0; i < W; i++)
            {
                const double r = 0.5 * diameter_m[i];
                const double h = 2.0 * head_compression_height_m[i];
                const double p = head_mass_density_kg_per_m3[i];
                head_mass_kg[i] = pi * r * r * h * p;
            }
        }

        /*             1        2
         * I = [ mp + --- mr ] r
         *             3
         */

        fn void calc_moments_of_inertia()
        {
            for(int i = 0; i < W; i++)
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
            for(int i = 0; i < W; i++)
            {
                const double Pg = chamber_static_pressure_pa[i];
                const double A = pi * diameter_m[i] * diameter_m[i];
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
            for(int i = 0; i < W; i++)
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

        fn void calc_volumetrics()
        {
            calc_thetas();
            calc_sin_cos();
            calc_positions();
            calc_volumes();
        }

        fn void calc()
        {
            calc_volumetrics();
            calc_masses();
            calc_moments_of_inertia();
            calc_gas_torques();
            calc_inertia_torques();
        }
    };

    template<typename... T>
    concept aggregate_of = (std::is_aggregate_v<T>&& ...);

    template<
        int W,
        int H,
        int PY,
        template<int> class P,
        template<int> class C,
        template<int> class S>
    requires(aggregate_of<crankshaft, P<W>, flow<H>, C<W>, S<W>>)
    struct as_engine : engine
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

        crankshaft crankshaft = {};
        P<W> pistons = {};
        C<W> inlet_cam = {};
        C<W> outlet_cam = {};
        S<W> sparkplugs = {};
        std::array<flow<H>, W> flow = {};
        diags back, front;

        void log_at(const unsigned x, const unsigned y)
        {
            if(x < W && y < H)
            {
                #define X(name) back[diags::channel::name].push_back(flow[x].name[y]);
                ENSIM_FLUIDS_LIST(X)
                #undef X
                if(y == PY)
                {
                    #define X(name) back[diags::channel::name].push_back(pistons.name[x]);
                    ENSIM_PISTONS_LIST(X)
                    #undef X
                }
            }
        }

        /*
         * State is broadcasted with a copy to prevent reference/pointer aliasing that
         * could interfere with vectorization.
         */

        fn void broadcast_states()
        {
            inlet_cam.crankshaft_theta_r = crankshaft.theta_r;
            outlet_cam.crankshaft_theta_r = crankshaft.theta_r;
            pistons.crankshaft_theta_r = crankshaft.theta_r;
            sparkplugs.crankshaft_theta_r = crankshaft.theta_r;
            pistons.crankshaft_angular_velocity_r_per_s = crankshaft.angular_velocity_r_per_s;
            for(int x = 0; x < W; x++)
            {
                flow[x].chamber_volume_m3[PY] = pistons.volumes_m3[x];
                flow[x].chamber_nozzle_open_ratio[PY - 1] = inlet_cam.open_ratio[x];
                flow[x].chamber_nozzle_open_ratio[PY + 0] = outlet_cam.open_ratio[x];
            }
            for(int x = 0; x < W; x++)
            {
                pistons.chamber_static_pressure_pa[x] = flow[x].chamber_static_pressure_pa[PY];
            }
        }

        void remember_volumes()
        {
            for(int x = 0; x < W; x++)
            {
                flow[x].chamber_prev_volume_m3 = flow[x].chamber_volume_m3;
            }
        }

        fn void reset_chambers()
        {
            for(int x = 0; x < W; x++)
            {
                flow[x].calc_chamber_ambients();
            }
        }

        void reset() override
        {
            pistons.calc_volumetrics();
            broadcast_states();
            remember_volumes();
            reset_chambers();
        }

        void run_crankshaft()
        {
            const bool otto_cycled = crankshaft.calc(0.0);
            if(otto_cycled)
            {
                std::swap(front, back);
                back.clear();
            }
        }

        void run_cams()
        {
            inlet_cam.calc();
            outlet_cam.calc();
        }

        void run_sparkplugs()
        {
            sparkplugs.calc();
        }

        void run_pistons()
        {
            pistons.calc();
        }

        void run_flows()
        {
            for(int x = 0; x < W; x++)
            {
                flow[x].calc();
            }
        }

        fn void run_injection()
        {
            for(int x = 0; x < W; x++)
            {
                if(sparkplugs.rising_edge[x])
                {
                    /* ignite */
                }
            }
        }

        void run(const unsigned steps, const unsigned x, const unsigned y) override
        {
            for(unsigned i = 0; i < steps; i++)
            {
                run_crankshaft();
                run_cams();
                run_sparkplugs();
                run_pistons();
                run_injection();
                run_flows();
                log_at(x, y);
                remember_volumes();
                broadcast_states();
            }
        }

        int width() override
        {
            return W;
        }

        int height() override
        {
            return H;
        }

        int piston_y() override
        {
            return PY;
        }

        size_t bytes() override
        {
            return sizeof *this;
        }

        const diags& get_diags() const override
        {
            return front;
        }

        std::vector<std::vector<double>> new_matrix(const int width, const int height)
        {
            std::vector<std::vector<double>> matrix;
            matrix.resize(height);
            for(auto& row : matrix)
            {
                row.resize(width);
            }
            return matrix;
        }

        std::vector<std::vector<double>> get_panics() override
        {
            std::vector<std::vector<double>> matrix = new_matrix(width(), height());
            for(int y = 0; y < H; y++)
            for(int x = 0; x < W; x++)
            {
                matrix[y][x] = flow[x].panic[y];
            }
            return matrix;
        }

        std::vector<std::vector<double>> get_port_open_ratios() override
        {
            std::vector<std::vector<double>> matrix = new_matrix(width(), height());
            for(int y = 0; y < H; y++)
            for(int x = 0; x < W; x++)
            {
                matrix[y][x] = flow[x].chamber_nozzle_open_ratio[y];
            }
            return matrix;
        }
    };

    struct inline8 : as_engine<8, 9, 4, inline_pistons, simple_cam, sparkplugs>
    {
        inline8()
        {
            this->crankshaft.mass_kg = 25.3;
            this->crankshaft.radius_m = 0.031;
            this->pistons.diameter_m.fill(0.085);
            this->pistons.crank_throw_length_m.fill(0.038);
            this->pistons.connecting_rod_length_m.fill(0.1);
            this->pistons.connecting_rod_mass_kg.fill(0.4);
            this->pistons.head_mass_density_kg_per_m3.fill(7800);
            this->pistons.head_compression_height_m.fill(0.025);
            this->pistons.head_clearance_height_m.fill(0.007);
            inlet_cam.ramp_theta_r.fill(pi);
            outlet_cam.ramp_theta_r.fill(pi);
            double theta0_r = 0.0;
            for(int i = 0; i < width(); i++)
            {
                this->pistons.theta0_r[i] = theta0_r;
                this->inlet_cam.engage_theta_r[i] = theta0_r + otto_intake_cycle_r;
                this->sparkplugs.engage_theta_r[i] = theta0_r + otto_combustion_cycle_r;
                this->outlet_cam.engage_theta_r[i] = theta0_r + otto_exhaust_cycle_r;
                theta0_r += otto_cycle_r / width();
            }
            for(auto& flow : this->flow)
            {
                flow.chamber_nozzle_open_ratio.fill(1.0);
                flow.chamber_volume_m3 = {
                    resevoir_volume_m3,
                    0.05,
                    0.05,
                    0.05,
                    0.05,
                    0.05,
                    0.05,
                    0.05,
                    resevoir_volume_m3,
                };
                flow.chamber_nozzle_flow_area_m2 = {
                    6e-3,
                    5e-3,
                    4e-3,
                    3e-3,
                    5e-3,
                    2e-3,
                    3e-3,
                    4e-3,
                };
            }
        }
    };

    std::unique_ptr<engine> new_engine(const engine::type type)
    {
        std::unique_ptr<engine> engine;
        switch(type)
        {
            default:
            case engine::type::inline8: engine = std::make_unique<inline8>(); break;
        }
        engine->reset();
        return engine;
    }
}
