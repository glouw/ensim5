#include "ensim5.hh"

#include <cmath>
#include <array>

#define fn __attribute__((used))

namespace ensim
{
    using std::sin, std::cos, std::sqrt, std::fmax, std::fmin, std::fmod, std::cbrt;

    void diags::clear()
    {
        for(auto& plot : plots)
        {
            plot.clear();
        }
    }

    std::vector<float>& diags::operator[](const int index)
    {
        return plots[index];
    }

    std::vector<float>& diags::operator[](const channel channel)
    {
        return plots[static_cast<int>(channel)];
    }

    template<std::floating_point T>
    class real
    {
        T value = {};
    public:
        constexpr real() = default;
        constexpr real(const double value): value(static_cast<T>(value)) {}
        constexpr real& operator-=(const real& other) { value -= other.value; return *this; }
        constexpr real& operator+=(const real& other) { value += other.value; return *this; }
        constexpr real& operator*=(const real& other) { value *= other.value; return *this; }
        constexpr operator T() const { return value; }
    };

    using val = real<double>;

    static constexpr int sample_rate_hz = 48000;
    static constexpr val dt_s = val(1.0) / sample_rate_hz;
    static constexpr val ambient_static_temperature_k = 293.0;
    static constexpr val ambient_static_pressure_pa = 101325.0;
    static constexpr val otto_cycle_r = 4.0 * M_PI;
    static constexpr val resevoir_volume_m3 = 1e9;

    consteval bool is_odd(const int value) { return value % 2 == 1; }
    consteval bool is_evn(const int value) { return value % 2 == 0; }

    template<int H>
    requires(is_odd(H))
    struct flow
    {
        static constexpr int N = H - 1;
        static_assert(is_evn(N));

        static constexpr val nozzle_gamma = 1.333;
        static constexpr val universal_gas_constant_j_per_mol_k = 8.3144598;
        static constexpr val cv_j_per_mol_k = universal_gas_constant_j_per_mol_k / (nozzle_gamma - 1.0);
        static constexpr val molar_mass_kg_per_mol = 0.023;
        static constexpr val specific_gas_constant_j_per_kg_k = universal_gas_constant_j_per_mol_k / molar_mass_kg_per_mol;

        using lane_t = std::array<val, H>;

        lane_t chamber_prev_volume_m3;
        lane_t chamber_volume_m3;
        lane_t chamber_nozzle_flow_area_m2;
        lane_t chamber_static_pressure_pa;
        lane_t chamber_static_temperature_k;
        lane_t chamber_mass_kg;
        lane_t nozzle_mach;
        lane_t nozzle_velocity_m_per_s;
        lane_t nozzle_static_density_kg_per_m3;
        lane_t nozzle_mass_flow_rate_kg_per_s;
        lane_t parcel_mass_kg;
        lane_t parcel_static_temperature_k;

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
                const val Ps = chamber_static_pressure_pa[i];
                const val V = chamber_volume_m3[i];
                const val Rs = specific_gas_constant_j_per_kg_k;
                const val Ts = chamber_static_temperature_k[i];
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
            for(int i = 0; i < N; i++)
            {
                const val m = chamber_mass_kg[i];
                const val Rs = specific_gas_constant_j_per_kg_k;
                const val Ts = chamber_static_temperature_k[i];
                const val V = chamber_volume_m3[i];
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
                const val X = (val(2.0) / (nozzle_gamma - val(1.0)));
                const val Pi = chamber_static_pressure_pa[i];
                const val Pj = chamber_static_pressure_pa[j];
                const val Pt = fmax(Pi, Pj);
                const val Ps = fmin(Pi, Pj);
                const val direction = Pi > Pj ? val(1.0) : val(-1.0);

                /*
                 *      y - 1    1.333 - 1                        1     1          __________
                 *      ----- = ---------- = 0.25                --- * ---        /    ______
                 *        y       1.333                0.25       2     2        /    /
                 *  term                         = term    = term          = _  / _  /  term
                 *                                                            \/   \/
                 */

                static_assert(nozzle_gamma == val(1.333), "nozzle_gamma (y) must be 1.333 to use manual pow optimization");
                const val Y = sqrt(sqrt(Pt / Ps));
                const val M0 = direction * sqrt(X * (Y - val(1.0)));
                const val M1 = fmax(M0, val(-1.0));
                const val M2 = fmin(M1, val(+1.0));
                nozzle_mach[i] = M2;
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
                const val Rs = specific_gas_constant_j_per_kg_k;
                const val Tt = chamber_static_temperature_k[i];
                const val M = nozzle_mach[i];
                const val X = nozzle_gamma * Rs * Tt;
                const val Y = val(0.5) * (nozzle_gamma - val(1.0)) * M * M;
                const val u = M * sqrt(X / (val(1.0) + Y));
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
                const val Pt = chamber_static_pressure_pa[i];
                const val Rs = specific_gas_constant_j_per_kg_k;
                const val Tt = chamber_static_temperature_k[i];
                const val M = nozzle_mach[i];
                const val X = Pt / (Rs * Tt);
                const val C = val(1.0) + val(0.5) * (nozzle_gamma - val(1.0)) * M * M;

                /*
                 *         1         1
                 *       ----- = --------- = 3.00
                 *       y - 1   1.333 - 1
                 *  term                         = term * term * term
                 */

                static_assert(nozzle_gamma == val(1.333), "nozzle_gamma (y) must be 1.333 to use manual pow optimization");
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
                const val ps = nozzle_static_density_kg_per_m3[i];
                const val A = chamber_nozzle_flow_area_m2[i];
                const val u = nozzle_velocity_m_per_s[i];
                const val mdot = ps * A * u;
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
                const val mdot = nozzle_mass_flow_rate_kg_per_s[i];
                parcel_mass_kg[i] = mdot * dt_s;

                /*
                 * Since this is an internal combustion engine and not a rocket engine,
                 * discard nozzle enthalpy and assume parcel temperature is equal to upstream
                 * chamber static temperature.
                 *
                 */

                const val Ti = chamber_static_temperature_k[i];
                const val Tj = chamber_static_temperature_k[j];
                parcel_static_temperature_k[i] = mdot > val(0.0) ? Ti : Tj;
            }
        }

        fn void calc_chamber_temperatures_from_compression()
        {
            /*                  y - 1
             *               V1
             * Ts2 = Ts1 * (----)
             *               V2
             */

            for(int i = 0; i < N; i++)
            {
                const val V1 = chamber_volume_m3[i];
                const val V2 = chamber_prev_volume_m3[i];
                const val dv = V1 / V2;

                /*                                     _____
                 *                         1          /
                 *      y = 1.333 - 1.0 = --- = _ 3  / term
                 * term                    3     \  /
                 *                                \/
                 */

                static_assert(nozzle_gamma == val(1.333), "nozzle_gamma (y) must be 1.333 to use manual pow optimization");
                chamber_static_temperature_k[i] *= cbrt(dv);
            }
        }

        /*
         * ---> is positive (+) convention
         *
         * +----+      +-----+       +----+     +----+
         * |    | m1   | m2  |  m3   |    |     |    |
         * |    | ---> |     |  ---> |    | ... |    |
         * |    | Ts1  | Ts2 |  Ts2  |    |     |    |
         * +----+      +-----+       +----+     +----+
         *
         *        m2 Ts2 + m1 Ts1 - m3 Ts3
         * Ts2 = --------------------------
         *              m2 + m1 - m3
         *
         */

        fn void calc_chamber_temperatures_from_mixing()
        {

            for(int j = 1; j < N; j++)
            {
                const int i = j - 1;
                const int k = j + 1;
                const val m1 = parcel_mass_kg[i];
                const val m2 = chamber_mass_kg[j];
                const val m3 = parcel_mass_kg[k];
                const val Ts1 = parcel_static_temperature_k[i];
                const val Ts2 = chamber_static_temperature_k[j];
                const val Ts3 = parcel_static_temperature_k[k];
                chamber_static_temperature_k[j] = (m2 * Ts2 + m1 * Ts1 - m3 * Ts3) / (m2 + m1 - m3);
            }
        }

        fn void calc_chamber_masses_from_transfer()
        {
            for(int i = 0; i < N; i++)
            {
                const float m = parcel_mass_kg[i];
                const int j = i + 1;
                chamber_mass_kg[i] -= m;
                chamber_mass_kg[j] += m;
            }

            /*
             * This safe guard effectively injects mass from outside the simulation
             * control volume to prevent chamber mass from approaching 0.0 and
             * introducing NaNs downstream.
             *
             */

            for(int i = 0; i < N; i++)
            {
                chamber_mass_kg[i] = fmax(chamber_mass_kg[i], 1e-5);
            }
        }

        fn void calc()
        {
            calc_chamber_pressures();
            calc_nozzle_machs();
            calc_nozzle_velocities();
            calc_nozzle_static_densities();
            calc_nozzle_mass_flow_rates();
            calc_nozzle_parcels();
            calc_chamber_temperatures_from_compression();
            calc_chamber_temperatures_from_mixing();
            calc_chamber_masses_from_transfer();

            /*
             * Nozzle equations assume upstream total conditions are equal to upstream static conditions.
             * Chamber momentum calculations can therefor be skipped.
             *
             */
        }
    };

    struct crankshaft
    {
        val mass_kg;
        val radius_m;
        val theta_r;
        val last_theta_r;
        val angular_velocity_r_per_s;
        val moment_of_inertia_kg_m2;

        /*
         * dw = dw/dt * dt
         *
         */

        fn void accelerate(const val angular_acceleration_r_per_s2)
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
            const val t0 = fmod(last_theta_r, otto_cycle_r);
            const val t1 = fmod(theta_r, otto_cycle_r);
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
            moment_of_inertia_kg_m2 = val(0.5) * mass_kg * radius_m * radius_m;
        }

        fn bool calc(const val angular_acceleration_r_per_s2)
        {
            accelerate(angular_acceleration_r_per_s2);
            turn();
            const bool cycled = otto_cycled();
            calc_moment_of_inertia();
            return cycled;
        }
    };

    template<int W, int PY>
    struct pistons
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

        using lane_t = std::array<val, W>;

        lane_t diameter_m;
        lane_t crank_throw_length_m;
        lane_t connecting_rod_length_m;
        lane_t connecting_rod_mass_kg;
        lane_t head_mass_density_kg_per_m3;
        lane_t head_compression_height_m;
        lane_t head_clearance_height_m;
        lane_t theta0_r;
        lane_t theta_r;
        lane_t sint;
        lane_t cost;
        lane_t pin_x_m;
        lane_t pin_y_m;
        lane_t bearing_x_m;
        lane_t bearing_y_m;
        lane_t volumes_m3;
        lane_t head_mass_kg;
        lane_t moment_of_inertia_kg_m2;
        lane_t gas_torque_n_m;
        lane_t inertia_torque_n_m;

        /*
         * t = t0 + t1
         */

        fn void calc_thetas(const crankshaft& __restrict crankshaft)
        {
            for(int i = 0; i < W; i++)
            {
                theta_r[i] = theta0_r[i] + crankshaft.theta_r;
            }
        }

        fn void calc_sin_cos()
        {
            for(int i = 0; i < W; i++)
            {
                const val t = theta_r[i];
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
                const val r = crank_throw_length_m[i];
                const val l = connecting_rod_length_m[i];
                const val t = theta_r[i];
                const val x = r * sint[i];
                const val y = r * cost[i];
                bearing_x_m[i] = x;
                bearing_y_m[i] = y;
                pin_x_m[i] = val(0.0);
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
                const val r = crank_throw_length_m[i];
                const val l = connecting_rod_length_m[i];
                const val cm = head_compression_height_m[i];
                const val cl = head_clearance_height_m[i];
                const val block_deck_surface_m = r + l + cm + cl;
                const val y = pin_y_m[i] + cm;
                const val radius = diameter_m[i] / val(2.0);
                const val h = block_deck_surface_m - y;
                volumes_m3[i] = val(M_PI) * radius * radius * h;
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
                const val r = val(0.5) * diameter_m[i];
                const val h = val(2.0) * head_compression_height_m[i];
                const val p = head_mass_density_kg_per_m3[i];
                head_mass_kg[i] = val(M_PI) * r * r * h * p;
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
                const val r = crank_throw_length_m[i];
                const val mp = head_mass_kg[i];
                const val mr = connecting_rod_mass_kg[i];
                moment_of_inertia_kg_m2[i] = (mp + (val(1.0 / 3.0) * mr)) * r * r;
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
            for(int i = 0; i < W; i++)
            {
                const val Pg = flow.chamber_static_pressure_pa[PY];
                const val A = val(M_PI) * diameter_m[i] * diameter_m[i];
                const val r = crank_throw_length_m[i];
                const val l = connecting_rod_length_m[i];
                const val t = theta_r[i];
                const val X = Pg * A * r * sint[i];
                const val Y = val(1.0) + (r / l) * cost[i];
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

        fn void calc_inertia_torques(const crankshaft& __restrict crankshaft)
        {
            for(int i = 0; i < W; i++)
            {
                const val r = crank_throw_length_m[i];
                const val l = connecting_rod_length_m[i];
                const val I = moment_of_inertia_kg_m2[i];
                const val w = crankshaft.angular_velocity_r_per_s;
                const val t = theta_r[i];
                const val rl = r / l;
                const val s = sint[i];
                const val c = cost[i];
                const val X = val(0.25) * rl * s;
                const val Y = s * c;
                const val Z = val(0.75) * rl * (val(3.0) * s - val(4.0) * s * s * s);
                inertia_torque_n_m[i] = I * w * w * (X - Y - Z);
            }
        }

        fn void calc_volumetrics(const crankshaft& __restrict crankshaft)
        {
            calc_thetas(crankshaft);
            calc_sin_cos();
            calc_positions();
            calc_volumes();
        }

        fn void calc(const auto& __restrict flow, const crankshaft& __restrict crankshaft)
        {
            calc_volumetrics(crankshaft);
            calc_masses();
            calc_moments_of_inertia();
            calc_gas_torques(flow);
            calc_inertia_torques(crankshaft);
        }
    };

    template<typename... T>
    concept aggregate_of = (std::is_aggregate_v<T> && ...);

    template<int W, int H, int PY>
    requires(aggregate_of<crankshaft, pistons<W, PY>, flow<H>>)
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
        pistons<W, PY> pistons = {};
        std::array<flow<H>, W> flow = {};
        diags back;
        diags front;

        void log(const int x, const int y)
        {
            back[diags::channel::chamber_volume_m3              ].push_back(flow[x].chamber_volume_m3              [y]);
            back[diags::channel::chamber_nozzle_flow_area_m2    ].push_back(flow[x].chamber_nozzle_flow_area_m2    [y]);
            back[diags::channel::chamber_static_pressure_pa     ].push_back(flow[x].chamber_static_pressure_pa     [y]);
            back[diags::channel::chamber_static_temperature_k   ].push_back(flow[x].chamber_static_temperature_k   [y]);
            back[diags::channel::chamber_mass_kg                ].push_back(flow[x].chamber_mass_kg                [y]);
            back[diags::channel::nozzle_mach                    ].push_back(flow[x].nozzle_mach                    [y]);
            back[diags::channel::nozzle_velocity_m_per_s        ].push_back(flow[x].nozzle_velocity_m_per_s        [y]);
            back[diags::channel::nozzle_static_density_kg_per_m3].push_back(flow[x].nozzle_static_density_kg_per_m3[y]);
            back[diags::channel::nozzle_mass_flow_rate_kg_per_s ].push_back(flow[x].nozzle_mass_flow_rate_kg_per_s [y]);
        }

        void relay_volume()
        {
            for(int x = 0; x < W; x++)
            {
                flow[x].chamber_volume_m3[PY] = pistons.volumes_m3[x];
            }
        }

        void remember_volume()
        {
            for(int x = 0; x < W; x++)
            {
                flow[x].chamber_prev_volume_m3 = flow[x].chamber_volume_m3;
            }
        }

        void reset() override
        {
            pistons.calc_volumetrics(crankshaft);
            relay_volume();
            remember_volume();
            for(int x = 0; x < W; x++)
            {
                flow[x].calc_chamber_ambients();
            }
        }

        bool in_bounds(const int x, const int y)
        {
            return x >= 0 && y >= 0 && x < W && y < H;
        }

        void run(const int steps, const int x, const int y) override
        {
            for(int i = 0; i < steps; i++)
            {
                remember_volume();
                crankshaft.angular_velocity_r_per_s = 500.0;
                const bool otto_cycled = crankshaft.calc(0.0);
                if(otto_cycled)
                {
                    std::swap(front, back);
                    back.clear();
                }
                for(int j = 0; j < W; j++)
                {
                    pistons.calc(flow[j], crankshaft);
                }
                relay_volume();
                for(int j = 0; j < W; j++)
                {
                    flow[j].calc();
                }
                if(in_bounds(x, y))
                {
                    log(x, y);
                }
            }
        }

        int get_w() override
        {
            return W;
        }

        int get_h() override
        {
            return H;
        }

        int get_y() override
        {
            return PY;
        }

        diags& get_diags() override
        {
            return front;
        }

        size_t get_size() override
        {
            return sizeof *this;
        }
    };

    struct inline8 : as_engine<8, 9, 4>
    {
        inline8()
        {
            this->crankshaft.mass_kg  = 0.1;
            this->crankshaft.radius_m = 0.25;
            this->pistons.diameter_m.fill(0.1);
            this->pistons.crank_throw_length_m.fill(0.1);
            this->pistons.connecting_rod_length_m.fill(0.1);
            this->pistons.connecting_rod_mass_kg.fill(0.1);
            this->pistons.head_mass_density_kg_per_m3.fill(0.1);
            this->pistons.head_compression_height_m.fill(0.1);
            this->pistons.head_clearance_height_m.fill(0.01);
            for(auto& flow : this->flow)
            {
                flow.chamber_volume_m3 = {
                    resevoir_volume_m3,
                    0.10,
                    0.10,
                    0.10,
                    0.10,
                    0.10,
                    0.10,
                    0.10,
                    resevoir_volume_m3,
                };
                flow.chamber_nozzle_flow_area_m2 = {
                    0.01,
                    0.01,
                    0.01,
                    0.01,
                    0.01,
                    0.01,
                    0.01,
                    0.01,
                    0.01,
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
