#include "ensim5.hh"

#include <cmath>
#include <array>

#define fn __attribute__((used))

using std::sin, std::cos, std::sqrt, std::fmax, std::fmin, std::fmod, std::cbrt;

template<std::floating_point T>
class real_c
{
    T value = {};
public:
    constexpr real_c() = default;
    constexpr real_c(const double value): value(static_cast<T>(value)) {}
    constexpr real_c& operator-=(const real_c other) { value -= other.value; return *this; }
    constexpr real_c& operator+=(const real_c other) { value += other.value; return *this; }
    constexpr real_c& operator*=(const real_c other) { value *= other.value; return *this; }
    constexpr operator T() const { return value; }
};

using val_t = real_c<float>;

static constexpr int g_sample_rate_hz = 48000;
static constexpr val_t g_dt_s = val_t(1.0) / g_sample_rate_hz;
static constexpr val_t g_ambient_static_temperature_k = 293.0;
static constexpr val_t g_ambient_static_pressure_pa = 101325.0;
static constexpr val_t g_otto_cycle_r = 4.0 * M_PI;

consteval bool is_odd(const int value) { return value % 2 == 1; }
consteval bool is_evn(const int value) { return value % 2 == 0; }

template<int H>
requires(is_odd(H))
struct flow_s
{
    static constexpr int N = H - 1;
    static_assert(is_evn(N));

    static constexpr val_t nozzle_gamma = 1.333;
    static constexpr val_t universal_gas_constant_j_per_mol_k = 8.3144598;
    static constexpr val_t cv_j_per_mol_k = universal_gas_constant_j_per_mol_k / (nozzle_gamma - 1.0);
    static constexpr val_t molar_mass_kg_per_mol = 0.023;
    static constexpr val_t specific_gas_constant_j_per_kg_k = universal_gas_constant_j_per_mol_k / molar_mass_kg_per_mol;

    using lane_t = std::array<val_t, H>;

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
            chamber_static_temperature_k[i] = g_ambient_static_temperature_k;
            chamber_static_pressure_pa[i] = g_ambient_static_pressure_pa;
        }
        for(int i = 0; i < H; i++)
        {
            const val_t Ps = chamber_static_pressure_pa[i];
            const val_t V = chamber_volume_m3[i];
            const val_t Rs = specific_gas_constant_j_per_kg_k;
            const val_t Ts = chamber_static_temperature_k[i];
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
            const val_t m = chamber_mass_kg[i];
            const val_t Rs = specific_gas_constant_j_per_kg_k;
            const val_t Ts = chamber_static_temperature_k[i];
            const val_t V = chamber_volume_m3[i];
            chamber_static_pressure_pa[i] = m * Rs * Ts / V;
        }
    }

    /*
     * Ignore dynamic conditions and assume Pt is upstream chamber static
     * pressure and Ps is downstream chamber static pressure.
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
            const val_t X = (val_t(2.0) / (nozzle_gamma - val_t(1.0)));
            const val_t Px = chamber_static_pressure_pa[i];
            const val_t Py = chamber_static_pressure_pa[j];
            const val_t Pt = fmax(Px, Py);
            const val_t Ps = fmin(Px, Py);
            const val_t direction = Px > Py ? val_t(1.0) : val_t(-1.0);

            /*
             *      y - 1    1.333 - 1                        1     1          __________
             *      ----- = ---------- = 0.25                --- * ---        /    ______
             *        y       1.333                0.25       2     2        /    /
             *  term                         = term    = term          = _  / _  /  term
             *                                                            \/   \/
             */

            static_assert(nozzle_gamma == val_t(1.333), "nozzle_gamma (y) must be 1.333 to use manual pow optimization");
            const val_t Y = sqrt(sqrt(Pt / Ps));
            const val_t M0 = direction * sqrt(X * (Y - val_t(1.0)));
            const val_t M1 = fmax(M0, val_t(-1.0));
            const val_t M2 = fmin(M1, val_t(+1.0));
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
            const val_t Rs = specific_gas_constant_j_per_kg_k;
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
        for(int i = 0; i < N; i++)
        {
            const val_t Pt = chamber_static_pressure_pa[i];
            const val_t Rs = specific_gas_constant_j_per_kg_k;
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

    /* .
     * m = ps A u
     *
     */

    fn void calc_nozzle_mass_flow_rates()
    {
        for(int i = 0; i < N; i++)
        {
            const val_t ps = nozzle_static_density_kg_per_m3[i];
            const val_t A = chamber_nozzle_flow_area_m2[i];
            const val_t u = nozzle_velocity_m_per_s[i];
            const val_t mdot = ps * A * u;
            nozzle_mass_flow_rate_kg_per_s[i] = mdot;
        }
    }

    fn void calc_nozzle_parcels()
    {
        for(int i = 0; i < N; i++)
        {
            parcel_mass_kg[i] = nozzle_mass_flow_rate_kg_per_s[i] * g_dt_s;
        }
        for(int i = 0; i < N; i++)
        {
            const int j = i + 1;
            const val_t Px = chamber_static_pressure_pa[i];
            const val_t Py = chamber_static_pressure_pa[j];
            const val_t Tx = chamber_static_temperature_k[i];
            const val_t Ty = chamber_static_temperature_k[j];
            parcel_static_temperature_k[i] = Px > Py ? Tx : Ty;
        }
    }

    fn void calc_chamber_temperatures()
    {
        /*
         * 1. Adiabatic Compression.
         *
         *                   y - 1
         *               V1
         * Ts2 = Ts1 * (----)
         *               V2
         */

        for(int i = 0; i < N; i++)
        {
            const val_t V1 = chamber_volume_m3[i];
            const val_t V2 = chamber_prev_volume_m3[i];
            const val_t dv = V1 / V2;

            /*                                     _____
             *                         1          /
             *      y = 1.333 - 1.0 = --- = _ 3  / term
             * term                    3     \  /
             *                                \/
             */

            static_assert(nozzle_gamma == val_t(1.333), "nozzle_gamma (y) must be 1.333 to use manual pow optimization");
            chamber_static_temperature_k[i] *= cbrt(dv);
        }

        /*
         * 2. Gas temperature mix (Simplified energy term).
         *
         *        (m2 - m3) * Ts2 + m1 * Ts1
         * Ts2 = ----------------------------
         *              m2 + m1 - m3
         */

        /*
         * +----+       +-----+       +----+ +----+
         * |    |  m1   | m2  |  m3   |    | |    |
         * |    | <---> |     | <---> |    | |    |
         * |    |  Ts1  | Ts2 |  Ts2  |    | |    |
         * +----+       +-----+       +----+ +----+
         *
         */

        for(int j = 1; j < N; j++)
        {
            const int i = j - 1;
            const int k = j + 1;
            const val_t m1 = parcel_mass_kg[i];
            const val_t m2 = chamber_mass_kg[j];
            const val_t m3 = parcel_mass_kg[k];
            const val_t Ts1 = parcel_static_temperature_k[i];
            const val_t Ts2 = chamber_static_temperature_k[j];
            const val_t Ts3 = parcel_static_temperature_k[k];
            chamber_static_temperature_k[j] = ((m2 - m3) * Ts2 + m1 * Ts1) / (m2 + m1 - m3);
        }

        /* 3. Gas mass transfer (Simplified conservation of mass). */

        for(int i = 0; i < N; i++)
        {
            const float m = parcel_mass_kg[i];
            const int j = i + 1;
            chamber_mass_kg[i] -= m;
            chamber_mass_kg[j] += m;
        }

        /* 4. Skip momentum since static conditions largely dominate dynamic conditions.
         *
         * */
    }

    fn void calc_flow()
    {
        calc_chamber_pressures();
        calc_nozzle_machs();
        calc_nozzle_velocities();
        calc_nozzle_static_densities();
        calc_nozzle_mass_flow_rates();
        calc_nozzle_parcels();
        calc_chamber_temperatures();
    }
};

struct crankshaft_s
{
    val_t mass_kg;
    val_t radius_m;
    val_t theta_r;
    val_t last_theta_r;
    val_t angular_velocity_r_per_s;
    val_t moment_of_inertia_kg_m2;

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
        last_theta_r = theta_r;
        theta_r += angular_velocity_r_per_s * g_dt_s;
    }

    fn bool otto_cycled()
    {
        const val_t t0 = fmod(last_theta_r, g_otto_cycle_r);
        const val_t t1 = fmod(theta_r, g_otto_cycle_r);
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
        moment_of_inertia_kg_m2 = val_t(0.5) * mass_kg * radius_m * radius_m;
    }

    fn bool calc_crankshaft(const val_t angular_acceleration_r_per_s2)
    {
        accelerate(angular_acceleration_r_per_s2);
        turn();
        const bool cycled = otto_cycled();
        calc_moment_of_inertia();
        return cycled;
    }
};

template<int W, int PY>
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

    using lane_t = std::array<val_t, W>;

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

    fn void calc_thetas(const crankshaft_s& __restrict crankshaft)
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
            const val_t t = theta_r[i];
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
            const val_t r = crank_throw_length_m[i];
            const val_t l = connecting_rod_length_m[i];
            const val_t t = theta_r[i];
            const val_t x = r * sint[i];
            const val_t y = r * cost[i];
            bearing_x_m[i] = x;
            bearing_y_m[i] = y;
            pin_x_m[i] = val_t(0.0);
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
        for(int i = 0; i < W; i++)
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
        for(int i = 0; i < W; i++)
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
        for(int i = 0; i < W; i++)
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
        for(int i = 0; i < W; i++)
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

template<typename... T>
concept aggregate_of = (std::is_aggregate_v<T> && ...);

template<int W, int H, int PY>
requires(aggregate_of<crankshaft_s, pistons_s<W, PY>, flow_s<H>>)
struct engine_s : public ensim5_s
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

    crankshaft_s crankshaft = {};
    pistons_s<W, PY> pistons = {};
    std::array<flow_s<H>, W> flow = {};
    ensim5_diag_s back, front;

    fn void log_engine(const int x, const int y)
    {
        back[ensim5_diag_e::chamber_volume_m3              ].push_back(flow[x].chamber_volume_m3              [y]);
        back[ensim5_diag_e::chamber_nozzle_flow_area_m2    ].push_back(flow[x].chamber_nozzle_flow_area_m2    [y]);
        back[ensim5_diag_e::chamber_static_pressure_pa     ].push_back(flow[x].chamber_static_pressure_pa     [y]);
        back[ensim5_diag_e::chamber_static_temperature_k   ].push_back(flow[x].chamber_static_temperature_k   [y]);
        back[ensim5_diag_e::chamber_mass_kg                ].push_back(flow[x].chamber_mass_kg                [y]);
        back[ensim5_diag_e::nozzle_mach                    ].push_back(flow[x].nozzle_mach                    [y]);
        back[ensim5_diag_e::nozzle_velocity_m_per_s        ].push_back(flow[x].nozzle_velocity_m_per_s        [y]);
        back[ensim5_diag_e::nozzle_static_density_kg_per_m3].push_back(flow[x].nozzle_static_density_kg_per_m3[y]);
        back[ensim5_diag_e::nozzle_mass_flow_rate_kg_per_s ].push_back(flow[x].nozzle_mass_flow_rate_kg_per_s [y]);
    }

    fn void relay_volume()
    {
        for(int x = 0; x < W; x++)
        {
            flow[x].chamber_volume_m3[PY] = pistons.volumes_m3[x];
        }
    }

    fn void remember_volume()
    {
        for(int x = 0; x < W; x++)
        {
            flow[x].chamber_prev_volume_m3 = flow[x].chamber_volume_m3;
        }
    }

    fn void reset_engine() override
    {
        pistons.calc_volumetrics(crankshaft);
        relay_volume();
        remember_volume();
        for(int x = 0; x < W; x++)
        {
            flow[x].calc_chamber_ambients();
        }
    }

    fn bool in_bounds(const int x, const int y)
    {
        return x >= 0 && y >= 0 && x < W && y < H;
    }

    fn void run_engine(const int steps, const int x, const int y) override
    {
        for(int i = 0; i < steps; i++)
        {
            crankshaft.angular_velocity_r_per_s = 500.0;
            const bool otto_cycled = crankshaft.calc_crankshaft(0.0);
            if(otto_cycled)
            {
                std::swap(front, back);
                back.clear();
            }
            for(int j = 0; j < W; j++)
            {
                pistons.calc_pistons(flow[j], crankshaft);
            }
            remember_volume();
            relay_volume();
            for(int j = 0; j < W; j++)
            {
                flow[j].calc_flow();
            }
            if(in_bounds(x, y))
            {
                log_engine(x, y);
            }
        }
    }

    fn int get_w() override { return W;  }
    fn int get_h() override { return H;  }
    fn int get_y() override { return PY; }
    fn ensim5_diag_s& get_diags() override { return front; }
};

std::unique_ptr<ensim5_s> new_ensim5_inline_8()
{
    auto engine = std::make_unique<engine_s<8, 9, 4>>();
    engine->crankshaft.mass_kg  = 0.1;
    engine->crankshaft.radius_m = 0.25;
    engine->pistons.diameter_m                  = { 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    engine->pistons.crank_throw_length_m        = { 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    engine->pistons.connecting_rod_length_m     = { 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    engine->pistons.connecting_rod_mass_kg      = { 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    engine->pistons.head_mass_density_kg_per_m3 = { 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    engine->pistons.head_compression_height_m   = { 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    engine->pistons.head_clearance_height_m     = { 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    engine->pistons.theta0_r                    = { 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 };
    engine->flow[0].chamber_volume_m3           = { 1000, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1000 };
    engine->flow[1].chamber_volume_m3           = { 1000, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1000 };
    engine->flow[2].chamber_volume_m3           = { 1000, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1000 };
    engine->flow[3].chamber_volume_m3           = { 1000, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1000 };
    engine->flow[4].chamber_volume_m3           = { 1000, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1000 };
    engine->flow[5].chamber_volume_m3           = { 1000, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1000 };
    engine->flow[6].chamber_volume_m3           = { 1000, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1000 };
    engine->flow[7].chamber_volume_m3           = { 1000, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1000 };

    engine->flow[0].chamber_nozzle_flow_area_m2 = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
    engine->flow[1].chamber_nozzle_flow_area_m2 = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
    engine->flow[2].chamber_nozzle_flow_area_m2 = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
    engine->flow[3].chamber_nozzle_flow_area_m2 = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
    engine->flow[4].chamber_nozzle_flow_area_m2 = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
    engine->flow[5].chamber_nozzle_flow_area_m2 = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
    engine->flow[6].chamber_nozzle_flow_area_m2 = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
    engine->flow[7].chamber_nozzle_flow_area_m2 = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
    engine->reset_engine();
    return engine;
}
