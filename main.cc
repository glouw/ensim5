#include <cmath>

namespace ensim5
{
using std::sin;
using std::cos;
using std::sqrt;
using std::fmax;
using std::fmin;

template<typename T> requires(std::same_as<T, float> || std::same_as<T, double>)
class dbl_c
{
    T _value;

public:
    constexpr dbl_c() = default;
    constexpr dbl_c(const double v): _value(static_cast<T>(v)) {}

    constexpr dbl_c& operator+=(const dbl_c other)
    {
        _value += other._value;
        return *this;
    }

    constexpr operator T() const
    {
        return _value;
    }
};

using val_t = dbl_c<float>;

constexpr size_t g_sample_rate_hz = 48000;
constexpr val_t g_dt_s = val_t(1.0) / g_sample_rate_hz;

template<size_t W> requires(W % 2 == 1)
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

    static constexpr val_t _nozzle_gamma = 1.333;
    static constexpr val_t _chamber_cv_n2_j_per_mol_k = 23.0;
    static constexpr val_t _chamber_cv_o2_j_per_mol_k = 22.5;
    static constexpr val_t _chamber_cv_ar_j_per_mol_k = 12.5;
    static constexpr val_t _chamber_cv_c8h18_j_per_mol_k = 250.0;
    static constexpr val_t _chamber_cv_co2_j_per_mol_k = 33.0;
    static constexpr val_t _chamber_cv_h2o_j_per_mol_k = 30.0;
    static constexpr val_t _chamber_molar_mass_kg_per_mol_c8h18 = 0.1142285200;
    static constexpr val_t _chamber_molar_mass_kg_per_mol_o2 = 0.0319988000;
    static constexpr val_t _chamber_molar_mass_kg_per_mol_n2 = 0.0280134000;
    static constexpr val_t _chamber_molar_mass_kg_per_mol_ar = 0.0399480000;
    static constexpr val_t _chamber_molar_mass_kg_per_mol_co2 = 0.0440095000;
    static constexpr val_t _chamber_molar_mass_kg_per_mol_h2o = 0.0180152800;
    static constexpr val_t _chamber_universal_gas_constant_j_per_mol_k = 8.3144598;

    using lane_t = val_t[W];
    using mask_t = bool[W];

    lane_t _chamber_nozzle_flow_area_m2;
    lane_t _chamber_volume_m3;
    lane_t _chamber_molar_mass_kg_per_mol;
    lane_t _chamber_cv_j_per_mol_k;
    lane_t _chamber_specific_gas_constant_j_per_kg_k;
    lane_t _chamber_static_pressure_pa;
    lane_t _chamber_moles;
    lane_t _chamber_total_cv_j_per_k;
    lane_t _chamber_mol_ratio_n2;
    lane_t _chamber_mol_ratio_o2;
    lane_t _chamber_mol_ratio_ar;
    lane_t _chamber_mol_ratio_c8h18;
    lane_t _chamber_mol_ratio_co2;
    lane_t _chamber_mol_ratio_h2o;
    lane_t _chamber_static_temperature_k;
    lane_t _chamber_mass_kg;

    lane_t _nozzle_mach;
    lane_t _nozzle_velocity_m_per_s;
    lane_t _nozzle_static_density_kg_per_m3;
    lane_t _nozzle_mass_flow_rate_kg_per_s;

    /* Parcels are tiny chamber masses determined by nozzle mass flow rate. */
    lane_t _parcel_molar_mass_kg_per_mol;
    lane_t _parcel_cv_j_per_mol_k;
    lane_t _parcel_moles;
    lane_t _parcel_total_cv_j_per_k;
    lane_t _parcel_mol_ratio_n2;
    lane_t _parcel_mol_ratio_o2;
    lane_t _parcel_mol_ratio_ar;
    lane_t _parcel_mol_ratio_c8h18;
    lane_t _parcel_mol_ratio_co2;
    lane_t _parcel_mol_ratio_h2o;
    lane_t _parcel_static_temperature_k;
    lane_t _parcel_mass_kg;
    mask_t _parcel_flow_right;

    val_t clamp(const val_t val, const val_t lower, const val_t upper)
    {
        return fmax(fmin(val, val_t(upper)), val_t(lower));
    }

    val_t calc_mix(const val_t v1, const val_t w1, const val_t v2, const val_t w2)
    {
        return (v1 * w1 + v2 * w2) / (w1 + w2);
    }

    /*
     *      ____
     *      \
     * cv = /___ mr[i] * M
     *         i
     */

    inline void calc_chamber_molar_masses()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            _chamber_molar_mass_kg_per_mol [i]
                = _chamber_mol_ratio_n2    [i] * _chamber_molar_mass_kg_per_mol_n2
                + _chamber_mol_ratio_o2    [i] * _chamber_molar_mass_kg_per_mol_o2
                + _chamber_mol_ratio_ar    [i] * _chamber_molar_mass_kg_per_mol_ar
                + _chamber_mol_ratio_c8h18 [i] * _chamber_molar_mass_kg_per_mol_c8h18
                + _chamber_mol_ratio_co2   [i] * _chamber_molar_mass_kg_per_mol_co2
                + _chamber_mol_ratio_h2o   [i] * _chamber_molar_mass_kg_per_mol_h2o;
        }
    }

    /*
     *      ____
     *      \
     * cv = /___ mr[i] * cv
     *         i
     */

    inline void calc_chamber_molar_cvs()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            _chamber_cv_j_per_mol_k[i]
                = _chamber_mol_ratio_n2    [i] * _chamber_cv_n2_j_per_mol_k
                + _chamber_mol_ratio_o2    [i] * _chamber_cv_o2_j_per_mol_k
                + _chamber_mol_ratio_ar    [i] * _chamber_cv_ar_j_per_mol_k
                + _chamber_mol_ratio_c8h18 [i] * _chamber_cv_c8h18_j_per_mol_k
                + _chamber_mol_ratio_co2   [i] * _chamber_cv_co2_j_per_mol_k
                + _chamber_mol_ratio_h2o   [i] * _chamber_cv_h2o_j_per_mol_k;
        }
    }

    /*
     *       R
     * Rs = ---
     *       M
     */

    inline void calc_chamber_specific_gas_constants()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            const val_t R = _chamber_universal_gas_constant_j_per_mol_k;
            const val_t M = _chamber_molar_mass_kg_per_mol[i];
            _chamber_specific_gas_constant_j_per_kg_k[i] = R / M;
        }
    }

    /*
     *      m * Rs * Ts
     * Ps = -----------
     *           V
     */

    inline void calc_chamber_static_pressures()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            const val_t m = _chamber_mass_kg[i];
            const val_t Rs = _chamber_specific_gas_constant_j_per_kg_k[i];
            const val_t Ts = _chamber_static_temperature_k[i];
            const val_t V = _chamber_volume_m3[i];
            _chamber_static_pressure_pa[i] = m * Rs * Ts / V;
        }
    }

    /*
     *      m
     * n = ---
     *      M
     */

    inline void calc_chamber_moles()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            const val_t m = _chamber_mass_kg[i];
            const val_t M = _chamber_molar_mass_kg_per_mol[i];
            _chamber_moles[i] = m / M;
        }
    }

    /*
     * cv(total) = n * cv
     */

    inline void calc_chamber_total_cvs()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            _chamber_total_cv_j_per_k[i] = _chamber_moles[i] * _chamber_cv_j_per_mol_k[i];
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

    inline void calc_nozzle_machs()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            const size_t j = i + 1;
            const val_t X = (val_t(2.0) / (_nozzle_gamma - val_t(1.0)));
            const val_t Px = _chamber_static_pressure_pa[i];
            const val_t Py = _chamber_static_pressure_pa[j];
            const val_t Pt = fmax(Px, Py);
            const val_t Ps = fmin(Px, Py);
            const val_t direction = Px > Py ? val_t(1.0) : val_t(-1.0);
            _parcel_flow_right[i] = direction > val_t(0.0);

            /*
             *      y - 1    1.333 - 1                        1     1          __________
             *      ----- = ---------- = 0.25                --- * ---        /    ______
             *        y       1.333                0.25       2     2        /    /
             *  term                         = term    = term          = _  / _  /  term
             *                                                            \/   \/
             */

            static_assert(_nozzle_gamma == val_t(1.333), "_nozzle_gamma (y) must be 1.333 to use manual pow optimization");
            const val_t Y = sqrt(sqrt(Pt / Ps));
            const val_t M = direction * sqrt(X * (Y - val_t(1.0)));
            _nozzle_mach[i] = clamp(M, val_t(-1.0), val_t(1.0));
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

    inline void calc_nozzle_velocities()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            const val_t Rs = _chamber_specific_gas_constant_j_per_kg_k[i];
            const val_t Tt = _chamber_static_temperature_k[i];
            const val_t M = _nozzle_mach[i];
            const val_t X = _nozzle_gamma * Rs * Tt;
            const val_t Y = val_t(0.5) * (_nozzle_gamma - val_t(1.0)) * M * M;
            const val_t u = M * sqrt(X / (val_t(1.0) + Y));
            _nozzle_velocity_m_per_s[i] = u;
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

    inline void calc_nozzle_static_densities()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            const val_t Pt = _chamber_static_pressure_pa[i];
            const val_t Rs = _chamber_specific_gas_constant_j_per_kg_k[i];
            const val_t Tt = _chamber_static_temperature_k[i];
            const val_t M = _nozzle_mach[i];
            const val_t X = Pt / (Rs * Tt);
            const val_t C = val_t(1.0) + val_t(0.5) * (_nozzle_gamma - val_t(1.0)) * M * M;

            /*
             *         1         1
             *       ----- = --------- = 3.00
             *       y - 1   1.333 - 1
             *  term                         = term * term * term
             */

            static_assert(_nozzle_gamma == val_t(1.333), "_nozzle_gamma (y) must be 1.333 to use manual pow optimization");
            _nozzle_static_density_kg_per_m3[i] = X / (C * C * C);
        }
    }

    inline void calc_parcel_property(lane_t& __restrict nozzle, const lane_t& __restrict chamber)
    {
        for(size_t i = 0; i < W - 1; i++) nozzle[i] = val_t(0.0);
        for(size_t i = 0; i < W - 1; i++) nozzle[i] += _parcel_flow_right[i] ? chamber[i + 0] : val_t(0.0);
        for(size_t i = 0; i < W - 1; i++) nozzle[i] += _parcel_flow_right[i] ? val_t(0.0) : chamber[i + 1];
    }

    inline void calc_parcel_properties()
    {
        calc_parcel_property(_parcel_molar_mass_kg_per_mol, _chamber_molar_mass_kg_per_mol);
        calc_parcel_property(_parcel_cv_j_per_mol_k,        _chamber_cv_j_per_mol_k);
        calc_parcel_property(_parcel_mol_ratio_n2,          _chamber_mol_ratio_n2);
        calc_parcel_property(_parcel_mol_ratio_o2,          _chamber_mol_ratio_o2);
        calc_parcel_property(_parcel_mol_ratio_ar,          _chamber_mol_ratio_ar);
        calc_parcel_property(_parcel_mol_ratio_c8h18,       _chamber_mol_ratio_c8h18);
        calc_parcel_property(_parcel_mol_ratio_co2,         _chamber_mol_ratio_co2);
        calc_parcel_property(_parcel_mol_ratio_h2o,         _chamber_mol_ratio_h2o);
        calc_parcel_property(_parcel_static_temperature_k,  _chamber_static_temperature_k);
    }

    /* .
     * m = ps A u
     *
     */

    inline void calc_nozzle_mass_flow_rates()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            const val_t ps = _nozzle_static_density_kg_per_m3[i];
            const val_t A = _chamber_nozzle_flow_area_m2[i];
            const val_t u = _nozzle_velocity_m_per_s[i];
            const val_t mdot = ps * A * u;
            _nozzle_mass_flow_rate_kg_per_s[i] = mdot;
            const val_t dm = mdot * g_dt_s;
            _parcel_mass_kg[i] = dm;
            _parcel_moles[i] = dm / _parcel_molar_mass_kg_per_mol[i];
        }
    }

    /*
     *        nr1 n1 + nr2 n2
     * nr1 = --------------------
     *            n1 + n2
     */

    inline void calc_chamber_mol_ratios()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            _chamber_mol_ratio_n2   [i] = calc_mix(_chamber_mol_ratio_n2   [i], _chamber_moles[i], _parcel_mol_ratio_n2   [i], _parcel_moles[i]);
            _chamber_mol_ratio_o2   [i] = calc_mix(_chamber_mol_ratio_o2   [i], _chamber_moles[i], _parcel_mol_ratio_o2   [i], _parcel_moles[i]);
            _chamber_mol_ratio_ar   [i] = calc_mix(_chamber_mol_ratio_ar   [i], _chamber_moles[i], _parcel_mol_ratio_ar   [i], _parcel_moles[i]);
            _chamber_mol_ratio_c8h18[i] = calc_mix(_chamber_mol_ratio_c8h18[i], _chamber_moles[i], _parcel_mol_ratio_c8h18[i], _parcel_moles[i]);
            _chamber_mol_ratio_co2  [i] = calc_mix(_chamber_mol_ratio_co2  [i], _chamber_moles[i], _parcel_mol_ratio_co2  [i], _parcel_moles[i]);
            _chamber_mol_ratio_h2o  [i] = calc_mix(_chamber_mol_ratio_h2o  [i], _chamber_moles[i], _parcel_mol_ratio_h2o  [i], _parcel_moles[i]);
        }
    }

    /*
     * cv(total) = n * cv
     */

    inline void calc_parcel_total_cvs()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            _parcel_total_cv_j_per_k[i] = _parcel_moles[i] * _parcel_cv_j_per_mol_k[i];
        }
    }

    /*
     *        n1 cv1 Ts1 + n2 cv2 Ts2
     * Ts1 = -------------------------
     *            n1 cv1 + n2 cv2
     */

    inline void calc_chamber_static_temperatures()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            _chamber_static_temperature_k[i] = calc_mix(
                _chamber_static_temperature_k[i],
                _chamber_total_cv_j_per_k[i],
                _parcel_static_temperature_k[i],
                _parcel_total_cv_j_per_k[i]
            );
        }
    }

    inline void calc_chamber_masses()
    {
        for(size_t i = 0; i < W - 1; i++)
        {
            _chamber_mass_kg[i] += _parcel_mass_kg[i];
        }
    }

    void calc_flow()
    {
        calc_chamber_molar_masses();
        calc_chamber_molar_cvs();
        calc_chamber_specific_gas_constants();
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
    val_t _theta_r;
    val_t _angular_velocity_r_per_s;
    val_t _mass_kg;
    val_t _radius_m;
    val_t _moment_of_inertia_kg_m2;

    /*
     * dw = dw/dt * dt
     *
     */

    inline void accelerate(const val_t angular_acceleration_r_per_s2)
    {
        _angular_velocity_r_per_s += angular_acceleration_r_per_s2 * g_dt_s;
    }

    /*
     * dth = dth/dt * dt
     *
     */

    inline void turn()
    {
        _theta_r += _angular_velocity_r_per_s * g_dt_s;
    }

    /*
     *      1         2
     * I = --- * m * r
     *      2
     *
     */

    inline void calc_moment_of_inertia()
    {
        _moment_of_inertia_kg_m2 = val_t(0.5) * _mass_kg * _radius_m * _radius_m;
    }

    void calc_crankshaft()
    {
        calc_moment_of_inertia();
    }
};

template<size_t H, size_t PistonX>
struct pistons_s
{
    using lane_t = val_t[H];

    lane_t _diameter_m;
    lane_t _pin_x_m;
    lane_t _pin_y_m;
    lane_t _bearing_x_m;
    lane_t _bearing_y_m;
    lane_t _theta0_r;
    lane_t _crank_throw_length_m;
    lane_t _connecting_rod_length_m;
    lane_t _connecting_rod_mass_kg;
    lane_t _head_mass_density_kg_per_m3;
    lane_t _head_compression_height_m;
    lane_t _head_clearance_height_m;
    lane_t _head_mass_kg;
    lane_t _theta_r;
    lane_t _gas_torque_n_m;
    lane_t _moment_of_inertia_kg_m2;
    lane_t _inertia_torque_n_m;
    lane_t _sin;
    lane_t _cos;

    /*
     * t = t0 + t1
     */

    inline void calc_thetas(const crankshaft_s& __restrict crankshaft)
    {
        for(size_t i = 0; i < H; i++)
        {
            _theta_r[i] = _theta0_r[i] + crankshaft._theta_r;
        }
    }

    inline void calc_sin_cos()
    {
        for(size_t i = 0; i < H; i++)
        {
            _sin[i] = sin(_theta_r[i]);
            _cos[i] = cos(_theta_r[i]);
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

    inline void calc_positions()
    {
        for(size_t i = 0; i < H; i++)
        {
            const val_t r = _crank_throw_length_m[i];
            const val_t l = _connecting_rod_length_m[i];
            const val_t t = _theta_r[i];
            const val_t x = r * _sin[i];
            const val_t y = r * _cos[i];
            _bearing_x_m[i] = x;
            _bearing_y_m[i] = y;
            _pin_x_m[i] = val_t(0.0);
            _pin_y_m[i] = y + sqrt(l * l - x * x);
        }
    }

    /*           2
     * M = pi * r  * h * p
     *
     */

    inline void calc_masses()
    {
        for(size_t i = 0; i < H; i++)
        {
            const val_t r = val_t(0.5) * _diameter_m[i];
            const val_t h = val_t(2.0) * _head_compression_height_m[i];
            const val_t p = _head_mass_density_kg_per_m3[i];
            _head_mass_kg[i] = val_t(M_PI) * r * r * h * p;
        }
    }

    /*             1        2
     * I = [ mp + --- mr ] r
     *             3
     */

    inline void calc_moments_of_inertia()
    {
        for(size_t i = 0; i < H; i++)
        {
            const val_t r = _crank_throw_length_m[i];
            const val_t mp = _head_mass_kg[i];
            const val_t mr = _connecting_rod_mass_kg[i];
            _moment_of_inertia_kg_m2[i] = (mp + (val_t(1.0 / 3.0) * mr)) * r * r;
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

    template<size_t W>
    inline void calc_gas_torques(const flow_s<W>& __restrict flow)
    {
        for(size_t i = 0; i < H; i++)
        {
            const val_t Pg = flow._chamber_static_pressure_pa[PistonX];
            const val_t A = val_t(M_PI) * _diameter_m[i] * _diameter_m[i];
            const val_t r = _crank_throw_length_m[i];
            const val_t l = _connecting_rod_length_m[i];
            const val_t t = _theta_r[i];
            const val_t X = Pg * A * r * _sin[i];
            const val_t Y = val_t(1.0) + (r / l) * _cos[i];
            _gas_torque_n_m[i] = X * Y;
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

    inline void calc_inertia_torques(const crankshaft_s& __restrict crankshaft)
    {
        for(size_t i = 0; i < H; i++)
        {
            const val_t r = _crank_throw_length_m[i];
            const val_t l = _connecting_rod_length_m[i];
            const val_t I = _moment_of_inertia_kg_m2[i];
            const val_t w = crankshaft._angular_velocity_r_per_s;
            const val_t t = _theta_r[i];
            const val_t rl = r / l;
            const val_t s = _sin[i];
            const val_t c = _cos[i];
            const val_t X = val_t(0.25) * rl * s;
            const val_t Y = s * c;
            const val_t Z = val_t(0.75) * rl * (val_t(3.0) * s - val_t(4.0) * s * s * s);
            _inertia_torque_n_m[i] = I * w * w * (X - Y - Z);
        }
    }

    template<size_t W>
    void calc_pistons(const flow_s<W>& __restrict flow, const crankshaft_s& __restrict crankshaft)
    {
        calc_thetas(crankshaft);
        calc_sin_cos();
        calc_positions();
        calc_masses();
        calc_moments_of_inertia();
        calc_gas_torques(flow);
        calc_inertia_torques(crankshaft);
    }
};

template struct flow_s<9>;
template struct pistons_s<8, 5>;

void run()
{
    crankshaft_s crankshaft;
    flow_s<9> flow;
    pistons_s<8, 5> pistons;
    crankshaft.calc_crankshaft();
    flow.calc_flow();
    pistons.calc_pistons(flow, crankshaft);
}
}
