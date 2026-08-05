#pragma once

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
            const real Pg = chamber_static_pressure_pa[i];
            const real A = g_pi_r * diameter_m[i] * diameter_m[i];
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
