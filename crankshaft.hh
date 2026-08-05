#pragma once

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
