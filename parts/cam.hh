#pragma once

template<size_t W>
struct simple_cam
{
    lane<W> engage_theta_r = {};
    lane<W> ramp_theta_r = {};
    lane<W> open_ratio = {};

    real crankshaft_theta_r = 0.0_r;

    /*      4            1      2      3
     * r = t  [ 35 - 84 t + 70 t - 20 t ]
     *
     */

    fn void calc_open_ratios()
    {
        for(size_t i = 0; i < W; i++)
        {
            real theta0_r = modulos(engage_theta_r[i], g_otto_cycle_r);
            if(theta0_r < 0.0_r)
            {
                theta0_r += g_otto_cycle_r;
            }
            real theta1_r = modulos(crankshaft_theta_r, g_otto_cycle_r);
            if(theta1_r < theta0_r)
            {
                theta1_r += g_otto_cycle_r;
            }
            const real open_r = theta1_r - theta0_r;
            const real t = open_r / ramp_theta_r[i];
            const real a = t * t * t * t;
            const real b = t * a;
            const real c = t * b;
            const real d = t * c;
            const real A = 35.0_r * a;
            const real B = 84.0_r * b;
            const real C = 70.0_r * c;
            const real D = 20.0_r * d;
            const real R = clamper(A - B + C - D, 0.0_r, 1.0_r);
            open_ratio[i] = theta1_r < theta0_r ? 0.0_r : R;
        }
    }

    void update()
    {
        calc_open_ratios();
    }
};
