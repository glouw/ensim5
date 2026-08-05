#pragma once

template<size_t W>
struct sparkplugs
{
    static constexpr real fire_delay_theta_r = 1e-1_r;
    lane<W> engage_theta_r = {};
    mask<W> prev_fired = {};
    mask<W> fired = {};
    mask<W> rising_edge = {};
    real crankshaft_theta_r = 0.0_r;

    fn void calc_fired()
    {
        for(size_t i = 0; i < W; i++)
        {
            prev_fired[i] = fired[i];
            const real theta0_r = modulos(crankshaft_theta_r, g_otto_cycle_r);
            const real theta1_r = engage_theta_r[i] >= g_otto_cycle_r ? (engage_theta_r[i] - g_otto_cycle_r) : engage_theta_r[i];
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
