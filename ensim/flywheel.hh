#pragma once

struct flywheel
{
    real mass_kg = 0.0_r;
    real radius_m = 0.0_r;
    real moment_of_inertia_kg_m2 = 0.0_r;

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

    void update()
    {
        calc_moment_of_inertia();
    }
};
