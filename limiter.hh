#pragma once

struct limiter
{
    real max_angular_velocity_r_per_s = 600.0_r;
    real crankshaft_angular_velocity_r_per_s = 0.0_r;
    real limit_time_s = 0.1;
    real cycles = 0;
    bool limiting = false;

    void update()
    {
        if(not limiting)
        {
            if(crankshaft_angular_velocity_r_per_s > max_angular_velocity_r_per_s)
            {
                limiting = true;
            }
        }
        if(limiting)
        {
            const real time_s = g_dt_s * cycles;
            if(time_s > limit_time_s)
            {
                limiting = false;
                cycles = 0;
            }
            cycles++;
        }
    }
};
