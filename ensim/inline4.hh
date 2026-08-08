#pragma once

struct inline4 : as_engine<
    /* W          */ 4,
    /* H          */ 9,
    /* THROTTLE_Y */ 2,
    /* PISTON_Y   */ 4,
    /* AUDIO_Y    */ 7,
    /* PIPE_CELLS */ 192,
    inline_pistons,
    basic_cams,
    basic_sparkplugs,
    lfr_pipe>
{
    inline4()
    {
        this->limiter.max_angular_velocity_r_per_s = 1500.0_r;
        this->limiter.limit_time_s = 0.040_r;
        this->flywheel.mass_kg = 9.0_r;
        this->flywheel.radius_m = 0.12_r;
        this->crankshaft.mass_kg = 12.5_r;
        this->crankshaft.radius_m = 0.045_r;
        this->crankshaft.angular_velocity_r_per_s = 500.0_r;
        this->pistons.diameter_m.fill(0.086_r);
        this->pistons.crank_throw_length_m.fill(0.043_r);
        this->pistons.connecting_rod_length_m.fill(0.145_r);
        this->pistons.connecting_rod_mass_kg.fill(0.45_r);
        this->pistons.head_mass_density_kg_per_m3.fill(2700.0_r);
        this->pistons.head_compression_height_m.fill(0.030_r);
        this->pistons.head_clearance_height_m.fill(0.007_r);
        this->pistons.friction_n_m_s2_per_r2.fill(0.000035_r);
        this->inlet_cam.ramp_theta_r.fill(g_pi_r * 0.8_r);
        this->outlet_cam.ramp_theta_r.fill(g_pi_r * 0.53_r);
        real theta0_r = 0.0_r;
        for(size_t i = 0; i < get_width(); i++)
        {
            this->pistons.theta0_r[i] = theta0_r;
            this->inlet_cam.engage_theta_r[i]  = theta0_r + g_otto_intake_cycle_r - g_pi_r / 8.0_r;
            this->sparkplugs.engage_theta_r[i] = theta0_r + g_otto_combustion_cycle_r - g_pi_r / 6.2_r;
            this->outlet_cam.engage_theta_r[i] = theta0_r + g_otto_exhaust_cycle_r + 0.8;
            theta0_r += g_otto_cycle_r / static_cast<real>(get_width());
        }
        for(auto& flow : this->flows)
        {
            flow.chamber_nozzle_open_ratio.fill(1.0_r);
            flow.chamber_nozzle_flow_area_m2 = {
                0.00250_r, /* ATMOSPHERIC SOURCE -> INTAKE */
                0.00120_r, /* INTAKE -> THROTTLE           */
                0.00030_r, /* THROTTLE -> RUNNER           */
                0.00060_r, /* RUNNER -> PISTON             */
                0.00060_r, /* PISTON -> RUNNER             */
                0.00110_r, /* RUNNER -> CHAMBER1           */
                0.00180_r, /* CHAMBER1 -> CHAMBER2         */
                0.00300_r, /* CHAMBER2 -> ATMOSPHERIC SINK */
            };
        }
        /*                             ATMOSPHERIC SOURCE    INTAKE   THROTTLE  RUNNER    PISTON     RUNNER    CHAMBER1  CHAMBER2  ATMOSPHERIC SINK     */
        flows[0].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.00020_r, 0.0006_r, 0.0006_r, g_resevoir_volume_m3 };
        flows[1].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.00048_r, 0.0006_r, 0.0006_r, g_resevoir_volume_m3 };
        flows[2].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.00043_r, 0.0006_r, 0.0006_r, g_resevoir_volume_m3 };
        flows[3].chamber_volume_m3 = { g_resevoir_volume_m3, 0.003_r, 0.0008_r, 0.0003_r, 0.00000_r, 0.00020_r, 0.0006_r, 0.0006_r, g_resevoir_volume_m3 };
        this->throttle.table = {
            0.00050_r,
            0.01000_r,
            0.10000_r,
            1.00000_r,
        };
        /* ALL CHAMBER2s CONNECT TO 1D CFD PIPE */
        const real pipe_length_m = 1.5_r;
        const real pipe_radius_m = 0.005_r;
        const real pipe_mic_position_ratio = 0.66_r;
        const real pipe_reflection_ratio = -0.33_r;
        const size_t pipe_solver_substeps = 12;
        this->pipe.configure(pipe_length_m, pipe_radius_m, pipe_mic_position_ratio, pipe_reflection_ratio, pipe_solver_substeps);
        this->dc.set_cutoff_frequency(120.0_r);
        this->gain.ratio = 0.0000002_r;
    }
};
