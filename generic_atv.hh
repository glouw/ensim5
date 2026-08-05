#pragma once

struct generic_atv : as_engine<4, 9, 2, 4, 5, inline_pistons, simple_cam, sparkplugs>
{
    generic_atv()
    {
        this->dc.set_cutoff_frequency(1000.0_r);
        this->pregain.ratio = 0.1_r;
        this->gain.ratio = 0.00001_r;
        this->limiter.max_angular_velocity_r_per_s = 900.0_r;
        this->limiter.limit_time_s = 0.1_r;
        this->flywheel.mass_kg = 3.0_r;
        this->flywheel.radius_m = 0.60_r;
        this->crankshaft.mass_kg = 7.3_r;
        this->crankshaft.radius_m = 0.200_r;
        this->crankshaft.angular_velocity_r_per_s = 250.0_r;
        this->pistons.diameter_m.fill(0.089_r);
        this->pistons.crank_throw_length_m.fill(0.035_r);
        this->pistons.connecting_rod_length_m.fill(0.120_r);
        this->pistons.connecting_rod_mass_kg.fill(0.55_r);
        this->pistons.head_mass_density_kg_per_m3.fill(2700.0_r);
        this->pistons.head_compression_height_m.fill(0.03_r);
        this->pistons.head_clearance_height_m.fill(0.006_r);
        this->pistons.friction_n_m_s2_per_r2.fill(0.00005_r);
        this->inlet_cam.ramp_theta_r.fill(g_pi_r * 0.8_r);
        this->outlet_cam.ramp_theta_r.fill(g_pi_r * 0.8_r);
        real theta0_r = 0.0_r;
        for(size_t i = 0; i < get_width(); i++)
        {
            this->pistons.theta0_r[i] = theta0_r;
            this->inlet_cam.engage_theta_r[i]  = theta0_r + g_otto_intake_cycle_r - g_pi_r / 8.0_r;
            this->sparkplugs.engage_theta_r[i] = theta0_r + g_otto_combustion_cycle_r - g_pi_r / 8.0_r;
            this->outlet_cam.engage_theta_r[i] = theta0_r + g_otto_exhaust_cycle_r - g_pi_r / 8.0_r;
            theta0_r += g_otto_cycle_r / static_cast<real>(get_width());
        }
        for(auto& flow : this->flows)
        {
            flow.chamber_nozzle_open_ratio.fill(1.0_r);
            flow.chamber_volume_m3 = {
                g_resevoir_volume_m3,
                0.00200_r, /* Intake   */
                0.00050_r, /* Throttle */
                0.00025_r, /* Runner   */
                0.00000_r, /* Piston   */
                0.00065_r, /* Runner   */
                0.00010_r, /* Exhaust1 */
                0.00013_r, /* Exhaust2 */
                g_resevoir_volume_m3,
            };
            flow.chamber_nozzle_flow_area_m2 = {
                0.00100_r, /* Atmosphere to Intake   */
                0.00080_r, /* Intake to Throttle     */
                0.00075_r, /* Throttle to Runner     */
                0.00055_r, /* Runner to Piston       */
                0.00040_r, /* Piston to Runner       */
                0.00120_r, /* Runner to Exhaust1     */
                0.00170_r, /* Exhaust1 to Exhaust2   */
                0.00250_r, /* Exhaust2 to Atmosphere */
            };
        }
        this->throttle.table = {
            0.00000_r,
            0.00015_r,
            0.00500_r,
            0.20000_r,
        };
    }
};
