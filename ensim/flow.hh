#pragma once

template<size_t H, size_t PY>
struct flow
{
    static constexpr size_t N = H - 1;
    static_assert(N % 2 == 0);

    lane<H> chamber_prev_volume_m3 = {};
    lane<H> chamber_volume_m3 = {};
    lane<H> chamber_nozzle_flow_area_m2 = {};
    lane<H> chamber_nozzle_real_flow_area_m2 = {};
    lane<H> chamber_nozzle_open_ratio = {};
    lane<H> chamber_static_pressure_pa = {};
    lane<H> chamber_dynamic_pressure_pa = {};
    lane<H> chamber_total_pressure_pa = {};
    lane<H> chamber_static_temperature_k = {};
    lane<H> chamber_dynamic_temperature_k = {};
    lane<H> chamber_total_temperature_k = {};
    lane<H> chamber_mass_kg = {};
    lane<H> chamber_bulk_momentum_kg_m_per_s = {};
    lane<H> nozzle_mach = {};
    lane<H> nozzle_velocity_m_per_s = {};
    lane<H> nozzle_static_density_kg_per_m3 = {};
    lane<H> nozzle_mass_flow_rate_kg_per_s = {};
    lane<H> parcel_mass_kg = {};
    lane<H> parcel_static_temperature_k = {};
    mask<H> panic = {};
    real piston_injection_enabled = 0.0_r;
    real piston_chamber_flame_height_m = 0.0_r;
    real piston_chamber_mass_burned_m3 = 0.0_r;
    real piston_chamber_radius_m = 0.0_r;
    bool piston_chamber_on_fire = false;

    /*
     *     Ps * V
     * m = -------
     *     Rs * Ts
     */

    fn void calc_chamber_ambients()
    {
        for(size_t i = 0; i < H; i++)
        {
            chamber_static_temperature_k[i] = g_ambient_temperature_k;
            chamber_static_pressure_pa[i] = g_ambient_pressure_pa;
        }
        for(size_t i = 0; i < H; i++)
        {
            const real Ps = chamber_static_pressure_pa[i];
            const real V = chamber_volume_m3[i];
            const real Rs = g_specific_gas_constant_j_per_kg_k;
            const real Ts = chamber_static_temperature_k[i];
            chamber_mass_kg[i] = (Ps * V) / (Rs * Ts);
        }
    }

    /*
     *       m * Rs * T
     * Ps = ------------
     *           V
     */

    fn void calc_chamber_static_pressures()
    {
        for(size_t i = 0; i < N; i++)
        {
            const real m = chamber_mass_kg[i];
            const real Rs = g_specific_gas_constant_j_per_kg_k;
            const real Ts = chamber_static_temperature_k[i];
            const real V = chamber_volume_m3[i];
            chamber_static_pressure_pa[i] = m * Rs * Ts / V;
        }
    }

    /*              ____________________
     *             /
     *            /            y - 1
     *           /             -----
     *          /                y
     *         /   2         Pt
     * M = _  /  ----- * [ (----) - 1 ]
     *      \/   y - 1       Ps
     */

    fn void calc_nozzle_machs()
    {
        for(size_t i = 0; i < N; i++)
        {
            const size_t j = i + 1;
            const real X = (2.0_r / (g_gamma - 1.0_r));
            const real Pi = chamber_total_pressure_pa[i];
            const real Pj = chamber_total_pressure_pa[j];
            const real Pt = fmax(Pi, Pj);
            const real Ps = fmin(Pi, Pj);
            const real direction = Pi > Pj ? 1.0_r : -1.0_r;

            /*
             *      y - 1                        3
             *      ----- = 0.3333... where y = ---
             *        y                          2
             * Term
             */

            static_assert(g_gamma == 3.0_r / 2.0_r);
            const real Y = cuberoot(Pt / Ps);
            const real M = direction * sqrt(X * (Y - 1.0_r));
            nozzle_mach[i] = M;
        }
    }

    /*
     *
     * Pt = Ps + Pd
     *
     */

    fn void calc_chamber_total_pressures()
    {
        for(size_t i = 0; i < H; i++)
        {
            const real Ps = chamber_static_pressure_pa[i];
            const real Pd = chamber_dynamic_pressure_pa[i];
            chamber_total_pressure_pa[i] = Ps + Pd;
        }
    }

    /*
     *
     * Tt = Ts + Td
     *
     */

    fn void calc_chamber_total_temperatures()
    {
        for(size_t i = 0; i < H; i++)
        {
            const real Ts = chamber_static_temperature_k[i];
            const real Td = chamber_dynamic_temperature_k[i];
            chamber_total_temperature_k[i] = Ts + Td;
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
        for(size_t i = 0; i < N; i++)
        {
            const real Rs = g_specific_gas_constant_j_per_kg_k;
            const real Tt = chamber_total_temperature_k[i];
            const real M = nozzle_mach[i];
            const real X = g_gamma * Rs * Tt;
            const real Y = 0.5_r * (g_gamma - 1.0_r) * M * M;
            const real u = M * sqrt(X / (1.0_r + Y));
            const real A = chamber_nozzle_real_flow_area_m2[i];
            const real mute = A == 0.0_r ? 0.0_r : 1.0_r;
            nozzle_velocity_m_per_s[i] = u * mute;
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
        for(size_t i = 0; i < N; i++)
        {
            const real Pt = chamber_total_pressure_pa[i];
            const real Rs = g_specific_gas_constant_j_per_kg_k;
            const real Tt = chamber_total_temperature_k[i];
            const real M = nozzle_mach[i];
            const real pt = Pt / (Rs * Tt);

            /*
             *        1                  3
             *      ----- = 2 where y = ---
             *      y - 1                2
             * Term
             */

            static_assert(g_gamma == 3.0_r / 2.0_r);
            const real C = 1.0_r + 0.5_r * (g_gamma - 1.0_r) * M * M;
            nozzle_static_density_kg_per_m3[i] = pt / (C * C);
        }
    }

    fn void calc_nozzle_real_flow_areas()
    {
        for(size_t i = 0; i < N; i++)
        {
            const real r = chamber_nozzle_open_ratio[i];
            const real A = chamber_nozzle_flow_area_m2[i];
            chamber_nozzle_real_flow_area_m2[i] = r * A;
        }
    }

    /* .
     * m = ps A u
     *
     */

    fn void calc_nozzle_mass_flow_rates()
    {
        for(size_t i = 0; i < N; i++)
        {
            const real ps = nozzle_static_density_kg_per_m3[i];
            const real A = chamber_nozzle_real_flow_area_m2[i];
            const real u = nozzle_velocity_m_per_s[i];
            const real mdot = ps * A * u;
            nozzle_mass_flow_rate_kg_per_s[i] = mdot;
        }
    }

    /*
     *      .
     * mp = m dt
     * Tsp = Ts,upstream
     *
     */

    fn void calc_nozzle_parcels()
    {
        for(size_t i = 0; i < N; i++)
        {
            const size_t j = i + 1;
            const real mdot = nozzle_mass_flow_rate_kg_per_s[i];
            const real dm = mdot * g_dt_s;
            parcel_mass_kg[i] = dm;
            const real Tsi = chamber_static_temperature_k[i];
            const real Tsj = chamber_static_temperature_k[j];
            parcel_static_temperature_k[i] = mdot > 0.0_r ? Tsi : Tsj;
        }
    }

    /*                  y - 1
     *               V1
     * Ts2 = Ts1 * (----)
     *               V2
     */

    fn void calc_compressions()
    {
        for(size_t i = 0; i < N; i++)
        {
            const real Ts1 = chamber_static_temperature_k[i];
            const real V1 = chamber_prev_volume_m3[i];
            const real V2 = chamber_volume_m3[i];
            const real dv = V1 / V2;

            /*
             *               1             3
             *      y - 1 = --- where y = ---
             * Term          2             2
             *
             */

            static_assert(g_gamma == 3.0_r / 2.0_r);
            chamber_static_temperature_k[i] = Ts1 * sqrt(dv);
        }
    }

    /*
     *       Ts m + Tp dm
     * Ts = --------------
     *          m + dm
     */

    fn void calc_forward_energy_transfers()
    {
        for(size_t i = 0; i < N; i++)
        {
            const size_t j = i + 1;
            const real dm = parcel_mass_kg[i];
            const real m = chamber_mass_kg[j];
            const real Tsp = parcel_static_temperature_k[i];
            const real Ts0 = chamber_static_temperature_k[j];
            const real Ts1 = (Ts0 * m + Tsp * dm) / (m + dm);
            chamber_static_temperature_k[j] = dm > 0.0_r ? Ts1 : Ts0;
        }
    }

    fn void calc_reverse_energy_transfers()
    {
        for(size_t i = N; i > 0; i--)
        {
            const real dm = parcel_mass_kg[i];
            const real m = chamber_mass_kg[i];
            const real Tsp = parcel_static_temperature_k[i];
            const real Ts0 = chamber_static_temperature_k[i];
            const real Ts1 = (Ts0 * m - Tsp * dm) / (m - dm);
            chamber_static_temperature_k[i] = dm < 0.0_r ? Ts1 : Ts0;
        }
    }

    /*                  y - 1
     *               m1
     * Ts2 = Ts1 * (----)
     *               m2
     */

    fn void calc_mass_transfers()
    {
        for(size_t i = 0; i < N; i++)
        {
            const size_t j = i + 1;
            const real mi = chamber_mass_kg[i];
            const real mj = chamber_mass_kg[j];
            const real dm = parcel_mass_kg[i];
            const real m0 = mi - dm;
            const real m1 = mj + dm;

            /*
             *               1             3
             *      y - 1 = --- where y = ---
             * Term          2             2
             *
             */

            static_assert(g_gamma == 3.0_r / 2.0_r);
            chamber_static_temperature_k[i] *= sqrt(m0 / mi);
            chamber_static_temperature_k[j] *= sqrt(m1 / mj);
            chamber_mass_kg[i] = m0;
            chamber_mass_kg[j] = m1;
        }

        for(size_t i = 0; i < N; i++)
        {
            const size_t j = i + 1;
            const real dm = parcel_mass_kg[i];
            const real u = nozzle_velocity_m_per_s[i];
            const real p = dm * u;
            chamber_bulk_momentum_kg_m_per_s[i] -= p;
            chamber_bulk_momentum_kg_m_per_s[j] += p;
        }

        /*               ___________
         *              /
         * pmax = m _  / y * Rs * Ts
         *           \/
         */

        for(size_t i = 0; i < N; i++)
        {
            const real Rs = g_specific_gas_constant_j_per_kg_k;
            const real Ts = chamber_static_temperature_k[i];
            const real m = chamber_mass_kg[i];
            const real p = chamber_bulk_momentum_kg_m_per_s[i];
            const real pmax = m * sqrt(g_gamma * Rs * Ts);
            chamber_bulk_momentum_kg_m_per_s[i] = clamper(p, -pmax, pmax);
        }
    }

    /*
     *      1     2
     * q = --- p u
     *      2
     */

    fn void calc_chamber_dynamic_pressures()
    {
        for(size_t i = 0; i < N; i++)
        {
            const real u = chamber_bulk_momentum_kg_m_per_s[i] / chamber_mass_kg[i];
            const real p = chamber_mass_kg[i] / chamber_volume_m3[i];
            const real q = 0.5_r * p * u * u;
            chamber_dynamic_pressure_pa[i] = q;
        }
    }

    /*         2
     *        u
     * Td = ------
     *       2 Cp
     */

    fn void calc_chamber_dynamic_temperatures()
    {
        for(size_t i = 0; i < N; i++)
        {
            const real u = chamber_bulk_momentum_kg_m_per_s[i] / chamber_mass_kg[i];
            const real Cv = g_cv_j_per_kg_k;
            const real Td = 0.5_r * u * u / Cv;
            chamber_dynamic_temperature_k[i] = Td;
        }
    }

    fn void ignite_piston_chamber()
    {
        piston_chamber_flame_height_m = 0.0_r;
        piston_chamber_mass_burned_m3 = 0.0_r;
        piston_chamber_on_fire = true;
    }

    fn void calc_combustion()
    {
        if(piston_chamber_on_fire)
        {
            if(not piston_injection_enabled)
            {
                piston_chamber_on_fire = false;
                return;
            }
            const real M = chamber_mass_kg[PY];
            const real V = chamber_volume_m3[PY];

            /*
             *              Ts    2
             *           [------]
             *             Ts0
             * S =  0.4 ----------------
             *              Ps    0.125
             *           [------]
             *             Ps0
             */

            const real Ts = chamber_static_temperature_k[PY];
            const real Ts0 = g_ambient_temperature_k;
            const real Ps = chamber_static_pressure_pa[PY];
            const real Ps0 = g_ambient_pressure_pa;
            const real Tr = Ts / Ts0;
            const real Pr = Ps/ Ps0;
            const real S = 0.4_r * Tr * Tr / sqrt(sqrt(sqrt(Pr)));

            /*
             *                 2
             * Vb = dh * pi * r
             *
             */

            const real randomness = 0.2_r;
            const real dh = S * g_dt_s;
            const real drh = dh * (1.0_r + randomness * frand());
            const real h1 = piston_chamber_flame_height_m;
            const real h2 = h1 + drh;
            const real r = piston_chamber_radius_m;
            const real Vb = (h2 - h1) * g_pi_r * r * r;

            /*
             *      M
             * p = ---
             *      V
             *
             * Mburned = Vb * p
             *
             */

            const real p = M / V;
            const real Mb = Vb * p;

            /*
             *             Q
             * Ts = Ts + ------
             *            M Cv
             *
             */

            const real TMb = piston_chamber_mass_burned_m3 + Mb;
            if(TMb / M < 1.0_r)
            {
                const real MFb = Mb / (1.0_r + g_stoich_air_fuel_ratio);
                const real Q = MFb * g_energy_octane_j_per_kg;
                const real Cv = g_cv_j_per_kg_k;
                const real dTs = Q / (M * Cv);
                chamber_static_temperature_k[PY] += dTs;
            }
            else
            {
                piston_chamber_on_fire = false;
            }
            piston_chamber_flame_height_m = h2;
            piston_chamber_mass_burned_m3 = TMb;
        }
    }

    fn void calc_panics()
    {
        for(size_t i = 0; i < N; i++)
        {
            panic[i] |= chamber_mass_kg[i] <= 0.0_r;
            panic[i] |= chamber_static_temperature_k[i] <= 0.0_r;
            panic[i] |= chamber_static_pressure_pa[i] <= 0.0_r;
        }
    }

    fn void update()
    {
        calc_chamber_dynamic_pressures();
        calc_chamber_dynamic_temperatures();
        calc_chamber_static_pressures();
        calc_chamber_total_pressures();
        calc_chamber_total_temperatures();
        calc_nozzle_machs();
        calc_nozzle_velocities();
        calc_nozzle_static_densities();
        calc_nozzle_real_flow_areas();
        calc_nozzle_mass_flow_rates();
        calc_nozzle_parcels();
        calc_compressions();
        calc_forward_energy_transfers();
        calc_reverse_energy_transfers();
        calc_mass_transfers();
        calc_combustion();
        calc_panics();
    }
};
