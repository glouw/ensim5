#pragma once

template<size_t W>
struct pipe
{
    static constexpr size_t L = 256;
    static constexpr size_t M = L - 1;
    static constexpr size_t substeps = 10;
    static constexpr real dt_s = g_dt_s / static_cast<real>(substeps);
    const real length_m = 1.50_r;
    const real radius_m = 0.01_r;
    const real dx_m = length_m / static_cast<real>(L);
    const real dx_volume_m3 = g_pi_r * radius_m * radius_m * dx_m;
    const real max_wave_speed_m_per_s = 800.0_r;
    const real cfl = max_wave_speed_m_per_s * dt_s / dx_m;
    const real max_cfl = 0.75_r;
    const real dx_dt = dx_m / dt_s;
    const real dt_dx = dt_s / dx_m;
    const real mic_position_ratio = 0.75_r;
    line pipe_pressure;

    virtual void calc_flux_faces() = 0;

    pipe()
    {
        reset();
    }

    pipe(const real length_m, const real radius_m, const real mic_position_ratio)
        : length_m(length_m)
        , radius_m(radius_m)
        , mic_position_ratio(mic_position_ratio)
    {
        reset();
    }

    lane<W> in_mass_flow_rate_kg_per_s = {};
    lane<W> in_velocity_m_per_s = {};
    lane<W> in_static_temperature_k = {};

    /*
     * +------+     +------+     +------+     +------+
     * |  U0  |     |  U1  |     |  U2  | ... |  UL  |
     * +------+     +------+     +------+     +------+
     * +------+     +------+     +------+     +------+
     * |  F0  |     |  F1  |     |  F2  | ... |  FL  |
     * +------+     +------+     +------+     +------+
     *        +-----+      +-----+      +-----+
     *        | Ff0 |      | Ff1 |      | FfM | M == L - 1
     *        +-----+      +-----+      +-----+
     */

    lane<L> U_r = {};
    lane<L> U_ru = {};
    lane<L> U_rE = {};
    lane<L> F_r = {};
    lane<L> F_ru = {};
    lane<L> F_rE = {};
    lane<M> Ff_r = {};
    lane<M> Ff_ru = {};
    lane<M> Ff_rE = {};

    void to_ambient(const size_t i)
    {
        const real r = g_ambient_density_kg_per_m3;
        const real u = 0.0_r;
        const real ru = r * u;
        const real Rs = g_specific_gas_constant_j_per_kg_k;
        const real Ts = g_ambient_temperature_k;
        const real rE = r * Rs * Ts / (g_gamma - 1.0_r) + ru * ru / (2.0_r * r);
        U_r[i] = r;
        U_ru[i] = ru;
        U_rE[i] = rE;
    }

    void reset()
    {
        if(cfl >= max_cfl)
        {
            throw std::runtime_error("cfl condition unmet: " + std::to_string(cfl) + " >= " + std::to_string(max_cfl));
        }
        for(size_t i = 0; i < L; i++)
        {
            to_ambient(i);
        }
    }

    void energize(const size_t i, const real mass_kg, const real momentum_kg_m_per_s, const real total_energy_j)
    {
        if(mass_kg != 0.0_r)
        {
            const real V = dx_volume_m3;
            const real r = mass_kg / V;
            const real p = momentum_kg_m_per_s;
            const real E = total_energy_j / mass_kg;
            U_r[i] += r;
            U_ru[i] += p / V;
            U_rE[i] += r * E;
        }
    }

    void carry(const size_t i)
    {
        U_r [i] = 1.0_r * U_r [i - 1];
        U_ru[i] = 1.0_r * U_ru[i - 1];
        U_rE[i] = 1.0_r * U_rE[i - 1];
    }

    void calc_boundary_conditions()
    {
        real mass_kg = 0.0_r;
        real momentum_kg_m_per_s = 0.0_r;
        real total_energy_j = 0.0_r;
        for(size_t i = 0; i < W; i++)
        {
            const real mdot = in_mass_flow_rate_kg_per_s[i];
            const real u = in_velocity_m_per_s[i];
            const real Ts = in_static_temperature_k[i];
            const real Cv = g_cv_j_per_kg_k;
            const real m = mdot * dt_s;
            const real E = Cv * Ts + 0.5_r * u * u;
            mass_kg += m;
            momentum_kg_m_per_s += m * u;
            total_energy_j += E * m;
        }
        to_ambient(0);
        carry(L - 1);
        energize(0, mass_kg, momentum_kg_m_per_s, total_energy_j);
    }

    /*
     *  F = [ rr \ ruu + Ps \ u(rE + Ps) ]
     */

    void calc_fluxes()
    {
        for(size_t i = 0; i < L; i++)
        {
            const real r = U_r[i];
            const real ru = U_ru[i];
            const real rE = U_rE[i];
            const real u = ru / r;
            const real E = rE / r;
            const real Ps = (g_gamma - 1.0_r) * r * (E - 0.5_r * u * u);
            F_r[i] = ru;
            F_ru[i]= ru * u + Ps;
            F_rE[i] = u * (rE + Ps);
        }
    }

    /*
     *          dt
     * U = U - ---- * [ Ffr - Ffl ]
     *          dx
     */

    void calc_cells()
    {
        for(size_t i = 1; i < L - 1; i++)
        {
            const size_t j = i - 1;
            U_r [i] -= dt_dx * (Ff_r [i] - Ff_r [j]);
            U_ru[i] -= dt_dx * (Ff_ru[i] - Ff_ru[j]);
            U_rE[i] -= dt_dx * (Ff_rE[i] - Ff_rE[j]);
        }
    }

    void update()
    {
        for(size_t i = 0; i < substeps; i++)
        {
            calc_boundary_conditions();
            calc_fluxes();
            calc_flux_faces();
            calc_cells();
        }
    }

    real calc_static_pressure(const size_t i)
    {
        const real r = U_r[i];
        const real ru = U_ru[i];
        const real rE = U_rE[i];
        const real u = ru / r;
        const real E = rE / r;
        const real Ps = (g_gamma - 1.0_r) * r * (E - 0.5_r * u * u);
        return Ps;
    }

    void calc_pipe_pressure()
    {
        pipe_pressure.clear();
        for(size_t i = 0; i < L; i++)
        {
            const real Ps = calc_static_pressure(i);
            pipe_pressure.push_back(Ps);
        }
    }

    real calc_audio_sample()
    {
        const size_t i = mic_position_ratio * (L - 1);
        return calc_static_pressure(i);
    }
};

template<size_t W>
struct lf_pipe : pipe<W>
{
    using pipe<W>::M;
    using pipe<W>::dx_dt;
    using pipe<W>::Ff_r;
    using pipe<W>::Ff_ru;
    using pipe<W>::Ff_rE;
    using pipe<W>::F_r;
    using pipe<W>::F_ru;
    using pipe<W>::F_rE;
    using pipe<W>::U_r;
    using pipe<W>::U_ru;
    using pipe<W>::U_rE;
    using pipe<W>::calc_static_pressure;

    /*
     *       1                 1
     * Ff = --- [ Fl + Fr ] - --- dx/dt * [ Ur - Ul ]
     *       2                 2
     */

    void calc_flux_faces() override
    {
        for(size_t i = 0; i < M; i++)
        {
            const size_t j = i + 1;
            Ff_r [i] = 0.5_r * ((F_r [i] + F_r [j]) - dx_dt * (U_r [j] - U_r [i]));
            Ff_ru[i] = 0.5_r * ((F_ru[i] + F_ru[j]) - dx_dt * (U_ru[j] - U_ru[i]));
            Ff_rE[i] = 0.5_r * ((F_rE[i] + F_rE[j]) - dx_dt * (U_rE[j] - U_rE[i]));
        }
    }
};

template<size_t W>
struct lfr_pipe : pipe<W>
{
    using pipe<W>::M;
    using pipe<W>::dx_dt;
    using pipe<W>::Ff_r;
    using pipe<W>::Ff_ru;
    using pipe<W>::Ff_rE;
    using pipe<W>::F_r;
    using pipe<W>::F_ru;
    using pipe<W>::F_rE;
    using pipe<W>::U_r;
    using pipe<W>::U_ru;
    using pipe<W>::U_rE;
    using pipe<W>::calc_static_pressure;

    void calc_flux_faces() override
    {
        for(size_t i = 0; i < M; i++)
        {
            const size_t j = i + 1;
            const real cl = sqrt(g_gamma * calc_static_pressure(i) / U_r[i]);
            const real cr = sqrt(g_gamma * calc_static_pressure(j) / U_r[j]);
            const real ul = U_ru[i] / U_r[i];
            const real ur = U_ru[j] / U_r[j];
            const real alpha = fmax(fabs(ul) + cl, fabs(ur) + cr);
            Ff_r [i] = 0.5_r * ((F_r [i] + F_r [j]) - alpha * (U_r [j] - U_r [i]));
            Ff_ru[i] = 0.5_r * ((F_ru[i] + F_ru[j]) - alpha * (U_ru[j] - U_ru[i]));
            Ff_rE[i] = 0.5_r * ((F_rE[i] + F_rE[j]) - alpha * (U_rE[j] - U_rE[i]));
        }
    }
};

template<size_t W>
struct hlle_pipe : pipe<W>
{
    using pipe<W>::M;
    using pipe<W>::dx_dt;
    using pipe<W>::Ff_r;
    using pipe<W>::Ff_ru;
    using pipe<W>::Ff_rE;
    using pipe<W>::F_r;
    using pipe<W>::F_ru;
    using pipe<W>::F_rE;
    using pipe<W>::U_r;
    using pipe<W>::U_ru;
    using pipe<W>::U_rE;
    using pipe<W>::calc_static_pressure;

    void calc_flux_faces() override
    {
        for(size_t i = 0; i < M; i++)
        {
            const size_t j = i + 1;
            const real pL = calc_static_pressure(i);
            const real pR = calc_static_pressure(j);
            const real cL = sqrt(g_gamma * pL / U_r[i]);
            const real cR = sqrt(g_gamma * pR / U_r[j]);
            const real uL = U_ru[i] / U_r[i];
            const real uR = U_ru[j] / U_r[j];
            const real SL = fmin(uL - cL, uR - cR);
            const real SR = fmax(uL + cL, uR + cR);
            if(SL >= 0.0_r)
            {
                Ff_r[i]  = F_r[i];
                Ff_ru[i] = F_ru[i];
                Ff_rE[i] = F_rE[i];
            }
            else
            if(SR <= 0.0_r)
            {
                Ff_r[i]  = F_r[j];
                Ff_ru[i] = F_ru[j];
                Ff_rE[i] = F_rE[j];
            }
            else
            {
                const real inv = 1.0_r / (SR - SL);
                Ff_r[i] = inv * (SR * F_r[i] - SL * F_r[j] + SL * SR * (U_r[j] - U_r[i]));
                Ff_ru[i] = inv * (SR * F_ru[i] - SL * F_ru[j] + SL * SR * (U_ru[j] - U_ru[i]));
                Ff_rE[i] = inv * (SR * F_rE[i] - SL * F_rE[j] + SL * SR * (U_rE[j] - U_rE[i]));
            }
        }
    }
};

template<size_t W>
struct hllc_pipe : pipe<W>
{
    using pipe<W>::M;
    using pipe<W>::dx_dt;
    using pipe<W>::Ff_r;
    using pipe<W>::Ff_ru;
    using pipe<W>::Ff_rE;
    using pipe<W>::F_r;
    using pipe<W>::F_ru;
    using pipe<W>::F_rE;
    using pipe<W>::U_r;
    using pipe<W>::U_ru;
    using pipe<W>::U_rE;
    using pipe<W>::calc_static_pressure;

    void calc_flux_faces() override
    {
        for(size_t i = 0; i < M; i++)
        {
            const size_t j = i + 1;
            const real rhoL = U_r[i];
            const real rhoR = U_r[j];
            const real uL = U_ru[i] / rhoL;
            const real uR = U_ru[j] / rhoR;
            const real pL = calc_static_pressure(i);
            const real pR = calc_static_pressure(j);
            const real cL = sqrt(g_gamma * pL / rhoL);
            const real cR = sqrt(g_gamma * pR / rhoR);
            const real SL = fmin(uL - cL, uR - cR);
            const real SR = fmax(uL + cL, uR + cR);
            const real SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / (rhoL * (SL - uL) - rhoR * (SR - uR));
            if(SL >= 0.0_r)
            {
                Ff_r[i]  = F_r[i];
                Ff_ru[i] = F_ru[i];
                Ff_rE[i] = F_rE[i];
            }
            else
            if(SM >= 0.0_r)
            {
                const real rhoSL = rhoL * (SL - uL) / (SL - SM);
                const real pStar = pL + rhoL * (SL - uL) * (SM - uL);
                const real E_L = U_rE[i] / rhoL;
                const real EStar = ((SL - uL) * E_L - pL * uL / rhoL + pStar * SM / rhoSL) / (SL - SM);
                const real Ustar_r  = rhoSL;
                const real Ustar_ru = rhoSL * SM;
                const real Ustar_rE = rhoSL * EStar;
                Ff_r[i]  = F_r[i] + SL * (Ustar_r - U_r[i]);
                Ff_ru[i] = F_ru[i] + SL * (Ustar_ru - U_ru[i]);
                Ff_rE[i] = F_rE[i] + SL * (Ustar_rE - U_rE[i]);
            }
            else
            if(SR > 0.0_r)
            {
                const real rhoSR = rhoR * (SR - uR) / (SR - SM);
                const real pStar = pR + rhoR * (SR - uR) * (SM - uR);
                const real E_R = U_rE[j] / rhoR;
                const real EStar = ((SR - uR) * E_R - pR * uR / rhoR + pStar * SM / rhoSR) / (SR - SM);
                const real Ustar_r  = rhoSR;
                const real Ustar_ru = rhoSR * SM;
                const real Ustar_rE = rhoSR * EStar;
                Ff_r[i] = F_r[j] + SR * (Ustar_r - U_r[j]);
                Ff_ru[i] = F_ru[j] + SR * (Ustar_ru - U_ru[j]);
                Ff_rE[i] = F_rE[j] + SR * (Ustar_rE - U_rE[j]);
            }
            else
            {
                Ff_r[i]  = F_r[j];
                Ff_ru[i] = F_ru[j];
                Ff_rE[i] = F_rE[j];
            }
        }
    }
};
