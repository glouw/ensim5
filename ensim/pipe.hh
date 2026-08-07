#pragma once

template<size_t W, size_t L>
struct pipe
{
    static constexpr size_t M = L - 1;
    size_t substeps = 0;
    real dt_s = 0.0_r;
    real length_m = 0.0_r;
    real radius_m = 0.0_r;
    real dx_m = 0.0_r;
    real dx_volume_m3 = 0.0_r;
    real cfl = 0.0_r;
    real dx_dt = 0.0_r;
    real dt_dx = 0.0_r;
    real mic_position_ratio = 0.0_r;
    real reflective_ratio = 0.0_r;
    line pipe_pressure_signal = {};

    void configure(
        const real length_m,
        const real radius_m,
        const real mic_position_ratio,
        const real reflective_ratio,
        const size_t substeps)
    {
        this->substeps = substeps;
        dt_s = g_dt_s / substeps;
        this->length_m = length_m;
        this->radius_m = radius_m;
        dx_m = length_m / static_cast<real>(L);
        dx_volume_m3 = g_pi_r * radius_m * radius_m * dx_m;
        dx_dt = dx_m / dt_s;
        dt_dx = dt_s / dx_m;
        this->mic_position_ratio = mic_position_ratio;
        this->reflective_ratio = reflective_ratio;
        reset();
    }

    /*
     * Pipe Junction - mix mass W flow rates, velocities, static temperatures streams.
     *
     * --- +-------+
     *  |  |  0    | ---+
     *  |  +-------+    |
     *  |  +-------+    |
     *  W  |  1    | ---+ --> U0 ... UL-1
     *  |  +-------+    |
     *  |     ...       |
     *  |  +-------+    |
     *  |  |  W-1  | ---+
     * --- +-------+
     */

    lane<W> in_mass_flow_rate_kg_per_s = {};
    lane<W> in_velocity_m_per_s = {};
    lane<W> in_static_temperature_k = {};

    /*
     * Cells. U is Conserved state: F is flux state: Ff is flux face state.
     *
     * +------+     +------+     +------+       +--------+
     * |  U0  |     |  U1  |     |  U2  | ..... |  UL-1  |
     * +------+     +------+     +------+       +--------+
     * +------+     +------+     +------+       +--------+
     * |  F0  |     |  F1  |     |  F2  | ..... |  FL-1  |
     * +------+     +------+     +------+       +--------+
     *        +-----+      +-----+      +-------+
     *        | Ff0 |      | Ff1 |      | FfM-1 | M = L-1
     *        +-----+      +-----+      +-------+
     *
     * |---------------------- L ------------------------|
     *
     */

    lane<L> U_r = {};
    lane<L> U_ru = {};
    lane<L> U_rEs = {};
    lane<L> F_r = {};
    lane<L> F_ru = {};
    lane<L> F_rEs = {};
    lane<M> Ff_r = {};
    lane<M> Ff_ru = {};
    lane<M> Ff_rEs = {};

    lane<L> static_pressure_pa = {};

    void to_ambient(const size_t i)
    {
        const real r = g_ambient_density_kg_per_m3;
        const real u = 0.0_r;
        const real ru = r * u;
        const real Rs = g_specific_gas_constant_j_per_kg_k;
        const real Ts = g_ambient_temperature_k;
        const real rEs = r * Rs * Ts / (g_gamma - 1.0_r) + ru * ru / (2.0_r * r);
        U_r[i] = r;
        U_ru[i] = ru;
        U_rEs[i] = rEs;
    }

    void reset()
    {
        if(cfl >= 1.0_r)
        {
            throw std::runtime_error("cfl condition unmet: " + std::to_string(cfl));
        }
        for(size_t i = 0; i < L; i++)
        {
            to_ambient(i);
        }
    }

    /*
     *                1   2
     *  Es = Cv Ts + --- u
     *                2
     *  m = m dt
     *  p = m u
     *  E = m Es
     */

    void inject(const size_t cell)
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
            const real Es = Cv * Ts + 0.5_r * u * u;
            mass_kg += m;
            momentum_kg_m_per_s += m * u;
            total_energy_j += m * Es;
        }
        if(mass_kg != 0.0_r)
        {
            const real V = dx_volume_m3;
            const real r = mass_kg / V;
            const real p = momentum_kg_m_per_s;
            const real Es = total_energy_j / mass_kg;
            U_r[cell] += r;
            U_ru[cell] += p / V;
            U_rEs[cell] += r * Es;
        }
    }

    void reflective_carry(const size_t cell)
    {
        U_r[cell] = U_r[cell - 1];
        U_ru[cell] = reflective_ratio * U_ru[cell - 1];
        U_rEs[cell] = U_rEs[cell - 1];
    }

    void calc_boundary_conditions()
    {
        to_ambient(0);
        inject(0);
        reflective_carry(L - 1);
    }

    /*
     *  F = [ rr \ ruu + Ps \ u(rEs + Ps) ]
     */

    void calc_fluxes()
    {
        for(size_t i = 0; i < L; i++)
        {
            const real r = U_r[i];
            const real ru = U_ru[i];
            const real rEs = U_rEs[i];
            const real u = ru / r;
            const real Es = rEs / r;
            const real Ps = (g_gamma - 1.0_r) * r * (Es - 0.5_r * u * u);
            F_r[i] = ru;
            F_ru[i] = ru * u + Ps;
            F_rEs[i] = u * (rEs + Ps);
        }
    }

    /*
     *          dt
     * U = U - ---- * [ Ffr - Ffl ]
     *          dx
     */

    void calc_conserved()
    {
        for(size_t i = 1; i < L - 1; i++)
        {
            const size_t j = i - 1;
            U_r  [i] -= dt_dx * (Ff_r  [i] - Ff_r  [j]);
            U_ru [i] -= dt_dx * (Ff_ru [i] - Ff_ru [j]);
            U_rEs[i] -= dt_dx * (Ff_rEs[i] - Ff_rEs[j]);
        }
    }

    /*
     *                            1   2
     * Ps = (y  - 1) * p * [ E - --- u ]
     *                            2
     */

    real calc_static_pressure(const size_t i)
    {
        const real r = U_r[i];
        const real ru = U_ru[i];
        const real rEs = U_rEs[i];
        const real u = ru / r;
        const real Es = rEs / r;
        const real Ps = (g_gamma - 1.0_r) * r * (Es - 0.5_r * u * u);
        return Ps;
    }

    void calc_static_pressures()
    {
        for(size_t i = 0; i < L; i++)
        {
            const real Ps = calc_static_pressure(i);
            static_pressure_pa[i] = Ps;
        }
    }

    virtual void calc_flux_faces() = 0;

    void update()
    {
        for(size_t i = 0; i < substeps; i++)
        {
            calc_boundary_conditions();
            calc_static_pressures();
            calc_fluxes();
            calc_flux_faces();
            calc_conserved();
        }
    }

    void gather_pipe_pressure_signal()
    {
        pipe_pressure_signal.clear();
        for(size_t i = 0; i < L; i++)
        {
            const real Ps = static_pressure_pa[i];
            pipe_pressure_signal.push_back(Ps);
        }
    }

    real calc_audio_sample()
    {
        const size_t i = mic_position_ratio * (L - 1);
        return static_pressure_pa[i];
    }
};

template<size_t W, size_t L>
struct lf_pipe : pipe<W, L>
{
    using pipe<W, L>::M;
    using pipe<W, L>::dx_dt;
    using pipe<W, L>::Ff_r;
    using pipe<W, L>::Ff_ru;
    using pipe<W, L>::Ff_rEs;
    using pipe<W, L>::F_r;
    using pipe<W, L>::F_ru;
    using pipe<W, L>::F_rEs;
    using pipe<W, L>::U_r;
    using pipe<W, L>::U_ru;
    using pipe<W, L>::U_rEs;

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
            Ff_r  [i] = 0.5_r * ((F_r  [i] + F_r  [j]) - dx_dt * (U_r  [j] - U_r  [i]));
            Ff_ru [i] = 0.5_r * ((F_ru [i] + F_ru [j]) - dx_dt * (U_ru [j] - U_ru [i]));
            Ff_rEs[i] = 0.5_r * ((F_rEs[i] + F_rEs[j]) - dx_dt * (U_rEs[j] - U_rEs[i]));
        }
    }
};

template<size_t W, size_t L>
struct lfr_pipe : pipe<W, L>
{
    using pipe<W, L>::M;
    using pipe<W, L>::dx_dt;
    using pipe<W, L>::Ff_r;
    using pipe<W, L>::Ff_ru;
    using pipe<W, L>::Ff_rEs;
    using pipe<W, L>::F_r;
    using pipe<W, L>::F_ru;
    using pipe<W, L>::F_rEs;
    using pipe<W, L>::U_r;
    using pipe<W, L>::U_ru;
    using pipe<W, L>::U_rEs;
    using pipe<W, L>::static_pressure_pa;

    lane<L> speed_of_sound_m_per_s = {};
    lane<L> local_speed_of_sound_m_per_s = {};
    lane<L> absolute_speed_of_sound_m_per_s = {};

    fn void calc_speed_of_sounds()
    {
        for(size_t i = 0; i < L; i++)
        {
            const real Ps = static_pressure_pa[i];
            const real C = sqrt(g_gamma * Ps / U_r[i]);
            speed_of_sound_m_per_s[i] = C;
        }
    }

    fn void calc_local_speed_of_sounds()
    {
        for(size_t i = 0; i < L; i++)
        {
            local_speed_of_sound_m_per_s[i] = U_ru[i] / U_r[i];
        }
    }

    fn void calc_absolute_speed_of_sounds()
    {
        for(size_t i = 0; i < L; i++)
        {
            const real U = local_speed_of_sound_m_per_s[i];
            const real C = speed_of_sound_m_per_s[i];
            absolute_speed_of_sound_m_per_s[i] = fabs(U) + C;
        }
    }

    /*
     *       1                 1
     * Ff = --- [ Fl + Fr ] - --- alpha * [ Ur - Ul ]
     *       2                 2
     */

    fn void calc_flux_faces() override
    {
        calc_speed_of_sounds();
        calc_local_speed_of_sounds();
        calc_absolute_speed_of_sounds();
        for(size_t i = 0; i < M; i++)
        {
            const size_t j = i + 1;
            const real Al = absolute_speed_of_sound_m_per_s[i];
            const real Ar = absolute_speed_of_sound_m_per_s[j];
            const real alpha = fmax(Al, Ar);
            Ff_r  [i] = 0.5_r * ((F_r  [i] + F_r  [j]) - alpha * (U_r  [j] - U_r  [i]));
            Ff_ru [i] = 0.5_r * ((F_ru [i] + F_ru [j]) - alpha * (U_ru [j] - U_ru [i]));
            Ff_rEs[i] = 0.5_r * ((F_rEs[i] + F_rEs[j]) - alpha * (U_rEs[j] - U_rEs[i]));
        }
    }
};
