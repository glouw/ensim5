#pragma once
/*
 * [ ]  ...  [ ] | <- Source
 * [ ]  ...  [ | | <- Intake
 * [ ]  ...  [ ] | <- Intake Manifold
 * [ ]  ...  [ ] | <- Intake Runner
 * [ ]  ...  [ ] H <- Piston (PY)
 * [ ]  ...  [ ] | <- Exhaust Runner
 * [ ]  ...  [ ] | <- Exhaust Manifold
 * [ ]  ...  [ ] | <- Exhaust
 * [ ]  ...  [ ] | <- Sink
 * +---- W ----+ +
 *
 */

template<
    size_t W,
    size_t H,
    size_t TY,
    size_t PY,
    size_t AY,
    template<size_t> class P,
    template<size_t> class C,
    template<size_t> class S>
struct as_engine : engine
{
    struct P<W> pistons;
    struct C<W> inlet_cam;
    struct C<W> outlet_cam;
    struct S<W> sparkplugs;
    std::array<struct flow<H, PY>, W> flows;
    struct limiter limiter;
    struct throttle throttle;
    struct crankshaft crankshaft;
    struct flywheel flywheel;
    struct dc_filter dc;
    struct gain_filter pregain;
    struct gain_filter gain;
    struct clamp_filter clamp;
    struct convolution_filter convolution;
    struct diags diags;
    struct lfr_pipe<W> pipe;
    line audio_signal;
    struct mailbox<W, H> mailbox;
    mutable std::vector<float> audio_data;
    line pipe_pressure;
    std::mutex swap_mutex;

    void log_step(const size_t x, const size_t y)
    {
        if(x < W and y < H)
        {
            #define X(name) diags.back[g_##name].push_back(flows[x].name[y]);
            FLUIDS(X)
            #undef X
            if(y == PY)
            {
                #define X(name) diags.back[g_##name].push_back(pistons.name[x]);
                PISTONS(X)
                #undef X
            }
        }
    }

    fn void broadcast(const real throttle_open_ratio, const bool injection_enabled)
    {
        /*
         * Crankshaft theta -> inlet/outlet cams + pistons + sparkplugs thetas.
         *
         */

        inlet_cam.crankshaft_theta_r = crankshaft.theta_r;
        outlet_cam.crankshaft_theta_r = crankshaft.theta_r;
        pistons.crankshaft_theta_r = crankshaft.theta_r;
        sparkplugs.crankshaft_theta_r = crankshaft.theta_r;
        pistons.crankshaft_angular_velocity_r_per_s = crankshaft.angular_velocity_r_per_s;
        limiter.crankshaft_angular_velocity_r_per_s = crankshaft.angular_velocity_r_per_s;

        /*
         * Cam open ratios -> chamber open ratios.
         *
         */

        for(size_t x = 0; x < W; x++)
        {
            flows[x].chamber_nozzle_open_ratio[PY - 1] = inlet_cam.open_ratio[x];
            flows[x].chamber_nozzle_open_ratio[PY + 0] = outlet_cam.open_ratio[x];
        }

        /*
         * Throttle open ratios -> chamber open ratios.
         *
         */

        for(size_t x = 0; x < W; x++)
        {
            const real open_ratio = throttle.lookup(throttle_open_ratio);
            flows[x].chamber_nozzle_open_ratio[TY] = open_ratio;
        }

        /*
         * Piston shapes <-> flow shapes.
         *
         */

        for(size_t x = 0; x < W; x++)
        {
            flows[x].piston_injection_enabled = injection_enabled;
            flows[x].piston_chamber_radius_m = pistons.diameter_m[x] / 2.0_r;
            flows[x].chamber_volume_m3[PY] = pistons.volumes_m3[x];
            pistons.chamber_static_pressure_pa[x] = flows[x].chamber_static_pressure_pa[PY];
        }

        /*
         * System torque + inertia to system angular acceleration.
         *
         *      t
         * a = ---
         *      I
         */

        real I = 0.0_r;
        for(size_t x = 0; x < W; x++)
        {
            I += pistons.moment_of_inertia_kg_m2[x];
        }
        I += flywheel.moment_of_inertia_kg_m2;
        I += crankshaft.moment_of_inertia_kg_m2;
        real t = 0.0_r;
        for(size_t x = 0; x < W; x++)
        {
            t += pistons.total_torque_n_m[x];
        }
        crankshaft.angular_acceleration_r_per_s2 = t / I;
        crankshaft.angular_velocity_r_per_s = fmax(crankshaft.angular_velocity_r_per_s, 0.0_r);
    }

    void remember_volumes()
    {
        for(size_t x = 0; x < W; x++)
        {
            flows[x].chamber_prev_volume_m3 = flows[x].chamber_volume_m3;
        }
    }

    fn void reset_chambers()
    {
        for(size_t x = 0; x < W; x++)
        {
            flows[x].calc_chamber_ambients();
        }
    }

    void reset() override
    {
        flywheel.update();
        crankshaft.update();
        pistons.calc_volumetrics();
        broadcast(0.0_r, false);
        remember_volumes();
        reset_chambers();
    }

    bool diags_swap()
    {
        if(swap_mutex.try_lock())
        {
            for(size_t i = 0; i < diags.front.size(); i++)
            {
                std::swap(diags.front[i], diags.back[i]);
            }
            for(auto& line : diags.back)
            {
                line.clear();
            }
            pipe.calc_pipe_pressure();
            std::swap(pipe_pressure, pipe.pipe_pressure);
            swap_mutex.unlock();
            return true;
        }
        else
        {
            /*
             * Discard back if front in use by renderer.
             * Addng more samples to back will distort diags oscilloscope trigger.
             */

            for(auto& line : diags.back)
            {
                line.clear();
            }
            return false;
        }
    }

    void sample_audio()
    {
        real x0 = pipe.calc_audio_sample();
        x0 = dc.filter(x0);
        x0 = pregain.filter(x0);
        x0 = convolution.filter(x0);
        x0 = gain.filter(x0);
        x0 = clamp.filter(x0);
        audio_signal.push_back(x0);
    }

    void update_pipe()
    {
        for(size_t x = 0; x < W; x++)
        {
            pipe.in_mass_flow_rate_kg_per_s[x] = flows[x].nozzle_mass_flow_rate_kg_per_s[AY];
            pipe.in_velocity_m_per_s[x] = flows[x].nozzle_velocity_m_per_s[AY];
            pipe.in_static_temperature_k[x] = flows[x].parcel_static_temperature_k[AY];
        }
        pipe.update();
    }

    void update_limiter()
    {
        limiter.update();
    }

    void update_flywheel()
    {
        flywheel.update();
    }

    bool update_crankshaft()
    {
        const bool otto_cycled = crankshaft.update();
        if(otto_cycled)
        {
            return diags_swap();
        }
        return false;
    }

    void update_cams()
    {
        inlet_cam.update();
        outlet_cam.update();
    }

    void update_sparkplugs()
    {
        sparkplugs.update();
    }

    void update_pistons()
    {
        pistons.update();
    }

    void update_ignition()
    {
        for(size_t x = 0; x < W; x++)
        {
            if(sparkplugs.rising_edge[x])
            {
                flows[x].ignite_piston_chamber();
            }
        }
    }

    void update_flows()
    {
        for(size_t x = 0; x < W; x++)
        {
            flows[x].update();
        }
    }

    void post_mailbox(const size_t swap_drops)
    {
        mailbox.engine_angular_velocity_r_per_s = crankshaft.angular_velocity_r_per_s;
        for(size_t y = 0; y < H; y++)
        for(size_t x = 0; x < W; x++)
        {
            mailbox.port_open_ratios[y][x] = flows[x].chamber_nozzle_open_ratio[y];
            mailbox.panics[y][x] = flows[x].panic[y];
        }
        mailbox.swap_drops += swap_drops;
    }

    void run(const size_t steps) override
    {
        const real throttle_open_ratio = mailbox.throttle_open_ratio;
        const size_t log_x = mailbox.log_x;
        const size_t log_y = mailbox.log_y;
        const bool injection_enabled = mailbox.injection_enabled;
        audio_signal.clear();
        audio_signal.reserve(steps);
        size_t swap_drops = 0;
        for(size_t step = 0; step < steps; step++)
        {
            update_limiter();
            update_flywheel();
            if(update_crankshaft())
            {
                swap_drops++;
            }
            update_cams();
            update_sparkplugs();
            update_pistons();
            update_ignition();
            update_flows();
            update_pipe();
            log_step(log_x, log_y);
            remember_volumes();
            const bool injection_overrided = injection_enabled && not limiter.limiting;
            broadcast(throttle_open_ratio, injection_overrided);
            sample_audio();
        }
        post_mailbox(swap_drops);
    }

    size_t get_width() const override
    {
        return W;
    }

    size_t get_height() const override
    {
        return H;
    }

    size_t get_piston_y() const override
    {
        return PY;
    }

    size_t get_audio_y() const override
    {
        return AY;
    }

    size_t get_throttle_y() const override
    {
        return TY;
    }

    size_t get_bytes() const override
    {
        return sizeof *this;
    }

    std::string_view get_signal_name(const size_t index) const override
    {
        return g_signal_names[index];
    }

    const atom& get_angular_velocity_r_per_s() const override
    {
        return mailbox.engine_angular_velocity_r_per_s;
    }

    const atom& get_port_open_ratio(const size_t x, const size_t y) const override
    {
        return mailbox.port_open_ratios[y][x];
    }

    const trip& get_panic(const size_t x, const size_t y) const override
    {
        return mailbox.panics[y][x];
    }

    size_t get_swap_drops() const override
    {
        return mailbox.swap_drops;
    }

    const line& get_signal(const size_t index) const override
    {
        return diags.front[index];
    }

    const line& get_static_temperature_signal_k() const override
    {
        return get_signal(g_chamber_static_temperature_k);
    }

    const line& get_static_pressure_signal_pa() const override
    {
        return get_signal(g_chamber_static_pressure_pa);
    }

    const line& get_volume_signal_m3() const override
    {
        return get_signal(g_chamber_volume_m3);
    }

    const line& get_audio_signal() const override
    {
        return audio_signal;
    }

    const line& get_impulse_signal() const override
    {
        return g_impulse;
    }

    const std::vector<float>& get_audio_data() const override
    {
        audio_data.clear();
        const line& audio_signal = get_audio_signal();
        for(const real& x : audio_signal)
        {
            audio_data.push_back(static_cast<float>(x));
        }
        return audio_data;
    }

    const line& get_pipe_pressure_signal() const override
    {
        return pipe_pressure;
    }

    void set_throttle_open_ratio(const real open_ratio) override
    {
        mailbox.throttle_open_ratio = open_ratio;
    }

    void set_injection_on() override
    {
        mailbox.injection_enabled = true;
    }

    void set_injection_off() override
    {
        mailbox.injection_enabled = false;
    }

    void set_logger(const size_t x, const size_t y) override
    {
        mailbox.log_x = x;
        mailbox.log_y = y;
    }

    void set_swap_lock_on() override
    {
        swap_mutex.lock();
    }

    void set_swap_lock_off() override
    {
        swap_mutex.unlock();
    }
};
