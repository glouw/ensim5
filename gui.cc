#include "ensim.hh"

#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <atomic>
#include <mutex>
#include <numeric>
#include <raylib.h>
#include <raymath.h>
#include <sstream>
#include <thread>
#include <vector>

static constexpr int g_fps = 60;
static constexpr int g_audio_bit_depth = 32;
static constexpr int g_audio_channels = 1;
static constexpr float g_margin_p = 8.0f;
static constexpr size_t g_xres_p = 1920;
static constexpr size_t g_yres_p = 1080;
static constexpr float g_xmid_p = g_xres_p / 2.0f;
static constexpr const char* g_name = "ensim5";
static constexpr float g_sidebar_width_p = g_xres_p / 5;
static constexpr size_t g_producer_block_size = 200;
static constexpr size_t g_mesh_slices = 6;
static constexpr float g_font_size = 1.0f;
static constexpr Vector3 g_render_scale = { 1.0f, 1.0f, 1.0f };
static constexpr size_t g_stream_buffer_size = 4096;
static constexpr size_t g_producer_queue_size = 4096;
static constexpr size_t g_history_size = 512;
static constexpr Color g_plot_signal_colors[] = { RED, BLUE, GREEN, PURPLE, BROWN, YELLOW };
static constexpr float g_grid_steps = 16.0f;
static constexpr float g_grid_step_size = 0.15f;
static constexpr size_t g_signal_samples = g_sidebar_width_p;
static constexpr size_t g_signal_default_samples = -1;
static constexpr float g_hud_gauge_width_p = 192.0f;

template<typename T, size_t N>
class ring
{
    std::array<T, N> self = {};
    size_t head = 0;
    size_t tail = 0;
    size_t elems = 0;

public:
    void clear()
    {
        head = tail = elems = 0;
    }

    size_t size() const
    {
        return elems;
    }

    size_t capacity() const
    {
        return N;
    }

    void push_back(const T& value)
    {
        self[tail++] = value;
        tail %= N;
        if(elems == N)
        {
            head++;
            head %= N;
        }
        else
        {
            elems++;
        }
    }

    bool empty() const
    {
        return elems == 0;
    }

    void pop_front()
    {
        if(empty())
        {
            throw std::runtime_error("ring empty!");
        }
        head++;
        head %= N;
        elems--;
    }

    size_t index(const size_t i) const
    {
        return (head + i) % N;
    }

    T& operator[](const size_t i)
    {
        return self[index(i)];
    }

    const T& operator[](const size_t i) const
    {
        return self[index(i)];
    }
};

class producer
{
    ensim::engine* engine = nullptr;
    ring<float, g_producer_queue_size> queue = {};
    ring<size_t, g_history_size> history = {};
    ring<size_t, g_history_size> snapshot = {};
    std::mutex engine_mutex = {};
    std::mutex queue_mutex = {};
    std::mutex history_mutex = {};
    std::atomic<bool> done = false;
    std::thread thread = std::thread(&producer::run, this);
    std::vector<float> ready = {};

public:
    std::vector<float>& consume(const size_t samples)
    {
        std::lock_guard lock1(queue_mutex);
        if(samples > queue.size())
        {
            throw std::runtime_error("producer underrun");
        }
        ready.clear();
        for(size_t i = 0; i < samples; i++)
        {
            ready.push_back(queue[0]);
            queue.pop_front();
        }
        std::lock_guard lock2(history_mutex);
        history.push_back(queue.size());
        return ready;
    }

    bool produce()
    {
        std::lock_guard lock(queue_mutex);
        if(queue.size() + g_producer_block_size < queue.capacity())
        {
            std::lock_guard lock(engine_mutex);
            if(engine)
            {
                engine->run(g_producer_block_size);
                const std::vector<float>& slice = engine->get_audio_signal();
                for(float value : slice)
                {
                    queue.push_back(value);
                }
                return true;
            }
        }
        return false;
    }

    void run()
    {
        while(not done)
        {
            if(not produce())
            {
                const auto delay = std::chrono::microseconds(50);
                std::this_thread::sleep_for(delay);
            }
        }
    }

    void set(ensim::engine* other)
    {
        std::lock_guard lock(engine_mutex);
        engine = other;
    }

    void stop()
    {
        done = true;
        thread.join();
    }

    const auto& get_history()
    {
        std::lock_guard lock1(history_mutex);
        snapshot = history;
        return snapshot;
    }
};

producer g_producer = {};

class audio
{
    AudioStream stream = {};

public:
    audio(const AudioCallback audio_callback)
    {
        InitAudioDevice();
        SetAudioStreamBufferSizeDefault(g_stream_buffer_size);
        stream = LoadAudioStream(ensim::g_sample_rate_hz, g_audio_bit_depth, g_audio_channels);
        SetAudioStreamCallback(stream, audio_callback);
        PlayAudioStream(stream);
    }

    ~audio()
    {
        UnloadAudioStream(stream);
        CloseAudioDevice();
    }
};

class widget
{
protected:
    Rectangle rectangle = {};

public:
    void place(const float x_p, const float y_p, const float w_p, const float h_p)
    {
        rectangle.x = x_p;
        rectangle.y = y_p;
        rectangle.width = w_p;
        rectangle.height = h_p;
    }

    virtual void draw() const = 0;
    virtual ~widget() = default;
};

template<typename T>
using sample = std::function<T()>;

class tachometer : public widget
{
    sample<float> value = {};
    sample<float> max = {};
    size_t ticks = 0;
    static constexpr float start_r = 3.0f * std::numbers::pi_v<float> / 4.0f;
    static constexpr float sweep_r = 3.0f * std::numbers::pi_v<float> / 2.0f;
    static constexpr float inner_radius_ratio = 0.82f;
    static constexpr float outer_radius_ratio = 0.92f;
    static constexpr float needle_radius_ratio = 0.80f;
    static constexpr float tick_radius_ratio = 0.70f;
    static constexpr float tick_line_width_p = 2.0f;
    static constexpr float needle_thickness = 3.0f;
    static constexpr float gauge_thickness = 5.0f;

public:
    tachometer(const sample<float> value, const sample<float> max, const size_t ticks)
        : value(value)
        , max(max)
        , ticks(ticks)
    {
    }

    Vector2 to_polar(const float x, const float y, const float angle_r, const float radius_p) const
    {
        return {
            x + std::cos(angle_r) * radius_p,
            y + std::sin(angle_r) * radius_p,
        };
    }

    void draw() const override
    {
        const float radius_p = rectangle.width / 2.0f;
        const Vector2 mid = {
            rectangle.x + radius_p,
            rectangle.y + radius_p,
        };
        DrawCircle(mid.x, mid.y, radius_p, BLACK);
        DrawCircleLines(mid.x, mid.y, radius_p, WHITE);
        for(size_t i = 0; i <= ticks; i++)
        {
            const float tick = static_cast<float>(i) / ticks;
            const float angle_r = start_r + tick * sweep_r;
            const float inner_radius_p = radius_p * inner_radius_ratio;
            const float outer_radius_p = radius_p * outer_radius_ratio;
            DrawLineEx(
                to_polar(mid.x, mid.y, angle_r, inner_radius_p),
                to_polar(mid.x, mid.y, angle_r, outer_radius_p),
                tick_line_width_p,
                WHITE
            );
            const int tick_value = max() * tick;
            const std::string text = std::to_string(tick_value);
            const float tick_radius_p = radius_p * tick_radius_ratio;
            const Vector2 tick_text = to_polar(mid.x, mid.y, angle_r, tick_radius_p);
            const int width_p = MeasureText(text.c_str(), 1);
            DrawText(
                text.c_str(),
                tick_text.x - width_p / 2.0,
                tick_text.y - 5.0f, // TODO: Get font height somehow
                g_font_size,
                WHITE
            );
        }
        const float ratio = value() / max();
        const float angle_r = start_r + ratio * sweep_r;
        const float needle_radius_p = radius_p * needle_radius_ratio;
        DrawLineEx(mid, to_polar(mid.x, mid.y, angle_r, needle_radius_p), needle_thickness, RED);
        DrawCircle(mid.x, mid.y, gauge_thickness, WHITE);
    }
};

template<typename T>
using signal = std::function<const T&(size_t)>;

template<typename T>
class plot : public widget
{
    std::string title = {};
    signal<T> x = {};
    signal<T> y = {};
    size_t samples = {};
    size_t signals = {};

public:
    plot(const std::string& title, const signal<T> x, const signal<T> y, const size_t samples, const size_t signals)
        : title(title)
        , x(x)
        , y(y)
        , samples(samples)
        , signals(signals) {}

    struct range
    {
        double xmin = std::numeric_limits<double>::max();
        double xmax = std::numeric_limits<double>::lowest();
        double ymin = std::numeric_limits<double>::max();
        double ymax = std::numeric_limits<double>::lowest();

        bool valid() const
        {
            return xmin <= xmax and ymin <= ymax;
        }
    };

    range get_size() const
    {
        range range = {};
        for(size_t signal = 0; signal < signals; signal++)
        {
            const T& sx = x(signal);
            const T& sy = y(signal);
            for(size_t i = 0; i < sx.size(); i++)
            {
                const double xv = sx[i];
                const double yv = sy[i];
                range.xmin = std::min(range.xmin, xv);
                range.xmax = std::max(range.xmax, xv);
                range.ymin = std::min(range.ymin, yv);
                range.ymax = std::max(range.ymax, yv);
            }
        }
        return range;
    }

    void draw() const override
    {
        const range range = get_size();
        if(range.valid())
        {
            const double xrange = range.xmax - range.xmin;
            const double yrange = range.ymax - range.ymin;
            std::vector<Vector2> points = {};
            const double x_p = rectangle.x + g_margin_p;
            const double y_p = rectangle.y + g_margin_p;
            size_t width = 0;
            for(size_t signal = 0; signal < signals; signal++)
            {
                const T& sx = x(signal);
                const T& sy = y(signal);
                const size_t size = sx.size();
                width = samples > size ? size : samples;
                points.resize(width);
                for (size_t at = 0; at < width; at++)
                {
                    const size_t i = at * (size - 1) / (width - 1);
                    const double nx_p = xrange > 0.0 ? (sx[i] - range.xmin) / xrange : 0.5;
                    const double ny_p = yrange > 0.0 ? (sy[i] - range.ymin) / yrange : 0.5;
                    const double w_p = -2.0 * g_margin_p + rectangle.width;
                    const double h_p = -2.0 * g_margin_p + rectangle.height;
                    points[at].x = x_p + nx_p * w_p;
                    points[at].y = y_p + (1.0 - ny_p) * h_p;
                }
                DrawLineStrip(points.data(), points.size(), g_plot_signal_colors[signal]);
            }
            const double precision = 6.0;
            const double div = range.ymax / range.ymin;
            std::ostringstream metrics;
            metrics
                << title
                << "\nmax: " << std::fixed << std::setprecision(precision) << range.ymax
                << "\nmin: " << std::fixed << std::setprecision(precision) << range.ymin
                << "\nrng: " << std::fixed << std::setprecision(precision) << range.ymax - range.ymin
                << "\ndiv: " << std::fixed << std::setprecision(precision) << div
                << "\nsamples: " << width;
            DrawText(metrics.str().data(), x_p, y_p, g_font_size, WHITE);
        }
    }
};

template<typename T>
const T& lingen(const size_t size)
{
    static T linear;
    linear.clear();
    for(size_t i = 0; i < size; i++)
    {
        linear.push_back(i);
    }
    return linear;
}

template<typename T>
std::unique_ptr<plot<T>> make_plot(
    const std::string& title,
    const signal<T> y,
    const size_t samples,
    const size_t signals = 1)
{
    const auto x = [y](const size_t signal)-> auto& { return lingen<T>(y(signal).size()); };
    return std::make_unique<plot<T>>(title, x, y, samples, signals);
}

template<typename T>
std::unique_ptr<plot<T>> make_plot(
    const std::string& title,
    const signal<T> x,
    const signal<T> y,
    const size_t samples,
    const size_t signals = 1)
{
    return std::make_unique<plot<T>>(title, x, y, samples, signals);
}

class widgets
{
    std::vector<std::unique_ptr<widget>> self = {};
    ensim::engine* engine = {};

public:
    auto begin()
    {
        return self.begin();
    }

    const auto begin() const
    {
        return self.begin();
    }

    const auto end() const
    {
        return self.end();
    }

    void set(ensim::engine* engine)
    {
        this->engine = engine;
        regen();
    }

    void regen_hud()
    {
        std::unique_ptr<widget> engine_angular_velocity = std::make_unique<tachometer>(
            [this]()-> auto { return engine->get_engine_angular_velocity_r_per_s(); },
            [this]()-> auto { return engine->get_limiter_angular_velocity_r_per_s(); },
            10
        );
        std::unique_ptr<widget> load_angular_velocity = std::make_unique<tachometer>(
            [this]()-> auto { return engine->get_load_angular_velocity_r_per_s(); },
            [this]()-> auto { return 300.0; },
            10
        );
        std::unique_ptr<widget> gear = std::make_unique<tachometer>(
            [this]()-> auto { return engine->get_gear(); },
            [this]()-> auto { return 8.0; },
            8
        );
        const float margin_p = 4.0 * g_margin_p;
        const float w_p = g_hud_gauge_width_p;
        engine_angular_velocity->place(
            g_xmid_p - margin_p - w_p,
            g_yres_p - w_p,
            w_p,
            w_p
        );
        load_angular_velocity->place(
            g_xmid_p + margin_p,
            g_yres_p - w_p,
            w_p,
            w_p
        );
        const float gear_w_p = g_hud_gauge_width_p / 1.75;
        gear->place(
            g_xmid_p - gear_w_p / 2.0,
            g_yres_p - gear_w_p,
            gear_w_p,
            gear_w_p
        );
        self.push_back(std::move(engine_angular_velocity));
        self.push_back(std::move(load_angular_velocity));
        self.push_back(std::move(gear));
    }

    void regen_right_sidebar()
    {
        std::vector<std::unique_ptr<widget>> right;
        const size_t size = engine->get_signal_count();
        for(size_t i = 0; i < size; i++)
        {
            right.push_back(
                make_plot<std::vector<double>>(
                    std::string(engine->get_signal_name(i)),
                    [this, i](auto)-> auto& { return engine->get_signal(i); },
                    g_signal_samples
                )
            );
        }
        float y_p = 0.0f;
        const float x_p = g_xres_p - g_sidebar_width_p;
        const float h_p = g_yres_p / size;
        for(auto& widget : right)
        {
            widget->place(x_p, y_p, g_sidebar_width_p, h_p);
            y_p += h_p;
            self.push_back(std::move(widget));
        }
    }

    void regen_left_sidebar()
    {
        std::vector<std::unique_ptr<widget>> left;
        left.push_back(
            make_plot<std::vector<float>>(
                "pipe pressure (pascals)",
                [this](auto signal)-> auto& {
                    return engine->get_pipe_pressure_signal(signal);
                },
                g_signal_default_samples,
                engine->get_pipe_count()
            )
        );
        left.push_back(
            make_plot<ring<size_t, g_history_size>>(
                "producer audio queue size",
                [](auto)-> auto& {
                    return g_producer.get_history();
                },
                g_signal_default_samples
            )
        );
        left.push_back(
            make_plot<std::vector<double>>(
                "pressure-volume (pascals-m3) graph",
                [this](auto)-> auto& {
                    return engine->get_volume_signal_m3();
                },
                [this](auto)-> auto& {
                    return engine->get_static_pressure_signal_pa();
                },
                g_signal_samples
            )
        );
        left.push_back(
            make_plot<std::vector<double>>(
                "temperature-volume (kelvin-m3) graph",
                [this](auto)-> auto& {
                    return engine->get_volume_signal_m3();
                },
                [this](auto)-> auto& {
                    return engine->get_static_temperature_signal_k();
                },
                g_signal_samples
            )
        );
        float y_p = 0.0f;
        const float x_p = 0.0f;
        const float h_p = g_yres_p / left.size();
        for(auto& widget : left)
        {
            widget->place(x_p, y_p, g_sidebar_width_p, h_p);
            y_p += h_p;
            self.push_back(std::move(widget));
        }
    }

    void regen()
    {
        self.clear();
        regen_left_sidebar();
        regen_right_sidebar();
        regen_hud();
    }
};

struct shape
{
    Model model = {};
    Color color = {};
    Vector3 rotation = {};
    Vector3 position = {};
    double theta_degrees = {};
};

class part
{
protected:
    float x = 0.0;
    float z = 0.0;
    size_t x_id = 0;
    size_t y_id = 0;

    part(const float x, const float z):
        x(x * g_grid_step_size),
        z(z * g_grid_step_size),
        x_id(x),
        y_id(z)
        {
        }

public:
    size_t get_x_id() { return x_id; };
    size_t get_y_id() { return y_id; };
    virtual std::span<const shape> get_shapes() const = 0;
    virtual void update() = 0;
    virtual ~part() = default;
};

class chamber : public part
{
    float radius_m = g_grid_step_size / 4.0f;
    float height_m = 0.0f;
    std::array<shape, 1> shapes = {};

public:
    chamber(
        const float x,
        const float z,
        const float volume_m3)
        : part(x, z)
        , height_m(volume_m3 / (std::numbers::pi_v<float> * radius_m * radius_m))
        {
            shapes[0].model = LoadModelFromMesh(GenMeshCylinder(radius_m, height_m, g_mesh_slices));
            shapes[0].color = DARKBLUE;
        }

    void update() override
    {
        shapes[0].position = { x, 0.0f, z };
    }

    std::span<const shape> get_shapes() const override
    {
        return shapes;
    }
};

class piston : public part
{
    float head_height_m = 0.0f;
    float conrod_width_m = 0.015f;
    float conrod_depth_m = 0.007f;
    float counter_weight_height_m = 0.0f;
    float counter_weight_depth_m = 2.0f * conrod_depth_m;
    float counter_weight_width_m = 2.0f * conrod_width_m;
    std::array<shape, 3> shapes = {};
    std::function<float()> get_pin_y_m = {};
    std::function<float()> get_pin_phi_r = {};
    std::function<float()> get_crank_theta_r = {};

public:
    piston(
        const float x,
        const float z,
        const float head_radius_m,
        const float head_height_m,
        const float conrod_height_m,
        const float crank_diameter_m,
        const std::function<float()> get_pin_y_m,
        const std::function<float()> get_pin_phi_r,
        const std::function<float()> get_crank_theta_r)
        : part(x, z)
        , head_height_m(head_height_m)
        , counter_weight_height_m(1.2f * crank_diameter_m)
        , get_pin_y_m(get_pin_y_m)
        , get_pin_phi_r(get_pin_phi_r)
        , get_crank_theta_r(get_crank_theta_r)
    {
        shapes[0].model = LoadModelFromMesh(GenMeshCylinder(head_radius_m, head_height_m, g_mesh_slices));
        shapes[1].model = LoadModelFromMesh(GenMeshCube(conrod_width_m, conrod_height_m, conrod_depth_m));
        shapes[2].model = LoadModelFromMesh(GenMeshCube(counter_weight_width_m, counter_weight_height_m, counter_weight_depth_m));

        shapes[0].color = DARKBLUE;
        shapes[1].color = DARKBLUE;
        shapes[2].color = DARKBLUE;

        const float rot90_r = 90.0f * DEG2RAD;
        shapes[0].model.transform = MatrixRotateY(rot90_r);
        shapes[1].model.transform = MatrixMultiply(MatrixRotateY(rot90_r), MatrixTranslate(0.0f, -conrod_height_m * 0.5f, 0.0f));
        shapes[2].model.transform = MatrixRotateY(rot90_r);
    }

    void update() override
    {
        shapes[0].position = { x, get_pin_y_m() - head_height_m / 2.0f, z };
        shapes[1].position = { x, get_pin_y_m(), z };
        shapes[2].position = { x - conrod_depth_m, 0.0f, z };
        shapes[1].rotation = { 1.0f, 0.0f, 0.0f };
        shapes[2].rotation = { 1.0f, 0.0f, 0.0f };
        shapes[1].theta_degrees = RAD2DEG * get_pin_phi_r();
        shapes[2].theta_degrees = RAD2DEG * get_crank_theta_r();
    }

    std::span<const shape> get_shapes() const override
    {
        return shapes;
    }

    ~piston()
    {
        for(const auto& shape : shapes)
        {
            UnloadModel(shape.model);
        }
    }
};

class parts
{
    std::vector<std::unique_ptr<part>> self = {};
    ensim::engine* engine = nullptr;

public:
    auto begin()
    {
        return self.begin();
    }

    const auto begin() const
    {
        return self.begin();
    }

    auto end()
    {
        return self.end();
    }

    const auto end() const
    {
        return self.end();
    }

    void set(ensim::engine* engine)
    {
        self.clear();
        this->engine = engine;
        const size_t w = engine->get_width();
        const size_t h = engine->get_height();
        for(size_t x = 0; x < w; x++)
        for(size_t y = 0; y < h; y++)
        {
            const float xx = x;
            const float zz = y;
            if(y == engine->get_source_y())
            {
                // Do not draw
            }
            else
            if(y == engine->get_sink_y())
            {
                // Do not draw
            }
            else
            if(y == engine->get_piston_y())
            {
                self.push_back(std::make_unique<piston>(
                    xx,
                    zz,
                    engine->get_piston_head_radius_m(x),
                    engine->get_piston_head_height_m(x),
                    engine->get_piston_connecting_rod_length_m(x),
                    engine->get_piston_crank_diameter_m(x),
                    [this, x]()-> float {
                        return this->engine->get_piston_pin_y_m(x);
                    },
                    [this, x]()-> float {
                        return this->engine->get_piston_pin_phi_r(x);
                    },
                    [this, x]()-> float {
                        return this->engine->get_piston_crank_theta_r(x);
                    }
                ));
            }
            else
            {
                self.push_back(std::make_unique<chamber>(xx, zz, engine->get_chamber_volume_m3(x, y)));
            }
        }
    }

    void update()
    {
        for(const auto& part : self)
        {
            part->update();
        }
    }
};

struct look
{
    float yaw_r = std::numbers::pi_v<float> / 4.0f;
    float pitch_r = atan(sqrt(0.5f));
    float x = 0.0f;
    float y = 0.0f;
    float distance = 2.0f;
    const float max_pitch_r = std::numbers::pi_v<float> / 2.01f;
    const float min_pitch_r = 0.0f;

    Camera3D camera = {
        .position = {},
        .target = {},
        .up = { 0.0f, 1.0f, 0.0f },
        .fovy = 45.0f,
        .projection = CAMERA_PERSPECTIVE,
    };

    void update()
    {
        pitch_r = std::clamp(pitch_r, min_pitch_r, max_pitch_r);
        const float xx = x * g_grid_step_size;
        const float yy = 0.0f;
        const float zz = y * g_grid_step_size;
        camera.target.x = xx;
        camera.target.y = yy;
        camera.target.z = zz;
        camera.position.x = xx + distance * sinf(yaw_r) * cosf(pitch_r);
        camera.position.y = yy + distance * sinf(pitch_r);
        camera.position.z = zz + distance * cosf(yaw_r) * cosf(pitch_r);
    }
};

class window
{
    ensim::engine* engine = {};
    look look = {};
    float throttle_ratio = 0.0f;

public:
    window()
    {
        SetConfigFlags(FLAG_FULLSCREEN_MODE | FLAG_VSYNC_HINT);
        InitWindow(g_xres_p, g_yres_p, g_name);
        SetTargetFPS(g_fps);
        HideCursor();
    }

    ~window()
    {
        CloseWindow();
    }

    void draw(const widgets& widgets)
    {
        for(const auto& widget : widgets)
        {
            widget->draw();
        }
    }

    void draw(const std::span<const shape> shapes, const bool highlight)
    {
        for(const auto& shape : shapes)
        {
            const Color color = highlight ? RED : LIGHTGRAY;
            DrawModelWiresEx(shape.model, shape.position, shape.rotation, shape.theta_degrees, g_render_scale, color);
        }
    }

    void draw(const parts& parts)
    {
        BeginMode3D(look.camera);
        DrawGrid(g_grid_steps, g_grid_step_size);
        for(const auto& part : parts)
        {
            const bool highlight = part->get_x_id() == look.x && part->get_y_id() == look.y;
            draw(part->get_shapes(), highlight);
        }
        EndMode3D();
    }

    void draw(const widgets& widgets, const parts& parts)
    {
        BeginDrawing();
        ClearBackground(BLACK);
        engine->set_swap_lock_on();
        draw(parts);
        draw(widgets);
        engine->set_swap_lock_off();
        EndDrawing();
    }

    void handle_input()
    {
        if(IsKeyPressed(KEY_ZERO))
        {
            engine->set_throttle_open_ratio(throttle_ratio = 0.00f);
            engine->set_injection_off();
        }
        if(IsKeyPressed(KEY_ONE))
        {
            engine->set_throttle_open_ratio(throttle_ratio = 0.00f);
            engine->set_injection_on();
        }
        if(IsKeyPressed(KEY_TWO))
        {
            engine->set_throttle_open_ratio(throttle_ratio = 0.33f);
            engine->set_injection_on();
        }
        if(IsKeyPressed(KEY_THREE))
        {
            engine->set_throttle_open_ratio(throttle_ratio = 0.66f);
            engine->set_injection_on();
        }
        if(IsKeyPressed(KEY_FOUR))
        {
            engine->set_throttle_open_ratio(throttle_ratio = 0.99f);
            engine->set_injection_on();
        }
        if(IsKeyPressed(KEY_PERIOD))
        {
            engine->set_throttle_open_ratio(0.0);
            engine->disengage_clutch();
            engine->increment_gear();
        }
        if(IsKeyPressed(KEY_COMMA))
        {
            engine->set_throttle_open_ratio(0.0);
            engine->disengage_clutch();
            engine->decrement_gear();
        }
        if(IsKeyReleased(KEY_PERIOD))
        {
            engine->set_throttle_open_ratio(throttle_ratio);
            engine->engage_clutch();
        }
        if(IsKeyReleased(KEY_COMMA))
        {
            engine->set_throttle_open_ratio(throttle_ratio);
            engine->engage_clutch();
        }
        if(IsKeyDown(KEY_Q))
        {
            look.distance += 0.05f;
        }
        if(IsKeyDown(KEY_E))
        {
            look.distance -= 0.05f;
        }
        if(IsKeyDown(KEY_H))
        {
            look.yaw_r -= 0.1f;
        }
        if(IsKeyDown(KEY_L))
        {
            look.yaw_r += 0.1f;
        }
        if(IsKeyDown(KEY_J))
        {
            look.pitch_r -= 0.1f;
        }
        if(IsKeyDown(KEY_K))
        {
            look.pitch_r += 0.1f;
        }
        int dx = 0;
        int dy = 0;
        if(IsKeyPressed(KEY_W))
        {
            dy = -1;
        }
        else
        if(IsKeyPressed(KEY_S))
        {
            dy = 1;
        }
        else
        if(IsKeyPressed(KEY_A))
        {
            dx = -1;
        }
        else
        if(IsKeyPressed(KEY_D))
        {
            dx = 1;
        }
        const int dirs = 4;
        const int dir = std::round(look.yaw_r / (std::numbers::pi_v<float> / 2.0f));
        const int id = ((dir % dirs) + dirs) % dirs;
        int mx = 0;
        int my = 0;
        switch(id)
        {
            case 0: mx =  dx; my =  dy; break;
            case 1: mx =  dy; my = -dx; break;
            case 2: mx = -dx; my = -dy; break;
            case 3: mx = -dy; my =  dx; break;
        }
        look.x += mx;
        look.y += my;
    }

    void handle_look()
    {
        look.update();
        engine->set_logger(look.x, look.y);
    }

    void loop(const widgets& widgets, parts& parts)
    {
        while(not done())
        {
            handle_input();
            handle_look();
            parts.update();
            draw(widgets, parts);
        }
    }

    void set(ensim::engine* engine)
    {
        this->engine = engine;
        look.x = engine->get_width() / 2;
        look.y = engine->get_piston_y();
    }

    bool done()
    {
        return WindowShouldClose();
    }
};

static void set(window& window, widgets& widgets, producer& producer, parts& parts, ensim::engine* engine)
{
    window.set(engine);
    widgets.set(engine);
    producer.set(engine);
    parts.set(engine);
}

static void audio_callback(void* const data, const unsigned frames)
{
    float* const buffer = static_cast<float* const>(data);
    const std::vector<float>& block = g_producer.consume(frames);
    std::copy(block.begin(), block.end(), buffer);
};

int main()
{
    std::unique_ptr<ensim::engine> engine = ensim::new_engine(ensim::type::inline8);
    window window;
    widgets widgets;
    audio audio(audio_callback);
    parts parts;
    set(window, widgets, g_producer, parts, engine.get());
    window.loop(widgets, parts);
    g_producer.stop();
}

// TODO
// * add semicircle for counter balance mass
