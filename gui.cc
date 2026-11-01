#include "ensim.hh"

#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <atomic>
#include <mutex>
#include <numeric>
#include <sstream>
#include <vector>

#include <raylib.h>
#include <raymath.h>

static constexpr int g_fps = 60;
static constexpr int g_panic_fps = g_fps / 10;
static constexpr float g_margin_p = 8.0f;
static constexpr size_t g_xres_p = 1600;
static constexpr size_t g_yres_p = 900;
static constexpr float g_xmid_p = g_xres_p / 2.0f;
static constexpr size_t g_mesh_slices = 10;
static constexpr float g_font_size = 1.0f;
static constexpr Color g_plot_signal_colors[] = { RED, BLUE, GREEN, PURPLE, BROWN, YELLOW };
static constexpr float g_grid_step_size = 0.125f;
static constexpr float g_rot90_r = 90.0f * DEG2RAD;

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

template<typename T>
using signal = std::function<const T&(size_t)>;

class tachometer: public widget
{
    static constexpr float start_r = 3.0f * std::numbers::pi_v<float> / 4.0f;
    static constexpr float sweep_r = 3.0f * std::numbers::pi_v<float> / 2.0f;
    static constexpr float inner_radius_ratio = 0.82f;
    static constexpr float outer_radius_ratio = 0.92f;
    static constexpr float needle_radius_ratio = 0.80f;
    static constexpr float tick_radius_ratio = 0.70f;
    static constexpr float tick_line_width_p = 2.0f;
    static constexpr float needle_thickness = 3.0f;
    static constexpr float gauge_thickness = 5.0f;
    sample<float> value = {};
    sample<float> max = {};
    size_t ticks = 0;

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
                tick_text.x - width_p / 2.0f,
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
class plot: public widget
{
    std::string title = {};
    signal<T> x = {};
    signal<T> y = {};
    size_t samples = {};
    size_t signals = {};
    bool clip_check = false;

public:
    plot(const std::string& title, const signal<T> x, const signal<T> y, const size_t samples, const size_t signals, const bool clip_check)
        : title(title)
        , x(x)
        , y(y)
        , samples(samples)
        , signals(signals)
        , clip_check(clip_check)
        {
        }

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

    range get_range() const
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
        const range range = get_range();
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
            const bool clipping = range.ymax >= 1.0 or range.ymin <= -1.0;
            DrawText(metrics.str().data(), x_p, y_p, g_font_size, clip_check and clipping ? RED : WHITE);
        }
    }
};

template <typename T>
concept Growable = requires(T self, const size_t i)
{
    self.clear();
    self.push_back(i);
};

template<Growable T>
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
    const size_t signals,
    const bool clip_check)
{
    const auto x = [y](const size_t signal)-> auto& {
        return lingen<T>(y(signal).size());
    };
    return std::make_unique<plot<T>>(title, x, y, samples, signals, clip_check);
}

template<typename T>
std::unique_ptr<plot<T>> make_plot(
    const std::string& title,
    const signal<T> x,
    const signal<T> y,
    const size_t samples,
    const size_t signals,
    const bool clip_check)
{
    return std::make_unique<plot<T>>(title, x, y, samples, signals, clip_check);
}

class widgets
{
    static constexpr float hud_gauge_width_p = 192.0f;
    static constexpr float sidebar_width_p = g_xres_p / 3;
    static constexpr size_t signal_samples = sidebar_width_p;
    static constexpr size_t tachometer_ticks = 10;
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
        this->engine = engine;
        regen();
    }

    void regen_hud()
    {
        std::unique_ptr<widget> engine_angular_velocity = std::make_unique<tachometer>(
            [this]()-> auto {
                return engine->get_engine_angular_velocity_r_per_s();
            },
            [this]()-> auto {
                return engine->get_limiter_angular_velocity_r_per_s();
            },
            tachometer_ticks
        );
        std::unique_ptr<widget> load_angular_velocity = std::make_unique<tachometer>(
            [this]()-> auto {
                return engine->get_load_angular_velocity_r_per_s();
            },
            [this]()-> auto {
                return engine->get_load_max_angular_velocity_r_per_s();
            },
            tachometer_ticks
        );
        std::unique_ptr<widget> gear = std::make_unique<tachometer>(
            [this]()-> auto {
                return engine->get_gear();
            },
            [this]()-> auto {
                return engine->get_max_gear();
            },
            engine->get_max_gear()
        );
        const float margin_p = 4.0f * g_margin_p;
        const float w_p = hud_gauge_width_p;
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
        const float gear_w_p = hud_gauge_width_p / 1.75f;
        gear->place(
            g_xmid_p - gear_w_p / 2.0f,
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
                    [this, i](auto)-> auto& {
                        return engine->get_signal(i);
                    },
                    signal_samples,
                    1,
                    false
                )
            );
        }
        float y_p = 0.0f;
        const float x_p = g_xres_p - sidebar_width_p;
        const float h_p = g_yres_p / size;
        for(auto& widget: right)
        {
            widget->place(x_p, y_p, sidebar_width_p, h_p);
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
                -1, engine->get_pipe_count(), false
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
                signal_samples, 1, false
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
                signal_samples, 1, false
            )
        );
        left.push_back(
            make_plot<std::vector<float>>(
                "audio signal",
                [this](auto)-> auto& {
                    return engine->get_audio_signal();
                },
                -1, 1, true
            )
        );
        left.push_back(
            make_plot<std::vector<size_t>>(
                "producer audio queue size",
                [](auto)-> auto& {
                    return ensim::get_buffer_size_history();
                },
                -1, 1, false
            )
        );
        left.push_back(
            make_plot<std::vector<float>>(
                "room impulse response",
                [this](auto)-> auto& {
                    return engine->get_impulse_signal();
                },
                -1, 1, false
            )
        );
        float y_p = 0.0f;
        const float x_p = 0.0f;
        const float h_p = g_yres_p / left.size();
        for(auto& widget: left)
        {
            widget->place(x_p, y_p, sidebar_width_p, h_p);
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

struct figure
{
    Model model = {};

    figure() = default;
    figure(const figure&) = delete;
    figure& operator=(const figure&) = delete;
    figure& operator=(figure&&) = delete;
    figure(figure&&) = delete;

    void load(const Mesh& mesh)
    {
        UnloadModel(model);
        model = LoadModelFromMesh(mesh);
    }

    ~figure()
    {
        UnloadModel(model);
    }
};

struct shape
{
    struct figure figure = {};
    Color color = {};
    Vector3 rotation = {};
    Vector3 position = {};
    float theta_degrees = {};
};

static float calc_cylinder_height_m(const float volume_m3, const float radius_m)
{
    return volume_m3 / (std::numbers::pi_v<float> * radius_m * radius_m);
}

static float calc_circle_radius_m(const float area_m2)
{
    return std::sqrt(area_m2 / std::numbers::pi_v<float>);
}

class part
{
protected:
    float x = 0.0;
    float z = 0.0;
    size_t x_id = 0;
    size_t y_id = 0;
    bool panicking = false;
    size_t panicks = 0;
    Vector3 origin = {};

    part(const float x, const float z)
        : x(x * g_grid_step_size)
        , z(z * g_grid_step_size)
        , x_id(x)
        , y_id(z)
        , origin(this->x, 0, this->z)
        {
        }

    void panic()
    {
        panicking ^= panicks ++ % g_panic_fps == 0;
    }

public:
    size_t get_x_id() const { return x_id; }
    size_t get_y_id() const { return y_id; }
    virtual std::span<const shape> get_shapes() const = 0;
    virtual void update() = 0;
    virtual ~part() = default;
};

class port: public part
{
    std::array<shape, 1> shapes = {};
    float radius_m = 0.0f;
    sample<float> get_open_ratio;

public:
    port(const float x,
         const float z,
         const float area_m2,
         const sample<float> get_open_ratio)
        : part(x, z)
        , radius_m(calc_circle_radius_m(area_m2))
        , get_open_ratio(get_open_ratio)
        {
            shapes[0].figure.load(GenMeshCylinder(radius_m, 0.0f, g_mesh_slices));
            shapes[0].figure.model.transform = MatrixRotateX(g_rot90_r);
        }

    void update() override
    {
        Vector3 position = origin;
        position.z += g_grid_step_size / 2.0f;
        position.y += radius_m;
        shapes[0].position = position;
        shapes[0].color = ColorLerp(RED, GREEN, get_open_ratio());
    }

    std::span<const shape> get_shapes() const override
    {
        return shapes;
    }
};

class chamber: public part
{
    static constexpr float side_m = g_grid_step_size / 2.0f;
    std::array<shape, 1> shapes = {};
    sample<bool> get_panic = {};

public:
    chamber(
        const float x,
        const float z,
        const float volume_m3,
        const sample<bool> get_panic)
        : part(x, z)
        , get_panic(get_panic)
        {
            const float height_m = volume_m3 / (side_m * side_m);
            shapes[0].figure.load(GenMeshCube(side_m, height_m, side_m));
            shapes[0].figure.model.transform = MatrixTranslate(0.0f, 0.5f * height_m, 0.0f);
        }

    void update() override
    {
        if(get_panic())
        {
            panic();
        }
        shapes[0].position = origin;
        shapes[0].color = panicking ? RED : WHITE;
    }

    std::span<const shape> get_shapes() const override
    {
        return shapes;
    }
};

class disk: public part
{
protected:
    static constexpr float density_kg_per_m3 = 7600.0f;
    std::array<shape, 1> shapes = {};
    sample<float> get_theta_r = {};

public:
    disk(
        const float x,
        const float z,
        const float mass_kg,
        const float radius_m,
        const sample<float> get_theta_r)
        : part(x, z)
        , get_theta_r(get_theta_r)
        {
            const float height_m = calc_cylinder_height_m(mass_kg / density_kg_per_m3, radius_m);
            shapes[0].figure.load(GenMeshCylinder(radius_m, height_m, g_mesh_slices));
            shapes[0].figure.model.transform = MatrixRotateX(g_rot90_r) * MatrixRotateY(g_rot90_r);
        }

    void update() override
    {
        shapes[0].position = origin;
        shapes[0].rotation = { 1.0f, 0.0f, 0.0f };
        shapes[0].theta_degrees = RAD2DEG * get_theta_r();
        shapes[0].color = WHITE;
    }

    std::span<const shape> get_shapes() const override
    {
        return shapes;
    }
};

class pop_disk: public disk
{
    static constexpr float pop_by = 0.25f * g_grid_step_size;
    sample<bool> pop = {};

public:
    pop_disk(
        const float x,
        const float z,
        const float mass_kg,
        const float radius_m,
        const sample<float> get_theta_r,
        const sample<bool> pop)
        : disk(x, z, mass_kg, radius_m, get_theta_r)
        , pop(pop)
        {
        }

    void update() override
    {
        disk::update();
        Vector3 position = origin;
        if(pop())
        {
            position.x -= pop_by;
        }
        shapes[0].position = position;
    }
};

class piston: public part
{
    float head_height_m = 0.0f;
    float conrod_width_m = 0.015f;
    float conrod_depth_m = 0.007f;
    float counter_weight_height_m = 0.0f;
    float counter_weight_depth_m = 2.0f * conrod_depth_m;
    float counter_weight_width_m = 2.0f * conrod_width_m;
    std::array<shape, 3> shapes = {};
    sample<float> get_pin_y_m = {};
    sample<float> get_pin_phi_r = {};
    sample<float> get_crank_theta_r = {};
    sample<bool> get_panic = {};

public:
    piston(
        const float x,
        const float z,
        const float head_radius_m,
        const float head_height_m,
        const float conrod_height_m,
        const float crank_diameter_m,
        const sample<float> get_pin_y_m,
        const sample<float> get_pin_phi_r,
        const sample<float> get_crank_theta_r,
        const sample<bool> get_panic)
        : part(x, z)
        , head_height_m(head_height_m)
        , counter_weight_height_m(1.2f * crank_diameter_m)
        , get_pin_y_m(get_pin_y_m)
        , get_pin_phi_r(get_pin_phi_r)
        , get_crank_theta_r(get_crank_theta_r)
        , get_panic(get_panic)
        {
            shapes[0].figure.load(GenMeshCylinder(head_radius_m, head_height_m, g_mesh_slices));
            shapes[1].figure.load(GenMeshCube(conrod_width_m, conrod_height_m, conrod_depth_m));
            shapes[2].figure.load(GenMeshCube(counter_weight_width_m, counter_weight_height_m, counter_weight_depth_m));
            shapes[0].figure.model.transform = MatrixRotateY(g_rot90_r);
            shapes[1].figure.model.transform = MatrixRotateY(g_rot90_r) * MatrixTranslate(0.0f, -0.5f * conrod_height_m, 0.0f);
            shapes[2].figure.model.transform = MatrixRotateY(g_rot90_r);
        }

    void update() override
    {
        if(get_panic())
        {
            panic();
        }
        shapes[0].position = { x, get_pin_y_m() - head_height_m / 2.0f, z };
        shapes[1].position = { x, get_pin_y_m(), z };
        shapes[2].position = { x - conrod_depth_m, 0.0f, z };
        shapes[1].rotation = { 1.0f, 0.0f, 0.0f };
        shapes[2].rotation = { 1.0f, 0.0f, 0.0f };
        shapes[1].theta_degrees = RAD2DEG * get_pin_phi_r();
        shapes[2].theta_degrees = RAD2DEG * get_crank_theta_r();
        shapes[0].color = panicking ? RED : WHITE;
        shapes[1].color = panicking ? RED : WHITE;
        shapes[2].color = panicking ? RED : WHITE;
    }

    std::span<const shape> get_shapes() const override
    {
        return shapes;
    }
};

class pipe: public part
{
    static constexpr float radius_m = 0.025f;
    static constexpr size_t capacity = 1024;
    std::array<shape, capacity> shapes = {};
    float y = 0;
    size_t cells = 0;
    signal<std::vector<float>> get_pipe_pressure = {};
    float min_pressure_pa = 0.0f;
    float max_pressure_pa = 0.0f;

public:
    pipe(
        const float x,
        const float y,
        const float z,
        const float length_m,
        const size_t cells,
        const signal<std::vector<float>> get_pipe_pressure,
        const float min_pressure_pa,
        const float max_pressure_pa)
        : part(x, z)
        , y(y)
        , cells(cells)
        , get_pipe_pressure(get_pipe_pressure)
        , min_pressure_pa(min_pressure_pa)
        , max_pressure_pa(max_pressure_pa)
        {
            const float yy = y * g_grid_step_size;
            const float cell_height_m = length_m / cells;
            for(size_t i = 0; i < cells; i++)
            {
                shapes[i].figure.load(GenMeshCylinder(radius_m, cell_height_m, g_mesh_slices));
                shapes[i].figure.model.transform = MatrixRotateZ(g_rot90_r) * MatrixTranslate(i * cell_height_m, yy, 0.0f);
            }
        }

    void update() override
    {
        const size_t y_id = y;
        const std::vector<float>& pipe_pressure = get_pipe_pressure(y_id);
        if(pipe_pressure.size() == cells)
        {
            size_t cell = 0;
            for(auto& shape: shapes)
            {
                const float cell_pressure_pa = pipe_pressure[cell];
                const float ratio = (cell_pressure_pa - min_pressure_pa) / (max_pressure_pa - min_pressure_pa);
                const Color color = g_plot_signal_colors[y_id];
                const float clamp = std::clamp(ratio, -1.0f, 1.0f);
                shape.color = ColorLerp(WHITE, color, clamp);
                shape.position = origin;
                cell++;
            }
        }
    }

    std::span<const shape> get_shapes() const override
    {
        const std::span<const shape> full = shapes;
        const std::span<const shape> less = full.first(cells);
        return less;
    }
};

class room: public part
{
    static constexpr float sphere_radius_m = 0.05;
    std::array<shape, 3> shapes = {};
    std::tuple<float,float,float> mic_position;
    std::tuple<float,float,float> impulse_position;

public:
    room(
        const float x,
        const float z,
        const std::tuple<float,float,float> dimensions,
        const std::tuple<float,float,float> mic_position,
        const std::tuple<float,float,float> impulse_position)
        : part(x, z)
        , mic_position(mic_position)
        , impulse_position(impulse_position)
        {
            {
                const auto [width_m, height_m, length_m] = dimensions;
                shapes[0].figure.load(GenMeshCube(width_m, height_m, length_m));
                shapes[0].figure.model.transform = MatrixTranslate(0.5 * width_m, 0.5 * height_m, 0.5 * length_m);
            }{
                shapes[1].figure.load(GenMeshSphere(sphere_radius_m, g_mesh_slices, g_mesh_slices));
            }{
                shapes[2].figure.load(GenMeshSphere(sphere_radius_m, g_mesh_slices, g_mesh_slices));
            }
        }

    void update() override
    {
        {
            shapes[0].position = origin;
            shapes[0].color = YELLOW;
        }{
            const auto [x_m, y_m, z_m] = mic_position;
            shapes[1].position = origin + Vector3(x_m, y_m, z_m);
            shapes[1].color = RED;
        }{
            const auto [x_m, y_m, z_m] = impulse_position;
            shapes[2].position = origin + Vector3(x_m, y_m, z_m);
            shapes[2].color = BLUE;
        }
    }

    std::span<const shape> get_shapes() const override
    {
        return shapes;
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
            if(y == engine->get_source_y() or y == engine->get_sink_y())
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
                    [this, x]()-> auto {
                        return this->engine->get_piston_pin_y_m(x);
                    },
                    [this, x]()-> auto {
                        return this->engine->get_piston_pin_phi_r(x);
                    },
                    [this, x]()-> auto {
                        return this->engine->get_piston_crank_theta_r(x);
                    },
                    [this, x, y]()-> auto {
                        return this->engine->get_panic(x, y);
                    }));
            }
            else
            {
                self.push_back(std::make_unique<chamber>(
                    xx,
                    zz,
                    engine->get_chamber_volume_m3(x, y),
                    [this, x, y]()-> auto {
                        return this->engine->get_panic(x, y);
                    }));
            }
            self.push_back(std::make_unique<port>(
                xx,
                zz,
                engine->get_port_flow_area_m2(x, y),
                [this, x, y]()-> auto {
                    return this->engine->get_port_open_ratio(x, y);
                }));
        }
        const float flywheel_x = w;
        const float flywheel_z = engine->get_piston_y();
        const float load_x = flywheel_x + 1;
        const float load_z = flywheel_z;
        const float pipe_z = h;
        self.push_back(std::make_unique<pop_disk>(
            flywheel_x,
            flywheel_z,
            engine->get_flywheel_mass_kg(),
            engine->get_flywheel_radius_m(),
            [this]()-> auto {
                return this->engine->get_flywheel_theta_r();
            },
            [this]()-> auto {
                return this->engine->get_clutch_disengaged();
            }));
        self.push_back(std::make_unique<disk>(
            load_x,
            load_z,
            engine->get_load_mass_kg(),
            engine->get_load_radius_m(),
            [this]()-> auto {
                return this->engine->get_load_theta_r();
            }));
        for(size_t y = 0; y < engine->get_pipe_count(); y++)
        {
            self.push_back(std::make_unique<pipe>(
                0.0f,
                y,
                pipe_z,
                engine->get_pipe_length_m(y),
                engine->get_pipe_cells(),
                [this](auto signal)-> const std::vector<float>& {
                    return this->engine->get_pipe_pressure_signal(signal);
                },
                 80000.0f, // TODO: PER engine
                140000.0f
            ));
        }
        const float room_x = w;
        const float room_z = h + 2;
        self.push_back(std::make_unique<room>(
            0,
            room_z,
            engine->get_room_dimensions_m(),
            engine->get_mic_position(),
            engine->get_impulse_position()
            ));
    }

    void update()
    {
        for(const auto& part: self)
        {
            part->update();
        }
    }
};

struct look
{
    static constexpr float max_pitch_r = std::numbers::pi_v<float> / 2.1f;
    static constexpr float min_pitch_r = 0.0f;
    static constexpr float min_distance = 0.1f;
    float yaw_r = std::numbers::pi_v<float> / 4.0f;
    float pitch_r = atan(sqrt(0.5f));
    float x = 0.0f;
    float y = 0.0f;
    float distance = 2.0f;

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
        distance = std::max(distance, min_distance);
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
    static constexpr float grid_steps = 60.0f;
    static constexpr Vector3 render_scale = { 1.0f, 1.0f, 1.0f };
    static constexpr const char* name = "ensim5";
    ensim::engine* engine = {};
    struct look look = {};
    float throttle_ratio = 0.0f;

public:
    window()
    {
        SetConfigFlags(FLAG_VSYNC_HINT);
        InitWindow(g_xres_p, g_yres_p, name);
        SetTargetFPS(g_fps);
        HideCursor();
        SetTraceLogLevel(LOG_WARNING);
    }

    ~window()
    {
        CloseWindow();
    }

    void draw(const widgets& widgets)
    {
        for(const auto& widget: widgets)
        {
            widget->draw();
        }
    }

    void draw(const part& part)
    {
        const auto shapes = part.get_shapes();
        for(const auto& shape: shapes)
        {
            const bool selected = part.get_x_id() == look.x and part.get_y_id() == look.y;
            const Color color = selected ? BLUE : shape.color;
            DrawModelWiresEx(shape.figure.model, shape.position, shape.rotation, shape.theta_degrees, render_scale, color);
        }
    }

    void draw(const parts& parts)
    {
        BeginMode3D(look.camera);
        DrawGrid(grid_steps, g_grid_step_size);
        for(const auto& part: parts)
        {
            draw(*part.get());
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
            engine->set_throttle(throttle_ratio = 0.00f);
            engine->set_injection_off();
        }
        if(IsKeyPressed(KEY_ONE))
        {
            engine->set_throttle(throttle_ratio = 0.00f);
            engine->set_injection_on();
        }
        if(IsKeyPressed(KEY_TWO))
        {
            engine->set_throttle(throttle_ratio = 0.33f);
            engine->set_injection_on();
        }
        if(IsKeyPressed(KEY_THREE))
        {
            engine->set_throttle(throttle_ratio = 0.66f);
            engine->set_injection_on();
        }
        if(IsKeyPressed(KEY_FOUR))
        {
            engine->set_throttle(throttle_ratio = 0.99f);
            engine->set_injection_on();
        }
        if(IsKeyPressed(KEY_PERIOD))
        {
            engine->set_throttle(0.0f);
            engine->disengage_clutch();
            engine->increment_gear();
        }
        if(IsKeyReleased(KEY_PERIOD))
        {
            engine->set_throttle(throttle_ratio);
            engine->engage_clutch();
        }
        if(IsKeyPressed(KEY_COMMA))
        {
            engine->set_throttle(0.99f);
            engine->disengage_clutch();
            engine->decrement_gear();
        }
        if(IsKeyReleased(KEY_COMMA))
        {
            engine->set_throttle(throttle_ratio);
            engine->engage_clutch();
        }
        if(IsKeyPressed(KEY_SPACE))
        {
            engine->engage_load_brake();
        }
        if(IsKeyReleased(KEY_SPACE))
        {
            engine->disengage_load_brake();
        }
        if(IsKeyDown(KEY_Q)) look.distance += 0.05f;
        if(IsKeyDown(KEY_E)) look.distance -= 0.05f;
        if(IsKeyDown(KEY_H)) look.yaw_r -= 0.1f;
        if(IsKeyDown(KEY_L)) look.yaw_r += 0.1f;
        if(IsKeyDown(KEY_J)) look.pitch_r -= 0.1f;
        if(IsKeyDown(KEY_K)) look.pitch_r += 0.1f;
        int dx = 0;
        int dy = 0;
        if(IsKeyPressed(KEY_W)) dy = -1;
        if(IsKeyPressed(KEY_S)) dy = +1;
        if(IsKeyPressed(KEY_A)) dx = -1;
        if(IsKeyPressed(KEY_D)) dx = +1;
        const int dirs = 4;
        const int dir = std::round(look.yaw_r / (std::numbers::pi_v<float> / 2.0f));
        const int id = ((dir % dirs) + dirs) % dirs;
        int mx = 0;
        int my = 0;
        switch(id)
        {
        case 0: mx =  dx; my = +dy; break;
        case 1: mx =  dy; my = -dx; break;
        case 2: mx = -dx; my = -dy; break;
        case 3: mx = -dy; my = +dx; break;
        }
        look.x += mx;
        look.y += my;
    }

    void handle_look()
    {
        look.update();
        engine->set_logger(look.x, look.y);
    }

    void handle_part_updates(parts& parts)
    {
        parts.update();
    }

    void loop(const widgets& widgets, parts& parts)
    {
        while(not done())
        {
            handle_input();
            handle_look();
            handle_part_updates(parts);
            draw(widgets, parts);
        }
    }

    void set(ensim::engine* other)
    {
        engine = other;
        look.x = engine->get_width() - 1;;
        look.y = engine->get_piston_y();
    }

    bool done()
    {
        return WindowShouldClose();
    }
};

class audio
{
    static constexpr int bit_depth = 32;
    static constexpr int channels = 1;
    static constexpr size_t buffer_size = 4096;
    AudioStream stream = {};

    static void callback(void* const data, const unsigned frames)
    {
        float* const buffer = static_cast<float*>(data);
        const std::vector<float>& block = ensim::consume(frames);
        std::copy(block.begin(), block.end(), buffer);
    };

public:
    audio()
    {
        InitAudioDevice();
        SetAudioStreamBufferSizeDefault(buffer_size);
        stream = LoadAudioStream(ensim::g_sample_rate_hz, bit_depth, channels);
        SetAudioStreamCallback(stream, callback);
        PlayAudioStream(stream);
    }

    ~audio()
    {
        UnloadAudioStream(stream);
        CloseAudioDevice();
    }
};

static void set(window& window, widgets& widgets, parts& parts, ensim::engine* engine)
{
    window.set(engine);
    widgets.set(engine);
    ensim::attach(engine);
    parts.set(engine);
}

int main()
{
    std::unique_ptr<ensim::engine> engine = ensim::new_engine(ensim::type::i8_2400cc);
    window window;
    widgets widgets;
    audio audio;
    parts parts;
    set(window, widgets, parts, engine.get());
    window.loop(widgets, parts);
    ensim::kill();
}

// TODO
// * add semicircle for counter balance mass
// * add flywheel
