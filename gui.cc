#include "ensim.hh"

#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <numeric>
#include <raylib.h>
#include <raymath.h>
#include <sstream>
#include <thread>
#include <vector>

static constexpr float g_margin_p = 8.0f;
static constexpr int g_xres_p = 1920;
static constexpr int g_yres_p = 1080;
static constexpr int g_fps = 60;
static constexpr const char* g_name = "ensim5";
static constexpr float g_sidebar_width_p = g_xres_p / 5;
static constexpr size_t g_producer_block_size = 200;
static constexpr int g_audio_bit_depth = 32;
static constexpr int g_audio_channels = 1;
static constexpr size_t g_mesh_slices = 8;
static constexpr float g_font_size = 1.0f;
static constexpr float g_plot_border_line_thickness = 1.0f;
static constexpr Vector3 g_render_scale = { 1.0f, 1.0f, 1.0f };
static constexpr int g_stream_buffer_size = 8192;
static constexpr Color g_plot_border_color = GRAY;
static constexpr Color g_plot_signal_colors[] = { RED, BLUE, GREEN, PURPLE, BROWN, YELLOW };
static constexpr size_t g_signal_stride = 8;
static constexpr float g_grid_steps = 32.0f;
static constexpr float g_grid_step_size = 0.15f;

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
    ring<float, 2048> queue = {};
    ring<size_t, 1024> history = {};
    ring<size_t, 1024> snapshot = {};
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
                const auto delay = std::chrono::microseconds(100);
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
using signal = std::function<const T&(size_t)>;

template<typename T>
class plot : public widget
{
    std::string title = {};
    signal<T> x = {};
    signal<T> y = {};
    size_t stride = {};
    size_t signals = {};

public:
    plot(const std::string& title, const signal<T> x, const signal<T> y, const size_t stride, const size_t signals)
        : title(title)
        , x(x)
        , y(y)
        , stride(stride)
        , signals(signals) {}

    struct range
    {
        double xmin = std::numeric_limits<double>::max();
        double xmax = std::numeric_limits<double>::lowest();
        double ymin = std::numeric_limits<double>::max();
        double ymax = std::numeric_limits<double>::lowest();

        bool valid() const
        {
            return xmin < xmax and ymin < ymax;
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
            size_t samples = 0;
            for(size_t signal = 0; signal < signals; signal++)
            {
                const T& sx = x(signal);
                const T& sy = y(signal);
                const size_t size = sx.size();
                samples = size / stride;
                points.resize(samples);
                size_t at = 0;
                for(size_t i = 0; i < size and at < samples; i += stride)
                {
                    const double nx_p = xrange > 0.0 ? (sx[i] - range.xmin) / xrange : 0.5;
                    const double ny_p = yrange > 0.0 ? (sy[i] - range.ymin) / yrange : 0.5;
                    const double w_p = -2.0 * g_margin_p + rectangle.width;
                    const double h_p = -2.0 * g_margin_p + rectangle.height;
                    points[at].x = x_p + (0.0 + nx_p) * w_p;
                    points[at].y = y_p + (1.0 - ny_p) * h_p;
                    at++;
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
                << "\nsamples: " << samples;
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
    const size_t stride,
    const size_t signals = 1)
{
    const auto x = [y](const size_t signal)-> auto& {
        return lingen<T>(y(signal).size());
    };
    return std::make_unique<plot<T>>(title, x, y, stride, signals);
}

template<typename T>
std::unique_ptr<plot<T>> make_plot(
    const std::string& title,
    const signal<T> x,
    const signal<T> y,
    const size_t stride,
    const size_t signals = 1)
{
    return std::make_unique<plot<T>>(title, x, y, stride, signals);
}

class sidebar
{
    std::vector<std::unique_ptr<widget>> widgets = {};
    ensim::engine* engine = {};

public:
    auto begin()
    {
        return widgets.begin();
    }

    auto end()
    {
        return widgets.end();
    }

    void set(ensim::engine* engine)
    {
        this->engine = engine;
        regen();
    }

    void regen_right()
    {
        std::vector<std::unique_ptr<widget>> right;
        const size_t size = engine->get_signal_count();
        for(size_t i = 0; i < size; i++)
        {
            right.push_back(
                make_plot<std::vector<double>>(
                    std::string(engine->get_signal_name(i)),
                    [this, i](auto)-> auto& { return engine->get_signal(i); },
                    g_signal_stride
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
            widgets.push_back(std::move(widget));
        }
    }

    void regen_left()
    {
        std::vector<std::unique_ptr<widget>> left;
        left.push_back(
            make_plot<std::vector<float>>(
                "pipe pressure (pascals)",
                [this](auto signal)-> auto& {
                    return engine->get_pipe_pressure_signal(signal);
                },
                1,
                engine->get_pipe_count()
            )
        );
        left.push_back(
            make_plot<ring<size_t, 1024>>(
                "producer audio queue size",
                [](auto)-> auto& {
                    return g_producer.get_history();
                },
                1
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
                g_signal_stride
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
                g_signal_stride
            )
        );
        float y_p = 0.0f;
        const float x_p = 0.0f;
        const float h_p = g_yres_p / left.size();
        for(auto& widget : left)
        {
            widget->place(x_p, y_p, g_sidebar_width_p, h_p);
            y_p += h_p;
            widgets.push_back(std::move(widget));
        }
    }

    void regen()
    {
        widgets.clear();
        regen_left();
        regen_right();
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
    part(const float x, const float z): x(x * g_grid_step_size), z(z * g_grid_step_size) {}

public:
    virtual std::span<const shape> get_shapes() const = 0;
    virtual void update() = 0;
    virtual ~part() = default;
};

class cylinder : public part
{
    float radius_m = g_grid_step_size / 4.0f;
    float height_m = 0.0f;
    std::array<shape, 1> shapes = {};

public:
    cylinder(
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
    float head_radius_m = 0.0f;
    float head_height_m = 0.0f;
    float conrod_height_m = 0.0f;
    float conrod_width_m = 0.015f;
    float conrod_depth_m = 0.007f;
    float counter_weight_height_m = 0.0f;
    float counter_weight_depth_m = 1.75f * conrod_depth_m;
    float counter_weight_width_m = 1.75f * conrod_width_m;
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
        , head_radius_m(head_radius_m)
        , head_height_m(head_height_m)
        , conrod_height_m(conrod_height_m)
        , counter_weight_height_m(1.5f * crank_diameter_m)
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

        const float theta_r = 90.0f * DEG2RAD;
        shapes[0].model.transform = MatrixRotateY(theta_r);
        shapes[1].model.transform = MatrixMultiply(MatrixTranslate(0.0f, -conrod_height_m * 0.5f, 0.0f), MatrixRotateY(theta_r));
        shapes[2].model.transform = MatrixRotateY(theta_r);
    }

    void update() override
    {
        shapes[0].position = { x, get_pin_y_m() - head_height_m / 2.0f, z };
        shapes[1].position = { x, get_pin_y_m(), z };
        shapes[2].position = { x - conrod_depth_m, 0.0f, z };
        shapes[1].rotation = { 1.0f, 0.0f, 0.0f };
        shapes[2].rotation = { 1.0f, 0.0f, 0.0f };
        shapes[1].theta_degrees = +RAD2DEG * get_pin_phi_r();
        shapes[2].theta_degrees = -RAD2DEG * get_crank_theta_r();
    }

    std::span<const shape> get_shapes() const override
    {
        return shapes;
    }

    ~piston()
    {
        for(auto& shape : shapes)
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

    auto end()
    {
        return self.end();
    }

    void set(ensim::engine* engine)
    {
        self.clear();
        this->engine = engine;
        for(size_t x = 0; x < engine->get_width(); x++)
        for(size_t y = 0; y < engine->get_height(); y++)
        {
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
                push(std::make_unique<piston>(
                    x,
                    y,
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
                push(std::make_unique<cylinder>(x, y, engine->get_chamber_volume_m3(x, y)));
            }
        }
    }

    void push(std::unique_ptr<part> part)
    {
        self.push_back(std::move(part));
    }

    void update()
    {
        for(auto& part : self)
        {
            part->update();
        }
    }
};

class window
{
    ensim::engine* engine = {};
    Camera3D camera = {};

public:
    window()
    {
        SetConfigFlags(FLAG_FULLSCREEN_MODE);
        InitWindow(g_xres_p, g_yres_p, g_name);
        SetTargetFPS(g_fps);
    }

    ~window()
    {
        CloseWindow();
    }

    void draw(sidebar& sidebar)
    {
        for(auto& widget : sidebar)
        {
            widget->draw();
        }
    }

    void draw(parts& parts)
    {
        BeginMode3D(camera);
        DrawGrid(16, g_grid_step_size);
        for(auto& part : parts)
        {
            const std::span<const shape> shapes = part->get_shapes();
            for(const auto& shape : shapes)
            {
                DrawModelWiresEx(shape.model, shape.position, shape.rotation, shape.theta_degrees, g_render_scale, shape.color);
            }
        }
        EndMode3D();
    }

    void loop(sidebar& sidebar, parts& parts)
    {
        while(not done())
        {
            parts.update();
            open();
            draw(parts);
            draw(sidebar);
            close();
        }
    }

    void set(ensim::engine* engine)
    {
        this->engine = engine;
        camera = {
            .position = { 1.0f, 0.0f, 0.0f },
            .target = { 0.0f, 0.0f, 0.0f },
            .up = { 0.0f, 1.0f, 0.0f },
            .fovy = 45.0f,
            .projection = CAMERA_PERSPECTIVE,
        };
    }

    void open()
    {
        BeginDrawing();
        ClearBackground(BLACK);
        UpdateCamera(&camera, CAMERA_THIRD_PERSON);
        engine->set_swap_lock_on();
    }

    void close()
    {
        engine->set_swap_lock_off();
        EndDrawing();
    }

    bool done()
    {
        return WindowShouldClose();
    }
};

static void set(window& window, sidebar& sidebar, producer& producer, parts& parts, ensim::engine* engine)
{
    window.set(engine);
    sidebar.set(engine);
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
    engine->set_logger(0, 4);
    window window;
    sidebar sidebar;
    audio audio(audio_callback);
    parts parts;
    set(window, sidebar, g_producer, parts, engine.get());
    window.loop(sidebar, parts);
    g_producer.stop();
}
