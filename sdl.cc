#include "ensim.hh"

#include <SDL3/SDL.h>
#include <algorithm>
#include <numbers>
#include <deque>
#include <thread>
#include <chrono>
#include <cmath>
#include <list>

static constexpr size_t g_signals_count = 9;

struct point
{
    SDL_FPoint self;
    uint32_t color;

    point(const int x_p, const int y_p, const uint32_t color)
        : self(x_p, y_p)
        , color(color)
        {
        }
};

struct points
{
    std::vector<SDL_FPoint> self;
    uint32_t color = 0x0;

    points(const uint32_t color):
        color(color)
        {
        }

    void append(const SDL_FPoint& point)
    {
        self.push_back(point);
    }

    void append(const double x, const double y)
    {
        const float xx = x;
        const float yy = y;
        const SDL_FPoint point = { xx, yy };
        append(point);
    }
};

struct rect
{
    SDL_FRect self;
    uint32_t color;
    double color_ratio;

    rect(const int x_p, const int y_p, const int w_p, const int h_p, const uint32_t color, const double color_ratio = 1.0)
        : self(x_p, y_p, w_p, h_p)
        , color(color)
        , color_ratio(color_ratio)
        {
        }

    rect(const rect& other, const uint32_t color)
    {
        *this = other;
        this->color = color;
    }

    SDL_FPoint project(const double x_ratio, const double y_ratio) const
    {
        SDL_FPoint out;
        out.x = self.x + x_ratio * self.w;
        out.y = self.y + (1.0 - y_ratio) * self.h;
        return out;
    }
};

struct circle
{
    SDL_FPoint self;
    double radius;
    uint32_t color;

    circle(const int x_p, const int y_p, const uint32_t color)
        : self(x_p, y_p)
        , color(color)
        {
        }

    circle(const rect& rect, const uint32_t color, const double border_ratio = 1.0)
    {
        self.x = rect.self.x + rect.self.w * 0.5;
        self.y = rect.self.y + rect.self.h * 0.5;
        radius = std::min(rect.self.w, rect.self.h) * 0.5 * border_ratio;
        this->color = color;
    }
};

struct message: point
{
    std::string string;

    message(const int x_p, const int y_p, const uint32_t color, const std::string& string)
        : point(x_p, y_p, color)
        , string(string)
        {
        }
};

struct sdl
{
    static constexpr int w_p = 1920;
    static constexpr int h_p = 1080;
    static constexpr uint32_t font_p = 8;
    static constexpr uint32_t line_p = 1.5 * font_p;
    static constexpr uint32_t grey   = 0xFFAAAAAA;
    static constexpr uint32_t white  = 0xFFFFFFFF;
    static constexpr uint32_t black  = 0xFF101010;
    static constexpr uint32_t green  = 0xFF00FF00;
    static constexpr uint32_t purple = 0xFFFF00FF;
    static constexpr uint32_t blue   = 0xFF00AAFF;
    static constexpr uint32_t orange = 0xFFFFAA00;
    static constexpr uint32_t red    = 0xFFFF2222;
    static constexpr uint32_t yellow = 0xFFFFFF00;

    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_AudioStream* audio_stream;
    SDL_AudioSpec audio_spec;

    sdl()
    {
        SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO);
        SDL_CreateWindowAndRenderer("ensim5", w_p, h_p, SDL_WINDOW_BORDERLESS, &window, &renderer);
        SDL_SetRenderVSync(renderer, true);
        audio_spec.channels = 1;
        audio_spec.format = SDL_AUDIO_F32;
        audio_spec.freq = ensim::g_sample_rate_hz;
        audio_stream = SDL_OpenAudioDeviceStream(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, &audio_spec, nullptr, nullptr);
        SDL_ResumeAudioStreamDevice(audio_stream);
    }

    ~sdl()
    {
        SDL_DestroyAudioStream(audio_stream);
        SDL_DestroyWindow(window);
        SDL_DestroyRenderer(renderer);
    }

    std::vector<SDL_Event> poll()
    {
        std::vector<SDL_Event> events;
        SDL_Event event;
        while(SDL_PollEvent(&event))
        {
            events.push_back(event);
        }
        return events;
    }

    void set_color(const uint32_t hex, const double ratio = 1.0)
    {
        const uint8_t a = uint8_t(hex >> 24) * ratio;
        const uint8_t r = uint8_t(hex >> 16) * ratio;
        const uint8_t g = uint8_t(hex >>  8) * ratio;
        const uint8_t b = uint8_t(hex >>  0) * ratio;
        SDL_SetRenderDrawColor(renderer, r, g, b, a);
    }

    void render()
    {
        SDL_RenderPresent(renderer);
    }

    void clear()
    {
        set_color(black);
        SDL_RenderClear(renderer);
    }

    void delay(const int ms)
    {
        SDL_Delay(ms);
    }

    void outline(const rect& rect)
    {
        set_color(rect.color);
        SDL_RenderRect(renderer, &rect.self);
    }

    void fill(const rect& rect)
    {
        set_color(rect.color, rect.color_ratio);
        SDL_RenderFillRect(renderer, &rect.self);
        set_color(white);
        SDL_RenderRect(renderer, &rect.self);
    }

    void write(const message& message, const bool center = false)
    {
        set_color(message.color);
        double x_p = message.self.x;
        double y_p = message.self.y;
        if(center)
        {
            x_p -= sdl::font_p / 2 * message.string.size();
            y_p -= sdl::font_p / 2;
        }
        SDL_RenderDebugText(renderer, x_p, y_p, message.string.data());
    }

    void write(const point& point, const std::vector<std::string>& strings)
    {
        size_t i = 0;
        for(const auto& x : strings)
        {
            const message message(point.self.x, point.self.y + i * sdl::line_p, point.color, x);
            write(message);
            i++;
        }
    }

    void draw_line(const point& from, const point& to)
    {
        set_color(from.color);
        SDL_RenderLine(renderer, from.self.x, from.self.y, to.self.x, to.self.y);
    }

    void draw_lines(const points& points)
    {
        const size_t size = points.self.size();
        if(size > 0)
        {
            set_color(points.color);
            SDL_RenderLines(renderer, points.self.data(), size);
        }
    }

    void draw_points(const points& points)
    {
        const size_t size = points.self.size();
        if(size > 0)
        {
            set_color(points.color);
            SDL_RenderPoints(renderer, points.self.data(), size);
        }
    }

    void draw_circle(const circle& circle)
    {
        points points(circle.color);
        int x = 0;
        int y = circle.radius;
        int d = 3 - 2 * circle.radius;
        while(x <= y)
        {
            const SDL_FPoint steps[] = {
                { circle.self.x + x, circle.self.y + y },
                { circle.self.x - x, circle.self.y + y },
                { circle.self.x + x, circle.self.y - y },
                { circle.self.x - x, circle.self.y - y },
                { circle.self.x + y, circle.self.y + x },
                { circle.self.x - y, circle.self.y + x },
                { circle.self.x + y, circle.self.y - x },
                { circle.self.x - y, circle.self.y - x },
            };
            for(const auto& step : steps)
            {
                points.append(step);
            }
            if(d < 0)
            {
                d += 4 * x + 6;
            }
            else
            {
                d += 4 * (x - y) + 10;
                y--;
            }
            x++;
        }
        draw_points(points);
    }

    int get_audio_buffer_size()
    {
        return SDL_GetAudioStreamQueued(audio_stream) / sizeof(double);
    }

    void buffer_audio(const std::vector<float>& audio)
    {
        SDL_PutAudioStreamData(audio_stream, audio.data(), audio.size() * sizeof(float));
    }
};

auto minmax(const auto& line)
{
    const auto [x, y] = std::minmax_element(line.begin(), line.end());
    const auto min = (x == line.end()) ? 0.0 : *x;
    const auto max = (y == line.end()) ? 0.0 : *y;
    return std::pair { min, max };
}

std::vector<float> lingen(const size_t size)
{
    std::vector<float > linear;
    for(size_t i = 0; i < size; i++)
    {
        linear.push_back(i);
    }
    return linear;
}

std::vector<float> normalize(const std::vector<float>& line)
{
    const auto [min, max] = minmax(line);
    std::vector<float> out = line;
    for(auto& x : out)
    {
        x /= max;
    }
    return out;
}

std::vector<float> downsample(const auto& line, const size_t size)
{
    std::vector<float> out;
    if(line.empty())
    {
        return out;
    }
    out.resize(size);
    for(size_t i = 0; i < size; i++)
    {
        const size_t j = i * line.size() / size;
        out[i] = line[j];
    }
    return out;
}

points project(const auto& xx, const auto& yy, const rect& rect)
{
    const auto [x_min, x_max] = minmax(xx);
    const auto [y_min, y_max] = minmax(yy);
    const size_t size = xx.size();
    points points(rect.color);
    points.self.reserve(size);
    for(size_t i = 0; i < size; i++)
    {
        const double dx = x_max - x_min;
        const double dy = y_max - y_min;
        const double x_ratio = (xx[i] - x_min) / dx;
        const double y_ratio = (yy[i] - y_min) / dy;
        const SDL_FPoint point = rect.project(x_ratio, y_ratio);
        points.append(point);
    }
    return points;
}

points project_1d(const auto& y, const rect& rect, const size_t size)
{
    const std::vector<float> temp = downsample(y, size);
    const std::vector<float> xx = lingen(temp.size());
    const std::vector<float> yy = normalize(temp);
    return project(xx, yy, rect);
}

points project_2d(const auto& x, const auto& y, const rect& rect, const size_t size)
{
    const std::vector<float> xx = normalize(downsample(x, size));
    const std::vector<float> yy = normalize(downsample(y, size));
    return project(xx, yy, rect);
}

struct cell
{
    virtual void draw(sdl& sdl) = 0;
    virtual ~cell() = default;
};

struct chamber: cell
{
    static constexpr int w_p = sdl::w_p / 32;
    static constexpr int h_p = sdl::h_p / g_signals_count;
    static constexpr int ws_p = w_p / 3;
    static constexpr int hs_p = w_p / 3;
    static constexpr uint32_t fill_color = sdl::black;
    static constexpr uint32_t throttle_color = sdl::orange;
    static constexpr uint32_t piston_color = sdl::purple;
    static constexpr uint32_t audio_color = sdl::blue;
    static constexpr uint32_t panic_color = sdl::red;
    static constexpr uint32_t select_color = sdl::yellow;
    int x;
    int y;
    int& x_select;
    int& y_select;
    ensim::engine& engine;

    chamber(const int x, const int y, int& x_select, int& y_select, ensim::engine& engine)
        : x(x)
        , y(y)
        , x_select(x_select)
        , y_select(y_select)
        , engine(engine)
        {
        }

    void draw(sdl& sdl) override
    {
        const int x_p = x * w_p;
        const int y_p = y * h_p;
        const rect container(x_p, y_p, w_p, h_p, fill_color);
        const int xs_p = x_p + (w_p - ws_p) / 2;
        const int ys_p = y_p + (h_p - hs_p) / 2;
        const int xf_p = x_p + w_p / 2;
        const int yf_p = y_p + h_p / 2;
        const int y_throttle = engine.get_throttle_y();
        const int y_piston = engine.get_piston_y();
        const int y_audio = engine.get_audio_y();
        sdl.fill(container);
        if(x == x_select and y == y_select)
        {
            sdl.fill(rect(xs_p, ys_p, ws_p, hs_p, select_color));
        }
        if(engine.get_panic(x, y))
        {
            sdl.fill(rect(xs_p, ys_p, ws_p, hs_p, panic_color));
        }
        if(y == y_throttle)
        {
            sdl.write(message(xf_p, yf_p, throttle_color, "T"), true);
        }
        if(y == y_piston)
        {
            sdl.write(message(xf_p, yf_p, piston_color, "P"), true);
        }
        if(y == y_audio)
        {
            sdl.write(message(xf_p, yf_p, audio_color, "A"), true);
        }
    }
};

struct frame: cell
{
    static constexpr uint32_t frame_color = sdl::white;

    void draw(sdl& sdl) override
    {
        const int x_p = 0;
        const int y_p = 0;
        const rect rect(x_p, y_p, sdl::w_p, sdl::h_p, frame_color);
        sdl.outline(rect);
    }
};

struct port: cell
{
    static constexpr int w_p = chamber::w_p / 4;
    static constexpr int h_p = chamber::w_p / 4;
    static constexpr uint32_t fill_color = sdl::green;

    int x;
    int y;
    ensim::engine& engine;

    port(const int x, const int y, ensim::engine& engine)
        : x(x)
        , y(y)
        , engine(engine)
        {
        }

    void draw(sdl& sdl) override
    {
        const int x_p = chamber::w_p * x + chamber::w_p / 2 - w_p / 2;
        const int y_p = chamber::h_p * y + chamber::h_p / 1 - h_p;
        const std::atomic<double>& ratio = engine.get_port_open_ratio(x, y);
        const rect rect(x_p, y_p, w_p, h_p, fill_color, ratio);
        sdl.fill(rect);
    }
};

struct plot: cell
{
    static constexpr size_t max_points = sdl::w_p;
    static constexpr int h_p = sdl::h_p / g_signals_count;
    static constexpr uint32_t signal_color = sdl::red;
    static constexpr uint32_t fill_color = sdl::black;
    static constexpr uint32_t text_color = sdl::white;
    static constexpr uint32_t zero_line_color = sdl::grey;
    int y;
    ensim::engine& engine;

    plot(const int y, ensim::engine& engine)
        : y(y)
        , engine(engine)
        {
        }

    void draw(sdl& sdl) override
    {
        const std::string_view name = engine.get_signal_name(y);
        const std::vector<double>& y_signal = engine.get_signal(y);
        const int x_p = engine.get_width() * chamber::w_p;
        const int y_p = y * h_p;
        const rect fill(x_p, y_p, sdl::w_p - x_p, h_p, fill_color);
        const rect signal(fill, signal_color);
        const points data = project_1d(y_signal, signal, max_points);
        const auto [min, max] = minmax(y_signal);
        const point font(
            x_p + sdl::line_p,
            y_p + sdl::line_p,
            text_color
        );
        const std::vector<std::string> strings = {
            std::string(name),
            "max " + std::to_string(max),
            "min " + std::to_string(min),
            "div " + (min ? std::to_string(max / min) : std::string("N/A")),
            "rng " + std::to_string(max - min),
        };
        const int yz_p = y_p + h_p * (max / (max - min));
        const point yz0_p(x_p, yz_p, zero_line_color);
        const point yz1_p(sdl::w_p, yz_p, zero_line_color);
        sdl.fill(fill);
        sdl.draw_line(yz0_p, yz1_p);
        sdl.draw_lines(data);
        sdl.write(font, strings);
    }
};

struct popup: cell
{
    static constexpr size_t max_points = 1024;
    static constexpr int w_p = sdl::h_p / 2;
    static constexpr int h_p = sdl::h_p / 2;
    static constexpr int dw_p = 2 * sdl::line_p;
    static constexpr int dh_p = 2 * sdl::line_p;
    int index;

    popup(const int index)
        : index(index)
        {
        }

    std::pair<int, int> calc_position() const
    {
        const int x_p = sdl::w_p - w_p - index * dw_p;
        const int y_p = index * dh_p;
        return { x_p, y_p };
    }

    rect calc_rect() const
    {
        const auto [x_p, y_p] = calc_position();
        return rect(x_p, y_p, w_p, h_p, sdl::black);
    }
};

struct plot_popup: popup
{
    static constexpr uint32_t signal_color = sdl::red;
    static constexpr uint32_t text_color = sdl::white;
    std::string name;
    const std::vector<double>& x_signal;
    const std::vector<double>& y_signal;

    plot_popup(const int index, const std::string& name, const std::vector<double>& x_signal, const std::vector<double>& y_signal)
        : popup(index)
        , name(name)
        , x_signal(x_signal)
        , y_signal(y_signal)
        {
        }

    void draw(sdl& sdl) override
    {
        const rect fill = calc_rect();
        const rect signal(fill, signal_color);
        sdl.fill(fill);
        const points data = project_2d(x_signal, y_signal, signal, max_points);
        const point font(fill.self.x + sdl::line_p, fill.self.y + sdl::line_p, text_color);
        const auto [x_min, x_max] = minmax(x_signal);
        const auto [y_min, y_max] = minmax(y_signal);
        const std::vector<std::string> strings = {
            name,
            "x min = " + std::to_string(x_min),
            "x max = " + std::to_string(x_max),
            "y max = " + std::to_string(y_max),
            "y min = " + std::to_string(y_min),
        };
        sdl.draw_lines(data);
        sdl.write(font, strings);
    }
};

struct audio_popup: popup
{
    static constexpr uint32_t signal_color = sdl::red;
    static constexpr uint32_t text_color = sdl::white;
    const std::string name;
    const std::vector<float>& audio_signal;

    audio_popup(const int index, const std::string& name, const std::vector<float>& audio_signal)
        : popup(index)
        , name(name)
        , audio_signal(audio_signal)
        {
        }

    void draw(sdl& sdl) override
    {
        const rect fill = calc_rect();
        const rect signal(fill, signal_color);
        const points data = project_1d(audio_signal, signal, max_points);
        const point font(fill.self.x + sdl::line_p, fill.self.y + sdl::line_p, text_color);
        const auto [y_min, y_max] = minmax(audio_signal);
        const double amplitude = y_max - y_min;
        const bool clipping = y_min <= -1.0 || y_max >= 1.0;
        const std::vector<std::string> strings = {
            name,
            "max = " + std::to_string(y_max),
            "min = " + std::to_string(y_min),
            "amplitude = " + std::to_string(amplitude),
            "samples = " + std::to_string(audio_signal.size()),
            std::string(clipping ? "CLIPPING!" : ""),
        };
        sdl.fill(fill);
        sdl.draw_lines(data);
        sdl.write(font, strings);
    }
};

struct pipe_popup: popup
{
    static constexpr uint32_t signal_colors[] = {
        sdl::green,
        sdl::purple,
        sdl::blue,
        sdl::orange,
        sdl::red,
        sdl::yellow,
    };
    const std::string name;
    const ensim::engine& engine;

    pipe_popup(const int index, const std::string& name, const ensim::engine& engine)
        : popup(index)
        , name(name)
        , engine(engine)
        {
        }

    void draw(sdl& sdl) override
    {
        const rect fill = calc_rect();
        sdl.fill(fill);
        const size_t pipes = engine.get_pipe_count();
        for(size_t i = 0; i < pipes; i++)
        {
            const std::span<const float> pipe_signal = engine.get_pipe_pressure_signal(i);
            const rect signal(fill, signal_colors[i]);
            const points data = project_1d(pipe_signal, signal, max_points);
            const auto [y_min, y_max] = minmax(pipe_signal);
            const double amplitude = y_max - y_min;
            sdl.draw_lines(data);
            if(i == 0)
            {
                const point font(fill.self.x + sdl::line_p, fill.self.y + sdl::line_p, signal_colors[i]);
                const std::vector<std::string> strings = {
                    name,
                    "max = " + std::to_string(y_max),
                    "min = " + std::to_string(y_min),
                    "amplitude = " + std::to_string(amplitude),
                    "samples = " + std::to_string(pipe_signal.size()),
                };
                sdl.write(font, strings);
            }
        }
    }
};

struct gauge_popup: popup
{
    static constexpr double tick_divisor = 50.0;
    static constexpr double start_theta_r = (4.0 / 3.0) * std::numbers::pi_v<double>;
    static constexpr double sweep_theta_r = (5.0 / 3.0) * std::numbers::pi_v<double>;
    static constexpr double outer_ratio = 0.80;
    static constexpr double inner_ratio = 0.05;
    static constexpr double ticks_ratio = 0.82;
    static constexpr double needle_ratio = 0.85;
    static constexpr uint32_t needle1_color = sdl::red;
    static constexpr uint32_t needle2_color = sdl::yellow;
    static constexpr uint32_t inner_color = sdl::grey;
    static constexpr uint32_t outer_color = sdl::grey;
    static constexpr uint32_t ticks_color = sdl::white;
    static constexpr uint32_t gear_color = sdl::green;
    static constexpr uint32_t text_color = sdl::white;
    const std::string name;
    ensim::engine& engine;
    const std::atomic<double>& w1;
    const std::atomic<double>& w2;
    const std::atomic<size_t>& gear;
    const std::atomic<double>& max;

    gauge_popup(const int index, const std::string& name, ensim::engine& engine)
        : popup(index)
        , name(name)
        , engine(engine)
        , w1(engine.get_engine_angular_velocity_r_per_s())
        , w2(engine.get_load_angular_velocity_r_per_s())
        , gear(engine.get_gear())
        , max(engine.get_limiter_angular_velocity_r_per_s())
        {
        }

    void draw(sdl& sdl) override
    {
        const rect rect = calc_rect();
        const point text(
            rect.self.x + sdl::line_p,
            rect.self.y + sdl::line_p,
            text_color
        );
        const std::vector<std::string> strings = {
            name,
            std::to_string(w1),
            std::to_string(w2),
        };
        const circle outer(rect, outer_color, outer_ratio);
        const circle inner(rect, inner_color, inner_ratio);
        const message above(inner.self.x, inner.self.y - outer.radius / 4, gear_color, std::to_string(gear));
        sdl.fill(rect);
        sdl.draw_circle(outer);
        sdl.draw_circle(inner);
        sdl.write(above, true);
        draw_ticks(sdl, outer);
        draw_needle(sdl, outer, w1, needle1_color);
        draw_needle(sdl, outer, w2, needle2_color);
        sdl.write(text, strings);
    }

    double to_angle(const double at) const
    {
        return start_theta_r - (at / max) * sweep_theta_r;
    }

    void draw_needle(sdl& sdl, const circle& outer, const double at, const uint32_t color) const
    {
        const point middle(outer.self.x, outer.self.y, color);
        const double angle_r = to_angle(at);
        const double radius = needle_ratio * outer.radius;
        const point tip(
            middle.self.x + std::cos(angle_r) * radius,
            middle.self.y - std::sin(angle_r) * radius,
            color
        );
        sdl.draw_line(middle, tip);
    }

    void draw_ticks(sdl& sdl, const circle& outer) const
    {
        const size_t needle_ticks = max / tick_divisor;
        const double step = max / needle_ticks;
        const double radius = outer.radius * ticks_ratio;
        for(size_t i = 0; i <= needle_ticks; i++)
        {
            const double tick = i * step;
            const double angle_r = to_angle(tick);
            const message message(
                outer.self.x + std::cos(angle_r) * radius,
                outer.self.y - std::sin(angle_r) * radius,
                ticks_color,
                std::to_string(static_cast<int>(tick))
            );
            sdl.write(message, true);
        }
    }
};

struct help_popup: popup
{
    static constexpr uint32_t text_color = sdl::white;
    ensim::engine& engine;

    help_popup(const int index, ensim::engine& engine)
        : popup(index)
        , engine(engine)
        {
        }

    void draw(sdl& sdl) override
    {
        const rect rect = calc_rect();
        const point point(
            rect.self.x + sdl::line_p,
            rect.self.y + sdl::line_p,
            text_color
        );
        const std::vector<std::string> strings = {
            "                   _           ______  ",
            "  ___  ____  _____(_)___ ___  / ____/  ",
            " / _ \\/ __ \\/ ___/ / __ `__ \\/___ \\",
            "/  __/ / / (__  ) / / / / / /___/ /    ",
            "\\___/_/ /_/____/_/_/ /_/ /_/_____/    ",
            "",
            "1,2,3,4 for throttle",
            "Q,E to cycle through these popups.",
            "W,A,S,D for chamber select.",
            "Render drops: " + std::to_string(engine.get_swap_drops()),
            "Engine bytes: " + std::to_string(engine.get_bytes()),
            "",
            "Copyright (C) 2026 Gustav Louw",
            "ensim.cc, ensim.hh, sdl.cc, raylib.cc",
            "Licensed under GNU AGPLv3.",
            "See: https://www.gnu.org/licenses/agpl-3.0.html",
        };
        sdl.fill(rect);
        sdl.write(point, strings);
    }
};

struct ui
{
    static constexpr size_t history_capacity = 512;
    std::vector<float> history = {};
    std::vector<std::unique_ptr<cell>> base = {};
    std::vector<std::unique_ptr<cell>> popups = {};
    ensim::engine& engine;
    int x_select = 0;
    int y_select = 0;
    bool done = false;
    double throttle = 0.0;

    ui(ensim::engine& engine)
        : engine(engine)
        , y_select(engine.get_piston_y())
        {
            engine.set_logger(x_select, y_select);
            if(engine.get_height() > g_signals_count)
            {
                throw std::runtime_error("max signals supported: " + std::to_string(g_signals_count));
            }
            for(size_t y = 0; y < engine.get_height(); y++)
            for(size_t x = 0; x < engine.get_width(); x++)
            {
                base.push_back(std::make_unique<chamber>(x, y, x_select, y_select, engine));
            }
            for(size_t y = 0; y < engine.get_height(); y++)
            {
                base.push_back(std::make_unique<plot>(y, engine));
            }
            for(size_t y = 0; y < engine.get_height(); y++)
            for(size_t x = 0; x < engine.get_width(); x++)
            {
                base.push_back(std::make_unique<port>(x, y, engine));
            }
            base.push_back(std::make_unique<frame>());
        }

    std::unique_ptr<popup> make_popup()
    {
        const std::vector<double>& volume = engine.get_volume_signal_m3();
        const std::vector<double>& temperature = engine.get_static_temperature_signal_k();
        const std::vector<double>& pressure = engine.get_static_pressure_signal_pa();
        const std::vector<float>& audio = engine.get_audio_signal();
        const std::vector<float>& impulse = engine.get_impulse_signal();
        const size_t next = popups.size() + 1;
        switch(next)
        {
        case 1: return std::make_unique<gauge_popup>(next, "angular velocity (r/s)", engine);
        case 2: return std::make_unique<audio_popup>(next, "audio", audio);
        case 3: return std::make_unique<audio_popup>(next, "audio buffer history", history);
        case 4: return std::make_unique<pipe_popup> (next, "pipe pressure", engine);
        case 5: return std::make_unique<plot_popup> (next, "static pressure (p) volume (m3) diagram", volume, pressure);
        case 6: return std::make_unique<plot_popup> (next, "static temperature (k) volume (m3) diagram", volume, temperature);
        case 7: return std::make_unique<audio_popup>(next, "impulse signal", impulse);
        case 8: return std::make_unique<help_popup> (next, engine);
        }
        return nullptr;
    }

    void push_popup()
    {
        if(std::unique_ptr<cell> popup = make_popup())
        {
            popups.push_back(std::move(popup));
        }
    }

    void pop_popup()
    {
        if(popups.size() > 0)
        {
            popups.pop_back();
        }
    }

    void draw(sdl& sdl)
    {
        engine.set_swap_lock_on();
        for(const auto& elem : base)
        {
            elem->draw(sdl);
        }
        for(const auto& popup : popups)
        {
            popup->draw(sdl);
        }
        engine.set_swap_lock_off();
    }

    void log(sdl& sdl)
    {
        if(history.size() == history_capacity)
        {
            history.erase(history.begin());
        }
        history.push_back(sdl.get_audio_buffer_size());
    }

    void poll(sdl& sdl)
    {
        const std::vector<SDL_Event> events = sdl.poll();
        for(const auto& event : events)
        {
            if(event.type == SDL_EVENT_QUIT)
            {
                done = true;
            }
            if(event.type == SDL_EVENT_KEY_DOWN)
            {
                if(event.key.key == SDLK_0)
                {
                    engine.set_throttle_open_ratio(throttle = 0.00);
                    engine.set_injection_off();
                }
                if(event.key.key == SDLK_1)
                {
                    engine.set_throttle_open_ratio(throttle = 0.00);
                    engine.set_injection_on();
                }
                if(event.key.key == SDLK_2)
                {
                    engine.set_throttle_open_ratio(throttle = 0.33);
                    engine.set_injection_on();
                }
                if(event.key.key == SDLK_3)
                {
                    engine.set_throttle_open_ratio(throttle = 0.66);
                    engine.set_injection_on();
                }
                if(event.key.key == SDLK_4)
                {
                    engine.set_throttle_open_ratio(throttle = 0.99);
                    engine.set_injection_on();
                }
                if(event.key.key == SDLK_J)
                {
                    engine.set_throttle_open_ratio(0.0);
                    engine.disengage_clutch();
                }
                if(event.key.key == SDLK_K)
                {
                    engine.set_throttle_open_ratio(0.0);
                    engine.disengage_clutch();
                }
                if(event.key.key == SDLK_W) y_select -= 1;
                if(event.key.key == SDLK_S) y_select += 1;
                if(event.key.key == SDLK_D) x_select += 1;
                if(event.key.key == SDLK_A) x_select -= 1;
                x_select %= engine.get_width();
                y_select %= engine.get_height();
                engine.set_logger(x_select, y_select);
                if(event.key.key == SDLK_E) pop_popup();
                if(event.key.key == SDLK_Q) push_popup();
            }
            if(event.type == SDL_EVENT_KEY_UP)
            {
                if(event.key.key == SDLK_J)
                {
                    engine.set_throttle_open_ratio(throttle);
                    engine.decrement_gear();
                    engine.engage_clutch();
                }
                if(event.key.key == SDLK_K)
                {
                    engine.set_throttle_open_ratio(throttle);
                    engine.increment_gear();
                    engine.engage_clutch();
                }
            }
        }
    }
};

int main(int argc, const char* const*)
{
    if(argc == 2)
    {
        auto engine = ensim::new_engine(ensim::type::inline8);
        engine->set_logger(0, 0);
        engine->run(ensim::g_sample_rate_hz);
        printf("%lu bytes\n", engine->get_bytes());
        return 0;
    }
    std::atomic<bool> done = false;
    auto engine = ensim::new_engine(ensim::type::trx450r);
    sdl sdl;
    std::jthread thread(
        [&done, &engine, &sdl]()
        {
            while(not done)
            {
                static constexpr int audio_buffer_minimum_size = 1024;
                if(sdl.get_audio_buffer_size() < audio_buffer_minimum_size)
                {
                    engine->run(200);
                    sdl.buffer_audio(engine->get_audio_signal());
                }
                else
                {
                    using namespace std::chrono_literals;
                    std::this_thread::sleep_for(1ms);
                }
            }
        }
    );
    ui ui(*engine);
    while(not ui.done)
    {
        sdl.clear();
        ui.draw(sdl);
        ui.log(sdl);
        ui.poll(sdl);
        sdl.render();
    }
    done = true;
}
