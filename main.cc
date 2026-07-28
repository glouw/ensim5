#include <SDL3/SDL.h>
#include <algorithm>
#include <numbers>
#include <cmath>
#include <list>

#include "ensim.hh"

using ensim::operator""_r;

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
    uint32_t color;

    points(const uint32_t color): color(color) {}

    void append(const SDL_FPoint& point)
    {
        self.push_back(point);
    }

    void append(const float x, const float y)
    {
        append({x, y});
    }
};

struct rect
{
    SDL_FRect self;
    uint32_t color;
    float color_ratio;

    rect(const int x_p, const int y_p, const int w_p, const int h_p, const uint32_t color, const float color_ratio = 1.0f)
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

    SDL_FPoint project(const float x_ratio, const float y_ratio) const
    {
        SDL_FPoint out;
        out.x = self.x + x_ratio * self.w;
        out.y = self.y + (1.0f - y_ratio) * self.h;
        return out;
    }
};

struct circle
{
    SDL_FPoint self;
    float radius;
    uint32_t color;

    circle(const int x_p, const int y_p, const uint32_t color)
        : self(x_p, y_p)
        , color(color)
        {
        }

    circle(const rect& rect, const uint32_t color, const float border_ratio = 1.0f)
    {
        self.x = rect.self.x + rect.self.w * 0.5f;
        self.y = rect.self.y + rect.self.h * 0.5f;
        radius = std::min(rect.self.w, rect.self.h) * 0.5f * border_ratio;
        this->color = color;
    }
};

struct message : point
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
    static constexpr uint32_t line_p = 16;
    static constexpr uint32_t font_p = 8;

    static constexpr uint32_t grey   = 0xFFAAAAAA;
    static constexpr uint32_t white  = 0xFFFFFFFF;
    static constexpr uint32_t black  = 0xFF101010;
    static constexpr uint32_t green  = 0xFF00FF00;
    static constexpr uint32_t purple = 0xFFFF00FF;
    static constexpr uint32_t blue   = 0xFF00AAFF;
    static constexpr uint32_t orange = 0xFFFFAA00;
    static constexpr uint32_t red    = 0xFFFF2222;
    static constexpr uint32_t yellow = 0xFFFFFF00;

    sdl()
    {
        SDL_Init(SDL_INIT_VIDEO);
        SDL_CreateWindowAndRenderer(nullptr, w_p, h_p, SDL_WINDOW_FULLSCREEN, &window, &renderer);
        SDL_SetRenderVSync(renderer, true);
    }

    ~sdl()
    {
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

    void set_color(const uint32_t hex, const float ratio = 1.0f)
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
        set_color(0x0);
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
        float x_p = message.self.x;
        float y_p = message.self.y;
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

private:
    SDL_Window* window;
    SDL_Renderer* renderer;
};

std::pair<float, float> minmax(const ensim::line& line)
{
    const auto [min_it, max_it] = std::minmax_element(line.begin(), line.end());
    const float min = (min_it == line.end()) ? 0.0f : *min_it;
    const float max = (max_it == line.end()) ? 0.0f : *max_it;
    return { min, max };
}

struct signals
{
    static constexpr int count = 9;

    ensim::line lingen(const size_t size) const
    {
        ensim::line linear;
        for(size_t i = 0; i < size; i++)
        {
            linear.push_back(i);
        }
        return linear;
    }

    ensim::line normalize(const ensim::line& line) const
    {
        const auto [min, max] = minmax(line);
        ensim::line out = line;
        for(auto& x : out)
        {
            x = (max == 0.0f) ? 1.0f : x / max;
        }
        return out;
    }

    ensim::line downsample(const ensim::line& line, const size_t size) const
    {
        ensim::line out;
        out.reserve(size);
        for(size_t i = 0; i < size; i++)
        {
            const size_t j = i * line.size() / size;
            out.push_back(line[j]);
        }
        return out;
    }

    points project(const ensim::line& xx, const ensim::line& yy, const rect& rect) const
    {
        const auto [x_min, x_max] = minmax(xx);
        const auto [y_min, y_max] = minmax(yy);
        const size_t size = xx.size();
        points points(rect.color);
        points.self.reserve(size);
        for(size_t i = 0; i < size; i++)
        {
            const float dx = x_max - x_min;
            const float dy = y_max - y_min;
            const float x_ratio = dx == 0.0f ? 1.0f : ((xx[i] - x_min) / dx);
            const float y_ratio = dy == 0.0f ? 1.0f : ((yy[i] - y_min) / dy);
            points.append(rect.project(x_ratio, y_ratio));
        }
        return points;
    }

    points project_1d(const ensim::line& signal, const rect& rect, const size_t size) const
    {
        const ensim::line xx = lingen(size);
        const ensim::line yy = normalize(downsample(signal, size));
        return project(xx, yy, rect);
    }

    points project_2d(const ensim::line& x_signal, const ensim::line& y_signal, const rect& rect, const size_t size) const
    {
        const ensim::line xx = normalize(downsample(x_signal, size));
        const ensim::line yy = normalize(downsample(y_signal, size));
        return project(xx, yy, rect);
    }
};

struct cell
{
    virtual void draw(sdl& sdl) = 0;
    virtual ~cell() = default;
};

struct chamber : cell
{
    static constexpr int w_p = sdl::w_p / 32;
    static constexpr int h_p = sdl::h_p / signals::count;
    static constexpr int ws_p = w_p / 3;
    static constexpr int hs_p = w_p / 3;

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
        sdl.fill(container);
        const int xs_p = x_p + (w_p - ws_p) / 2;
        const int ys_p = y_p + (h_p - hs_p) / 2;
        const int xf_p = x_p + w_p / 2;
        const int yf_p = y_p + h_p / 2;
        const int y_throttle = engine.get_throttle_y();
        const int y_piston = engine.get_piston_y();
        const int y_audio = engine.get_audio_y();
        if(engine.get_panic_status(x, y))
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
        if(x == x_select and y == y_select)
        {
            sdl.fill(rect(xs_p, ys_p, ws_p, hs_p, select_color));
        }
    }

private:
    int x;
    int y;
    int& x_select;
    int& y_select;
    ensim::engine& engine;
    static constexpr uint32_t fill_color = sdl::black;
    static constexpr uint32_t throttle_color = sdl::orange;
    static constexpr uint32_t piston_color = sdl::purple;
    static constexpr uint32_t audio_color = sdl::blue;
    static constexpr uint32_t panic_color = sdl::red;
    static constexpr uint32_t select_color = sdl::yellow;
};

struct frame : cell
{
    void draw(sdl& sdl) override
    {
        const int x_p = 0;
        const int y_p = 0;
        const rect rect(x_p, y_p, sdl::w_p, sdl::h_p, frame_color);
        sdl.outline(rect);
    }

private:
    static constexpr uint32_t frame_color = sdl::white;
};

struct port : cell
{
    static constexpr int w_p = chamber::w_p / 4;
    static constexpr int h_p = chamber::w_p / 4;

    port(int x, int y, ensim::engine& engine)
        : x(x)
        , y(y)
        , engine(engine)
        {
        }

    void draw(sdl& sdl) override
    {
        const int x_p = chamber::w_p * x + chamber::w_p / 2 - w_p / 2;
        const int y_p = chamber::h_p * y + chamber::h_p / 1 - h_p;
        const float ratio = engine.get_port_open_ratio(x, y);
        const rect rect(x_p, y_p, w_p, h_p, fill_color, ratio);
        sdl.fill(rect);
    }

private:
    int x;
    int y;
    ensim::engine& engine;
    static constexpr uint32_t fill_color = sdl::green;
};

struct plot : signals, cell
{
    static constexpr size_t max_points = sdl::w_p;
    static constexpr int h_p = sdl::h_p / signals::count;

    plot(const int y, ensim::engine& engine)
        : y(y)
        , engine(engine)
        {
        }

    void draw(sdl& sdl) override
    {
        const std::string_view name = engine.get_signal_name(y);
        const ensim::line& y_signal = engine.get_signal(y);
        if(y_signal.empty())
        {
            return;
        }
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
        };
        const int yz_p = y_p + h_p * (max / (max - min));
        const point yz0_p(x_p, yz_p, zero_line_color);
        const point yz1_p(sdl::w_p, yz_p, zero_line_color);
        sdl.fill(fill);
        sdl.draw_line(yz0_p, yz1_p);
        sdl.draw_lines(data);
        sdl.write(font, strings);
    }

private:
    int y;
    ensim::engine& engine;
    static constexpr uint32_t signal_color = sdl::red;
    static constexpr uint32_t fill_color = sdl::black;
    static constexpr uint32_t text_color = sdl::white;
    static constexpr uint32_t zero_line_color = sdl::grey;
};

struct popup : cell
{
    static constexpr size_t max_points = 1024;
    static constexpr int w_p = sdl::h_p / 2;
    static constexpr int h_p = sdl::h_p / 2;
    static constexpr int dw_p = 2 * sdl::line_p;
    static constexpr int dh_p = 2 * sdl::line_p;

    popup(const int index)
        : index(index)
        {
        }

protected:
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

private:
    int index;
};

struct plot_popup : signals, popup
{
    plot_popup(const int index, const std::string& name, const ensim::line& x_signal, const ensim::line& y_signal)
        : popup(index)
        , name(name)
        , x_signal(x_signal)
        , y_signal(y_signal)
        {
        }

    void draw(sdl& sdl) override
    {
        if(x_signal.empty() or y_signal.empty())
        {
            return;
        }
        const rect fill = calc_rect();
        const rect signal(fill, signal_color);
        const points data = project_2d(x_signal, y_signal, signal, max_points);
        const point font(fill.self.x + sdl::line_p, fill.self.y + sdl::line_p, text_color);
        const auto [x_min, x_max] = minmax(x_signal);
        const auto [y_min, y_max] = minmax(y_signal);
        const std::vector<std::string> strings = {
            name,
            "x min = " + std::to_string(x_min),
            "x max = " + std::to_string(x_max),
            "y min = " + std::to_string(y_min),
            "y max = " + std::to_string(y_max),
        };
        sdl.fill(fill);
        sdl.draw_lines(data);
        sdl.write(font, strings);
    }

private:
    std::string name;
    const ensim::line& x_signal;
    const ensim::line& y_signal;
    static constexpr uint32_t signal_color = sdl::red;
    static constexpr uint32_t text_color = sdl::white;
};

struct audio_popup : signals, popup
{
    audio_popup(const int index, const ensim::line& audio_signal)
        : popup(index)
        , audio_signal(audio_signal)
        {
        }

    void draw(sdl& sdl) override
    {
        if(audio_signal.empty())
        {
            return;
        }
        const rect fill = calc_rect();
        const rect signal(fill, signal_color);
        const points data = project_1d(audio_signal, signal, max_points);
        const point font(fill.self.x + sdl::line_p, fill.self.y + sdl::line_p, text_color);
        const auto [y_min, y_max] = minmax(audio_signal);
        const std::vector<std::string> strings = {
            name,
            "min = " + std::to_string(y_min),
            "max = " + std::to_string(y_max),
            "samples = " + std::to_string(audio_signal.size()),
        };
        sdl.fill(fill);
        sdl.draw_lines(data);
        sdl.write(font, strings);
    }

private:
    const ensim::line& audio_signal;
    static constexpr std::string name = "audio_signal";
    static constexpr uint32_t signal_color = sdl::red;
    static constexpr uint32_t text_color = sdl::white;
};

struct gauge_popup : popup
{
    gauge_popup(const int index, const std::string& name, const ensim::real& value, const float max_value, const size_t needle_ticks)
        : popup(index)
        , name(name)
        , value(value)
        , max_value(max_value)
        , needle_ticks(needle_ticks)
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
        const float at = value;
        const std::vector<std::string> strings = {
            name,
            std::to_string(at),
            "max: " + std::to_string(max_value),
        };
        const circle outer(rect, outer_color, outer_ratio);
        const circle inner(rect, inner_color, inner_ratio);
        sdl.fill(rect);
        sdl.draw_circle(outer);
        sdl.draw_circle(inner);
        draw_ticks(sdl, outer);
        draw_needle(sdl, outer, at);
        sdl.write(text, strings);
    }

private:

    float to_angle(const float at) const
    {
        return start_theta_r - (at / max_value) * sweep_theta_r;
    }

    void draw_needle(sdl& sdl, const circle& outer, const float at) const
    {
        const point middle(outer.self.x, outer.self.y, needle_color);
        const float angle_r = to_angle(at);
        const float radius = needle_ratio * outer.radius;
        const point tip(
            middle.self.x + std::cos(angle_r) * radius,
            middle.self.y - std::sin(angle_r) * radius,
            needle_color
        );
        sdl.draw_line(middle, tip);
    }

    void draw_ticks(sdl& sdl, const circle& outer) const
    {
        const float step = max_value / needle_ticks;
        const float radius = outer.radius * ticks_ratio;
        for(size_t i = 0; i <= needle_ticks; i++)
        {
            const float tick = i * step;
            const float angle_r = to_angle(tick);
            const message message(
                outer.self.x + std::cos(angle_r) * radius,
                outer.self.y - std::sin(angle_r) * radius,
                ticks_color,
                std::to_string(static_cast<int>(tick))
            );
            sdl.write(message, true);
        }
    }

    const std::string name;
    const ensim::real& value;
    const float max_value;
    const size_t needle_ticks;
    static constexpr float start_theta_r = (4.0f / 3.0f) * std::numbers::pi_v<float>;
    static constexpr float sweep_theta_r = (5.0f / 3.0f) * std::numbers::pi_v<float>;
    static constexpr float outer_ratio = 0.75f;
    static constexpr float inner_ratio = 0.05f;
    static constexpr float ticks_ratio = 0.82f;
    static constexpr float needle_ratio = 0.85f;
    static constexpr uint32_t needle_color = sdl::red;
    static constexpr uint32_t inner_color = sdl::grey;
    static constexpr uint32_t outer_color = sdl::grey;
    static constexpr uint32_t ticks_color = sdl::white;
    static constexpr uint32_t text_color = sdl::white;
};

struct help_popup : popup
{
    help_popup(const int index)
        : popup(index)
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
            "help!",
            "Q,E to cycle through these popups.",
            "W,A,S,D for chamber select.",
            "AGPL V3."
        };
        sdl.fill(rect);
        sdl.write(point, strings);
    }

private:
    static constexpr uint32_t text_color = sdl::white;
};

struct ui
{
    ui(ensim::engine& engine)
        : y_select(engine.get_piston_y())
    {
        if(engine.get_width() >= signals::count)
        {
            throw std::runtime_error("max signals supported: " + std::to_string(signals::count));
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

    std::unique_ptr<popup> make_popup(ensim::engine& engine)
    {
        const size_t next = popups.size() + 1;
        switch(next)
        {
        case 1:
            return std::make_unique<gauge_popup>(
                next,
                "crankshaft_angular_velocity_r_per_s",
                engine.get_crankshaft_angular_velocity_r_per_s(),
                2000.0,
                20
            );
        case 2:
            return std::make_unique<audio_popup>(
                next,
                engine.get_audio_signal()
            );
        case 3:
            return std::make_unique<plot_popup>(
                next,
                "pressure_volume_diagram",
                engine.get_volume_signal_m3(),
                engine.get_pressure_signal_pa()
            );
        case 4:
            return std::make_unique<plot_popup>(
                next,
                "temperature_volume_diagram",
                engine.get_volume_signal_m3(),
                engine.get_temperature_signal_k()
            );
        case 5:
            return std::make_unique<help_popup>(
                next
            );
        }
        return nullptr;
    }

    void push_popup(ensim::engine& engine)
    {
        if(std::unique_ptr<cell> popup = make_popup(engine))
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
        for(const auto& elem : base)
        {
            elem->draw(sdl);
        }
        for(const auto& popup : popups)
        {
            popup->draw(sdl);
        }
    }

    void poll(sdl& sdl, ensim::engine& engine)
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
                if(event.key.key == SDLK_1)
                {
                    engine.set_throttle_open_ratio(0.00_r);
                    engine.set_injection_off();
                }
                if(event.key.key == SDLK_2)
                {
                    engine.set_throttle_open_ratio(0.001_r);
                    engine.set_injection_on();
                }
                if(event.key.key == SDLK_3)
                {
                    engine.set_throttle_open_ratio(0.01_r);
                    engine.set_injection_on();
                }
                if(event.key.key == SDLK_4)
                {
                    engine.set_throttle_open_ratio(0.1_r);
                    engine.set_injection_on();
                }
                if(event.key.key == SDLK_5)
                {
                    engine.set_throttle_open_ratio(0.5_r);
                    engine.set_injection_on();
                }
                if(event.key.key == SDLK_W) y_select -= 1;
                if(event.key.key == SDLK_S) y_select += 1;
                if(event.key.key == SDLK_D) x_select += 1;
                if(event.key.key == SDLK_A) x_select -= 1;
                if(event.key.key == SDLK_E) pop_popup();
                if(event.key.key == SDLK_Q) push_popup(engine);
                x_select %= engine.get_width();
                y_select %= engine.get_height();
            }
        }
    }

private:
    std::vector<std::unique_ptr<cell>> base;
    std::vector<std::unique_ptr<cell>> popups;

public:
    int x_select = 0;
    int y_select = 0;
    bool done = false;
};

int main(int argc, const char* const*)
{
    auto engine = ensim::new_engine(ensim::type::inline8);
    if(argc == 2)
    {
        engine->run(44800);
    }
    else
    {
        ui ui(*engine);
        sdl sdl;
        while(not ui.done)
        {
            engine->run(800, ui.x_select, ui.y_select);
            sdl.clear();
            ui.draw(sdl);
            ui.poll(sdl, *engine);
            sdl.render();
        }
    }
}
