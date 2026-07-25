#include <SDL3/SDL.h>

#include "ensim.hh"

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
};

struct rect
{
    SDL_FRect self;
    uint32_t color;
    float ratio;

    rect(const int x_p, const int y_p, const int w_p, const int h_p,
         const uint32_t color, const float ratio = 1.0f)
        : self(x_p, y_p, w_p, h_p)
        , color(color)
        , ratio(ratio)
        {
        }

    SDL_FPoint project(const float x_ratio, const float y_ratio) const
    {
        SDL_FPoint out;
        out.x = self.x + x_ratio * self.w;
        out.y = self.y + (1.0f - y_ratio) * self.h;
        return out;
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
    static constexpr uint32_t white = 0xFFFFFFFF;
    static constexpr uint32_t black = 0xFF000000;
    static constexpr uint32_t green = 0xFF00FF00;
    static constexpr uint32_t red   = 0xFFFF0000;
    static constexpr uint32_t font_p = 16;

    SDL_Window* window;
    SDL_Renderer* renderer;

    sdl()
    {
        SDL_Init(SDL_INIT_VIDEO);
        SDL_CreateWindowAndRenderer(nullptr, w_p, h_p, SDL_WINDOW_FULLSCREEN, &window, &renderer);
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
        set_color(rect.color, rect.ratio);
        SDL_RenderFillRect(renderer, &rect.self);
        set_color(white);
        SDL_RenderRect(renderer, &rect.self);
    }

    void write(const message& message)
    {
        set_color(message.color);
        SDL_RenderDebugText(renderer, message.self.x, message.self.y, message.string.data());
    }

    void draw(const points& points)
    {
        const size_t size = points.self.size();
        if(size > 0)
        {
            set_color(points.color);
            SDL_RenderLines(renderer, points.self.data(), size);
        }
    }
};

std::pair<ensim::real, ensim::real> minmax(const ensim::line& line)
{
    const auto [min_it, max_it] = std::minmax_element(line.begin(), line.end());
    const ensim::real min = (min_it == line.end()) ? ensim::real() : *min_it;
    const ensim::real max = (max_it == line.end()) ? ensim::real() : *max_it;
    return { min, max };
}

static constexpr int g_signals = 9;

struct cell
{
    virtual void draw(sdl& sdl) = 0;
    virtual ~cell() = default;
};

struct chamber : cell
{
    static constexpr int w_p = sdl::w_p / 32;
    static constexpr int h_p = sdl::h_p / g_signals;

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
        const rect container(x_p, y_p, w_p, h_p, sdl::black);
        sdl.fill(container);
        const int ws_p = w_p / 2;
        const int hs_p = w_p / 2;
        const int xs_p = x_p + (w_p - ws_p) / 2;
        const int ys_p = y_p + (h_p - hs_p) / 2;
        if(engine.get_panic_status(x, y))
        {
            const rect select(xs_p, ys_p, ws_p, hs_p, sdl::red);
            sdl.fill(select);
        }
        if(x == x_select and y == y_select)
        {
            const rect select(xs_p, ys_p, ws_p, hs_p, sdl::green);
            sdl.fill(select);
        }
    }
};

struct frame : cell
{
    void draw(sdl& sdl) override
    {
        const int x_p = 0;
        const int y_p = 0;
        const rect rect(x_p, y_p, sdl::w_p, sdl::h_p, sdl::white);
        sdl.outline(rect);
    }
};

ensim::line lingen(const size_t size)
{
    ensim::line linear;
    for(size_t i = 0; i < size; i++)
    {
        linear.push_back(i);
    }
    return linear;
}

ensim::line normalize(const ensim::line& line)
{
    const auto [min, max] = minmax(line);
    ensim::line out = line;
    for(auto& x : out)
    {
        x /= max;
    }
    return out;
}

ensim::line downsample(const ensim::line& line, const size_t size)
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

points project(const ensim::line& xx, const ensim::line& yy, const rect& rect, const uint32_t color)
{
    const auto [x_min, x_max] = minmax(xx);
    const auto [y_min, y_max] = minmax(yy);
    const size_t size = xx.size();
    points points(color);
    points.self.reserve(size);
    for(size_t i = 0; i < size; i++)
    {
        const float x_ratio = (xx[i] - x_min) / (x_max - x_min);
        const float y_ratio = (yy[i] - y_min) / (y_max - y_min);
        points.self.push_back(rect.project(x_ratio, y_ratio));
    }
    return points;
}

struct plot : cell
{
    static constexpr size_t max_points = 512;
    static constexpr int h_p = sdl::h_p / g_signals;

    int y;
    ensim::engine& engine;

    plot(const int y, ensim::engine& engine)
        : y(y)
        , engine(engine)
        {
        }

    void draw(sdl& sdl) override
    {
        const int x_p = engine.get_w() * chamber::w_p;
        const int y_p = y * h_p;
        const rect rect(x_p, y_p, sdl::w_p - x_p, h_p, sdl::black);
        const std::string_view name = engine.get_signal_name(y);
        const ensim::line& y_signal = engine.get_signal(y);
        if(y_signal.empty())
        {
            return;
        }
        const ensim::line xx = lingen(max_points);
        const ensim::line yy = normalize(downsample(y_signal, max_points));
        const points points = project(xx, yy, rect, sdl::red);
        const auto [min, max] = minmax(y_signal);
        const ensim::real div = max / min;
        const int xm_p = x_p + sdl::font_p;
        const int ym_p = y_p + sdl::font_p;
        const std::array<message, 4> messages = {
            message(xm_p, ym_p + 0 * sdl::font_p, sdl::white, std::string(name)),
            message(xm_p, ym_p + 1 * sdl::font_p, sdl::white, "max " + std::to_string(max)),
            message(xm_p, ym_p + 2 * sdl::font_p, sdl::white, "min " + std::to_string(min)),
            message(xm_p, ym_p + 3 * sdl::font_p, sdl::white, "div " + std::to_string(div)),
        };
        sdl.fill(rect);
        sdl.draw(points);
        for(const auto& x : messages)
        {
            sdl.write(x);
        }
    }
};

struct popup : cell
{
    static constexpr size_t max_points = 1024;
    static constexpr int w_p = sdl::h_p / 2;
    static constexpr int h_p = sdl::h_p / 2;
    static constexpr int dw_p = w_p / 8;
    static constexpr int dh_p = h_p / 8;

    int index;
    std::string name;
    const ensim::line& x_signal;
    const ensim::line& y_signal;

    popup(const int index, const std::string& name, const ensim::line& x_signal, const ensim::line& y_signal)
        : index(index)
        , name(name)
        , x_signal(x_signal)
        , y_signal(y_signal)
        {
        }

    void draw(sdl& sdl) override
    {
        const int x_p = sdl::w_p - w_p - index * dw_p;
        const int y_p = index * dh_p;
        const rect rect(x_p, y_p, w_p, h_p, sdl::black);
        if(x_signal.empty() or y_signal.empty())
        {
            return;
        }
        const auto [x_min, x_max] = minmax(x_signal);
        const auto [y_min, y_max] = minmax(y_signal);
        const ensim::line xx = normalize(downsample(x_signal, max_points));
        const ensim::line yy = normalize(downsample(y_signal, max_points));
        const points points = project(xx, yy, rect, sdl::red);
        const int xm_p = x_p + sdl::font_p;
        const int ym_p = y_p + sdl::font_p;
        const std::array<message, 5> messages = {
            message(xm_p, ym_p + 0 * sdl::font_p, sdl::white, name),
            message(xm_p, ym_p + 1 * sdl::font_p, sdl::white, "x_min = " + std::to_string(x_min)),
            message(xm_p, ym_p + 2 * sdl::font_p, sdl::white, "x_max = " + std::to_string(x_max)),
            message(xm_p, ym_p + 3 * sdl::font_p, sdl::white, "y_min = " + std::to_string(y_min)),
            message(xm_p, ym_p + 4 * sdl::font_p, sdl::white, "y_max = " + std::to_string(y_max)),
        };
        sdl.fill(rect);
        sdl.draw(points);
        for(const auto& x : messages)
        {
            sdl.write(x);
        }
    }
};

struct port : cell
{
    static constexpr int w_p = chamber::w_p / 4;
    static constexpr int h_p = chamber::w_p / 4;

    int x;
    int y;
    ensim::engine& engine;

    port(int x, int y, ensim::engine& engine): x(x), y(y), engine(engine) {}

    void draw(sdl& sdl) override
    {
        const int x_p = chamber::w_p * x + chamber::w_p / 2 - w_p / 2;
        const int y_p = chamber::h_p * y + chamber::h_p / 1 - h_p;
        const float ratio = engine.get_port_open_ratio(x, y);
        const rect rect(x_p, y_p, w_p, h_p, sdl::green, ratio);
        sdl.fill(rect);
    }
};

struct ui
{
    int x_select = 0;
    int y_select = 0;
    bool done = false;

    ui(ensim::engine& engine)
        : engine_w(engine.get_w())
        , engine_h(engine.get_h())
    {
        if(engine_w >= g_signals)
        {
            throw std::runtime_error("max signals supported: " + std::to_string(g_signals));
        }
        for(size_t y = 0; y < engine_h; y++)
        for(size_t x = 0; x < engine_w; x++)
        {
            base.push_back(std::make_unique<chamber>(x, y, x_select, y_select, engine));
        }
        for(size_t y = 0; y < engine_h; y++)
        {
            base.push_back(std::make_unique<plot>(y, engine));
        }
        for(size_t y = 0; y < engine_h; y++)
        for(size_t x = 0; x < engine_w; x++)
        {
            base.push_back(std::make_unique<port>(x, y, engine));
        }
        base.push_back(std::make_unique<frame>());
    }

    void push_popup(ensim::engine& engine)
    {
        const int next = popups.size() + 1;
        if(next == 1)
        {
            popups.push_back(
                std::make_unique<popup>(
                    next,
                    "pressure_volume_diagram",
                    engine.get_volume_signal_m3(),
                    engine.get_static_pressure_signal_pa()
                )
            );
        }
        if(next == 2)
        {
            popups.push_back(
                std::make_unique<popup>(
                    next,
                    "temperature_volume_diagram",
                    engine.get_volume_signal_m3(),
                    engine.get_static_temperature_signal_k()
                )
            );
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
                if(event.key.key == SDLK_W) y_select -= 1;
                if(event.key.key == SDLK_S) y_select += 1;
                if(event.key.key == SDLK_D) x_select += 1;
                if(event.key.key == SDLK_A) x_select -= 1;
                if(event.key.key == SDLK_E) pop_popup();
                if(event.key.key == SDLK_Q) push_popup(engine);
                x_select %= engine_w;
                y_select %= engine_h;
            }
        }
    }

private:
    size_t engine_w;
    size_t engine_h;
    std::vector<std::unique_ptr<cell>> base;
    std::vector<std::unique_ptr<cell>> popups;
};

int main(int argc, char**)
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
            engine->run(512, ui.x_select, ui.y_select);
            sdl.clear();
            ui.draw(sdl);
            ui.poll(sdl, *engine);
            sdl.render();
        }
    }
}
