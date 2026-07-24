#include <SDL3/SDL.h>
#include <list>
#include <algorithm>
#include <memory>
#include <format>

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

    std::list<SDL_Event> poll()
    {
        std::list<SDL_Event> events;
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

struct plot : cell
{
    static constexpr int h_p = sdl::h_p / g_signals;

    int y;
    ensim::engine& engine;

    plot(int y, ensim::engine& engine)
        : y(y)
        , engine(engine)
        {
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

    ensim::line downsample(const ensim::line& line)
    {
        static constexpr size_t max = 1024;
        const size_t size = line.size();
        if(size <= max)
        {
            return line;
        }
        ensim::line out;
        out.reserve(max);
        const size_t step = size / max;
        for(size_t i = 0, j = 0; i < size; i += step)
        {
            out.push_back(line[i]);
            if(++j == max)
            {
                break;
            }
        }
        return out;
    }

    points project(const ensim::line& line, const rect& rect, const uint32_t color)
    {
        const auto [min, max] = minmax(line);
        const size_t size = line.size();
        points points(color);
        points.self.reserve(size);
        for(size_t i = 0; i < size; i++)
        {
            const float x_ratio = float(i) / (size - 1);
            const float y_ratio = (line[i] - min) / (max - min);
            points.self.push_back(rect.project(x_ratio, y_ratio));
        }
        return points;
    }

    void draw(sdl& sdl) override
    {
        const int x_p = engine.get_w() * chamber::w_p;
        const int y_p = y * h_p;
        const rect rect(x_p, y_p, sdl::w_p - x_p, h_p, sdl::black);
        const std::string_view name = engine.get_signal_name(y);
        const ensim::line& raw = engine.get_signal(y);
        if(not raw.empty())
        {
            const ensim::line signal = normalize(downsample(raw));
            const points points = project(signal, rect, sdl::red);
            const auto [min, max] = minmax(raw);
            const ensim::real div = max / min;
            const int xm_p = x_p + sdl::font_p;
            const int ym_p = y_p + sdl::font_p;
            const std::array<message, 4> messages = {
                message(xm_p, ym_p + 0 * sdl::font_p, sdl::white, std::string(name)),
                message(xm_p, ym_p + 1 * sdl::font_p, sdl::white, std::format("max {:.6f}", max)),
                message(xm_p, ym_p + 2 * sdl::font_p, sdl::white, std::format("min {:.6f}", min)),
                message(xm_p, ym_p + 3 * sdl::font_p, sdl::white, std::format("div {:.6f}", div)),
            };
            sdl.fill(rect);
            sdl.draw(points);
            for(const auto& x : messages)
            {
                sdl.write(x);
            }
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
    std::list<std::unique_ptr<cell>> data;
    int x_select = 0;
    int y_select = 0;
    bool done = false;
    size_t engine_w;
    size_t engine_h;

    ui(ensim::engine& engine)
        : engine_w(engine.get_w())
        , engine_h(engine.get_h())
    {
        if(engine_w >= g_signals)
        {
            throw std::runtime_error(std::format("max {} signals supported", g_signals));
        }
        for(size_t y = 0; y < engine_h; y++)
        for(size_t x = 0; x < engine_w; x++)
        {
            append(std::make_unique<chamber>(x, y, x_select, y_select, engine));
        }
        for(size_t y = 0; y < engine_h; y++)
        {
            append(std::make_unique<plot>(y, engine));
        }
        for(size_t y = 0; y < engine_h; y++)
        for(size_t x = 0; x < engine_w; x++)
        {
            append(std::make_unique<port>(x, y, engine));
        }
        append(std::make_unique<frame>());
    }

    void append(std::unique_ptr<cell> cell)
    {
        data.push_back(std::move(cell));
    }

    void draw(sdl& sdl)
    {
        for(const auto& x : data)
        {
            x->draw(sdl);
        }
    }

    void poll(sdl& sdl)
    {
        const std::list<SDL_Event> events = sdl.poll();
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
                x_select %= engine_w;
                y_select %= engine_h;
            }
        }
    }
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
            ui.poll(sdl);
            sdl.render();
        }
    }
}
