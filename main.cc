#include <SDL3/SDL.h>
#include <list>
#include <memory>
#include "ensim.hh"

struct sdl
{
    static constexpr int w_p = 1920;
    static constexpr int h_p = 1080;
    static constexpr uint32_t white = 0xFFFFFFFF;
    static constexpr uint32_t black = 0xFF000000;
    static constexpr uint32_t green = 0xFF00FF00;
    static constexpr uint32_t font_p = 8;

    SDL_Window* window;
    SDL_Renderer* renderer;

    sdl()
    {
        SDL_CreateWindowAndRenderer(nullptr, w_p, h_p, SDL_WINDOW_FULLSCREEN, &window, &renderer);
    }

    void set_color(const uint32_t hex, float ratio = 1.0f)
    {
        const int a = ratio * (hex >> 24);
        const int r = ratio * (hex >> 16);
        const int g = ratio * (hex >>  8);
        const int b = ratio * (hex >>  0);
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

    void delay(int ms)
    {
        SDL_Delay(ms);
    }

    void outline(const SDL_FRect& rect, const uint32_t color)
    {
        set_color(color);
        SDL_RenderRect(renderer, &rect);
    }

    void fill(const SDL_FRect& rect, const uint32_t color, const float ratio = 1.0f)
    {
        set_color(color, ratio);
        SDL_RenderFillRect(renderer, &rect);
        outline(rect, white);
    }

    void write(const SDL_FPoint& point, const std::string_view& message)
    {
        SDL_RenderDebugText(renderer, point.x, point.y, message.data());
    }

    ~sdl()
    {
        SDL_DestroyWindow(window);
        SDL_DestroyRenderer(renderer);
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
    static constexpr int h_p = sdl::h_p / ensim::channels;

    int x;
    int y;

    chamber(int x, int y): x(x), y(y) {}

    void draw(sdl& sdl) override
    {
        SDL_FRect rect;
        rect.x = w_p * x;
        rect.y = h_p * y;
        rect.w = w_p;
        rect.h = h_p;
        sdl.fill(rect, sdl::black);
    }
};

struct frame : cell
{
    void draw(sdl& sdl) override
    {
        SDL_FRect rect;
        rect.x = 0;
        rect.y = 0;
        rect.w = sdl::w_p;
        rect.h = sdl::h_p;
        sdl.outline(rect, sdl::white);
    }
};

struct plot : cell
{
    static constexpr int h_p = sdl::h_p / ensim::channels;

    int y;
    int x_p;
    int w_p;
    std::string_view name;

    plot(int y, int engine_w, const std::string_view& name)
        : y(y)
        , x_p(engine_w * chamber::w_p)
        , w_p(sdl::w_p - x_p)
        , name(name)
        {
        }

    void draw(sdl& sdl) override
    {
        SDL_FRect rect;
        rect.x = x_p;
        rect.y = h_p * y;
        rect.w = w_p;
        rect.h = h_p;
        sdl.fill(rect, sdl::black);

        SDL_FPoint point;
        point.x = rect.x + sdl::font_p;
        point.y = rect.y + sdl::font_p;
        sdl.write(point, name.data());
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
        ensim::grid grid = engine.get_port_open_ratios();
        SDL_FRect rect;
        rect.x = chamber::w_p * x + chamber::w_p / 2 - w_p / 2;
        rect.y = chamber::h_p * y + chamber::h_p / 1 - h_p / 2;
        rect.w = w_p;
        rect.h = h_p;
        sdl.fill(rect, sdl::black, grid[y][x]);
    }
};

struct ui
{
    std::list<std::unique_ptr<cell>> data;

    ui(ensim::engine& engine)
    {
        for(size_t x = 0; x < engine.get_w(); x++)
        for(size_t y = 0; y < engine.get_h(); y++)
        {
            append(std::make_unique<chamber>(x, y));
        }
        for(size_t y = 0; y < engine.get_h(); y++)
        {
            append(std::make_unique<plot>(y, engine.get_w(), ensim::diags::name[y]));
        }
        for(size_t x = 0; x < engine.get_w(); x++)
        for(size_t y = 0; y < engine.get_h(); y++)
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
};


int main()
{
    std::unique_ptr<ensim::engine> engine = ensim::new_engine(ensim::engine::type::inline8);
    ui ui(*engine);
    sdl sdl;
    sdl.clear();
    ui.draw(sdl);
    sdl.render();
    sdl.delay(5000);
}
