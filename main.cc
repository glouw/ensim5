#include <SDL3/SDL.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>
#include "ensim.hh"

struct sdl_s
{
    static constexpr unsigned black = 0x000000FF;
    static constexpr unsigned white = 0xFFFFFFFF;
    static constexpr unsigned grey  = 0xAAAAAAFF;
    static constexpr unsigned red   = 0xFF0000FF;
    static constexpr unsigned green = 0x00FF00FF;
    static constexpr unsigned blue  = 0x0000FFFF;
    static constexpr int xres_p = 1920;
    static constexpr int yres_p = 1080;
    const int engine_w;
    const int engine_h;
    const int engine_y;
    const float chamber_w_p = xres_p / 32.0f;
    const float chamber_h_p = yres_p / engine_h;
    const float port_w_p = chamber_w_p / 4.0f;
    const float port_h_p = chamber_w_p / 4.0f;
    const float plot_x_p = chamber_w_p * engine_w;
    const float plot_y_p = 0;
    const float plot_w_p = xres_p - plot_x_p;
    const float plot_h_p = chamber_h_p;
    const float margin_p = 8;
    const float font_p = 16;
    float mouse_x_p;
    float mouse_y_p;
    int chamber_select_x = 0;
    int chamber_select_y = 0;
    bool quit = false;
    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_Cursor* cursor;

    sdl_s(const int engine_w, const int engine_h, const int engine_y)
        : engine_w(engine_w)
        , engine_h(engine_h)
        , engine_y(engine_y)
    {
        SDL_Init(SDL_INIT_VIDEO);
        window = SDL_CreateWindow("ensim", xres_p, yres_p, SDL_WINDOW_FULLSCREEN);
        renderer = SDL_CreateRenderer(window, nullptr);
        cursor = SDL_CreateSystemCursor(SDL_SYSTEM_CURSOR_CROSSHAIR);
        SDL_SetCursor(cursor);
    }

    void poll_quit()
    {
        SDL_Event event;
        while(SDL_PollEvent(&event))
        {
            quit = event.type == SDL_EVENT_QUIT;
        }
    }

    void set_color(const unsigned color, const float brightness = 1.0f)
    {
        SDL_SetRenderDrawColor(
            renderer,
            brightness * (0xFF & (color >> 24)),
            brightness * (0xFF & (color >> 16)),
            brightness * (0xFF & (color >>  8)),
            brightness * (0xFF & (color >>  0))
        );
    }

    void clear_screen()
    {
        set_color(black);
        SDL_RenderClear(renderer);
    }

    void draw_chambers(const std::vector<std::vector<float>>& open_ratios)
    {
        for(int y = 0; y < engine_h; y++)
        for(int x = 0; x < engine_w; x++)
        {
            const SDL_FRect chamber = {
                .x = chamber_w_p * x,
                .y = chamber_h_p * y,
                .w = chamber_w_p,
                .h = chamber_h_p,
            };
            set_color(grey);
            SDL_RenderRect(renderer, &chamber);
        }
        for(int y = 0; y < engine_h; y++)
        for(int x = 0; x < engine_w; x++)
        {
            const SDL_FRect chamber = {
                .x = chamber_w_p * x,
                .y = chamber_h_p * y,
                .w = chamber_w_p,
                .h = chamber_h_p,
            };
            const SDL_FRect port = {
                .x = chamber.x + chamber.w / 2.0f - port_w_p / 2.0f,
                .y = chamber.y + chamber.h - port_h_p / 2.0f,
                .w = port_w_p,
                .h = port_h_p,
            };
            set_color(red, open_ratios[y][x]);
            SDL_RenderFillRect(renderer, &port);
        }
    }

    void select_chamber()
    {
        unsigned buttons = SDL_GetMouseState(&mouse_x_p, &mouse_y_p);
        if(buttons & SDL_BUTTON_LEFT)
        {
            chamber_select_x = std::floor(mouse_x_p / chamber_w_p);
            chamber_select_y = std::floor(mouse_y_p / chamber_h_p);
            chamber_select_x = std::min(chamber_select_x, engine_w - 1);
            chamber_select_y = std::min(chamber_select_y, engine_h - 1);
        }
    }

    void draw_chamber_selection()
    {
        const float select_x_p = chamber_w_p * chamber_select_x;
        const float select_y_p = chamber_h_p * chamber_select_y;
        const SDL_FRect rect = {
            .x = select_x_p,
            .y = select_y_p,
            .w = chamber_w_p,
            .h = chamber_h_p,
        };
        chamber_select_y == engine_y ? set_color(red) : set_color(blue);
        SDL_RenderFillRect(renderer, &rect);
    }

    void draw_plots(ensim::diags& diags)
    {
        if(mouse_x_p > plot_x_p)
        {
            set_color(grey);
            SDL_RenderLine(renderer, mouse_x_p, 0.0f, mouse_x_p, yres_p);
            const float ratio = (mouse_x_p - plot_x_p) / plot_w_p;
            const float radians = 4.0f * M_PI * ratio;
            const std::string unit = std::to_string(radians) + "rad";
            set_color(white);
            SDL_RenderDebugText(renderer, mouse_x_p + font_p, mouse_y_p, unit.data());
        }

        for(int y = 0; y < diags.channels; y++)
        {
            const std::vector<float>& data = diags[y];
            if(data.empty())
            {
                continue;
            }
            const SDL_FRect rect = {
                .x = plot_x_p,
                .y = plot_y_p + y * plot_h_p,
                .w = plot_w_p,
                .h = plot_h_p,
            };
            const float min_val = *std::min_element(data.begin(), data.end());
            const float max_val = *std::max_element(data.begin(), data.end());

            /*
             * Draw zero line.
             *
             */

            if(max_val > min_val)
            {
                const float t = std::clamp((0.0f - min_val) / (max_val - min_val), 0.0f, 1.0f);
                const float zero_y = rect.y + rect.h * (1.0f - t);
                set_color(red);
                SDL_RenderLine(renderer, rect.x, zero_y, rect.x + rect.w, zero_y);
            }

            /*
             * Draw points.
             *
             */

            const float range = max_val - min_val;
            const size_t points = data.size();
            const size_t max_points = plot_w_p;
            const size_t point_count = std::min(points, max_points);
            set_color(green);
            for(size_t i = 0; i < point_count; i++)
            {
                size_t index = i;
                if(points > max_points)
                {
                    index = (i * (points - 1)) / (point_count - 1);
                }
                float normalized = 0.0f;
                if(range != 0.0f)
                {
                    normalized = ((data[index] - min_val) / range) * 2.0f - 1.0f;
                }
                float x0 = rect.x;
                if(point_count > 1)
                {
                    x0 += ((float) i / (point_count - 1)) * rect.w;
                }
                const float y0 = rect.y + rect.h * 0.5f - normalized * rect.h * 0.5f;
                SDL_RenderPoint(renderer, x0, y0);
            }

            /*
             * Draw bounding box.
             *
             */

            set_color(grey);
            SDL_RenderRect(renderer, &rect);

            /*
             * Draw plot goodies.
             *
             */

            const float div_val = max_val / min_val;
            const float index = data.size() * (mouse_x_p - plot_x_p) / plot_w_p;
            set_color(white);
            const float xx = rect.x + margin_p;
            const float yy = rect.y + margin_p;
            const std::string max_string = "max: " + std::to_string(max_val);
            const std::string min_string = "min: " + std::to_string(min_val);
            const std::string sel_string = "sel: " + (mouse_x_p < plot_x_p ? "?" : std::to_string(data[index]));
            const std::string div_string = "div: " + std::to_string(div_val);
            SDL_RenderDebugText(renderer, xx, yy + 0 * font_p, ensim::diags::name[y]);
            SDL_RenderDebugText(renderer, xx, yy + 1 * font_p, max_string.data());
            SDL_RenderDebugText(renderer, xx, yy + 2 * font_p, min_string.data());
            SDL_RenderDebugText(renderer, xx, yy + 3 * font_p, sel_string.data());
            SDL_RenderDebugText(renderer, xx, yy + 4 * font_p, div_string.data());
        }
    }

    void render()
    {
        SDL_RenderPresent(renderer);
    }

    ~sdl_s()
    {
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_DestroyCursor(cursor);
        SDL_Quit();
    }
};

int main()
{
    std::unique_ptr<ensim::engine> engine = ensim::new_engine(ensim::engine::type::inline8);
#ifdef ENSIM_PERF
    auto t0 = std::chrono::high_resolution_clock::now();
    engine->run(44800);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    std::cout << "Watermark Percent: " << dt.count() / 1000.0 << std::endl;
#else
    sdl_s sdl(engine->width(), engine->height(), engine->piston_y());
    while(!sdl.quit)
    {
        sdl.poll_quit();
        sdl.clear_screen();
        sdl.draw_chambers(engine->get_port_open_ratios());
        sdl.select_chamber();
        sdl.draw_chamber_selection();
        sdl.draw_plots(engine->get_diags());
        sdl.render();
        engine->run(512, sdl.chamber_select_x, sdl.chamber_select_y);
    }
    std::cout << engine->bytes() << std::endl;
#endif
}
