#include <SDL3/SDL.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>
#include "ensim.hh"

struct sdl_s
{
    const int xres_p = 1920;
    const int yres_p = 1080;
    const int engine_w;
    const int engine_h;
    const int engine_y;
    const int plots = 9;
    const float chamber_w_p = xres_p / 32.0f;
    const float chamber_h_p = yres_p / engine_h;
    const float plot_x_p = chamber_w_p * engine_w;
    const float plot_y_p = 0;
    const float plot_w_p = xres_p - plot_x_p;
    const float plot_h_p = chamber_h_p;
    const float margin_p = 8;
    const float font_p = 16;
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
        if(engine_h != plots)
        {
            std::cerr << "demo only compatible with engine heights of 9" << std::endl;
            std::terminate();
        }
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

    void clear_screen()
    {
        SDL_SetRenderDrawColor(renderer, 0x00, 0x00, 0x00, 0xFF);
        SDL_RenderClear(renderer);
    }

    void draw_chambers()
    {
        SDL_SetRenderDrawColor(renderer, 0xAA, 0xAA, 0xAA, 0xFF);
        for(int x = 0; x < engine_w; x++)
        for(int y = 0; y < engine_h; y++)
        {
            const SDL_FRect rect = {
                .x = chamber_w_p * x,
                .y = chamber_h_p * y,
                .w = chamber_w_p,
                .h = chamber_h_p,
            };
            SDL_RenderRect(renderer, &rect);
        }
    }

    void select_chamber()
    {
        float mouse_x_p;
        float mouse_y_p;
        SDL_GetMouseState(&mouse_x_p, &mouse_y_p);
        chamber_select_x = std::floor(mouse_x_p / chamber_w_p);
        chamber_select_y = std::floor(mouse_y_p / chamber_h_p);
        chamber_select_x = std::min(chamber_select_x, engine_w - 1);
        chamber_select_y = std::min(chamber_select_y, engine_h - 1);
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
        chamber_select_y == engine_y
            ? SDL_SetRenderDrawColor(renderer, 0xAA, 0x00, 0x00, 0xFF)
            : SDL_SetRenderDrawColor(renderer, 0x00, 0x00, 0xAA, 0xFF);
        SDL_RenderFillRect(renderer, &rect);
    }

    void draw_plots(ensim::diags& diags)
    {
        for(int y = 0; y < plots; y++)
        {
            std::vector<float> data = diags[y];
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
                SDL_SetRenderDrawColor(renderer, 0xFF, 0x00, 0x00, 0xFF);
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
            SDL_SetRenderDrawColor(renderer, 0x00, 0xFF, 0x00, 0xFF);
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

            SDL_SetRenderDrawColor(renderer, 0xAA, 0xAA, 0xAA, 0xFF);
            SDL_RenderRect(renderer, &rect);
            SDL_SetRenderDrawColor(renderer, 0xFF, 0xFF, 0xFF, 0xFF);
            const float xx = rect.x + margin_p;
            const float yy = rect.y + margin_p;
            const std::string max_string = "max: " + std::to_string(max_val);
            const std::string min_string = "min: " + std::to_string(min_val);
            SDL_RenderDebugText(renderer, xx, yy + 0 * font_p, ensim::diags::name[y]);
            SDL_RenderDebugText(renderer, xx, yy + 1 * font_p, max_string.data());
            SDL_RenderDebugText(renderer, xx, yy + 2 * font_p, min_string.data());
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
    sdl_s sdl(engine->get_w(), engine->get_h(), engine->get_y());
    while(!sdl.quit)
    {
        sdl.poll_quit();
        sdl.clear_screen();
        sdl.draw_chambers();
        sdl.select_chamber();
        sdl.draw_chamber_selection();
        sdl.draw_plots(engine->get_diags());
        sdl.render();
        engine->run(512, sdl.chamber_select_x, sdl.chamber_select_y);
    }
    std::cout << engine->get_size() << std::endl;
#endif
}
