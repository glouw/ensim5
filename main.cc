#include <SDL3/SDL.h>
#include <iostream>
#include <string_view>
#include <cmath>
#include <algorithm>
#include "ensim.hh"

struct sdl_s
{
    static constexpr unsigned black  = 0x000000FF;
    static constexpr unsigned white  = 0xFFFFFFFF;
    static constexpr unsigned grey   = 0xAAAAAAFF;
    static constexpr unsigned red    = 0xFF0000FF;
    static constexpr unsigned yellow = 0xFFFF00FF;
    static constexpr unsigned green  = 0x00FF00FF;
    static constexpr unsigned blue   = 0x0000FFFF;
    static constexpr int xres_p = 1920;
    static constexpr int yres_p = 1080;
    const int engine_w;
    const int engine_h;
    const int engine_y;
    const double chamber_w_p = xres_p / 32.0;
    const double chamber_h_p = yres_p / engine_h;
    const double port_w_p = chamber_w_p / 5.0;
    const double port_h_p = chamber_w_p / 5.0;
    const double plot_x_p = chamber_w_p * engine_w;
    const double plot_y_p = 0;
    const double plot_w_p = xres_p - plot_x_p;
    const double plot_h_p = chamber_h_p;
    const double pv_window_w_p = 0.75 * yres_p;
    const double pv_window_h_p = 0.75 * yres_p;
    const double pv_window_x_p = xres_p / 2.0 - pv_window_w_p / 2.0;
    const double pv_window_y_p = yres_p / 2.0 - pv_window_h_p / 2.0;
    const double margin_p = 8;
    const double font_p = 16;
    const double border_p = 18;
    float mouse_x_p;
    float mouse_y_p;
    int chamber_select_x = 0;
    int chamber_select_y = engine_y;
    bool quit = false;
    bool pressure_volume_window = false;
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
            if(event.type == SDL_EVENT_QUIT)
            {
                quit = true;
            }
            if(event.type == SDL_EVENT_KEY_UP)
            {
                if(event.key.key == SDLK_Q)
                {
                    pressure_volume_window ^= true;
                }
            }
            if(event.type == SDL_EVENT_KEY_DOWN)
            {
                if(event.key.key == SDLK_W) { chamber_select_y -= 1; }
                if(event.key.key == SDLK_S) { chamber_select_y += 1; }
                if(event.key.key == SDLK_D) { chamber_select_x += 1; }
                if(event.key.key == SDLK_A) { chamber_select_x -= 1; }
                chamber_select_y = std::clamp(chamber_select_y, 0, engine_h - 1);
                chamber_select_x = std::clamp(chamber_select_x, 0, engine_w - 1);
            }
        }
    }

    void set_color(const unsigned color, const double brightness = 1.0)
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

    void draw_box(const SDL_FRect& rect, const unsigned color, const double brightness = 1.0)
    {
        SDL_FRect outline = rect;
        outline.w += 1;
        outline.h += 1;
        set_color(color, brightness);
        SDL_RenderFillRect(renderer, &rect);
        set_color(white);
        SDL_RenderRect(renderer, &outline);
    }

    void draw_chambers(
        const std::vector<std::vector<double>>& open_ratios,
        const std::vector<std::vector<double>>& panics)
    {
        /*
         * Chamber boxes.
         *
         */

        for(int y = 0; y < engine_h; y++)
        for(int x = 0; x < engine_w; x++)
        {
            SDL_FRect chamber;
            chamber.x = chamber_w_p * x;
            chamber.y = chamber_h_p * y;
            chamber.w = chamber_w_p;
            chamber.h = chamber_h_p;
            draw_box(chamber, black);
            if(panics[y][x])
            {
                chamber.x += border_p;
                chamber.y += border_p;
                chamber.w -= 2.0 * border_p;
                chamber.h -= 2.0 * border_p;
                set_color(red);
                SDL_RenderFillRect(renderer, &chamber);
            }
        }

        /*
         * Port boxes.
         *
         */

        for(int y = 0; y < engine_h; y++)
        for(int x = 0; x < engine_w; x++)
        {
            SDL_FRect chamber;
            chamber.x = chamber_w_p * x;
            chamber.y = chamber_h_p * y;
            chamber.w = chamber_w_p;
            chamber.h = chamber_h_p;
            SDL_FRect port;
            port.x = chamber.x + chamber.w / 2.0 - port_w_p / 2.0;
            port.y = chamber.y + chamber.h - port_h_p;
            port.w = port_w_p;
            port.h = port_h_p;
            draw_box(port, green, open_ratios[y][x]);
        }
    }

    void select_chamber()
    {
        const unsigned buttons = SDL_GetMouseState(&mouse_x_p, &mouse_y_p);
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
        const double select_x_p = chamber_w_p * chamber_select_x;
        const double select_y_p = chamber_h_p * chamber_select_y;
        SDL_FRect rect;
        rect.x = select_x_p + border_p;
        rect.y = select_y_p + border_p;
        rect.w = chamber_w_p - 2.0 * border_p;
        rect.h = chamber_h_p - 2.0 * border_p;
        draw_box(rect, chamber_select_y == engine_y ? yellow : green);
    }

    void draw_plots(const ensim::diags& diags)
    {
        for(int y = 0; y < diags.channels; y++)
        {
            const std::vector<double>& data = diags[y];
            if(data.empty())
            {
                continue;
            }
            SDL_FRect rect;
            rect.x = plot_x_p;
            rect.y = plot_y_p + y * plot_h_p;
            rect.w = plot_w_p;
            rect.h = plot_h_p;
            const double min_val = *std::min_element(data.begin(), data.end());
            const double max_val = *std::max_element(data.begin(), data.end());

            /*
             * Draw bounding box.
             *
             */

            draw_box(rect, black);

            /*
             * Draw zero line.
             *
             */

            if(max_val > min_val)
            {
                const double t = std::clamp((0.0 - min_val) / (max_val - min_val), 0.0, 1.0);
                const double zero_y = rect.y + rect.h * (1.0 - t);
                set_color(red);
                SDL_RenderLine(renderer, rect.x, zero_y, rect.x + rect.w, zero_y);
            }

            /*
             * Draw points.
             *
             */

            const double range = max_val - min_val;
            const size_t points = data.size();
            const size_t max_points = plot_w_p;
            const size_t point_count = std::min(points, max_points);
            std::vector<SDL_FPoint> plot_points;
            plot_points.reserve(point_count);
            for(size_t i = 0; i < point_count; i++)
            {
                size_t index = i;
                if(points > max_points)
                {
                    index = (i * (points - 1)) / (point_count - 1);
                }
                double normalized = 0.0;
                if(range != 0.0)
                {
                    normalized = ((data[index] - min_val) / range) * 2.0 - 1.0;
                }
                double x0 = rect.x;
                if(point_count > 1)
                {
                    x0 += ((double)i / (point_count - 1)) * rect.w;
                }
                const double y0 = rect.y + rect.h * 0.5 - normalized * rect.h * 0.5;
                plot_points.push_back({
                    (float) x0,
                    (float) y0,
                });
            }
            set_color(green);
            for(size_t i = 1; i < plot_points.size(); i++)
            {
                const size_t j = i - 1;
                SDL_RenderLine(renderer, plot_points[j].x, plot_points[j].y, plot_points[i].x, plot_points[i].y);
            }

            /*
             * Draw plot goodies.
             *
             */

            const double div_val = max_val / min_val;
            const double index = data.size() * (mouse_x_p - plot_x_p) / plot_w_p;
            set_color(white);
            const double xx = rect.x + margin_p;
            const double yy = rect.y + margin_p;
            const std::string max_string = "max: " + std::to_string(max_val);
            const std::string min_string = "min: " + std::to_string(min_val);
            const std::string sel_string = "sel: " + (mouse_x_p < plot_x_p ? "?" : std::to_string(data[index]));
            const std::string div_string = "div: " + std::to_string(div_val);
            const std::string siz_string = "siz: " + std::to_string(point_count);
            SDL_RenderDebugText(renderer, xx, yy + 0 * font_p, ensim::diags::name[y]);
            SDL_RenderDebugText(renderer, xx, yy + 1 * font_p, max_string.data());
            SDL_RenderDebugText(renderer, xx, yy + 2 * font_p, min_string.data());
            SDL_RenderDebugText(renderer, xx, yy + 3 * font_p, sel_string.data());
            SDL_RenderDebugText(renderer, xx, yy + 4 * font_p, siz_string.data());
        }
        if(mouse_x_p > plot_x_p)
        {
            set_color(white);
            SDL_RenderLine(renderer, mouse_x_p, 0.0, mouse_x_p, yres_p);
            const double ratio = (mouse_x_p - plot_x_p) / plot_w_p;
            const double radians = 4.0 * M_PI * ratio;
            const std::string unit = std::to_string(radians) + "rad";
            set_color(white);
            SDL_RenderDebugText(renderer, mouse_x_p + font_p, mouse_y_p, unit.data());
        }
    }

    void draw_pressure_volume(const ensim::diags& diags)
    {
        if(pressure_volume_window)
        {
            /*
             * Draw bounding box.
             *
             */

            SDL_FRect rect;
            rect.x = pv_window_x_p;
            rect.y = pv_window_y_p;
            rect.w = pv_window_w_p;
            rect.h = pv_window_h_p;
            draw_box(rect, black);

            /*
             * Draw points.
             *
             */

            const std::vector<double>& static_pressure_pa = diags[ensim::diags::channel::chamber_static_pressure_pa];
            const std::vector<double>& volume_m3 = diags[ensim::diags::channel::chamber_volume_m3];
            const size_t size = static_pressure_pa.size();
            if(size == 0)
            {
                return;
            }
            const double v_min = *std::min_element(volume_m3.begin(), volume_m3.end());
            const double v_max = *std::max_element(volume_m3.begin(), volume_m3.end());
            const double p_min = *std::min_element(static_pressure_pa.begin(), static_pressure_pa.end());
            const double p_max = *std::max_element(static_pressure_pa.begin(), static_pressure_pa.end());
            const double v_range = (v_max - v_min > 0.0) ? (v_max - v_min) : 1.0;
            const double p_range = (p_max - p_min > 0.0) ? (p_max - p_min) : 1.0;
            const double margin_ratio = 0.05;
            const double margin_x_p = rect.w * margin_ratio;
            const double margin_y_p = rect.h * margin_ratio;
            const double rect_x_p = rect.x + margin_x_p;
            const double rect_y_p = rect.y + margin_y_p;
            const double rect_w_p = rect.w - 2.0 * margin_x_p;
            const double rect_h_p = rect.h - 2.0 * margin_y_p;
            std::vector<SDL_FPoint> points;
            points.reserve(size);
            for(size_t i = 0; i < size; ++i)
            {
                const double x = rect_x_p + ((volume_m3[i] - v_min) / v_range) * rect_w_p;
                const double y = rect_y_p + rect_h_p - ((static_pressure_pa[i] - p_min) / p_range) * rect_h_p;
                points.push_back({
                    (float) x,
                    (float) y,
                });
            }
            set_color(green);
            SDL_RenderPoints(renderer, points.data(), points.size());
            set_color(white);
            SDL_RenderDebugText(renderer, rect.x + font_p, rect.y + font_p, "pv-curve (Q to close)");
        }
    }

    void draw_border()
    {
        const double x0 = 0.0;
        const double y0 = 0.0;
        const double x1 = xres_p - 1.0;
        const double y1 = yres_p - 1.0;
        SDL_RenderLine(renderer, x0, y0, x1, y0);
        SDL_RenderLine(renderer, x0, y1, x1, y1);
        SDL_RenderLine(renderer, x0, y0, x0, y1);
        SDL_RenderLine(renderer, x1, y0, x1, y1);
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

int main(int argc, char* argv[])
{
    std::vector<std::string_view> args(argv, argv + argc);
    std::unique_ptr<ensim::engine> engine = ensim::new_engine(ensim::engine::type::inline8);
    auto perf = std::find(args.begin(), args.end(), "--perf");
    if(perf != args.end())
    {
        engine->run(44800);
        std::cout << engine->bytes() << std::endl;
    }
    else
    {
        sdl_s sdl(engine->width(), engine->height(), engine->piston_y());
        while(!sdl.quit)
        {
            sdl.poll_quit();
            sdl.clear_screen();
            sdl.draw_chambers(engine->get_port_open_ratios(), engine->get_panics());
            sdl.select_chamber();
            sdl.draw_chamber_selection();
            sdl.draw_plots(engine->get_diags());
            sdl.draw_pressure_volume(engine->get_diags());
            sdl.draw_border();
            sdl.render();
            engine->run(512, sdl.chamber_select_x, sdl.chamber_select_y);
        }
    }
}
