#pragma once

std::unique_ptr<engine> new_engine(const type type)
{
    std::unique_ptr<engine> engine;
    if(type == type::inline4)
    {
        engine = std::make_unique<ensim::inline4>();
    }
    engine->reset();
    return engine;
}
