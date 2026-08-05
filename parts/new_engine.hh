#pragma once

std::unique_ptr<engine> new_engine(const type type)
{
    std::unique_ptr<engine> engine;
    if(type == type::generic_atv)
    {
        engine = std::make_unique<generic_atv>();
    }
    engine->reset();
    return engine;
}
