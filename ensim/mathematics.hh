#pragma once

using std::sin;
using std::cos;
using std::fmax;
using std::fmin;
using std::log;
using std::sqrt;
using std::trunc;
using std::exp;

fn constexpr real clamper(const real value, const real lower, const real upper)
{
    return fmax(fmin(value, upper), lower);
}

fn constexpr real modulos(const real value, const real by)
{
    return value - trunc(value / by) * by;
}

fn real cuberoot(const real x)
{
    return exp(log(x) / 3.0_r);
}

fn real frand()
{
    const real random = 2.0_r * rand() / static_cast<real>(RAND_MAX);
    return random - 1.0_r;
}
