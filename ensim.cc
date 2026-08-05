#include "ensim.hh"

#include <array>
#include <exception>
#include <numbers>
#include <cmath>
#include <mutex>

#define fn __attribute__((used))

namespace ensim
{
#include "parts/lane.hh"
#include "parts/const.hh"
#include "parts/mathematics.hh"
#include "parts/flow.hh"
#include "parts/sparkplugs.hh"
#include "parts/cam.hh"
#include "parts/flywheel.hh"
#include "parts/pistons.hh"
#include "parts/throttle.hh"
#include "parts/limiter.hh"
#include "parts/crankshaft.hh"
#include "parts/pipe.hh"
#include "parts/diags.hh"
#include "parts/filter.hh"
#include "parts/mailbox.hh"
#include "parts/as_engine.hh"
#include "parts/generic_atv.hh"
#include "parts/new_engine.hh"
}
