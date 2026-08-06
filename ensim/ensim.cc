#include "ensim.hh"

#include <array>
#include <exception>
#include <numbers>
#include <cmath>
#include <mutex>

#define fn __attribute__((used))

namespace ensim
{
#include "lane.hh"
#include "const.hh"
#include "mathematics.hh"
#include "flow.hh"
#include "sparkplugs.hh"
#include "cam.hh"
#include "flywheel.hh"
#include "pistons.hh"
#include "throttle.hh"
#include "limiter.hh"
#include "crankshaft.hh"
#include "pipe.hh"
#include "diags.hh"
#include "filter.hh"
#include "mailbox.hh"
#include "as_engine.hh"
#include "inline4.hh"
#include "new_engine.hh"
}
