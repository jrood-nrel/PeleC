#ifndef TRANSPORTTYPES_H
#define TRANSPORTTYPES_H

namespace pele::physics {
namespace transport {
struct ConstTransport;
struct SimpleTransport;
} // namespace transport
#ifdef USE_CONSTANT_TRANSPORT
using TransportType = transport::ConstTransport;
#elif USE_SIMPLE_TRANSPORT
using TransportType = transport::SimpleTransport;
#else
static_assert(false, "Invalid Transport specified");
#endif
} // namespace pele::physics
#endif
