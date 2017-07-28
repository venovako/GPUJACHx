#ifndef HZ_HPP
#define HZ_HPP

#ifndef HZ_MAX_DEVICES
#define HZ_MAX_DEVICES 1u
#else // HZ_MAX_DEVICES
#error HZ_MAX_DEVICES not definable externally
#endif // !HZ_MAX_DEVICES

#ifndef HZ_MAX_LEVELS
#define HZ_MAX_LEVELS 2u
#else // HZ_MAX_LEVELS
#error HZ_MAX_LEVELS not definable externally
#endif // !HZ_MAX_LEVELS

#include "HZ_L.hpp"

#endif // !HZ_HPP
