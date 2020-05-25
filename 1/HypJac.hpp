#ifndef HYPJAC_HPP
#define HYPJAC_HPP

#ifndef HYPJAC_MAX_DEVICES
#define HYPJAC_MAX_DEVICES 1u
#else /* HYPJAC_MAX_DEVICES */
#error HYPJAC_MAX_DEVICES not definable externally
#endif /* ?HYPJAC_MAX_DEVICES */

#ifndef HYPJAC_MAX_LEVELS
#define HYPJAC_MAX_LEVELS 2u
#else /* HYPJAC_MAX_LEVELS */
#error HYPJAC_MAX_LEVELS not definable externally
#endif /* ?HYPJAC_MAX_LEVELS */

#include "HypJacL.hpp"

#endif /* !HYPJAC_HPP */
