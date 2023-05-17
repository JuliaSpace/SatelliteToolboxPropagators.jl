module SatelliteToolboxPropagators

using Crayons
using Dates
using Reexport
using SatelliteToolboxSgp4

@reexport using SatelliteToolboxBase
@reexport using SatelliteToolboxSgp4

############################################################################################
#                                           API
############################################################################################

# Orbit propagators API.
include("./api/Propagators.jl")

using .Propagators
export Propagators, OrbitPropagator

############################################################################################
#                                          Types
############################################################################################

include("./types.jl")

############################################################################################
#                                        Constants
############################################################################################

# Escape sequences related to the crayons.
const _D = Crayon(reset = true)
const _B = crayon"bold"

############################################################################################
#                                         Includes
############################################################################################

include("./api/j2.jl")
include("./api/j2osc.jl")
include("./api/j4.jl")
include("./api/sgp4.jl")

include("./propagators/j2.jl")
include("./propagators/j2osc.jl")
include("./propagators/j4.jl")

end # module SatelliteToolboxPropagators
