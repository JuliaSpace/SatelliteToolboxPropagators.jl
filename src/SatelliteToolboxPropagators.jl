module SatelliteToolboxPropagators

using Crayons
using Dates
using Reexport
using SatelliteToolboxSgp4
using SatelliteToolboxTle

@reexport using SatelliteToolboxBase

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
include("./propagators/j2.jl")

end # module SatelliteToolboxPropagators
