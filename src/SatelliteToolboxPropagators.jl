module SatelliteToolboxPropagators

using Dates
using LinearAlgebra
using Printf
using StyledStrings

using ForwardDiff
using Reexport
using SatelliteToolboxOrbitDataMessages
using StaticArrays

@reexport using SatelliteToolboxBase
@reexport using SatelliteToolboxSgp4

import PrecompileTools

############################################################################################
#                                           API                                            #
############################################################################################

# Orbit propagators API.
include("./api/Propagators.jl")

using .Propagators
export Propagators, OrbitPropagator

############################################################################################
#                                          Types                                           #
############################################################################################

include("./types.jl")

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./api/j2.jl")
include("./api/j2osc.jl")
include("./api/j4.jl")
include("./api/j4osc.jl")
include("./api/sgp4.jl")
include("./api/twobody.jl")

include("./propagators/osculating.jl")
include("./propagators/j2.jl")
include("./propagators/j2osc.jl")
include("./propagators/j4.jl")
include("./propagators/j4osc.jl")
include("./propagators/twobody.jl")
include("./propagators/fit.jl")

include("./show.jl")

include("./precompile.jl")

end # module SatelliteToolboxPropagators
