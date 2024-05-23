module SatelliteToolboxPropagators

using Crayons
using Dates
using Printf
using LinearAlgebra
using Reexport
using SatelliteToolboxSgp4
using StaticArrays

@reexport using SatelliteToolboxBase
@reexport using SatelliteToolboxSgp4

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
#                                        Constants                                         #
############################################################################################

# Escape sequences related to the crayons.
const _D = Crayon(reset = true)
const _B = crayon"bold"
const _Y = crayon"bold yellow"

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./api/j2.jl")
include("./api/j2osc.jl")
include("./api/j4.jl")
include("./api/j4osc.jl")
include("./api/sgp4.jl")
include("./api/twobody.jl")

include("./propagators/j2.jl")
include("./propagators/j2osc.jl")
include("./propagators/j4.jl")
include("./propagators/j4osc.jl")
include("./propagators/twobody.jl")

include("./precompile.jl")

end # module SatelliteToolboxPropagators
