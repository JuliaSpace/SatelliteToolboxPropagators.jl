## Description #############################################################################
#
# Define the function `copy` for the structures created here.
#
############################################################################################

function Base.copy(orbp::OrbitPropagatorSgp4{Tepoch, T}) where {Tepoch <: Number, T <: Number} 
    return OrbitPropagatorSgp4(copy(orbp.sgp4d))
end
