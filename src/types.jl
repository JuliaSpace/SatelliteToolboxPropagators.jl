# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Definition of types and structures.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export J2PropagatorConstants, J2Propagator
export OrbitPropagator
export OrbirPropagatorJ2

############################################################################################
#                                   J2 Orbit Propagator
############################################################################################

"""
    J2PropagatorConstants{T}

Constants for the J2 orbit propagator.

# Fields

- `R0::T`: Earth equatorial radius [m].
- `μm::T`: √(GM / R0^3) [er/s]^(3/2).
- `J2::T`: The second gravitational zonal harmonic of the Earth.
"""
struct J2PropagatorConstants{T<:Number}
    R0::T
    μm::T
    J2::T
end

"""
    mutable struct J2Propagator{Tepoch<:Number, T<:Number}

J2 orbit propagator structure.
"""
mutable struct J2Propagator{Tepoch<:Number, T<:Number}
    orb₀::KeplerianElements{Tepoch, T} # ............ Initial mean orbit elements [SI units]
    orbk::KeplerianElements{Tepoch, T} # ............ Current mean orbit elements [SI units]
    dn_o2::T                           # ... First time derivative of mean motion [rad / s²]
    ddn_o6::T                          # .. Second time derivative of mean motion [rad / s³]
    j2c::J2PropagatorConstants{T}      # .............................. Propagator constants
    Δt::T                              # ..... Timespan from the initial elements' epoch [s]

    # Auxiliary Variables
    # ======================================================================================

    al₀::T    # ............................... Initial mean normalized semi-major axis [er]
    M₀::T     # ............................................ Initial mean mean anomaly [rad]
    δa::T     # ................................... Semi-major axis time derivative [er / s]
    δe::T     # ....................................... Eccentricity time derivative [1 / s]
    δΩ::T     # ............................................. RAAN time derivative [rad / s]
    δω::T     # .............................. Argument of perigee time derivative [rad / s]
    n̄::T      # ............................................... "Mean" mean motion [rad / s]
end

############################################################################################
#                                           API
############################################################################################

"""
    abstract type OrbitPropagator{Tepoch<:Number, T<:Number}

Abstract type for the orbit propagators.
"""
abstract type OrbitPropagator{Tepoch<:Number, T<:Number} end

#                                   J2 Orbit Propagator
# ==========================================================================================

"""
    OrbitPropagatorJ2{Tepoch, T} <: OrbitPropagator{Tepoch, T}

J2 orbit propagator.

# Fields

- `j2d`: Structure that stores the J2 orbit propagator data (see [`J2Propagator`](@ref)).
"""
struct OrbitPropagatorJ2{Tepoch<:Number, T<:Number} <: OrbitPropagator{Tepoch, T}
    j2d::J2Propagator{Tepoch, T}
end
