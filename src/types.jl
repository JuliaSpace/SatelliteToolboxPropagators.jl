# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Definition of types and structures.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export J2PropagatorConstants, J2Propagator, J2OsculatingPropagator
export J4PropagatorConstants, J4Propagator, J4OsculatingPropagator
export TwoBodyPropagator
export OrbitPropagatorJ2
export OrbitPropagatorJ2Osculating
export OrbitPropagatorJ4
export OrbitPropagatorJ4Osculating
export OrbitPropagatorSgp4
export OrbitPropagatorTwoBody

############################################################################################
#                                   J2 Orbit Propagator
############################################################################################

"""
    struct J2PropagatorConstants{T<:Number}

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

    # Constructors
    # ======================================================================================

    J2Propagator{Tepoch, T}(args...) where {Tepoch<:Number, T<:Number} = new(args...)
    J2Propagator{Tepoch, T}() where {Tepoch<:Number, T<:Number} = new()
end

############################################################################################
#                              J2 Osculating Orbit Propagator
############################################################################################

"""
    mutable struct J2OsculatingPropagator{Tepoch<:Number, T<:Number}

J2 osculating orbit propagator structure.
"""
mutable struct J2OsculatingPropagator{Tepoch<:Number, T<:Number}
    # J2 orbit propagator to propagate the mean elements.
    j2d::J2Propagator{Tepoch, T}

    # Propagation time from epoch.
    Δt::T

    # Current osculating Keplerian elements.
    orbk::KeplerianElements{Tepoch, T}

    # Constructors
    # ======================================================================================

    J2OsculatingPropagator{Tepoch, T}(args...) where {Tepoch<:Number, T<:Number} = new(args...)
    J2OsculatingPropagator{Tepoch, T}() where {Tepoch<:Number, T<:Number} = new()
end

############################################################################################
#                                   J4 Orbit Propagator
############################################################################################

"""
    struct J4PropagatorConstants{T<:Number}

Constants for the J4 orbit propagator.

# Fields

- `R0::T`: Earth equatorial radius [m].
- `μm::T`: √(GM / R0^3) [er/s]^(3/2).
- `J2::T`: The second gravitational zonal harmonic of the Earth.
- `J4::T`: The fourth gravitational zonal harmonic of the Earth.
"""
struct J4PropagatorConstants{T<:Number}
    R0::T
    μm::T
    J2::T
    J4::T
end

"""
    J4Propagator{Tepoch, T}

J4 orbit propagator structure.
"""
mutable struct J4Propagator{Tepoch, T}
    orb₀::KeplerianElements{Tepoch, T} # ............ Initial mean orbit elements [SI units]
    orbk::KeplerianElements{Tepoch, T} # ............ Current mean orbit elements [SI units]
    dn_o2::T                           # ... First time derivative of mean motion [rad / s²]
    ddn_o6::T                          # .. Second time derivative of mean motion [rad / s³]
    j4c::J4PropagatorConstants{T}      # .............................. Propagator constants
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

    # Constructors
    # ======================================================================================

    J4Propagator{Tepoch, T}(args...) where {Tepoch<:Number, T<:Number} = new(args...)
    J4Propagator{Tepoch, T}() where {Tepoch<:Number, T<:Number} = new()
end

############################################################################################
#                              J4 Osculating Orbit Propagator
############################################################################################

"""
    mutable struct J4OsculatingPropagator{Tepoch<:Number, T<:Number}

J4 osculating orbit propagator structure.
"""
mutable struct J4OsculatingPropagator{Tepoch<:Number, T<:Number}
    # J4 orbit propagator to propagate the mean elements.
    j4d::J4Propagator{Tepoch, T}

    # Propagation time from epoch.
    Δt::T

    # Current osculating Keplerian elements.
    orbk::KeplerianElements{Tepoch, T}

    # Constructors
    # ======================================================================================

    J4OsculatingPropagator{Tepoch, T}(args...) where {Tepoch<:Number, T<:Number} = new(args...)
    J4OsculatingPropagator{Tepoch, T}() where {Tepoch<:Number, T<:Number} = new()
end

############################################################################################
#                                   Two Body Propagator
############################################################################################

"""
    mutable struct TwoBodyPropagator{Tepoch<:Number, T<:Number}

Two body orbit propagator structure.
"""
mutable struct TwoBodyPropagator{Tepoch<:Number, T<:Number}
    orb₀::KeplerianElements{Tepoch, T} # ............ Initial mean orbit elements [SI units]
    orbk::KeplerianElements{Tepoch, T} # ............ Current mean orbit elements [SI units]
    μ::T                               # . Central body std. gravitational parameter [m³/s²]
    Δt::T                              # ..... Timespan from the initial elements' epoch [s]

    # Auxiliary Variables
    # ======================================================================================

    M₀::T  # ............................................... Initial mean mean anomaly [rad]
    n₀::T  # ........................................................  Mean motion [rad / s]

    # Constructors
    # ======================================================================================

    TwoBodyPropagator{Tepoch, T}(args...) where {Tepoch<:Number, T<:Number} = new(args...)
    TwoBodyPropagator{Tepoch, T}() where {Tepoch<:Number, T<:Number} = new()
end

############################################################################################
#                                           API
############################################################################################

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

#                              J2 Osculating Orbit Propagator
# ==========================================================================================

"""
    OrbitPropagatorJ2Osculating{Tepoch, T} <: OrbitPropagator{Tepoch, T}

J2 osculating orbit propagator.

# Fields

- `j2oscd`: Structure that stores the J2 osculating orbit propagator data (see
    [`J2OsculatingPropagator`](@ref)).
"""
struct OrbitPropagatorJ2Osculating{Tepoch<:Number, T<:Number} <: OrbitPropagator{Tepoch, T}
    j2oscd::J2OsculatingPropagator{Tepoch, T}
end

#                                   J4 Orbit Propagator
# ==========================================================================================

"""
    OrbitPropagatorJ4{Tepoch, T} <: OrbitPropagator{Tepoch, T}

J4 orbit propagator.

# Fields

- `j4d`: Structure that stores the J4 orbit propagator data (see [`J4Propagator`](@ref)).
"""
struct OrbitPropagatorJ4{Tepoch<:Number, T<:Number} <: OrbitPropagator{Tepoch, T}
    j4d::J4Propagator{Tepoch, T}
end

#                              J4 Osculating Orbit Propagator
# ==========================================================================================

"""
    OrbitPropagatorJ4Osculating{Tepoch, T} <: OrbitPropagator{Tepoch, T}

J4 osculating orbit propagator.

# Fields

- `j4oscd`: Structure that stores the J4 osculating orbit propagator data (see
    [`J4OsculatingPropagator`](@ref)).
"""
struct OrbitPropagatorJ4Osculating{Tepoch<:Number, T<:Number} <: OrbitPropagator{Tepoch, T}
    j4oscd::J4OsculatingPropagator{Tepoch, T}
end

#                                  SGP4 Orbit Propagator
# ==========================================================================================

"""
    OrbitPropagatorSgp4{Tepoch, T} <: OrbitPropagator{Tepoch, T}

SGP4 orbit propagator.

# Fields

- `sgp4d`: Structure that stores the SGP4 orbit propagator data.
"""
struct OrbitPropagatorSgp4{Tepoch<:Number, T<:Number} <: OrbitPropagator{Tepoch, T}
    sgp4d::Sgp4Propagator{Tepoch, T}
end

#                                  Two Orbit Propagator
# ==========================================================================================

"""
    OrbitPropagatorTwoBody{Tepoch, T} <: OrbitPropagator{Tepoch, T}

Two body orbit propagator.

# Fields

- `tbd`: Structure that stores the two body orbit propagator data (see
    [`TwoBodyPropagator`](@ref)).
"""
struct OrbitPropagatorTwoBody{Tepoch<:Number, T<:Number} <: OrbitPropagator{Tepoch, T}
    tbd::TwoBodyPropagator{Tepoch, T}
end
