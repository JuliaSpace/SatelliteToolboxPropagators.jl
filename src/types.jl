## Description #############################################################################
#
#   Definition of types and structures.
#
############################################################################################

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
#                                   J2 Orbit Propagator                                    #
############################################################################################

"""
    struct J2PropagatorConstants{T<:Number}

Constants for the J2 orbit propagator.

# Fields

- `R0::T`: Earth equatorial radius [m].
- `μm::T`: Standard gravitational parameter normalized by the Earth equatorial radius,
    √(GM / R0³) [rad / s].
- `J2::T`: The second gravitational zonal harmonic of the Earth.
"""
struct J2PropagatorConstants{T <: Number}
    R0::T
    μm::T
    J2::T
end

function J2PropagatorConstants{T}(j2c::J2PropagatorConstants) where {T <: Number}
    return J2PropagatorConstants{T}(T(j2c.R0), T(j2c.μm), T(j2c.J2))
end

function Base.convert(
    ::Type{J2PropagatorConstants{T}}, j2c::J2PropagatorConstants
) where {T <: Number}
    return J2PropagatorConstants{T}(j2c)
end

"""
    mutable struct J2Propagator{Tepoch<:Number, T<:Number}

J2 orbit propagator structure.
"""
mutable struct J2Propagator{Tepoch <: Number, T <: Number}
    # Initial mean orbit elements [SI units].
    orb₀::KeplerianElements{MeanAnomaly, Tepoch, T}

    # Current mean orbit elements [SI units].
    orbk::KeplerianElements{MeanAnomaly, Tepoch, T}

    # Propagator constants.
    j2c::J2PropagatorConstants{T}

    # Timespan from the initial elements' epoch [s].
    Δt::T

    # == Auxiliary Variables ===============================================================

    ∂Ω::T # ................................................. RAAN time derivative [rad / s]
    ∂ω::T # .................................. Argument of perigee time derivative [rad / s]
    n̄::T  # ................................................ Perturbed mean motion [rad / s]

    # == Constructors ======================================================================

    J2Propagator{Tepoch, T}(args...) where {Tepoch <: Number, T <: Number} = new(args...)
    J2Propagator{Tepoch, T}() where {Tepoch <: Number, T <: Number} = new()
end

############################################################################################
#                              J2 Osculating Orbit Propagator                              #
############################################################################################

"""
    mutable struct J2OsculatingPropagator{Tepoch<:Number, T<:Number}

J2 osculating orbit propagator structure.
"""
mutable struct J2OsculatingPropagator{Tepoch <: Number, T <: Number}
    # J2 orbit propagator to propagate the mean elements.
    j2d::J2Propagator{Tepoch, T}

    # Propagation time from epoch.
    Δt::T

    # Current osculating Keplerian elements.
    orbk::KeplerianElements{TrueAnomaly, Tepoch, T}

    # == Constructors ======================================================================

    J2OsculatingPropagator{Tepoch, T}(args...) where {Tepoch <: Number, T <: Number} =
        new(args...)
    J2OsculatingPropagator{Tepoch, T}() where {Tepoch <: Number, T <: Number} = new()
end

############################################################################################
#                                   J4 Orbit Propagator                                    #
############################################################################################

"""
    struct J4PropagatorConstants{T<:Number}

Constants for the J4 orbit propagator.

# Fields

- `R0::T`: Earth equatorial radius [m].
- `μm::T`: Standard gravitational parameter normalized by the Earth equatorial radius,
    √(GM / R0³) [rad / s].
- `J2::T`: The second gravitational zonal harmonic of the Earth.
- `J4::T`: The fourth gravitational zonal harmonic of the Earth.
"""
struct J4PropagatorConstants{T <: Number}
    R0::T
    μm::T
    J2::T
    J4::T
end

function J4PropagatorConstants{T}(j4c::J4PropagatorConstants) where {T <: Number}
    return J4PropagatorConstants{T}(T(j4c.R0), T(j4c.μm), T(j4c.J2), T(j4c.J4))
end

function Base.convert(
    ::Type{J4PropagatorConstants{T}}, j4c::J4PropagatorConstants
) where {T <: Number}
    return J4PropagatorConstants{T}(j4c)
end

"""
    mutable struct J4Propagator{Tepoch<:Number, T<:Number}

J4 orbit propagator structure.
"""
mutable struct J4Propagator{Tepoch <: Number, T <: Number}
    # Initial mean orbit elements [SI units].
    orb₀::KeplerianElements{MeanAnomaly, Tepoch, T}

    # Current mean orbit elements [SI units].
    orbk::KeplerianElements{MeanAnomaly, Tepoch, T}

    # Propagator constants.
    j4c::J4PropagatorConstants{T}

    # Timespan from the initial elements' epoch [s].
    Δt::T

    # == Auxiliary Variables ===============================================================

    ∂Ω::T # ................................................. RAAN time derivative [rad / s]
    ∂ω::T # .................................. Argument of perigee time derivative [rad / s]
    n̄::T  # ................................................ Perturbed mean motion [rad / s]

    # == Constructors ======================================================================

    J4Propagator{Tepoch, T}(args...) where {Tepoch <: Number, T <: Number} = new(args...)
    J4Propagator{Tepoch, T}() where {Tepoch <: Number, T <: Number} = new()
end

############################################################################################
#                              J4 Osculating Orbit Propagator                              #
############################################################################################

"""
    mutable struct J4OsculatingPropagator{Tepoch<:Number, T<:Number}

J4 osculating orbit propagator structure.
"""
mutable struct J4OsculatingPropagator{Tepoch <: Number, T <: Number}
    # J4 orbit propagator to propagate the mean elements.
    j4d::J4Propagator{Tepoch, T}

    # Propagation time from epoch.
    Δt::T

    # Current osculating Keplerian elements.
    orbk::KeplerianElements{TrueAnomaly, Tepoch, T}

    # == Constructors ======================================================================

    J4OsculatingPropagator{Tepoch, T}(args...) where {Tepoch <: Number, T <: Number} =
        new(args...)
    J4OsculatingPropagator{Tepoch, T}() where {Tepoch <: Number, T <: Number} = new()
end

############################################################################################
#                                   Two Body Propagator                                    #
############################################################################################

"""
    mutable struct TwoBodyPropagator{Tepoch<:Number, T<:Number}

Two body orbit propagator structure.
"""
mutable struct TwoBodyPropagator{Tepoch <: Number, T <: Number}
    # Initial mean orbit elements [SI units].
    orb₀::KeplerianElements{MeanAnomaly, Tepoch, T}

    # Current mean orbit elements [SI units].
    orbk::KeplerianElements{MeanAnomaly, Tepoch, T}

    # Central body std. gravitational parameter [m³/s²].
    μ::T

    # Timespan from the initial elements' epoch [s].
    Δt::T

    # == Auxiliary Variables ===============================================================

    n₀::T  # ........................................................  Mean motion [rad / s]

    # == Constructors ======================================================================

    TwoBodyPropagator{Tepoch, T}(args...) where {Tepoch <: Number, T <: Number} =
        new(args...)
    TwoBodyPropagator{Tepoch, T}() where {Tepoch <: Number, T <: Number} = new()
end

############################################################################################
#                                           API                                            #
############################################################################################

# == J2 Orbit Propagator ===================================================================

"""
    OrbitPropagatorJ2{Tepoch, T} <: OrbitPropagator{Tepoch, T}

J2 orbit propagator.

# Fields

- `j2d`: Structure that stores the J2 orbit propagator data (see [`J2Propagator`](@ref)).
"""
struct OrbitPropagatorJ2{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}
    j2d::J2Propagator{Tepoch, T}
end

# == J2 Osculating Orbit Propagator ========================================================

"""
    OrbitPropagatorJ2Osculating{Tepoch, T} <: OrbitPropagator{Tepoch, T}

J2 osculating orbit propagator.

# Fields

- `j2oscd`: Structure that stores the J2 osculating orbit propagator data (see
    [`J2OsculatingPropagator`](@ref)).
"""
struct OrbitPropagatorJ2Osculating{Tepoch <: Number, T <: Number} <:
       OrbitPropagator{Tepoch, T}
    j2oscd::J2OsculatingPropagator{Tepoch, T}
end

# == J4 Orbit Propagator ===================================================================

"""
    OrbitPropagatorJ4{Tepoch, T} <: OrbitPropagator{Tepoch, T}

J4 orbit propagator.

# Fields

- `j4d`: Structure that stores the J4 orbit propagator data (see [`J4Propagator`](@ref)).
"""
struct OrbitPropagatorJ4{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}
    j4d::J4Propagator{Tepoch, T}
end

# == J4 Osculating Orbit Propagator ========================================================

"""
    OrbitPropagatorJ4Osculating{Tepoch, T} <: OrbitPropagator{Tepoch, T}

J4 osculating orbit propagator.

# Fields

- `j4oscd`: Structure that stores the J4 osculating orbit propagator data (see
    [`J4OsculatingPropagator`](@ref)).
"""
struct OrbitPropagatorJ4Osculating{Tepoch <: Number, T <: Number} <:
       OrbitPropagator{Tepoch, T}
    j4oscd::J4OsculatingPropagator{Tepoch, T}
end

# == SGP4 Orbit Propagator =================================================================

"""
    OrbitPropagatorSgp4{Tepoch, T} <: OrbitPropagator{Tepoch, T}

SGP4 orbit propagator.

# Fields

- `sgp4d`: Structure that stores the SGP4 orbit propagator data.
"""
struct OrbitPropagatorSgp4{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}
    sgp4d::Sgp4Propagator{Tepoch, T}
end

# == Two Body Orbit Propagator =================================================

"""
    OrbitPropagatorTwoBody{Tepoch, T} <: OrbitPropagator{Tepoch, T}

Two body orbit propagator.

# Fields

- `tbd`: Structure that stores the two body orbit propagator data (see
    [`TwoBodyPropagator`](@ref)).
"""
struct OrbitPropagatorTwoBody{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}
    tbd::TwoBodyPropagator{Tepoch, T}
end

############################################################################################
#                                        Julia API                                         #
############################################################################################

# Union of the low-level propagator structures defined in this package.
const _PropagatorData{Tepoch, T} = Union{
    J2Propagator{Tepoch, T},
    J2OsculatingPropagator{Tepoch, T},
    J4Propagator{Tepoch, T},
    J4OsculatingPropagator{Tepoch, T},
    TwoBodyPropagator{Tepoch, T},
}

"""
    Base.copy(pd::P) where {P <: _PropagatorData} -> P

Create a copy of the propagator structure `pd`. The fields that are propagator structures
themselves, such as the J2 propagator inside the J2 osculating propagator, are copied
recursively so that the copy can be propagated independently of `pd`.
"""
function Base.copy(pd::P) where {P <: _PropagatorData}
    return P(ntuple(i -> _copy_field(getfield(pd, i)), Val(fieldcount(P)))...)
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _copy_field(x) -> typeof(x)

Return the value stored in a propagator field when copying the structure. Nested propagator
structures are copied, whereas every other field is immutable and is returned as is.
"""
_copy_field(x) = x
_copy_field(pd::_PropagatorData) = copy(pd)
