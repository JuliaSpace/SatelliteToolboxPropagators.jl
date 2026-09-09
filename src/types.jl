## Description #############################################################################
#
# Definition of types and structures.
#
############################################################################################

export J2PropagatorConstants, J2Propagator, J2OsculatingPropagator
export J4PropagatorConstants, J4Propagator, J4OsculatingPropagator
export TwoBodyPropagator
export MeanElementsFitDivergenceError
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
    struct J2PropagatorConstants{T <: Number}

Constants for the J2 orbit propagator.

# Fields

- `R0::T`: Earth equatorial radius [m].
- `μm::T`: Standard gravitational parameter normalized by the Earth equatorial radius,
    √(GM / R0³) [rad / s].
- `J2::T`: The second gravitational zonal harmonic of the Earth [-].
"""
struct J2PropagatorConstants{T <: Number}
    R0::T
    μm::T
    J2::T
end

"""
    J2PropagatorConstants{T}(
        j2c::J2PropagatorConstants
    ) where {T <: Number} -> J2PropagatorConstants{T}

Create a copy of the J2 propagator constants `j2c` with the element type `T`.
"""
function J2PropagatorConstants{T}(j2c::J2PropagatorConstants) where {T <: Number}
    return J2PropagatorConstants{T}(T(j2c.R0), T(j2c.μm), T(j2c.J2))
end

function Base.convert(
    ::Type{J2PropagatorConstants{T}},
    j2c::J2PropagatorConstants
) where {T <: Number}
    return J2PropagatorConstants{T}(j2c)
end

"""
    mutable struct J2Propagator{Tepoch <: Number, T <: Number}

J2 orbit propagator structure, where `Tepoch` is the type of the epoch and `T` is the type
of the elements and of the internal variables. The fields are set by `j2_init!` and updated
by `j2!`, and must not be modified directly by the callers, except for `j2c`, which must be
assigned before initializing a structure created with the empty constructor.

# Fields

- `orb₀::KeplerianElements{MeanAnomaly, Tepoch, T}`: Initial mean orbit elements [SI units].
- `orbk::KeplerianElements{MeanAnomaly, Tepoch, T}`: Current mean orbit elements [SI units].
- `j2c::J2PropagatorConstants{T}`: Propagator constants.
- `Δt::T`: Timespan from the initial elements' epoch [s], which is `NaN` until the structure
    is initialized.
- `∂Ω::T`: RAAN time derivative [rad / s].
- `∂ω::T`: Argument of perigee time derivative [rad / s].
- `n̄::T`: Perturbed mean motion [rad / s].

# Extended help

## Printing

`show(io, pd)` prints the compact form: the type with its parameters and the epoch of the
initial mean elements as a Julian Day and as a date. `show(io, MIME("text/plain"), pd)`
prints a tree with the sections holding the initial mean elements and their epoch, the
secular rates, the constants, and the last propagation instant, one field per line with its
unit, with the labels in bold and the units dimmed if `io` supports color. The body can be
printed under another header with `SatelliteToolboxBase.print_tree_body`. An
uninitialized structure prints the status `not initialized` instead.
"""
mutable struct J2Propagator{Tepoch <: Number, T <: Number}
    orb₀::KeplerianElements{MeanAnomaly, Tepoch, T}
    orbk::KeplerianElements{MeanAnomaly, Tepoch, T}
    j2c::J2PropagatorConstants{T}
    Δt::T

    # == Auxiliary Variables ===============================================================

    ∂Ω::T
    ∂ω::T
    n̄::T

    # == Constructors ======================================================================

    J2Propagator{Tepoch, T}(args...) where {Tepoch <: Number, T <: Number} = new(args...)
    function J2Propagator{Tepoch, T}() where {Tepoch <: Number, T <: Number}
        pd    = new()
        pd.Δt = _uninitialized_instant(T)
        return pd
    end
end

############################################################################################
#                              J2 Osculating Orbit Propagator                              #
############################################################################################

"""
    mutable struct J2OsculatingPropagator{Tepoch <: Number, T <: Number}

J2 osculating orbit propagator structure, where `Tepoch` is the type of the epoch and `T` is
the type of the elements and of the internal variables. The fields are set by `j2osc_init!`
and updated by `j2osc!`, and must not be modified directly by the callers, except for `j2d`,
which must be assigned before initializing a structure created with the empty constructor.

# Fields

- `j2d::J2Propagator{Tepoch, T}`: J2 orbit propagator that propagates the mean elements.
- `Δt::T`: Timespan from the initial elements' epoch [s], which is `NaN` until the structure
    is initialized.
- `orbk::KeplerianElements{TrueAnomaly, Tepoch, T}`: Current osculating orbit elements [SI
    units].

# Extended help

## Printing

`show(io, pd)` prints the compact form: the type with its parameters and the epoch of the
initial mean elements as a Julian Day and as a date. `show(io, MIME("text/plain"), pd)`
prints a tree with the sections holding the initial mean elements and their epoch, the
secular rates, the constants, and the last propagation instant, one field per line with its
unit, with the labels in bold and the units dimmed if `io` supports color. The body can be
printed under another header with `SatelliteToolboxBase.print_tree_body`. An
uninitialized structure prints the status `not initialized` instead.
"""
mutable struct J2OsculatingPropagator{Tepoch <: Number, T <: Number}
    j2d::J2Propagator{Tepoch, T}
    Δt::T
    orbk::KeplerianElements{TrueAnomaly, Tepoch, T}

    # == Constructors ======================================================================

    J2OsculatingPropagator{Tepoch, T}(args...) where {Tepoch <: Number, T <: Number} =
        new(args...)
    function J2OsculatingPropagator{Tepoch, T}() where {Tepoch <: Number, T <: Number}
        pd    = new()
        pd.Δt = _uninitialized_instant(T)
        return pd
    end
end

############################################################################################
#                                   J4 Orbit Propagator                                    #
############################################################################################

"""
    struct J4PropagatorConstants{T <: Number}

Constants for the J4 orbit propagator.

# Fields

- `R0::T`: Earth equatorial radius [m].
- `μm::T`: Standard gravitational parameter normalized by the Earth equatorial radius,
    √(GM / R0³) [rad / s].
- `J2::T`: The second gravitational zonal harmonic of the Earth [-].
- `J4::T`: The fourth gravitational zonal harmonic of the Earth [-].
"""
struct J4PropagatorConstants{T <: Number}
    R0::T
    μm::T
    J2::T
    J4::T
end

"""
    J4PropagatorConstants{T}(
        j4c::J4PropagatorConstants
    ) where {T <: Number} -> J4PropagatorConstants{T}

Create a copy of the J4 propagator constants `j4c` with the element type `T`.
"""
function J4PropagatorConstants{T}(j4c::J4PropagatorConstants) where {T <: Number}
    return J4PropagatorConstants{T}(T(j4c.R0), T(j4c.μm), T(j4c.J2), T(j4c.J4))
end

function Base.convert(
    ::Type{J4PropagatorConstants{T}}, j4c::J4PropagatorConstants
) where {T <: Number}
    return J4PropagatorConstants{T}(j4c)
end

"""
    mutable struct J4Propagator{Tepoch <: Number, T <: Number}

J4 orbit propagator structure, where `Tepoch` is the type of the epoch and `T` is the type
of the elements and of the internal variables. The fields are set by `j4_init!` and updated
by `j4!`, and must not be modified directly by the callers, except for `j4c`, which must be
assigned before initializing a structure created with the empty constructor.

# Fields

- `orb₀::KeplerianElements{MeanAnomaly, Tepoch, T}`: Initial mean orbit elements [SI units].
- `orbk::KeplerianElements{MeanAnomaly, Tepoch, T}`: Current mean orbit elements [SI units].
- `j4c::J4PropagatorConstants{T}`: Propagator constants.
- `Δt::T`: Timespan from the initial elements' epoch [s], which is `NaN` until the structure
    is initialized.
- `∂Ω::T`: RAAN time derivative [rad / s].
- `∂ω::T`: Argument of perigee time derivative [rad / s].
- `n̄::T`: Perturbed mean motion [rad / s].

# Extended help

## Printing

`show(io, pd)` prints the compact form: the type with its parameters and the epoch of the
initial mean elements as a Julian Day and as a date. `show(io, MIME("text/plain"), pd)`
prints a tree with the sections holding the initial mean elements and their epoch, the
secular rates, the constants, and the last propagation instant, one field per line with its
unit, with the labels in bold and the units dimmed if `io` supports color. The body can be
printed under another header with `SatelliteToolboxBase.print_tree_body`. An
uninitialized structure prints the status `not initialized` instead.
"""
mutable struct J4Propagator{Tepoch <: Number, T <: Number}
    orb₀::KeplerianElements{MeanAnomaly, Tepoch, T}
    orbk::KeplerianElements{MeanAnomaly, Tepoch, T}
    j4c::J4PropagatorConstants{T}
    Δt::T

    # == Auxiliary Variables ===============================================================

    ∂Ω::T
    ∂ω::T
    n̄::T

    # == Constructors ======================================================================

    J4Propagator{Tepoch, T}(args...) where {Tepoch <: Number, T <: Number} = new(args...)
    function J4Propagator{Tepoch, T}() where {Tepoch <: Number, T <: Number}
        pd    = new()
        pd.Δt = _uninitialized_instant(T)
        return pd
    end
end

############################################################################################
#                              J4 Osculating Orbit Propagator                              #
############################################################################################

"""
    mutable struct J4OsculatingPropagator{Tepoch <: Number, T <: Number}

J4 osculating orbit propagator structure, where `Tepoch` is the type of the epoch and `T` is
the type of the elements and of the internal variables. The fields are set by `j4osc_init!`
and updated by `j4osc!`, and must not be modified directly by the callers, except for `j4d`,
which must be assigned before initializing a structure created with the empty constructor.

# Fields

- `j4d::J4Propagator{Tepoch, T}`: J4 orbit propagator that propagates the mean elements.
- `Δt::T`: Timespan from the initial elements' epoch [s], which is `NaN` until the structure
    is initialized.
- `orbk::KeplerianElements{TrueAnomaly, Tepoch, T}`: Current osculating orbit elements [SI
    units].

# Extended help

## Printing

`show(io, pd)` prints the compact form: the type with its parameters and the epoch of the
initial mean elements as a Julian Day and as a date. `show(io, MIME("text/plain"), pd)`
prints a tree with the sections holding the initial mean elements and their epoch, the
secular rates, the constants, and the last propagation instant, one field per line with its
unit, with the labels in bold and the units dimmed if `io` supports color. The body can be
printed under another header with `SatelliteToolboxBase.print_tree_body`. An
uninitialized structure prints the status `not initialized` instead.
"""
mutable struct J4OsculatingPropagator{Tepoch <: Number, T <: Number}
    j4d::J4Propagator{Tepoch, T}
    Δt::T
    orbk::KeplerianElements{TrueAnomaly, Tepoch, T}

    # == Constructors ======================================================================

    J4OsculatingPropagator{Tepoch, T}(args...) where {Tepoch <: Number, T <: Number} =
        new(args...)
    function J4OsculatingPropagator{Tepoch, T}() where {Tepoch <: Number, T <: Number}
        pd    = new()
        pd.Δt = _uninitialized_instant(T)
        return pd
    end
end

############################################################################################
#                                   Two-Body Propagator                                    #
############################################################################################

"""
    mutable struct TwoBodyPropagator{Tepoch <: Number, T <: Number}

Two-body orbit propagator structure, where `Tepoch` is the type of the epoch and `T` is the
type of the elements and of the internal variables. The fields are set by `twobody_init!`
and updated by `twobody!`, and must not be modified directly by the callers, except for `μ`,
which must be assigned before initializing a structure created with the empty constructor.

# Fields

- `orb₀::KeplerianElements{MeanAnomaly, Tepoch, T}`: Initial mean orbit elements [SI units].
- `orbk::KeplerianElements{MeanAnomaly, Tepoch, T}`: Current mean orbit elements [SI units].
- `μ::T`: Standard gravitational parameter of the central body [m³ / s²].
- `Δt::T`: Timespan from the initial elements' epoch [s], which is `NaN` until the structure
    is initialized.
- `n₀::T`: Mean motion [rad / s].

# Extended help

## Printing

`show(io, pd)` prints the compact form: the type with its parameters and the epoch of the
initial mean elements as a Julian Day and as a date. `show(io, MIME("text/plain"), pd)`
prints a tree with the sections holding the initial mean elements and their epoch, the
mean motion, the constants, and the last propagation instant, one field per line with its
unit, with the labels in bold and the units dimmed if `io` supports color. The body can be
printed under another header with `SatelliteToolboxBase.print_tree_body`. An
uninitialized structure prints the status `not initialized` instead.
"""
mutable struct TwoBodyPropagator{Tepoch <: Number, T <: Number}
    orb₀::KeplerianElements{MeanAnomaly, Tepoch, T}
    orbk::KeplerianElements{MeanAnomaly, Tepoch, T}
    μ::T
    Δt::T

    # == Auxiliary Variables ===============================================================

    n₀::T

    # == Constructors ======================================================================

    TwoBodyPropagator{Tepoch, T}(args...) where {Tepoch <: Number, T <: Number} =
        new(args...)
    function TwoBodyPropagator{Tepoch, T}() where {Tepoch <: Number, T <: Number}
        pd    = new()
        pd.Δt = _uninitialized_instant(T)
        return pd
    end
end

############################################################################################
#                                        Exceptions                                        #
############################################################################################

"""
    struct MeanElementsFitDivergenceError <: Exception

Exception thrown when the least-square iterations used to fit the mean elements of the
analytical propagators diverge.

# Fields

- `iteration::Int`: Iteration in which the divergence was detected.
- `residue::Float64`: Total RMSE of the residue in that iteration.
"""
struct MeanElementsFitDivergenceError <: Exception
    iteration::Int
    residue::Float64
end

function Base.showerror(io::IO, e::MeanElementsFitDivergenceError)
    print(
        io,
        "MeanElementsFitDivergenceError: The least-square iterations diverged at ",
        "iteration ",
        e.iteration,
        " with a total RMSE of ",
        e.residue,
        ".",
    )
    return nothing
end

############################################################################################
#                                           API                                            #
############################################################################################

# == J2 Orbit Propagator ===================================================================

"""
    struct OrbitPropagatorJ2{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}

J2 orbit propagator for the `Propagators` API.

# Fields

- `j2d::J2Propagator{Tepoch, T}`: Structure that stores the J2 orbit propagator data (see
    [`J2Propagator`](@ref)).

# Extended help

## Printing

`show(io, orbp)` prints the compact form: the type with its parameters, the propagator
name, and the epoch as a Julian Day and as a date. `show(io, MIME("text/plain"), orbp)`
prints the same header followed by the body of the rich representation of `j2d`.
"""
struct OrbitPropagatorJ2{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}
    j2d::J2Propagator{Tepoch, T}
end

# == J2 Osculating Orbit Propagator ========================================================

"""
    struct OrbitPropagatorJ2Osculating{
        Tepoch <: Number,
        T <: Number
    } <: OrbitPropagator{Tepoch, T}

J2 osculating orbit propagator for the `Propagators` API.

# Fields

- `j2oscd::J2OsculatingPropagator{Tepoch, T}`: Structure that stores the J2 osculating orbit
    propagator data (see [`J2OsculatingPropagator`](@ref)).

# Extended help

## Printing

`show(io, orbp)` prints the compact form: the type with its parameters, the propagator
name, and the epoch as a Julian Day and as a date. `show(io, MIME("text/plain"), orbp)`
prints the same header followed by the body of the rich representation of `j2oscd`.
"""
struct OrbitPropagatorJ2Osculating{Tepoch <: Number, T <: Number} <:
       OrbitPropagator{Tepoch, T}
    j2oscd::J2OsculatingPropagator{Tepoch, T}
end

# == J4 Orbit Propagator ===================================================================

"""
    struct OrbitPropagatorJ4{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}

J4 orbit propagator for the `Propagators` API.

# Fields

- `j4d::J4Propagator{Tepoch, T}`: Structure that stores the J4 orbit propagator data (see
    [`J4Propagator`](@ref)).

# Extended help

## Printing

`show(io, orbp)` prints the compact form: the type with its parameters, the propagator
name, and the epoch as a Julian Day and as a date. `show(io, MIME("text/plain"), orbp)`
prints the same header followed by the body of the rich representation of `j4d`.
"""
struct OrbitPropagatorJ4{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}
    j4d::J4Propagator{Tepoch, T}
end

# == J4 Osculating Orbit Propagator ========================================================

"""
    struct OrbitPropagatorJ4Osculating{
        Tepoch <: Number,
        T <: Number
    } <: OrbitPropagator{Tepoch, T}

J4 osculating orbit propagator for the `Propagators` API.

# Fields

- `j4oscd::J4OsculatingPropagator{Tepoch, T}`: Structure that stores the J4 osculating orbit
    propagator data (see [`J4OsculatingPropagator`](@ref)).

# Extended help

## Printing

`show(io, orbp)` prints the compact form: the type with its parameters, the propagator
name, and the epoch as a Julian Day and as a date. `show(io, MIME("text/plain"), orbp)`
prints the same header followed by the body of the rich representation of `j4oscd`.
"""
struct OrbitPropagatorJ4Osculating{Tepoch <: Number, T <: Number} <:
       OrbitPropagator{Tepoch, T}
    j4oscd::J4OsculatingPropagator{Tepoch, T}
end

# == SGP4 Orbit Propagator =================================================================

"""
    struct OrbitPropagatorSgp4{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}

SGP4 orbit propagator for the `Propagators` API.

# Fields

- `sgp4d::Sgp4Propagator{Tepoch, T}`: Structure that stores the SGP4 orbit propagator data
    (see `Sgp4Propagator` in **SatelliteToolboxSgp4.jl**).

# Extended help

## Printing

`show(io, orbp)` prints the compact form: the type with its parameters, the propagator
name, and the epoch as a Julian Day and as a date. `show(io, MIME("text/plain"), orbp)`
prints the same header followed by the body of the rich representation of `sgp4d`.
"""
struct OrbitPropagatorSgp4{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}
    sgp4d::Sgp4Propagator{Tepoch, T}
end

# == Two-Body Orbit Propagator =============================================================

"""
    struct OrbitPropagatorTwoBody{
        Tepoch <: Number,
        T <: Number
    } <: OrbitPropagator{Tepoch, T}

Two-body orbit propagator for the `Propagators` API.

# Fields

- `tbd::TwoBodyPropagator{Tepoch, T}`: Structure that stores the two-body orbit propagator
    data (see [`TwoBodyPropagator`](@ref)).

# Extended help

## Printing

`show(io, orbp)` prints the compact form: the type with its parameters, the propagator
name, and the epoch as a Julian Day and as a date. `show(io, MIME("text/plain"), orbp)`
prints the same header followed by the body of the rich representation of `tbd`.
"""
struct OrbitPropagatorTwoBody{Tepoch <: Number, T <: Number} <: OrbitPropagator{Tepoch, T}
    tbd::TwoBodyPropagator{Tepoch, T}
end

############################################################################################
#                                        Julia API                                         #
############################################################################################

const PropagatorData{Tepoch, T} = Union{
    J2Propagator{Tepoch, T},
    J2OsculatingPropagator{Tepoch, T},
    J4Propagator{Tepoch, T},
    J4OsculatingPropagator{Tepoch, T},
    TwoBodyPropagator{Tepoch, T},
}

function Base.copy(pd::P) where {P <: PropagatorData}
    return P(ntuple(i -> _copy_field(getfield(pd, i)), Val(fieldcount(P)))...)
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _uninitialized_instant(::Type{T}) where {T <: Number} -> T

Return the propagation instant stored in the field `Δt` by the empty constructors of the
propagator structures, `NaN` converted to `T`, which marks a structure as not initialized
until an initialization function assigns its fields. Hence, `T` must be able to represent
`NaN`, as every floating-point type and the dual numbers do.
"""
_uninitialized_instant(::Type{T}) where {T <: Number} = T(NaN)

"""
    _copy_field(x) -> typeof(x)

Return the value stored in a propagator field when copying the structure. Nested propagator
structures are copied, whereas every other field is immutable and is returned as is.
"""
_copy_field(x) = x
_copy_field(pd::PropagatorData) = copy(pd)

"""
    _propagator_eltype(
        ::Type{Tconstants},
        ::Type{Tkepler}
    ) where {Tconstants <: Number, Tkepler <: Number} -> Type

Return the element type used by a propagator initialized with constants of type
`Tconstants` and Keplerian elements of type `Tkepler`. If `Tkepler` is an `AbstractFloat`,
the constants type is used, so that, e.g., `Float32` constants select a `Float32`
propagation regardless of the elements type. Otherwise, the types are promoted, which keeps
the propagation differentiable when the elements are `ForwardDiff.Dual` numbers.
"""
function _propagator_eltype(
    ::Type{Tconstants}, ::Type{Tkepler}
) where {Tconstants <: Number, Tkepler <: AbstractFloat}
    return Tconstants
end

function _propagator_eltype(
    ::Type{Tconstants}, ::Type{Tkepler}
) where {Tconstants <: Number, Tkepler <: Number}
    return promote_type(Tconstants, Tkepler)
end
