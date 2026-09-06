## Description #############################################################################
#
#   J4 osculating orbit propagator algorithm.
#
#   This algorithm propagates the orbit considering the secular perturbations from
#   the terms J2, J2², and J4 and the short-period perturbations from the J2
#   gravitational term only. The algorithm is based on Kwok version as indicated in
#   [1, p. 708-710].
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorn, CA, USA.
#
############################################################################################

export j4osc_init, j4osc_init!, j4osc, j4osc!
export fit_j4osc_mean_elements, fit_j4osc_mean_elements!
export update_j4osc_mean_elements_epoch, update_j4osc_mean_elements_epoch!

############################################################################################
#                                        Functions                                         #
############################################################################################

"""
    j4osc_init(orb₀::KeplerianElements; kwargs...) -> J4OsculatingPropagator

Create and initialize the J4 osculating orbit propagator structure using the mean Keplerian
elements `orb₀`.

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j4c`.

# Keywords

- `j4c::J4PropagatorConstants`: J4 orbit propagator constants (see
    [`J4PropagatorConstants`](@ref)).
    (**Default**: `J4C_EGM2008`)
"""
function j4osc_init(
    orb₀::KeplerianElements{Tanomaly, Tepoch, Tkepler};
    j4c::J4PropagatorConstants{Tj4c} = J4C_EGM2008,
) where {Tanomaly <: AbstractAnomaly, Tepoch <: Number, Tkepler <: Number, Tj4c <: Number}
    T = _propagator_eltype(Tj4c, Tkepler)

    # Allocate the J4 propagator structure that will propagate the mean elements and assign
    # the constants, which are used in the initialization.
    j4d = J4Propagator{Tepoch, T}()
    j4d.j4c = convert(J4PropagatorConstants{T}, j4c)

    # Allocate the J4 osculating propagator structure.
    j4oscd = J4OsculatingPropagator{Tepoch, T}()
    j4oscd.j4d = j4d

    # Initialize the propagator and return.
    j4osc_init!(j4oscd, orb₀)

    return j4oscd
end

"""
    j4osc_init!(j4oscd::J4OsculatingPropagator, orb₀::KeplerianElements) -> Nothing

Initialize the J4 osculating orbit propagator structure `j4oscd` using the mean Keplerian
elements `orb₀`.

!!! warning

    The propagation constants `j4c::J4PropagatorConstants` in `j4oscd.j4d` will not be
    changed. Hence, they must be initialized.
"""
function j4osc_init!(j4oscd::J4OsculatingPropagator, orb₀::KeplerianElements)
    _j4osc_init!(j4oscd, orb₀)

    # Call the propagation one time to update the osculating elements.
    j4osc!(j4oscd, 0)

    return nothing
end

"""
    j4osc(
        Δt::Number,
        orb₀::KeplerianElements;
        kwargs...
    ) -> SVector{3, T}, SVector{3, T}, J4OsculatingPropagator

Initialize the J4 osculating propagator structure using the input elements `orb₀` and
propagate the orbit until the time Δt [s].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j4c`.

# Keywords

- `j4c::J4PropagatorConstants{T}`: J4 orbit propagator constants (see
  [`J4PropagatorConstants`](@ref)).
    (**Default**: `J4C_EGM2008`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`J4OsculatingPropagator`](@ref): Structure with the initialized propagator.

# Remarks

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j4osc(
    Δt::Number, orb₀::KeplerianElements; j4c::J4PropagatorConstants = J4C_EGM2008
)
    j4oscd = j4osc_init(orb₀; j4c = j4c)
    r_i, v_i = j4osc!(j4oscd, Δt)
    return r_i, v_i, j4oscd
end

"""
    j4osc!(
        j4oscd::J4OsculatingPropagator{Tepoch, T},
        t::Number
    ) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit defined in `j4oscd` (see [`J4OsculatingPropagator`](@ref)) to `t` [s]
after the epoch of the input mean elements in `j4oscd`.

!!! note

    The internal values in `j4oscd` will be modified.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

# Remarks

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j4osc!(
    j4oscd::J4OsculatingPropagator{Tepoch, T}, t::Number
) where {Tepoch <: Number, T <: Number}
    # First, we need to propagate the mean elements since they are necessary to compute the
    # short-periodic perturbations.
    j4d = j4oscd.j4d
    mean_orbk = _j4_mean_elements!(j4d, t)

    # Unpack the propagator constants.
    j4c = j4d.j4c
    R₀  = j4c.R0
    μm  = j4c.μm
    J₂  = j4c.J2

    orbk = _osculating_elements(mean_orbk, R₀, μm, J₂)

    # Compute the position and velocity considering the osculating elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the J4 orbit propagator structure.
    j4oscd.Δt   = T(t)
    j4oscd.orbk = orbk

    return r_i_k, v_i_k
end

"""
    fit_j4osc_mean_elements(
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {
        Tjd <: Number,
        Tv <: AbstractVector
    } -> KeplerianElements{MeanAnomaly, Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the J4 osculating orbit propagator using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    This algorithm version will allocate a new J4 osculating propagator with the default
    constants `J4C_EGM2008`. If another set of constants are required, use the function
    [`fit_j4osc_mean_elements!`](@ref) instead.

# Keywords

- `atol::Number`: Tolerance for the residual absolute value. If the residual is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default**: 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residuals. If the
    relative difference between the residuals in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default**: 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default**: nothing)
- `jacobian_method::Union{FiniteDiffJacobian, ForwardDiffJacobian}`: Method used to compute
    the Jacobian matrix. Use `FiniteDiffJacobian()` for finite differences or
    `ForwardDiffJacobian()` for `ForwardDiff.jl` automatic differentiation.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. Only used with
    `FiniteDiffJacobian()`.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`. Only used with `FiniteDiffJacobian()`.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default**: 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default**: vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default**: `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{MeanAnomaly, Float64, Float64}`: Fitted Keplerian elements.
- `SMatrix{6, 6, Float64}`: Final covariance matrix of the least-square algorithm.

# Examples

```julia-repl
julia> vr_i = [
           [-6792.402703741442, 2192.6458461287293, 0.18851758695295118] .* 1000,
           [-6357.88873265975, 2391.9476768911686, 2181.838771262736] .* 1000
       ];

julia> vv_i = [
           [0.3445760107690598, 1.0395135806993514, 7.393686131436984] .* 1000,
           [2.5285015912807003, 0.27812476784300005, 7.030323100703928] .* 1000
       ];

julia> vjd = [
           2.46002818657856e6,
           2.460028190050782e6
       ];

julia> orb, P = fit_j4osc_mean_elements(vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J4 osculating propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          4          1.68763e-05           0.00259611              2.59616          1.08586e-09 %

(KeplerianElements{MeanAnomaly, Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.9999427681054915 -0.004901179110064787 … -0.00012021347573371965 -0.00011550159594742305; -0.004901179104382454 1.0007325998821814 … 0.0032661543927207434 3.8582051095353965e-5; … ; -0.00012021347571372953 0.0032661543927214936 … 2.204822799517442e-5 5.3982379473473e-8; -0.00011550159594878459 3.8582051094144375e-5 … 5.3982379468632544e-8 2.1657607939466597e-5])

julia> orb
KeplerianElements{MeanAnomaly, Float64, Float64}:
  Epoch              : 2.46003e6 (2023-03-24T16:33:40.388)
  Semi-Major Axis    : 7135.792461 km
  Eccentricity       : 0.001353135365
  Inclination        : 98.4304116°
  RA of Asc. Node    : 162.1131631°
  Arg. of Pericenter : 64.96868276°
  Mean Anomaly       : 313.1552992°
```
"""
function fit_j4osc_mean_elements(
    vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...
) where {Tjd <: Number, Tv <: AbstractVector}
    # Allocate the J4 propagator structure that will propagate the mean elements.
    j4d = J4Propagator{Float64, Float64}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = J4C_EGM2008

    # Allocate the J4 osculating propagator structure.
    j4oscd = J4OsculatingPropagator{Float64, Float64}()
    j4oscd.j4d = j4d

    return fit_j4osc_mean_elements!(j4oscd, vjd, vr_i, vv_i; kwargs...)
end

"""
    fit_j4osc_mean_elements!(
        j4oscd::J4OsculatingPropagator{Tepoch, T},
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {
        T <: Number,
        Tepoch <: Number,
        Tjd <: Number,
        Tv <: AbstractVector
    } -> KeplerianElements{MeanAnomaly, Tepoch, T}, SMatrix{6, 6, T}

Fit a set of mean Keplerian elements for the J4 osculating orbit propagator `j4oscd` using
the osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    The J4 osculating orbit propagator `j4oscd` will be initialized with the Keplerian
    elements returned by the function.

# Keywords

- `atol::Number`: Tolerance for the residual absolute value. If the residual is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default**: 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residuals. If the
    relative difference between the residuals in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default**: 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements
    fitting process. If it is `nothing`, the algorithm will obtain an initial estimate from
    the osculating elements in `vr_i` and `vv_i`.
    (**Default**: nothing)
- `jacobian_method::Union{FiniteDiffJacobian, ForwardDiffJacobian}`: Method used to compute
    the Jacobian matrix. Use `FiniteDiffJacobian()` for finite differences or
    `ForwardDiffJacobian()` for `ForwardDiff.jl` automatic differentiation.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix. Only used with
    `FiniteDiffJacobian()`.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until its absolute value is higher than
    `jacobian_perturbation_tol`. Only used with `FiniteDiffJacobian()`.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default**: 50)
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements.
    (**Default**: vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default**: `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{MeanAnomaly, Tepoch, T}`: Fitted Keplerian elements.
- `SMatrix{6, 6, T}`: Final covariance matrix of the least-square algorithm.

# Examples

```julia-repl
# Allocate a new J4 osculating orbit propagator using a dummy set of Keplerian elements.
julia> j4oscd = j4osc_init(KeplerianElements(0.0, 7000e3, 0, 0, 0, 0, 0));

julia> vr_i = [
           [-6792.402703741442, 2192.6458461287293, 0.18851758695295118] .* 1000,
           [-6357.88873265975, 2391.9476768911686, 2181.838771262736] .* 1000
       ];

julia> vv_i = [
           [0.3445760107690598, 1.0395135806993514, 7.393686131436984] .* 1000,
           [2.5285015912807003, 0.27812476784300005, 7.030323100703928] .* 1000
       ];

julia> vjd = [
           2.46002818657856e6,
           2.460028190050782e6
       ];

julia> orb, P = fit_j4osc_mean_elements!(j4oscd, vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J4 osculating propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          4          1.68763e-05           0.00259611              2.59616          1.08586e-09 %

(KeplerianElements{MeanAnomaly, Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.9999427681054915 -0.004901179110064787 … -0.00012021347573371965 -0.00011550159594742305; -0.004901179104382454 1.0007325998821814 … 0.0032661543927207434 3.8582051095353965e-5; … ; -0.00012021347571372953 0.0032661543927214936 … 2.204822799517442e-5 5.3982379473473e-8; -0.00011550159594878459 3.8582051094144375e-5 … 5.3982379468632544e-8 2.1657607939466597e-5])

julia> orb
KeplerianElements{MeanAnomaly, Float64, Float64}:
  Epoch              : 2.46003e6 (2023-03-24T16:33:40.388)
  Semi-Major Axis    : 7135.792461 km
  Eccentricity       : 0.001353135365
  Inclination        : 98.4304116°
  RA of Asc. Node    : 162.1131631°
  Arg. of Pericenter : 64.96868276°
  Mean Anomaly       : 313.1552992°
```
"""
function fit_j4osc_mean_elements!(
    j4oscd::J4OsculatingPropagator{Tepoch, T},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tepoch <: Number, T <: Number, Tjd <: Number, Tv <: AbstractVector}
    return _fit_mean_elements!(j4oscd, vjd, vr_i, vv_i; kwargs...)
end

"""
    update_j4osc_mean_elements_epoch(
        orb::KeplerianElements,
        new_epoch::Union{Number, DateTime}
    ) -> KeplerianElements{MeanAnomaly}

Update the epoch of the mean elements `orb` using a J4 osculating orbit propagator to
`new_epoch`, which can be represented by a Julian Day or a `DateTime`.

!!! note

    This algorithm version will allocate a new J4 osculating propagator with the default
    constants `J4C_EGM2008`. If another set of constants are required, use the function
    [`update_j4osc_mean_elements_epoch!`](@ref) instead.

# Examples

```julia-repl
julia> orb = KeplerianElements(
           DateTime("2023-01-01") |> datetime2julian,
           7190.982e3,
           0.001111,
           98.405 |> deg2rad,
           90     |> deg2rad,
           200    |> deg2rad,
           45     |> deg2rad
       )
KeplerianElements{TrueAnomaly, Float64, Float64}:
  Epoch              : 2.45995e6 (2023-01-01T00:00:00)
  Semi-Major Axis    : 7190.982 km
  Eccentricity       : 0.001111
  Inclination        : 98.405°
  RA of Asc. Node    : 90.0°
  Arg. of Pericenter : 200.0°
  True Anomaly       : 45.0°

julia> update_j4osc_mean_elements_epoch(orb, DateTime("2023-01-02"))
KeplerianElements{MeanAnomaly, Float64, Float64}:
  Epoch              : 2.45995e6 (2023-01-02T00:00:00)
  Semi-Major Axis    : 7190.982 km
  Eccentricity       : 0.001111
  Inclination        : 98.405°
  RA of Asc. Node    : 90.95551368°
  Arg. of Pericenter : 197.0785362°
  Mean Anomaly       : 127.191195°
```
"""
function update_j4osc_mean_elements_epoch(
    orb::KeplerianElements{Tanomaly, Tepoch, T}, new_epoch::Union{Number, DateTime}
) where {Tanomaly <: AbstractAnomaly, Tepoch <: Number, T <: Number}
    # Allocate the J4 propagator structure that will propagate the mean elements.
    j4d = J4Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j4d.j4c = J4C_EGM2008

    # Allocate the J4 osculating propagator structure.
    j4oscd = J4OsculatingPropagator{Tepoch, T}()
    j4oscd.j4d = j4d

    return update_j4osc_mean_elements_epoch!(j4oscd, orb, new_epoch)
end

"""
    update_j4osc_mean_elements_epoch!(
        j4oscd::J4OsculatingPropagator,
        orb::KeplerianElements,
        new_epoch::Union{Number, DateTime}
    ) -> KeplerianElements{MeanAnomaly}

Update the epoch of the mean elements `orb` using the propagator `j4oscd` to `new_epoch`,
which can be represented by a Julian Day or a `DateTime`.

!!! note

    The J4 osculating orbit propagator `j4oscd` will be initialized with the Keplerian
    elements returned by the function.

# Examples

```julia-repl
julia> orb = KeplerianElements(
           DateTime("2023-01-01") |> datetime2julian,
           7190.982e3,
           0.001111,
           98.405 |> deg2rad,
           90     |> deg2rad,
           200    |> deg2rad,
           45     |> deg2rad
       )
KeplerianElements{TrueAnomaly, Float64, Float64}:
  Epoch              : 2.45995e6 (2023-01-01T00:00:00)
  Semi-Major Axis    : 7190.982 km
  Eccentricity       : 0.001111
  Inclination        : 98.405°
  RA of Asc. Node    : 90.0°
  Arg. of Pericenter : 200.0°
  True Anomaly       : 45.0°

# Allocate a new J4 osculating orbit propagator using the created Keplerian elements. Notice
# that any set of Keplerian elements can be used here.
julia> j4oscd = j4osc_init(orb);

julia> update_j4osc_mean_elements_epoch!(j4oscd, orb, DateTime("2023-01-02"))
KeplerianElements{MeanAnomaly, Float64, Float64}:
  Epoch              : 2.45995e6 (2023-01-02T00:00:00)
  Semi-Major Axis    : 7190.982 km
  Eccentricity       : 0.001111
  Inclination        : 98.405°
  RA of Asc. Node    : 90.95551368°
  Arg. of Pericenter : 197.0785362°
  Mean Anomaly       : 127.191195°
```
"""
function update_j4osc_mean_elements_epoch!(
    j4oscd::J4OsculatingPropagator, orb::KeplerianElements, new_epoch::DateTime
)
    dt = datetime2julian(new_epoch)
    return update_j4osc_mean_elements_epoch!(j4oscd, orb, dt)
end

function update_j4osc_mean_elements_epoch!(
    j4oscd::J4OsculatingPropagator, orb::KeplerianElements, new_epoch::Number
)
    orb = _update_mean_elements_epoch!(j4oscd, orb, new_epoch)

    # The generic algorithm initializes the J4 propagator only. Hence, we must compute the
    # osculating elements at the new epoch to complete the initialization.
    j4osc!(j4oscd, 0)

    return orb
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _similar_propagator(
        j4oscd::J4OsculatingPropagator{Tepoch},
        ::Type{T}
    ) where {Tepoch <: Number, T <: Number} -> J4OsculatingPropagator{Tepoch, T}

Create an uninitialized J4 osculating propagator with the same constants as `j4oscd`
converted to the element type `T`.
"""
function _similar_propagator(
    j4oscd::J4OsculatingPropagator{Tepoch}, ::Type{T}
) where {Tepoch <: Number, T <: Number}
    new_j4oscd = J4OsculatingPropagator{Tepoch, T}()
    new_j4oscd.j4d = _similar_propagator(j4oscd.j4d, T)

    return new_j4oscd
end

"""
    _j4osc_init!(j4oscd::J4OsculatingPropagator, orb₀::KeplerianElements) -> Nothing

Initialize the J4 osculating orbit propagator `j4oscd` with the mean Keplerian elements
`orb₀` [SI units] without computing the osculating elements at the initial instant. The
callers that propagate the orbit right afterwards use this function, since the osculating
elements would be overwritten anyway.
"""
function _j4osc_init!(j4oscd::J4OsculatingPropagator, orb₀::KeplerianElements)
    # Initialize the J4 propagator that will propagate the mean elements.
    j4_init!(j4oscd.j4d, orb₀)

    return nothing
end
