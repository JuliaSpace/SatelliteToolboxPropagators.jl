## Description #############################################################################
#
# J2 osculating orbit propagator algorithm.
#
# This algorithm propagates the orbit considering the secular and short-period
# perturbations introduced by the J2 gravitational term. The algorithm is based on Kwok
# version as indicated in [1, p. 708-710].
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#     Press, Hawthorn, CA, USA.
#
############################################################################################

export j2osc_init, j2osc_init!, j2osc, j2osc!
export fit_j2osc_mean_elements, fit_j2osc_mean_elements!
export update_j2osc_mean_elements_epoch, update_j2osc_mean_elements_epoch!

############################################################################################
#                                        Functions                                         #
############################################################################################

"""
    j2osc_init(orb₀::KeplerianElements; kwargs...) -> J2OsculatingPropagator

Create and initialize the J2 osculating orbit propagator structure using the mean Keplerian
elements `orb₀` [SI units].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.

# Keywords

- `j2c::J2PropagatorConstants`: J2 orbit propagator constants (see
    [`J2PropagatorConstants`](@ref)).
    (**Default**: `J2C_EGM2008`)
"""
function j2osc_init(
    orb₀::KeplerianElements{Tanomaly, Tepoch, Tkepler};
    j2c::J2PropagatorConstants{Tj2c} = J2C_EGM2008,
) where {Tanomaly <: AbstractAnomaly, Tepoch <: Number, Tkepler <: Number, Tj2c <: Number}
    T = _propagator_eltype(Tj2c, Tkepler)

    # Allocate the J2 propagator structure that will propagate the mean elements and assign
    # the constants, which are used in the initialization.
    j2d = J2Propagator{Tepoch, T}()
    j2d.j2c = convert(J2PropagatorConstants{T}, j2c)

    # Allocate the J2 osculating propagator structure.
    j2oscd = J2OsculatingPropagator{Tepoch, T}()
    j2oscd.j2d = j2d

    # Initialize the propagator and return.
    j2osc_init!(j2oscd, orb₀)

    return j2oscd
end

"""
    j2osc_init!(j2oscd::J2OsculatingPropagator, orb₀::KeplerianElements) -> Nothing

Initialize the J2 osculating orbit propagator structure `j2oscd` using the mean Keplerian
elements `orb₀` [SI units].

!!! warning

    The propagation constants `j2c::J2PropagatorConstants` in `j2oscd.j2d` will not be
    changed. Hence, they must be initialized.
"""
function j2osc_init!(j2oscd::J2OsculatingPropagator, orb₀::KeplerianElements)
    _j2osc_init!(j2oscd, orb₀)

    # Call the propagation one time to update the osculating elements.
    j2osc!(j2oscd, 0)

    return nothing
end

"""
    j2osc(
        Δt::Number,
        orb₀::KeplerianElements;
        kwargs...
    ) -> SVector{3, T}, SVector{3, T}, J2OsculatingPropagator

Initialize the J2 osculating propagator structure using the input elements `orb₀` [SI units]
and propagate the orbit until the time Δt [s].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `j2c`.

# Keywords

- `j2c::J2PropagatorConstants{T}`: J2 orbit propagator constants (see
    [`J2PropagatorConstants`](@ref)).
    (**Default**: `J2C_EGM2008`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`J2OsculatingPropagator`](@ref): Structure with the initialized propagator.

# Remarks

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters. Notice that the perturbation theory requires an inertial
frame with true equator.
"""
function j2osc(
    Δt::Number, orb₀::KeplerianElements; j2c::J2PropagatorConstants = J2C_EGM2008
)
    j2oscd = j2osc_init(orb₀; j2c = j2c)
    r_i, v_i = j2osc!(j2oscd, Δt)
    return r_i, v_i, j2oscd
end

"""
    j2osc!(
        j2oscd::J2OsculatingPropagator{Tepoch, T},
        t::Number
    ) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit defined in `j2oscd` (see [`J2OsculatingPropagator`](@ref)) to `t` [s]
after the epoch of the input mean elements in `j2oscd`.

!!! note

    The internal values in `j2oscd` will be modified.

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
function j2osc!(
    j2oscd::J2OsculatingPropagator{Tepoch, T}, t::Number
) where {Tepoch <: Number, T <: Number}
    # First, we need to propagate the mean elements since they are necessary to compute the
    # short-periodic perturbations.
    j2d = j2oscd.j2d
    mean_orbk = _j2_mean_elements!(j2d, t)

    # Unpack the propagator constants.
    j2c = j2d.j2c
    R₀  = j2c.R0
    μm  = j2c.μm
    J₂  = j2c.J2

    orbk = _osculating_elements(mean_orbk, R₀, μm, J₂)

    # Compute the position and velocity considering the osculating elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the J2 orbit propagator structure.
    j2oscd.Δt   = T(t)
    j2oscd.orbk = orbk

    return r_i_k, v_i_k
end

"""
    fit_j2osc_mean_elements(
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {
        Tjd <: Number,
        Tv <: AbstractVector
    } -> KeplerianElements{MeanAnomaly, Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the J2 osculating orbit propagator using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    This algorithm version will allocate a new J2 osculating propagator with the default
    constants `J2C_EGM2008`. If another set of constants are required, use the function
    [`fit_j2osc_mean_elements!`](@ref) instead.

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

julia> orb, P = fit_j2osc_mean_elements(vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J2 osculating propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          4          1.69072e-05           0.00260088              2.60093         -2.38817e-10 %

(KeplerianElements{MeanAnomaly, Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.9999427807079662 -0.004901152245568996 … -0.00012021282775592569 -0.00011550583198160739; -0.004901152248513471 1.000732580069942 … 0.003266156907020634 3.8567169210589895e-5; … ; -0.00012021282776408784 0.0032661569070203842 … 2.2048267812548562e-5 5.3827666735694705e-8; -0.00011550583198608087 3.856716921258199e-5 … 5.382766674268047e-8 2.165742075903671e-5])

julia> orb
KeplerianElements{MeanAnomaly, Float64, Float64}:
              Epoch :    2.46003e6 (2023-03-24T16:33:40.388)
    Semi-major axis : 7135.8        km
       Eccentricity :    0.00135383
        Inclination :   98.4304     °
               RAAN :  162.113      °
  Arg. of Periapsis :   64.9256     °
       Mean Anomaly :  313.198      °
```
"""
function fit_j2osc_mean_elements(
    vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...
) where {Tjd <: Number, Tv <: AbstractVector}
    # Allocate the J2 propagator structure that will propagate the mean elements.
    j2d = J2Propagator{Float64, Float64}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = J2C_EGM2008

    # Allocate the J2 osculating propagator structure.
    j2oscd = J2OsculatingPropagator{Float64, Float64}()
    j2oscd.j2d = j2d

    return fit_j2osc_mean_elements!(j2oscd, vjd, vr_i, vv_i; kwargs...)
end

"""
    fit_j2osc_mean_elements!(
        j2oscd::J2OsculatingPropagator{Tepoch, T},
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

Fit a set of mean Keplerian elements for the J2 osculating orbit propagator `j2oscd` using
the osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    The J2 osculating orbit propagator `j2oscd` will be initialized with the Keplerian
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
# Allocate a new J2 osculating orbit propagator using a dummy set of Keplerian elements.
julia> j2oscd = j2osc_init(KeplerianElements(0.0, 7000e3, 0, 0, 0, 0, 0));

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

julia> orb, P = fit_j2osc_mean_elements!(j2oscd, vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the J2 osculating propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          4          1.69072e-05           0.00260088              2.60093         -2.38817e-10 %

(KeplerianElements{MeanAnomaly, Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.9999427807079662 -0.004901152245568996 … -0.00012021282775592569 -0.00011550583198160739; -0.004901152248513471 1.000732580069942 … 0.003266156907020634 3.8567169210589895e-5; … ; -0.00012021282776408784 0.0032661569070203842 … 2.2048267812548562e-5 5.3827666735694705e-8; -0.00011550583198608087 3.856716921258199e-5 … 5.382766674268047e-8 2.165742075903671e-5])

julia> orb
KeplerianElements{MeanAnomaly, Float64, Float64}:
              Epoch :    2.46003e6 (2023-03-24T16:33:40.388)
    Semi-major axis : 7135.8        km
       Eccentricity :    0.00135383
        Inclination :   98.4304     °
               RAAN :  162.113      °
  Arg. of Periapsis :   64.9256     °
       Mean Anomaly :  313.198      °
```
"""
function fit_j2osc_mean_elements!(
    j2oscd::J2OsculatingPropagator{Tepoch, T},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tepoch <: Number, T <: Number, Tjd <: Number, Tv <: AbstractVector}
    return _fit_mean_elements!(j2oscd, vjd, vr_i, vv_i; kwargs...)
end

"""
    update_j2osc_mean_elements_epoch(
        orb::KeplerianElements,
        new_epoch::Union{Number, DateTime}
    ) -> KeplerianElements{MeanAnomaly}

Update the epoch of the mean elements `orb` using a J2 osculating orbit propagator to
`new_epoch`, which can be represented by a Julian Day or a `DateTime`.

!!! note

    This algorithm version will allocate a new J2 osculating propagator with the default
    constants `J2C_EGM2008`. If another set of constants are required, use the function
    [`update_j2osc_mean_elements_epoch!`](@ref) instead.

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
              Epoch :    2.45995e6 (2023-01-01T00:00:00)
    Semi-major axis : 7190.98     km
       Eccentricity :    0.001111
        Inclination :   98.405    °
               RAAN :   90.0      °
  Arg. of Periapsis :  200.0      °
       True Anomaly :   45.0      °

julia> update_j2osc_mean_elements_epoch(orb, DateTime("2023-01-02"))
KeplerianElements{MeanAnomaly, Float64, Float64}:
              Epoch :    2.45995e6 (2023-01-02T00:00:00)
    Semi-major axis : 7190.98     km
       Eccentricity :    0.001111
        Inclination :   98.405    °
               RAAN :   90.9565   °
  Arg. of Periapsis :  197.078    °
       Mean Anomaly :  127.189    °
```
"""
function update_j2osc_mean_elements_epoch(
    orb::KeplerianElements{Tanomaly, Tepoch, T}, new_epoch::Union{Number, DateTime}
) where {Tanomaly <: AbstractAnomaly, Tepoch <: Number, T <: Number}
    # Allocate the J2 propagator structure that will propagate the mean elements.
    j2d = J2Propagator{Tepoch, T}()

    # Assign the constants, which are used in the initialization.
    j2d.j2c = J2C_EGM2008

    # Allocate the J2 osculating propagator structure.
    j2oscd = J2OsculatingPropagator{Tepoch, T}()
    j2oscd.j2d = j2d

    return update_j2osc_mean_elements_epoch!(j2oscd, orb, new_epoch)
end

"""
    update_j2osc_mean_elements_epoch!(
        j2oscd::J2OsculatingPropagator,
        orb::KeplerianElements,
        new_epoch::Union{Number, DateTime}
    ) -> KeplerianElements{MeanAnomaly}

Update the epoch of the mean elements `orb` using the propagator `j2oscd` to `new_epoch`,
which can be represented by a Julian Day or a `DateTime`.

!!! note

    The J2 osculating orbit propagator `j2oscd` will be initialized with the Keplerian
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
              Epoch :    2.45995e6 (2023-01-01T00:00:00)
    Semi-major axis : 7190.98     km
       Eccentricity :    0.001111
        Inclination :   98.405    °
               RAAN :   90.0      °
  Arg. of Periapsis :  200.0      °
       True Anomaly :   45.0      °

# Allocate a new J2 osculating orbit propagator using the created Keplerian elements. Notice
# that any set of Keplerian elements can be used here.
julia> j2oscd = j2osc_init(orb);

julia> update_j2osc_mean_elements_epoch!(j2oscd, orb, DateTime("2023-01-02"))
KeplerianElements{MeanAnomaly, Float64, Float64}:
              Epoch :    2.45995e6 (2023-01-02T00:00:00)
    Semi-major axis : 7190.98     km
       Eccentricity :    0.001111
        Inclination :   98.405    °
               RAAN :   90.9565   °
  Arg. of Periapsis :  197.078    °
       Mean Anomaly :  127.189    °
```
"""
function update_j2osc_mean_elements_epoch!(
    j2oscd::J2OsculatingPropagator, orb::KeplerianElements, new_epoch::DateTime
)
    dt = datetime2julian(new_epoch)
    return update_j2osc_mean_elements_epoch!(j2oscd, orb, dt)
end

function update_j2osc_mean_elements_epoch!(
    j2oscd::J2OsculatingPropagator, orb::KeplerianElements, new_epoch::Number
)
    orb = _update_mean_elements_epoch!(j2oscd, orb, new_epoch)

    # The generic algorithm initializes the J2 propagator only. Hence, we must compute the
    # osculating elements at the new epoch to complete the initialization.
    j2osc!(j2oscd, 0)

    return orb
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _similar_propagator(
        j2oscd::J2OsculatingPropagator{Tepoch},
        ::Type{T}
    ) where {Tepoch <: Number, T <: Number} -> J2OsculatingPropagator{Tepoch, T}

Create an uninitialized J2 osculating propagator with the same constants as `j2oscd`
converted to the element type `T`.
"""
function _similar_propagator(
    j2oscd::J2OsculatingPropagator{Tepoch}, ::Type{T}
) where {Tepoch <: Number, T <: Number}
    new_j2oscd = J2OsculatingPropagator{Tepoch, T}()
    new_j2oscd.j2d = _similar_propagator(j2oscd.j2d, T)

    return new_j2oscd
end

"""
    _j2osc_init!(j2oscd::J2OsculatingPropagator, orb₀::KeplerianElements) -> Nothing

Initialize the J2 osculating orbit propagator `j2oscd` with the mean Keplerian elements
`orb₀` [SI units] without computing the osculating elements at the initial instant. The
callers that propagate the orbit right afterwards use this function, since the osculating
elements would be overwritten anyway.
"""
function _j2osc_init!(j2oscd::J2OsculatingPropagator, orb₀::KeplerianElements)
    # Initialize the J2 propagator that will propagate the mean elements.
    j2_init!(j2oscd.j2d, orb₀)

    return nothing
end
