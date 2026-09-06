## Description #############################################################################
#
#   Two-Body orbit propagator.
#
#   This algorithm considers a perfect Keplerian orbit. In other words, no perturbation is
#   considered during the propagation and the Earth is modeled as a perfect sphere.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#     Press, Hawthorn, CA, USA.
#
############################################################################################

export TBC_M0, TBC_M0_F32
export twobody_init, twobody_init!, twobody, twobody!
export fit_twobody_mean_elements, fit_twobody_mean_elements!
export update_twobody_mean_elements_epoch, update_twobody_mean_elements_epoch!

############################################################################################
#                                        Constants                                         #
############################################################################################

# Earth's standard gravitational parameter [m³/s²]
const TBC_M0     = 3.986004415e14
const TBC_M0_F32 = 3.986004415f14

############################################################################################
#                                        Functions                                         #
############################################################################################

"""
    twobody_init(orb₀::KeplerianElements; kwargs...) -> TwoBodyPropagator

Create and initialize the two-body propagator structure using the mean Keplerian elements
`orb₀`.

!!! note

    The type used in the propagation will be the same as used to define the standard
    gravitational parameter `m0`.

# Keywords

- `m0::T`: Standard gravitational parameter of the central body [m³ / s²].
    (**Default**: `TBC_M0`)
"""
function twobody_init(
    orb₀::KeplerianElements{Tanomaly, Tepoch, Tkepler}; m0::Tm0 = TBC_M0
) where {Tanomaly <: AbstractAnomaly, Tepoch <: Number, Tkepler <: Number, Tm0 <: Number}
    T = _propagator_eltype(Tm0, Tkepler)

    # Allocate the propagator structure and assign the constant, which is used in the
    # initialization.
    tbd = TwoBodyPropagator{Tepoch, T}()
    tbd.μ = T(m0)

    # Initialize the propagator and return.
    twobody_init!(tbd, orb₀)

    return tbd
end

"""
    twobody_init!(tbd::TwoBodyPropagator, orb₀::KeplerianElements) -> Nothing

Initialize the two-body propagator structure `tbd` using the mean Keplerian elements `orb₀`.

!!! warning

    The propagation constant `μ::Number` in `tbd` will not be changed. Hence, it must be
    initialized.
"""
function twobody_init!(
    tbd::TwoBodyPropagator{Tepoch, T}, orb₀::KeplerianElements
) where {Tepoch <: Number, T <: Number}
    # Unpack elements.
    a₀ = T(orb₀.semi_major_axis)
    e₀ = T(orb₀.eccentricity)

    # Without this check, the user would get a `DomainError` from an internal square root,
    # or silently wrong results, instead of a message pointing at the offending element.
    if !(0 <= e₀ < 1)
        throw(
            ArgumentError("The eccentricity must be in the interval [0, 1), but it is $e₀.")
        )
    end

    if a₀ * (1 - e₀) <= 0
        throw(
            ArgumentError(
                "The perigee radius must be positive, but the semi-major axis is $a₀ m and the eccentricity is $e₀.",
            ),
        )
    end

    # Make sure the Keplerian elements use the mean anomaly. The conversion requires a
    # valid eccentricity, so it must happen after the checks above.
    ke₀ = convert(KeplerianElements{MeanAnomaly, Tepoch, T}, orb₀)

    # Compute the mean motion using the semi-major axis.
    n₀ = √(tbd.μ / a₀^3)

    # Create and return the two-body orbit propagator structure.
    tbd.orb₀ = ke₀
    tbd.orbk = ke₀
    tbd.Δt   = 0
    tbd.n₀   = n₀

    return nothing
end

"""
    twobody(
        Δt::Number,
        orb₀::KeplerianElements;
        kwargs...
    ) -> SVector{3, T}, SVector{3, T}, TwoBodyPropagator

Initialize the two-body propagator structure using the input elements `orb₀` and propagate
the orbit until the time Δt [s].

!!! note

    The type used in the propagation will be the same as used to define the standard
    gravitational parameter `m0`.

# Keywords

- `m0::T`: Standard gravitational parameter of the central body [m³ / s²].
    (**Default**: `TBC_M0`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`TwoBodyPropagator`](@ref): Structure with the initialized propagator.

# Remarks

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters.
"""
function twobody(Δt::Number, orb₀::KeplerianElements; m0::T = TBC_M0) where {T <: Number}
    tbd = twobody_init(orb₀; m0 = m0)
    r_i, v_i = twobody!(tbd, Δt)
    return r_i, v_i, tbd
end

"""
    twobody!(
        tbd::TwoBodyPropagator{Tepoch, T},
        t::Number
    ) where {Tepoch, T} -> SVector{3, T}, SVector{3, T}

Propagate the orbit defined in `tbd` (see [`TwoBodyPropagator`](@ref)) to `t` [s] after the
epoch of the input mean elements in `tbd`.

!!! note

    The internal values in `tbd` will be modified.

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.

# Remarks

The inertial frame in which the output is represented depends on which frame was used to
generate the orbit parameters.
"""
function twobody!(
    tbd::TwoBodyPropagator{Tepoch, T}, t::Number
) where {Tepoch <: Number, T <: Number}
    # Unpack.
    orb₀ = tbd.orb₀
    a₀   = orb₀.semi_major_axis
    e₀   = orb₀.eccentricity
    i₀   = orb₀.inclination
    Ω₀   = orb₀.raan
    ω₀   = orb₀.argument_of_periapsis
    M₀   = mean_anomaly(orb₀)

    # Time elapsed since epoch.
    epoch = orb₀.epoch
    Δt    = T(t)

    # Propagate the orbital elements.
    M_k = mod(M₀ + tbd.n₀ * Δt, T(2π))

    # Assemble the current mean elements.
    orbk = KeplerianElements{MeanAnomaly}(
        epoch + Tepoch(t) / 86400, a₀, e₀, i₀, Ω₀, ω₀, M_k
    )

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(orbk)

    # Update the propagator structure.
    tbd.Δt   = Δt
    tbd.orbk = orbk

    # Return the position and velocity vector represented in the inertial
    # reference frame.
    return r_i_k, v_i_k
end

"""
    fit_twobody_mean_elements(
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {
        Tjd <: Number,
        Tv <: AbstractVector
    } -> KeplerianElements{MeanAnomaly, Float64, Float64}, SMatrix{6, 6, Float64}

Fit a set of mean Keplerian elements for the two-body orbit propagator using the osculating
elements represented by a set of position vectors `vr_i` [m] and a set of velocity vectors
`vv_i` [m / s] represented in an inertial reference frame at instants in the array `vjd`
[Julian Day].

!!! note

    This algorithm version will allocate a new two-body propagator with the default
    gravitational parameter `TBC_M0`. If another value is required, use the function
    [`fit_twobody_mean_elements!`](@ref) instead.

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
    (**Default**: `nothing`)
- `jacobian_method::Union{FiniteDiffJacobian, ForwardDiffJacobian}`: Method used to compute
    the Jacobian matrix. Use `FiniteDiffJacobian()` for finite differences or
    `ForwardDiffJacobian()` for **ForwardDiff.jl** automatic differentiation.
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
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements [Julian Day].
    (**Default**: `vjd[end]`)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: `true`)
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

julia> orb, P = fit_twobody_mean_elements(vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the two-body propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          3          1.67397e-06           0.00142531              1.42531         -4.87682e-05 %

(KeplerianElements{MeanAnomaly, Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.9999772153857773 3.4367217753253005e-7 … -9.849257816722872e-5 -6.851478727900722e-5; 3.4367203183548976e-7 0.9999780011041396 … 0.0032585770446388936 2.506098261672399e-5; … ; -9.849257816766779e-5 0.0032585770446391334 … 2.199935616160074e-5 1.1527517043362011e-7; -6.851478727602333e-5 2.5060982618551414e-5 … 1.152751704386378e-7 2.19648292850571e-5])

julia> orb
KeplerianElements{MeanAnomaly, Float64, Float64}:
  Epoch             : 2.46003e6 (2023-03-24T16:33:40.388)
  Semi-Major Axis   : 7139.282635 km
  Eccentricity      : 0.001327104531
  Inclination       : 98.42538016°
  RA of Asc. Node   : 162.1096997°
  Arg. of Periapsis : 79.45743611°
  Mean Anomaly      : 298.683562°
```
"""
function fit_twobody_mean_elements(
    vjd::AbstractVector{Tjd}, vr_i::AbstractVector{Tv}, vv_i::AbstractVector{Tv}; kwargs...
) where {Tjd <: Number, Tv <: AbstractVector}
    # Allocate the two-body propagator structure that will propagate the mean elements.
    tbd = TwoBodyPropagator{Float64, Float64}()

    # Assign the constant, which is used in the initialization.
    tbd.μ = TBC_M0

    return fit_twobody_mean_elements!(tbd, vjd, vr_i, vv_i; kwargs...)
end

"""
    fit_twobody_mean_elements!(
        tbd::TwoBodyPropagator{Tepoch, T},
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {
        Tepoch <: Number,
        T <: Number,
        Tjd <: Number,
        Tv <: AbstractVector
    } -> KeplerianElements{MeanAnomaly, Tepoch, T}, SMatrix{6, 6, T}

Fit a set of mean Keplerian elements for the two-body orbit propagator `tbd` using the
osculating elements represented by a set of position vectors `vr_i` [m] and a set of
velocity vectors `vv_i` [m / s] represented in an inertial reference frame at instants in
the array `vjd` [Julian Day].

!!! note

    The two-body orbit propagator `tbd` will be initialized with the Keplerian elements
    returned by the function.

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
    (**Default**: `nothing`)
- `jacobian_method::Union{FiniteDiffJacobian, ForwardDiffJacobian}`: Method used to compute
    the Jacobian matrix. Use `FiniteDiffJacobian()` for finite differences or
    `ForwardDiffJacobian()` for **ForwardDiff.jl** automatic differentiation.
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
- `mean_elements_epoch::Number`: Epoch for the fitted mean elements [Julian Day].
    (**Default**: `vjd[end]`)
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default**: `true`)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default**: `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{MeanAnomaly, Tepoch, T}`: Fitted Keplerian elements.
- `SMatrix{6, 6, T}`: Final covariance matrix of the least-square algorithm.

# Examples

```julia-repl
# Allocate a new two-body orbit propagator using a dummy set of Keplerian elements.
julia> tbd = twobody_init(KeplerianElements(0.0, 7000e3, 0, 0, 0, 0, 0));

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

julia> orb, P = fit_twobody_mean_elements!(tbd, vjd, vr_i, vv_i)
ACTION:   Fitting the mean elements for the two-body propagator.
           Iteration        Position RMSE        Velocity RMSE           Total RMSE       RMSE Variation
                                     [km]             [km / s]                  [ ]
PROGRESS:          3          1.67397e-06           0.00142531              1.42531         -4.87682e-05 %

(KeplerianElements{MeanAnomaly, Float64, Float64}: Epoch = 2.46003e6 (2023-03-24T16:33:40.388), [0.9999772153857773 3.4367217753253005e-7 … -9.849257816722872e-5 -6.851478727900722e-5; 3.4367203183548976e-7 0.9999780011041396 … 0.0032585770446388936 2.506098261672399e-5; … ; -9.849257816766779e-5 0.0032585770446391334 … 2.199935616160074e-5 1.1527517043362011e-7; -6.851478727602333e-5 2.5060982618551414e-5 … 1.152751704386378e-7 2.19648292850571e-5])

julia> orb
KeplerianElements{MeanAnomaly, Float64, Float64}:
  Epoch             : 2.46003e6 (2023-03-24T16:33:40.388)
  Semi-Major Axis   : 7139.282635 km
  Eccentricity      : 0.001327104531
  Inclination       : 98.42538016°
  RA of Asc. Node   : 162.1096997°
  Arg. of Periapsis : 79.45743611°
  Mean Anomaly      : 298.683562°
```
"""
function fit_twobody_mean_elements!(
    tbd::TwoBodyPropagator{Tepoch, T},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    kwargs...,
) where {Tepoch <: Number, T <: Number, Tjd <: Number, Tv <: AbstractVector}
    return _fit_mean_elements!(tbd, vjd, vr_i, vv_i; kwargs...)
end

"""
    update_twobody_mean_elements_epoch(
        orb::KeplerianElements,
        new_epoch::Union{Number, DateTime}
    ) -> KeplerianElements{MeanAnomaly}

Update the epoch of the mean elements `orb` using a two-body orbit propagator to
`new_epoch`, which can be represented by a Julian Day or a `DateTime`.

!!! note

    This algorithm version will allocate a new two-body propagator with the default
    gravitational parameter `TBC_M0`. If another value is required, use the function
    [`update_twobody_mean_elements_epoch!`](@ref) instead.

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
  Epoch             : 2.45995e6 (2023-01-01T00:00:00)
  Semi-Major Axis   : 7190.982 km
  Eccentricity      : 0.001111
  Inclination       : 98.405°
  RA of Asc. Node   : 90.0°
  Arg. of Periapsis : 200.0°
  True Anomaly      : 45.0°

julia> update_twobody_mean_elements_epoch(orb, DateTime("2023-01-02"))
KeplerianElements{MeanAnomaly, Float64, Float64}:
  Epoch             : 2.45995e6 (2023-01-02T00:00:00)
  Semi-Major Axis   : 7190.982 km
  Eccentricity      : 0.001111
  Inclination       : 98.405°
  RA of Asc. Node   : 90.0°
  Arg. of Periapsis : 200.0°
  Mean Anomaly      : 130.2533301°
```
"""
function update_twobody_mean_elements_epoch(
    orb::KeplerianElements{Tanomaly, Tepoch, T}, new_epoch::Union{Number, DateTime}
) where {Tanomaly <: AbstractAnomaly, Tepoch <: Number, T <: Number}
    # Allocate the two-body propagator structure that will propagate the mean elements.
    tbd = TwoBodyPropagator{Tepoch, T}()

    # Assign the constant, which is used in the initialization.
    tbd.μ = TBC_M0

    return update_twobody_mean_elements_epoch!(tbd, orb, new_epoch)
end

"""
    update_twobody_mean_elements_epoch!(
        tbd::TwoBodyPropagator,
        orb::KeplerianElements,
        new_epoch::Union{Number, DateTime}
    ) -> KeplerianElements{MeanAnomaly}

Update the epoch of the mean elements `orb` using the propagator `tbd` to `new_epoch`, which
can be represented by a Julian Day or a `DateTime`.

!!! note

    The two-body orbit propagator `tbd` will be initialized with the Keplerian elements
    returned by the function.

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
  Epoch             : 2.45995e6 (2023-01-01T00:00:00)
  Semi-Major Axis   : 7190.982 km
  Eccentricity      : 0.001111
  Inclination       : 98.405°
  RA of Asc. Node   : 90.0°
  Arg. of Periapsis : 200.0°
  True Anomaly      : 45.0°

# Allocate a new two-body orbit propagator using the created Keplerian elements. Notice that
# any set of Keplerian elements can be used here.
julia> tbd = twobody_init(orb);

julia> update_twobody_mean_elements_epoch!(tbd, orb, DateTime("2023-01-02"))
KeplerianElements{MeanAnomaly, Float64, Float64}:
  Epoch             : 2.45995e6 (2023-01-02T00:00:00)
  Semi-Major Axis   : 7190.982 km
  Eccentricity      : 0.001111
  Inclination       : 98.405°
  RA of Asc. Node   : 90.0°
  Arg. of Periapsis : 200.0°
  Mean Anomaly      : 130.2533301°
```
"""
function update_twobody_mean_elements_epoch!(
    tbd::TwoBodyPropagator, orb::KeplerianElements, new_epoch::DateTime
)
    dt = datetime2julian(new_epoch)
    return update_twobody_mean_elements_epoch!(tbd, orb, dt)
end

function update_twobody_mean_elements_epoch!(
    tbd::TwoBodyPropagator, orb::KeplerianElements, new_epoch::Number
)
    return _update_mean_elements_epoch!(tbd, orb, new_epoch)
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _similar_propagator(
        tbd::TwoBodyPropagator{Tepoch},
        ::Type{T}
    ) where {Tepoch <: Number, T <: Number} -> TwoBodyPropagator{Tepoch, T}

Create an uninitialized two-body propagator with the same gravitational parameter as `tbd`
converted to the element type `T`.
"""
function _similar_propagator(
    tbd::TwoBodyPropagator{Tepoch}, ::Type{T}
) where {Tepoch <: Number, T <: Number}
    new_tbd = TwoBodyPropagator{Tepoch, T}()
    new_tbd.μ = T(tbd.μ)

    return new_tbd
end
