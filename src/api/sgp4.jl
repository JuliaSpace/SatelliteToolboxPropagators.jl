## Description #############################################################################
#
#  API implementation for SGP4 orbit propagator.
#
############################################################################################

Propagators.epoch(orbp::OrbitPropagatorSgp4)         = orbp.sgp4d.epoch
Propagators.last_instant(orbp::OrbitPropagatorSgp4)  = orbp.sgp4d.Δt * 60
Propagators.name(orbp::OrbitPropagatorSgp4)          = "SGP4 Orbit Propagator"

function Propagators.mean_elements(orbp::OrbitPropagatorSgp4)
    # We need to copy the propagator to avoid modifying it.
    sgp4d = deepcopy(orbp.sgp4d)
    sgp4c = sgp4d.sgp4c

    # First, we need to create a TLE based on the initial parameters.
    dt  = julian2datetime(Propagators.epoch(orbp))
    dt₀ = DateTime(Year(dt))

    dt_year    = year(dt)
    epoch_year = dt_year < 1980 ? dt_year - 1900 : dt_year - 2000
    epoch_day  = (dt - dt₀).value / 1000 / 86400 + 1

    tle = TLE(
        epoch_year          = epoch_year,
        epoch_day           = epoch_day,
        bstar               = sgp4d.bstar,
        inclination         = sgp4d.i₀ |> rad2deg,
        raan                = sgp4d.Ω₀ |> rad2deg,
        eccentricity        = sgp4d.e₀,
        argument_of_perigee = sgp4d.ω₀ |> rad2deg,
        mean_anomaly        = sgp4d.M₀ |> rad2deg,
        mean_motion         = 720 * sgp4d.n₀ / π,
    )

    # Now, we update the TLE epoch to the current propagation instant.
    new_epoch = Propagators.epoch(orbp) + Propagators.last_instant(orbp) / 86400
    updated_tle = update_sgp4_tle_epoch!(sgp4d, tle, new_epoch; verbose = false)

    # Initialize the propagator to obtain the mean elements.
    sgp4_init!(sgp4d, updated_tle)

    # Create and return the Keplerian elements.
    return KeplerianElements(
        new_epoch,
        (sgp4c.XKE / sgp4d.n₀)^(2 / 3) * (1000 * sgp4c.R0),
        sgp4d.e₀,
        sgp4d.i₀,
        sgp4d.Ω₀,
        sgp4d.ω₀,
        mean_to_true_anomaly(sgp4d.e₀, sgp4d.M₀)
    )
end

"""
    Propagators.fit_mean_elements(::Val{:SGP4}, vjd::AbstractVector{Tjd}, vr_teme::AbstractVector{Tv}, vv_teme::AbstractVector{Tv}; kwargs...) -> KeplerianElements{Float64, Float64}, SMatrix{6, 6, Float64}

Fit a Two-Line Element set (`TLE`) for the SGP4 orbit propagator using the osculating
elements represented by a set of position vectors `vr_teme` [m] and a set of velocity
vectors `vv_teme` [m / s] represented in the True-Equator, Mean-Equinox reference frame
(TEME) at instants in the array `vjd` [Julian Day].

This algorithm was based on **[1]**.

!!! note

    This algorithm version will allocate a new SGP4 propagator with the default constants
    `sgp4c_wgs84`. If another set of constants are required, use the function
    [`Propagators.fit_mean_elements!`](@ref) instead.

# Keywords

- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default** = 2e-4)
- `estimate_bstar::Bool`: If `true`, the algorithm will try to estimate the B* parameter.
    Otherwise, it will be set to 0 or to the value in initial guess (see section  **Initial
    Guess**).
    (**Default** = true)
- `initial_guess::Union{Nothing, AbstractVector, TLE}`: Initial guess for the TLE fitting
    process. If it is `nothing`, the algorithm will obtain an initial estimate from the
    osculating elements in `vr_teme` and `vv_teme`. For more information, see the section
    **Initial Guess**.
    (**Default** = nothing)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix.
    (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until it absolute value is higher than
    `jacobian_perturbation_tol`.
    (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted TLE.
    (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default** = `@SVector(ones(Bool, 6))`)
- `classification::Char`: Satellite classification character for the output TLE.
    (**Default** = 'U')
- `element_set_number::Int`: Element set number for the output TLE.
    (**Default** = 0)
- `international_designator::String`: International designator string for the output TLE.
    (**Default** = "999999")
- `name::String`: Satellite name for the output TLE.
    (**Default** = "UNDEFINED")
- `revolution_number::Int`: Revolution number for the output TLE.
    (**Default** = 0)
- `satellite_number::Int`: Satellite number for the output TLE.
    (**Default** = 9999)

# Returns

- `TLE`: The fitted TLE.
- `SMatrix{7, 7, T}`: Final covariance matrix of the least-square algorithm.

# Initial Guess

This algorithm uses a least-square algorithm to fit a TLE based on a set of osculating state
vectors. Since the system is chaotic, a good initial guess is paramount for algorithm
convergence. We can provide an initial guess using the keyword `initial_guess`.

If `initial_guess` is a `TLE`, we update the TLE epoch using the function
[`update_sgp4_tle_epoch!`](@ref) to the desired one in `mean_elements_epoch`. Afterward, we
use this new TLE as the initial guess.

If `initial_guess` is an `AbstractVector`, we use this vector as the initial mean state
vector for the algorithm. It must contain 7 elements as follows:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

If `initial_guess` is `nothing`, the algorithm takes the closest osculating state vector to
the `mean_elements_epoch` and uses it as the initial mean state vector. In this case, the
epoch is set to the same epoch of the osculating data in `vjd`. When the fitted TLE is
obtained, the algorithm uses the function [`update_sgp4_tle_epoch!`](@ref) to change its
epoch to `mean_elements_epoch`.

!!! note

    If `initial_guess` is not `nothing`, the B* initial estimate is obtained from the TLE or
    the state vector. Hence, if `estimate_bstar` is `false`, it will be kept constant with
    this initial value.

# References

- **[1]** Vallado, D. A., Crawford, P (2008). SGP4 Orbit Determination. American Institute
    of Aeronautics ans Astronautics.
"""
function Propagators.fit_mean_elements(
    ::Val{:SGP4},
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv};
    kwargs...
) where {Tjd<:Number, Tv<:AbstractVector}
    return fit_sgp4_tle(vjd, vr_teme ./ 1000, vv_teme ./ 1000; kwargs...)
end

"""
    Propagators.fit_mean_elements!(orbp::OrbitPropagatorSgp4, vjd::AbstractVector{Tjd}, vr_teme::AbstractVector{Tv}, vv_teme::AbstractVector{Tv}; kwargs...) where {Tjd<:Number, Tv<:AbstractVector}

Fit a Two-Line Element set (`TLE`) for the SGP4 orbit propagator `orbp` using the
osculating elements represented by a set of position vectors `vr_teme` [m] and a set of
velocity vectors `vv_teme` [m / s] represented in the True-Equator, Mean-Equinox reference
frame (TEME) at instants in the array `vjd` [Julian Day].

This algorithm was based on **[1]**.

!!! note

    The SGP4 orbit propagator `orbp` will be initialized with the TLE returned by the
    function.

# Keywords

- `atol::Number`: Tolerance for the residue absolute value. If the residue is lower than
    `atol` at any iteration, the computation loop stops.
    (**Default** = 2e-4)
- `rtol::Number`: Tolerance for the relative difference between the residues. If the
    relative difference between the residues in two consecutive iterations is lower than
    `rtol`, the computation loop stops.
    (**Default** = 2e-4)
- `estimate_bstar::Bool`: If `true`, the algorithm will try to estimate the B* parameter.
    Otherwise, it will be set to 0 or to the value in initial guess (see section  **Initial
    Guess**).
    (**Default** = true)
- `initial_guess::Union{Nothing, AbstractVector, TLE}`: Initial guess for the TLE fitting
    process. If it is `nothing`, the algorithm will obtain an initial estimate from the
    osculating elements in `vr_teme` and `vv_teme`. For more information, see the section
    **Initial Guess**.
    (**Default** = nothing)
- `jacobian_perturbation::Number`: Initial state perturbation to compute the
    finite-difference when calculating the Jacobian matrix.
    (**Default** = 1e-3)
- `jacobian_perturbation_tol::Number`: Tolerance to accept the perturbation when calculating
    the Jacobian matrix. If the computed perturbation is lower than
    `jacobian_perturbation_tol`, we increase it until it absolute value is higher than
    `jacobian_perturbation_tol`.
    (**Default** = 1e-7)
- `max_iterations::Int`: Maximum number of iterations allowed for the least-square fitting.
    (**Default** = 50)
- `mean_elements_epoch::Number`: Epoch for the fitted TLE.
    (**Default** = vjd[end])
- `verbose::Bool`: If `true`, the algorithm prints debugging information to `stdout`.
    (**Default** = true)
- `weight_vector::AbstractVector`: Vector with the measurements weights for the least-square
    algorithm. We assemble the weight matrix `W` as a diagonal matrix with the elements in
    `weight_vector` at its diagonal.
    (**Default** = `@SVector(ones(Bool, 6))`)
- `classification::Char`: Satellite classification character for the output TLE.
    (**Default** = 'U')
- `element_set_number::Int`: Element set number for the output TLE.
    (**Default** = 0)
- `international_designator::String`: International designator string for the output TLE.
    (**Default** = "999999")
- `name::String`: Satellite name for the output TLE.
    (**Default** = "UNDEFINED")
- `revolution_number::Int`: Revolution number for the output TLE.
    (**Default** = 0)
- `satellite_number::Int`: Satellite number for the output TLE.
    (**Default** = 9999)

# Returns

- `TLE`: The fitted TLE.
- `SMatrix{7, 7, T}`: Final covariance matrix of the least-square algorithm.

# Initial Guess

This algorithm uses a least-square algorithm to fit a TLE based on a set of osculating state
vectors. Since the system is chaotic, a good initial guess is paramount for algorithm
convergence. We can provide an initial guess using the keyword `initial_guess`.

If `initial_guess` is a `TLE`, we update the TLE epoch using the function
[`update_sgp4_tle_epoch!`](@ref) to the desired one in `mean_elements_epoch`. Afterward, we
use this new TLE as the initial guess.

If `initial_guess` is an `AbstractVector`, we use this vector as the initial mean state
vector for the algorithm. It must contain 7 elements as follows:

    ┌                                    ┐
    │ IDs 1 to 3: Mean position [km]     │
    │ IDs 4 to 6: Mean velocity [km / s] │
    │ ID  7:      Bstar         [1 / er] │
    └                                    ┘

If `initial_guess` is `nothing`, the algorithm takes the closest osculating state vector to
the `mean_elements_epoch` and uses it as the initial mean state vector. In this case, the
epoch is set to the same epoch of the osculating data in `vjd`. When the fitted TLE is
obtained, the algorithm uses the function [`update_sgp4_tle_epoch!`](@ref) to change its
epoch to `mean_elements_epoch`.

!!! note

    If `initial_guess` is not `nothing`, the B* initial estimate is obtained from the TLE or
    the state vector. Hence, if `estimate_bstar` is `false`, it will be kept constant with
    this initial value.

# References

- **[1]** Vallado, D. A., Crawford, P (2008). SGP4 Orbit Determination. American Institute
    of Aeronautics ans Astronautics.
"""
function Propagators.fit_mean_elements!(
    orbp::OrbitPropagatorSgp4,
    vjd::AbstractVector{Tjd},
    vr_teme::AbstractVector{Tv},
    vv_teme::AbstractVector{Tv};
    kwargs...
) where {Tjd<:Number, Tv<:AbstractVector}
    return fit_sgp4_tle!(orbp.sgp4d, vjd, vr_teme ./ 1000, vv_teme ./ 1000; kwargs...)
end

"""
    Propagators.init(Val(:SGP4), epoch::Number, n₀::Number, e₀::Number, i₀::Number, Ω₀::Number, ω₀::Number, M₀::Number, bstar::Number; kwargs...) -> OrbitPropagatorSgp4
    Propagators.init(Val(:SGP4), tle::TLE; kwargs...) -> OrbitPropagatorSgp4

Create and initialize the SGP4 orbit propagator structure using the initial orbit specified
by the arguments.

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `sgp4c`.

# Arguments

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/s].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).
- `tle::TLE`: Two-line elements used for the initialization.

# Keywords

- `sgp4c::Sgp4Constants`: SGP4 orbit propagator constants (see [`Sgp4Constants`](@ref)).
    (**Default** = `sgp4c_wgs84`)
"""
function Propagators.init(::Val{:SGP4}, tle::TLE; sgp4c::Sgp4Constants = sgp4c_wgs84)
    sgp4d = sgp4_init(tle; sgp4c = sgp4c)
    return OrbitPropagatorSgp4(sgp4d)
end

function Propagators.init(
    ::Val{:SGP4},
    epoch::Number,
    n₀::Number,
    e₀::Number,
    i₀::Number,
    Ω₀::Number,
    ω₀::Number,
    M₀::Number,
    bstar::Number;
    sgp4c::Sgp4Constants = sgp4c_wgs84
)
    sgp4d = sgp4_init(epoch, 60n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar; sgp4c = sgp4c)
    return OrbitPropagatorSgp4(sgp4d)
end

"""
    Propagators.init!(orbp::OrbitPropagatorSgp4, epoch::Number, n₀::Number, e₀::Number, i₀::Number, Ω₀::Number, ω₀::Number, M₀::Number, bstar::Number; kwargs...) -> Nothing
    Propagators.init!(orbp::OrbitPropagatorSgp4, tle::TLE; kwargs...) -> Nothing

Initialize the SGP4 orbit propagator structure `orbp` using the initial orbit specified
by the arguments.

!!! warning

    The propagation constants `sgp4c::Sgp4Constants` in `orbp.sgp4d` will not be changed.
    Hence, they must be initialized.

# Arguments

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/s].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).
- `tle::TLE`: Two-line elements used for the initialization.
"""
function Propagators.init!(orbp::OrbitPropagatorSgp4, tle::TLE)
    sgp4_init!(orbp.sgp4d, tle)
    return nothing
end

function Propagators.init!(
    orbp::OrbitPropagatorSgp4,
    epoch::Number,
    n₀::Number,
    e₀::Number,
    i₀::Number,
    Ω₀::Number,
    ω₀::Number,
    M₀::Number,
    bstar::Number
)
    sgp4_init!(orbp.sgp4d, epoch, 60n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar)
    return nothing
end

"""
    Propagators.propagate(Val(:SGP4), Δt::Number, epoch::Number, n₀::Number, e₀::Number, i₀::Number, Ω₀::Number, ω₀::Number, M₀::Number, bstar::Number; kwargs...)
    Propagators.propagate(Val(:SGP4), Δt::Number, tle::TLE; kwargs...)

Initialize the SGP4 propagator structure using the input arguments and propagate the orbit
until the time Δt [s].

!!! note

    The type used in the propagation will be the same as used to define the constants in the
    structure `sgp4c`.

# Arguments

- `epoch::Number`: Epoch of the orbital elements [Julian Day].
- `n₀::Number`: SGP type "mean" mean motion at epoch [rad/s].
- `e₀::Number`: "Mean" eccentricity at epoch.
- `i₀::Number`: "Mean" inclination at epoch [rad].
- `Ω₀::Number`: "Mean" longitude of the ascending node at epoch [rad].
- `ω₀::Number`: "Mean" argument of perigee at epoch [rad].
- `M₀::Number`: "Mean" mean anomaly at epoch [rad].
- `bstar::Number`: Drag parameter (B*).
- `tle::TLE`: Two-line elements used for the initialization.

# Keywords

- `sgp4c::Sgp4Constants{T}`: SGP4 orbit propagator constants (see [`Sgp4Constants`](@ref)).
    (**Default** = `sgp4c_wgs84`)

# Returns

- `SVector{3, T}`: Position vector [m] represented in the inertial frame at propagation
    instant.
- `SVector{3, T}`: Velocity vector [m / s] represented in the inertial frame at propagation
    instant.
- [`OrbitPropagatorSgp4`](@ref): Structure with the initialized parameters.
"""
function Propagators.propagate(
    ::Val{:SGP4},
    Δt::Number,
    tle::TLE;
    sgp4c::Sgp4Constants = sgp4c_wgs84
)
    r_i, v_i, sgp4d = sgp4(Δt / 60, tle; sgp4c = sgp4c)
    return 1000r_i, 1000v_i, OrbitPropagatorSgp4(sgp4d)
end

function Propagators.propagate(
    ::Val{:SGP4},
    Δt::Number,
    epoch::Number,
    n₀::Number,
    e₀::Number,
    i₀::Number,
    Ω₀::Number,
    ω₀::Number,
    M₀::Number,
    bstar::Number;
    sgp4c::Sgp4Constants = sgp4c_wgs84
)
    r_i, v_i, sgp4d = sgp4(
        Δt / 60,
        epoch,
        60n₀,
        e₀,
        i₀,
        Ω₀,
        ω₀,
        M₀,
        bstar;
        sgp4c = sgp4c
    )
    return 1000r_i, 1000v_i, OrbitPropagatorSgp4(sgp4d)
end

function Propagators.propagate!(orbp::OrbitPropagatorSgp4, t::Number)
    # Auxiliary variables.
    sgp4d = orbp.sgp4d

    # Propagate the orbit.
    r_i, v_i = sgp4!(sgp4d, t / 60)

    return 1000r_i, 1000v_i
end
