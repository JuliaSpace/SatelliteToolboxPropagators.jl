## Description #############################################################################
#
# Generic least-square algorithm to fit a set of mean elements from osculating state
# vectors.
#
# The analytical propagators share the same algorithm. They only differ in how the
# propagator is initialized and propagated, which the small set of methods below provides
# for each of them.
#
############################################################################################

############################################################################################
#                                        Constants                                         #
############################################################################################

# The decorated strings printed by the fitting algorithm are rendered once here with
# **StyledStrings**, in the plain and in the colored versions, so that the algorithm only
# prints plain strings. Otherwise, printing an annotated string inside the algorithm adds
# many allocation sites to it, even though they are only reachable when `verbose` is `true`.
const _FIT_ACTION_TAG = (
    "ACTION:",
    sprint(
        print, styled"{(foreground=yellow,weight=bold):ACTION:}"; context = :color => true
    ),
)

const _FIT_PROGRESS_TAG = (
    "PROGRESS:", sprint(print, styled"{bold:PROGRESS:}"; context = :color => true)
)

const _FIT_HEADER = let
    header = @sprintf(
        "%10s %20s %20s %20s %20s",
        "Iteration",
        "Position RMSE",
        "Velocity RMSE",
        "Total RMSE",
        "RMSE Variation"
    )
    (
        header,
        sprint(
            print,
            styled"{(foreground=yellow,weight=bold):$header}";
            context = :color => true,
        ),
    )
end

const _FIT_UNITS = let
    units = @sprintf("%10s %20s %20s %20s", "", "[km]", "[km / s]", "[ ]")
    (units, sprint(print, styled"{bold:$units}"; context = :color => true))
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# Union of the propagator structures that support fitting a set of mean elements from
# osculating state vectors using `_fit_mean_elements!`.
const AbstractMeanElementsPropagator{Tepoch, T} = Union{
    J2Propagator{Tepoch, T},
    J2OsculatingPropagator{Tepoch, T},
    J4Propagator{Tepoch, T},
    J4OsculatingPropagator{Tepoch, T},
    TwoBodyPropagator{Tepoch, T},
}

# == Propagator-Specific Operations ========================================================
#
# The generic algorithm only needs to know how to name, initialize, and propagate each
# structure, and which structure stores the mean elements. Everything else is shared.

"""
    _propagator_name(pd::AbstractMeanElementsPropagator) -> String

Return the name of the propagator `pd` used in the messages of the fitting algorithm.
"""
_propagator_name(::J2Propagator)           = "J2"
_propagator_name(::J2OsculatingPropagator) = "J2 osculating"
_propagator_name(::J4Propagator)           = "J4"
_propagator_name(::J4OsculatingPropagator) = "J4 osculating"
_propagator_name(::TwoBodyPropagator)      = "two-body"

"""
    _mean_elements_init!(
        pd::AbstractMeanElementsPropagator,
        orb::KeplerianElements
    ) -> Nothing

Initialize the propagator `pd` with the mean elements `orb` [SI units]. The osculating
propagators use the internal initialization, which skips the propagation to the initial
instant, because the algorithm always propagates right afterwards.
"""
_mean_elements_init!(pd::J2Propagator, orb)           = j2_init!(pd, orb)
_mean_elements_init!(pd::J2OsculatingPropagator, orb) = _j2osc_init!(pd, orb)
_mean_elements_init!(pd::J4Propagator, orb)           = j4_init!(pd, orb)
_mean_elements_init!(pd::J4OsculatingPropagator, orb) = _j4osc_init!(pd, orb)
_mean_elements_init!(pd::TwoBodyPropagator, orb)      = twobody_init!(pd, orb)

"""
    _mean_elements_propagate!(
        pd::AbstractMeanElementsPropagator{Tepoch, T},
        Δt::Number
    ) where {Tepoch <: Number, T <: Number} -> SVector{3, T}, SVector{3, T}

Propagate the orbit using `pd` to `Δt` [s] after the epoch of its initial mean elements,
returning the position [m] and the velocity [m / s] vectors in the inertial frame.
"""
_mean_elements_propagate!(pd::J2Propagator, Δt)           = j2!(pd, Δt)
_mean_elements_propagate!(pd::J2OsculatingPropagator, Δt) = j2osc!(pd, Δt)
_mean_elements_propagate!(pd::J4Propagator, Δt)           = j4!(pd, Δt)
_mean_elements_propagate!(pd::J4OsculatingPropagator, Δt) = j4osc!(pd, Δt)
_mean_elements_propagate!(pd::TwoBodyPropagator, Δt)      = twobody!(pd, Δt)

"""
    _mean_elements_propagator(
        pd::AbstractMeanElementsPropagator
    ) -> Union{J2Propagator, J4Propagator, TwoBodyPropagator}

Return the structure inside `pd` that stores the initial (`orb₀`) and the current (`orbk`)
mean elements, which is `pd` itself for the propagators of mean elements.
"""
_mean_elements_propagator(pd::J2Propagator)           = pd
_mean_elements_propagator(pd::J2OsculatingPropagator) = pd.j2d
_mean_elements_propagator(pd::J4Propagator)           = pd
_mean_elements_propagator(pd::J4OsculatingPropagator) = pd.j4d
_mean_elements_propagator(pd::TwoBodyPropagator)      = pd

# == Algorithm =============================================================================

"""
    _fit_mean_elements!(
        pd::AbstractMeanElementsPropagator{Tepoch, T},
        vjd::AbstractVector{Tjd},
        vr_i::AbstractVector{Tv},
        vv_i::AbstractVector{Tv};
        kwargs...
    ) where {
        Tepoch <: Number,
        T <: Number,
        Tjd <: Number,
        Tv <: AbstractVector
    } -> KeplerianElements{MeanAnomaly, Tepoch, T}, SMatrix{6, 6, T}, NamedTuple

Fit a set of mean Keplerian elements for the propagator `pd` using the osculating elements
represented by a set of position vectors `vr_i` [m] and a set of velocity vectors `vv_i`
[m / s] represented in an inertial reference frame at instants in the array `vjd` [Julian
Day]. The algorithm is an iterative weighted least-square that propagates the current
estimate of the mean state vector to every measurement instant and corrects it using the
Jacobian of the propagation. `pd` is left initialized with the fitted elements. The fitting
fails if the residual diverges.

# Keywords

- `atol::Number`: Tolerance for the residual absolute value that stops the iterations.
    (**Default**: 2e-4)
- `rtol::Number`: Tolerance for the relative residual variation that stops the iterations.
    (**Default**: 2e-4)
- `initial_guess::Union{Nothing, KeplerianElements}`: Initial guess for the mean elements.
    If it is `nothing`, the closest measurement to `mean_elements_epoch` is used.
    (**Default**: `nothing`)
- `jacobian_method::AbstractJacobianMethod`: Method used to compute the Jacobian matrix.
    (**Default**: `FiniteDiffJacobian()`)
- `jacobian_perturbation::Number`: Initial relative state perturbation of the
    finite-difference Jacobian.
    (**Default**: 1e-3)
- `jacobian_perturbation_tol::Number`: Minimum absolute perturbation of the
    finite-difference Jacobian.
    (**Default**: 1e-7)
- `max_iterations::Int`: Maximum number of iterations.
    (**Default**: 50)
- `mean_elements_epoch::Union{Number, DateTime}`: Epoch of the fitted mean elements,
    represented by a Julian Day [UTC] or a `DateTime` [UTC].
    (**Default**: `vjd[end]`)
- `verbose::Bool`: If `true`, the algorithm prints its progress to `stdout`.
    (**Default**: `true`)
- `weight_vector::AbstractVector`: Diagonal of the weight matrix of the least-square
    algorithm, with six elements.
    (**Default**: `@SVector(ones(Bool, 6))`)

# Returns

- `KeplerianElements{MeanAnomaly, Tepoch, T}`: Fitted mean Keplerian elements [SI units].
- `SMatrix{6, 6, T}`: Final covariance matrix of the least-square algorithm.
- `NamedTuple`: Statistics of the least-square algorithm with the following fields:
    - `converged::Bool`: `true` if the iterations stopped because the residue was lower
        than `atol` or its relative variation was lower than `rtol`, or `false` if they
        stopped by reaching `max_iterations`.
    - `iterations::Int`: Number of iterations performed.
    - `position_rmse::T`: RMSE of the position residue in the last iteration [m].
    - `velocity_rmse::T`: RMSE of the velocity residue in the last iteration [m / s].
    - `total_rmse::T`: Weighted RMSE of the residue in the last iteration.

    The statistics refer to the fitting of the mean elements. If their epoch is updated
    afterward to match `mean_elements_epoch`, the statistics of that update are not
    returned.

# Extended help

## Throws

- `ArgumentError`: If `vjd`, `vr_i`, and `vv_i` do not have the same length, if
    `weight_vector` does not have six elements, or if `max_iterations` is lower than 1.
- `MeanElementsFitDivergenceError`: If the residual increases for three consecutive
    iterations and exceeds 5e11, indicating that the iterations diverged.
"""
function _fit_mean_elements!(
    pd::AbstractMeanElementsPropagator{Tepoch, T},
    vjd::AbstractVector{Tjd},
    vr_i::AbstractVector{Tv},
    vv_i::AbstractVector{Tv};
    atol::Number                                     = 2e-4,
    rtol::Number                                     = 2e-4,
    initial_guess::Union{Nothing, KeplerianElements} = nothing,
    jacobian_method::AbstractJacobianMethod          = FiniteDiffJacobian(),
    jacobian_perturbation::Number                    = 1e-3,
    jacobian_perturbation_tol::Number                = 1e-7,
    max_iterations::Int                              = 50,
    mean_elements_epoch::Union{Number, DateTime}     = vjd[end],
    verbose::Bool                                    = true,
    weight_vector::AbstractVector                    = @SVector(ones(Bool, 6)),
) where {T <: Number, Tepoch <: Number, Tjd <: Number, Tv <: AbstractVector}
    # Number of available measurements.
    num_measurements = length(vjd)

    # Check the inputs.
    length(vr_i) != num_measurements &&
        throw(ArgumentError("The number of elements in `vjd` and `vr_i` must be the same."))

    length(vv_i) != num_measurements &&
        throw(ArgumentError("The number of elements in `vjd` and `vv_i` must be the same."))

    if length(weight_vector) != 6
        throw(ArgumentError("The weight vector must have 6 elements."))
    end

    if max_iterations < 1
        throw(ArgumentError("The maximum number of iterations must be at least 1."))
    end

    # Desired epoch of the fitted mean elements [Julian Day].
    desired_epoch = _julian_day(mean_elements_epoch)

    # Check once if `stdout` supports colors, which selects the decorated strings printed by
    # the algorithm when `verbose` is `true`.
    has_color = get(stdout, :color, false)::Bool

    # Assemble the weight vector. Since the weight matrix is diagonal, we store only the
    # diagonal to improve performance by avoiding the Diagonal wrapper.
    W = @SVector T[
        weight_vector[begin],
        weight_vector[1 + begin],
        weight_vector[2 + begin],
        weight_vector[3 + begin],
        weight_vector[4 + begin],
        weight_vector[5 + begin],
    ]

    # Initial guess of the mean elements.
    #
    # NOTE: x₁ is the previous estimate and x₂ is the current estimate.
    if !isnothing(initial_guess)
        epoch = Tepoch(desired_epoch)

        # First, we need to update the mean elements to the desired epoch.
        verbose && _fit_print_action(
            has_color,
            "Updating the epoch of the initial mean elements guess to match the desired one.",
        )
        orb = _update_mean_elements_epoch!(pd, initial_guess, epoch)

        r_i, v_i = kepler_to_rv(orb)
        x₁ = SVector{6, T}(r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3])
    else
        # In this case, we must find the closest osculating vector to the desired epoch.
        id = firstindex(vjd)
        v  = abs(vjd[id] - desired_epoch)

        for k in eachindex(vjd)
            vk = abs(vjd[k] - desired_epoch)
            if vk < v
                id = k
                v  = vk
            end
        end

        epoch = Tepoch(vjd[id])
        r_i   = vr_i[id]
        v_i   = vv_i[id]
        x₁    = SVector{6, T}(r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3])
    end

    x₂ = x₁

    # Number of states in the input vector.
    num_states = 6

    # Statistics returned after the iterations.
    converged  = false
    iterations = 0
    σ_i        = T(0)
    σp_i       = T(0)
    σv_i       = T(0)

    # Variable to store the last residue.
    σ_i_₁ = T(0)

    # Variable to store how many iterations the residue increased. This is used to account
    # for divergence.
    Δd = 0

    # Header.
    if verbose
        _fit_print_action(
            has_color,
            "Fitting the mean elements for the $(_propagator_name(pd)) propagator.",
        )

        _fit_print_header(has_color)
    end

    # We need a reference to the covariance inverse because we will invert it and return
    # after the iterations.
    ΣJ′WJ = @SMatrix zeros(T, num_states, num_states)

    pd_ad = jacobian_method isa ForwardDiffJacobian ? _create_ad_propagator(pd) : nothing

    # The finite-difference Jacobian overwrites the propagator it uses. Giving it its own
    # propagator lets us initialize `pd` once per iteration instead of once per
    # measurement.
    pd_fd = jacobian_method isa FiniteDiffJacobian ? _create_fd_propagator(pd) : nothing

    # Loop until the maximum allowed iteration. The measurement loop only accesses the input
    # vectors with indices derived from `num_measurements`, whose consistency was checked
    # above, and fixed indices of static arrays. Hence, we can skip the bounds checking and
    # avoid copies when slicing the static arrays.
    @inbounds @views for it in 1:max_iterations
        x₁ = x₂
        iterations = it

        # Variables to store the summations to compute the least square fitting algorithm.
        ΣJ′WJ = @SMatrix zeros(T, num_states, num_states)
        ΣJ′Wb = @SVector zeros(T, num_states)

        # Variables to store the RMS errors in this iteration.
        σ_i  = T(0)
        σp_i = T(0)
        σv_i = T(0)

        # The estimated mean elements do not change within an iteration, so the propagator
        # only needs to be initialized once, outside the measurement loop.
        orb = rv_to_kepler(x₁[SOneTo(3)], x₁[StaticArrays.SUnitRange(4, 6)], epoch)
        _mean_elements_init!(pd, orb)

        for k in 1:num_measurements
            # Obtain the measured ephemerides.
            r_i = vr_i[k - 1 + begin]
            v_i = vv_i[k - 1 + begin]

            y = SVector{6, T}(
                r_i[begin],
                r_i[1 + begin],
                r_i[2 + begin],
                v_i[begin],
                v_i[1 + begin],
                v_i[2 + begin],
            )

            # Obtain the propagation time for this measurement.
            Δt = (vjd[k - 1 + begin] - epoch) * 86400

            # Propagate the orbit.
            r̂_i, v̂_i = _mean_elements_propagate!(pd, Δt)
            ŷ = vcat(r̂_i, v̂_i)

            # Compute the residue.
            b = y - ŷ

            # Compute the Jacobian in-place.
            J = _mean_elements_jacobian(
                jacobian_method,
                pd,
                Δt,
                x₁,
                ŷ;
                perturbation     = jacobian_perturbation,
                perturbation_tol = jacobian_perturbation_tol,
                pd_ad            = pd_ad,
                pd_fd            = pd_fd,
            )

            # Accumulation.
            ΣJ′WJ += J' * (W .* J)
            ΣJ′Wb += J' * (W .* b)
            σ_i   += dot(b, W .* b)
            σp_i  += dot(b[SOneTo(3)], b[SOneTo(3)])
            σv_i  += dot(b[StaticArrays.SUnitRange(4, 6)], b[StaticArrays.SUnitRange(4, 6)])
        end

        # Normalize and compute the RMS errors.
        σ_i  = √(σ_i / num_measurements)
        σp_i = √(σp_i / num_measurements)
        σv_i = √(σv_i / num_measurements)

        # Update the estimate.
        δx = ΣJ′WJ \ ΣJ′Wb

        # Limit the correction to avoid divergence. Each component is limited against the
        # norm of its position or velocity block instead of its own magnitude. Otherwise, a
        # component that is exactly zero, e.g. the z-position of an equatorial initial
        # guess, would have its correction clamped to zero, freezing the component at zero
        # for all iterations, and a component much smaller than its block norm would only
        # be able to change by a small fraction of its magnitude per iteration, slowing
        # down the convergence.
        norm_r₁ = norm(x₁[SOneTo(3)])
        norm_v₁ = norm(x₁[StaticArrays.SUnitRange(4, 6)])

        for i in 1:num_states
            threshold = T(0.1)
            scale = i <= 3 ? norm_r₁ : norm_v₁

            if abs(δx[i]) > threshold * scale
                δx = setindex(δx, threshold * scale * sign(δx[i]), i)
            end
        end

        x₂ = x₁ + δx

        # We cannot compute the RMSE variation in the first iteration.
        if it == 1
            verbose && _fit_print_progress(
                has_color,
                @sprintf(
                    "%10d %20g %20g %20g %20s", it, σp_i / 1000, σv_i / 1000, σ_i, "---"
                )
            )

        else
            # Compute the RMSE variation.
            Δσ = (σ_i - σ_i_₁) / σ_i_₁

            verbose && _fit_print_progress(
                has_color,
                @sprintf(
                    "%10d %20g %20g %20g %20g %%",
                    it,
                    σp_i / 1000,
                    σv_i / 1000,
                    σ_i,
                    100 * Δσ
                )
            )

            # Check if the RMSE is increasing.
            if σ_i < σ_i_₁
                Δd = 0
            else
                Δd += 1
            end

            # If the RMSE increased by three iterations and its value is higher than 5e11,
            # we abort because the iterations are diverging.
            ((Δd ≥ 3) && (σ_i > 5e11)) && throw(MeanElementsFitDivergenceError(it, σ_i))

            # Check if the condition to stop has been reached.
            if (abs(Δσ) < rtol) || (σ_i < atol)
                converged = true
                break
            end
        end

        σ_i_₁ = σ_i
    end

    verbose && println()

    # Obtain the mean elements, which store the mean anomaly.
    orb = convert(
        KeplerianElements{MeanAnomaly, Tepoch, T},
        rv_to_kepler(x₂[SOneTo(3)], x₂[StaticArrays.SUnitRange(4, 6)], epoch),
    )

    # Update the epoch of the fitted mean elements to match the desired one.
    if abs(epoch - desired_epoch) > 0.001 / 86400
        verbose && _fit_print_action(
            has_color,
            "Updating the epoch of the fitted mean elements to match the desired one.",
        )
        orb = _update_mean_elements_epoch!(pd, orb, desired_epoch)
    end

    # Initialize the propagator with the mean elements.
    _mean_elements_init!(pd, orb)

    # Compute the final covariance.
    P = pinv(ΣJ′WJ)

    # Assemble the statistics of the fit.
    stats = (;
        converged,
        iterations,
        position_rmse = σp_i,
        velocity_rmse = σv_i,
        total_rmse    = σ_i,
    )

    return orb, P, stats
end

# == Helpers ===============================================================================

"""
    _fit_print_action(has_color::Bool, msg::AbstractString) -> Nothing

Print to `stdout` the action message `msg` of the fitting algorithm, prefixed by an
`ACTION:` tag, which is highlighted if `has_color` is `true`.
"""
# The helper is not inlined so that its allocation sites, which are only reachable when the
# algorithm is verbose, are counted once regardless of the number of call sites.
@noinline function _fit_print_action(has_color::Bool, msg::AbstractString)
    println(_fit_decorated(_FIT_ACTION_TAG, has_color), "   ", msg)
    return nothing
end

"""
    _fit_print_header(has_color::Bool) -> Nothing

Print to `stdout` the header of the progress table of the fitting algorithm, followed by an
empty line that the first progress line overwrites. The header is decorated if `has_color`
is `true`.
"""
# The helper is not inlined so that its allocation sites, which are only reachable when the
# algorithm is verbose, are counted once regardless of the number of call sites.
@noinline function _fit_print_header(has_color::Bool)
    print(
        "          ",
        _fit_decorated(_FIT_HEADER, has_color),
        "\n          ",
        _fit_decorated(_FIT_UNITS, has_color),
        "\n\n",
    )
    return nothing
end

"""
    _fit_print_progress(has_color::Bool, msg::AbstractString) -> Nothing

Print to `stdout` the progress line `msg` of the fitting algorithm, prefixed by a
`PROGRESS:` tag, which is highlighted if `has_color` is `true`. The previous line is erased
first, so consecutive calls update the same terminal line.
"""
# The helper is not inlined so that its allocation sites, which are only reachable when the
# algorithm is verbose, are counted once regardless of the number of call sites.
@noinline function _fit_print_progress(has_color::Bool, msg::AbstractString)
    print("\x1b[A\x1b[2K\r", _fit_decorated(_FIT_PROGRESS_TAG, has_color), " ", msg, "\n")
    return nothing
end

"""
    _fit_decorated(versions::Tuple{String, String}, has_color::Bool) -> String

Return the colored version of a string printed by the fitting algorithm, stored as the
second element of `versions`, if `has_color` is `true`. Otherwise, return the plain version
stored as its first element.
"""
function _fit_decorated(versions::Tuple{String, String}, has_color::Bool)
    return versions[has_color ? 2 : 1]
end

"""
    _julian_day(epoch::Union{Number, DateTime}) -> Number

Return the `epoch` as a Julian Day, converting it if it is a `DateTime` [UTC].
"""
_julian_day(epoch::Number)   = epoch
_julian_day(epoch::DateTime) = datetime2julian(epoch)

"""
    _dual_type(::Type{T}) where {T <: Number} -> Type

Return the `ForwardDiff.Dual` type with six partial derivatives used to differentiate a
propagation with element type `T` with respect to the six components of the state vector.
"""
function _dual_type(::Type{T}) where {T <: Number}
    return ForwardDiff.Dual{ForwardDiff.Tag{Nothing, T}, T, 6}
end

"""
    _create_fd_propagator(
        pd::AbstractMeanElementsPropagator{Tepoch, T}
    ) where {Tepoch <: Number, T <: Number} -> typeof(pd)

Create an uninitialized propagator with the same constants as `pd` that the
finite-difference Jacobian uses as scratch space, so it does not clobber the propagator kept
initialized by the fitting loop.
"""
function _create_fd_propagator(
    pd::AbstractMeanElementsPropagator{Tepoch, T}
) where {Tepoch <: Number, T <: Number}
    return _similar_propagator(pd, T)
end

"""
    _create_ad_propagator(
        pd::AbstractMeanElementsPropagator{Tepoch, T}
    ) where {Tepoch <: Number, T <: Number} -> AbstractMeanElementsPropagator{Tepoch, D}

Create an uninitialized propagator with the same constants as `pd` but with the element type
`D = _dual_type(T)`, which the ForwardDiff Jacobian uses to propagate dual numbers.
"""
function _create_ad_propagator(
    pd::AbstractMeanElementsPropagator{Tepoch, T}
) where {Tepoch <: Number, T <: Number}
    return _similar_propagator(pd, _dual_type(T))
end

"""
    _mean_elements_epoch(
        pd::AbstractMeanElementsPropagator{Tepoch, T}
    ) where {Tepoch <: Number, T <: Number} -> Tepoch

Return the epoch [Julian Day] of the initial mean elements stored in `pd`.
"""
function _mean_elements_epoch(pd::AbstractMeanElementsPropagator)
    return _mean_elements_propagator(pd).orb₀.epoch
end

"""
    _update_mean_elements_epoch!(
        pd::AbstractMeanElementsPropagator{Tepoch, T},
        orb::KeplerianElements,
        new_epoch::Number
    ) where {Tepoch <: Number, T <: Number} -> KeplerianElements{MeanAnomaly, Tepoch, T}

Update the epoch of the mean elements `orb` [SI units] to `new_epoch` [Julian Day] by
propagating them with `pd`, and initialize `pd` with the returned elements. The osculating
propagators are initialized without computing the osculating elements at the new epoch; the
public functions `update_*_mean_elements_epoch!` complete that step.
"""
function _update_mean_elements_epoch!(
    pd::AbstractMeanElementsPropagator, orb::KeplerianElements, new_epoch::Number
)
    # First, we need to initialize the propagator with the mean elements.
    _mean_elements_init!(pd, orb)

    # Now, we just need to propagate the orbit to the desired instant and obtain the mean
    # elements from the propagator structure.
    Δt = (new_epoch - _mean_elements_epoch(pd)) * 86400
    _mean_elements_propagate!(pd, Δt)
    orbk = _mean_elements_propagator(pd).orbk

    # Finally, we initialize the propagator with the new set of mean elements.
    _mean_elements_init!(pd, orbk)

    return orbk
end

"""
    _mean_elements_jacobian(
        ::FiniteDiffJacobian,
        pd::AbstractMeanElementsPropagator{Tepoch, T},
        Δt::Number,
        x₁::SVector{6, T},
        y₁::SVector{6, T};
        kwargs...
    ) where {Tepoch <: Number, T <: Number} -> SMatrix{6, 6, T}

    _mean_elements_jacobian(
        ::ForwardDiffJacobian,
        pd::AbstractMeanElementsPropagator{Tepoch, T},
        Δt::Number,
        x₁::SVector{6, T},
        y₁::SVector{6, T};
        kwargs...
    ) where {Tepoch <: Number, T <: Number} -> SMatrix{6, 6, T}

Compute the Jacobian of the state vector `y₁` [m; m / s] propagated by `Δt` [s] with the
propagator `pd` with respect to the initial state vector `x₁` [m; m / s] at the epoch of the
initial mean elements in `pd`. The first argument selects the method: finite differences,
which perturbs each component of `x₁` and propagates the orbit again, or forward-mode
automatic differentiation with **ForwardDiff.jl**, which propagates dual numbers once.

!!! note

    The finite-difference method reinitializes `pd` for each perturbed state unless the
    keyword `pd_fd` is provided. The ForwardDiff method never modifies `pd`.

# Keywords

- `perturbation::Number`: Initial relative state perturbation to compute the finite
    differences. Only used with `FiniteDiffJacobian`.
    (**Default**: `T(1e-3)`)
- `perturbation_tol::Number`: Tolerance to accept the perturbation. If the computed
    perturbation is lower than `perturbation_tol`, we increase it until its absolute value
    is higher than `perturbation_tol`. Only used with `FiniteDiffJacobian`.
    (**Default**: `T(1e-7)`)
- `pd_ad::Union{Nothing, AbstractMeanElementsPropagator}`: Propagator with the dual element
    type used by the ForwardDiff method, as created by `_create_ad_propagator`. If it is
    `nothing`, a new one is allocated.
    (**Default**: `nothing`)
- `pd_fd::Union{Nothing, AbstractMeanElementsPropagator}`: Scratch propagator overwritten by
    the finite-difference method, as created by `_create_fd_propagator`. If it is `nothing`,
    `pd` itself is used and left initialized with the last perturbed state.
    (**Default**: `nothing`)
"""
function _mean_elements_jacobian(
    ::FiniteDiffJacobian,
    pd::AbstractMeanElementsPropagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{6, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7),
    pd_ad::Union{Nothing, AbstractMeanElementsPropagator} = nothing,
    pd_fd::Union{Nothing, AbstractMeanElementsPropagator} = nothing,
) where {Tepoch <: Number, T <: Number}
    # The perturbed propagations below overwrite the propagator. Use the scratch propagator
    # when the caller provides one, so it can keep `pd` initialized across the measurements
    # instead of reinitializing it for each one.
    epoch = _mean_elements_epoch(pd)
    fd    = isnothing(pd_fd) ? pd : pd_fd

    J  = MMatrix{6, 6, T}(undef)
    x₂ = x₁

    # Convert the perturbation parameters to the element type once. Otherwise, `ϵ` would be
    # assigned values of two different types, and the Jacobian columns would be computed in
    # the promoted type instead of in `T`.
    ϵ₀    = T(perturbation)
    ϵ_tol = T(perturbation_tol)

    # The loop below only accesses fixed indices of 6-element static arrays. Hence, we can
    # skip the bounds checking.
    @inbounds for j in 1:6
        α = x₂[j]
        ϵ = α * ϵ₀

        for _ in 1:5
            abs(ϵ) > ϵ_tol && break
            ϵ *= T(1.4)
        end

        if abs(ϵ) < ϵ_tol
            ϵ = signbit(α) ? -ϵ_tol : ϵ_tol
        end

        α += ϵ
        x₂ = setindex(x₂, α, j)

        orb = rv_to_kepler(x₂[SOneTo(3)], x₂[StaticArrays.SUnitRange(4, 6)], epoch)
        _mean_elements_init!(fd, orb)
        r_i, v_i = _mean_elements_propagate!(fd, Δt)
        y₂ = vcat(r_i, v_i)

        J[:, j] .= (y₂ .- y₁) ./ ϵ
        x₂ = setindex(x₂, x₁[j], j)
    end

    return SMatrix{6, 6, T}(J)
end

function _mean_elements_jacobian(
    ::ForwardDiffJacobian,
    pd::AbstractMeanElementsPropagator{Tepoch, T},
    Δt::Number,
    x₁::SVector{6, T},
    y₁::SVector{6, T};
    perturbation::Number = T(1e-3),
    perturbation_tol::Number = T(1e-7),
    pd_ad::Union{Nothing, AbstractMeanElementsPropagator} = nothing,
    pd_fd::Union{Nothing, AbstractMeanElementsPropagator} = nothing,
) where {Tepoch <: Number, T <: Number}
    epoch = _mean_elements_epoch(pd)
    ad    = isnothing(pd_ad) ? _create_ad_propagator(pd) : pd_ad

    # Seed the dual numbers so that the k-th partial derivative is taken with respect to the
    # k-th component of the state vector.
    D      = _dual_type(T)
    seeds  = ntuple(i -> ForwardDiff.Partials(ntuple(j -> T(i == j), Val(6))), Val(6))
    x_dual = SVector{6, D}(ntuple(i -> D(x₁[i], seeds[i]), Val(6)))

    orb = rv_to_kepler(x_dual[SOneTo(3)], x_dual[StaticArrays.SUnitRange(4, 6)], epoch)
    _mean_elements_init!(ad, orb)
    r_i, v_i = _mean_elements_propagate!(ad, Δt)
    y_dual   = vcat(r_i, v_i)

    # Assemble the Jacobian column by column from the partial derivatives.
    return SMatrix{6, 6, T}(
        ntuple(k -> ForwardDiff.partials(y_dual[mod1(k, 6)], cld(k, 6)), Val(36))
    )
end
