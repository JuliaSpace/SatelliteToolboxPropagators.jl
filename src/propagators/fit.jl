## Description #############################################################################
#
# Generic least-square algorithm to fit a set of mean elements from osculating state
# vectors.
#
# The four analytical propagators share the same algorithm. They only differ in how the
# propagator is initialized, propagated, and differentiated, which the small set of methods
# below provides for each of them.
#
############################################################################################

# Union of the propagator structures that support fitting a set of mean elements from
# osculating state vectors using `_fit_mean_elements!`.
const AbstractMeanElementsPropagator{Tepoch, T} = Union{
    J2Propagator{Tepoch, T},
    J2OsculatingPropagator{Tepoch, T},
    J4Propagator{Tepoch, T},
    J4OsculatingPropagator{Tepoch, T},
}

# == Propagator-Specific Operations ========================================================

_propagator_name(::J2Propagator)           = "J2"
_propagator_name(::J2OsculatingPropagator) = "J2 osculating"
_propagator_name(::J4Propagator)           = "J4"
_propagator_name(::J4OsculatingPropagator) = "J4 osculating"

# The osculating propagators use the internal initialization, which skips the propagation to
# the initial instant, because the algorithm always propagates right afterwards.
_mean_elements_init!(pd::J2Propagator, orb)           = j2_init!(pd, orb)
_mean_elements_init!(pd::J2OsculatingPropagator, orb) = _j2osc_init!(pd, orb)
_mean_elements_init!(pd::J4Propagator, orb)           = j4_init!(pd, orb)
_mean_elements_init!(pd::J4OsculatingPropagator, orb) = _j4osc_init!(pd, orb)

_mean_elements_propagate!(pd::J2Propagator, Δt)           = j2!(pd, Δt)
_mean_elements_propagate!(pd::J2OsculatingPropagator, Δt) = j2osc!(pd, Δt)
_mean_elements_propagate!(pd::J4Propagator, Δt)           = j4!(pd, Δt)
_mean_elements_propagate!(pd::J4OsculatingPropagator, Δt) = j4osc!(pd, Δt)

_create_ad_propagator(pd::J2Propagator)           = _create_j2_ad_propagator(pd)
_create_ad_propagator(pd::J2OsculatingPropagator) = _create_j2osc_ad_propagator(pd)
_create_ad_propagator(pd::J4Propagator)           = _create_j4_ad_propagator(pd)
_create_ad_propagator(pd::J4OsculatingPropagator) = _create_j4osc_ad_propagator(pd)

_create_fd_propagator(pd::J2Propagator)           = _create_j2_fd_propagator(pd)
_create_fd_propagator(pd::J2OsculatingPropagator) = _create_j2osc_fd_propagator(pd)
_create_fd_propagator(pd::J4Propagator)           = _create_j4_fd_propagator(pd)
_create_fd_propagator(pd::J4OsculatingPropagator) = _create_j4osc_fd_propagator(pd)

_update_mean_elements_epoch!(pd::J2Propagator, orb, epoch) =
    update_j2_mean_elements_epoch!(pd, orb, epoch)
_update_mean_elements_epoch!(pd::J2OsculatingPropagator, orb, epoch) =
    update_j2osc_mean_elements_epoch!(pd, orb, epoch)
_update_mean_elements_epoch!(pd::J4Propagator, orb, epoch) =
    update_j4_mean_elements_epoch!(pd, orb, epoch)
_update_mean_elements_epoch!(pd::J4OsculatingPropagator, orb, epoch) =
    update_j4osc_mean_elements_epoch!(pd, orb, epoch)

_mean_elements_jacobian(m, pd::J2Propagator, args...; kwargs...) =
    _j2_jacobian(m, pd, args...; kwargs...)
_mean_elements_jacobian(m, pd::J2OsculatingPropagator, args...; kwargs...) =
    _j2osc_jacobian(m, pd, args...; kwargs...)
_mean_elements_jacobian(m, pd::J4Propagator, args...; kwargs...) =
    _j4_jacobian(m, pd, args...; kwargs...)
_mean_elements_jacobian(m, pd::J4OsculatingPropagator, args...; kwargs...) =
    _j4osc_jacobian(m, pd, args...; kwargs...)

# == Algorithm =============================================================================

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
    mean_elements_epoch::Number                      = vjd[end],
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
        epoch = Tepoch(mean_elements_epoch)

        # First, we need to update the mean elements to the desired epoch.
        verbose && _fit_print_action(
            "Updating the epoch of the initial mean elements guess to match the desired one.",
        )
        orb = _update_mean_elements_epoch!(pd, initial_guess, epoch)

        r_i, v_i = kepler_to_rv(orb)
        x₁ = SVector{6, T}(r_i[1], r_i[2], r_i[3], v_i[1], v_i[2], v_i[3])
    else
        # In this case, we must find the closest osculating vector to the desired epoch.
        id = firstindex(vjd)
        v  = abs(vjd[id] - mean_elements_epoch)

        for k in eachindex(vjd)
            vk = abs(vjd[k] - mean_elements_epoch)
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

    # Variable to store the last residue.
    σ_i_₁ = T(0)

    # Variable to store how many iterations the residue increased. This is used to account
    # for divergence.
    Δd = 0

    # Header.
    if verbose
        _fit_print_action(
            "Fitting the mean elements for the $(_propagator_name(pd)) propagator."
        )

        header = @sprintf(
            "%10s %20s %20s %20s %20s",
            "Iteration",
            "Position RMSE",
            "Velocity RMSE",
            "Total RMSE",
            "RMSE Variation"
        )
        units = @sprintf("%10s %20s %20s %20s %20s", "", "[km]", "[km / s]", "[ ]", "")

        println("          ", styled"{(foreground=yellow,weight=bold):$header}")
        println("          ", styled"{bold:$units}")
        println()
    end

    # We need a reference to the covariance inverse because we will invert it and return
    # after the iterations.
    ΣJ′WJ = @SMatrix zeros(T, num_states, num_states)

    pd_ad = jacobian_method isa ForwardDiffJacobian ? _create_ad_propagator(pd) : nothing

    # The finite-difference Jacobian overwrites the propagator it uses. Giving it its own
    # propagator lets us initialize `pd` once per iteration instead of once per
    # measurement.
    pd_fd = jacobian_method isa FiniteDiffJacobian ? _create_fd_propagator(pd) : nothing

    # Loop until the maximum allowed iteration.
    @inbounds @views for it in 1:max_iterations
        x₁ = x₂

        # Variables to store the summations to compute the least square fitting algorithm.
        ΣJ′WJ = @SMatrix zeros(T, num_states, num_states)
        ΣJ′Wb = @SVector zeros(T, num_states)

        # Variable to store the RMS errors in this iteration.
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
                @sprintf(
                    "%10d %20g %20g %20g %20s", it, σp_i / 1000, σv_i / 1000, σ_i, "---"
                )
            )

        else
            # Compute the RMSE variation.
            Δσ = (σ_i - σ_i_₁) / σ_i_₁

            verbose && _fit_print_progress(
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
            ((Δd ≥ 3) && (σ_i > 5e11)) && error("The iterations diverged!")

            # Check if the condition to stop has been reached.
            ((abs(Δσ) < rtol) || (σ_i < atol) || (it ≥ max_iterations)) && break
        end

        σ_i_₁ = σ_i
    end

    verbose && println()

    # Obtain the mean elements.
    orb = rv_to_kepler(x₂[SOneTo(3)], x₂[StaticArrays.SUnitRange(4, 6)], epoch)

    # Update the epoch of the fitted mean elements to match the desired one.
    if abs(epoch - mean_elements_epoch) > 0.001 / 86400
        verbose && _fit_print_action(
            "Updating the epoch of the fitted mean elements to match the desired one."
        )
        orb = convert(
            typeof(orb), _update_mean_elements_epoch!(pd, orb, mean_elements_epoch)
        )
    end

    # Initialize the propagator with the mean elements.
    _mean_elements_init!(pd, orb)

    # Compute the final covariance.
    P = pinv(ΣJ′WJ)

    # Return the mean elements and the covariance.
    return orb, P
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _fit_print_action(msg::AbstractString) -> Nothing

Print to `stdout` the action message `msg` of the fitting algorithm, prefixed by a
highlighted `ACTION:` tag. The decorations are only emitted if `stdout` supports colors.
"""
function _fit_print_action(msg::AbstractString)
    println(styled"{(foreground=yellow,weight=bold):ACTION:}   ", msg)
    return nothing
end

"""
    _fit_print_progress(msg::AbstractString) -> Nothing

Print to `stdout` the progress line `msg` of the fitting algorithm, prefixed by a highlighted
`PROGRESS:` tag. The previous line is erased first, so consecutive calls update the same
terminal line. The decorations are only emitted if `stdout` supports colors.
"""
function _fit_print_progress(msg::AbstractString)
    print("\x1b[A\x1b[2K\r", styled"{bold:PROGRESS:} ", msg, "\n")
    return nothing
end
