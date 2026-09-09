## Description #############################################################################
#
# Functions to print the structures of the orbit propagators implemented in this package.
#
# The representations follow the layout of SatelliteToolboxBase.jl: the compact form prints
# the type with its parameters and the epoch, whereas the rich form is a tree with the
# sections holding the initial mean elements and their epoch, the secular rates, the
# constants, and the last propagation instant. The `OrbitPropagator` wrappers of the API
# print the same tree under their own header in `src/api/Propagators.jl`.
#
############################################################################################

############################################################################################
#                                      Compact Format                                      #
############################################################################################

function Base.show(io::IO, pd::PropagatorData)
    name = type_name(pd)

    if !_is_initialized(pd)
        print(io, name, " (not initialized)")
        return nothing
    end

    print_compact(io, name, _initial_epoch(pd))

    return nothing
end

############################################################################################
#                                       Rich Format                                        #
############################################################################################

function Base.show(io::IO, ::MIME"text/plain", pd::PropagatorData)
    print_tree(io, type_name(pd), pd)
    return nothing
end

# The body of the rich representation is overloaded so that the API wrappers can print it
# under their own header.
function print_tree_body(io::IO, pd::PropagatorData)
    if !_is_initialized(pd)
        print_status(io, "not initialized")
        return nothing
    end

    sections = _sections(pd)
    push!(sections, "Propagation" => _propagation_fields(_last_instant(pd), "s"))

    print_tree_body(io, PrintedField[], sections)

    return nothing
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _constants(
        pd::PropagatorData
    ) -> Union{J2PropagatorConstants, J4PropagatorConstants, Number}

Return the constants of the initialized propagator structure `pd`: the structure with the
gravitational constants of the J2 and J4 propagators, or the standard gravitational
parameter [m³ / s²] of the two-body propagator.
"""
_constants(pd::J2Propagator)           = pd.j2c
_constants(pd::J2OsculatingPropagator) = _constants(pd.j2d)
_constants(pd::J4Propagator)           = pd.j4c
_constants(pd::J4OsculatingPropagator) = _constants(pd.j4d)
_constants(pd::TwoBodyPropagator)      = pd.μ

"""
    _constants_fields(constants) -> Vector{PrintedField}

Return the fields that print the propagator `constants`, which are a
`J2PropagatorConstants`, a `J4PropagatorConstants`, or the standard gravitational parameter
[m³ / s²] of the two-body propagator. The equatorial radius is printed in kilometers.
"""
function _constants_fields(j2c::J2PropagatorConstants)
    return PrintedField[
        ("R₀", format_value(j2c.R0 / 1000), "km"),
        ("μm", format_value(j2c.μm),        "rad/s"),
        ("J₂", format_value(j2c.J2),        ""),
    ]
end

function _constants_fields(j4c::J4PropagatorConstants)
    return PrintedField[
        ("R₀", format_value(j4c.R0 / 1000), "km"),
        ("μm", format_value(j4c.μm),        "rad/s"),
        ("J₂", format_value(j4c.J2),        ""),
        ("J₄", format_value(j4c.J4),        ""),
    ]
end

function _constants_fields(μ::Number)
    return PrintedField[("μ", format_value(μ), "m³/s²")]
end

"""
    _initial_elements(pd::PropagatorData) -> KeplerianElements{MeanAnomaly}

Return the initial mean elements stored in the initialized propagator structure `pd`. The
osculating propagators return the elements of the propagator they wrap.
"""
_initial_elements(pd::J2Propagator)           = pd.orb₀
_initial_elements(pd::J2OsculatingPropagator) = _initial_elements(pd.j2d)
_initial_elements(pd::J4Propagator)           = pd.orb₀
_initial_elements(pd::J4OsculatingPropagator) = _initial_elements(pd.j4d)
_initial_elements(pd::TwoBodyPropagator)      = pd.orb₀

"""
    _initial_epoch(pd::PropagatorData) -> Number

Return the epoch [Julian Day] of the initial mean elements stored in the initialized
propagator structure `pd`.
"""
_initial_epoch(pd::PropagatorData) = _initial_elements(pd).epoch

"""
    _propagation_fields(Δt::Number, unit::String) -> Vector{PrintedField}

Return the fields of the section that prints the last propagation instant `Δt`, measured
from the epoch of the initial mean elements in the `unit` given as a string.
"""
function _propagation_fields(Δt::Number, unit::String)
    return PrintedField[("Last Instant", format_value(Δt), unit)]
end

"""
    _mean_elements_fields(orb₀::KeplerianElements{MeanAnomaly}) -> Vector{PrintedField}

Return the fields that print the initial mean elements `orb₀` and their epoch in the rich
representation of a propagator structure. The semi-major axis is printed in kilometers and
the angles in degrees.
"""
function _mean_elements_fields(orb₀::KeplerianElements{MeanAnomaly})
    return PrintedField[
        ("Epoch",             epoch_string(orb₀.epoch),                          ""),
        ("Semi-Major Axis",   format_value(orb₀.semi_major_axis / 1000),         "km"),
        ("Eccentricity",      format_value(orb₀.eccentricity),                   ""),
        ("Inclination",       format_value(rad2deg(orb₀.inclination)),           "°"),
        ("RA of Asc. Node",   format_value(rad2deg(orb₀.raan)),                  "°"),
        ("Arg. of Periapsis", format_value(rad2deg(orb₀.argument_of_periapsis)), "°"),
        ("Mean Anomaly",      format_value(rad2deg(orb₀.anomaly)),               "°"),
    ]
end

"""
    _secular_rates_fields(n̄::Number[, ∂Ω::Number, ∂ω::Number]) -> Vector{PrintedField}

Return the fields that print the perturbed mean motion `n̄` [rad / s] and, if provided, the
RAAN rate `∂Ω` [rad / s] and the argument of periapsis rate `∂ω` [rad / s] in the rich
representation of a propagator structure. The mean motion is printed in revolutions per day
and the rates in degrees per day.
"""
function _secular_rates_fields(n̄::Number)
    return PrintedField[("Mean Motion", format_value(86400 * n̄ / 2π), "rev/day")]
end

function _secular_rates_fields(n̄::Number, ∂Ω::Number, ∂ω::Number)
    return PrintedField[
        ("Mean Motion",            format_value(86400 * n̄ / 2π),     "rev/day"),
        ("RAAN Rate",              format_value(86400 * rad2deg(∂Ω)), "°/day"),
        ("Arg. of Periapsis Rate", format_value(86400 * rad2deg(∂ω)), "°/day"),
    ]
end

"""
    _sections(pd::PropagatorData) -> Vector{PrintedSection}

Return the sections of the rich representation of the initialized propagator structure
`pd`: the initial mean elements with their epoch, the secular rates, and the constants. The
last propagation instant is appended by the caller. The osculating propagators print the
sections of the propagator they wrap, since the short-period corrections do not add any
parameter.
"""
function _sections(pd::Union{J2Propagator, J4Propagator})
    return PrintedSection[
        "Mean Elements" => _mean_elements_fields(pd.orb₀),
        "Secular Rates" => _secular_rates_fields(pd.n̄, pd.∂Ω, pd.∂ω),
        "Constants"     => _constants_fields(_constants(pd)),
    ]
end

function _sections(pd::TwoBodyPropagator)
    return PrintedSection[
        "Mean Elements" => _mean_elements_fields(pd.orb₀),
        "Secular Rates" => _secular_rates_fields(pd.n₀),
        "Constants"     => _constants_fields(pd.μ),
    ]
end

_sections(pd::J2OsculatingPropagator) = _sections(pd.j2d)
_sections(pd::J4OsculatingPropagator) = _sections(pd.j4d)
