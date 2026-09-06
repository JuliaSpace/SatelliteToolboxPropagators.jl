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
    name = Propagators._type_name(pd)

    if !_is_initialized(pd)
        print(io, name, " (not initialized)")
        return nothing
    end

    SatelliteToolboxBase.print_compact(io, name, _initial_epoch(pd))

    return nothing
end

############################################################################################
#                                       Rich Format                                        #
############################################################################################

function Base.show(io::IO, ::MIME"text/plain", pd::PropagatorData)
    SatelliteToolboxBase.print_tree(io, Propagators._type_name(pd), pd)
    return nothing
end

# The body of the rich representation is overloaded so that the API wrappers can print it
# under their own header.
function SatelliteToolboxBase.print_tree_body(io::IO, pd::PropagatorData)
    if !_is_initialized(pd)
        fields   = SatelliteToolboxBase.PrintedField[("Status", "not initialized", "")]
        sections = SatelliteToolboxBase.PrintedSection[]
        SatelliteToolboxBase.print_tree_body(io, fields, sections)
        return nothing
    end

    sections = _sections(pd)
    push!(sections, "Propagation" => _propagation_fields(pd.Δt, "s"))

    SatelliteToolboxBase.print_tree_body(io, SatelliteToolboxBase.PrintedField[], sections)

    return nothing
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _constants_fields(constants) -> Vector{SatelliteToolboxBase.PrintedField}

Return the fields that print the propagator `constants`, which are a
`J2PropagatorConstants`, a `J4PropagatorConstants`, or the standard gravitational parameter
[m³ / s²] of the two-body propagator. The equatorial radius is printed in kilometers.
"""
function _constants_fields(j2c::J2PropagatorConstants)
    format_value = SatelliteToolboxBase.format_value

    return SatelliteToolboxBase.PrintedField[
        ("R₀", format_value(j2c.R0 / 1000), "km"),
        ("μm", format_value(j2c.μm),        "rad/s"),
        ("J₂", format_value(j2c.J2),        ""),
    ]
end

function _constants_fields(j4c::J4PropagatorConstants)
    format_value = SatelliteToolboxBase.format_value

    return SatelliteToolboxBase.PrintedField[
        ("R₀", format_value(j4c.R0 / 1000), "km"),
        ("μm", format_value(j4c.μm),        "rad/s"),
        ("J₂", format_value(j4c.J2),        ""),
        ("J₄", format_value(j4c.J4),        ""),
    ]
end

function _constants_fields(μ::Number)
    return SatelliteToolboxBase.PrintedField[
        ("μ", SatelliteToolboxBase.format_value(μ), "m³/s²"),
    ]
end

"""
    _initial_epoch(pd::PropagatorData) -> Number

Return the epoch [Julian Day] of the initial mean elements stored in the initialized
propagator structure `pd`.
"""
_initial_epoch(pd::J2Propagator)           = pd.orb₀.epoch
_initial_epoch(pd::J2OsculatingPropagator) = _initial_epoch(pd.j2d)
_initial_epoch(pd::J4Propagator)           = pd.orb₀.epoch
_initial_epoch(pd::J4OsculatingPropagator) = _initial_epoch(pd.j4d)
_initial_epoch(pd::TwoBodyPropagator)      = pd.orb₀.epoch

"""
    _propagation_fields(
        Δt::Number,
        unit::String
    ) -> Vector{SatelliteToolboxBase.PrintedField}

Return the fields of the section that prints the last propagation instant `Δt`, measured
from the epoch of the initial mean elements in the `unit` given as a string.
"""
function _propagation_fields(Δt::Number, unit::String)
    return SatelliteToolboxBase.PrintedField[
        ("Last Instant", SatelliteToolboxBase.format_value(Δt), unit),
    ]
end

"""
    _is_initialized(pd::PropagatorData) -> Bool

Return whether the propagator structure `pd` has been initialized. The empty constructors
set the field `Δt` to `NaN`, which is replaced by the initialization functions.
"""
_is_initialized(pd::PropagatorData) = !isnan(pd.Δt)

"""
    _mean_elements_fields(
        orb₀::KeplerianElements{MeanAnomaly}
    ) -> Vector{SatelliteToolboxBase.PrintedField}

Return the fields that print the initial mean elements `orb₀` and their epoch in the rich
representation of a propagator structure. The semi-major axis is printed in kilometers and
the angles in degrees.
"""
function _mean_elements_fields(orb₀::KeplerianElements{MeanAnomaly})
    format_value = SatelliteToolboxBase.format_value

    return SatelliteToolboxBase.PrintedField[
        ("Epoch",             SatelliteToolboxBase.epoch_string(orb₀.epoch),     ""),
        ("Semi-Major Axis",   format_value(orb₀.semi_major_axis / 1000),         "km"),
        ("Eccentricity",      format_value(orb₀.eccentricity),                   ""),
        ("Inclination",       format_value(rad2deg(orb₀.inclination)),           "°"),
        ("RA of Asc. Node",   format_value(rad2deg(orb₀.raan)),                  "°"),
        ("Arg. of Periapsis", format_value(rad2deg(orb₀.argument_of_periapsis)), "°"),
        ("Mean Anomaly",      format_value(rad2deg(orb₀.anomaly)),               "°"),
    ]
end

"""
    _secular_rates_fields(
        n̄::Number[, ∂Ω::Number, ∂ω::Number]
    ) -> Vector{SatelliteToolboxBase.PrintedField}

Return the fields that print the perturbed mean motion `n̄` [rad / s] and, if provided, the
RAAN rate `∂Ω` [rad / s] and the argument of pericenter rate `∂ω` [rad / s] in the rich
representation of a propagator structure. The mean motion is printed in revolutions per day
and the rates in degrees per day.
"""
function _secular_rates_fields(n̄::Number)
    return SatelliteToolboxBase.PrintedField[
        ("Mean Motion", SatelliteToolboxBase.format_value(86400 * n̄ / 2π), "rev/day"),
    ]
end

function _secular_rates_fields(n̄::Number, ∂Ω::Number, ∂ω::Number)
    format_value = SatelliteToolboxBase.format_value

    return SatelliteToolboxBase.PrintedField[
        ("Mean Motion",            format_value(86400 * n̄ / 2π),     "rev/day"),
        ("RAAN Rate",              format_value(86400 * rad2deg(∂Ω)), "°/day"),
        ("Arg. of Periapsis Rate", format_value(86400 * rad2deg(∂ω)), "°/day"),
    ]
end

"""
    _sections(pd::PropagatorData) -> Vector{SatelliteToolboxBase.PrintedSection}

Return the sections of the rich representation of the initialized propagator structure
`pd`: the initial mean elements with their epoch, the secular rates, and the constants. The
last propagation instant is appended by the caller. The osculating propagators print the
sections of the propagator they wrap, since the short-period corrections do not add any
parameter.
"""
function _sections(pd::J2Propagator)
    return SatelliteToolboxBase.PrintedSection[
        "Mean Elements" => _mean_elements_fields(pd.orb₀),
        "Secular Rates" => _secular_rates_fields(pd.n̄, pd.∂Ω, pd.∂ω),
        "Constants"     => _constants_fields(pd.j2c),
    ]
end

function _sections(pd::J4Propagator)
    return SatelliteToolboxBase.PrintedSection[
        "Mean Elements" => _mean_elements_fields(pd.orb₀),
        "Secular Rates" => _secular_rates_fields(pd.n̄, pd.∂Ω, pd.∂ω),
        "Constants"     => _constants_fields(pd.j4c),
    ]
end

function _sections(pd::TwoBodyPropagator)
    return SatelliteToolboxBase.PrintedSection[
        "Mean Elements" => _mean_elements_fields(pd.orb₀),
        "Secular Rates" => _secular_rates_fields(pd.n₀),
        "Constants"     => _constants_fields(pd.μ),
    ]
end

_sections(pd::J2OsculatingPropagator) = _sections(pd.j2d)
_sections(pd::J4OsculatingPropagator) = _sections(pd.j4d)
