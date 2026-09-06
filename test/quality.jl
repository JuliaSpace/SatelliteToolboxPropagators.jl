## Description #############################################################################
#
# Quality tests using Aqua.jl and JET.jl.
#
############################################################################################

@testset "Aqua.jl" begin
    # The method ambiguities are not checked recursively because the dependencies define
    # ambiguous methods that this package cannot fix. The compat bounds of the extras are
    # not checked because the test dependencies are bounded in the same `[compat]` section.
    Aqua.test_all(
        SatelliteToolboxPropagators;
        ambiguities = (recursive = false),
        deps_compat = (check_extras = false),
    )
end

@testset "JET.jl" begin
    JET.test_package(
        SatelliteToolboxPropagators;
        toplevel_logger = nothing,
        target_modules = (@__MODULE__,),
    )
end
