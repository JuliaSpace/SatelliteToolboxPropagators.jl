## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

@testset "Aqua.jl" begin
    Aqua.test_all(
        SatelliteToolboxPropagators;
        ambiguities  = (recursive = false),
        deps_compat  = (check_extras = false),
    )
end

if VERSION >= v"1.12"
    @warn "JET.jl test skipped on Julia 1.12+ due to MethodTableView incompatibility"
else
    @testset "JET Testing" begin
        rep = JET.test_package(SatelliteToolboxPropagators; toplevel_logger=nothing, target_modules=(@__MODULE__,))
    end
end

@testset "Allocation Check" begin

    # == Two-Body ==========================================================================
    if VERSION >= v"1.12"
        @warn "Allocation Check skipped on Julia 1.12+ as it is falsely flagging rem2pi internals"
    else
        @testset "twobody_init!" begin
            @test length(
                check_allocs(
                    (tbd, orb₀) -> twobody_init!(tbd, orb₀),
                    (TwoBodyPropagator{Float64, Float64}, KeplerianElements{Float64, Float64})
                )
            ) == 0
        end

        @testset "twobody!" begin
            @test length(
                check_allocs(
                    (tbd, t) -> twobody!(tbd, t),
                    (TwoBodyPropagator{Float64, Float64}, Float64)
                )
            ) == 0
        end

        # == J2 ================================================================================

        @testset "j2_init!" begin
            @test length(
                check_allocs(
                    (j2d, orb₀) -> j2_init!(j2d, orb₀),
                    (J2Propagator{Float64, Float64}, KeplerianElements{Float64, Float64})
                )
            ) == 0
        end

        @testset "j2!" begin
            @test length(
                check_allocs(
                    (j2d, t) -> j2!(j2d, t),
                    (J2Propagator{Float64, Float64}, Float64)
                )
            ) == 0
        end

        # == J2 Osculating =====================================================================

        @testset "j2osc_init!" begin
            @test length(
                check_allocs(
                    (j2oscd, orb₀) -> j2osc_init!(j2oscd, orb₀),
                    (J2OsculatingPropagator{Float64, Float64}, KeplerianElements{Float64, Float64})
                )
            ) == 0
        end

        @testset "j2osc!" begin
            @test length(
                check_allocs(
                    (j2oscd, t) -> j2osc!(j2oscd, t),
                    (J2OsculatingPropagator{Float64, Float64}, Float64)
                )
            ) == 0
        end

        # == J4 ================================================================================

        @testset "j4_init!" begin
            @test length(
                check_allocs(
                    (j4d, orb₀) -> j4_init!(j4d, orb₀),
                    (J4Propagator{Float64, Float64}, KeplerianElements{Float64, Float64})
                )
            ) == 0
        end

        @testset "j4!" begin
            @test length(
                check_allocs(
                    (j4d, t) -> j4!(j4d, t),
                    (J4Propagator{Float64, Float64}, Float64)
                )
            ) == 0
        end

        # == J4 Osculating =====================================================================

        @testset "j4osc_init!" begin
            @test length(
                check_allocs(
                    (j4oscd, orb₀) -> j4osc_init!(j4oscd, orb₀),
                    (J4OsculatingPropagator{Float64, Float64}, KeplerianElements{Float64, Float64})
                )
            ) == 0
        end

        @testset "j4osc!" begin
            @test length(
                check_allocs(
                    (j4oscd, t) -> j4osc!(j4oscd, t),
                    (J4OsculatingPropagator{Float64, Float64}, Float64)
                )
            ) == 0
        end
    end
end
