# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests related to the orbit propagator API.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

@testset "Broadcast" begin
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)
    jd₁ = date_to_jd(2023, 1, 5, 0, 0, 0)

    # Float64
    # ======================================================================================

    let T = Float64
        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |> deg2rad,
            T(100)  |> deg2rad,
            T(200)  |> deg2rad,
            T(45)   |> deg2rad
        )

        orbp = Propagators.init(Val(:J2), orb; j2c = j2c_egm08)
        ret  = Propagators.propagate!.(orbp, 1:1:100)

        @test length(ret) == 100
        @test ret isa Vector{Tuple{SVector{3, Float64}, SVector{3, Float64}}}

        for k in 1:100
            r, v = Propagators.propagate!(orbp, k)
            @test r ≈ ret[k][1]
            @test v ≈ ret[k][2]
        end
    end

    # Float32
    # ======================================================================================

    let T = Float32
        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |> deg2rad,
            T(100)  |> deg2rad,
            T(200)  |> deg2rad,
            T(45)   |> deg2rad
        )

        orbp = Propagators.init(Val(:J2), orb; j2c = j2c_egm08_f32)
        ret  = Propagators.propagate!.(orbp, 1:1:100)

        @test length(ret) == 100
        @test ret isa Vector{Tuple{SVector{3, Float32}, SVector{3, Float32}}}

        for k in 1:100
            r, v = Propagators.propagate!(orbp, k)
            @test r ≈ ret[k][1]
            @test v ≈ ret[k][2]
        end
    end
end

@testset "Show" begin
    let T = Float64
        jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)
        dt₀ = date_to_jd(2023, 1, 1, 0, 0, 0) |> julian2datetime
        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |> deg2rad,
            T(100)  |> deg2rad,
            T(200)  |> deg2rad,
            T(45)   |> deg2rad
        )

        orbp = Propagators.init(Val(:J2), orb; j2c = j2c_egm08)

        expected = "J2 Orbit Propagator (Epoch = $(string(dt₀)), Δt = 0.0 s)"
        result = sprint(show, orbp)
        @test result == expected

        expected = """
        OrbitPropagatorJ2{Float64, Float64}:
           Propagator name : J2 Orbit Propagator
          Propagator epoch : $(string(dt₀))
          Last propagation : $(string(dt₀))"""
        result = sprint(show, MIME("text/plain"), orbp)

        @test result == expected
    end
end
