## Description #############################################################################
#
#   Tests related to the orbit propagator API.
#
############################################################################################

struct DummyPropagator{Tepoch, T} <: OrbitPropagator{Tepoch, T} end

@testset "Broadcast" verbose = true begin
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)
    jd₁ = date_to_jd(2023, 1, 5, 0, 0, 0)

    # == Float64 ===========================================================================

    @testset "Float64" begin
        T = Float64

        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |> deg2rad,
            T(100)  |> deg2rad,
            T(200)  |> deg2rad,
            T(45)   |> deg2rad
        )

        orbp = Propagators.init(Val(:J2), orb; j2c = j2c_egm2008)
        ret  = Propagators.propagate!.(orbp, 1:1:100)

        @test length(ret) == 100
        @test ret isa Vector{Tuple{SVector{3, Float64}, SVector{3, Float64}}}

        for k in 1:100
            r, v = Propagators.propagate!(orbp, k)
            @test r ≈ ret[k][1]
            @test v ≈ ret[k][2]
        end
    end

    # == Float32 ===========================================================================

    @testset "Float32" begin
        T = Float32

        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |> deg2rad,
            T(100)  |> deg2rad,
            T(200)  |> deg2rad,
            T(45)   |> deg2rad
        )

        orbp = Propagators.init(Val(:J2), orb; j2c = j2c_egm2008_f32)
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

@testset "DateTime Support" verbose = true begin
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

    for (T, j2c) in ((Float64, j2c_egm2008), (Float32, j2c_egm2008_f32))
        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |> deg2rad,
            T(100)  |> deg2rad,
            T(200)  |> deg2rad,
            T(45)   |> deg2rad
        )

        orbp_ref = Propagators.init(Val(:J2), orb; j2c = j2c)

        # == propagate_to_epoch ============================================================

        r, v, orbp = Propagators.propagate_to_epoch(
            Val(:J2),
            DateTime("2024-01-01"),
            orb;
            j2c = j2c
        )

        r_ref, v_ref = Propagators.propagate_to_epoch!(orbp_ref, date_to_jd(2024, 1, 1))

        @test orbp isa typeof(orbp_ref)
        @test r isa SVector{3, T}
        @test v isa SVector{3, T}
        @test r ≈ r_ref
        @test v ≈ v_ref
        @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

        # == propagate_to_epoch! ===========================================================

        orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

        r, v = Propagators.propagate_to_epoch!(orbp, DateTime("2024-01-01"))

        r_ref, v_ref = Propagators.propagate_to_epoch!(orbp_ref, date_to_jd(2024, 1, 1))

        @test r isa SVector{3, T}
        @test v isa SVector{3, T}
        @test r ≈ r_ref
        @test v ≈ v_ref
        @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)
    end
end

@testset "Multi-thread Propagation" verbose = true begin
    for (T, j2c) in ((Float64, j2c_egm2008), (Float32, j2c_egm2008_f32))
        @testset "$T" begin
            jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)
            jd₁ = date_to_jd(2023, 1, 5, 0, 0, 0)

            orb = KeplerianElements(
                jd₀,
                T(8000e3),
                T(0.015),
                T(28.5) |> deg2rad,
                T(100)  |> deg2rad,
                T(200)  |> deg2rad,
                T(45)   |> deg2rad
            )

            orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

            # == propagate! ================================================================

            ret  = Propagators.propagate!.(orbp, 1:1:100)
            r, v = Propagators.propagate!(orbp, 1:1:100)

            @test length(r) == 100
            @test length(v) == 100
            @test r isa Vector{SVector{3, T}}
            @test v isa Vector{SVector{3, T}}
            @test Propagators.last_instant(orbp) == 100.0

            for k in 1:100
                @test ret[k][1] == r[k]
                @test ret[k][2] == v[k]
            end

            # == propagate_to_epoch! =======================================================

            vjd  = collect(jd₀:0.1:jd₁)
            ret  = Propagators.propagate_to_epoch!.(orbp, vjd)
            r, v = Propagators.propagate_to_epoch!(orbp, vjd)

            @test length(r) == 41
            @test length(v) == 41
            @test r isa Vector{SVector{3, T}}
            @test v isa Vector{SVector{3, T}}
            @test Propagators.last_instant(orbp) == 345600.0

            for k in 1:41
                @test ret[k][1] == r[k]
                @test ret[k][2] == v[k]
            end

            # -- DateTime Support ----------------------------------------------------------

            vdt  = [DateTime(2024, 1, i) for i in 1:30]
            ret  = Propagators.propagate_to_epoch!.(orbp, vdt)
            r, v = Propagators.propagate_to_epoch!(orbp, vdt)

            @test length(r) == 30
            @test length(v) == 30
            @test r isa Vector{SVector{3, T}}
            @test v isa Vector{SVector{3, T}}
            @test Propagators.last_instant(orbp) == 3.40416e7

            for k in 1:30
                @test ret[k][1] == r[k]
                @test ret[k][2] == v[k]
            end

            # == Simultaneous Initialization and Propagation ===============================

            ret = Propagators.propagate!.(orbp, 1:1:100)
            r, v, orbp = Propagators.propagate(Val(:J2), 1:1:100, orb; j2c = j2c)

            @test length(r) == 100
            @test length(v) == 100
            @test r isa Vector{SVector{3, T}}
            @test v isa Vector{SVector{3, T}}
            @test orbp isa OrbitPropagatorJ2{Float64, T}
            @test Propagators.last_instant(orbp) == 100.0

            for k in 1:100
                @test ret[k][1] == r[k]
                @test ret[k][2] == v[k]
            end

        end
    end
end

@testset "Show" verbose = true begin
    T = Float64

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

    # == J2 Orbit Propagator ===============================================================

    orbp = Propagators.init(Val(:J2), orb)

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

    # == J4 Orbit Propagator ===============================================================

    orbp = Propagators.init(Val(:J4), orb)

    expected = "J4 Orbit Propagator (Epoch = $(string(dt₀)), Δt = 0.0 s)"
    result = sprint(show, orbp)
    @test result == expected

    expected = """
        OrbitPropagatorJ4{Float64, Float64}:
           Propagator name : J4 Orbit Propagator
          Propagator epoch : $(string(dt₀))
          Last propagation : $(string(dt₀))"""
    result = sprint(show, MIME("text/plain"), orbp)

    @test result == expected
end

@testset "Default Functions in the API" begin
    orbp = DummyPropagator{Float64, Float64}()
    @test Propagators.name(orbp) == "DummyPropagator{Float64, Float64}"
    @test Propagators.mean_elements(orbp) === nothing
end
