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
            T(100) |> deg2rad,
            T(200) |> deg2rad,
            T(45) |> deg2rad,
        )

        orbp = Propagators.init(Val(:J2), orb; j2c = J2C_EGM2008)
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
            T(100) |> deg2rad,
            T(200) |> deg2rad,
            T(45) |> deg2rad,
        )

        orbp = Propagators.init(Val(:J2), orb; j2c = J2C_EGM2008_F32)
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

@testset "Dates Support" verbose = true begin
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

    for (T, j2c) in ((Float64, J2C_EGM2008), (Float32, J2C_EGM2008_F32))
        @testset "$T" begin
            orb = KeplerianElements(
                jd₀,
                T(8000e3),
                T(0.015),
                T(28.5) |> deg2rad,
                T(100) |> deg2rad,
                T(200) |> deg2rad,
                T(45) |> deg2rad,
            )

            orbp_ref = Propagators.init(Val(:J2), orb; j2c = j2c)

            # == propagate =================================================================

            r, v, orbp = Propagators.propagate(
                Val(:J2), Dates.Minute(1) + Dates.Second(1), orb; j2c = j2c
            )

            r_ref, v_ref = Propagators.propagate!(orbp_ref, 61)

            @test orbp isa typeof(orbp_ref)
            @test r isa SVector{3, T}
            @test v isa SVector{3, T}
            @test r ≈ r_ref
            @test v ≈ v_ref
            @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

            # == propagate! ================================================================

            orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

            r, v = Propagators.propagate!(orbp, Dates.Minute(1) + Dates.Second(1))

            r_ref, v_ref = Propagators.propagate!(orbp_ref, 61)

            @test orbp isa typeof(orbp_ref)
            @test r isa SVector{3, T}
            @test v isa SVector{3, T}
            @test r ≈ r_ref
            @test v ≈ v_ref
            @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

            # == propagate_to_epoch ========================================================

            r, v, orbp = Propagators.propagate_to_epoch(
                Val(:J2), DateTime("2024-01-01"), orb; j2c = j2c
            )

            r_ref, v_ref = Propagators.propagate_to_epoch!(orbp_ref, date_to_jd(2024, 1, 1))

            @test orbp isa typeof(orbp_ref)
            @test r isa SVector{3, T}
            @test v isa SVector{3, T}
            @test r ≈ r_ref
            @test v ≈ v_ref
            @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

            # == propagate_to_epoch! =======================================================

            orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

            r, v = Propagators.propagate_to_epoch!(orbp, DateTime("2024-01-01"))

            r_ref, v_ref = Propagators.propagate_to_epoch!(orbp_ref, date_to_jd(2024, 1, 1))

            @test r isa SVector{3, T}
            @test v isa SVector{3, T}
            @test r ≈ r_ref
            @test v ≈ v_ref
            @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

            # == step! =====================================================================

            orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

            r, v = Propagators.step!(orbp, Dates.Day(365))

            r_ref, v_ref = Propagators.propagate_to_epoch!(orbp_ref, date_to_jd(2024, 1, 1))

            @test r isa SVector{3, T}
            @test v isa SVector{3, T}
            @test r ≈ r_ref
            @test v ≈ v_ref
            @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)
        end
    end
end

@testset "Multi-thread Propagation" verbose = true begin
    for (T, j2c) in ((Float64, J2C_EGM2008), (Float32, J2C_EGM2008_F32))
        @testset "$T" begin
            jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)
            jd₁ = date_to_jd(2023, 1, 5, 0, 0, 0)

            orb = KeplerianElements(
                jd₀,
                T(8000e3),
                T(0.015),
                T(28.5) |> deg2rad,
                T(100) |> deg2rad,
                T(200) |> deg2rad,
                T(45) |> deg2rad,
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

            # -- Dates Support -------------------------------------------------------------

            vp   = [Dates.Second(i) for i in 1:1:100]
            ret  = Propagators.propagate!.(orbp, 1:1:100)
            r, v = Propagators.propagate!(orbp, vp)

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

            # -- Dates Support -------------------------------------------------------------

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

            vjd = collect(jd₀:0.1:jd₁)
            ret = Propagators.propagate_to_epoch!.(orbp, vjd)
            r, v, orbp = Propagators.propagate_to_epoch(Val(:J2), vjd, orb; j2c = j2c)

            @test length(r) == 41
            @test length(v) == 41
            @test r isa Vector{SVector{3, T}}
            @test v isa Vector{SVector{3, T}}
            @test orbp isa OrbitPropagatorJ2{Float64, T}
            @test Propagators.last_instant(orbp) == 345600.0

            for k in 1:41
                @test ret[k][1] == r[k]
                @test ret[k][2] == v[k]
            end

            # -- Dates Support -------------------------------------------------------------

            vp = [Dates.Second(i) for i in 1:1:100]
            ret = Propagators.propagate!.(orbp, 1:1:100)
            r, v, orbp = Propagators.propagate(Val(:J2), vp, orb; j2c = j2c)

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

            vdt = [DateTime(2024, 1, i) for i in 1:30]
            ret = Propagators.propagate_to_epoch!.(orbp, vdt)
            r, v, orbp = Propagators.propagate_to_epoch(Val(:J2), vdt, orb; j2c = j2c)

            @test length(r) == 30
            @test length(v) == 30
            @test r isa Vector{SVector{3, T}}
            @test v isa Vector{SVector{3, T}}
            @test orbp isa OrbitPropagatorJ2{Float64, T}
            @test Propagators.last_instant(orbp) == 3.40416e7

            for k in 1:30
                @test ret[k][1] == r[k]
                @test ret[k][2] == v[k]
            end
        end
    end
end

@testset "Sink Options" verbose = true begin
    @testset "Julian Day" verbose = true begin
        jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

        for (T, j2c) in ((Float64, J2C_EGM2008), (Float32, J2C_EGM2008_F32))
            @testset "$T" begin
                orb = KeplerianElements(
                    jd₀,
                    T(8000e3),
                    T(0.015),
                    T(28.5) |> deg2rad,
                    T(100) |> deg2rad,
                    T(200) |> deg2rad,
                    T(45) |> deg2rad,
                )

                orbp_ref = Propagators.init(Val(:J2), orb; j2c = j2c)

                # == propagate =============================================================

                r_ref, v_ref = Propagators.propagate!(orbp_ref, 61)

                r, v, orbp = Propagators.propagate(Tuple, Val(:J2), 61, orb; j2c = j2c)

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                sv, orbp = Propagators.propagate(
                    OrbitStateVector, Val(:J2), 61, orb; j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                vr_ref, vv_ref = Propagators.propagate!(orbp_ref, [60, 61])

                vr, vv, orbp = Propagators.propagate(
                    Tuple, Val(:J2), [60, 61], orb; j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test vr isa Vector{SVector{3, T}}
                @test vv isa Vector{SVector{3, T}}
                @test vr ≈ vr_ref
                @test vv ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                vsv, orbp = Propagators.propagate(
                    OrbitStateVector, Val(:J2), [60, 61], orb; j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == propagate! ============================================================

                r_ref, v_ref = Propagators.propagate!(orbp_ref, 61)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                r, v = Propagators.propagate!(orbp, 61, Tuple)

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                sv = Propagators.propagate!(orbp, 61, OrbitStateVector)

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                vr_ref, vv_ref = Propagators.propagate!(orbp_ref, [60, 61])

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                r, v = Propagators.propagate!(orbp, [60, 61], Tuple)

                @test orbp isa typeof(orbp_ref)
                @test vr isa Vector{SVector{3, T}}
                @test vv isa Vector{SVector{3, T}}
                @test vr ≈ vr_ref
                @test vv ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                vsv = Propagators.propagate!(orbp, [60, 61], OrbitStateVector)

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == propagate_to_epoch ====================================================

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref, date_to_jd(2024, 1, 1)
                )

                r, v, orbp = Propagators.propagate_to_epoch(
                    Tuple, Val(:J2), date_to_jd(2024, 1, 1), orb; j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                sv, orbp = Propagators.propagate_to_epoch(
                    OrbitStateVector, Val(:J2), date_to_jd(2024, 1, 1), orb; j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                vr_ref, vv_ref = Propagators.propagate_to_epoch!(
                    orbp_ref, [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)]
                )

                vr, vv, orbp = Propagators.propagate_to_epoch(
                    Tuple,
                    Val(:J2),
                    [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)],
                    orb;
                    j2c = j2c,
                )

                @test orbp isa typeof(orbp_ref)
                @test vr isa Vector{SVector{3, T}}
                @test vv isa Vector{SVector{3, T}}
                @test vr ≈ vr_ref
                @test vv ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                vsv, orbp = Propagators.propagate_to_epoch(
                    OrbitStateVector,
                    Val(:J2),
                    [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)],
                    orb;
                    j2c = j2c,
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == propagate_to_epoch! ===================================================

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref, date_to_jd(2024, 1, 1)
                )

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                r, v = Propagators.propagate_to_epoch!(orbp, date_to_jd(2024, 1, 1), Tuple)

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                sv = Propagators.propagate_to_epoch!(
                    orbp, date_to_jd(2024, 1, 1), OrbitStateVector
                )

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref, [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)]
                )

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                vr, vv = Propagators.propagate_to_epoch!(
                    orbp, [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)], Tuple
                )

                @test orbp isa typeof(orbp_ref)
                @test vr isa Vector{SVector{3, T}}
                @test vv isa Vector{SVector{3, T}}
                @test vr ≈ vr_ref
                @test vv ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                vsv = Propagators.propagate_to_epoch!(
                    orbp, [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)], OrbitStateVector
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == step! =================================================================

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref, date_to_jd(2024, 1, 1)
                )

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                r, v = Propagators.step!(orbp, 365 * 86400, Tuple)

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                sv = Propagators.step!(orbp, 365 * 86400, OrbitStateVector)

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)
            end
        end
    end

    @testset "Dates Support" verbose = true begin
        jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

        for (T, j2c) in ((Float64, J2C_EGM2008), (Float32, J2C_EGM2008_F32))
            @testset "$T" begin
                orb = KeplerianElements(
                    jd₀,
                    T(8000e3),
                    T(0.015),
                    T(28.5) |> deg2rad,
                    T(100) |> deg2rad,
                    T(200) |> deg2rad,
                    T(45) |> deg2rad,
                )

                orbp_ref = Propagators.init(Val(:J2), orb; j2c = j2c)

                # == propagate =============================================================

                r_ref, v_ref = Propagators.propagate!(orbp_ref, 61)

                r, v, orbp = Propagators.propagate(
                    Tuple, Val(:J2), Dates.Minute(1) + Dates.Second(1), orb; j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                sv, orbp = Propagators.propagate(
                    OrbitStateVector,
                    Val(:J2),
                    Dates.Minute(1) + Dates.Second(1),
                    orb;
                    j2c = j2c,
                )

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                vr_ref, vv_ref = Propagators.propagate!(orbp_ref, [60, 61])

                vr, vv, orbp = Propagators.propagate(
                    Tuple,
                    Val(:J2),
                    [Dates.Minute(1) + Dates.Second(0), Dates.Minute(1) + Dates.Second(1)],
                    orb;
                    j2c = j2c,
                )

                @test orbp isa typeof(orbp_ref)
                @test vr isa Vector{SVector{3, T}}
                @test vv isa Vector{SVector{3, T}}
                @test vr ≈ vr_ref
                @test vv ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                vsv, orbp = Propagators.propagate(
                    OrbitStateVector,
                    Val(:J2),
                    [Dates.Minute(1) + Dates.Second(0), Dates.Minute(1) + Dates.Second(1)],
                    orb;
                    j2c = j2c,
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == propagate! ============================================================

                r_ref, v_ref = Propagators.propagate!(orbp_ref, 61)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                r, v = Propagators.propagate!(
                    orbp, Dates.Minute(1) + Dates.Second(1), Tuple
                )

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                sv = Propagators.propagate!(
                    orbp, Dates.Minute(1) + Dates.Second(1), OrbitStateVector
                )

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                vr_ref, vv_ref = Propagators.propagate!(orbp_ref, [60, 61])

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                vr, vv = Propagators.propagate!(
                    orbp,
                    [Dates.Minute(1) + Dates.Second(0), Dates.Minute(1) + Dates.Second(1)],
                    Tuple,
                )

                @test orbp isa typeof(orbp_ref)
                @test vr isa Vector{SVector{3, T}}
                @test vv isa Vector{SVector{3, T}}
                @test vr ≈ vr_ref
                @test vv ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                vsv = Propagators.propagate!(
                    orbp,
                    [Dates.Minute(1) + Dates.Second(0), Dates.Minute(1) + Dates.Second(1)],
                    OrbitStateVector,
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == propagate_to_epoch ========================================================

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref, date_to_jd(2024, 1, 1)
                )

                r, v, orbp = Propagators.propagate_to_epoch(
                    Tuple, Val(:J2), DateTime("2024-01-01"), orb; j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                sv, orbp = Propagators.propagate_to_epoch(
                    OrbitStateVector, Val(:J2), DateTime("2024-01-01"), orb; j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                vr_ref, vv_ref = Propagators.propagate_to_epoch!(
                    orbp_ref, [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)]
                )

                vr, vv, orbp = Propagators.propagate_to_epoch(
                    Tuple,
                    Val(:J2),
                    [DateTime("2024-01-01"), DateTime("2024-01-02")],
                    orb;
                    j2c = j2c,
                )

                @test orbp isa typeof(orbp_ref)
                @test vr isa Vector{SVector{3, T}}
                @test vv isa Vector{SVector{3, T}}
                @test vr ≈ vr_ref
                @test vv ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                vsv, orbp = Propagators.propagate_to_epoch(
                    OrbitStateVector,
                    Val(:J2),
                    [DateTime("2024-01-01"), DateTime("2024-01-02")],
                    orb;
                    j2c = j2c,
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == propagate_to_epoch! =======================================================

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref, date_to_jd(2024, 1, 1)
                )

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                r, v = Propagators.propagate_to_epoch!(orbp, DateTime("2024-01-01"), Tuple)

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                sv = Propagators.propagate_to_epoch!(
                    orbp, DateTime("2024-01-01"), OrbitStateVector
                )

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref, [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)]
                )

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                vr, vv = Propagators.propagate_to_epoch!(
                    orbp, [DateTime("2024-01-01"), DateTime("2024-01-02")], Tuple
                )

                @test orbp isa typeof(orbp_ref)
                @test vr isa Vector{SVector{3, T}}
                @test vv isa Vector{SVector{3, T}}
                @test vr ≈ vr_ref
                @test vv ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                vsv = Propagators.propagate_to_epoch!(
                    orbp, [DateTime("2024-01-01"), DateTime("2024-01-02")], OrbitStateVector
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == step! =================================================================

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref, date_to_jd(2024, 1, 1)
                )

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                r, v = Propagators.step!(orbp, Dates.Day(365), Tuple)

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                sv = Propagators.step!(orbp, Dates.Day(365), OrbitStateVector)

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)
            end
        end
    end
end

@testset "Fitting Mean Elements Using OrbitStateVector" verbose = true begin
    for (T, j2c) in ((Float64, J2C_EGM2008), (Float32, J2C_EGM2008_F32))
        @testset "$T" begin
            jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

            orb_input = KeplerianElements(
                jd₀,
                T(8000e3),
                T(0.015),
                T(28.5) |> deg2rad,
                T(100) |> deg2rad,
                T(200) |> deg2rad,
                T(45) |> deg2rad,
            )

            orbp_input = Propagators.init(Val(:J2osc), orb_input; j2c = j2c)

            # Create the osculating state vectors that we will use to fit the mean elements
            # and compare with the reference Keplerian elements.
            vsv = Propagators.propagate!(orbp_input, 0:10:12_000, OrbitStateVector)

            # == Fit the Orbit =============================================================

            # Create a new, uninitialized propagator.
            j2d = J2Propagator{Float64, T}()
            j2d.j2c = j2c
            j2oscd = J2OsculatingPropagator{Float64, T}()
            j2oscd.j2d = j2d
            orbp = OrbitPropagatorJ2Osculating(j2oscd)

            orb, _ = redirect_stdout(devnull) do
                Propagators.fit_mean_elements!(
                    orbp, vsv; max_iterations = 10, mean_elements_epoch = jd₀
                )
            end

            @test orb.t ≈ orb_input.t
            @test orb.a ≈ orb_input.a
            @test orb.e ≈ orb_input.e
            @test orb.i ≈ orb_input.i
            @test orb.Ω ≈ orb_input.Ω
            @test orb.ω ≈ orb_input.ω
            @test orb.f ≈ orb_input.f atol = (T === Float32 ? 1e-4 : 1e-5)

            # Fit the mean elements without explicitly initializing the propagator.
            orb, _ = redirect_stdout(devnull) do
                Propagators.fit_mean_elements(
                    Val(:J2osc), vsv; max_iterations = 10, mean_elements_epoch = jd₀
                )
            end

            @test orb.t ≈ orb_input.t
            @test orb.a ≈ orb_input.a
            @test orb.e ≈ orb_input.e
            @test orb.i ≈ orb_input.i
            @test orb.Ω ≈ orb_input.Ω
            @test orb.ω ≈ orb_input.ω
            @test orb.f ≈ orb_input.f atol = (T === Float32 ? 1e-4 : 1e-5)
        end
    end
end

@testset "Default Functions in the API" begin
    orbp = DummyPropagator{Float64, Float64}()
    @test Propagators.name(orbp) == "DummyPropagator{Float64, Float64}"
    @test Propagators.mean_elements(orbp) === nothing
end

@testset "Number of Tasks in the Vectorized Propagation" verbose = true begin
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

    orb = KeplerianElements(
        jd₀,
        8000e3,
        0.015,
        28.5 |> deg2rad,
        100.0 |> deg2rad,
        200.0 |> deg2rad,
        45.0 |> deg2rad,
    )

    vt  = collect(0.0:10:200)
    vjd = jd₀ .+ vt ./ 86400

    # Reference obtained propagating each instant separately.
    orbp = Propagators.init(Val(:J2), orb)
    vr_ref = [Propagators.propagate!(orbp, t)[1] for t in vt]
    vv_ref = [Propagators.propagate!(orbp, t)[2] for t in vt]

    # The result must not depend on the number of tasks. `ntasks` lower than 1 must fall
    # back to a sequential propagation instead of leaving the output uninitialized, and
    # `ntasks` higher than the number of propagated instants must not lead to tasks writing
    # concurrently to the same output elements.
    @testset "propagate! and propagate_to_epoch!" begin
        for ntasks in (-1, 0, 1, 2, 3, 5, 128)
            orbp = Propagators.init(Val(:J2), orb)
            vr, vv = Propagators.propagate!(orbp, vt; ntasks = ntasks)

            @test vr == vr_ref
            @test vv == vv_ref

            # The propagator must be left at the last requested instant.
            @test Propagators.last_instant(orbp) == last(vt)

            orbp = Propagators.init(Val(:J2), orb)
            vr, vv = Propagators.propagate_to_epoch!(orbp, vjd; ntasks = ntasks)

            @test length(vr) == length(vt)
            @test length(vv) == length(vt)
        end
    end

    # The vector must be propagated correctly regardless of its length, which exercises the
    # cases in which there are fewer instants to partition than tasks.
    @testset "Short Time Vectors" begin
        for len in 1:5, ntasks in (0, 1, 8)
            orbp = Propagators.init(Val(:J2), orb)
            vr, vv = Propagators.propagate!(orbp, vt[1:len]; ntasks = ntasks)

            @test vr == vr_ref[1:len]
            @test vv == vv_ref[1:len]
        end
    end

    # `ntasks` is a propagation keyword and must not be forwarded to the initialization
    # function, which does not accept it.
    @testset "Non-mutating Functions" begin
        for ntasks in (0, 1, 4)
            vr, vv, orbp = Propagators.propagate(Val(:J2), vt, orb; ntasks = ntasks)

            @test vr == vr_ref
            @test vv == vv_ref
            @test orbp isa OrbitPropagatorJ2

            vsv, orbp = Propagators.propagate(
                OrbitStateVector, Val(:J2), vt, orb; ntasks = ntasks
            )

            @test length(vsv) == length(vt)
            @test [sv.r for sv in vsv] == vr_ref

            vr, vv, orbp = Propagators.propagate_to_epoch(
                Val(:J2), vjd, orb; ntasks = ntasks
            )

            @test length(vr) == length(vt)

            vsv, orbp = Propagators.propagate_to_epoch(
                OrbitStateVector, Val(:J2), vjd, orb; ntasks = ntasks
            )

            @test length(vsv) == length(vt)
            @test [sv.t for sv in vsv] == vjd
        end
    end
end

@testset "Epoch Type When Fitting Mean Elements" verbose = true begin
    # The fitted mean elements must keep the epoch in `Tepoch`, which is the whole point of
    # having a separate epoch type. Converting it to the element type truncates the Julian
    # Day and makes the return type differ from the documented one.
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

    orb = KeplerianElements(
        jd₀,
        Float32(7130.982e3),
        Float32(0.001111),
        Float32(98.405 |> deg2rad),
        Float32(90.0 |> deg2rad),
        Float32(200.0 |> deg2rad),
        Float32(45.0 |> deg2rad),
    )

    vt  = collect(0.0:60:6000)
    vjd = jd₀ .+ vt ./ 86400

    @testset "$prop" for (prop, kwargs) in (
        (:J2, (; j2c = J2C_EGM2008_F32)),
        (:J2osc, (; j2c = J2C_EGM2008_F32)),
        (:J4, (; j4c = J4C_EGM2008_F32)),
        (:J4osc, (; j4c = J4C_EGM2008_F32)),
        (:TwoBody, (; m0 = TBC_M0_F32)),
    )
        orbp = Propagators.init(Val(prop), orb; kwargs...)
        ret  = [Propagators.propagate!(orbp, t) for t in vt]

        vr_i = [SVector{3, Float32}(r) for r in first.(ret)]
        vv_i = [SVector{3, Float32}(v) for v in last.(ret)]

        orbk, P, stats = Propagators.fit_mean_elements!(
            orbp,
            vjd,
            vr_i,
            vv_i;
            mean_elements_epoch = vjd[begin],
            verbose             = false,
        )

        @test orbk isa KeplerianElements{MeanAnomaly, Float64, Float32}
        @test P isa SMatrix{6, 6, Float32}
        @test stats.converged isa Bool
        @test stats.iterations isa Int
        @test stats.position_rmse isa Float32
        @test stats.velocity_rmse isa Float32
        @test stats.total_rmse isa Float32
        @test orbk.t == vjd[begin]
    end
end

@testset "Initialization With Invalid Orbit Elements" verbose = true begin
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

    orbit(a, e) = KeplerianElements(jd₀, a, e, 10.0 |> deg2rad, 0.0, 0.0, 0.0)

    # The propagators implement theories that are only valid for elliptical orbits. An
    # invalid element must produce a message pointing at it instead of a `DomainError`
    # raised by an internal square root.
    @testset "$prop" for prop in (:J2, :J2osc, :J4, :J4osc, :TwoBody)
        @test_throws ArgumentError Propagators.init(Val(prop), orbit(7000e3, 1.0))
        @test_throws ArgumentError Propagators.init(Val(prop), orbit(7000e3, 1.5))
        @test_throws ArgumentError Propagators.init(Val(prop), orbit(7000e3, -0.1))
        @test_throws ArgumentError Propagators.init(Val(prop), orbit(-7000e3, 0.0))

        # A valid orbit must still be accepted.
        @test Propagators.init(Val(prop), orbit(7000e3, 0.01)) isa OrbitPropagator
    end
end

@testset "Osculating Elements Are Wrapped" verbose = true begin
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

    # The osculating angles must be wrapped to [0, 2π), like the mean elements. Otherwise,
    # the argument of perigee can be negative or larger than 2π, since it is obtained from
    # the difference between the argument of latitude and the true anomaly.
    orb = KeplerianElements(
        jd₀,
        7190.982e3,
        0.001111,
        98.405 |> deg2rad,
        350.0 |> deg2rad,
        359.0 |> deg2rad,
        359.0 |> deg2rad,
    )

    @testset "$prop" for prop in (:J2osc, :J4osc)
        orbp = Propagators.init(Val(prop), orb)

        for t in 0:97:6000
            Propagators.propagate!(orbp, Float64(t))
            orbk = Propagators.mean_elements(orbp)

            @test 0 <= orbk.Ω < 2π
            @test 0 <= orbk.ω < 2π
            @test 0 <= orbk.f < 2π
        end
    end
end

@testset "Osculating Correction Does Not Overflow in Float32" verbose = true begin
    # The short-period correction to the radial rate used to be scaled by `√(p^5)`, which
    # overflows in `Float32` for orbits above roughly 50 900 km. Dividing by the resulting
    # `Inf` silently discarded the whole correction instead of raising an error.
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

    @testset "$prop" for (prop, kwargs_32, kwargs_64) in (
        (:J2osc, (; j2c = J2C_EGM2008_F32), (; j2c = J2C_EGM2008)),
        (:J4osc, (; j4c = J4C_EGM2008_F32), (; j4c = J4C_EGM2008)),
    )
        for a in (6.0e7, 1.0e8, 3.844e8)
            orb_32 = KeplerianElements(
                jd₀, Float32(a), 0.01f0, Float32(10 |> deg2rad), 0.0f0, 0.0f0, 0.0f0
            )

            orb_64 = KeplerianElements(jd₀, a, 0.01, 10 |> deg2rad, 0.0, 0.0, 0.0)

            orbp_32 = Propagators.init(Val(prop), orb_32; kwargs_32...)
            orbp_64 = Propagators.init(Val(prop), orb_64; kwargs_64...)

            r_32, v_32 = Propagators.propagate!(orbp_32, 3600.0f0)
            r_64, v_64 = Propagators.propagate!(orbp_64, 3600.0)

            @test all(isfinite, r_32)
            @test all(isfinite, v_32)

            # The `Float32` result must agree with the `Float64` one to within the `Float32`
            # resolution, which is not the case if the correction is discarded.
            @test maximum(abs.(v_32 .- v_64)) < 1e-3 * maximum(abs.(v_64))
        end
    end
end

@testset "Allocations When Fitting From Non-Static Vectors" verbose = true begin
    # The docstring examples pass the measurements as `Vector{Vector{Float64}}`. Building
    # the measurement vector with `vcat` allocated a new array for every measurement of
    # every iteration, which `test/performance.jl` never caught because it only exercises
    # `Vector{SVector{3, Float64}}`.
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

    orb = KeplerianElements(
        jd₀,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        90.0 |> deg2rad,
        200.0 |> deg2rad,
        45.0 |> deg2rad,
    )

    vt  = collect(0.0:60.0:6000.0)
    vjd = jd₀ .+ vt ./ 86400

    @testset "$prop" for (prop, fit!, build) in (
        (:J2, fit_j2_mean_elements!, () -> begin
            d = J2Propagator{Float64, Float64}()
            d.j2c = J2C_EGM2008
            d
        end),
        (
            :J2osc,
            fit_j2osc_mean_elements!,
            () -> begin
                d = J2OsculatingPropagator{Float64, Float64}()
                d.j2d = J2Propagator{Float64, Float64}()
                d.j2d.j2c = J2C_EGM2008
                d
            end,
        ),
        (:J4, fit_j4_mean_elements!, () -> begin
            d = J4Propagator{Float64, Float64}()
            d.j4c = J4C_EGM2008
            d
        end),
        (
            :J4osc,
            fit_j4osc_mean_elements!,
            () -> begin
                d = J4OsculatingPropagator{Float64, Float64}()
                d.j4d = J4Propagator{Float64, Float64}()
                d.j4d.j4c = J4C_EGM2008
                d
            end,
        ),
        (
            :TwoBody,
            fit_twobody_mean_elements!,
            () -> begin
                d = TwoBodyPropagator{Float64, Float64}()
                d.μ = TBC_M0
                d
            end,
        ),
    )
        orbp = Propagators.init(Val(prop), orb)
        ret  = [Propagators.propagate!(orbp, t) for t in vt]

        vr_i = [collect(r) for r in first.(ret)]
        vv_i = [collect(v) for v in last.(ret)]

        f() = fit!(
            build(),
            vjd,
            vr_i,
            vv_i;
            mean_elements_epoch = vjd[begin],
            jacobian_method     = FiniteDiffJacobian(),
            verbose             = false,
        )

        # Compile before measuring.
        f()

        # This input used to allocate between 34 kB and 57 kB, all of it proportional to
        # the number of measurements. The limit is well above the roughly 1 kB allocated
        # now, but far below the previous figures.
        @test (@allocated f()) < 10_000
    end
end
