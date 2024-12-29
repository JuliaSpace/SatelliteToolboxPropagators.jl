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

@testset "Dates Support" verbose = true begin
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

    for (T, j2c) in ((Float64, j2c_egm2008), (Float32, j2c_egm2008_f32))
        @testset "$T" begin
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

            # == propagate =================================================================

            r, v, orbp = Propagators.propagate(
                Val(:J2),
                Dates.Minute(1) + Dates.Second(1),
                orb;
                j2c = j2c
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

            vjd  = collect(jd₀:0.1:jd₁)
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

@testset "Sink Options"  verbose = true begin
    @testset "Julian Day" verbose = true begin
        jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)

        for (T, j2c) in ((Float64, j2c_egm2008), (Float32, j2c_egm2008_f32))
            @testset "$T" begin
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

                # == propagate =============================================================

                r_ref, v_ref = Propagators.propagate!(orbp_ref, 61)

                r, v, orbp = Propagators.propagate(
                    Tuple,
                    Val(:J2),
                    61,
                    orb;
                    j2c = j2c
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
                    61,
                    orb;
                    j2c = j2c
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
                    [60, 61],
                    orb;
                    j2c = j2c
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
                    [60, 61],
                    orb;
                    j2c = j2c
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
                    orbp_ref,
                    date_to_jd(2024, 1, 1)
                )

                r, v, orbp = Propagators.propagate_to_epoch(
                    Tuple,
                    Val(:J2),
                    date_to_jd(2024, 1, 1),
                    orb;
                    j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                sv, orbp = Propagators.propagate_to_epoch(
                    OrbitStateVector,
                    Val(:J2),
                    date_to_jd(2024, 1, 1),
                    orb;
                    j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                vr_ref, vv_ref = Propagators.propagate_to_epoch!(
                    orbp_ref,
                    [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)]
                )

                vr, vv, orbp = Propagators.propagate_to_epoch(
                    Tuple,
                    Val(:J2),
                    [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)],
                    orb;
                    j2c = j2c
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
                    j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == propagate_to_epoch! ===================================================

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref,
                    date_to_jd(2024, 1, 1)
                )

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                r, v = Propagators.propagate_to_epoch!(
                    orbp,
                    date_to_jd(2024, 1, 1),
                    Tuple
                )

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                sv = Propagators.propagate_to_epoch!(
                    orbp,
                    date_to_jd(2024, 1, 1),
                    OrbitStateVector
                )

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref,
                    [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)]
                )

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                vr, vv = Propagators.propagate_to_epoch!(
                    orbp,
                    [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)],
                    Tuple
                )

                @test orbp isa typeof(orbp_ref)
                @test vr isa Vector{SVector{3, T}}
                @test vv isa Vector{SVector{3, T}}
                @test vr ≈ vr_ref
                @test vv ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                vsv = Propagators.propagate_to_epoch!(
                    orbp,
                    [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)],
                    OrbitStateVector
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == step! =================================================================

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref,
                    date_to_jd(2024, 1, 1)
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

        for (T, j2c) in ((Float64, j2c_egm2008), (Float32, j2c_egm2008_f32))
            @testset "$T" begin
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

                # == propagate =============================================================

                r_ref, v_ref = Propagators.propagate!(orbp_ref, 61)

                r, v, orbp = Propagators.propagate(
                    Tuple,
                    Val(:J2),
                    Dates.Minute(1) + Dates.Second(1),
                    orb;
                    j2c = j2c
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
                    j2c = j2c
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
                    j2c = j2c
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
                    j2c = j2c
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
                    orbp,
                    Dates.Minute(1) + Dates.Second(1),
                    Tuple
                )

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                sv = Propagators.propagate!(
                    orbp,
                    Dates.Minute(1) + Dates.Second(1),
                    OrbitStateVector
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
                    Tuple
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
                    OrbitStateVector
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == propagate_to_epoch ========================================================

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref,
                    date_to_jd(2024, 1, 1)
                )

                r, v, orbp = Propagators.propagate_to_epoch(
                    Tuple,
                    Val(:J2),
                    DateTime("2024-01-01"),
                    orb;
                    j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                sv, orbp = Propagators.propagate_to_epoch(
                    OrbitStateVector,
                    Val(:J2),
                    DateTime("2024-01-01"),
                    orb;
                    j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                vr_ref, vv_ref = Propagators.propagate_to_epoch!(
                    orbp_ref,
                    [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)]
                )

                vr, vv, orbp = Propagators.propagate_to_epoch(
                    Tuple,
                    Val(:J2),
                    [DateTime("2024-01-01"), DateTime("2024-01-02")],
                    orb;
                    j2c = j2c
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
                    j2c = j2c
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == propagate_to_epoch! =======================================================

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref,
                    date_to_jd(2024, 1, 1)
                )

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                r, v = Propagators.propagate_to_epoch!(
                    orbp,
                    DateTime("2024-01-01"),
                    Tuple
                )

                @test orbp isa typeof(orbp_ref)
                @test r isa SVector{3, T}
                @test v isa SVector{3, T}
                @test r ≈ r_ref
                @test v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                sv = Propagators.propagate_to_epoch!(
                    orbp,
                    DateTime("2024-01-01"),
                    OrbitStateVector
                )

                @test orbp isa typeof(orbp_ref)
                @test sv isa OrbitStateVector{Float64, T}
                @test sv.r ≈ r_ref
                @test sv.v ≈ v_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # -- Array -----------------------------------------------------------------

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref,
                    [date_to_jd(2024, 1, 1), date_to_jd(2024, 1, 2)]
                )

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                vr, vv = Propagators.propagate_to_epoch!(
                    orbp,
                    [DateTime("2024-01-01"), DateTime("2024-01-02")],
                    Tuple
                )

                @test orbp isa typeof(orbp_ref)
                @test vr isa Vector{SVector{3, T}}
                @test vv isa Vector{SVector{3, T}}
                @test vr ≈ vr_ref
                @test vv ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                orbp = Propagators.init(Val(:J2), orb; j2c = j2c)

                vsv = Propagators.propagate_to_epoch!(
                    orbp,
                    [DateTime("2024-01-01"), DateTime("2024-01-02")],
                    OrbitStateVector
                )

                @test orbp isa typeof(orbp_ref)
                @test vsv isa Vector{OrbitStateVector{Float64, T}}
                @test map(x -> x.r, vsv) ≈ vr_ref
                @test map(x -> x.v, vsv) ≈ vv_ref
                @test Propagators.last_instant(orbp) == Propagators.last_instant(orbp_ref)

                # == step! =================================================================

                r_ref, v_ref = Propagators.propagate_to_epoch!(
                    orbp_ref,
                    date_to_jd(2024, 1, 1)
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
