using Test

using Dates
using StaticArrays
using SatelliteToolboxPropagators

@testset "J2 Orbit Propagator" begin
    include("./j2.jl")
end

@testset "J4 Orbit Propagator" begin
    include("./j4.jl")
end

@testset "SGP4 Orbit Propagator" begin
    include("./sgp4.jl")
end

@testset "Orbit Propagator API" begin
    include("./api.jl")
end
