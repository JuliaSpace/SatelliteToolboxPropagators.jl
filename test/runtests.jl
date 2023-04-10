using Test

using Dates
using StaticArrays
using SatelliteToolboxPropagators

@testset "J2 Orbit Propagator" begin
    include("./j2.jl")
end

@testset "Orbit Propagator API" begin
    include("./api.jl")
end
