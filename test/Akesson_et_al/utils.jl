cd(@__DIR__)
using Test
using UnPack
using RCall
using Random
using Distributions
using LinearAlgebra
using ForwardDiff
using OrdinaryDiffEq
using EcoEvoModelZoo

test_rcall = false
test_plot = false

test_rcall ? (using RCall) : nothing
test_plot ? (using PyPlot) : nothing

@testset "generate_network" begin
    SR = 10
    SC = 10
    @test size(generate_network(SR, SC)) == (SR + SC, SR + SC)
end

@testset "organize_data" begin
    println("no test implemented")
end

@testset "smoothstep" begin
    @test smoothstep(-1.0) == 0
    @test smoothstep(2.0) == 1
    @test 0 < smoothstep(0.5) < 1
end

@testset "Landscape" begin
    land = Landscape(10)
    @test land.mig isa Matrix
end


@testset "Temperature" begin
    temp = Temperature()
    land = Landscape(1)
    t = 50
    @test temp.(land.x, t) isa Vector
    land = Landscape(10)
    @test temp.(land.x, t) isa Vector

end

@testset "funcresp" begin
    SR = 10
    SC = 10
    n = rand(SC + SR) # Vector of population densities of all species in a given patch
    W = generate_network(SR, SC)
    Th = ones(SC + SR) # Vector of handling times (with dummy values for resource species)
    arate = ones(SC + SR)
    F = zeros(SR + SC, SR + SC)
    p = Dict{String, Any}()
    @pack! p = W, Th, arate
    trophic = Trophic{true}()
    funcresp!(F, n, p, trophic)
    @test F isa Matrix
    @test !any(isnan.(F))
    @test !any(isinf.(F))
end

if test_rcall
    @testset "init_params" begin

        # default parameters
        pars, ic = init_params()

        # default R parameters
        R"source('/Users/victorboussange/ETHZ/projects/piecewise-inference/code/model/spatial_ecoevo-1.0.0/ecoevo_norun.R')"
        @rget SR SC S L rho kappa a eta eps W venv vmat s nmin aw bw Tmax Tmin Th arate Cmax Cmin tE d mig model ninit muinit
        ic = [ninit; muinit] # merge initial conditions into a vector
        S = S |> Int
        SR = SR |> Int
        SC = SC |> Int
        L = L |> Int
        x = collect(0:L-1) ./ (Float64(L) - 1.0)

        # coerce parameters into a dictionary
        parsR = Dict{String,Any}()
        @pack! parsR = SR, SC, S, L, rho, kappa, a, eta, eps, W, venv, vmat, s, nmin, aw, bw, Tmax, Tmin, Th, arate, Cmax, Cmin, tE, d, mig, model, x

        @test pars[:S] == parsR["S"]
        @test pars[:SR] == parsR["SR"]
        @test pars[:SC] == parsR["SC"]
        @test pars[:L] == parsR["L"]
        @test size(pars[:rho]) == size(parsR["rho"])
        # uncomplete
    end

    @testset "Compare `eqs` evaluation against R script" begin
        R"source('/Users/victorboussange/ETHZ/projects/piecewise-inference/code/model/spatial_ecoevo-1.0.0/ecoevo_norun.R')"
        @rget SR SC S L rho kappa a eta eps W venv vmat s nmin aw bw Tmax Tmin Th arate Cmax Cmin tE d mig model ninit muinit
        ic = [ninit; muinit] # merge initial conditions into a vector
        S = S |> Int
        SR = SR |> Int
        SC = SC |> Int
        L = L |> Int
        x = collect(0:L-1) ./ (Float64(L) - 1.0)

        # coerce parameters into a dictionary
        pars = Dict{String,Any}()
        @pack! pars = SR, SC, S, L, rho, kappa, a, eta, eps, W, venv, vmat, s, nmin, aw, bw, Tmax, Tmin, Th, arate, Cmax, Cmin, tE, d, mig, model, x

        dudt = similar(ic)

        eqs!(dudt, ic, pars, 10.0)

        # for comparing, one needs to reorder
        dudt_rlike = [dudt[1:SC+SR, :][:]; dudt[SC+SR+1:end, :][:]]

        # test with R
        dudtR = R"eqs(1., ic, pars)"[1]

        @test all(isapprox.(dudt_rlike, dudtR, atol=1e-3))
    end
end