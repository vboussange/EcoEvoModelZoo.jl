using EcoEvoModelZoo
using ParametricModels
using Random; Random.seed!(1234)
using LinearAlgebra
using UnPack
using Test
using OrdinaryDiffEq
using Statistics
false ? (using PyPlot) : nothing
# Set a random seed for reproduceable behaviour

@testset "McCann model" begin
    alg = BS3()
    abstol = 1e-6
    reltol = 1e-6
    tspan = (0., 800)
    tsteps = 550:4:800
    p_true = (x_c = [0.4], 
            x_p = [0.08], 
            y_c = [2.01], 
            y_p = [5.00], 
            R_0 = [0.16129], 
            C_0 = [0.5])

    u0_true = [0.5,0.8,0.5]


    model = EcosystemModelMcCann(ModelParams(;p = p_true,
                                            tspan,
                                            u0 = u0_true,
                                            alg,
                                            reltol,
                                            abstol,
                                            saveat = tsteps,
                                            verbose = false, # suppresses warnings for maxiters
                                            maxiters = 50_000,
                                            ))
    sol = simulate(model, u0 = u0_true)
    @test sol.retcode == :Success
end

@testset "McCann model with omnivory" begin
    alg = BS3()
    abstol = 1e-6
    reltol = 1e-6
    tspan = (0., 800)
    tsteps = 550:4:800
    p_true = (x_c = [0.4], 
            x_p = [0.08], 
            y_c = [2.01], 
            y_pr = [2.00], 
            y_pc = [5.0], 
            R_0 = [0.16129], 
            R_02 = [ 0.5], 
            C_0 = [0.5],
            ω =[ 0.4])

    u0_true = [0.5,0.8,0.5]


    model = EcosystemModelOmnivory(ModelParams(;p = p_true,
                                            tspan,
                                            u0 = u0_true,
                                            alg,
                                            reltol,
                                            abstol,
                                            saveat = tsteps,
                                            verbose = false, # suppresses warnings for maxiters
                                            maxiters = 50_000,
                                            ))
    sol = simulate(model, u0 = u0_true)
    @test sol.retcode == :Success
end