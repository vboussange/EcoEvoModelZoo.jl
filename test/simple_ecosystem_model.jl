using EcoEvoModelZoo
using ParametricModels
# using Random; Random.seed!(1234)
using LinearAlgebra
using UnPack
using Test
using OrdinaryDiffEq
using Statistics
using Graphs
# Set a random seed for reproduceable behaviour


@testset "5 species ecosystem model" begin
        alg = BS3()
        abstol = 1e-6
        reltol = 1e-6
        tspan = (0., 100)
        tsteps = (tspan[1]):1.:(tspan[end])
        N = 5 # number of compartment

        # constructing simple chain foodweb
        foodweb = DiGraph(N)
        [add_edge!(foodweb,i=>i-1) for i in 2:N]

        p_true = (r = vcat(1., - 0.3 * ones(N-1)), 
                K = ones(N), 
                A = diagm(vcat(1, zeros(N-1))), 
                ϵ = vcat(0., ones(N-1)), 
                q = vcat(0., ones(N-1)),
                H = vcat(0., ones(N-1)), 
                W = adjacency_matrix(foodweb))

        u0_true = rand(N)


        model = SimpleEcosystemModel(;mp = ModelParams(;p = p_true,
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
        @test SciMLBase.successful_retcode(sol.retcode)
        # plotting
        # using PythonCall; plt = pyimport("matplotlib.pyplot")

        # fig, ax = plt.subplots(1)
        # for i in 1:N
        #         ax.plot(Array(sol)[i,:], label = "Species $i")
        # end
        # fig.legend()
        # display(fig)
end
# Testing with hastings
@testset "3 species ecosystem model" begin
        alg = BS3()
        abstol = 1e-6
        reltol = 1e-6
        tspan = (0., 800)
        tsteps = 550:4:800
        N = 3 # number of compartment

        # constructing simple chain foodweb
        foodweb = DiGraph(N)
        [add_edge!(foodweb,i=>i-1) for i in 2:N]

        p_mccann = (x_c = [0.4], 
                x_p = [0.08], 
                y_c = [2.01], 
                y_p = [5.00], 
                R_0 = [0.16129], 
                C_0 = [0.5])

        @unpack x_c, x_p, y_c, y_p, R_0, C_0 = p_mccann
        p = (r = vcat(1., -x_c, -x_p) , 
                K = ones(N), 
                A = diagm(vcat(1, zeros(N-1))), 
                ϵ = vcat(0., ones(N-1)), 
                q = vcat(0., x_c .* y_c ./ R_0, x_p .* y_p ./ C_0),
                H = vcat(0., 1. ./ (x_c .* y_c), 1. ./ (x_p .* y_p)), 
                W = adjacency_matrix(foodweb))

        u0_true = [0.5,0.8,0.5]


        model = SimpleEcosystemModel(;mp = ModelParams(;p,
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
        @test SciMLBase.successful_retcode(sol.retcode)

        # plotting
        # using PythonCall; plt = pyimport("matplotlib.pyplot")
        # fig, ax = plt.subplots(1)
        # for i in 1:N
        #         ax.plot(Array(sol)[i,:], label = "Species $i")
        # end
        # fig.legend()
        # display(fig)

        # comparing
        model_mccann = EcosystemModelMcCann(ModelParams(;p = p_mccann,
                                                tspan,
                                                u0 = u0_true,
                                                alg,
                                                reltol,
                                                abstol,
                                                saveat = tsteps,
                                                verbose = false, # suppresses warnings for maxiters
                                                maxiters = 50_000,
                                                ))
        sol_mccann = simulate(model_mccann, u0 = u0_true)

        @test all(Array(sol_mccann) .≈ Array(sol))
end
