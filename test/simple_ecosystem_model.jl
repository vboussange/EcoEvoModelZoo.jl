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


@testset "3 species ecosystem model with omnivory, " begin
        # https://royalsocietypublishing.org/doi/10.1098/rspb.1997.0172
        using SimpleWeightedGraphs
        alg = BS3()
        abstol = 1e-6
        reltol = 1e-6
        tspan = (0., 800)
        tsteps = (tspan[1]):1.:(tspan[end])
        N = 3 # number of compartment

        # constructing simple chain foodweb
        # Population vector is organised in such a way
        # N = [R1, C1, P, R2, C2]


        p_mccann = (x_c = 0.4, 
                        x_p = 0.08, 
                        y_c = 2.01, 
                        y_pr = 2.00, 
                        y_pc = 5.0, 
                        R_0 = 0.16129, 
                        R_02 = 0.5, 
                        C_0 = 0.5,
                        ω = 0.2)

        @unpack x_c, x_p, y_c, y_pc, y_pr, R_0, R_02, C_0, ω = p_mccann

        foodweb = SimpleWeightedDiGraph(N)
        add_edge!(foodweb, 2, 1, 1.) # C to R
        add_edge!(foodweb, 3, 2, 1-ω) # P to C
        add_edge!(foodweb, 3, 1, ω) # P to R

        H = sparse(zeros(N,N))
        H[2,1] = 1 / (x_c * y_c); H[3,1] = 1 / (x_p * y_pr); H[3,2] = 1 / (x_p * y_pc)
        q = sparse(zeros(N,N))
        q[2,1] = x_c * y_c / R_0; q[3,1] = x_p * y_pr / R_02; q[3,2] = x_p * y_pc / C_0
        
        p = (r = vcat(1., -x_c, -x_p) , 
                K = ones(N), 
                A = diagm(vcat(1, zeros(N-1))), 
                ϵ = ones(N), 
                q = q,
                H = H, 
                W = adjacency_matrix(foodweb))

        u0_true = [0.5,0.8,0.5]

        # testing feeding
        F = EcoEvoModelZoo.default_feeding(u0_true, p, 0.)
        R, C, P = u0_true
        F_32 = (1. - ω[]) * x_p[] * y_pc[] / (ω[] * R + (1. - ω[]) * C + C_0[]) 
        F_31 = ω[] * x_p[] * y_pr[] / (ω[] * R + (1. - ω[]) * C + R_02[])
        @test F_32 ≈ F[3,2]
        @test F_31 ≈ F[3,1]

        # simulating model
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

        model_mccann = EcosystemModelOmnivory(ModelParams(;p = p_mccann,
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



@testset "5 species ecosystem mode by post et al" begin
        # https://www.jstor.org/stable/177129
        using SimpleWeightedGraphs
        alg = BS3()
        abstol = 1e-6
        reltol = 1e-6
        tspan = (0., 600)
        tsteps = (tspan[1]):1.:(tspan[end])
        N = 5 # number of compartment

        p_mccann = (x_c = 0.15, 
                x_p = 0.08, 
                y_c = 2.3, 
                y_p = 1.7, 
                R_0 = 0.25, 
                C_0 = 0.5,
                ω = 0.2)

        @unpack x_c, x_p, y_c, y_p, R_0, C_0, ω = p_mccann

        # constructing simple chain foodweb
        # Population vector is organised in such a way
        # N = [R1, C1, P, R2, C2]
        foodweb = SimpleWeightedDiGraph(N)
        add_edge!(foodweb, 2, 1, 1.) # C1 to R1
        add_edge!(foodweb, 5, 4, 1.) # C2 to R2
        add_edge!(foodweb, 3, 2, ω) # P to C1
        add_edge!(foodweb, 3, 5, 1-ω) # P to C2

        H = sparse(zeros(N,N))
        H[2,1] = 1 / (x_c * y_c); H[5,4] = 1 / (x_c * y_c); 
        H[3,2] = 1 / (x_p * y_p); H[3,5] = 1 / (x_p * y_p)

        q = sparse(zeros(N,N))
        q[2,1] = x_c * y_c / R_0; q[5,4] = x_c * y_c / R_0;
        q[3,2] = x_p * y_p / C_0; q[3,5] = x_p * y_p / C_0

        p = (r = vcat(1., -x_c, -x_p, 1., -x_c,), 
                K = ones(N), 
                A = diagm(vcat(1,0,0,1,0)), 
                ϵ = ones(N), 
                q = q,
                H = H, 
                W = adjacency_matrix(foodweb))

        u0_true = rand(N)


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
        using PythonCall; plt = pyimport("matplotlib.pyplot")
        fig, ax = plt.subplots(1)
        for i in 1:N
                ax.plot(Array(sol)[i,:], label = "Species $i")
        end
        # ax.set_yscale("log")
        fig.legend()
        display(fig)
end

# To implement: https://www.google.com/search?client=safari&rls=en&q=Coupled+predator%E2%80%93prey+oscillations+in+a+chaotic+food+web&ie=UTF-8&oe=UTF-8