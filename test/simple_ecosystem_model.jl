using EcoEvoModelZoo
using ParametricModels
using Random; Random.seed!(1234)
using LinearAlgebra
using UnPack
using Test
using OrdinaryDiffEq
using Statistics
using Graphs
false ? (using PyPlot) : nothing
# Set a random seed for reproduceable behaviour

# TODO: this is work in progress and should be tested

alg = BS3()
abstol = 1e-6
reltol = 1e-6
tspan = (0., 800)
tsteps = (tspan[1]):1.:(tspan[end])
N = 5 # number of compartment

foodweb = Graph(N)

p_true = (r = vcat(1, zeros(N-1)), 
        K = vcat(1, zeros(N-1)), 
        A = diagm(vcat(1, zeros(N-1))), 
        ϵ = vcat(0., ones(N-1)), 
        q = vcat(0., ones(N-1)),
        H = vcat(0., ones(N-1)), 
        W = adjacency_matrix(foodweb))

u0_true = rand(N)


model = SimpleEcosystemModel(ModelParams(;p = p_true,
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
