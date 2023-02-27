using ParametricModels
using UnPack
using DocStringExtensions
using Statistics
using ComponentArrays
using Distributions
using OrdinaryDiffEq
using UnPack
using Test

## Defining phenotypic space and spatial graph
M = 7
g = star_graph(M)

rS = 1f0
dS = 0.02f0
phen_space = collect(range(-rS,rS,step=dS)) #grid

##
# Defining birth and competition functions
soptim = 0.5f0 * [-1,1,-1,1,-1,1,-1]
birth_fn(x, s, p) = max(0f0, 1f0 - (soptim[x] - s)^2)
competition_fn(u_xs2, x, s1, s2, p) = u_xs2 ./ p.K

## Defining parameters
σ_mu = 5f-2;
mu = 0.1f0
m = 0.1
K = 1.

p = ComponentArray(σ_mu = σ_mu,
                    mu = mu,
                    m = m,
                    K = K)


## rest of the simulation
tend = 1000f0
tspan = (0f0,tend)
tsteps = (tspan[1]):1.:(tspan[end])
u0 = vcat([K .* pdf.(Normal(so,σ_mu),phen_space') for so in soptim]...)

@testset "EcoEvoGraph" begin 
    mp = ModelParams(;p,
                    tspan,
                    u0,
                    alg=Tsit5(),
                    saveat = tsteps)

    model = EcoEvoGraph(mp, g, phen_space, birth_fn, competition_fn)
    sol = simulate(model)
    @test SciMLBase.successful_retcode(sol.retcode)
end