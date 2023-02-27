using DiffEqOperators
using Graphs

# this function is called to evaluate the convolution of local u and alpha evaluated at phen_space[i]
function int_competition_fn(model, u, x, p) where T <: AbstractArray
    @unpack phen_space = model
    dS = phen_space[2]-phen_space[1]
    C = 0.5f0 * (competition_fn(u[1], x, phen_space[1], p) + competition_fn(u[end], x, phen_space[end], p))
    C += sum(competition_fn.(u[2:end-1], x, phen_space[2:end-1], Ref(p)))
    return C*dS
end

function rw_laplacian_matrix(g::AbstractGraph)
    L = laplacian_matrix(g,dir=:out) ./ outdegree(g)
    replace!(L,NaN=>0)
end

"""
    $SIGNATURES

This model is inspired from [Boussange & Pellissier. 2022](https://www.nature.com/articles/s42003-022-03595-3)

# Example
```julia
# Testing model
using ParametricModels
using UnPack
using DocStringExtensions
using Statistics
using ComponentArrays
using Distributions
using OrdinaryDiffEq
using UnPack
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
competition_fn(u_xs, x, s, p) = u_xs ./ p.K

## Defining parameters
σ_mu = 5f-2;
mu = 0.1f0
m = 0.1
K = 1.

p = ComponentArray(σ_mu = [σ_mu],
                    mu = [mu],
                    m = [m],
                    K = [K])


## rest of the simulation
tend = 1000f0
tspan = (0f0,tend)
tsteps = (tspan[1]):1.:(tspan[end])
u0 = vcat([K .* pdf.(Normal(so,σ_mu),phen_space') for so in soptim]...)

mp = ModelParams(;p,
                tspan,
                u0,
                alg=Tsit5(),
                saveat = tsteps)

model = EcoEvoGraph(mp, g, phen_space, birth_fn, competition_fn)

simulate(model)
```
"""
Base.@kwdef struct EcoEvoGraph{MP,G,PS,DS,DX,BFN,CFN} <: AbstractModel
    mp::MP
    g::G # graph
    phen_space::PS
    Δ_s::DS
    Δ_x::DX
    birth_fn::BFN # number of resources
    competition_fn::CFN # supply rates
end

function (model::EcoEvoGraph)(du,u,p,t)
    @unpack phen_space, Δ_s, Δ_x, g, birth_fn, competition_fn = model
    @unpack m, σ_mu, mu = p

    u[u[:,:] .< eps()] .= 0
    B_mat = birth_fn.(1:nv(g), phen_space', Ref(p))

    for i in 1:nv(g) # iterating over graph nodes
        u_i = @view u[i,:]
        C_i = int_competition_fn(model, u_i, i, p)
        # here there is no b (1-m) because the -m is included in the laplacian
        du[i,:] .= u_i .* (B_mat[i,:] .- C_i )
    end

    lux = m .* Δ_x' * (u .* B_mat)
    lus = σ_mu.^2/2 .* mu .* (u .* B_mat) * Δ_s 
    @. du = du - lux + lus
    return du
end

function EcoEvoGraph(mp::ModelParams, g, phen_space, birth_fn, competition_fn)
    dS = phen_space[2]-phen_space[1]
    N = length(phen_space)
    Δ_s = Array(CenteredDifference{1}(2, 2, dS, N))[:,2:end-1]
    Δ_x = rw_laplacian_matrix(g)
    EcoEvoGraph(mp, g, phen_space, Δ_s, Δ_x, birth_fn, competition_fn)
end