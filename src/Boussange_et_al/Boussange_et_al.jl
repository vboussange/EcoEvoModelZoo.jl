using ParametricModels
using UnPack
using DocStringExtensions
using Statistics
using Graphs

# this function is called to evaluate the convolution of local u and alpha evaluated at S[i]
function int_u(u::T, p) where T <: AbstractArray
    @unpack N,S,dS = p
    C = 0.5f0 * (u[1] + u[N])
    C += sum(u[2:N-1])
    return C*dS
end

"""
    $SIGNATURES

This model is inspired from [Boussange & Pellissier. 2022](https://www.nature.com/articles/s42003-022-03595-3)

# Example
```julia

## Defining parameters
mu= 0.1f0
K1 = 1f0 #scaled to 1 for decreasing numerical errors
M = 7;
dS = 0.02f0;
rS = 1f0;
σ_mu = 5f-2;
c = 1f0 / K1;
tend = 1000f0
tspan = (0f0,tend)

## rest of the simulation
g = star_graph(M)
S = collect(range(-rS,rS,step=dS)) #grid
N = length(S)
Δ_s = Array(σ_mu^2/2 *mu * CenteredDifference{1}(2, 2, dS, N))[:,2:end-1]
X = 1:M
soptim = 0.5f0 * Float32[-1,1,-1,1,-1,1,-1]
u0 = vcat([K1 .* pdf.(Normal(so,σ_mu),S') for so in soptim]...)

B(x,s) = max(0f0,1f0 - (soptim[x] - s)^2)
```
"""
Base.@kwdef struct EcoEvoGraph{MP,G,BFN,CFN} <: AbstractModel
    mp::MP
    g::G # graph
    birth_fn::BFN # number of resources
    competition_fn::CFN # supply rates
end

function (m::EcoEvoGraph)(du,u,p,t)
    @unpack N, M, S, dS, X, Δ_s, Δ_x = m
    @unpack m = p

    u[u[:,:] .< eps()] .= 0
    B_i = B.(1:length(g), S')
    Δ_x = m * rw_laplacian_matrix(g)

    for i in 1:M
        u_i = @view u[i,:]
        C_i = int_u((u_i), p)
        # here there is no b (1-m) because the -m is included in the laplacian
        du[i,:] .= u_i .* (B_i[i,:] .- C_i / K1  )
    end
    # @show size(u)
    # @show size(Δ_x')
    # @show size(B_i)
    # @show size(u .* B_i)
    mul!(lux, Δ_x', (u .* B_i))
    mul!(lus, u .* B_i, Δ_s )
    @. du = du - lux + lus
    return du
end

## Parameters used
mu= 0.1f0
K1 = 1f0 #scaled to 1 for decreasing numerical errors
M = 7;
dS = 0.02f0;
rS = 1f0;
σ_mu = 5f-2;
c = 1f0 / K1;
tend = 1000f0
tspan = (0f0,tend)

## rest of the simulation
g = star_graph(M)
S = collect(range(-rS,rS,step=dS)) #grid
N = length(S)
Δ_s = Array(σ_mu^2/2 *mu * CenteredDifference{1}(2, 2, dS, N))[:,2:end-1]
X = 1:M
soptim = 0.5f0 * Float32[-1,1,-1,1,-1,1,-1]
u0 = vcat([K1 .* pdf.(Normal(so,σ_mu),S') for so in soptim]...)

B(x,s) = max(0f0,1f0 - (soptim[x] - s)^2)

lux = similar(u0);
lus = similar(u0);

# non linear term, i.e. the logistic summand of the IDE


p_default = Dict{String,Any}()
@pack! p_default = M,S,dS,N,Δ_s,σ_mu,X

    #test