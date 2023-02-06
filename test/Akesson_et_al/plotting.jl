#=
This script tests the plotting functions in `plotting_Akesson.jl`

(Test by hand)
=#
cd(@__DIR__)
using Test
using UnPack
using RCall
using Random
using Distributions
using LinearAlgebra
using ForwardDiff
using EcoEvoModelZoo
using PythonCall; plt = pyimport("matplotlib")

R"source('./ecoevo_norun.R')"
@rget SR SC S L rho kappa a eta eps W venv vmat s nmin aw bw Tmax Tmin Th arate Cmax Cmin tE d mig model ninit muinit
ic = [ninit; muinit] # merge initial conditions into a vector
S = S |> Int
SR = SR |> Int
SC = SC |> Int
L = L |> Int

# coerce parameters into a dictionary
pars = Dict{Symbol, Any}()
@pack! pars = SR, SC, S, L, rho, kappa, a, eta, eps, W, venv, vmat, s, nmin, aw, bw, Tmax, Tmin, Th, arate, Cmax, Cmin, tE, d, mig, model

pars = NamedTuple([pair for pair in pars])

# integrate prob with Julia
tspan = (tstart,0)
saveat = tspan[1]:200:tspan[2]
prob = ODEProblem(eqs!, ic, tspan, pars)

println("Integrating with Julia")
@time before_cc = solve(prob, 
                        alg=Tsit5(),
                        saveat = saveat)

u = Array(before_cc)

#=
Ploting N
=#
fig, ax = subplots(1)
ts = 1:length(saveat)
l = 25 # patch index
plotting_N_through_time(ax, u, ts, l, pars)
display(fig)

#=
Ploting mu
=#
fig, ax = subplots(1)
l = 10 # patch index
plotting_mu_through_time(ax, u, ts, l, pars)
display(fig)

#=
Ploting local phenotypic distribution
=#
fig, ax = subplots(1)
l = 10 # patch index
species = 1
t = length(saveat)
plotting_distribution(ax, u, t, species, l, pars)
display(fig)


#=
Plotting global phenotypic distribution
=#
fig, ax = subplots(1)
t = 1
# t = length(saveat)
# plotting resources
for species in sample(1:SR, 5, replace=false)
    plotting_distribution_global(ax, u, t, species, pars)
    display(fig)
end