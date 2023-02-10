using ParametricModels
using UnPack
using DocStringExtensions
using Statistics
using StatsBase

struct WidthGrowth{T} end # width growth type
struct Trophic{T} end # trophic type
struct Competition{T} end # competition type

"""
    generate_network(SR::Int, SC::Int)

return matrix W[i,j], which is nonzero if consumer i eats resource j
SR is the number of resource, SC the number of consumer species.

The bipartite network is generated as follows. 
First, both resources and consumers are labeled consecutively, 
based on their initial temperature adaptations: 
resource 1 / consumer 1 are the most cold-adapted, 
and resource S / consumer S the most warm-adapted. 
Next, we always put a feeding link between consumer i and resource i. 
Finally, each consumer is randomly linked to four other resource species. 
This yields a feeding network where every consumer is 
connected to five resources altogether.

Note: it seems that the definition of W_ij is that it determines which resource i is eat by consumer j
"""
function generate_network(SR::Int, SC::Int)
    w = zeros(SR+SC, SR+SC) # initialize adjacency matrix
    for i in 1:SR
        # determine which resources each consumer eats: it must eat
        # the one with matching trait value, plus a fixed number of
        # randomly assigned ones (in this case 4 more, for 5 resources per consumer)
        indices = sort([i; sample((1:SC)[.!(i .== (1:SC))], 4, replace = false)])
        w[i+SR, indices] .= 1
    end
    # omega = [] # initialize matrix of consumption efforts
    # rsum = sum(w, dims=1) # omega[i,j] is the proportion of i's consumption rate,
    # for i in 1:(SR+SC) 
    #     push!(omega, rsum) # targeted at consuming j
    #     omega[omega .!= 0] = 1 ./ omega[omega .!= 0] # if not 0, set to proportion
    #     W = w .* omega # only the product of w and omega is used
    # end
    omega = repeat(sum(w, dims=2), outer = (1,SR + SC))
    omega[omega .!= 0] .= 1. ./ omega[omega .!= 0]
    return w .* omega
end


# using DataFrames
# """
#     organize_data(dat, times, pars)

# This code is probably not working
# """
# # This code is probably not working
# function organize_data(dat, times, pars)
#     dat = DataFrame(dat) # convert to data frame 
#     dat = filter(dat, dat.time .∈ times) # only keep specified time points
#     names!(dat, Symbol("time")) # name the first column "time"
#     index = 2 # keep track of which column we are naming with this counter
#     for k in 1:pars.L
#         for i in 1:pars.S
#             # name columns for densities
#             # naming convention: "type_species_patch" - type is either m (trait), or n (density)
#             names!(dat, index, "n_$i_$k") 
#             index += 1
#         end
#     end
#     for k in 1:pars.L
#         for i in 1:pars.S
#             # name columns for trait values
#             # (same naming convention)
#             names!(dat, index, "m_$i_$k")
#             index += 1
#         end
#     end
#     dat = pivot_longer(dat, cols=2:size(dat, 2), names_to="variable", values_to="v") # normalize table by collapsing columns into a key-value column pair
#     dat = separate(dat, :variable, ["type", "species", "patch"], sep="_") # split "variable" into value type (density or trait), species, and patch
#     dat = mutate(dat, :species => parse.(Int, :species), :patch => parse.(Int, :patch)) # convert species & patch from string ("1","2",...) to integer (1,2,...)
#     dat = pivot_wider(dat, names_from="type", values_from="v") # split trait and abundance values into two columns
#     dat = mutate(dat, :tl => ifelse.(:species .> SR, "C", "R")) # trophic level (tl): species with index greater than SR are consumers ("C"), the rest are resources ("R")
#     return dat # return tidy table
# end


# from C++ code
"""
Apply twice continuously differentiable smoothed step function to a number x
Input: x: Distance from pole, measured in units of the pole-to-equator distance
Output: 0 if x < 0; 10*x^3-15*x^4+6*x^5 if 0 <= x <= 1; otherwise 1
"""
function smoothstep(x)
    y = 0.0
    if x < 0.0
        y = 0.0
    elseif x > 1.0
        y = 1.0
    else
        y = x^3 * (10.0 + x * (-15.0 + 6.0 * x))
    end
    return y
end

"""
    $SIGNATURES

Temperature as a function of space (x), time (t), and some climate parameters
"""
Base.@kwdef struct Temperature
    tE::Float64 = 300.0 # time at which climate change stops (assuming it starts at t = 0),
    Cmax::Float64 = 9.66 # projected temperature increase at poles
    Cmin::Float64 = 1.26 # projected temperature increase at equator
    Tmax::Float64 = 25.0 # initial mean temperature at equator
    Tmin::Float64 = -10.0 # initial mean temperature at poles
end

function (temp_fun::Temperature)(x, t)
    @unpack tE, Cmax, Cmin, Tmax, Tmin = temp_fun
    mytemp = (Tmax-Tmin)*x+Tmin+((Cmin-Cmax)*x+Cmax)*smoothstep(t/tE)
    return mytemp
end

"""
    $SIGNATURES
Type II functional response
Input:
- n: Vector of population densities of all species in a given patch
- Th: Vector of handling times (with dummy values for resource species)
- arate: Vector of attack rates (with dummy values for resource species)
- W: Adjacency matrix of trophic network; W(i,j)=1 if i eats j and 0 otherwise
Output:
- A matrix F(i,j), the feeding rate of consumer i on resource j
"""
function funcresp!(F, n, p, ::Trophic{true})
    @unpack arate, Th, W = p
    S = length(n)
    for i=1:S
        Wn = 0.0
        for j=1:S
            Wn += W[i, j]*n[j]
        end
        for j=1:S
            F[i, j] = arate[i]*W[i, j]*n[j]/(1.0 + arate[i]*Th[i]*Wn)
        end
    end
end

function funcresp!(F, n, p, ::Trophic{false})
    return nothing
end

"""
    $SIGNATURES

Returns landscape parameters
"""
struct Landscape
    L::Int64 # number of patches
    x::Vector{Float64} # latitude for each patch
    mig::Matrix{Int64} # dispersal matrix
end

function Landscape(L)
    # dispersal matrix
    if L > 1
        x = collect(0:L-1) ./ (Float64(L) - 1.0)
    else
        x = [0.5]
    end
    mig = zeros(L, L) # initialize dispersal matrix
    for k in 2:L
        mig[k-1, k] = 1 # each species can only migrate to the two
    end
    mig = mig + mig' # nearest-neighbor patches
    Landscape(L, x, mig)
end