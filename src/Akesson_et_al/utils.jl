using ParametricModels
using UnPack
using DocStringExtensions
using Statistics
using StatsBase

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


using DataFrames
"""
    organize_data(dat, times, pars)

This code is probably not working
"""
# This code is probably not working
function organize_data(dat, times, pars)
    dat = DataFrame(dat) # convert to data frame 
    dat = filter(dat, dat.time .∈ times) # only keep specified time points
    names!(dat, Symbol("time")) # name the first column "time"
    index = 2 # keep track of which column we are naming with this counter
    for k in 1:pars.L
        for i in 1:pars.S
            # name columns for densities
            # naming convention: "type_species_patch" - type is either m (trait), or n (density)
            names!(dat, index, "n_$i_$k") 
            index += 1
        end
    end
    for k in 1:pars.L
        for i in 1:pars.S
            # name columns for trait values
            # (same naming convention)
            names!(dat, index, "m_$i_$k")
            index += 1
        end
    end
    dat = pivot_longer(dat, cols=2:size(dat, 2), names_to="variable", values_to="v") # normalize table by collapsing columns into a key-value column pair
    dat = separate(dat, :variable, ["type", "species", "patch"], sep="_") # split "variable" into value type (density or trait), species, and patch
    dat = mutate(dat, :species => parse.(Int, :species), :patch => parse.(Int, :patch)) # convert species & patch from string ("1","2",...) to integer (1,2,...)
    dat = pivot_wider(dat, names_from="type", values_from="v") # split trait and abundance values into two columns
    dat = mutate(dat, :tl => ifelse.(:species .> SR, "C", "R")) # trophic level (tl): species with index greater than SR are consumers ("C"), the rest are resources ("R")
    return dat # return tidy table
end


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
Temperature as a function of space (x), time (t), and some climate parameters
"""
function Temp(x, t, tE, Cmax, Cmin, Tmax, Tmin)
    T = (Tmax-Tmin)*x+Tmin+((Cmin-Cmax)*x+Cmax)*smoothstep(t/tE)
    return T
end

"""
funcresp(n::Vector{<:Real}, Th::Vector{<:Real}, arate::Vector{<:Real}, W::Matrix{<:Real})
Type II functional response
Input:
- n: Vector of population densities of all species in a given patch
- Th: Vector of handling times (with dummy values for resource species)
- arate: Vector of attack rates (with dummy values for resource species)
- W: Adjacency matrix of trophic network; W(i,j)=1 if i eats j and 0 otherwise
Output:
- A matrix F(i,j), the feeding rate of consumer i on resource j
"""
function funcresp!(F, n::Vector{<:Real}, Th::Vector{<:Real}, arate::Vector{<:Real}, W::Matrix{<:Real})
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
    return F
end

"""
    eqs(u, pars, t)

Right-hand side of dynamical equations
 Input:
 - time: Time at which function is evaluated (explicit time-dependence)
 - u: Vector of u variables, with 2*S*L entries, where S is the number
 of species and L the number of patches. The first S*L entries are the
 densities, the second S*L entries are the trait means.
 - pars: Model parameters, given as members of a list
 Output:
 - The derivatives of the densities and trait means, as a vector in a list */
"""
function eqs!(dudt, u, pars, t)
    # Parameters
    @unpack S, SC, SR, L, eta, nmin, venv, tE, Cmax, Cmin, Tmax, Tmin, aw, bw, kappa, d, s, Th, rho, arate, eps, vmat, W, mig, a, model, x = pars
    V = s # convention from Akesson

    # Variables
    sumgr = 0.0
    summig = 0.0
    bsumgr = 0.0
    bsummig = 0.0

    # Assign matrices required; calculate local temperatures
    # TODO: for optimization, those matrices should be stored
    F = zeros(eltype(u), S, S)
    alpha = zeros(eltype(u), S, S)
    beta = zeros(eltype(u), S, S)
    n = @view u[1:SC+SR,:]
    m = @view u[SC+SR+1:end,:]
    dndt = @view dudt[1:SC+SR,:]
    dmdt = @view dudt[SC+SR+1:end,:]
    # initialising vectors
    dndt .= 0.
    dmdt .= 0.

    T = Temp.(x, t, tE, Cmax, Cmin, Tmax, Tmin) # Vector of temperatures

    # Assign competition coeffs alpha_ij^k and selection pressures beta_ij^k
    for k = 1:L
        # If we have temperature-dependent competition:
        if (model == "Tdep") || (model == "Tdep_trophic")
            for i = 1:(SR-1)
                alpha[i,i] = eta / sqrt(2.0 * V[i] + 2.0 * V[i] + eta^2)
                for j = (i+1):SR
                    Omega = 2.0 * V[i] + 2.0 * V[j] + eta^2
                    dm = m[j,k] - m[i,k]
                    alpha[i,j] = eta * exp(-dm*dm/Omega)/sqrt(Omega)
                    alpha[j,i] = alpha[i,j]
                    beta[i,j] = 2.0 * V[i] * alpha[i,j] * dm / Omega
                    beta[j,i] = - beta[i,j] * V[j]/V[i]
                end
            end
            alpha[SR, SR] = eta/sqrt(2.0*V[SR]+2.0*V[SR]+eta^2)
        else # If no temperature-dependent competition, it's much simpler:
            alpha .= a
        end

        funcresp!(F, n[:, k], Th, arate, W) # Feeding rate of species i on j in patch k
        
        # For debugging
        # if k == 1
        #     @show F[1,1]
        #     @show alpha[1,1]
        #     @show beta[1,3]
        #     # @show x
        #     @show n[1,3]
        # end

        # Loop over species
        for i = 1:S
            sumgr = 0.0
            bsumgr = 0.0
            # Species interaction terms in density and then trait evolution equations
            for j = 1:S
                sumgr += -n[i, k] * alpha[i, j] * n[j, k] + eps[i] * n[i, k] * F[i, j] - n[j, k] * F[j, i]
                bsumgr += beta[i, j] * n[j, k]
            end
            summig = 0.0
            bsummig = 0.0
            if L > 1 # otherwise, no spatial structure
                # Dispersal terms in density and then trait evolution equations
                for l in 1:L
                    summig += mig[k, l] * n[i, l] - n[i, k] * mig[l, k]
                    bsummig += mig[k, l] * n[i, l] * (m[i, l]-m[i, k]) / (n[i, k]+1.0e-10)
                end
                # Growth terms in the equations
                summig *= d[i]
                bsummig *= d[i]
            end
            w = bw - aw * m[i, k]
            sw = w^2 + V[i]
            ef = rho[i] * exp(-(T[k]-m[i, k]) * (T[k]-m[i, k])/(2.0*sw)) / sqrt(sw)
            b = ef - kappa
            g = ef * V[i] * (T[k]-m[i, k]) / sw
            q = vmat[i, k] * smoothstep(n[i, k]/nmin)

            h2 = q ./ (q .+ venv) # Heritability

            # Assign calculated rates to vector of derivatives for output
            dndt[i, k] = (n[i, k]*b+sumgr)*smoothstep(n[i, k]/1.0e-6) + summig
            dmdt[i, k] = h2 * (g - bsumgr + bsummig)
            
            # For debugging
            # if i == 1 && k == 1
            #     @show bsummig
            #     @show bsumgr
            #     @show summig
            #     @show sumgr
            # end

            # Periodic boundary conditions
            if L > 1 # otherwise, no spatial structure
                if (k == 1)
                    dndt[i, k] += d[i]*(mig[1, 2]*n[i, 2] - mig[2, 1] * n[i, 1])
                    dmdt[i, k] += d[i] * h2 * mig[1, 2] * n[i, 2] * (m[i, 2]-m[i, 1]) / (n[i, 1]+1.0e-10)
                elseif (k == L)
                    # update dudt for density
                    dndt[i, k]  += d[i] * (mig[k, k-1] * n[i, k-1] - mig[k-1, k] * n[i,k])
                    # update dudt for trait evolution
                    dmdt[i, k] += d[i] * h2 * mig[k,k-1] * n[i,k-1] * (m[i,k-1] - m[i,k]) / (n[i,k] + 1.0e-10)
                end
            end
        end
    end
end