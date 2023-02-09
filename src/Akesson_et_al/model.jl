using ParametricModels
using UnPack
using DocStringExtensions

struct WidthGrowth{T} end # width growth type
struct Trophic{T} end # trophic type
struct Competition{T} end # competition type

"""
    $SIGNATURES

This model is inspired from [Akesson et al. 2021](https://www.nature.com/articles/s41467-021-24977-x).

Specialized variants are provided in other `AbstractModel`s.
"""
Base.@kwdef struct AkessonModel{MP,WG,CP,Tr,Fr,A,B} <: AbstractModel
    mp::MP # standard model parameters
    width_growth::WG#width growth type
    competition::CP #competition type
    trophic::Tr #whether there is a trophic level or not
    feeding_rates::Fr # stored matrix for efficiency
    alpha::A # stored matrix for efficiency
    beta::B # stored matrix for efficiency
    L::Int64 # number of patches
    SR::Int64 # number of resources
    SC::Int64 # number of consumers
end

function init_akesson_model(mp::ModelParams; 
                            width_growth = WidthGrowth{:TraitDep}(), 
                            competition = Competition{:TraitDep}(), 
                            trophic::Tr = Trophic{false}(), 
                            L = 50, 
                            SR = 50, 
                            SC = 0, 
                            kwargs...) where Tr
    if Tr <: Trophic{true}
        @assert SC > 0
    else
        @assert SC == 0
    end
    S = SR + SC # set S to be the total number of species
    # initialising parameters
    pars, u0 = init_params_akesson_model(L, SR, SC; width_growth, competition, trophic, kwargs...)
    pars = NamedTuple([pair for pair in pars])
    mp = ParametricModels.remake(mp, u0 = u0, p = pars)
    # initializing matrices
    feeding_rates = zeros(eltype(u0), S, S)
    alpha = zeros(eltype(u0), S, S)
    beta = zeros(eltype(u0), S, S)
    AkessonModel(mp, width_growth, competition, trophic, feeding_rates, alpha, beta, L, SR, SC)
end

get_nb_species(em::AkessonModel) = em.SR + em.SC

function (em::AkessonModel)(du, u, pars, t)
     # Parameters
     @unpack nmin, venv, tE, Cmax, Cmin, Tmax, Tmin, kappa, d, s, rho, eps, vmat, mig, x = pars
     @unpack width_growth, competition, feeding_rates, alpha, beta, L, SR, SC = em
     S = get_nb_species(em)
     V = s # convention from Akesson
 
     # Variables
     sumgr = 0.0
     summig = 0.0
     bsumgr = 0.0
     bsummig = 0.0
 
     # Assign matrices required; calculate local temperatures

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
        _update_alpha_beta!(em, p)
        funcresp!(feeding_rates, p, trophic) # Feeding rate of species i on j in patch k
 
         # Loop over species
         for i = 1:S
             sumgr = 0.0
             bsumgr = 0.0
             # Species interaction terms in density and then trait evolution equations
             for j = 1:S
                 sumgr += -n[i, k] * alpha[i, j] * n[j, k] + eps[i] * n[i, k] * feeding_rates[i, j] - n[j, k] * feeding_rates[j, i]
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
             w = _get_width_growth_curve(p, m, width_growth)
             sw = w^2 + V[i]
             ef = rho[i] * exp(-(T[k]-m[i, k]) * (T[k]-m[i, k])/(2.0*sw)) / sqrt(sw)
             b = ef - kappa
             g = ef * V[i] * (T[k]-m[i, k]) / sw
             q = vmat[i, k] * smoothstep(n[i, k]/nmin)
 
             h2 = q ./ (q .+ venv) # Heritability
 
             # Assign calculated rates to vector of derivatives for output
             dndt[i, k] = (n[i, k]*b+sumgr)*smoothstep(n[i, k]/1.0e-6) + summig
             dmdt[i, k] = h2 * (g - bsumgr + bsummig)
 
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
  
 function _get_width_growth_curve(p, m, ::WidthGrowth{:TraitDep})
     @unpack bw, aw = p
     w = bw .- aw .* m
     return w
 end
 
 function _get_width_growth_curve(p, m, ::WidthGrowth{:Standard})
     @unpack w = p
     return w
 end
 
 function _update_alpha_beta!(em::AkessonModel{MP,WG,CP}, p) where {MP,WG,CP <: Competition{:TraitDep}}
     @unpack eta, V = p
     @unpack alpha, beta = em
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
     alpha[SR, SR] = eta/sqrt(2.0*V[SR]+2.0*V[SR]+eta^2);
 end
 
 function _update_alpha_beta!(em::AkessonModel{MP,WG,CP}, p)  where {MP,WG,CP <: Competition{:Standard}}
    @unpack alpha, beta = em
     alpha .= p.a
     beta .= 0.
 end