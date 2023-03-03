using ParametricModels
using UnPack
using DocStringExtensions

"""
    $SIGNATURES

This model is inspired from [Akesson et al.
2021](https://www.nature.com/articles/s41467-021-24977-x).

Specialized variants are provided in other `AbstractModel`s. Best constructure
with `init_akesson_model`.
"""
Base.@kwdef struct AkessonModel{MP,WG,CP,Tr,Fr,A,B} <: AbstractModel
    mp::MP # standard model parameters
    width_growth::WG#width growth type
    competition::CP #competition type
    trophic::Tr #whether there is a trophic level or not
    # TODO: matrix stored for efficiency should be deleted, as they impede the model to be differentiable
    feeding_rates::Fr # stored matrix for efficiency
    alpha::A # stored matrix for efficiency
    beta::B # stored matrix for efficiency
    landscape::Landscape
    temp::Temperature
    SR::Int64 # number of resources
    SC::Int64 # number of consumers
end


"""
$SIGNATURES

# Arguments
- `width_growth::WidthGrowth`: an instance of the WidthGrowth struct representing the width growth function.
- `competition::Competition`: an instance of the Competition struct representing the competition function.
- `trophic::Tr`: an instance of the Trophic struct representing the trophic interaction. Defaults to Trophic{false}() representing a single trophic level.
- `SR::Int`: the number of species in each patch. Defaults to 50.
- `SC::Int`: the number of patches. Defaults to 0 (only one patch).
- `kwargs`: additional keyword arguments to pass to the Base.ODEProblem constructor.

# Example

Here we simulate Akesson model in the most simple version with one trophic
level, one patch, and where competition is temperature dependent. The code
starts by setting the simulation time frame, choosing a differential equation
solver algorithm, and defining the tolerance levels for the solver. The model is
then initialized with various parameters, including the landscape size, carrying
capacity, and temperature dependency of competition. The temperature over a
period of 6500 years is plotted using matplotlib. The model is then simulated
before and after climate change, with the population density and population mean
trait (optimal temperature) plotted for each time step.

```julia
using EcoEvoModelZoo
using PythonCall
using OrdinaryDiffEq
using ParametricModels
plt = pyimport("matplotlib.pyplot")

tstart = -4000 # starting time (relative to start of climate change at t = 0)
tend = 2500 # time at which integration ends
alg = Tsit5()
abstol = 1e-6
reltol = 1e-6
# integrate prob with Julia
tspan = (tstart, 0)
tsteps = tspan[1]:200:tspan[2]
L = 10
SR = 50
temp = Temperature()
land = Landscape(L)

model = init_akesson_model(;mp = ModelParams(;tspan, alg, reltol, abstol, saveat = tsteps),
                            land = Landscape(L),
                            temp,
                            SR = SR,
                            width_growth = WidthGrowth{:TraitDep}(), 
                            competition = Competition{:TraitDep}(), 
                            trophic= Trophic{false}(), 
                            )

# Plotting temperature over 6500 years, from -4000 to 2500
fig, ax = plt.subplots(1)
ts = -100:1:400
p = get_p(model)
temp_ts = [temp.(land.x, t)[1] for t in ts]
ax.plot(ts, temp_ts)
display(fig)

#################
### BEFORE Climate Change ###
#################
sol = simulate(model)
ode_data = Array(sol)

# initial temperature
temp.(land.x, 100.)

# plotting evolution before CC
fig, ax = plt.subplots(1)
l = 1 # patch index
plotting_N_through_time(ax, ode_data, tsteps, l, model)
# ax.set_yscale("log")
display(fig)


#################
### After Climate Change  ###
#################
tspan = (0, 300)
step = 1
tsteps = tspan[1]:1.:tspan[2]
sol = simulate(model, tspan = tspan, saveat = tsteps, u0 = ode_data[:,:,end])
ode_data = Array(sol)

# plotting evolution after CC
fig, axs = plt.subplots(1, 2)
ts = 1:length(tsteps)
l = 1 # patch index
plotting_N_through_time(axs[0], ode_data, tsteps, l, model)
axs[0].set_ylabel("Population density")
# axs[0].set_yscale("log")
plotting_mu_through_time(axs[1], ode_data, tsteps, l, model)
axs[1].set_ylabel("Population mean trait\n(Optimal temperature)")
fig.tight_layout()
```

"""
function init_akesson_model(;mp = ModelParams(),
                            land = Landscape(50),
                            temp = Temperature(),
                            width_growth = WidthGrowth{:TraitDep}(), 
                            competition = Competition{:TraitDep}(), 
                            trophic::Tr = Trophic{false}(), 
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
    if isnothing(mp.p) && isnothing(mp.u0)
        pars, u0 = init_params_akesson_model(land, temp, SR, SC; width_growth, competition, trophic, kwargs...)
        pars = NamedTuple([pair for pair in pars])
        mp = ParametricModels.remake(mp, u0 = u0, p = pars)
    end
    # initializing matrices
    feeding_rates = zeros(eltype(mp.u0), S, S)
    alpha = zeros(eltype(mp.u0), S, S)
    beta = zeros(eltype(mp.u0), S, S)
    AkessonModel(mp, width_growth, competition, trophic, feeding_rates, alpha, beta, land, temp, SR, SC)
end

get_nb_species(em::AkessonModel) = em.SR + em.SC

function (em::AkessonModel)(dudt, u, pars, t)
    # Parameters
    @unpack nmin, venv, kappa, d, v, rho = pars
    @unpack width_growth, competition, trophic, feeding_rates, alpha, beta, SR, SC, landscape, temp = em
    @unpack mig, x, L = landscape

    S = get_nb_species(em)

    # internal variable declaration
    sumgr = 0.0
    summig = 0.0
    bsumgr = 0.0
    bsummig = 0.0
    V = v .+ venv # species' total phenotypic variances


    # Assign matrices required; calculate local temperatures

    n = @view u[1:SC+SR,:]
    m = @view u[SC+SR+1:end,:]
    dndt = @view dudt[1:SC+SR,:]
    dmdt = @view dudt[SC+SR+1:end,:]
    # initialising vectors
    dndt .= 0.
    dmdt .= 0.

    T = temp.(x, t) # Vector of temperatures

    # Assign competition coeffs alpha_ij^k and selection pressures beta_ij^k
    for k = 1:L
        nk = @view n[:, k]
        mk = @view m[:, k]
        # _update_alpha_beta!(mk, em, V, pars)
        funcresp!(feeding_rates, nk, pars, trophic) # Feeding rate of species i on j in patch k
        wk = _get_width_growth_curve(pars, mk, width_growth)
        sw = wk.^2 .+ V

        # Loop over species
        for i = 1:S
            sumgr = 0.0
            bsumgr = 0.0
            # Species interaction terms in density and then trait evolution equations
            for j = 1:S
                alpha_ij, beta_ij = get_alpha_beta(mk[i], mk[j], V[i], V[j], pars, em)
                sumgr += -n[i, k] * alpha_ij * n[j, k] + _feeding(nk, feeding_rates, i, j, pars, trophic)
                bsumgr += beta_ij * n[j, k]
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

            ef = rho[i] * exp(-(T[k]-m[i, k]) * (T[k]-m[i, k])/(2.0*sw[i])) / sqrt(sw[i])
            b = ef - kappa[]
            g = ef * V[i] * (T[k]-m[i, k]) / sw[i]
            q = v[i] * smoothstep(n[i, k]/nmin)

            h2 = q / (q + venv[]) # Heritability

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
  
 function _get_width_growth_curve(p, mk, ::WidthGrowth{:TraitDep})
     @unpack bw, aw = p
     w = bw .- aw .* mk
     return w
 end
 
 function _get_width_growth_curve(p, mk, ::WidthGrowth{:Standard})
     @unpack w = p
     return w
 end
 
#  TODO: not used anymore
#  function _update_alpha_beta!(mk, em::AkessonModel{MP,WG,CP,Tr,Fr,A,B}, V, p) where {MP,WG,Tr,Fr,A,B, CP <: Competition{:TraitDep}}
#      @unpack eta = p
#      @unpack alpha, beta, SR = em
#      for i = 1:(SR-1)
#          alpha[i,i] = eta[] / sqrt(2.0 * V[i] + 2.0 * V[i] + eta[]^2)
#          for j = (i+1):SR
#              Omega = 2.0 * V[i] + 2.0 * V[j] + eta[]^2
#              dm = mk[j] - mk[i]
#              alpha[i,j] = eta[] * exp(-dm*dm/Omega)/sqrt(Omega)
#              alpha[j,i] = alpha[i,j]
#              beta[i,j] = 2.0 * V[i] * alpha[i,j] * dm / Omega
#              beta[j,i] = - beta[i,j] * V[j]/V[i]
#          end
#      end
#      alpha[SR, SR] = eta[]/sqrt(2.0*V[SR]+2.0*V[SR]+eta[]^2);
#  end
 
#  function _update_alpha_beta!(mk, em::AkessonModel{MP,WG,CP,Tr,Fr,A,B}, V, p)  where {MP,WG,Tr,Fr,A,B, CP <: Competition{:Standard}}
#     @unpack alpha, beta = em
#     alpha .= p.a
#     beta .= 0.
#  end


 function get_alpha_beta(mki, mkj, Vi, Vj, p, em::AkessonModel{MP,WG,CP,Tr,Fr,A,B}) where {MP,WG,Tr,Fr,A,B, CP <: Competition{:TraitDep}}
    @unpack eta = p
    @unpack SR = em
    Omega_ij = 2.0 .* Vi .+ 2.0 .* Vj .+ eta[]^2
    dm = mkj .- mki
    alpha_ij = eta[] * exp.(- dm^2 / Omega_ij) / sqrt.(Omega_ij)
    beta_ij = 2.0 * Vi * alpha_ij * dm / Omega_ij
    return alpha_ij, beta_ij
 end

 function get_alpha_beta(mki, mkj, Vi, Vj, p, em::AkessonModel{MP,WG,CP,Tr,Fr,A,B}) where {MP,WG,Tr,Fr,A,B, CP <: Competition{:Standard}}
    return p.a, 0.
 end

 function _feeding(nk, feeding_rates, i, j, p, ::Trophic{true})
    @unpack eps = p
    eps[i] * nk[i] * feeding_rates[i, j] - nk[j] * feeding_rates[j, i]
 end

 function _feeding(n, feeding_rates, i, j, p, ::Trophic{false})
    0.
 end