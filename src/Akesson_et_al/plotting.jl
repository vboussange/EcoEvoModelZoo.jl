using DocStringExtensions
"""
$SIGNATURES

Here we plot both consumers and resources
"""
function plotting_N_through_time(ax, u, ts, l, pars)
    @unpack SR, SC = pars
    @assert size(u, 1) == 2 * (SR + SC) # u stores N and μ
    _cmap = PyPlot.cm.get_cmap("tab20", SR+SC) # only for cb
    color_palette = [_cmap(i) for i in 1:SR+SC]
    for i in 1:SR+SC
        ax.plot(ts, u[i,l,:], c = color_palette[i])
    end
    ax.set_xlabel("time")
    ax.set_ylabel("density")
    ax.set_title("Patch $l")

end

"""
$SIGNATURES

- `ts`: should be indices of `saveat` and not real times
Here we plot both consumers and resources
"""
function plotting_mu_through_time(ax, u, ts, l, pars)
    @unpack SR, SC = pars
    @assert size(u, 1) == 2 * (SR + SC) # u stores N and μ
    _cmap = PyPlot.cm.get_cmap("tab20", SR+SC) # only for cb
    color_palette = [_cmap(i) for i in 1:SR+SC]
    for i in 1:SR+SC
        ax.plot(ts, u[SR+SC + i,l,:], c = color_palette[i])
    end
    ax.set_xlabel("time")
    ax.set_ylabel("mean trait")
    ax.set_title("Patch $l")
end


"""
$SIGNATURES

Here we the phenotypic distribution of each species, where local distributions have been aggregated
"""
function plotting_distribution_global(ax, u, t, i, pars)
    for l in 1:L #looping through space
        plotting_distribution(ax, u, t, i, l, pars)
    end
end

"""
$SIGNATURES

- `t`: should be index of `saveat` and not real times
Here we plot both consumers and resources
"""
function plotting_distribution(ax, u, t, i, l, pars)
    @unpack SR, SC, L, s, Tmin, Tmax = pars
    _cmap = PyPlot.cm.get_cmap("tab20", SR+SC) # only for cb
    color_palette = [_cmap(i) for i in 1:SR+SC]
    xs = range(Tmin - 5., Tmax + 5, length=150)
    Ni = u[i, l, t]
    mui = u[i+SR+SC, l, t]
    ax.fill_between(xs, Ni * pdf.(Normal(mui, s[i]), xs), color = color_palette[i], alpha = 0.4)
end