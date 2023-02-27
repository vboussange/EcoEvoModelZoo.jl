module EcoEvoModelZoo
    # Note: we may want to have all these models in submodules
    # which would allow to isolate similar utility functions etcs, defined
    # idiosyncratically
    include("Huisman_et_al.jl")
    export ResourceCompetition, ResourceCompetitionSmoothMin, μsmooth, μ
    include("McCann_et_al.jl")
    export EcosystemModelMcCann, EcosystemModelOmnivory
    include("Akesson_et_al/utils.jl")
    include("Akesson_et_al/model.jl")
    include("Akesson_et_al/init_params.jl")
    include("Boussange_et_al.jl")
    using Requires
    function __init__()
        @require PythonCall="6099a3de-0909-46bc-b1f4-468b9a2dfc0d" begin
            using PythonCall
            include("Akesson_et_al/plotting.jl")
            export plotting_N_through_time, plotting_mu_through_time, plotting_distribution, plotting_distribution_global
        end
    end
    export generate_network, organize_data, smoothstep, Temperature, Landscape, funcresp!, init_akesson_model, AkessonModel, WidthGrowth, Competition, Trophic
end
