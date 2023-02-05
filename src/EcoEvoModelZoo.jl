module EcoEvoModelZoo
    include("Huisman_et_al.jl")
    export ResourceCompetition, ResourceCompetitionSmoothMin, μsmooth, μ
    include("McCann_et_al.jl")
    export EcosystemModelMcCann, EcosystemModelOmnivory
    include("Akesson_et_al.jl")
    export generate_network, organize_data, smoothstep, Temp, funcresp!
end
