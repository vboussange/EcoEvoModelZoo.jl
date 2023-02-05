module EcoEvoModelZoo

    include("Huisman_et_al.jl")
    export ResourceCompetition, ResourceCompetitionSmoothMin, μsmooth, μ
    include("McCann_et_al.jl")
    export EcosystemModelMcCann, EcosystemModelOmnivory
end
