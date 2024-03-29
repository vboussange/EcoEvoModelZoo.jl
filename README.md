# EcoEvoModelZoo

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://vboussange.github.io/EcoEvoModelZoo.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://vboussange.github.io/EcoEvoModelZoo.jl/dev/)
[![Build Status](https://github.com/vboussange/EcoEvoModelZoo.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vboussange/EcoEvoModelZoo.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/vboussange/EcoEvoModelZoo.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/vboussange/EcoEvoModelZoo.jl)

![](docs/src/time_series_5_species_ecosyste_model.png)

**EcoEvoModelZoo.jl** is a Julia package that offers a collection of eco-evolutionary models written in a consistent style using [ParametricModels.jl](https://github.com/vboussange/ParametricModels.jl). Although these models are inspired from theoretical modelling works, they have been carefully implemented to be compatible with automatic differentiation and can be smoothly integrated with Machine Learning frameworks such as [Flux.jl](https://github.com/FluxML/Flux.jl) and [PiecewiseInference.jl](https://github.com/vboussange/PiecewiseInference.jl). This allows for the testing, refinement and validation of the models using real data from empirical systems.

The long-term goal of EcoEvoModelZoo.jl is to provide high-performance eco-evolutionary models that combine the strengths of theoretical and data-driven models. Our vision is to use EcoEvoModelZoo.jl models in conjunction with ML techniques to enhance our understanding of living systems and provide robust biodiversity forecasts.

## Model zoo
### Spatial eco-evolutionary models
- `EcoEvoGraph`: Eco-evolutionary model on spatial graph, based on [Boussange & Pellissier. 2022](https://www.nature.com/articles/s42003-022-03595-3). 
- `AkessonModel`: Multi-trophic eco-evolutionary model on lattices, based on [Akesson et al. 2021](https://www.nature.com/articles/s41467-021-24977-x).

### 3 species ecosystem model
- `EcosystemModelMcCann`: 3-species chaotic model based on [McCann 1994](http://doi.wiley.com/10.2307/1939558).
- `EcosystemModelOmnivory`: omnivory variant of McCann 1994, based on [McCann 1997](10.1098/rspb.1997.0172).

### N-species ecosystemo model
- `SimpleEcosystemModel`: generalized version of `EcosystemModelMcCann` with N-species. Can be used to reproduce e.g. the 5-species ecosystem model from [Post et al. 2000](https://www.jstor.org/stable/177129?seq=2)

### Resource competition model
- `ResourceCompetition`: plankton ecosystem model based on [Huisman et al. 1999 Nature](http://www.nature.com/articles/46540).
- `ResourceCompetitionSmoothMin`: variant of `ResourceCompetition` where Leibig's law is replaced by imperfect substituable resources (smooth minimum).


## Contributing
Please reach out if you are interested in this project!
<!-- We encourage contributions of new models and documentation. -->

### Share a new model
If you want to share a new model, we recommend following these guidelines:

- Models should be in a separate .jl file and include a test file. They should be of type AbstractModel from ParametricModels.jl.
- Models should include documentation to explain the purpose of the model, how to use it, and any resulting outcomes.
- Code should be concise, neat, and self-explanatory, with minimal boilerplate.
- Please ensure that it is compatible with automatic differentiation.

### Create or improve documentation
You can contribute by:
- Adding or improving documentation for existing models.
- Writing detailed tutorials for a model

## Acknowledgements
The idea of `EcoEvoModelZoo` was inspired by the [machine learning model zoo](https://github.com/FluxML/model-zoo) of Flux.jl.
