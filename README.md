## Purpose of StatisticalRethingTuring.jl

This `project` contains Julia versions of selected `code snippets` and `mcmc models` contained in the R package "rethinking" associated with the book [Statistical Rethinking](https://xcelab.net/rm/statistical-rethinking/) by Richard McElreath.

As stated many times by the author in his [online lectures](https://www.youtube.com/watch?v=ENxTrFf9a7c&list=PLDcUM9US4XdNM4Edgs7weiyIguLSToZRI), StatisticalRethinking is a hands-on course. This project is intended to assist with that aspect of learning the key ideas in StatisticalRethinking.

This project uses Turing as the underlying mcmc implementation.

## Usage

StatisticalRethinkingTuring.jl is a DrWatson project, with some added/re-purposed subdirectories:

1. `models`, contains most used Turing models,
2. `notebooks`, used to store Pluto notebooks and
3. `exercises`, can be used to store the exercises (not stored in the StatisticalRethinkingStan.jl repository)

The `data` directory is only used for locally generated data, exercises, etc.

All example data files are stored and maintained in StatisticalRethinking.jl and can be accessed via `srdatadir()`. 

This leads to a typical set of opening lines in each script:
```
using DrWatson
@quickactivate "StatisticalRethinkingTuring"
using StatisticalRethinking
using Turing                # If Turing is used
include(srcdir("quap.jl"))  # For Turing quap()

# To access e.g. the Howell1.csv data file:
d = CSV.read(srdatadir() * "/Howell1.csv", DataFrame)
d2 = d[d.age .>= 18, :]
```

To (locally) reproduce and use this project, do the following:

0. Download this code base.
1. Move to the downloaded directory.
2. Open a Julia console and, to run the first script, do:
   ```
   julia> include(joinpath(scriptsdir(), "00", "clip-00-01-03.jl")
   ```

This assumes your Julia setup includes `Pkg` and `DrWatson`, activates project `StatisticalrethinkingTuring` and everything should work out of the box.

For the notebooks you'll need to install Pluto.jl and PlutoUI.jl.

## Setup

All R snippets (fragments) have been organized in clips. Clips are named as `clip-cc-fs-ls[s|t|d].jl` where

`cc`      : Chapter number
`fs`      : First snippet in clip
`ls`      : Last snippet in clip
`[s|t|d]` : Mcmc flavor used (s : Stan, t : Turing)

Note: `d` is reserved for a combination Soss/DynamicHMC.

Scripts containing the clips are stored by chapter. In some chapter directories special introductory scripts have been included or scripts that generate figures in the book. These figures are stored by chapter in the `plots` directory.

A similar structure is used for models and Pluto notebooks.

## Status

StatisticalRethinkingTuring.jl is compatible with the 2nd edition of the book.

Expanded coverage of chapters 7 and beyond of the book will likely happen while working on StatisticalRethinkingTuring.jl.

 StructuralCausalModels.jl is included as en experimental dependency in the otherwise stripped down StatisticalRethinking.jl v3.0.0 package.

Any feedback is appreciated. Please open an issue.

## Acknowledgements

This repository and format is derived from work by Karajan, previous versions of StatisticalRethinking.jl and many other contributors.

The huge progress made by the Turing.jl team over the last 2 years, the availability of Julia `projects` in addition to Julia `packages` and the novel approach to notebooks in Pluto.jl were a few of the ideas that triggered exploring a new setup for the StatisticalRethinkingJulia.

## Versions

### Version 0.1.0 (in preparation, expected Oct 2020)

1. Initial version

