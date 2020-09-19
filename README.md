## Purpose of StatisticalRethingTuring.jl

As stated many times by the author in his [online lectures](https://www.youtube.com/watch?v=ENxTrFf9a7c&list=PLDcUM9US4XdNM4Edgs7weiyIguLSToZRI), StatisticalRethinking is a hands-on course. This project is intended to assist with the hands-on aspect of learning the key ideas in StatisticalRethinking. 

StatisticalRethinkingTuring is a Julia project that uses Pluto notebooks for this purpose. Each notebook demonstrates Julia versions of `code snippets` and `mcmc models` contained in the R package "rethinking" associated with the book [Statistical Rethinking](https://xcelab.net/rm/statistical-rethinking/) by Richard McElreath.

This Julia project uses Turing as the underlying mcmc implementation.

## Usage

StatisticalRethinkingTuring.jl is a DrWatson project, with some added/re-purposed subdirectories:

1. `models`, which contains Turing model scripts if needed repeatedly,
2. `notebooks`, used to store the Pluto notebooks and
3. `exercises`, can be used to store the exercises (not stored in the StatisticalRethinkingTuring.jl repository)

The `data` directory, in DrWatson accessible through `datadir()`, is only used for locally generated data, exercises, etc. All "rethinking" data files are stored and maintained in StatisticalRethinking.jl and can be accessed via `sr_datadir(...)`. 

This leads to a typical set of opening lines in each notebook:
```
using Pkg, DrWatson

# Note: Below sequence is important. First activate the project
# followed by `using` or `import` statements. Pretty much all
# scripts use StatisticalRethinking. If mcmc sampling is
# needed, it must be loaded before StatisticalRethinking:

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking

# To access e.g. the Howell1.csv data file:
d = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
d2 = d[d.age .>= 18, :]
```

To (locally) reproduce and use this project, do the following:

1. Download this [project](https://github.com/StatisticalRethinkingJulia/StatisticalRethinkingTuring.jl) from Github.
2. Move to the downloaded directory.
3. Start a Pluto notebook server.
4. Open a notebook in a browser.

This assumes your Julia setup includes `Pkg`, `DrWatson`, `Pluto`, `PlutoUI` and `Turing`.

Step 3 usually can be done by:
```
$ julia

julia> using Pluto
julia> Pluto.run()
```

By default the Pluto server uses port 1234. In your browser go to
`http://localhost:1234`.

Each notebook will activate the project `StatisticalrethinkingStan`.

## Setup

All R snippets (fragments) have been organized in clips. Each clip is a notebook. Clips are named as `clip-cc-fs-ls[s|t|d].jl` where

* `cc`      : Chapter number
* `fs`      : First snippet in clip
* `ls`      : Last snippet in cli
* `[s|sl|t|d|m]` : Mcmc flavor used (s : Stan, t : Turing)

Note: `d` is reserved for a combination Soss/DynamicHMC, `sl` is reserved for Stan models using the `logpdf` formulation and `m` is reserved for Mamba.

The notebooks containing the clips are stored by chapter.

Special introductory notebooks have been included in `notebooks/intros`, e.g.
`intro-R-users/broadcasting.jl` and `intro-R-users/distributions.jl`.

In addition to clips, in the early notebook chapters (0-3) it is shown how to create the figures in the book, e.g. `Fig2.5t.jl` in `notebooks/chapter/02`.

## Status

StatisticalRethinkingTuring.jl is compatible with the 2nd edition of the book.

StructuralCausalModels.jl is included as en experimental dependency in the StatisticalRethinking.jl v3 package.

Any feedback is appreciated. Please open an issue.

## Acknowledgements

This repository and format is derived from work by Karajan, previous versions of StatisticalRethinking.jl and many other contributors.

The availability of DynamicHMC, the huge progress made by the Turing.jl team over the last 2 years, the availability of Julia `projects` in addition to Julia `packages` and the novel approach to notebooks in Pluto.jl were a few of the ideas that triggered exploring a new setup for the StatisticalRethinkingJulia.

## Versions

### Version 0.1.0 (in preparation, expected Oct 2020)

1. Initial version

