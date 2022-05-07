## Purpose of SR2TuringPluto.jl


As stated many times by the author in his [online lectures](https://www.youtube.com/watch?v=ENxTrFf9a7c&list=PLDcUM9US4XdNM4Edgs7weiyIguLSToZRI), StatisticalRethinking is a hands-on course. This project is intended to assist with the hands-on aspect of learning the key ideas in StatisticalRethinking. 

SR2TuringPluto is a Julia project that uses Pluto notebooks for this purpose. Each notebook demonstrates Julia versions of `code snippets` and `mcmc models` contained in the R package "rethinking" associated with the book [Statistical Rethinking](https://xcelab.net/rm/statistical-rethinking/) by Richard McElreath.

This Julia project uses Turing as the underlying mcmc implementation.  A companion project ( [SR2StanPluto.jl](https://github.com/StatisticalRethinkingJulia/SR2StanPluto.jl) ) uses Stan.

**Note: A new version (v5.0.0) is under development. This is a breaking version. Probably safe to stick with tagged version 4.0.3 until then.**

## Installation

To (locally) reproduce and use this project, do the following:

1. Download this [project](https://github.com/StatisticalRethinkingJulia/SR2TuringPluto.jl) from Github and move to the downloaded directory, e.g.:

```
$ git clone https://github.com/StatisticalRethinkingJulia/SR2TuringPluto.jl SR2TuringPluto
$ cd SR2TuringPluto
```

If you want a specific tagged version, use:

```
$ git tag -l # To see available tags, followed by:
$ git checkout tags/<tag_name> # or simply:
$ git checkout v4.0.3
```

and in the Julia REPL:

```
julia> ]                           # Actvate Pkg mode
(@v1.8) pkg> activate .            # Activate pkg in .
(SR2TuringPluto) pkg> instantiate  # Install in pkg environment
(SR2TuringPluto) pkg> <delete>     # Exit package mode
```

If above procedure fails, if present, try to delete the Manifest.toml file and repeat above steps. **These steps are only needed the first time.**

The next step assumes your Julia setup includes at least `Pkg`, `DrWatson`, `Pluto` and `PlutoUI`.

2. Start a Pluto notebook server.
```
$ julia

julia> using Pluto
julia> Pluto.run()
```

3. A Pluto page should open in a browser.

## Usage

Note: *SR2TuringPluto v4 requires StatisticalRethinking.jl v 4.*

Select a notebook in the `open a file` entry box, e.g. type `./` and step to `./notebooks/TuringGuide.jl`.

SR2TuringPluto.jl is a DrWatson project, with some added/re-purposed subdirectories,

The `data` directory, in DrWatson accessible through `datadir()`, can be used for locally generated data, exercises, etc. All "rethinking" data files are stored and maintained in StatisticalRethinking.jl and can be accessed via `sr_datadir(...)`.

The scripts in the `scripts` subdirectory are directly generated from the notebooks and thus adhere to Pluto's programming restrictions.

This leads to a typical set of opening lines in each notebook:
```
using Pkg, DrWatson

# Note: Below sequence is important. First activate the project
# followed by `using` or `import` statements. Pretty much all
# scripts use StatisticalRethinking. If mcmc sampling is
# needed, it must be loaded before StatisticalRethinking:

using Turing
# more using lines, e.g.
using CSV, DataFrames, Distributions
using StatisticalRethinking
using StatisticalRethinkingPlots, StatsPlots, Plots

# To access e.g. the Howell1.csv data file:
df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
df = df[df.age .>= 18, :]
```

## Naming of models and results:

1. ppl5_1 or m5_1    : Turing model
1. m5_1t             : The instantiated Turing model (includes data)

Chain(s):

3. chns5_1t          : MCMCChains object (4000 samples from 4 chains)

Results as a DataFrame:

4. prior5_1t_df      : Prior samples (DataFrame)
5. post5_1t_df       : Posterior samples (DataFrame)
6. quap5_1t_df       : MAP approximation to posterior samples (DataFrame)
7. pred5_1t_df       : Posterior predictions (DataFrame)

As before, the `t` at the end of the model number indicates Turing.

## Status

SR2TuringPluto.jl is compatible with the 2nd edition of the book.

StructuralCausalModels.jl and ParetoSmoothedImportanceSampling.jl are included as experimental dependencies in the StatisticalRethinking.jl v3 package. Definitely work in progress!

Any feedback is appreciated. Please open an issue.

## Acknowledgements

Of course, without the excellent textbook by Richard McElreath, this package would not have been possible. The author has also been supportive of this work and gave permission to use the datasets.

This repository and format is derived from work by Max Lapan, Karajan, previous and Stan versions of StatisticalRethinking.jl and many other contributors.

## Versions

### Version 5.0.0 (Under development, will take time)

1. Complete overhaul.likely using Makie.jl, Graphs.jl and more.
2. Larger notebooks.
3. Dropped additional figures in early chapters.    r

### Version 4.0.5

1. Minor updates to chapters 1-5.

### Version 4.0.4

1. Package updates.

### Version 4.0.0-4.0.3

1. Switch to StatisticalRethinking v4.
2. Switch to notebooks under development by Max Lapan. Notebooks are being converted to Pluto (vs. Jupyter notebooks).

### versions 2 & 3

1. Many additions for 2nd edition of Statistical Rethinking book.
2. Version 3 uses many ideas proposed in Karajan's scripts

### Version 1.0.0 (in preparation, expected late Nov 2020)

1. Initial version

