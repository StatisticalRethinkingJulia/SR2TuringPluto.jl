## Purpose of SR2TuringPluto.jl

As stated many times by the author in his [online lectures](https://www.youtube.com/watch?v=ENxTrFf9a7c&list=PLDcUM9US4XdNM4Edgs7weiyIguLSToZRI), StatisticalRethinking is a hands-on course. This project is intended to assist with the hands-on aspect of learning the key ideas in StatisticalRethinking. 

SR2TuringPluto is a Julia project that uses Pluto notebooks for this purpose. Each notebook demonstrates Julia versions of `code snippets` and `mcmc models` contained in the R package "rethinking" associated with the book [Statistical Rethinking](https://xcelab.net/rm/statistical-rethinking/) by Richard McElreath.

If you prefer to work with scripts instead of notebooks, a utility in the `src` directory is provided (`generate_scripts.jl`) to create scripts from all notebooks and store those in a newly created `scripts` directory. Note that this is a simple tool and will override all files in the `scripts` directory. For exploration purposes I suggest to move some of those scripts to e.g. the `research` directory.

This Julia project uses Turing as the underlying mcmc implementation.  A companion project ( [StatisticalRethinkingStan.jl](https://github.com/StatisticalRethinkingJulia/StatisticalRethinkingStan.jl) ) uses Stan.

## Installation

To (locally) reproduce and use this project, do the following:

1. Download this [project](https://github.com/StatisticalRethinkingJulia/SR2TuringPluto.jl) from Github and move to the downloaded directory, e.g.:

```
$ git clone https://github.com/StatisticalRethinkingJulia/SR2TuringPluto.jl
$ cd SR2TuringPluto.jl
# Note that Turing.jl is usally only supported in released version of Julia.
```

If you want a specific tagged version, use:

```
$ git tag -l # To see available tags, followed by:
$ git checkout tags/<tag_name> # or simply:
$ git checkout v1.1.1
```

and in the Julia REPL:

```
julia> ]                                        # Actvate Pkg mode
(@v1.6) pkg> activate .                         # Activate pkg in .
(SR2TuringPluto) pkg> instantiate  # Install in pkg environment
(SR2TuringPluto) pkg> <delete>     # Exit package mode
```

If above procedure fails, if present, try to delete the Manifest.toml file and repeat above steps. As mentioned above, these steps are only needed the first time.

The next step assumes your Julia setup includes `Pkg`, `DrWatson`, `Pluto` and `PlutoUI`.

2. Start a Pluto notebook server.
```
$ julia

julia> using Pluto
julia> Pluto.run()
```

3. A Pluto page should open in a browser.

## Usage

Note: *SR2TuringPluto v4 requires StatisticalRethinking.jl v 4.*

Select a notebook in the `open a file` entry box, e.g. type `./` and step to `./notebooks/00/clip-00-01-03t.jl`.

SR2TuringPluto.jl is a DrWatson project, with some added/re-purposed subdirectories:

1. `models`, which contains a subset of the Turing language models,
2. `notebooks`, used to store the Pluto notebooks,
3. `scripts`, Julia scrips generated from the notebooks.

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

## Naming conventions

All R snippets (fragments) have been organized in clips. Each clip is a notebook. 

Clips are named as `clip-cc-fs-ls[s|t|d].jl` where

* `cc`      : Chapter number
* `fs`      : First snippet in clip
* `ls`      : Last snippet in cli
* `[s|sl|t|d|m]` : Mcmc flavor used (s : Stan, t : Turing)

Note: `d` is reserved for a combination Soss/DynamicHMC, `sl` is reserved for Stan models using the `logpdf` formulation and `m` is reserved for Mamba.

The notebooks containing the clips are stored by chapter.  In addition to clips, in the early notebook chapters (0-3) it is also shown how to create some of the figures in the book, e.g. `Fig2.5t.jl` in `notebooks/chapter/02`.

Special introductory notebooks have been included in `notebooks/intros`, e.g.
in subdirectories `intro-R-users`, `intro-pluto` and `intro-turing`. They are intended to illustrate ways of using Julia and Pluto and of basic patterns to work with Turing models.

Additional introductory notebooks showing Julia and statistics ( based on the [Statistics with Julia](https://statisticswithjulia.org/index.html) book ) can be found in [StatisticsWithJuliaPlutoNotebooks](https://github.com/StatisticalRethinkingJulia/StatisticsWithJuliaPlutoNotebooks.jl).

One goal for the changes in StatisticalRethinking v3 was to make it easier to compare and mix and match results from different mcmc implementations. Hence consistent naming of models and results is important. The models and the results of simulations are stored as follows:

Models and results:

0. @model            : Turing model
1. m5_1t             : The sampled StanSample model
2. q5_1st            : Stan quap model (NamedTuple similar to Turing)

Draws:

3. chns5_1t          : MCMCChains object (4000 samples from 4 chains)
4. part5_1t          : Stan samples (Particles notation)
5. quap5_1t          : Quap samples (Particles notation)
6. nt5_1t            : NamedTuple with samples values
7. ka5_1t            : KeyedArray object (see AxisArrays.jl)
8. da5_1t            : DimArray object (see DimensionalData.jl)
9. st5_1t            : StanTable 0bject

The default for `read_samples(m1_1s)` is a StanTable chains object.

Results as a DataFrame:

10. prior5_1t_df      : Prior samples (DataFrame)
11. post5_1t_df       : Posterior samples (DataFrame)
12. quap5_1t_df       : Quap approximation to posterior samples (DataFrame)
13. pred5_1t_df       : Posterior predictions (DataFrame)

As before, the `t` at the end of the model number indicates Turing.

## Status

SR2TuringPluto.jl is compatible with the 2nd edition of the book.

StructuralCausalModels.jl and ParetoSmoothedImportanceSampling.jl are included as experimental dependencies in the StatisticalRethinking.jl v3 package. Definitely work in progress!

Any feedback is appreciated. Please open an issue.

## Acknowledgements

Of course, without the excellent textbook by Richard McElreath, this package would not have been possible. The author has also been supportive of this work and gave permission to use the datasets.

This repository and format is derived from work by Karajan, previous versions of StatisticalRethinking.jl and many other contributors.

## Versions

### Version 4.0.0 (Under preparation!)

1. Switch to StatisticalRethinking v4.
2. Switch to notebooks under development by Max Lapan. Notebooks are being converted to Pluto (vs. Jupyter notebooks).

### versions 2 & 3

1. Many additions for 2nd edition of Statistical Rethinking book.
2. Version 3 uses many ideas proposed in Karajan's scripts

### Version 1.0.0 (in preparation, expected late Nov 2020)

1. Initial version

