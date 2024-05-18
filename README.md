## Note

After many years I have decided to step away from my work with Stan and Julia. My plan is to be around until the end of 2024 for support if someone decides to step in and take over further development and maintenance work.

At the end of 2024 I'll archive the different packages and projects included in the Github organisations StanJulia, StatisticalRethingJulia and RegressionAndOtherStoriesJulia if no one is interested (and time-wise able!) to take on this work.

I have thoroughly enjoyed working on both Julia and Stan and see both projects mature during the last 15 or so years. And I will always be grateful for the many folks who have helped me on numerous occasions. Both the Julia and the Stan community are awesome to work with! Thanks a lot!

## Purpose of SR2TuringPluto.jl

As stated many times by the author in his [online lectures](https://www.youtube.com/playlist?list=PLDcUM9US4XdMROZ57-OIRtIK0aOynbgZN), StatisticalRethinking is a hands-on course. This project is intended to assist with the hands-on aspect of learning the key ideas in Statistical Rethinking. 

SR2TuringPluto is a Julia project that uses Pluto notebooks for this purpose. Max Lapan has a [version using Jupyter](https://github.com/Shmuma/rethinking-2ed-julia). Many of the Pluto notebooks have been derived from Max Lapan's work!

This Julia project uses Turing as the underlying mcmc implementation.  A companion project ( [SR2StanPluto.jl](https://github.com/StatisticalRethinkingJulia/SR2StanPluto.jl) ) uses Stan.

Each notebook demonstrates Julia versions of `code snippets` and `mcmc models` contained in the R package "rethinking" associated with the book [Statistical Rethinking](https://xcelab.net/rm/statistical-rethinking/) by Richard McElreath. 

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
$ git checkout v4.5.0
```

The next step assumes your `base` Julia setup includes at least `Pkg` and `Pluto`.

2. Start a Pluto notebook server in the Julia REPL:
```
$ julia

julia> using Pluto
julia> Pluto.run()
```

3. A Pluto page should open in a browser.

## Usage

Select a notebook in the `open a file` entry box, e.g. type `./` and step to `./notebooks/Chapter_00.jl`.

All "rethinking" data files are stored and maintained in StatisticalRethinking.jl and can be accessed via `sr_datadir(...)`.

This leads to a typical set of opening lines in each notebook:
```
using Pkg

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

1. ppl5_1            : Turing model
1. m5_1t             : The instantiated Turing model (includes data)

Chain(s):

3. chns5_1t          : MCMCChains object (4000 samples from 4 chains)

Results as a DataFrame:

4. prior5_1t_df      : Prior samples (DataFrame)
5. post5_1t_df       : Posterior samples (DataFrame)
6. quap5_1t_df       : MAP approximation to posterior samples (DataFrame)
7. pred5_1t_df       : Posterior predictions (DataFrame)

As before, the `t` at the end of the model number indicates Turing.

**Note**: Naming is not yet consistent through all notebooks. Work in progress!

## Status

SR2TuringPluto.jl is compatible with the 2nd edition of the book.

StructuralCausalModels.jl and ParetoSmoothedImportanceSampling.jl are included as experimental dependencies in the StatisticalRethinking.jl package. Definitely work in progress!

Max Lapan added a package Dagitty.jl which covers options similar as available in StructuralCausalModels.jl. There is also a new package, ParetoSmooth.jl which overlaps with ParetoSmoothedImportanceSampling.jl.
As terminology differs from the terminology used in the Statistical Rethinking book, I have not used this package in the notebooks (yet?).

Any feedback is appreciated. Please open an issue.

## Acknowledgements

Of course, without the excellent textbook by Richard McElreath, this package would not have been possible. The author has also been supportive of this work and gave permission to use the datasets.

This repository is derived from work by Max Lapan, Karajan, previous and current Stan versions of StatisticalRethinking.jl. It has been improved through comments and suggestions by many other contributors.

## Versions

### Version 5.0.0 (Under development, will take time)

1. Complete overhaul. Likely using Makie.jl, Graphs.jl and more.

### Version 4.5.0

1. Adapted Max Lapan's chapters 5 to 14 to an (initial) Pluto format.

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

