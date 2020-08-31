# Clip-04-01-05.jl

using DrWatson
@quickactivate "StatisticalRethinkingTuring"
using StatisticalRethinking

# %% 4.1

pos = [sum(rand(Uniform(-1, 1), 16)) for _ in 1:1000]

# %% 4.2

prod(1 .+ rand(Uniform(0.0, 0.1), 12))

# %% 4.3

growth = [prod(1 .+ rand(Uniform(0.0, 0.1), 12)) for _ in 1:10_000]
density(growth)

# %% 4.4

big = [prod(1 .+ rand(Uniform(0.0, 0.5), 12)) for _ in 1:10_000]
small = [prod(1 .+ rand(Uniform(0.0, 0.01), 12)) for _ in 1:10_000]

# %% 4.5

log_big = [prod(1 .+ log.(rand(Uniform(0.0, 0.5), 12))) for _ in 1:10_000]
density(log_big)
density(log_big, xscale = :log10)

# End of clip-04-01-05.jl
