# Clip-04-06.jl

using DrWatson
@quickactivate "StatisticalRethinkingTuring"
using StatisticalRethinking

# %% 4.6
w = 6
n = 9
p_grid = range(0, 1, length = 100)
posterior = @. pdf(Binomial(n, p_grid), w) * pdf(Uniform(0, 1), p_grid)
posterior ./= sum(posterior)

plot(p_grid, posterior)

# End of clip-04-06.jl
