# Clip-03-02-19.jl

using DrWatson
@quickactivate "StstisticalRethinkingTuring"
using StatisticalRethinking

# %% 3.20 - 3.23

pdf.(Binomial(2, 0.7), 0:2)

rand(Binomial(2, 0.7), 1)

rand(Binomial(2, 0.7), 10)

N = 10_000
dummy_w = rand(Binomial(2, 0.7), N)
counts(dummy_w) ./ N |> display   # either
countmap(dummy_w)                 # or, depending on what you want

# %% 3.24, 3.25

dummy_w = rand(Binomial(9, 0.7), N)
histogram(dummy_w)

w = rand(Binomial(9, 0.6), N)

# %% 3.26

samples = let

    # %% 3.2 - 3.5
    
    n = 1000
    p_grid = range(0, 1, length = n)
    prob_p = ones(n)
    prob_data = @. pdf(Binomial(9, p_grid), 6)
    posterior = prob_data .* prob_p
    posterior ./= sum(posterior)
    weights = pweights(posterior)
    sample(p_grid, weights, 10_000)
end

# Version 1: Imagine this as having all the binomial distributions from the p of the
# samples added up. The random numbers are drawn from the resulting distribution.
mix = MixtureModel(Binomial.(9, samples))
w1 = rand(mix, 10_000)

# Version 2: From each p in the sample you make a binomial distribution and draw 1 value
# from it. Then you put them all into one big array.
w2 = rand.(Binomial.(9, samples), 1)
w2 = vcat(w2...)

# These are two different approaches but give the same result for lots of values drawn.
# However, if you are only drawing a few values (less than samples) the first version
# is preferrable since you can still use all the information in the samples.

histogram(w1, normalize = :probability, alpha = 0.6)
histogram!(w2, normalize = :probability, alpha = 0.6)

# End of clip-03-02-19.jl

