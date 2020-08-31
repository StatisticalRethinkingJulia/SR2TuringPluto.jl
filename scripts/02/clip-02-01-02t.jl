# Clip-02-01-02.jl

using DrWatson
@quickactivate "StatisticalRethinkingTuring"
using StatisticalRethinking

# snippet 2.1
ways = [0, 3, 8, 9, 0]
ways = ways / sum(ways)

# snippet 2.2

#=
With Distributions.jl working with distributions is a little different
that with R. Instead of having `rbinom`, `dbinom` and `pbinom` we just
got the `Binomial` distribution which to work with.
=#

@show d = Binomial(9, 0.5)             # Binomial distribution
d.n |> display
@show rand(d)                          # singe random draw
@show rand(d, 10)                      # 10 random draws
@show pdf(d, 6)                        # probability density of getting a 6
@show cdf(d, 6)                        # cumulative probability of getting 6

@show pdf(Binomial(9, 0.5), 6)         # probability density of getting a 6

# End of clip-02-01-02.jl
