
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Clip-02-01-02t.jl"

md"## snippet 2.1"

begin
	ways = [0, 3, 8, 9, 0]
	ways = ways / sum(ways)
end

md"## snippet 2.2"

md"##### With Distributions.jl working with distributions is a little different from R. Instead of having `rbinom`, `dbinom` and `pbinom` we just got the `Binomial` distribution which to work with."

d = Binomial(9, 0.5)             # Binomial distribution

d.n                              # Number of trials parameter

rand(d)                          # singe random draw

rand(d, 10)                      # 10 random draws

pdf(d, 6)                        # probability density of getting a 6

cdf(d, 6)                        # cumulative probability of getting 6

pdf(Binomial(9, 0.5), 6)         # probability density of getting a 6

md"## End of clip-02-01-02t.jl"

