
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## m0.1t.jl"

md"##### Define a simple Normal model with unknown mean and variance."

@model gdemo(x, y) = begin
  s ~ InverseGamma(2, 3)
  m ~ Normal(0, sqrt(s))
  x ~ Normal(m, sqrt(s))
  y ~ Normal(m, sqrt(s))
end

md"#####  Run sampler, collect results."

chns = mapreduce(c -> sample(gdemo(1.5, 2), NUTS(0.65), 2000), chainscat, 1:4)

md"##### Plot the chains."

plot(chns; seriestype=:traceplot)

plot(chns; seriestype=:density)

md"## End m0.1t.jl"

