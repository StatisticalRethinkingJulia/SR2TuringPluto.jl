# m0.1t.jl

using DrWatson
@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking

# Define a simple Normal model with unknown mean and variance.

@model gdemo(x, y) = begin
  s ~ InverseGamma(2, 3)
  m ~ Normal(0, sqrt(s))
  x ~ Normal(m, sqrt(s))
  y ~ Normal(m, sqrt(s))
end

#  Run sampler, collect results

chns = mapreduce(c -> sample(gdemo(1.5, 2), NUTS(0.65), 2000),
  chainscat, 1:4);

# Summarise results

#chns |> display

# Plot and save results if in ./dev

plot(chns)

# End m0.1t.jl