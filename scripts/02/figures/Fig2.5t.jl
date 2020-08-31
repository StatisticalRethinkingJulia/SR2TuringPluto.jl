# Clip-02-00.1.jl (Fig 2.5)

cd(@__DIR__)
using DrWatson
@quickactivate "StatisticalRethinkingTuring"
using StatisticalRethinking
using Turing

#=
This clip is only intended to generate Fig 2.5 (fig-00.png).
It is not intended to show how to use Stan (yet)!
=#

# 1. Create a Turing model for tossing the globe:

@model globetosses(n, k) = begin
  p ~ Uniform(0, 1)
  k ~ Binomial(n, p)
end

# n will go from 1:9

p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 9)
dens = Vector{DataFrame}(undef, 10)

k = [1,0,1,1,1,0,1,0,1]       # Sequence observed
x = range(0, stop=9, length=10)

for n in 1:9
  p[n] = plot(xlims=(0.0, 1.0), ylims=(0.0, 3.0), leg=false)

  chn = sample(globetosses(n, sum(k[1:n])), NUTS(0.65), 1000)

  # Summarise results (currently requires the master branch from MCMCChains)

  chn |> display
  dfs = DataFrame(chn)

  if n == 1
    hline!([1.0], line=(:dash))
  else
    density!(dens[n][:, :p], line=(:dash))
  end
  density!(dfs[:, :p])
  dens[n+1] = dfs

end

plot(p..., layout=(3, 3))
savefig(plotsdir("02", "figures", "Fig2.5.png"))

# End of clip-02-00.1.jl (Fig 2.5)
