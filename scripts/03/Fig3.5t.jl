
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Fig3.5t.jl"

begin
  N = 100000
  d = rand(Binomial(9, 0.7), N);
  histogram(d; normalize=:probability, 
    bins=-0.5:1:9.5, leg=false, xticks=0:9, bar_width=0.2)
end

md"## End of Fig3.5t.jl"

