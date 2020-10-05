
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Fig3.1t.jl"

begin
	N = 10000
	p_grid = range(0, stop=1, length=N)
	prior = ones(length(p_grid))
	likelihood = pdf.(Binomial.(9, p_grid), 6)
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
	samples = sample(p_grid, Weights(posterior), length(p_grid));
end;

begin
	p1 = scatter(samples, ylim=(0, 1), xlab="Sample number",
	  ylab="Proportion water(p)", leg=false)
	p2 = density(samples, xlim=(0.0, 1.0), ylim=(0.0, 3.0),
	  xlab="Proportion water (p)",
	  ylab="Density", leg=false)
end;

plot(p1, p2, layout=(1,2))

md"## End of Fig3.1t.jl"

