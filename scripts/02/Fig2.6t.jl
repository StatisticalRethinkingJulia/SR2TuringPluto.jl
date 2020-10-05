
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

md"## Fig 2.6t.jl"

md"##### Define a grid."

begin
	N = 201
	p_grid = range( 0 , stop=1 , length=N )
end

md"##### Define three priors."

begin
	prior = []
	append!(prior, [pdf.(Uniform(0, 1), p_grid)])
	append!(prior, [[p < 0.5 ? 0 : 1 for p in p_grid]])
	append!(prior, [[exp( -5*abs( p - 0.5 ) ) for p in p_grid]])
end;

likelihood = pdf.(Binomial.(9, p_grid), 6)

p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 9);

for i in 1:3
  j = (i-1)*3 + 1
  p[j] = plot(p_grid, prior[i], leg=false, ylims=(0, 1), title="Prior")
  p[j+1] = plot(p_grid, likelihood, leg=false, title="Likelihood")
  p[j+2] = plot(p_grid, prior[i].*likelihood, leg=false, title="Posterior")
end

plot(p..., layout=(3, 3))

md"## End of Fig 2.6t.jl"

