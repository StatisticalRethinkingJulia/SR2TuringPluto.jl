
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-02-06-07.jl"

md"## snippet 2.6"

@model globethrowing(w, l) = begin
    p ~ Uniform(0, 1)
    w ~ Binomial(w + l, p)
end

m = globethrowing(6, 3);

r = quap(m)

begin
	p_grid = range(0, 1, length = 20)
	prior = ones(20)
	likelihood = pdf.(Binomial.(9, p_grid), 6)
	posterior = likelihood .* prior
	posterior ./= sum(posterior)
	
	x = 0:0.01:1
	w = 6
	l = 3
	
	f = plot( x, pdf.(Normal(mean(r.coef.p), √r.vcov[1]), x ), lab="quap")
	plot!( x, pdf.(Beta( w+1 , l+1 ) , x ), lab="exact", leg=:topleft, title="n = $(w+l)")
end

md"## End clip-02-06-07t.jl"

