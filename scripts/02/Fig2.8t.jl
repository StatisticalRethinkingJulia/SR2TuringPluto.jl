
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Fig2.8t.jl"

@model globethrowing(w, l) = begin
    p ~ Uniform(0, 1)
    w ~ Binomial(w + l, p)
end

begin
	f = Vector{Plots.Plot{Plots.GRBackend}}(undef, 3)
	x = 0:0.01:1

	for (j, i) in enumerate([1, 2, 4])
		w = i * 6
		n = i * 9

		p_grid = range(0, stop=1, length=1000);
		prior = ones(length(p_grid));
		likelihood = [pdf(Binomial(n, p), w) for p in p_grid];
		posterior = likelihood .* prior;
		posterior = posterior / sum(posterior);

		N = 10000
		samples = sample(p_grid, Weights(posterior), N);

		# Analytical calculation

		f[j] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(0.0, 1.0), 
			lab="exact", leg=:topleft, title="n = $n")

		# Quadratic approximation using Turing.jl quap()

		m = globethrowing(w, n-w)
		r = quap(m)
		plot!( f[j], x, pdf.(Normal([r.coef.p][1], √collect(reshape(r.vcov, 1, 1))[1]) , x ), lab="quap")
		
	end
end

plot(f..., layout=(1, 3))

md"## End of Fig2.8t.jl"

