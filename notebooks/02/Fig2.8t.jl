### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 9a7e5880-f770-11ea-0a57-5daf86ffbbf1
using Pkg, DrWatson

# ╔═╡ 9a7e9840-f770-11ea-3bc2-e765554021a1
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ 5998e9d4-f770-11ea-3578-ebf8a2c7bc58
md"## Fig2.8t.jl"

# ╔═╡ dda50c28-f9e9-11ea-0778-350b0c22a6d9
@model globethrowing(w, l) = begin
    p ~ Uniform(0, 1)
    w ~ Binomial(w + l, p)
end

# ╔═╡ 9a7f1360-f770-11ea-0c32-ff097eff9095
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

# ╔═╡ 9a8cc690-f770-11ea-22bd-95a007507a6c
plot(f..., layout=(1, 3))

# ╔═╡ 9a8d70ae-f770-11ea-2df0-fdc0f8c29935
md"## End of Fig2.8t.jl"

# ╔═╡ Cell order:
# ╟─5998e9d4-f770-11ea-3578-ebf8a2c7bc58
# ╠═9a7e5880-f770-11ea-0a57-5daf86ffbbf1
# ╠═9a7e9840-f770-11ea-3bc2-e765554021a1
# ╠═dda50c28-f9e9-11ea-0778-350b0c22a6d9
# ╠═9a7f1360-f770-11ea-0c32-ff097eff9095
# ╠═9a8cc690-f770-11ea-22bd-95a007507a6c
# ╟─9a8d70ae-f770-11ea-2df0-fdc0f8c29935
