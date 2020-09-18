### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ de60e348-f76f-11ea-0234-7b891b4f2f05
using Pkg, DrWatson

# ╔═╡ de612664-f76f-11ea-32f9-493c9b688096
begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

# ╔═╡ aa07d020-f76f-11ea-2b25-d71cde16c90b
md"## Fig2.7t.jl"

# ╔═╡ de619a2c-f76f-11ea-371d-277e335c97b4
begin
	p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 3)
	N = [5, 20, 50]

	for i in 1:3            # Different priors
		p_grid = range( 0 , stop=1 , length=N[i] )
		prior = pdf.(Uniform(0, 1), p_grid)
		likelihood = [pdf.(Binomial(9, p), 6) for p in p_grid]
		post = (prior .* likelihood) / sum(prior .* likelihood)
		p[i] = plot(p_grid, post, leg=false, title="$(N[i]) points")
		p[i] = scatter!(p_grid, post, leg=false)
	end
end

# ╔═╡ de698b38-f76f-11ea-3820-bf20ccc49d9b
plot(p..., layout=(1, 3))

# ╔═╡ de6a09be-f76f-11ea-1fd7-990715a3e6d7
md"## End of Fig2.7t.jl"

# ╔═╡ Cell order:
# ╟─aa07d020-f76f-11ea-2b25-d71cde16c90b
# ╠═de60e348-f76f-11ea-0234-7b891b4f2f05
# ╠═de612664-f76f-11ea-32f9-493c9b688096
# ╠═de619a2c-f76f-11ea-371d-277e335c97b4
# ╠═de698b38-f76f-11ea-3820-bf20ccc49d9b
# ╟─de6a09be-f76f-11ea-1fd7-990715a3e6d7
