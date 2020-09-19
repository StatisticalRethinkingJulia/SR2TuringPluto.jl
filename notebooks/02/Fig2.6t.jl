### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 4eeb8e68-f925-11ea-010f-af9a137e7d8d
using Pkg, DrWatson

# ╔═╡ 4eebc90a-f925-11ea-3722-2328667bfd1b
begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
end

# ╔═╡ 0fbb31d0-f925-11ea-32ec-7dcffe6e6bde
md"## Fig 2.6t.jl"

# ╔═╡ 4ef73c18-f925-11ea-0289-fd955d24279b
md"##### Define a grid."

# ╔═╡ 4efb4ef2-f925-11ea-2ce0-39385c5362a5
begin
	N = 201
	p_grid = range( 0 , stop=1 , length=N )
end

# ╔═╡ 4f07c86c-f925-11ea-027e-15950c91a091
md"##### Define three priors."

# ╔═╡ 4f0870a0-f925-11ea-2ea9-afb62e0ac34d
begin
	prior = []
	append!(prior, [pdf.(Uniform(0, 1), p_grid)])
	append!(prior, [[p < 0.5 ? 0 : 1 for p in p_grid]])
	append!(prior, [[exp( -5*abs( p - 0.5 ) ) for p in p_grid]])
end;

# ╔═╡ 4f13916a-f925-11ea-1d83-592895f10cb4
likelihood = pdf.(Binomial.(9, p_grid), 6)

# ╔═╡ 4f1b63fe-f925-11ea-380d-fdff386cb1f9
p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 9);

# ╔═╡ 4f1ebcf0-f925-11ea-0e7e-f3e89dc72554
for i in 1:3
  j = (i-1)*3 + 1
  p[j] = plot(p_grid, prior[i], leg=false, ylims=(0, 1), title="Prior")
  p[j+1] = plot(p_grid, likelihood, leg=false, title="Likelihood")
  p[j+2] = plot(p_grid, prior[i].*likelihood, leg=false, title="Posterior")
end

# ╔═╡ 4f2bb452-f925-11ea-3d32-9f5a9a47485d
plot(p..., layout=(3, 3))

# ╔═╡ 4f387372-f925-11ea-37c5-1f17b6a64d1e
md"## End of Fig 2.6t.jl"

# ╔═╡ Cell order:
# ╟─0fbb31d0-f925-11ea-32ec-7dcffe6e6bde
# ╠═4eeb8e68-f925-11ea-010f-af9a137e7d8d
# ╠═4eebc90a-f925-11ea-3722-2328667bfd1b
# ╟─4ef73c18-f925-11ea-0289-fd955d24279b
# ╠═4efb4ef2-f925-11ea-2ce0-39385c5362a5
# ╟─4f07c86c-f925-11ea-027e-15950c91a091
# ╠═4f0870a0-f925-11ea-2ea9-afb62e0ac34d
# ╠═4f13916a-f925-11ea-1d83-592895f10cb4
# ╠═4f1b63fe-f925-11ea-380d-fdff386cb1f9
# ╠═4f1ebcf0-f925-11ea-0e7e-f3e89dc72554
# ╠═4f2bb452-f925-11ea-3d32-9f5a9a47485d
# ╟─4f387372-f925-11ea-37c5-1f17b6a64d1e
