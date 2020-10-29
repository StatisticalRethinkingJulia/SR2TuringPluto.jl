### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 5168c034-1879-11eb-2df4-e74ba7880da1
using Pkg, DrWatson

# ╔═╡ 5168ed8e-1879-11eb-31cd-d3a727fff08d
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing 
	using StatisticalRethinking
	Turing.turnprogress(false);
end;

# ╔═╡ e6220984-1878-11eb-275c-db7fe829037a
md"## Intro clip 1_priors.jl"

# ╔═╡ 51695ecc-1879-11eb-2fd1-b1d0abfb4b8f
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)
	Text(precis(df; io=String))
end

# ╔═╡ 5174d414-1879-11eb-28ed-b3f36412a25f
@model function m4_3(weights, heights)
    a ~ Normal(mean(weights), 50)
    b ~ Normal(0, 5)
    σ ~ Uniform(0, 20)
    μ = a .+ b .* weights
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

# ╔═╡ 11151d8c-1916-11eb-2d70-69180f7be854
begin
	m4_3t = m4_3(df.weight, df.height)
	priors4_3t = sample(m4_3t, Prior(), 50)
end

# ╔═╡ 1f9c6586-1916-11eb-055d-9dd844c9b0a7
begin
	priors4_3t_df = DataFrame(priors4_3t)
	Text(precis(priors4_3t_df; io=String))
end

# ╔═╡ f1943d70-1916-11eb-2cc5-5b6afaea5fbd
begin
	plot(xlims=(25, 70), ylims=(-10, 300), title="Possible prior regression lines for this model", leg=false)
	for row in eachrow(priors4_3t_df)
		plot!(df.weight, row.a .+ row.b * df.weight)
	end
	scatter!(df.weight, df.height)
end

# ╔═╡ 0365d298-187f-11eb-1dc1-b7f61a13c7b3
md"## End of intro clip 1_priors.jl"

# ╔═╡ Cell order:
# ╟─e6220984-1878-11eb-275c-db7fe829037a
# ╠═5168c034-1879-11eb-2df4-e74ba7880da1
# ╠═5168ed8e-1879-11eb-31cd-d3a727fff08d
# ╠═51695ecc-1879-11eb-2fd1-b1d0abfb4b8f
# ╠═5174d414-1879-11eb-28ed-b3f36412a25f
# ╠═11151d8c-1916-11eb-2d70-69180f7be854
# ╠═1f9c6586-1916-11eb-055d-9dd844c9b0a7
# ╠═f1943d70-1916-11eb-2cc5-5b6afaea5fbd
# ╟─0365d298-187f-11eb-1dc1-b7f61a13c7b3
