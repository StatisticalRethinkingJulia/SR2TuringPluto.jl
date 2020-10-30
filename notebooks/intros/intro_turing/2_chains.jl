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
md"## Intro clip 2_chains.jl"

# ╔═╡ 51695ecc-1879-11eb-2fd1-b1d0abfb4b8f
begin
	N = 200
	df = DataFrame()
	df.weight = rand(Uniform(20, 90), N)
	f(x) = mean(df.weight) + 1.6x + rand(Normal(0, 10))
	df.height = f.(df.weight)
	Text(precis(df; io=String))
end

# ╔═╡ 5174d414-1879-11eb-28ed-b3f36412a25f
@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* weights
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

# ╔═╡ 517544f8-1879-11eb-3559-81508191d427
begin
	m4_3t = m4_3(df.weight, df.height)
    sampler=NUTS(0.65)
	nsamples=2000
	nchains=4
    chns4_3t = mapreduce(c -> sample(m4_3t, sampler, nsamples), chainscat, 1:nchains)
end

# ╔═╡ 2c0b2e42-187b-11eb-0c00-43aeac8b5045
plot(chns4_3t; seriestype=:traceplot)

# ╔═╡ 9110be74-187b-11eb-34e8-c1d4da25b29e
plot(chns4_3t; seriestype=:density)

# ╔═╡ b1690f50-19fa-11eb-15c7-87b89cf8ffa5
md"##### For additional diagnostics, see the [MCMCChains README document](https://github.com/TuringLang/MCMCChains.jl/blob/master/README.md)"

# ╔═╡ 0365d298-187f-11eb-1dc1-b7f61a13c7b3
md"## End of Intro clip 2_chains.jl"

# ╔═╡ Cell order:
# ╟─e6220984-1878-11eb-275c-db7fe829037a
# ╠═5168c034-1879-11eb-2df4-e74ba7880da1
# ╠═5168ed8e-1879-11eb-31cd-d3a727fff08d
# ╠═51695ecc-1879-11eb-2fd1-b1d0abfb4b8f
# ╠═5174d414-1879-11eb-28ed-b3f36412a25f
# ╠═517544f8-1879-11eb-3559-81508191d427
# ╠═2c0b2e42-187b-11eb-0c00-43aeac8b5045
# ╠═9110be74-187b-11eb-34e8-c1d4da25b29e
# ╟─b1690f50-19fa-11eb-15c7-87b89cf8ffa5
# ╟─0365d298-187f-11eb-1dc1-b7f61a13c7b3
