### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ f3a59724-187d-11eb-2f7c-2df2f3339a6e
using Pkg, DrWatson

# ╔═╡ f3a5da02-187d-11eb-24ed-a7c460c0219a
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ f3b3b5aa-187c-11eb-2e75-07bc5598cd16
md"## Turing intro clip 3_quap.jl"

# ╔═╡ f3a66514-187d-11eb-22eb-47f74d2724cd
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)     # either
end;

# ╔═╡ 2f0ce4de-187e-11eb-296e-8b4e16814350
@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* weights
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

# ╔═╡ 33b27276-187e-11eb-26a8-1bb7ad288539
begin
	m4_3t = m4_3(df.weight, df.height)
    sampler=NUTS(0.65)
	nsamples=2000
	nchains=4
    chns4_3t = mapreduce(c -> sample(m4_3t, sampler, nsamples), chainscat, 1:nchains)
end

# ╔═╡ 7779262e-187e-11eb-0441-37f564373c1e
q4_3t = quap(m4_3t)

# ╔═╡ 3efadbfe-1882-11eb-0900-b17b00bf8f31
begin
	scatter(df.weight, df.height, leg=false)
	quap4_3t = DataFrame(rand(q4_3t.distr, 10_000)', q4_3t.params)
	a_map = mean(quap4_3t.a)
	b_map = mean(quap4_3t.b)
	plot!(x, a_map .+ b_map .* x)
end

# ╔═╡ Cell order:
# ╟─f3b3b5aa-187c-11eb-2e75-07bc5598cd16
# ╠═f3a59724-187d-11eb-2f7c-2df2f3339a6e
# ╠═f3a5da02-187d-11eb-24ed-a7c460c0219a
# ╠═f3a66514-187d-11eb-22eb-47f74d2724cd
# ╠═2f0ce4de-187e-11eb-296e-8b4e16814350
# ╠═33b27276-187e-11eb-26a8-1bb7ad288539
# ╠═7779262e-187e-11eb-0441-37f564373c1e
# ╠═3efadbfe-1882-11eb-0900-b17b00bf8f31
