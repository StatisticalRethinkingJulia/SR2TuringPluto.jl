### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ ac31c1f4-0709-11eb-09ab-41e27b54a47b
using Pkg, DrWatson

# ╔═╡ a42e0abc-0709-11eb-27eb-412153f97b8f
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ b2616f98-0709-11eb-32a2-514954fe4fa5
md"## Clip-04-42-47t.jl"

# ╔═╡ a131667e-0709-11eb-28be-7906b12ad346
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)     # either
end;

# ╔═╡ caba36be-070a-11eb-3536-3d17bead1bcd
md"### snippet 4.42"

# ╔═╡ d2b9c2ee-070a-11eb-09f8-3390f676f682
@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* (weights .- x̄)
    heights ~ MvNormal(μ, σ)
end

# ╔═╡ d2ba4cfa-070a-11eb-394f-57a4000ca479
m4_3t = m4_3(df.weight, df.height)

# ╔═╡ d2bae2e4-070a-11eb-30fc-638d5366be34
quap4_3t = quap(m4_3t, NelderMead())

# ╔═╡ d2cb37f4-070a-11eb-2aee-af17621d1840
md"### snippet 4.43"

# ╔═╡ d2dc3644-070a-11eb-0cd6-19687545a003
@model function m4_3l(weights, heights)
    a ~ Normal(178, 20)
    log_b ~ Normal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ exp(log_b) .* (weights .- mean(weights))
    heights .~ Normal.(μ, σ)
end

# ╔═╡ d2e087b2-070a-11eb-3089-7d6511ec927f
m4_3tl = m4_3l(df.weight, df.height)

# ╔═╡ d2e9d2c4-070a-11eb-2130-6397de3611e4
quap4_3tl = quap(m4_3tl, NelderMead())

# ╔═╡ d2f3e980-070a-11eb-0ac4-a9f0c01b3453
md"### snippets 4.44, 4.45"

# ╔═╡ d30746de-070a-11eb-2f08-d76d77bf0853
round.(quap4_3t.vcov, digits = 3)

# ╔═╡ d308808c-070a-11eb-080c-27e7c48d4fbe
md"### snippet 4.46"

# ╔═╡ d3222082-070a-11eb-2804-b9341057fa9a
begin
	scatter(df.weight, df.height, leg=false)
	post = DataFrame(rand(quap4_3t.distr, 10_000)', quap4_3t.params)
	a_map = mean(post.a)
	b_map = mean(post.b)
	plot!(x, a_map .+ b_map .* (x .- x̄))
end

# ╔═╡ d327a962-070a-11eb-10af-0d4f260b4aef
md"### snippets 4.47"

# ╔═╡ d33379f4-070a-11eb-14b7-39e2ba66bd11
Text(precis(post; io=String))

# ╔═╡ c45fa368-0709-11eb-3013-83becc440a56
md"## End of clip-04-42-47t.jl"

# ╔═╡ Cell order:
# ╠═b2616f98-0709-11eb-32a2-514954fe4fa5
# ╠═ac31c1f4-0709-11eb-09ab-41e27b54a47b
# ╠═a42e0abc-0709-11eb-27eb-412153f97b8f
# ╠═a131667e-0709-11eb-28be-7906b12ad346
# ╟─caba36be-070a-11eb-3536-3d17bead1bcd
# ╠═d2b9c2ee-070a-11eb-09f8-3390f676f682
# ╠═d2ba4cfa-070a-11eb-394f-57a4000ca479
# ╠═d2bae2e4-070a-11eb-30fc-638d5366be34
# ╟─d2cb37f4-070a-11eb-2aee-af17621d1840
# ╠═d2dc3644-070a-11eb-0cd6-19687545a003
# ╠═d2e087b2-070a-11eb-3089-7d6511ec927f
# ╠═d2e9d2c4-070a-11eb-2130-6397de3611e4
# ╟─d2f3e980-070a-11eb-0ac4-a9f0c01b3453
# ╠═d30746de-070a-11eb-2f08-d76d77bf0853
# ╟─d308808c-070a-11eb-080c-27e7c48d4fbe
# ╠═d3222082-070a-11eb-2804-b9341057fa9a
# ╟─d327a962-070a-11eb-10af-0d4f260b4aef
# ╠═d33379f4-070a-11eb-14b7-39e2ba66bd11
# ╟─c45fa368-0709-11eb-3013-83becc440a56
