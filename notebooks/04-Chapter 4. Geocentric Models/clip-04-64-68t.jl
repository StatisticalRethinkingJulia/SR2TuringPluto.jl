### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 81cdc1aa-04fd-11eb-3235-bff0e3f52535
using Pkg, DrWatson

# ╔═╡ 81cdee8c-04fd-11eb-29f5-85927396b6dc
begin
	@quickactivate "SR2TuringPluto"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ 5e06a11c-04f6-11eb-0116-1fd51d614f99
md"## clip-04-64-79t"

# ╔═╡ b8b8de5e-0512-11eb-0932-830c3055e00a
#include(srcdir("tools.jl"))

md"### snippet 4.64"

# ╔═╡ 81ce6d32-04fd-11eb-2c4f-77bdf6a3f9e2
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	scale!(df, [:weight])
end;

# ╔═╡ 81e16d54-04fd-11eb-26e6-5b62d66d3d34
f_parabola(weight_s, a, b1, b2) = a + b1 * weight_s + b2 * weight_s^2

# ╔═╡ 81ee40ec-04fd-11eb-1c71-6b4674695b56
@model function parabola(weight_s, heights)
    a ~ Normal(178, 20)
    b1 ~ LogNormal(0, 1)
    b2 ~ Normal(0, 1)
    σ ~ Uniform(0, 50)
    μ = f_parabola.(weight_s, a, b1, b2)
    heights ~ MvNormal(μ, σ)
end

# ╔═╡ 81f015b6-04fd-11eb-23c4-3f3f8a7201c8
begin
	quap4_5t = quap(parabola(df.weight_s, df.height), NelderMead())
	quap4_5t.coef
end

# ╔═╡ 960b36fe-0573-11eb-0e51-9f7bc0d175a5
quap4_5t.vcov

# ╔═╡ 81f6618c-04fd-11eb-27e7-dffa73914969
md"### snippet 4.67"

# ╔═╡ 3228440c-0527-11eb-00f0-9b38ae90372c
begin
	weight_seq_1 = range(-2.2, 2, length = 30)
	post_1 = DataFrame(rand(quap4_5t.distr, 1_000)', quap4_5t.params)
	mu_1 = f_parabola.(weight_seq_1, post_1.a', post_1.b1', post_1.b2')
	sim_1 = rand.(Normal.(mu_1, post_1.σ'))
	mu_1 = meanlowerupper(mu_1)
	sim_1 = meanlowerupper(sim_1)
end;

# ╔═╡ 82285bce-04fd-11eb-0bfd-27cbf60511ee
md"### snippet 4.68"

# ╔═╡ 93badee6-0527-11eb-16bc-b1c6541276ba
begin
	scatter(df.weight_s, df.height, ms = 3, alpha = 0.7, legend = false)
	plot!(weight_seq_1, mu_1.mean, ribbon = (mu_1.mean .- mu_1.lower, mu_1.upper .- mu_1.mean))
	plot!(weight_seq_1, sim_1.lower, fillrange = sim_1.upper, alpha = 0.3, linealpha = 0.0, c = 2)
end

# ╔═╡ 4f386d5c-0573-11eb-084c-1758217d06f9
md"## End of clip-04-64-68t.jl"

# ╔═╡ Cell order:
# ╟─5e06a11c-04f6-11eb-0116-1fd51d614f99
# ╠═81cdc1aa-04fd-11eb-3235-bff0e3f52535
# ╠═81cdee8c-04fd-11eb-29f5-85927396b6dc
# ╟─b8b8de5e-0512-11eb-0932-830c3055e00a
# ╠═81ce6d32-04fd-11eb-2c4f-77bdf6a3f9e2
# ╠═81e16d54-04fd-11eb-26e6-5b62d66d3d34
# ╠═81ee40ec-04fd-11eb-1c71-6b4674695b56
# ╠═81f015b6-04fd-11eb-23c4-3f3f8a7201c8
# ╠═960b36fe-0573-11eb-0e51-9f7bc0d175a5
# ╟─81f6618c-04fd-11eb-27e7-dffa73914969
# ╠═3228440c-0527-11eb-00f0-9b38ae90372c
# ╠═82285bce-04fd-11eb-0bfd-27cbf60511ee
# ╠═93badee6-0527-11eb-16bc-b1c6541276ba
# ╟─4f386d5c-0573-11eb-084c-1758217d06f9
