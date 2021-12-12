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
md"## clip-04-69-71t"

# ╔═╡ b8b8de5e-0512-11eb-0932-830c3055e00a
#include(srcdir("tools.jl"))

md"### snippet 4.64"

# ╔═╡ 81ce6d32-04fd-11eb-2c4f-77bdf6a3f9e2
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	scale!(df, [:weight])
end;

# ╔═╡ 82407c7c-04fd-11eb-2bba-3783e79b0671
md"### snippet 4.69"

# ╔═╡ da70abc2-0527-11eb-31fa-d5b029193ba5
f_cube(weight_s, a, b1, b2, b3) = a + b1 * weight_s + b2 * weight_s^2 + b3 * weight_s^3

# ╔═╡ 82483e30-04fd-11eb-1290-b98ae5d06a45
@model function cube(weight_s, heights)
    a ~ Normal(178, 20)
    b1 ~ LogNormal(0, 1)
    b2 ~ Normal(0, 10)
    b3 ~ Normal(0, 10)
    σ ~ Uniform(0, 50)
    μ = f_cube.(weight_s, a, b1, b2, b3)
    heights ~ MvNormal(μ, σ)
end

# ╔═╡ 8250f05c-04fd-11eb-0370-899604fffdfd
begin
	quap4_6t = quap(cube(df.weight_s, df.height), NelderMead())
	quap4_6t.coef
end

# ╔═╡ 825e8410-04fd-11eb-3e36-6d3bac027215
md"### snippets 4.70, 4.71"

# ╔═╡ 423cd608-0529-11eb-126e-7fde3de2165b
begin
	weight_seq_2 = range(-2.2, 2, length = 30)
	post_2 = DataFrame(rand(quap4_6t.distr, 1_000)', quap4_6t.params)
	mu_2 = f_cube.(weight_seq_2, post_2.a', post_2.b1', post_2.b2', post_2.b3') |> meanlowerupper
	sim_2 = rand.(Normal.(mu_2.raw, post_2.σ')) |> meanlowerupper
	weight_seq_rescaled = weight_seq_2 .* std(df.weight) .+ mean(df.weight)
	scatter(df.weight, df.height, ms = 3, alpha = 0.7, legend = false)
	plot!(weight_seq_rescaled, mu_2.mean, ribbon = (mu_2.mean .- mu_2.lower, mu_2.upper .- mu_2.mean))
	plot!(weight_seq_rescaled, sim_2.lower, fillrange = sim_2.upper, alpha = 0.3, la = 0.0, c = 2)
end

# ╔═╡ 4f386d5c-0573-11eb-084c-1758217d06f9
md"## End of clip-04-69-71t.jl"

# ╔═╡ Cell order:
# ╟─5e06a11c-04f6-11eb-0116-1fd51d614f99
# ╠═81cdc1aa-04fd-11eb-3235-bff0e3f52535
# ╠═81cdee8c-04fd-11eb-29f5-85927396b6dc
# ╟─b8b8de5e-0512-11eb-0932-830c3055e00a
# ╠═81ce6d32-04fd-11eb-2c4f-77bdf6a3f9e2
# ╟─82407c7c-04fd-11eb-2bba-3783e79b0671
# ╠═da70abc2-0527-11eb-31fa-d5b029193ba5
# ╠═82483e30-04fd-11eb-1290-b98ae5d06a45
# ╠═8250f05c-04fd-11eb-0370-899604fffdfd
# ╟─825e8410-04fd-11eb-3e36-6d3bac027215
# ╠═423cd608-0529-11eb-126e-7fde3de2165b
# ╟─4f386d5c-0573-11eb-084c-1758217d06f9
