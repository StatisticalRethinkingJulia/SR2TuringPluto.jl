### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 81cdc1aa-04fd-11eb-3235-bff0e3f52535
using Pkg, DrWatson

# ╔═╡ 81cdee8c-04fd-11eb-29f5-85927396b6dc
begin
	@quickactivate "StatisticalRethinkingTuring"
	#using BSplines: BSplineBasis, basismatrix
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

# ╔═╡ 82a6656e-04fd-11eb-08c9-3bd9d2d8af20
md"### snippet 4.72"

# ╔═╡ 8ad73da4-0529-11eb-3e51-910e74043bfe
df2 = CSV.read(sr_datadir("cherry_blossoms.csv"), DataFrame; missingstring = "NA");

# ╔═╡ af9e141e-0529-11eb-3c28-b986f5939349
scatter(df2.year, df2.doy, leg=false)

# ╔═╡ 82ba1b72-04fd-11eb-33ed-0db4307c66fc
md"### snippet 4.73"

# ╔═╡ d9b4a9fa-0529-11eb-21b5-2da845b95fe9
begin
	df3 = dropmissing(df2, :doy)
	num_knots = 15
	knot_list = quantile(df3.year, range(0, 1, length = num_knots))
	basis = BSplineBasis(4, knot_list)
	B = basismatrix(basis, df3.year)
end;

# ╔═╡ b0b1e776-0575-11eb-093b-f1cb2d9c4b2e
begin
	plot(legend = false, xlabel = "year", ylabel = "basis value")
	for y in eachcol(B)
		plot!(df3.year, y)
	end
	plot!()
end

# ╔═╡ 832a285e-04fd-11eb-061b-49a40155abe1
md"## snippet 4.76"

# ╔═╡ cdb02de2-0575-11eb-121e-df6ecd57ba7c
@model function spline_1(D, B = B)
    α ~ Normal(100, 10)
    w ~ filldist(Normal(0, 10), size(B, 2))
    σ ~ Exponential(1)
    μ = α .+ B * w
    D ~ MvNormal(μ, σ)
    return μ
end

# ╔═╡ 8343fe0a-04fd-11eb-2a3e-79daed2e6b98
md"### snippet 4.77"

# ╔═╡ 1233198e-0576-11eb-3c17-672034589062
begin
	quap4_7t = quap(spline_1(df3.doy))
	quap4_7t.coef
end

# ╔═╡ 833793e0-04fd-11eb-0706-7567eada0aed
begin
	w_str = ["w[$i]" for i in 1:length(basis)]
	post_3 = DataFrame(rand(quap4_7t.distr, 1000)', ["α"; w_str; "σ"])
	w_3 = mean.(eachcol(post_3[:, w_str]))              # either
	w_3 = [mean(post_3[:, col]) for col in w_str]       # or
	plot(legend = false, xlabel = "year", ylabel = "basis * weight")
	for y in eachcol(B .* w_3')
		plot!(df3.year, y)
	end
	plot!()
end	

# ╔═╡ 83b54bb4-04fd-11eb-3208-4162cb3eabc7
md"### snippet 4.78"

# ╔═╡ d9fbed6e-0572-11eb-29bd-1d30ae286045
begin
	mu_3 = post_3.α' .+ B * Array(post_3[!, w_str])'
	mu_3 = meanlowerupper(mu_3)
	scatter(df3.year, df3.doy, alpha = 0.3)
	plot!(df3.year, mu_3.mean, ribbon = (mu_3.mean .- mu_3.lower, mu_3.upper .- mu_3.mean))
end

# ╔═╡ 83f2d7f4-04fd-11eb-3114-d3faa79abb1a
md"### snippet 4.79"

# ╔═╡ 21d2a6a0-0573-11eb-1043-adeebe1e8fde
@model function spline_2(D, B = B)
    α ~ Normal(100, 10)
    w ~ filldist(Normal(0, 10), size(B, 2))
    σ ~ Exponential(1)
    μ = [α + sum(Brow .* w) for Brow in eachrow(B)]
    D ~ MvNormal(μ, σ)
end

# ╔═╡ 840196ae-04fd-11eb-1c8a-3511b18e2835
begin
	quap4_7_altt = quap(spline_2(df3.doy))
	quap4_7_altt.coef
end

# ╔═╡ 4f386d5c-0573-11eb-084c-1758217d06f9
md"## End of clip-04-64-79t.jl"

# ╔═╡ Cell order:
# ╠═5e06a11c-04f6-11eb-0116-1fd51d614f99
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
# ╟─82407c7c-04fd-11eb-2bba-3783e79b0671
# ╠═da70abc2-0527-11eb-31fa-d5b029193ba5
# ╠═82483e30-04fd-11eb-1290-b98ae5d06a45
# ╠═8250f05c-04fd-11eb-0370-899604fffdfd
# ╟─825e8410-04fd-11eb-3e36-6d3bac027215
# ╠═423cd608-0529-11eb-126e-7fde3de2165b
# ╟─82a6656e-04fd-11eb-08c9-3bd9d2d8af20
# ╠═8ad73da4-0529-11eb-3e51-910e74043bfe
# ╠═af9e141e-0529-11eb-3c28-b986f5939349
# ╟─82ba1b72-04fd-11eb-33ed-0db4307c66fc
# ╠═d9b4a9fa-0529-11eb-21b5-2da845b95fe9
# ╠═b0b1e776-0575-11eb-093b-f1cb2d9c4b2e
# ╟─832a285e-04fd-11eb-061b-49a40155abe1
# ╠═cdb02de2-0575-11eb-121e-df6ecd57ba7c
# ╟─8343fe0a-04fd-11eb-2a3e-79daed2e6b98
# ╠═1233198e-0576-11eb-3c17-672034589062
# ╠═833793e0-04fd-11eb-0706-7567eada0aed
# ╠═83b54bb4-04fd-11eb-3208-4162cb3eabc7
# ╠═d9fbed6e-0572-11eb-29bd-1d30ae286045
# ╟─83f2d7f4-04fd-11eb-3114-d3faa79abb1a
# ╠═21d2a6a0-0573-11eb-1043-adeebe1e8fde
# ╠═840196ae-04fd-11eb-1c8a-3511b18e2835
# ╟─4f386d5c-0573-11eb-084c-1758217d06f9
