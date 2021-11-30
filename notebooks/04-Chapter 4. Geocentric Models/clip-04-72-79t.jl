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
md"## clip-04-72-79t"

# ╔═╡ 82a6656e-04fd-11eb-08c9-3bd9d2d8af20
md"### snippet 4.72"

# ╔═╡ 8ad73da4-0529-11eb-3e51-910e74043bfe
df = CSV.read(sr_datadir("cherry_blossoms.csv"), DataFrame; missingstring = "NA");

# ╔═╡ af9e141e-0529-11eb-3c28-b986f5939349
scatter(df.year, df.doy, leg=false)

# ╔═╡ 82ba1b72-04fd-11eb-33ed-0db4307c66fc
md"### snippet 4.73"

# ╔═╡ d9b4a9fa-0529-11eb-21b5-2da845b95fe9
begin
	dropmissing!(df, :doy)
	num_knots = 15
	knot_list = quantile(df.year, range(0, 1, length = num_knots))
	basis = BSplineBasis(4, knot_list)
	B = basismatrix(basis, df.year)
end;

# ╔═╡ b0b1e776-0575-11eb-093b-f1cb2d9c4b2e
begin
	plot(legend = false, xlabel = "year", ylabel = "basis value")
	for y in eachcol(B)
		plot!(df.year, y)
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
	quap4_7t = quap(spline_1(df.doy))
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
	scatter(df.year, df.doy, alpha = 0.3)
	plot!(df.year, mu_3.mean, ribbon = (mu_3.mean .- mu_3.lower, mu_3.upper .- mu_3.mean))
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
	quap_alt_4_7t = quap(spline_2(df.doy))
	quap_alt_4_7t.coef
end

# ╔═╡ 4f386d5c-0573-11eb-084c-1758217d06f9
md"## End of clip-04-72-79t.jl"

# ╔═╡ Cell order:
# ╟─5e06a11c-04f6-11eb-0116-1fd51d614f99
# ╠═81cdc1aa-04fd-11eb-3235-bff0e3f52535
# ╠═81cdee8c-04fd-11eb-29f5-85927396b6dc
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
# ╟─83b54bb4-04fd-11eb-3208-4162cb3eabc7
# ╠═d9fbed6e-0572-11eb-29bd-1d30ae286045
# ╟─83f2d7f4-04fd-11eb-3114-d3faa79abb1a
# ╠═21d2a6a0-0573-11eb-1043-adeebe1e8fde
# ╠═840196ae-04fd-11eb-1c8a-3511b18e2835
# ╟─4f386d5c-0573-11eb-084c-1758217d06f9
