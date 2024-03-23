### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 8483abbb-011f-4e36-976b-d04967d86af0
using Pkg, DrWatson

# ╔═╡ 772c0ef3-aeb2-4c8f-aac9-fbd189770956
begin
	using Distributions
	using Optim
	using CSV
	using DataFrames
	using LinearAlgebra
	using Logging
	using Random
	using StatsBase
	using Turing
	using StatisticalRethinking
	using StatisticalRethinkingPlots
	using StatsPlots, LaTeXStrings
	using BSplines
end

# ╔═╡ 9cda4135-14bf-4a75-bce1-e0317c682d40
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 5%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ 74a4750e-98a8-48af-a789-2109e2c5e01f
md"#### Setting default attributes for plots."

# ╔═╡ dbad5708-44ff-488e-a6ef-654c9252f940
begin
	default(label=false)
	Logging.disable_logging(Logging.Warn);
end;

# ╔═╡ 13285996-780c-4079-9991-8010e46ad582
md"### Code 4.72"

# ╔═╡ 5bb57e92-1f89-40d7-b025-8294ef3bec09
begin
	d = CSV.read(sr_datadir("cherry_blossoms.csv"), missingstring="NA", DataFrame)
	describe(d)
end

# ╔═╡ bf352765-ed56-436a-8dca-333435ec78eb
md"### Code 4.73"

# ╔═╡ 755df2d9-894a-4ee6-9cd4-32971672f025
begin
	d2 = d[completecases(d[!,[:doy]]),:]
	d2 = disallowmissing(d2[!,[:year,:doy]])
	num_knots = 15
	knots_list = quantile(d2.year, range(0, 1; length=num_knots));
end

# ╔═╡ 881fb044-cdfe-46d4-9843-e28072b85387
md"### Code 4.74"

# ╔═╡ 20f677f0-56dc-4df2-9227-d3be826348e8
basis = BSplineBasis(3, knots_list)

# ╔═╡ d4c60a5d-cd5e-4e54-a496-38b204701af2
md"### Code 4.75"

# ╔═╡ c212dff1-5b2e-48a8-8906-9632fd1fc38f
begin
	p1 = plot(basis)
	scatter!(knots_list, repeat([1], num_knots); xlab="year", ylab="basis",
		legend=false)
end

# ╔═╡ 59024873-f622-4a7d-9763-bda8b4cd85f8
md"### Code 4.76"

# ╔═╡ dcd41a9e-e657-45a5-b0e2-8714d65b8aa8
md"##### This way of calucalting bsplines is slightly slower, than shown in the book (with pre-calculated matrix) but it is much cleaner in my perspective."

# ╔═╡ 35c9eb0e-1251-4359-8139-2b5351d696e4
md"##### You can do comparison yourself by precalculating spline matrix outside of the model and do matrix multiplication in the model instead of spline evialutaion. Example of doing this is at code block 4.79."

# ╔═╡ dc16d5af-6526-4043-81d0-8e6147d4bfdd
@model function model_splines(year, doy)
    w ~ MvNormal(zeros(length(basis)), 1)
    a ~ Normal(100, 10)
    s = Spline(basis, w)
    μ = a .+ s.(year)
    σ ~ Exponential(1)
    doy ~ MvNormal(μ, σ)
end

# ╔═╡ 7836a7a4-15cf-4e4e-9fe3-8577682e037e
m4_7 = sample(model_splines(d2.year, d2.doy), NUTS(0.65; init_ϵ = 9.765625e-5),
	1000)

# ╔═╡ 861b03a9-f0a5-45ab-aceb-ec7c445031ff
md"### Code 4.77"

# ╔═╡ 9b9513f5-7a81-405d-864f-d2fc23a015e5
begin
	post = DataFrame(m4_7)

	# convert columns w[*] into single column w
	
	w_df = DataFrames.select(post, r"w")
	post = DataFrames.select(post, Not(r"w"))
	post[!,:w] = Vector.(eachrow(w_df))

	# vector of 16 average w values
	
	w_mean = mean.(eachcol(w_df))
	p2 = plot(basis .* w_mean)
	scatter!(knots_list, repeat([1], num_knots); xlab="year", ylab="basis × weight")
end

# ╔═╡ e2dcf8e1-309b-401a-9512-432e599fa438
md"### Code 4.78"

# ╔═╡ 780ba069-4e8b-4832-a0e0-e5cef802fa4a
md"#### Explicit link logic."

# ╔═╡ ec818ff5-79db-4911-b965-3bb7fb429960
begin
	μ = [
		row.a .+ Spline(basis, row.w).(d2.year)
		for row ∈ eachrow(post)
	]
	μ = hcat(μ...);

	μ_PI = PI.(eachrow(μ))
	μ_PI = vcat(μ_PI'...);

	p3 = @df d2 scatter(:year, :doy; alpha=0.3)
	μ_mean = mean.(eachrow(μ_PI))
	plot!(d2.year, [μ_mean, μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.3,
		alpha=0)
end

# ╔═╡ 33bcafa1-c96a-488c-8f6b-6398ebf0588a
plot(p1, p2, p3; layout=(3, 1))

# ╔═╡ ab0f4159-075f-4b05-9188-322a5e3bfe5a
md"### Code 4.79"

# ╔═╡ bf04f6d3-925c-47be-90ff-b8a5c3ccaa05
md"#### How to build the model with explicit spline matrix calculation."

# ╔═╡ 662d9645-50b7-4904-bffb-6939af0e11cf
let
	basis = BSplineBasis(3, knots_list)

	# list of splines with 1 only at corresponding basis index
	splines = [
		Spline(basis, [float(idx == knot) for idx ∈ 1:length(basis)])
		for knot ∈ 1:length(basis)
	]

	# calculate each spline for every year. Resulting matrix B is 827x16
	B = [
		map(s -> s(year), splines)
		for year in d2.year
	]
	B = vcat(B'...);


	# do not need years parameter anymore, all the information is in B matrix
	@model function model_splines_matrix(doy)
		w ~ MvNormal(zeros(length(basis)), 1)
		a ~ Normal(100, 10)
		μ = a .+ B * w
		σ ~ Exponential(1)
		doy ~ MvNormal(μ, σ)
	end

	m4_7alt = sample(model_splines_matrix(d2.doy), NUTS(0.65;
			init_ϵ = 0.0001953125), 1000)
end

# ╔═╡ Cell order:
# ╠═9cda4135-14bf-4a75-bce1-e0317c682d40
# ╠═8483abbb-011f-4e36-976b-d04967d86af0
# ╠═772c0ef3-aeb2-4c8f-aac9-fbd189770956
# ╟─74a4750e-98a8-48af-a789-2109e2c5e01f
# ╠═dbad5708-44ff-488e-a6ef-654c9252f940
# ╟─13285996-780c-4079-9991-8010e46ad582
# ╠═5bb57e92-1f89-40d7-b025-8294ef3bec09
# ╟─bf352765-ed56-436a-8dca-333435ec78eb
# ╠═755df2d9-894a-4ee6-9cd4-32971672f025
# ╟─881fb044-cdfe-46d4-9843-e28072b85387
# ╠═20f677f0-56dc-4df2-9227-d3be826348e8
# ╟─d4c60a5d-cd5e-4e54-a496-38b204701af2
# ╠═c212dff1-5b2e-48a8-8906-9632fd1fc38f
# ╟─59024873-f622-4a7d-9763-bda8b4cd85f8
# ╟─dcd41a9e-e657-45a5-b0e2-8714d65b8aa8
# ╟─35c9eb0e-1251-4359-8139-2b5351d696e4
# ╠═dc16d5af-6526-4043-81d0-8e6147d4bfdd
# ╠═7836a7a4-15cf-4e4e-9fe3-8577682e037e
# ╟─861b03a9-f0a5-45ab-aceb-ec7c445031ff
# ╠═9b9513f5-7a81-405d-864f-d2fc23a015e5
# ╟─e2dcf8e1-309b-401a-9512-432e599fa438
# ╟─780ba069-4e8b-4832-a0e0-e5cef802fa4a
# ╠═ec818ff5-79db-4911-b965-3bb7fb429960
# ╠═33bcafa1-c96a-488c-8f6b-6398ebf0588a
# ╟─ab0f4159-075f-4b05-9188-322a5e3bfe5a
# ╟─bf04f6d3-925c-47be-90ff-b8a5c3ccaa05
# ╠═662d9645-50b7-4904-bffb-6939af0e11cf
