### A Pluto.jl notebook ###
# v0.19.3

html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

using Markdown
using InteractiveUtils

# ╔═╡ 6bd4ea23-0113-4c32-9120-85aedf3cbae9
using Pkg, DrWatson

# ╔═╡ ee02bdae-99a6-417b-8dad-e9647ce30ede
begin
	using Distributions
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

# ╔═╡ 37c493f6-21b6-49bb-9c8e-5e776a298895
md"#### Setting default attributes for plots."

# ╔═╡ f152cab6-6054-4ec1-bdd3-524398f08789
begin
	default(label=false)
	Logging.disable_logging(Logging.Warn);
end

# ╔═╡ d65bc366-66f7-43e2-8e22-6e60a113d71b
md"## 4.4 Linear predictions."

# ╔═╡ 677f355b-b62b-4de3-ab85-29f76f695418
md"### Code 4.37"

# ╔═╡ 7cf4d33a-03d7-4026-9a43-fdbeadc9207f
begin
	howell1_1 = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	howell1_2 = howell1_1[howell1_1.age .>= 18,:];
end

# ╔═╡ 96d82ebd-0a20-48a9-aee0-a249846d03ed
md"#### Fancy way of doing `scatter(howell1_2.weight, howell1_2.height)`."

# ╔═╡ ca806716-1aba-47cc-bf38-e28cd162a389
@df howell1_2 scatter(:weight, :height)

# ╔═╡ d3d34239-5dac-44e5-a054-c5dfc7d81e4f
md"### Code 4.38"

# ╔═╡ 966d9cea-b4e0-4a8e-a6ea-cc5658694b55
let
	N = 100
	global a = rand(Normal(178, 20), N)
	global b = rand(Normal(0, 10), N);
end

# ╔═╡ c3175f2f-e75a-4f3e-a352-6c4bbdc7ed45
md"### Code 4.39"

# ╔═╡ 3b0d7fe2-e606-440f-95be-9a27d32a168f
let
	p = hline([0, 272]; ylims=(-100, 400), xlabel="weight", ylabel="height")
	title!(L"\beta \sim \mathcal{N}(\mu=0,\sigma=10)")

	x_mean = mean(howell1_2.weight)
	xlims = extrema(howell1_2.weight)  # getting min and max in one pass

	for (α, β) ∈ zip(a, b)
		plot!(x -> α + β * (x - x_mean); xlims=xlims, c=:black, alpha=0.3)
	end
	p
end

# ╔═╡ 8d18813b-8d61-4e62-8891-270116773439
md"### Code 4.40"

# ╔═╡ 1ffaa331-c317-4e71-ac50-d7fecc4ca3f5
let
	b = rand(LogNormal(0, 1), 10_000)
	density(b, xlims=(0, 5), bandwidth=0.1)
end

# ╔═╡ 9406a78d-d2d9-4459-9ad7-49acb6c205c3
md"### Code 4.41"

# ╔═╡ 9f7e22c5-8c3b-4244-a087-971e1baadcf7
let
	Random.seed!(2971)
	N = 100
	global a1 = rand(Normal(178, 20), N)
	global b1 = rand(LogNormal(0, 1), N);
end

# ╔═╡ e31fe2c2-626f-4d16-a428-64b31a2859b7
md"### Code 4.42"

# ╔═╡ c61a2ddc-9a56-4b1e-b248-b4b0e47bcce1
xbar = mean(howell1_2.weight)

# ╔═╡ 07d34cef-f7c7-4f3c-8cdd-014093ce3dc2
begin
	@model function height_regr_model(weight, height)
		a ~ Normal(178, 20)
		b ~ LogNormal(0, 1)
		μ = @. a + b * (weight - xbar)
		σ ~ Uniform(0, 50)
		height ~ MvNormal(μ, σ)
	end

	m4_3 = sample(height_regr_model(howell1_2.weight, howell1_2.height),
		NUTS(), 1000)
	
	m4_3 = resetrange(m4_3);
end

# ╔═╡ 9ec5ee01-c0a0-4518-861e-757fb73ef5fa
md"### Code 4.43"

# ╔═╡ f1273dda-bb4b-4520-9e34-189ca9c8e1a0
let
	@model function height_regr_model_exp1(weight, height)
		a ~ Normal(178, 20)
		log_b ~ Normal(0, 1)
		μ = @. a + exp(log_b) * (weight - xbar)
		σ ~ Uniform(0, 50)
		height ~ MvNormal(μ, σ)
	end
	
	m4_3b = sample(height_regr_model_exp1(howell1_2.weight, howell1_2.height),
		NUTS(), 1000)
end

# ╔═╡ 89e79d9c-3d41-4ffa-835c-1bbc8e0332f9
md"### Code 4.44"

# ╔═╡ c2183caa-8f80-48ee-a658-e42b9e1ca073
begin
	m4_3_df = DataFrame(m4_3)
	PRECIS(m4_3_df)
end

# ╔═╡ a6aedc54-d5f5-4057-b485-dfd94a3310ff
md"### Code 4.45"

# ╔═╡ 9dc7eec7-c865-43ea-b786-080d8c104e3c
round.(cov(Matrix(m4_3_df)), digits=3)

# ╔═╡ 5fcdf2a6-00f3-4fa3-88af-fb1fa2de2ac9
md"### Code 4.46"

# ╔═╡ b0362507-0bac-4c0f-b1e0-8d59bdd6faa7
let
	p = @df howell1_2 scatter(:weight, :height; alpha=0.3)

	chain = resetrange(m4_3)
	samples = sample(chain, 1000)

	a_map = mean(samples[:a])
	b_map = mean(samples[:b])
	plot!(x -> a_map + b_map*(x-xbar))
end

# ╔═╡ 7b6f2ce5-7df6-40b8-9f31-5c69b173c30d
md"### Code 4.47"

# ╔═╡ f6485d48-4603-4262-aa88-e36de40fb2f8
let
	post = sample(m4_3, 1000)
	post_df = DataFrame(post)
	post_df[1:5,:]
end

# ╔═╡ c1a68767-ea57-450a-95aa-eaef5260a3c8
md"### Code 4.48"

# ╔═╡ 8008a0b2-98fe-4352-8f9b-bd4a5c5fd888
begin
	N = 10
	dN = howell1_2[1:N,:]

	@model function height_regr_model_N(weight, height)
		a ~ Normal(178, 20)
		b ~ LogNormal(0, 1)
		m_weight = mean(weight)
		μ = @. a + b * (weight - m_weight)
		σ ~ Uniform(0, 50)
		height ~ MvNormal(μ, σ)
	end

	mN = sample(height_regr_model_N(dN.weight, dN.height), NUTS(), 1000)
	mN = resetrange(mN);
end

# ╔═╡ a7ae1560-b590-4b48-822c-93a23bdd5d25
md"### Code 4.49"

# ╔═╡ 9c919e61-05af-40ae-961c-7eeffa4e18df
begin
	post = sample(mN, 20)
	post_df = DataFrame(post);

	xlims = extrema(howell1_2.weight)
	ylims = extrema(howell1_2.height)
	p = @df dN scatter(:weight, :height; xlims=xlims, ylims=ylims)
	title!("N = $N"; xlab="weight", ylab="height")

	x_mean = mean(dN.weight)
	for (a, b) ∈ zip(post_df.a, post_df.b)
		plot!(x -> a + b * (x-x_mean); c="black", alpha=0.3)
	end
	p
end

# ╔═╡ beda3a79-ef9a-4312-92e4-fed560470fff
md"### Code 4.50"

# ╔═╡ 60c5c43f-3a3b-40b2-8145-9c93c2735e27
let
	post = sample(m4_3, 1000)
	post_df = DataFrame(post)
	global μ_at_50 = @. post_df.a + post_df.b * (50 - xbar);
end

# ╔═╡ 58f72009-4b72-428d-90cb-4e3251d8fccb
md"### Code 4.51"

# ╔═╡ 75dec612-305d-4326-8744-a782e3eb2191
density(μ_at_50; lw=2, xlab=L"\mu|weight=50")

# ╔═╡ 63f6ad59-47f2-4877-8d17-19210778f033
md"### Code 4.52"

# ╔═╡ c1835019-5507-4e9e-9dc6-b26782a4d61d
PI(μ_at_50, perc_prob=0.89)

# ╔═╡ 1b909d62-15c8-424d-8200-1c09facc569d
md"### Code 4.53"

# ╔═╡ 31533bcb-7c41-4471-bf2d-a73430dd690d
let
	μ = StatisticalRethinking.link(post_df, [:a :b], howell1_2.weight, xbar);
	μ = hcat(μ...);
	Base.size(μ), μ[1:5,1]
end

# ╔═╡ a4252214-f620-4a50-8d94-f8aaec103661
md"### Code 4.54"

# ╔═╡ 2eb9634e-8b14-4ad6-ad41-f8898b946a89
begin
	weight_seq = 25:70
	μ = StatisticalRethinking.link(post_df, [:a :b], weight_seq, xbar);
	μ = hcat(μ...);
	Base.size(μ), μ[1:5,1]
end

# ╔═╡ 5a984e54-88db-4d81-8505-0a960707a341
md"### Code 4.55"

# ╔═╡ 95de9e9e-de83-4a05-993a-f560955af5ca
let
	p = plot()
	for i in 1:size(μ, 1)
		scatter!(weight_seq, μ[i,:]; c=:blue, alpha=0.2)
	end
	p
end

# ╔═╡ 2dc235c9-3b61-48e1-9eb1-af8e68d93944
md"### Code 4.56"

# ╔═╡ 5241b732-abc6-4cfb-8a19-8dd47d072929
let
	μ_mean = mean.(eachcol(μ))
	global μ_PI = PI.(eachcol(μ))
	μ_PI = vcat(μ_PI'...);

	# Code 4.57

	@df howell1_2 scatter(:weight, :height; alpha=0.2, xlab="weight", ylab="height")
	plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.3)
end

# ╔═╡ 5eb9ddc6-1f6b-492c-9286-45b096a0a618
md"### Code 4.58"

# ╔═╡ c86f1d9e-514d-43dd-abe3-b64c52d1c7da
let
	post = sample(m4_3, 1000)
	post = DataFrame(post)

	weight_seq = 25:70
	μ = map(w -> post.a + post.b * (w - xbar), weight_seq)
	μ = hcat(μ...)
	global μ_mean = mean.(eachcol(μ))
	μ_CI = PI.(eachcol(μ));
	# -

	# Code 4.59

	sim_height = simulate(post, [:a, :b, :σ], weight_seq .- xbar);
	Base.size(sim_height), sim_height[1:5,1]

	# Code 4.60

	height_PI = PI.(eachcol(sim_height))
	height_PI = vcat(height_PI'...);

	# Code 4.61

	@df howell1_2 scatter(:weight, :height; alpha=0.2, xlab="weight", ylab="height")
	plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.3)
	plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=height_PI, fillalpha=0.3)

	# Code 4.62

	post = sample(m4_3, 10_000)
	post = DataFrame(post)
	sim_height = simulate(post, [:a, :b, :σ], weight_seq .- xbar)
	height_PI = PI.(eachcol(sim_height))
	height_PI = vcat(height_PI'...);

	@df howell1_2 scatter(:weight, :height; alpha=0.2, xlab="weight", ylab="height")
	plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.3)
	plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=height_PI, fillalpha=0.3)
end

# ╔═╡ 3b6dc9df-940e-43cd-bb05-6f63c9bbb9ad
md"### Code 4.63"

# ╔═╡ a511472d-4efc-4aef-9e6a-e2684b41a899
let
	post = sample(m4_3, 1000)
	post = DataFrame(post)

	sim_height = [
		[
			rand(Normal(a + b * (w - xbar), σ))
			for (a, b, σ) ∈ zip(post.a, post.b, post.σ)
		]
		for w ∈ weight_seq
	]
	sim_height = hcat(sim_height...)

	height_PI = PI.(eachcol(sim_height));
	height_PI = vcat(height_PI'...);
	# -

	@df howell1_2 scatter(:weight, :height; alpha=0.2, xlab="weight", ylab="height")
	plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=μ_PI, fillalpha=0.3)
	plot!(weight_seq, [μ_mean μ_mean]; c=:black, fillrange=height_PI, fillalpha=0.3)
end

# ╔═╡ 1b8d53be-7082-4879-a1c9-d47ec9c49402
md"## 4.5 Curves from lines."

# ╔═╡ 3296c41c-45a8-4a8f-8bb5-f7cfe842e950
md"### Code 4.64"

# ╔═╡ 7a66883d-cc3f-43dd-b34d-349c0d4421bf
scatter(howell1_1.weight, howell1_1.height; alpha=0.3)

# ╔═╡ 9eaf1c56-3c01-4d67-a9ed-1ffc259c0876
md"### Code 4.65"

# ╔═╡ cfd2c877-f181-478c-8b2d-fd46ac1a50d3
let
	howell1_1[!, :weight_s] = standardize(ZScoreTransform, howell1_1.weight)
	howell1_1[!, :weight_s2] = howell1_1.weight_s.^2;

	@model function height_regr_m2(weight_s, weight_s2, height)
		a ~ Normal(178, 20)
		b1 ~ LogNormal(0, 1)
		b2 ~ Normal(0, 1)
		μ = @. a + b1 * weight_s + b2 * weight_s2
		σ ~ Uniform(0, 50)
		height ~ MvNormal(μ, σ)
	end

	m4_5 = sample(height_regr_m2(howell1_1.weight_s, howell1_1.weight_s2,
		howell1_1.height), NUTS(), 1000)
	m4_5 = resetrange(m4_5)
	global m4_5_df = DataFrame(m4_5);
end

# ╔═╡ 9d10ad06-61f9-40dd-9e4b-95791d80ebd9
md"### Code 4.66"

# ╔═╡ 0a384315-7a21-4441-9e46-c23fc0cdd7f2
PRECIS(m4_5_df)

# ╔═╡ c6db5270-a8b6-4e78-a697-6e213d1d194e
md"### Code 4.67"

# ╔═╡ dd42d76b-0fd6-432b-bb59-20545bea9107
let
	df = sample(m4_5_df, 1000)
	weight_seq = range(-2.2, 2; length=30)

	# explicit logic of link
	
	mu = [
		df.a + df.b1 * w_s + df.b2 * w_s^2
		for w_s ∈ weight_seq
	]

	mu = hcat(mu...)
	mu_mean = mean.(eachcol(mu))
	mu_PI = PI.(eachcol(mu))
	mu_PI = vcat(mu_PI'...)

	# explicit logic of sim
	
	sim_height = [
		[
			rand(Normal(row.a + row.b1 * w_s + row.b2 * w_s^2, row.σ))
			for row ∈ eachrow(df)
		]
		for w_s ∈ weight_seq
	]
	sim_height = hcat(sim_height...);

	height_PI = PI.(eachcol(sim_height))
	height_PI = vcat(height_PI'...);

	# Code 4.68

	global p_square = @df howell1_1 scatter(:weight_s, :height; alpha=0.3, title="Square poly")
	plot!(weight_seq, mu_mean; c=:black)
	plot!(weight_seq, [mu_mean mu_mean]; c=:black, fillrange=mu_PI, fillalpha=0.3)
	plot!(weight_seq, [mu_mean mu_mean]; c=:black, fillrange=height_PI,
		fillalpha=0.3)
end

# ╔═╡ 4be78229-d1f6-4b4a-a838-f4f848ce116b
md"### Code 4.69"

# ╔═╡ 925c4772-1445-46dc-9a6d-6a33da9a3eb7
names(howell1_1)

# ╔═╡ 4dcbc6d2-29e1-48b6-9d86-287f64a589a0
begin
	scale!(howell1_1, :weight)
	howell1_1[!, :weight_s3] = howell1_1.weight_s.^3
end

# ╔═╡ a44126dc-9050-41c0-a82b-d1a918c3b5a2
@model function height_regr_m3(weight_s, weight_s2, weight_s3, height)
    a ~ Normal(178, 20)
    b1 ~ LogNormal(0, 1)
    b2 ~ Normal(0, 1)
    b3 ~ Normal(0, 1)
    μ = @. a + b1 * weight_s + b2 * weight_s2 + b3 * weight_s3
    σ ~ Uniform(0, 50)
    height ~ MvNormal(μ, σ)
end

# ╔═╡ a97a5c36-d29c-40d3-b10a-d744a1b71b18
let
	m4_6 = sample(height_regr_m3(howell1_1.weight_s, howell1_1.weight_s2,
			howell1_1.weight_s3, howell1_1.height), NUTS(), 1000)
	m4_6 = resetrange(m4_6)
	global m4_6_df = DataFrame(m4_6)
	PRECIS(m4_6_df)
end

# ╔═╡ 6a23ad50-9e8a-4375-a0b2-e71102dd7730
let
	df = sample(m4_6_df, 1000)
	weight_seq = range(-2.2, 2; length=30)

	# explicit logic of link
	
	mu = [
		df.a + df.b1 * w_s + df.b2 * w_s^2 + df.b3 * w_s^3
		for w_s ∈ weight_seq
	]

	mu = hcat(mu...)
	global mu_mean = mean.(eachcol(mu))
	global mu_PI = PI.(eachcol(mu))
	mu_PI = vcat(mu_PI'...)

	# explicit logic of sim
	
	global sim_height = [
		[
			rand(Normal(row.a + row.b1 * w_s + row.b2 * w_s^2 + row.b3 * w_s^3,
					row.σ)) for row ∈ eachrow(df)
		]
		for w_s ∈ weight_seq
	]
	sim_height = hcat(sim_height...);

	height_PI = PI.(eachcol(sim_height))
	height_PI = vcat(height_PI'...);

	p_cubic = @df howell1_1 scatter(:weight_s, :height; alpha=0.3,
		title="Cubic poly")
	plot!(weight_seq, mu_mean; c=:black)
	plot!(weight_seq, [mu_mean mu_mean]; c=:black, fillrange=mu_PI, fillalpha=0.3)
	plot!(weight_seq, [mu_mean mu_mean]; c=:black, fillrange=height_PI,
		fillalpha=0.3)

	plot(p_square, p_cubic; layout=(1, 2))
end

# ╔═╡ a9c133a1-6db1-4881-8096-b24ebddcdf56
md"### Codes 4.70 and 4.71"

# ╔═╡ 16876c69-3b39-4b31-b9c6-1d08b8a86f49
# Looks like Julia plots don't support change of ticks proposed in the book.
# But much more natural way will be to remap values we're plotting back to
# the original scale.

# Example of this is below.

let
	weight_seq = range(-2.2, 2; length=30)
	μ = mean(howell1_1.weight)
	σ = std(howell1_1.weight)
	w = @. howell1_1.weight_s * σ + μ
	scatter(w, howell1_1.height; alpha=0.3)
	height_PI = PI.(eachcol(sim_height))
	height_PI = vcat(height_PI'...);

	w_s = @. weight_seq * σ + μ
	plot!(w_s, mu_mean; c=:black)
	plot!(w_s, [mu_mean mu_mean]; c=:black, fillrange=mu_PI, fillalpha=0.3)
	plot!(w_s, [mu_mean mu_mean]; c=:black, fillrange=height_PI, fillalpha=0.3)
end

# ╔═╡ Cell order:
# ╠═6bd4ea23-0113-4c32-9120-85aedf3cbae9
# ╠═ee02bdae-99a6-417b-8dad-e9647ce30ede
# ╠═37c493f6-21b6-49bb-9c8e-5e776a298895
# ╠═f152cab6-6054-4ec1-bdd3-524398f08789
# ╟─d65bc366-66f7-43e2-8e22-6e60a113d71b
# ╟─677f355b-b62b-4de3-ab85-29f76f695418
# ╠═7cf4d33a-03d7-4026-9a43-fdbeadc9207f
# ╟─96d82ebd-0a20-48a9-aee0-a249846d03ed
# ╠═ca806716-1aba-47cc-bf38-e28cd162a389
# ╟─d3d34239-5dac-44e5-a054-c5dfc7d81e4f
# ╠═966d9cea-b4e0-4a8e-a6ea-cc5658694b55
# ╟─c3175f2f-e75a-4f3e-a352-6c4bbdc7ed45
# ╠═3b0d7fe2-e606-440f-95be-9a27d32a168f
# ╟─8d18813b-8d61-4e62-8891-270116773439
# ╠═1ffaa331-c317-4e71-ac50-d7fecc4ca3f5
# ╟─9406a78d-d2d9-4459-9ad7-49acb6c205c3
# ╠═9f7e22c5-8c3b-4244-a087-971e1baadcf7
# ╟─e31fe2c2-626f-4d16-a428-64b31a2859b7
# ╠═c61a2ddc-9a56-4b1e-b248-b4b0e47bcce1
# ╠═07d34cef-f7c7-4f3c-8cdd-014093ce3dc2
# ╟─9ec5ee01-c0a0-4518-861e-757fb73ef5fa
# ╠═f1273dda-bb4b-4520-9e34-189ca9c8e1a0
# ╟─89e79d9c-3d41-4ffa-835c-1bbc8e0332f9
# ╠═c2183caa-8f80-48ee-a658-e42b9e1ca073
# ╟─a6aedc54-d5f5-4057-b485-dfd94a3310ff
# ╠═9dc7eec7-c865-43ea-b786-080d8c104e3c
# ╟─5fcdf2a6-00f3-4fa3-88af-fb1fa2de2ac9
# ╠═b0362507-0bac-4c0f-b1e0-8d59bdd6faa7
# ╟─7b6f2ce5-7df6-40b8-9f31-5c69b173c30d
# ╠═f6485d48-4603-4262-aa88-e36de40fb2f8
# ╟─c1a68767-ea57-450a-95aa-eaef5260a3c8
# ╠═8008a0b2-98fe-4352-8f9b-bd4a5c5fd888
# ╟─a7ae1560-b590-4b48-822c-93a23bdd5d25
# ╠═9c919e61-05af-40ae-961c-7eeffa4e18df
# ╟─beda3a79-ef9a-4312-92e4-fed560470fff
# ╠═60c5c43f-3a3b-40b2-8145-9c93c2735e27
# ╠═58f72009-4b72-428d-90cb-4e3251d8fccb
# ╠═75dec612-305d-4326-8744-a782e3eb2191
# ╟─63f6ad59-47f2-4877-8d17-19210778f033
# ╠═c1835019-5507-4e9e-9dc6-b26782a4d61d
# ╟─1b909d62-15c8-424d-8200-1c09facc569d
# ╠═31533bcb-7c41-4471-bf2d-a73430dd690d
# ╟─a4252214-f620-4a50-8d94-f8aaec103661
# ╠═2eb9634e-8b14-4ad6-ad41-f8898b946a89
# ╟─5a984e54-88db-4d81-8505-0a960707a341
# ╠═95de9e9e-de83-4a05-993a-f560955af5ca
# ╟─2dc235c9-3b61-48e1-9eb1-af8e68d93944
# ╠═5241b732-abc6-4cfb-8a19-8dd47d072929
# ╟─5eb9ddc6-1f6b-492c-9286-45b096a0a618
# ╠═c86f1d9e-514d-43dd-abe3-b64c52d1c7da
# ╟─3b6dc9df-940e-43cd-bb05-6f63c9bbb9ad
# ╠═a511472d-4efc-4aef-9e6a-e2684b41a899
# ╟─1b8d53be-7082-4879-a1c9-d47ec9c49402
# ╠═3296c41c-45a8-4a8f-8bb5-f7cfe842e950
# ╠═7a66883d-cc3f-43dd-b34d-349c0d4421bf
# ╟─9eaf1c56-3c01-4d67-a9ed-1ffc259c0876
# ╠═cfd2c877-f181-478c-8b2d-fd46ac1a50d3
# ╟─9d10ad06-61f9-40dd-9e4b-95791d80ebd9
# ╠═0a384315-7a21-4441-9e46-c23fc0cdd7f2
# ╟─c6db5270-a8b6-4e78-a697-6e213d1d194e
# ╠═dd42d76b-0fd6-432b-bb59-20545bea9107
# ╟─4be78229-d1f6-4b4a-a838-f4f848ce116b
# ╠═925c4772-1445-46dc-9a6d-6a33da9a3eb7
# ╠═4dcbc6d2-29e1-48b6-9d86-287f64a589a0
# ╠═a44126dc-9050-41c0-a82b-d1a918c3b5a2
# ╠═a97a5c36-d29c-40d3-b10a-d744a1b71b18
# ╠═6a23ad50-9e8a-4375-a0b2-e71102dd7730
# ╟─a9c133a1-6db1-4881-8096-b24ebddcdf56
# ╠═16876c69-3b39-4b31-b9c6-1d08b8a86f49
