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

# ╔═╡ bd02802f-e9c2-47d3-9e1f-22e9a7f4db79
using Pkg, DrWatson

# ╔═╡ 1eb249a5-7e42-4f65-be60-ccc4b5db9f66
begin
    using Distributions
    using StatsPlots
    using StatsBase
    using LaTeXStrings
    using CSV
    using DataFrames
    using LinearAlgebra
    using Logging
    using Random
    using Turing
    using Dagitty
    using StatisticalRethinking
    using StatisticalRethinkingPlots
end

# ╔═╡ 26e632d5-a797-4fb0-b492-70541a272680
md"### Set defaults for plotting and logging."

# ╔═╡ b29185e1-513e-4152-b9cf-b2d62621b21a
begin
	default(label=false)
	Logging.disable_logging(Logging.Warn);
end

# ╔═╡ a9d8b978-a482-4bf7-a31b-f7f3381db19e
md"## 5.3 Categorical variables."

# ╔═╡ e4d4dc29-e5c8-4666-956c-1a2a936c842d
md"### Code 5.45"

# ╔═╡ 90f4730d-36f0-4284-9390-e1066451286a
begin
	d = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	describe(d)
end

# ╔═╡ 9d8c7435-4323-4833-9bfe-1ee4d25e1aba
md"### Code 5.46"

# ╔═╡ b9170040-46e8-478e-9254-79e8d15cca3e
begin
	cnt = 10_000
	μ_female = rand(Normal(178, 20), cnt)
	μ_male = rand(Normal(178, 20), cnt) + rand(Normal(0, 10), cnt)
	PRECIS(DataFrame(:μ_female => μ_female, :μ_male => μ_male))
end

# ╔═╡ 05bf01b6-62b9-44d7-a599-532b2d60ccda
md"### Code 5.47"

# ╔═╡ 2c5475d6-4623-469c-8e22-a19d2a95f58c
begin
	d[!,:sex] = ifelse.(d.male .== 1, 2, 1)
	describe(d.sex)
end

# ╔═╡ a904a158-aba0-45ad-bb17-1a4678a0f651
md"### Code 5.48"

# ╔═╡ 10fa3e5f-85c2-4f4b-a47b-7b7702d26d4f
@model function model_m5_8(sex, height)
    σ ~ Uniform(0, 50)
    a ~ MvNormal([178, 178], 20)
    height ~ MvNormal(a[sex], σ)
end

# ╔═╡ 9a3630ee-a02c-4edf-bd6d-abe27834ec58
begin
	m5_8 = sample(model_m5_8(d.sex, d.height), NUTS(), 1000)
	m5_8_df = DataFrame(m5_8)
	PRECIS(m5_8_df)
end

# ╔═╡ 665efae8-c778-4cb0-a366-42cd377baaaa
md"### Code 5.49"

# ╔═╡ b3c04688-8ba8-4c38-89e4-94800f87af73
begin
	m5_8_df[!,:diff_fm] = m5_8_df[:,"a[1]"] - m5_8_df[:,"a[2]"]
	PRECIS(m5_8_df)
end

# ╔═╡ 128b10c0-b48c-46cc-90d5-e983e7e221f7
md"### Code 5.50"

# ╔═╡ 5f6050b7-2225-4446-8e1b-5586cc339910
begin
	d2 = CSV.read(sr_datadir("milk.csv"), DataFrame)

	# get rid of dots in column names
	rename!(n -> replace(n, "." => "_"), d2)

	levels(d2.clade)
end

# ╔═╡ 04602151-0f64-4956-8a60-0316dc0b0afc
md"### Code 5.51"

# ╔═╡ bd105f62-fec0-4875-8b96-7d607a01409f
d2[!,:clade_id] = indexin(d2.clade, levels(d2.clade));

# ╔═╡ da2626f9-94bf-4c51-b45e-c346530f3ad4
md"### Code 5.52"

# ╔═╡ a93f2020-596f-4758-ba03-b2542c2da03f
begin
	d2[!,:K] = standardize(ZScoreTransform, d2.kcal_per_g);
	clade_counts = maximum(levels(d2.clade_id))
end

# ╔═╡ 0fa9e67a-dfe8-4ba7-8067-05fe4258664f
@model function model_m5_9(clade_id, K)
    clade_μ = zeros(clade_counts)
    a ~ MvNormal(clade_μ, 0.5)
    σ ~ Exponential(1)
    K ~ MvNormal(a[clade_id], σ)
end

# ╔═╡ edc28ea5-d80a-4d59-8773-c3a410355fab
begin
	m5_9 = sample(model_m5_9(d2.clade_id, d2.K), NUTS(), 1000)
	m5_9_df = DataFrame(m5_9)
	# get rid of square brackets in a params
	rename!(n -> replace(n, r"\[|\]" => ""), m5_9_df)

	pars = [:a1, :a2, :a3, :a4]
	p_names = map(v -> "$(v[1]): $(v[2])", zip(pars, levels(d2.clade)))

	coeftab_plot(m5_9_df; pars=pars, pars_names=p_names, xlab="expected kcal (std)")
end

# ╔═╡ d5fb78fc-0462-470f-b10c-88a07951d3d4
md"### Code 5.53"

# ╔═╡ b131a042-f090-49b0-89e9-53a327bf830d
md"#### It took me a while to find a seed which make Slytherin to stand out. It is just a seed, not the model property!"

# ╔═╡ a60b5835-3cfc-49b6-a9ae-4ec72c0498b9
begin
	Random.seed!(31)
	d2[!,:house] = sample(1:4, nrow(d2));
end

# ╔═╡ 85eacd46-854c-46f8-b435-6c4b6d101d6e
md"### Code 5.54"

# ╔═╡ 26a02cf1-3ca7-4056-b3eb-199d4dad6fb7
house_counts = maximum(levels(d2.house))

# ╔═╡ fc65973a-8857-4ad4-9722-7ee0387390d8
@model function model_m5_10(clade_id, house, K)
    clade_μ = zeros(clade_counts)
    house_μ = zeros(house_counts)
    a ~ MvNormal(clade_μ, 0.5)
    h ~ MvNormal(house_μ, 0.5)
    σ ~ Exponential(1)
    μ = a[clade_id] .+ h[house]
    K ~ MvNormal(μ, σ)
end

# ╔═╡ 19d257f9-bc0b-47e6-8915-c448995a2018
begin
	m5_10 = sample(model_m5_10(d2.clade_id, d2.house, d2.K), NUTS(), 1000)
	m5_10_df = DataFrame(m5_10)
	PRECIS(m5_10_df)
end

# ╔═╡ Cell order:
# ╠═bd02802f-e9c2-47d3-9e1f-22e9a7f4db79
# ╠═1eb249a5-7e42-4f65-be60-ccc4b5db9f66
# ╠═26e632d5-a797-4fb0-b492-70541a272680
# ╠═b29185e1-513e-4152-b9cf-b2d62621b21a
# ╟─a9d8b978-a482-4bf7-a31b-f7f3381db19e
# ╟─e4d4dc29-e5c8-4666-956c-1a2a936c842d
# ╠═90f4730d-36f0-4284-9390-e1066451286a
# ╟─9d8c7435-4323-4833-9bfe-1ee4d25e1aba
# ╠═b9170040-46e8-478e-9254-79e8d15cca3e
# ╟─05bf01b6-62b9-44d7-a599-532b2d60ccda
# ╠═2c5475d6-4623-469c-8e22-a19d2a95f58c
# ╟─a904a158-aba0-45ad-bb17-1a4678a0f651
# ╠═10fa3e5f-85c2-4f4b-a47b-7b7702d26d4f
# ╠═9a3630ee-a02c-4edf-bd6d-abe27834ec58
# ╟─665efae8-c778-4cb0-a366-42cd377baaaa
# ╠═b3c04688-8ba8-4c38-89e4-94800f87af73
# ╟─128b10c0-b48c-46cc-90d5-e983e7e221f7
# ╠═5f6050b7-2225-4446-8e1b-5586cc339910
# ╟─04602151-0f64-4956-8a60-0316dc0b0afc
# ╠═bd105f62-fec0-4875-8b96-7d607a01409f
# ╟─da2626f9-94bf-4c51-b45e-c346530f3ad4
# ╠═a93f2020-596f-4758-ba03-b2542c2da03f
# ╠═0fa9e67a-dfe8-4ba7-8067-05fe4258664f
# ╠═edc28ea5-d80a-4d59-8773-c3a410355fab
# ╟─d5fb78fc-0462-470f-b10c-88a07951d3d4
# ╟─b131a042-f090-49b0-89e9-53a327bf830d
# ╠═a60b5835-3cfc-49b6-a9ae-4ec72c0498b9
# ╟─85eacd46-854c-46f8-b435-6c4b6d101d6e
# ╠═26a02cf1-3ca7-4056-b3eb-199d4dad6fb7
# ╠═fc65973a-8857-4ad4-9722-7ee0387390d8
# ╠═19d257f9-bc0b-47e6-8915-c448995a2018
