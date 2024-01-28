### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 3f51597b-b20c-4927-afb0-a39f795ad838
using Pkg, DrWatson

# ╔═╡ 6ce23b09-8255-491d-8a61-c6076f1be17d
begin
	using Distributed
	using GLM
	using CSV
	using Random
	using StatsBase
	using DataFrames
	using Optim
	using Dagitty
	using Turing
	using StatsPlots
	using StatisticalRethinking
	using StatisticalRethinkingPlots
	using Logging
	using StatsFuns
end

# ╔═╡ 555738fb-d2e2-44a0-bf4a-c660887352d7
begin
	default(labels=false)
	Logging.disable_logging(Logging.Warn);
end;

# ╔═╡ 327b639b-7627-451f-acd0-b3337131e9a9
begin
	sppnames = ["afarensis", "africanus", "habilis", "boisei",
		"rudolfensis", "ergaster", "sapiens"]
	brainvolcc = [438, 452, 612, 521, 752, 871, 1350]
	masskg = [37.0, 35.5, 34.5, 41.5, 55.5, 61.0, 53.5]
	d = DataFrame(:species => sppnames, :brain => brainvolcc, :mass => masskg);

	d[!,:mass_std] = (d.mass .- mean(d.mass))./std(d.mass)
	d[!,:brain_std] = d.brain ./ maximum(d.brain)
	d
end

# ╔═╡ cbc761fd-05fd-4cd6-ae0b-fa4e8489acaf
@model function model_m7_1(mass_std, brain_std)
    a ~ Normal(0.5, 1)
    b ~ Normal(0, 10)
    μ = @. a + b*mass_std
    log_σ ~ Normal(0, 1)
    brain_std ~ MvNormal(μ, exp(log_σ))
    return(μ)
end

# ╔═╡ 3eeeceba-9016-4aae-a9e5-8505bb898078
begin
	m7_1_ch = sample(model_m7_1(d.mass_std, d.brain_std), NUTS(),
		, 1000)
	m7_1 = DataFrame(m7_1_ch)
	describe(m7_1)
end

# ╔═╡ 7c55e830-8f3e-4840-84da-1226bb1972a8
map7_1_estimate = optimize(model_m7_1(d.mass_std, d.brain_std), MAP())

# ╔═╡ 6f8bf6b9-ad7a-44cc-9922-1d8f30b53b8a
begin
	map7_1_ch = sample(model_m7_1(d.mass_std, d.brain_std), Prior(),
		1_000, init_theta = map7_1_estimate.values.array)
	map7_1 = DataFrame(m7_1_ch)
	precis(map7_1)
end

# ╔═╡ eda8e448-c699-4fe1-9d21-2bf1b19dbe70
plot(m7_1_ch)

# ╔═╡ Cell order:
# ╠═3f51597b-b20c-4927-afb0-a39f795ad838
# ╠═6ce23b09-8255-491d-8a61-c6076f1be17d
# ╠═555738fb-d2e2-44a0-bf4a-c660887352d7
# ╠═327b639b-7627-451f-acd0-b3337131e9a9
# ╠═cbc761fd-05fd-4cd6-ae0b-fa4e8489acaf
# ╠═3eeeceba-9016-4aae-a9e5-8505bb898078
# ╠═7c55e830-8f3e-4840-84da-1226bb1972a8
# ╠═6f8bf6b9-ad7a-44cc-9922-1d8f30b53b8a
# ╠═eda8e448-c699-4fe1-9d21-2bf1b19dbe70
