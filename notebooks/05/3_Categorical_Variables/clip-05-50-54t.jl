### A Pluto.jl notebook ###
# v0.12.3

using Markdown
using InteractiveUtils

# ╔═╡ 504df9a8-0aff-11eb-13b1-6d8a03e123aa
using Pkg, DrWatson

# ╔═╡ 515f176e-0aff-11eb-342c-9f6e5c550932
begin
	@quickactivate "SR2TuringPluto"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ 425b24f6-0aff-11eb-124e-8f7e7dc68c84
md"## Clip-05-50-54t.jl"

# ╔═╡ 2998907c-0aff-11eb-279d-67723fa39fb3
md"### snippet 5.50"

# ╔═╡ a86c956c-0b00-11eb-0e54-a7fa8640895e
begin
	df = CSV.read(sr_datadir("milk.csv"), DataFrame; missingstring = "NA")
	df.clade |> unique |> sort
end;

# ╔═╡ a86d54ca-0b00-11eb-217f-798f08670430
md"### snippet 5.51"

# ╔═╡ a87df190-0b00-11eb-3747-c35abd52dda3
md"##### You can't just turn Strings into Integers in Julia but hashing them should give the same result"

# ╔═╡ a87e7c9e-0b00-11eb-12f2-53d3c65e745c
df.clade_id = Int.(indexin(df.clade, unique(df.clade)))

# ╔═╡ a88a6312-0b00-11eb-1238-0d168fc40850
md"### snippet 5.52"

# ╔═╡ a88af58e-0b00-11eb-3729-4d96fd03ce11
df.K = zscore(df.kcal_per_g)

# ╔═╡ a896bdae-0b00-11eb-09c8-4185f6e08bea
md"##### Clade DynamicPPL model"

# ╔═╡ a8989c8c-0b00-11eb-2964-8bb8b153bf44
@model function ppl5_9(clade_id, K)
    σ ~ Exponential(1)
    α ~ filldist(Normal(0, 0.5), length(unique(clade_id)))
    μ = α[clade_id]
    K ~ MvNormal(μ, σ)
end

# ╔═╡ a8a04290-0b00-11eb-3ea7-afbefe340644
begin
	m5_9t = ppl5_9(df.clade_id, df.K)
	quap5_9t = quap(m5_9t)
	post5_9t = DataFrame(rand(quap5_9t.distr, 1000)', quap5_9t.params)
	Text(precis(post5_9t; io=String))
end

# ╔═╡ a8a809bc-0b00-11eb-3d71-9fdca69bcc10
# TODO: plot

md"### snippet 5.54"

# ╔═╡ a8b74526-0b00-11eb-2e7e-43231317823c
df.house = rand(1:4, nrow(df));

# ╔═╡ a8ba5926-0b00-11eb-17d3-ed90f3311578
@model function ppl5_10(clade_id, house, K)
    σ ~ Exponential(1)
    a ~ filldist(Normal(0, 0.5), length(unique(clade_id)))
    h ~ filldist(Normal(0, 0.5), length(unique(house)))
    μ = a[clade_id] + h[house]
    K ~ MvNormal(μ, σ)
end

# ╔═╡ a8c94a96-0b00-11eb-0290-ef5f4f51848e
begin
	m5_10t = ppl5_10(df.clade_id, df.house, df.K)
	quap5_10t = quap(m5_10t)
	Text(precis(DataFrame(rand(quap5_10t.distr, 1000)', quap5_10t.params); io=String))
end

# ╔═╡ a8cc6bb8-0b00-11eb-1f59-099af38de865
md"## End of clip-05-45-54.jl"

# ╔═╡ Cell order:
# ╟─425b24f6-0aff-11eb-124e-8f7e7dc68c84
# ╠═504df9a8-0aff-11eb-13b1-6d8a03e123aa
# ╠═515f176e-0aff-11eb-342c-9f6e5c550932
# ╟─2998907c-0aff-11eb-279d-67723fa39fb3
# ╠═a86c956c-0b00-11eb-0e54-a7fa8640895e
# ╟─a86d54ca-0b00-11eb-217f-798f08670430
# ╟─a87df190-0b00-11eb-3747-c35abd52dda3
# ╠═a87e7c9e-0b00-11eb-12f2-53d3c65e745c
# ╟─a88a6312-0b00-11eb-1238-0d168fc40850
# ╠═a88af58e-0b00-11eb-3729-4d96fd03ce11
# ╟─a896bdae-0b00-11eb-09c8-4185f6e08bea
# ╠═a8989c8c-0b00-11eb-2964-8bb8b153bf44
# ╠═a8a04290-0b00-11eb-3ea7-afbefe340644
# ╠═a8a809bc-0b00-11eb-3d71-9fdca69bcc10
# ╠═a8b74526-0b00-11eb-2e7e-43231317823c
# ╠═a8ba5926-0b00-11eb-17d3-ed90f3311578
# ╠═a8c94a96-0b00-11eb-0290-ef5f4f51848e
# ╟─a8cc6bb8-0b00-11eb-1f59-099af38de865
