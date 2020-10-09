### A Pluto.jl notebook ###
# v0.12.3

using Markdown
using InteractiveUtils

# ╔═╡ 846cc37e-0a35-11eb-24ed-7f513230036a
using Pkg, DrWatson

# ╔═╡ 846cf9b6-0a35-11eb-3d6a-91ca6fd4e0d2
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ 8067a2c6-0a35-11eb-14cf-074eb476bf71
md"## Clip-05-45-79.jl"

# ╔═╡ 846d894e-0a35-11eb-365c-b95a0dc20b9f
md"### snippet 5.45"

# ╔═╡ 847e685e-0a35-11eb-3a5d-91ac3078947e
df = CSV.read(sr_datadir("Howell1.csv"), DataFrame);

# ╔═╡ 847ee568-0a35-11eb-1a4f-79f62939c240
md"### snippet 5.46"

# ╔═╡ 848e8644-0a35-11eb-2963-0d76c7c22c0b
begin
	μ_female = rand(Normal(178, 20), 10_000)
	μ_male = rand(Normal(178, 20), 10_000) .+ rand(Normal(0, 10), 10_000)
	Text(precis(DataFrame(; μ_female, μ_male); io=String))
end

# ╔═╡ 8496ffe0-0a35-11eb-2c3e-4159fcafa062
md"### snippet 5.47"

# ╔═╡ 849b6e0e-0a35-11eb-1451-279722a36141
df.sex = ifelse.(df.male .== 1, 2, 1);

# ╔═╡ 5da0339a-0a36-11eb-22d1-d9b7aea3644d
md"### snippet 5.48"

# ╔═╡ 6e0203fa-0a36-11eb-0eaa-357f2828ee4f
@model function m5_8(sex, height)
    σ ~ Uniform(0, 50)
    a ~ filldist(Normal(178, 20), 2)
    μ = a[sex]
    height ~ MvNormal(μ, σ)
end

# ╔═╡ 6e0238ac-0a36-11eb-2aad-316bea996d20
begin
	m5_8t = m5_8(df.sex, df.height)
	quap5_8t = quap(m5_8t, NelderMead())
	Text(precis(quap5_8t))
end

# ╔═╡ 6e02e748-0a36-11eb-313f-19358a19a2e1
md"### snippet 5.49"

# ╔═╡ 6e14cf6c-0a36-11eb-1b7c-61ee035c7e1d
begin
	post5_8t = DataFrame(rand(quap5_8t.distr, 1000)', quap5_8t.params)
	post5_8t.diff_fm = post5_8t[:, "a[1]"] .- post5_8t[:, "a[2]"]
	Text(precis(post5_8t; io=String))
end

# ╔═╡ c079b54e-09a8-11eb-1692-7b02d73f6cba
# %% 5.50
d = CSV.read(datadir("exp_raw/milk.csv"), DataFrame; missingstring = "NA")
d.clade |> unique |> sort

# %% 5.51
# You can't just turn Strings into Integers in Julia but hashing them should give the same
# result
d.clade_id = Int.(indexin(d.clade, unique(d.clade)))

# %% 5.52
d.K = zscore(d.kcal_per_g)

@model function clade(clade_id, K)
    σ ~ Exponential(1)
    α ~ filldist(Normal(0, 0.5), length(unique(clade_id)))
    μ = α[clade_id]
    K ~ MvNormal(μ, σ)
end

q5_9 = quap(clade(d.clade_id, d.K))
post = DataFrame(rand(q5_9.distr, 1000)', q5_9.params)
precis(post)

# TODO: plot

# %% 5.54
d.house = rand(1:4, nrow(d))

@model function clade_house(clade_id, house, K)
    σ ~ Exponential(1)
    a ~ filldist(Normal(0, 0.5), length(unique(clade_id)))
    h ~ filldist(Normal(0, 0.5), length(unique(house)))
    μ = a[clade_id] + h[house]
    K ~ MvNormal(μ, σ)
end

q5_10 = quap(clade_house(d.clade_id, d.house, d.K))
DataFrame(rand(q5_10.distr, 1000)', q5_10.params) |> precis

# End of clip-05-45-54.jl


# ╔═╡ Cell order:
# ╟─8067a2c6-0a35-11eb-14cf-074eb476bf71
# ╠═846cc37e-0a35-11eb-24ed-7f513230036a
# ╠═846cf9b6-0a35-11eb-3d6a-91ca6fd4e0d2
# ╟─846d894e-0a35-11eb-365c-b95a0dc20b9f
# ╠═847e685e-0a35-11eb-3a5d-91ac3078947e
# ╟─847ee568-0a35-11eb-1a4f-79f62939c240
# ╠═848e8644-0a35-11eb-2963-0d76c7c22c0b
# ╟─8496ffe0-0a35-11eb-2c3e-4159fcafa062
# ╠═849b6e0e-0a35-11eb-1451-279722a36141
# ╟─5da0339a-0a36-11eb-22d1-d9b7aea3644d
# ╠═6e0203fa-0a36-11eb-0eaa-357f2828ee4f
# ╠═6e0238ac-0a36-11eb-2aad-316bea996d20
# ╟─6e02e748-0a36-11eb-313f-19358a19a2e1
# ╠═6e14cf6c-0a36-11eb-1b7c-61ee035c7e1d
# ╠═c079b54e-09a8-11eb-1692-7b02d73f6cba
