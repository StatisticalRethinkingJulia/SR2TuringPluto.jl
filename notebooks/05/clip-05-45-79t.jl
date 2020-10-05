# Clip-05-45-79.jl

using DrWatson
@quickactivate "StatReth"
using Turing

include(srcdir("quap.jl"))

# %% 5.45
d = CSV.read(datadir("exp_raw/Howell_1.csv"), DataFrame)

# %% 5.46
μ_female = rand(Normal(178, 20), 10_000)
μ_male = rand(Normal(178, 20), 10_000) .+ rand(Normal(0, 10), 10_000)
DataFrame((; μ_female, μ_male)) |> precis

# %% 5.47
d.sex = ifelse.(d.male .== 1, 2, 1)

# %% 5.48
@model function categ(sex, height)
    σ ~ Uniform(0, 50)
    a ~ filldist(Normal(178, 20), 2)
    μ = a[sex]
    height ~ MvNormal(μ, σ)
end

q5_8 = quap(categ(d.sex, d.height), NelderMead())
# precis

# %% 5.49
post = DataFrame(rand(q5_8.distr, 1000)', q5_8.params)
post.diff_fm = post[:, "a[1]"] .- post[:, "a[2]"]
precis(post)

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
