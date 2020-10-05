
using DrWatson
@quickactivate "StatReth"
using Turing

include(srcdir("quap.jl"))

d = CSV.read(datadir("exp_raw/Howell_1.csv"), DataFrame)

μ_female = rand(Normal(178, 20), 10_000)
μ_male = rand(Normal(178, 20), 10_000) .+ rand(Normal(0, 10), 10_000)
DataFrame((; μ_female, μ_male)) |> precis

d.sex = ifelse.(d.male .== 1, 2, 1)

@model function categ(sex, height)
    σ ~ Uniform(0, 50)
    a ~ filldist(Normal(178, 20), 2)
    μ = a[sex]
    height ~ MvNormal(μ, σ)
end

q5_8 = quap(categ(d.sex, d.height), NelderMead())

post = DataFrame(rand(q5_8.distr, 1000)', q5_8.params)
post.diff_fm = post[:, "a[1]"] .- post[:, "a[2]"]
precis(post)

d = CSV.read(datadir("exp_raw/milk.csv"), DataFrame; missingstring = "NA")
d.clade |> unique |> sort

d.clade_id = Int.(indexin(d.clade, unique(d.clade)))

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

