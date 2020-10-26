using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing 
using StatisticalRethinking
Turing.turnprogress(false);

N = 200

df = DataFrame()
df.weight = rand(Uniform(20, 90), N)

f(x) = mean(df.weight) + 1.6x + rand(Normal(0, 10))
df.height = f.(df.weight)
scatter(df.weight, df.height, lab="observations")

@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* weights
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

m4_3t = m4_3(df.weight, df.height)
chns4_3t = sample(m4_3t, NUTS(100, 0.65), N);
chns4_3t |> display

g(x) = mean(chns4_3t[:a]) + mean(chns4_3t[:b]) * x
x = range(minimum(df.weight), stop=maximum(df.weight), length=N)
plot!(x, g.(x), lab="mean", leg=:bottomright)

K = Int(N/2)
pred4_3t = DataFrame()
wvec = [20, 40, 50, 60, 80]
pred4_3t.weight = repeat(wvec, inner=K)

x_test = pred4_3t.weight

m_test = m4_3(x_test, Vector{Union{Missing, Float64}}(undef, length(x_test)))
#m_test - m4_3(x_test, missing)

ranges = UnitRange{Int64}[]
for (i,j) in enumerate(1:length(wvec))
    append!(ranges, [((i-1)*K + 1):(j * K)])
end

predictions = Turing.Inference.predict(m_test, chns4_3t)
pred4_3t.height = [predictions[1, name, 1] for name in names(predictions)]
height2 = 
    [chns4_3t[:a] .+ chns4_3t[:b] .* x_test[ranges[i]]  for i in 1:length(wvec)]
summary_df = DataFrame()
summary_df.weight = [mean(pred4_3t.weight[ranges[i]]) for i in 1:length(wvec)]
summary_df.mu_height = [mean(pred4_3t.height[ranges[i]]) for i in 1:length(wvec)]
summary_df.sd_height = [std(pred4_3t.height[ranges[i]]) for i in 1:length(wvec)]
summary_df.mu_height2 = [mean(height2[i]) for i in 1:length(wvec)]
summary_df.sd_height2 = [std(height2[i]) for i in 1:length(wvec)]
summary_df |> display

scatter!(summary_df.weight, summary_df.mu_height2,
    markersize=10,lab="Predictions 2")

gui()
