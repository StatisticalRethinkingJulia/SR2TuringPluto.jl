using Pkg, DrWatson

begin
  @quickactivate "StatisticalRethinkingTuring"
  using Turing
  using StatisticalRethinking
end

begin
  df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
  df = df[df.age .>= 18, :]
  x̄ = mean(df.weight)
  x = range(minimum(df.weight), maximum(df.weight), length = 100)     # either
end;

@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* (weights .- x̄)
    heights ~ MvNormal(μ, σ)
end

m4_3t = m4_3(df.weight, df.height)

q4_3t = quap(m4_3t, NelderMead())

chns4_3t = sample(m4_3t, NUTS(100, 0.65), 200);

f(x) = mean(chns4_3t[:a]) + mean(chns4_3t[:b]) * x + 0.1 * randn()
x_test = x; y_test = f.(x_test)
m_test = m4_3(x_test, Vector{Union{Missing, Float64}}(undef, length(y_test)));
#predictions = Turing.Inference.predict(m_test, chns4_3t)
