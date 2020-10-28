
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing 
	using StatisticalRethinking
	Turing.turnprogress(false);
end

md"## Predict.jl"

begin
	N = 200
	df = DataFrame()
	df.weight = rand(Uniform(20, 90), N)
	f(x) = mean(df.weight) + 1.6x + rand(Normal(0, 10))
	df.height = f.(df.weight)
	Text(precis(df; io=String))
end

@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* weights
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

begin
	m4_3t = m4_3(df.weight, df.height)
	chns4_3t = sample(m4_3t, NUTS(100, 0.65), N);
end

md"##### Use DataFrame for easy handling."

begin
	chns4_3t_df = DataFrame(Array(chns4_3t, :parameters), names(chns4_3t, :parameters))
	Text(precis(chns4_3t_df; io=String))
end

Particles(chns4_3t[[:a, :b, :σ]])

Particles(chns4_3t_df)

md"##### Prediction part."

begin
	x_test = [20, 40, 50, 60, 80]
	m_test = m4_3(x_test, fill(missing, length(x_test)))
	predictions = Turing.Inference.predict(m_test, chns4_3t)
end

md"##### Convert to a DataFrame for easy handling."

begin
	pred_df = DataFrame(predictions)[:, 3:end]
	Text(precis(pred_df; io=String))
end

md"#### Get the mean predicted values."

begin
	pred_summary_df = DataFrame()
	pred_summary_df.parameters = names(predictions)
	pred_summary_df.mean = mean(predictions)[:, 2]
	pred_summary_df.std = vec(std(Array(group(predictions, :heights)); dims = 1))
	pred_summary_df.lower = hpd(predictions; alpha=0.11)[:, 2]
	pred_summary_df.upper = hpd(predictions; alpha=0.11)[:, 3]
	pred_summary_df
end

begin
	scatter(df.weight, df.height, lab="Observations",leg=:bottomright)
	scatter!(x_test, pred_summary_df.mean,
    markersize=5, color=:red, lab="Predictions at x_test")
	x = range(minimum(df.weight), stop=maximum(df.weight), length=N)
	g(x) = mean(chns4_3t[:a]) + mean(chns4_3t[:b]) * x
	plot!(x, g.(x), lab="Regression line", leg=:bottomright)
end

md"## End of predict.jl"

