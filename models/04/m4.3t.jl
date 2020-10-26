using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
	Turing.turnprogress(false)
end

md"## Model m4.3t"

begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)
end;

Text(precis(df; io=String))

@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* (weights .- x̄)
    for i in eachindex(heights)
        heights[i] ~ Normal(μ[i], σ)
    end
end

m4_3t = m4_3(df.weight, df.height)

q4_3t = quap(m4_3t, NelderMead())

round.(q4_3t.vcov, digits = 3)

begin
	scatter(df.weight, df.height, leg=false)
	quap4_3t = DataFrame(rand(q4_3t.distr, 10_000)', q4_3t.params)
	a_map = mean(quap4_3t.a)
	b_map = mean(quap4_3t.b)
	plot!(x, a_map .+ b_map .* (x .- x̄))
end

Text(precis(quap4_3t; io=String))

Text(precis(m4_3t; io=String))

begin
	chns4_3t = sample(m4_3t, NUTS(), 1000)
	Particles(chns4_3t[[:a, :b, :σ]])
end

f(x) = mean(chns4_3t[:a]) +  mean(chns4_3t[:b]) * x  + 3 * randn();

plot(x, f.(x))

begin
	Δ = 5.0
	x_test = [45.0, 50.0, 55.0, 60.0]
	y_test = f.(x_test)
	m_train = m4_3(x_test, y_test)
end

m_test = m4_3(x_test, missing)

predictions = Turing.Inference.predict(m_test, chns4_3t)

md"## End of m4.3t.jl"
