### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 22398270-13a7-11eb-04c8-8b0995208358
using Pkg, DrWatson

# ╔═╡ 2239c3fc-13a7-11eb-1a00-4f5870c3b7ab
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
	Turing.turnprogress(false)
end

# ╔═╡ 5a9929f0-13a6-11eb-18fd-6bd484be1583
md"## Model m4.3t"

# ╔═╡ 223a5448-13a7-11eb-33d5-bf1dff954915
begin
	df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
	df = df[df.age .>= 18, :]
	x̄ = mean(df.weight)
	x = range(minimum(df.weight), maximum(df.weight), length = 100)
end;

# ╔═╡ 2247ea56-13a7-11eb-1e8b-0388b6f47c98
Text(precis(df; io=String))

# ╔═╡ 224c0bc0-13a7-11eb-3153-cf29938ecffb
@model function m4_3(weights, heights)
    a ~ Normal(178, 20)
    b ~ LogNormal(0, 1)
    σ ~ Uniform(0, 50)
    μ = a .+ b .* (weights .- x̄)
    heights ~ MvNormal(μ, σ)
end

# ╔═╡ 2253f11e-13a7-11eb-1c93-039f8e3bf927
m4_3t = m4_3(df.weight, df.height)

# ╔═╡ 22548624-13a7-11eb-24d0-f90ba17217d8
q4_3t = quap(m4_3t, NelderMead())

# ╔═╡ 225e5a0a-13a7-11eb-140f-bde0811dd1e1
round.(q4_3t.vcov, digits = 3)

# ╔═╡ 2265cfe2-13a7-11eb-3714-41524ab70ec5
begin
	scatter(df.weight, df.height, leg=false)
	quap4_3t = DataFrame(rand(q4_3t.distr, 10_000)', q4_3t.params)
	a_map = mean(quap4_3t.a)
	b_map = mean(quap4_3t.b)
	plot!(x, a_map .+ b_map .* (x .- x̄))
end

# ╔═╡ 22667578-13a7-11eb-2eaf-a3ef29e736b1
Text(precis(quap4_3t; io=String))

# ╔═╡ 06969514-13a8-11eb-31d5-836a86a69d0e
Text(precis(m4_3t; io=String))

# ╔═╡ e088c452-13a7-11eb-20e2-1bda82ac31f1
begin
	chns4_3t = sample(m4_3t, NUTS(), 1000)
	Particles(chns4_3t[[:a, :b, :σ]])
end

# ╔═╡ 000e6282-13b7-11eb-0d30-7fc6978b6e06
f(x) = mean(chns4_3t[:a]) +  mean(chns4_3t[:b]) * x  + 3 * randn();

# ╔═╡ 8c5c4506-13b7-11eb-28d7-53d9c6062c90
plot(x, f.(x))

# ╔═╡ d806ade6-13b4-11eb-02b3-f9a7bc51a2b6
begin
	Δ = 5.0
	x_test = [45.0, 50.0, 55.0, 60.0]
	y_test = f.(x_test)
	m_train = m4_3(x_test, y_test)
end

# ╔═╡ 31a07b4e-13b9-11eb-0f6c-17e99de04f60
m_test = m4_3(x_test, Vector{Union{Missing, Float64}}(undef, length(y_test)))

# ╔═╡ c9280bb0-13b6-11eb-1c8c-690a8151861f
predictions = Turing.Inference.predict(m_test, chns4_3t)

# ╔═╡ 2277eb66-13a7-11eb-200c-ebf1cc181a92
md"## End of m4.3t.jl"

# ╔═╡ Cell order:
# ╟─5a9929f0-13a6-11eb-18fd-6bd484be1583
# ╠═22398270-13a7-11eb-04c8-8b0995208358
# ╠═2239c3fc-13a7-11eb-1a00-4f5870c3b7ab
# ╠═223a5448-13a7-11eb-33d5-bf1dff954915
# ╠═2247ea56-13a7-11eb-1e8b-0388b6f47c98
# ╠═224c0bc0-13a7-11eb-3153-cf29938ecffb
# ╠═2253f11e-13a7-11eb-1c93-039f8e3bf927
# ╠═22548624-13a7-11eb-24d0-f90ba17217d8
# ╠═225e5a0a-13a7-11eb-140f-bde0811dd1e1
# ╠═2265cfe2-13a7-11eb-3714-41524ab70ec5
# ╠═22667578-13a7-11eb-2eaf-a3ef29e736b1
# ╠═06969514-13a8-11eb-31d5-836a86a69d0e
# ╠═e088c452-13a7-11eb-20e2-1bda82ac31f1
# ╠═000e6282-13b7-11eb-0d30-7fc6978b6e06
# ╠═8c5c4506-13b7-11eb-28d7-53d9c6062c90
# ╟─d806ade6-13b4-11eb-02b3-f9a7bc51a2b6
# ╠═31a07b4e-13b9-11eb-0f6c-17e99de04f60
# ╠═c9280bb0-13b6-11eb-1c8c-690a8151861f
# ╟─2277eb66-13a7-11eb-200c-ebf1cc181a92
