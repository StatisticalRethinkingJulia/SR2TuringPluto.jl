### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 0831131a-077d-11eb-2412-fbd0e7cc3436
using Pkg, DrWatson

# ╔═╡ 08314e98-077d-11eb-1d86-f9c9dd52e2c3
begin
	@quickactivate "SR2TuringPluto"
	using Turing
	using StatisticalRethinking
end

# ╔═╡ 041b4048-077d-11eb-20b8-6f8b8ee72626
md"## Clip-05-01-18t.jl"

# ╔═╡ bbbbd388-3ada-4108-bd5e-d9bc107f9d01
versioninfo()

# ╔═╡ 08320f9a-077d-11eb-2181-9bacc9603e02
begin
	df = CSV.read(sr_datadir("WaffleDivorce.csv"), DataFrame)
	df.D = zscore(df.Divorce)
	df.M = zscore(df.Marriage)
	df.A = zscore(df.MedianAgeMarriage)
end

# ╔═╡ 084594e8-077d-11eb-082b-6f19dc8a8e0f
std(df.MedianAgeMarriage)

# ╔═╡ 084663b4-077d-11eb-3586-73882ea2ceec
@model function m5_1_A(A, D)
    a ~ Normal(0, 0.2)
    bA ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bA)
    D ~ MvNormal(μ, σ)
end

# ╔═╡ 0857f156-077d-11eb-0829-799d64e773d2
begin
	m5_1_At = m5_1_A(df.A, df.D)
	prior5_1_At = sample(m5_1_At, Prior(), 50) |> DataFrame
	prior5_1_At = prior5_1_At[:, [:a, :bA, :σ]]
	Text(precis(prior5_1_At; io=String))
end

# ╔═╡ 35f34180-07d0-11eb-2c53-d1a2454235d5
begin
	x = -2:0.1:2
	plot(;dpi = 460)
	for r in eachrow(prior5_1_At)
		p = lin(r.a, x, r.bA)
		plot!(x, p, color = :black, alpha = 0.4)
	end
	plot!(legend = false)

	quap5_1_At = quap(m5_1_At)
	dfa5_1_At = DataFrame(rand(quap5_1_At.distr, 1000)', quap5_1_At.params)

	A_seq = range(-3, 3.2, length = 30)
	mu5_1_At = lin(dfa5_1_At.a', A_seq, dfa5_1_At.bA') |> meanlowerupper

	scatter(df.A, df.D, alpha = 0.4, legend = false)
	plot!(A_seq, mu5_1_At.mean, ribbon =
		(mu5_1_At.mean .- mu5_1_At.lower, mu5_1_At.upper .- mu5_1_At.mean);dpi = 460)
	vline!([0])
end

# ╔═╡ ffe4b816-07d5-11eb-351d-2fa77b983973
@model function m5_2_M(M, D)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, M, bM)
    D ~ MvNormal(μ, σ)
end

# ╔═╡ 0d7bab3c-07d7-11eb-3406-f3eea1352ba4
begin
	m5_2_Mt = m5_2_M(df.M, df.D)
	quap5_2_Mt = quap(m5_2_Mt)
end

# ╔═╡ 0d7c05c8-07d7-11eb-2253-437844d8bb9b
begin
	dfa5_2_Mt = DataFrame(rand(quap5_2_Mt.distr, 1000)', quap5_2_Mt.params)
	M_seq = range(-3, 3.2, length = 30)
	mu5_2_Mt = lin(dfa5_2_Mt.a', M_seq, dfa5_2_Mt.bM') |> meanlowerupper
	scatter(df.M, df.D, alpha = 0.4, legend = false)
	plot!(M_seq, mu5_2_Mt.mean, ribbon = 
		(mu5_2_Mt.mean .- mu5_2_Mt.lower, mu5_2_Mt.upper .- mu5_2_Mt.mean))
	vline!([0];dpi = 460)
end

# ╔═╡ 0db57a6a-07d7-11eb-2536-914facad3260
@model function m5_3_A_M(A, M, D)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    bA ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bA, M, bM)
    D ~ MvNormal(μ, σ)
end

# ╔═╡ 0dc1c694-07d7-11eb-1b23-fff550d68b74
begin
	m5_3_A_Mt = m5_3_A_M(df.A, df.M, df.D)
	quap5_3_A_Mt = quap(m5_3_A_Mt)
end

# ╔═╡ fa977d88-08d8-11eb-1ba3-3dbba5377f58
begin
	dfa5_3_A_Mt = DataFrame(rand(quap5_3_A_Mt.distr, 1000)', quap5_3_A_Mt.params)
	mu5_3_A_Mt = lin(dfa5_3_A_Mt.a', M_seq, dfa5_3_A_Mt.bM') |> meanlowerupper
	scatter(df.M, df.D, alpha = 0.4, legend = false)
	plot!(M_seq, mu5_3_A_Mt.mean, ribbon = 
		(mu5_3_A_Mt.mean .- mu5_3_A_Mt.lower, mu5_3_A_Mt.upper .- mu5_3_A_Mt.mean))
	vline!([0];dpi = 460)
end

# ╔═╡ 86a824ae-0991-11eb-3ee5-8ffb37cce3ec
begin
	quap5_1_3, fig = plotcoef(
		[m5_1_At, m5_2_Mt, m5_3_A_Mt], 
		[:a, :bM, :bA, :σ];
		func=quap
	)
	plot(fig)
end

# ╔═╡ 15658b5c-0993-11eb-21a9-0fc78a4a6fdf
quap5_1_3

# ╔═╡ 11f58eb0-08da-11eb-0f7f-239cec1bfdd5
begin
	age = randn(50)
	mar = rand.(Normal.(-age))
	div = rand.(Normal.(age))
end;

# ╔═╡ 156073d8-08da-11eb-18d2-ef100b4f33b0
@model function m5_4_AM(A, M)
    a ~ Normal(0, 0.2)
    bAM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bAM)
    M ~ MvNormal(μ, σ)
end

# ╔═╡ 1560bce4-08da-11eb-10b7-cf52ec6e7f33
begin
	sort!(df, :A)  # for nice looking graphs
	quap5_4_AMt = quap(m5_4_AM(df.A, df.M))
	post5_4_AMt = DataFrame(rand(quap5_4_AMt.distr, 1000)', quap5_4_AMt.params)
	mu5_4_AMt = lin(post5_4_AMt.a', df.A, post5_4_AMt.bAM') |> meanlowerupper
	resid = df.M .- mu5_4_AMt.mean

	scatter(df.A, df.M, legend = false)
	plot!(df.A, mu5_4_AMt.mean, ribbon =
		(mu5_4_AMt.mean .- mu5_4_AMt.lower, mu5_4_AMt.upper .- mu5_4_AMt.mean))

	post5_3_A_Mt = DataFrame(rand(quap5_3_A_Mt.distr, 1000)', quap5_3_A_Mt.params)
	mu5_3_A_M2t = lin(post5_3_A_Mt.a', df.A, post5_3_A_Mt.bA', df.M, post5_3_A_Mt.bM') |> meanlowerupper

	scatter(df.D, mu5_3_A_M2t.mean, 
		yerror = (mu5_3_A_M2t.mean .- mu5_3_A_M2t.lower, mu5_3_A_M2t.upper .- mu5_3_A_M2t.mean), 
		legend = false)
	plot!(identity, range(extrema(df.D)..., length = 10))
	plot!(xlabel = "Observed divorce", ylabel = "Predicted divorce")
end;

# ╔═╡ 701051c2-08da-11eb-33cb-41036a9f2724
begin
	xy = [df.D mu5_3_A_M2t.mean]
	for loc in ("ID", "UT", "ME", "RI")
		coord = xy[df.Loc .== loc, :]
		annotate!(coord..., loc, :bottom)
	end
	plot!(;dpi = 460)
end

# ╔═╡ c13c8b1e-08db-11eb-2725-b521e559429c
begin
	x_real = randn(100)
	x_spur = rand.(Normal.(x_real))
	y = rand.(Normal.(x_real))
	d = DataFrame((; y, x_real, x_spur))  # the semicolon turns it into a NamedTuple which then
										  # gives the DataFrame the names of the columns

	@df d corrplot([:y :x_real :x_spur])  # either
	corrplot(Matrix(d), label = names(d);dpi = 460) # or
end

# ╔═╡ 11d3317c-08df-11eb-1843-13d7b9b390ac
md"## End of clip-05-01-18t.jl"

# ╔═╡ Cell order:
# ╟─041b4048-077d-11eb-20b8-6f8b8ee72626
# ╠═bbbbd388-3ada-4108-bd5e-d9bc107f9d01
# ╠═0831131a-077d-11eb-2412-fbd0e7cc3436
# ╠═08314e98-077d-11eb-1d86-f9c9dd52e2c3
# ╠═08320f9a-077d-11eb-2181-9bacc9603e02
# ╠═084594e8-077d-11eb-082b-6f19dc8a8e0f
# ╠═084663b4-077d-11eb-3586-73882ea2ceec
# ╠═0857f156-077d-11eb-0829-799d64e773d2
# ╠═35f34180-07d0-11eb-2c53-d1a2454235d5
# ╠═ffe4b816-07d5-11eb-351d-2fa77b983973
# ╠═0d7bab3c-07d7-11eb-3406-f3eea1352ba4
# ╠═0d7c05c8-07d7-11eb-2253-437844d8bb9b
# ╠═0db57a6a-07d7-11eb-2536-914facad3260
# ╠═0dc1c694-07d7-11eb-1b23-fff550d68b74
# ╠═fa977d88-08d8-11eb-1ba3-3dbba5377f58
# ╠═86a824ae-0991-11eb-3ee5-8ffb37cce3ec
# ╠═15658b5c-0993-11eb-21a9-0fc78a4a6fdf
# ╠═11f58eb0-08da-11eb-0f7f-239cec1bfdd5
# ╠═156073d8-08da-11eb-18d2-ef100b4f33b0
# ╠═1560bce4-08da-11eb-10b7-cf52ec6e7f33
# ╠═701051c2-08da-11eb-33cb-41036a9f2724
# ╠═c13c8b1e-08db-11eb-2725-b521e559429c
# ╟─11d3317c-08df-11eb-1843-13d7b9b390ac
