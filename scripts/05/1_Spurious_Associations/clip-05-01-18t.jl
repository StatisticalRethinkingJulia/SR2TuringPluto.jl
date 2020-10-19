
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Clip-05-01-18t.jl"

begin
	df = CSV.read(sr_datadir("WaffleDivorce.csv"), DataFrame)
	df.D = zscore(df.Divorce)
	df.M = zscore(df.Marriage)
	df.A = zscore(df.MedianAgeMarriage)
end

std(df.MedianAgeMarriage)

@model function m5_1_A(A, D)
    a ~ Normal(0, 0.2)
    bA ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bA)
    D ~ MvNormal(μ, σ)
end

begin
	m5_1_At = m5_1_A(df.A, df.D)
	prior5_1_At = sample(m5_1_At, Prior(), 50) |> DataFrame
	prior5_1_At = prior5_1_At[:, [:a, :bA, :σ]]
	Text(precis(prior5_1_At; io=String))
end

begin
	x = -2:0.1:2
	plot()
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
		(mu5_1_At.mean .- mu5_1_At.lower, mu5_1_At.upper .- mu5_1_At.mean))
	vline!([0])
end

@model function m5_2_M(M, D)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, M, bM)
    D ~ MvNormal(μ, σ)
end

begin
	m5_2_Mt = m5_2_M(df.M, df.D)
	quap5_2_Mt = quap(m5_2_Mt)
end

begin
	dfa5_2_Mt = DataFrame(rand(quap5_2_Mt.distr, 1000)', quap5_2_Mt.params)
	M_seq = range(-3, 3.2, length = 30)
	mu5_2_Mt = lin(dfa5_2_Mt.a', M_seq, dfa5_2_Mt.bM') |> meanlowerupper
	scatter(df.M, df.D, alpha = 0.4, legend = false)
	plot!(M_seq, mu5_2_Mt.mean, ribbon = 
		(mu5_2_Mt.mean .- mu5_2_Mt.lower, mu5_2_Mt.upper .- mu5_2_Mt.mean))
	vline!([0])
end

@model function m5_3_A_M(A, M, D)
    a ~ Normal(0, 0.2)
    bM ~ Normal(0, 0.5)
    bA ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bA, M, bM)
    D ~ MvNormal(μ, σ)
end

begin
	m5_3_A_Mt = m5_3_A_M(df.A, df.M, df.D)
	quap5_3_A_Mt = quap(m5_3_A_Mt)
end

begin
	dfa5_3_A_Mt = DataFrame(rand(quap5_3_A_Mt.distr, 1000)', quap5_3_A_Mt.params)
	mu5_3_A_Mt = lin(dfa5_3_A_Mt.a', M_seq, dfa5_3_A_Mt.bM') |> meanlowerupper
	scatter(df.M, df.D, alpha = 0.4, legend = false)
	plot!(M_seq, mu5_3_A_Mt.mean, ribbon = 
		(mu5_3_A_Mt.mean .- mu5_3_A_Mt.lower, mu5_3_A_Mt.upper .- mu5_3_A_Mt.mean))
	vline!([0])
end

begin
	quap5_1_3, fig = plotcoef(
		[m5_1_At, m5_2_Mt, m5_3_A_Mt], 
		[:a, :bM, :bA, :σ];
		func=quap
	)
	plot(fig)
end

quap5_1_3

begin
	age = randn(50)
	mar = rand.(Normal.(-age))
	div = rand.(Normal.(age))
end;

@model function m5_4_AM(A, M)
    a ~ Normal(0, 0.2)
    bAM ~ Normal(0, 0.5)
    σ ~ Exponential(1)
    μ = lin(a, A, bAM)
    M ~ MvNormal(μ, σ)
end

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

begin
	xy = [df.D mu5_3_A_M2t.mean]
	for loc in ("ID", "UT", "ME", "RI")
		coord = xy[df.Loc .== loc, :]
		annotate!(coord..., loc, :bottom)
	end
	plot!()
end

begin
	x_real = randn(100)
	x_spur = rand.(Normal.(x_real))
	y = rand.(Normal.(x_real))
	d = DataFrame((; y, x_real, x_spur))  # the semicolon turns it into a NamedTuple which then
										  # gives the DataFrame the names of the columns

	@df d corrplot([:y :x_real :x_spur])  # either
	corrplot(Matrix(d), label = names(d)) # or
end

md"## End of clip-05-01-18t.jl"

