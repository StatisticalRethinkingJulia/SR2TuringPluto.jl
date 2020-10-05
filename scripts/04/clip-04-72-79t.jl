
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## clip-04-72-79t"

md"### snippet 4.72"

df = CSV.read(sr_datadir("cherry_blossoms.csv"), DataFrame; missingstring = "NA");

scatter(df.year, df.doy, leg=false)

md"### snippet 4.73"

begin
	dropmissing!(df, :doy)
	num_knots = 15
	knot_list = quantile(df.year, range(0, 1, length = num_knots))
	basis = BSplineBasis(4, knot_list)
	B = basismatrix(basis, df.year)
end;

begin
	plot(legend = false, xlabel = "year", ylabel = "basis value")
	for y in eachcol(B)
		plot!(df.year, y)
	end
	plot!()
end

md"## snippet 4.76"

@model function spline_1(D, B = B)
    α ~ Normal(100, 10)
    w ~ filldist(Normal(0, 10), size(B, 2))
    σ ~ Exponential(1)
    μ = α .+ B * w
    D ~ MvNormal(μ, σ)
    return μ
end

md"### snippet 4.77"

begin
	quap4_7t = quap(spline_1(df.doy))
	quap4_7t.coef
end

begin
	w_str = ["w[$i]" for i in 1:length(basis)]
	post_3 = DataFrame(rand(quap4_7t.distr, 1000)', ["α"; w_str; "σ"])
	w_3 = mean.(eachcol(post_3[:, w_str]))              # either
	w_3 = [mean(post_3[:, col]) for col in w_str]       # or
	plot(legend = false, xlabel = "year", ylabel = "basis * weight")
	for y in eachcol(B .* w_3')
		plot!(df3.year, y)
	end
	plot!()
end	

md"### snippet 4.78"

begin
	mu_3 = post_3.α' .+ B * Array(post_3[!, w_str])'
	mu_3 = meanlowerupper(mu_3)
	scatter(df.year, df.doy, alpha = 0.3)
	plot!(df.year, mu_3.mean, ribbon = (mu_3.mean .- mu_3.lower, mu_3.upper .- mu_3.mean))
end

md"### snippet 4.79"

@model function spline_2(D, B = B)
    α ~ Normal(100, 10)
    w ~ filldist(Normal(0, 10), size(B, 2))
    σ ~ Exponential(1)
    μ = [α + sum(Brow .* w) for Brow in eachrow(B)]
    D ~ MvNormal(μ, σ)
end

begin
	quap_alt_4_7t = quap(spline_2(df.doy))
	quap_alt_4_7t.coef
end

md"## End of clip-04-72-79t.jl"

