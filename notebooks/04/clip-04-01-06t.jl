### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 3bd0847c-f2b4-11ea-35d8-c7657a170cf9
using DrWatson

# ╔═╡ 3bd0c52c-f2b4-11ea-09e6-05dbd556433f
begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
	using PlutoUI
end

# ╔═╡ 28342aca-f2b1-11ea-342d-95590e306ff4
md"## Clip-04-01-06s.jl"

# ╔═╡ 3bd18872-f2b4-11ea-2271-470ff5d51737
md"### snippet 4.1"

# ╔═╡ 3be20a12-f2b4-11ea-2179-2995ca45302f
md"###### No attempt has been made to condense this too fewer lines of code."

# ╔═╡ 4114e7b0-f85f-11ea-248a-e1c6723d2f1a
@bind noofwalks Slider(5:200, default=9)

# ╔═╡ 3be2d1ea-f2b4-11ea-09ce-1db161f69c7e
begin
	noofsteps = 20;
	pos = Array{Float64, 2}(rand(Uniform(-1, 1), noofsteps, noofwalks));
	pos[1, :] = zeros(noofwalks);
	csum = cumsum(pos, dims=1);
	mx = minimum(csum) * 0.9
end;

# ╔═╡ 3bf08696-f2b4-11ea-0fae-815209e20f46
md"###### Plot and annotate the random walks."

# ╔═╡ 3bfe421a-f2b4-11ea-12c7-f7f775c210a9
md"###### Generate 3 plots of densities at 3 different step numbers (4, 8 and 16)."

# ╔═╡ 3bff112a-f2b4-11ea-338a-791bc65b719f
begin
	f = Plots.font("DejaVu Sans", 6)
	p1 = plot(csum, leg=false, title="Random walks ($(noofwalks))")
	plot!(p1, csum[:, Int(floor(noofwalks/2))], leg=false, title="Random walks ($(noofwalks))", 				color=:black)
	plot!(p1, [5], seriestype="vline")
	annotate!(5, mx, text("step 4", f, :left))
	plot!(p1, [9], seriestype="vline")
	annotate!(9, mx, text("step 8", f, :left))
	plot!(p1, [17], seriestype="vline")
	annotate!(17, mx, text("step 16", f, :left))

	p2 = Vector{Plots.Plot{Plots.GRBackend}}(undef, 3);
	plt = 1
	for step in [4, 8, 16]
		indx = step + 1 								# We added the first line of zeros
		global plt
	  	fitl = fit_mle(Normal, csum[indx, :])
	  	lx = (fitl.μ-4*fitl.σ):0.01:(fitl.μ+4*fitl.σ)
	  	p2[plt] = density(csum[indx, :], legend=false, title="$(step) steps")
	 	plot!( p2[plt], lx, pdf.(Normal( fitl.μ , fitl.σ ) , lx ), fill=(0, .5,:orange))
	  	plt += 1
	end
	p3 = plot(p2..., layout=(1, 3))
	plot(p1, p3, layout=(2,1))
end

# ╔═╡ 3c10cb0e-f2b4-11ea-210f-1f6986fbb7d8
md"## snippet 4.2"

# ╔═╡ 3c1b683e-f2b4-11ea-206b-b9178296ec19
prod(1 .+ rand(Uniform(0, 0.1), 12))

# ╔═╡ 3c248840-f2b4-11ea-1397-9f2c313a2676
md"## snippet 4.3"

# ╔═╡ 3c26427c-f2b4-11ea-23f9-6303c62de7f6
begin
	growth = [prod(1 .+ rand(Uniform(0, 0.1), 12)) for i in 1:10000];
	fit2 = fit_mle(Normal, growth)
	plot(Normal(fit2.μ , fit2.σ ), fill=(0, .5,:orange), lab="Normal distribution")
	density!(growth, lab="'sample' distribution")
end

# ╔═╡ 3c37aa58-f2b4-11ea-2e87-a11117582011
md"## snippet 4.4"

# ╔═╡ 3c390506-f2b4-11ea-3232-296f6952aaa7
begin
	big = [prod(1 .+ rand(Uniform(0, 0.5), 12)) for i in 1:10000];
	small = [prod(1 .+ rand(Uniform(0, 0.01), 12)) for i in 1:10000];
	fitb = fit_mle(Normal, big)
	fits = fit_mle(Normal, small)
	p5 = plot(Normal(fitb.μ , fitb.σ ), lab="Big normal distribution", fill=(0, .5,:orange))
	p4 = plot(Normal(fits.μ , fits.σ ), lab="Small normal distribution", fill=(0, .5,:orange))
	density!(p5, big, lab="'big' distribution")
	density!(p4, small, lab="'small' distribution")
	plot(p5, p4, layout=(1, 2))
end

# ╔═╡ 3c4245ee-f2b4-11ea-2add-750014544e83
md"## snippet 4.5"

# ╔═╡ 3c4be516-f2b4-11ea-086b-f1461fa2e55e
begin
	log_big = [log(prod(1 .+ rand(Uniform(0, 0.5), 12))) for i in 1:10000];
	fit3 = fit_mle(Normal, log_big)
	plot(Normal(fit3.μ , fit3.σ ), fill=(0, .5,:orange), lab="Normal distribution")
	density!(log_big, lab="'sample' distribution")
end

# ╔═╡ 3c567b9a-f2b4-11ea-2798-b5b3bffc3f84
md"## snippet 4.6"

# ╔═╡ 3c60487a-f2b4-11ea-32f7-abda1ad883cf
md"###### Grid of 1001 steps."

# ╔═╡ 3c6a0c78-f2b4-11ea-12a9-a506e7a7cacb
p_grid = range(0, step=0.001, stop=1);

# ╔═╡ 3c7435ce-f2b4-11ea-36a4-e72a5024cfa3
md"###### All priors = 1.0."

# ╔═╡ 3c7f4458-f2b4-11ea-3006-fb69892e25da
prior = ones(length(p_grid));

# ╔═╡ 3c89c428-f2b4-11ea-3216-7dab773dfe79
md"###### Binomial pdf."

# ╔═╡ 3c9436d8-f2b4-11ea-3b46-fbfefa79d535
likelihood = [pdf(Binomial(9, p), 6) for p in p_grid];

# ╔═╡ 3ca753ee-f2b4-11ea-21f5-93f4fee09437
md"###### A Uniform prior has been used, unstandardized posterior is equal to likelihood."

# ╔═╡ 3cab5208-f2b4-11ea-0cff-ddab2410c607
posterior = likelihood .* prior;

# ╔═╡ 3cb6dac6-f2b4-11ea-0c94-a109c74c2d22
md"###### Scale posterior such that they become probabilities."

# ╔═╡ 3cc23a06-f2b4-11ea-3dd3-c5ca6ebe71c3
posterior2 = posterior / sum(posterior);

# ╔═╡ 3cd03f70-f2b4-11ea-1877-11b0102993c8
md"###### Sample using the computed posterior values as weights. In this example we keep the number of samples equal to the length of p_grid, but that is not required."

# ╔═╡ 3cdd29a6-f2b4-11ea-071d-fb6328c01d5f
begin
	samples = sample(p_grid, Weights(posterior), length(p_grid));
	p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 2)
	p[1] = scatter(1:length(p_grid), samples, markersize = 2, ylim=(0.0, 1.3), lab="Draws")
end

# ╔═╡ 3ce9eac4-f2b4-11ea-0bf6-9f130255ec5a
md"###### Analytical calculation."

# ╔═╡ 3cf665a6-f2b4-11ea-29ad-19ecaf1f423a
begin
	w = 6
	n = 9
	x = 0:0.01:1
	p[2] = density(samples, ylim=(0.0, 5.0), lab="Sample density")
	p[2] = plot!( x, pdf.(Beta( w+1 , n-w+1 ) , x ), lab="Conjugate solution")
end

# ╔═╡ 3d03bb5c-f2b4-11ea-315e-6fac666895eb
md"###### Quadratic approximation."

# ╔═╡ 3d11472c-f2b4-11ea-3bf6-9d517b13b291
begin
	plot!( p[2], x, pdf.(Normal( 0.67 , 0.16 ) , x ), lab="Normal approximation", fill=(0, .5,:orange))
	plot(p..., layout=(1, 2))
end

# ╔═╡ 3d1db2f0-f2b4-11ea-3bde-29b8f0be3087
md"## End of clip-04-01-06s.jl"

# ╔═╡ Cell order:
# ╟─28342aca-f2b1-11ea-342d-95590e306ff4
# ╠═3bd0847c-f2b4-11ea-35d8-c7657a170cf9
# ╠═3bd0c52c-f2b4-11ea-09e6-05dbd556433f
# ╟─3bd18872-f2b4-11ea-2271-470ff5d51737
# ╟─3be20a12-f2b4-11ea-2179-2995ca45302f
# ╠═4114e7b0-f85f-11ea-248a-e1c6723d2f1a
# ╠═3be2d1ea-f2b4-11ea-09ce-1db161f69c7e
# ╟─3bf08696-f2b4-11ea-0fae-815209e20f46
# ╟─3bfe421a-f2b4-11ea-12c7-f7f775c210a9
# ╠═3bff112a-f2b4-11ea-338a-791bc65b719f
# ╟─3c10cb0e-f2b4-11ea-210f-1f6986fbb7d8
# ╠═3c1b683e-f2b4-11ea-206b-b9178296ec19
# ╠═3c248840-f2b4-11ea-1397-9f2c313a2676
# ╠═3c26427c-f2b4-11ea-23f9-6303c62de7f6
# ╠═3c37aa58-f2b4-11ea-2e87-a11117582011
# ╠═3c390506-f2b4-11ea-3232-296f6952aaa7
# ╠═3c4245ee-f2b4-11ea-2add-750014544e83
# ╠═3c4be516-f2b4-11ea-086b-f1461fa2e55e
# ╟─3c567b9a-f2b4-11ea-2798-b5b3bffc3f84
# ╟─3c60487a-f2b4-11ea-32f7-abda1ad883cf
# ╠═3c6a0c78-f2b4-11ea-12a9-a506e7a7cacb
# ╟─3c7435ce-f2b4-11ea-36a4-e72a5024cfa3
# ╠═3c7f4458-f2b4-11ea-3006-fb69892e25da
# ╟─3c89c428-f2b4-11ea-3216-7dab773dfe79
# ╠═3c9436d8-f2b4-11ea-3b46-fbfefa79d535
# ╟─3ca753ee-f2b4-11ea-21f5-93f4fee09437
# ╠═3cab5208-f2b4-11ea-0cff-ddab2410c607
# ╠═3cb6dac6-f2b4-11ea-0c94-a109c74c2d22
# ╠═3cc23a06-f2b4-11ea-3dd3-c5ca6ebe71c3
# ╠═3cd03f70-f2b4-11ea-1877-11b0102993c8
# ╠═3cdd29a6-f2b4-11ea-071d-fb6328c01d5f
# ╟─3ce9eac4-f2b4-11ea-0bf6-9f130255ec5a
# ╠═3cf665a6-f2b4-11ea-29ad-19ecaf1f423a
# ╟─3d03bb5c-f2b4-11ea-315e-6fac666895eb
# ╠═3d11472c-f2b4-11ea-3bf6-9d517b13b291
# ╟─3d1db2f0-f2b4-11ea-3bde-29b8f0be3087
