
using Markdown
using InteractiveUtils

macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

using DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StatisticalRethinking
	using PlutoUI
end

md"## Clip-04-01-06s.jl"

md"### snippet 4.1"

md"###### No attempt has been made to condense this too fewer lines of code."

@bind noofwalks Slider(5:200, default=9)

begin
	noofsteps = 20;
	pos = Array{Float64, 2}(rand(Uniform(-1, 1), noofsteps, noofwalks));
	pos[1, :] = zeros(noofwalks);
	csum = cumsum(pos, dims=1);
	mx = minimum(csum) * 0.9
end;

md"###### Plot and annotate the random walks."

md"###### Generate 3 plots of densities at 3 different step numbers (4, 8 and 16)."

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

md"## snippet 4.2"

prod(1 .+ rand(Uniform(0, 0.1), 12))

md"## snippet 4.3"

begin
	growth = [prod(1 .+ rand(Uniform(0, 0.1), 12)) for i in 1:10000];
	fit2 = fit_mle(Normal, growth)
	plot(Normal(fit2.μ , fit2.σ ), fill=(0, .5,:orange), lab="Normal distribution")
	density!(growth, lab="'sample' distribution")
end

md"## snippet 4.4"

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

md"## snippet 4.5"

begin
	log_big = [log(prod(1 .+ rand(Uniform(0, 0.5), 12))) for i in 1:10000];
	fit3 = fit_mle(Normal, log_big)
	plot(Normal(fit3.μ , fit3.σ ), fill=(0, .5,:orange), lab="Normal distribution")
	density!(log_big, lab="'sample' distribution")
end

md"## snippet 4.6"

md"###### Grid of 1001 steps."

p_grid = range(0, step=0.001, stop=1);

md"###### All priors = 1.0."

prior = ones(length(p_grid));

md"###### Binomial pdf."

likelihood = [pdf(Binomial(9, p), 6) for p in p_grid];

md"###### A Uniform prior has been used, unstandardized posterior is equal to likelihood."

posterior = likelihood .* prior;

md"###### Scale posterior such that they become probabilities."

posterior2 = posterior / sum(posterior);

md"###### Sample using the computed posterior values as weights. In this example we keep the number of samples equal to the length of p_grid, but that is not required."

begin
	samples = sample(p_grid, Weights(posterior), length(p_grid));
	p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 2)
	p[1] = scatter(1:length(p_grid), samples, markersize = 2, ylim=(0.0, 1.3), lab="Draws")
end

md"###### Analytical calculation."

begin
	w = 6
	n = 9
	x = 0:0.01:1
	p[2] = density(samples, ylim=(0.0, 5.0), lab="Sample density")
	p[2] = plot!( x, pdf.(Beta( w+1 , n-w+1 ) , x ), lab="Conjugate solution")
end

md"###### Quadratic approximation."

begin
	plot!( p[2], x, pdf.(Normal( 0.67 , 0.16 ) , x ), lab="Normal approximation", fill=(0, .5,:orange))
	plot(p..., layout=(1, 2))
end

md"## End of clip-04-01-06s.jl"

