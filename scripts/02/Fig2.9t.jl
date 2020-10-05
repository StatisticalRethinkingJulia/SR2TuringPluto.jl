
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
end

md"## Fig2.9t.jl"

md"##### Comparisons between different quadratic approximations."

md"##### Grid of 1001 steps."

p_grid = range(0, step=0.001, stop=1);

md"##### Set all priors = 1.0."

prior = ones(length(p_grid));

md"##### Binomial pdf."

likelihood = pdf.(Binomial.(9, p_grid), 6);

md"##### As Uniform prior has been used, unstandardized posterior is equal to likelihood."

begin
	posterior = likelihood .* prior
	
	# Scale posterior such that they become probabilities."
	
	posterior = posterior / sum(posterior)
end;

md"##### Sample using the computed posterior values as weights."

begin
	N = 10000
	samples = sample(p_grid, Weights(posterior), N);
	chn = MCMCChains.Chains(reshape(samples, N, 1, 1), ["toss"]);
end;

md"##### Describe the chain."

chn

md"##### Compute the MAP (maximum a posteriori) estimate."

function loglik(x)
  ll = 0.0
  ll += log.(pdf.(Beta(1, 1), x[1]))
  ll += sum(log.(pdf.(Binomial(9, x[1]), repeat([6], 1))))
  -ll
end

opt = optimize(loglik, 0.0, 1.0)

qmap = Optim.minimizer(opt)

md"##### Fit Optim quadratic approximation."

quapfit = [qmap[1], std(samples, mean=qmap[1])]

f = Vector{Plots.Plot{Plots.GRBackend}}(undef, 4);

md"##### Fit Turing quap approximation."

begin
	w = 6
	l = 3
	n = 9
	x = 0:0.01:1
end

@model globethrowing(w, l) = begin
    p ~ Uniform(0, 1)
    w ~ Binomial(w + l, p)
end

begin
	m = globethrowing(6, 3)
	r = quap(m)
end

md"##### Setup the plots."

begin
	f[1] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0), 
	  lab="Conjugate solution", leg=:topleft)
	density!(f[1], samples, lab="Sample density")

	# Quadratic approximation using Optim

	f[2] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0),
	  lab="Conjugate solution", leg=:topleft)
	plot!( f[2], x, pdf.(Normal( quapfit[1], quapfit[2] ) , x ),
	  lab="Optim logpdf approx.")

	# Quadratic approximation using Turing.jl quap()

	f[3] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0),
	  lab="Conjugate solution", leg=:topleft)
	plot!( f[3], x, pdf.(Normal([r.coef.p][1], √collect(reshape(r.vcov, 1, 1))[1]) , x ), 
		lab="Turing quap approx.")

	# MLE approximation

	f[4] = plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), xlims=(-0.5, 1.0),
	  lab="Conjugate solution", leg=:topleft)
	fi = fit(Normal, samples)
	plot!(f[4], x, pdf.(Normal( fi.μ , fi.σ ) , x ), lab="Normal MLE approx.")
end;

plot(f..., layout=(2, 2))

md"## End of Fig2.9t.jl"

