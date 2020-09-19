### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 263e3c4c-f826-11ea-246e-db6675dc4719
using Pkg, DrWatson

# ╔═╡ 263e7f5c-f826-11ea-1794-2daccdb57f16
begin
	@quickactivate "StatisticalRethinkingTuring"
	using StanSample
	using StatisticalRethinking
end

# ╔═╡ d5b82444-f824-11ea-015a-9b79e6e41731
md"## Clip-04-07-15s.jl"

# ╔═╡ 263f6d54-f826-11ea-292b-517d8fe29e05
md"### snippet 4.7"

# ╔═╡ 26472062-f826-11ea-262e-2b65a759d719
d = CSV.read(sr_datadir("Howell1.csv"), DataFrame);

# ╔═╡ 265325d8-f826-11ea-25c5-51983ce8514f
md"### snippet 4.8"

# ╔═╡ 265d006c-f826-11ea-05b8-016f5f9e15da
md"##### Show a summary of the  DataFrame."

# ╔═╡ 265e60ce-f826-11ea-38f7-cfeb5f92ea85
Particles(d)

# ╔═╡ 266ac76a-f826-11ea-2d0e-03c06ea11c56
md"### snippet 4.9"

# ╔═╡ 8efcd538-f863-11ea-1677-111b9890814e
Text(precis(d; io=String))

# ╔═╡ 266b8554-f826-11ea-26c9-ffc3b0ab57d6
md"##### Compare with describe():"

# ╔═╡ 267b9d68-f826-11ea-0d91-e117f7f4b66e
describe(d, :all)

# ╔═╡ 267d089e-f826-11ea-1788-a3cf9d04a192
md"### snippet 4.10"

# ╔═╡ 268c8e40-f826-11ea-2499-132b929e08ac
d.height

# ╔═╡ 268d5b18-f826-11ea-23e1-fdd6f004869c
md"### snippet 4.11"

# ╔═╡ 26950052-f826-11ea-014b-f1b92eb7b7fd
md"##### Adults only."

# ╔═╡ 26a496fa-f826-11ea-36eb-2bebacba65e7
begin
	df = filter(row -> row[:age] >= 18, d);
	Particles(df)
end

# ╔═╡ 26a5dfa8-f826-11ea-0c29-d57b72e30315
Text(precis(df; io=String))

# ╔═╡ 26ae2848-f826-11ea-00f0-77f23ed3ccf3
md"##### Our model:"

# ╔═╡ 26bc65ae-f826-11ea-1a73-035c3f4beff9
m4_1_book = "
  height ~ Normal(μ, σ) # likelihood
  μ ~ Normal(178,20) # prior
  σ ~ Uniform(0, 50) # prior
";

# ╔═╡ 26bfc9cc-f826-11ea-2378-3da2914964d0
md"##### Plot the prior densities."

# ╔═╡ 26c847aa-f826-11ea-2051-93748f099c87
p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 4);

# ╔═╡ 26d1052a-f826-11ea-14ec-db479de75f85
md"### snippet 4.12"

# ╔═╡ 26da7a74-f826-11ea-04de-c1af74ca5f52
md"##### μ prior."

# ╔═╡ 26e489d8-f826-11ea-2a8e-63dc68ebdfcb
d1 = Normal(178, 20)

# ╔═╡ 26ee436a-f826-11ea-2b51-553109976eeb
p[1] = plot(100:250, [pdf(d1, μ) for μ in 100:250],
	xlab="mu",
	ylab="density",
	lab="Prior on mu");

# ╔═╡ 26f91aa6-f826-11ea-1b96-2be9514b1f4c
md"### snippet 4.13"

# ╔═╡ 2702ab5e-f826-11ea-007a-878156afed4e
md"##### Show σ  prior."

# ╔═╡ 270cde42-f826-11ea-25f3-71a99eae18e2
begin
	d2 = Uniform(0, 50)
	p[2] = plot(-5:0.1:55, [pdf(d2, σ) for σ in 0-5:0.1:55],
		xlab="sigma",
		ylab="density",
		lab="Prior on sigma")
end;

# ╔═╡ 2717e74c-f826-11ea-26d5-c5862cdeef06
md"### snippet 4.14."

# ╔═╡ 27222c68-f826-11ea-0c72-392e7c56f262
begin
	sample_mu_20 = rand(d1, 10000)
	sample_sigma = rand(d2, 10000)

	d3 = Normal(178, 100)
	sample_mu_100 = rand(d3, 10000)

	prior_height_20 = [rand(Normal(sample_mu_20[i], sample_sigma[i]), 1)[1] for i in 1:10000]
	p[3] = density(prior_height_20,
		xlab="height",
		ylab="density",
		lab="Prior predictive height")
end;

# ╔═╡ 272c9e2e-f826-11ea-221d-afb68d7bd316
begin
	prior_height_100 = [rand(Normal(sample_mu_100[i], sample_sigma[i]), 1)[1] for i in 1:10000]
	p[4] = density(prior_height_100,
		xlab="height",
		ylab="density",
		lab="Prior predictive mu")
end;

# ╔═╡ 2737cdaa-f826-11ea-2c70-c5c3d16d5a36
plot(p..., layout=(2,2))

# ╔═╡ 274bc54e-f826-11ea-026e-9f45eeb84c1b
md"##### Store in a DataFrame."

# ╔═╡ 2755d188-f826-11ea-0f5f-e11a4ffec021
df2 = DataFrame(
	mu_20 = sample_mu_20,
	mu_100 = sample_mu_100,
	sigma=sample_sigma,
	prior_height_20=prior_height_20,
	prior_height_100=prior_height_100);

# ╔═╡ 2763c7ac-f826-11ea-13ce-61b6a5ee2379
precis(df2)

# ╔═╡ 276f9e94-f826-11ea-2421-218735527bc2
md"##### On to Stan."

# ╔═╡ 277b474c-f826-11ea-1816-e518c1ed26eb
md"##### Recall our model:"

# ╔═╡ 278ac238-f826-11ea-0f32-799e4c8f10f5
m4_1_rethinking = "
  # Priors
  μ ~ Normal(178,20)
  σ ~ Uniform(0, 50)

  # Likelihood of data
  height ~ Normal(μ, σ)
";

# ╔═╡ 2797fe8e-f826-11ea-3f7b-b50d4cb27e92
m4_1 = "
// Inferring the mean and std
data {
  int N;
  real<lower=0> h[N];
}
parameters {
  real<lower=0> sigma;
  real<lower=0> sigma_prior;
  real<lower=0,upper=250> mu;
  real<lower=0,upper=250> mu_prior;

}
model {
  // Priors for mu
  mu ~ normal(178, 20);
  mu_prior ~ normal(178, 20);

  // Priors for sigma
  sigma ~ uniform( 0 , 50 );
  sigma_prior ~ uniform( 0 , 50 );

  // Observed heights
  h ~ normal(mu, sigma);
}
";

# ╔═╡ 27a52404-f826-11ea-0d0d-8530041392b5
md"##### Create a StanSample SampleModel:"

# ╔═╡ 27b38be8-f826-11ea-19d6-fd7d10c22e6f
m4_1s = SampleModel("heights", m4_1);

# ╔═╡ 27c90360-f826-11ea-0d57-a10c128fca32
md"##### Package the data:"

# ╔═╡ 27cfa486-f826-11ea-0083-3fcaa67feb12
heightsdata = Dict("N" => length(df[:, :height]), "h" => df[:, :height]);

# ╔═╡ 27dbd2a6-f826-11ea-2e23-59d7c6d2aac9
md"##### Run Stan's cmdstan:"

# ╔═╡ 27e9d8ea-f826-11ea-3cbe-31d9a7852345
rc = stan_sample(m4_1s, data=heightsdata);

# ╔═╡ 27f78c12-f826-11ea-1cd1-eda91d47196f
md"##### Check if sampling went ok:"

# ╔═╡ 280705b8-f826-11ea-1d51-251fdd5eff2e
md"##### Read in the samples and show a chain summary."

# ╔═╡ dc3a5122-f82c-11ea-1eff-c745d65ab11b
success(rc) && (chn = read_samples(m4_1s; output_format=:mcmcchains))

# ╔═╡ dc3af532-f82c-11ea-3212-f1b3c852513b
md"##### Plot the sampling trace."

# ╔═╡ 044a19a2-f866-11ea-2b89-a51866d89a50
plot(chn, seriestype = :traceplot)

# ╔═╡ dc476376-f82c-11ea-10ac-97fcb8c78627
md"##### Plot the density of posterior draws."

# ╔═╡ 25318e22-f866-11ea-015b-d736c83ebfaa
plot(chn, seriestype = :density)

# ╔═╡ 2814f2fc-f826-11ea-3fbc-0541fe904b97
md"## End of clip-04-07-15s.jl"

# ╔═╡ Cell order:
# ╠═d5b82444-f824-11ea-015a-9b79e6e41731
# ╠═263e3c4c-f826-11ea-246e-db6675dc4719
# ╠═263e7f5c-f826-11ea-1794-2daccdb57f16
# ╟─263f6d54-f826-11ea-292b-517d8fe29e05
# ╠═26472062-f826-11ea-262e-2b65a759d719
# ╟─265325d8-f826-11ea-25c5-51983ce8514f
# ╟─265d006c-f826-11ea-05b8-016f5f9e15da
# ╠═265e60ce-f826-11ea-38f7-cfeb5f92ea85
# ╟─266ac76a-f826-11ea-2d0e-03c06ea11c56
# ╠═8efcd538-f863-11ea-1677-111b9890814e
# ╟─266b8554-f826-11ea-26c9-ffc3b0ab57d6
# ╠═267b9d68-f826-11ea-0d91-e117f7f4b66e
# ╟─267d089e-f826-11ea-1788-a3cf9d04a192
# ╠═268c8e40-f826-11ea-2499-132b929e08ac
# ╟─268d5b18-f826-11ea-23e1-fdd6f004869c
# ╟─26950052-f826-11ea-014b-f1b92eb7b7fd
# ╠═26a496fa-f826-11ea-36eb-2bebacba65e7
# ╠═26a5dfa8-f826-11ea-0c29-d57b72e30315
# ╟─26ae2848-f826-11ea-00f0-77f23ed3ccf3
# ╠═26bc65ae-f826-11ea-1a73-035c3f4beff9
# ╟─26bfc9cc-f826-11ea-2378-3da2914964d0
# ╠═26c847aa-f826-11ea-2051-93748f099c87
# ╟─26d1052a-f826-11ea-14ec-db479de75f85
# ╟─26da7a74-f826-11ea-04de-c1af74ca5f52
# ╠═26e489d8-f826-11ea-2a8e-63dc68ebdfcb
# ╠═26ee436a-f826-11ea-2b51-553109976eeb
# ╟─26f91aa6-f826-11ea-1b96-2be9514b1f4c
# ╟─2702ab5e-f826-11ea-007a-878156afed4e
# ╠═270cde42-f826-11ea-25f3-71a99eae18e2
# ╟─2717e74c-f826-11ea-26d5-c5862cdeef06
# ╠═27222c68-f826-11ea-0c72-392e7c56f262
# ╠═272c9e2e-f826-11ea-221d-afb68d7bd316
# ╠═2737cdaa-f826-11ea-2c70-c5c3d16d5a36
# ╠═274bc54e-f826-11ea-026e-9f45eeb84c1b
# ╠═2755d188-f826-11ea-0f5f-e11a4ffec021
# ╠═2763c7ac-f826-11ea-13ce-61b6a5ee2379
# ╠═276f9e94-f826-11ea-2421-218735527bc2
# ╟─277b474c-f826-11ea-1816-e518c1ed26eb
# ╠═278ac238-f826-11ea-0f32-799e4c8f10f5
# ╠═2797fe8e-f826-11ea-3f7b-b50d4cb27e92
# ╟─27a52404-f826-11ea-0d0d-8530041392b5
# ╠═27b38be8-f826-11ea-19d6-fd7d10c22e6f
# ╟─27c90360-f826-11ea-0d57-a10c128fca32
# ╠═27cfa486-f826-11ea-0083-3fcaa67feb12
# ╟─27dbd2a6-f826-11ea-2e23-59d7c6d2aac9
# ╠═27e9d8ea-f826-11ea-3cbe-31d9a7852345
# ╟─27f78c12-f826-11ea-1cd1-eda91d47196f
# ╟─280705b8-f826-11ea-1d51-251fdd5eff2e
# ╠═dc3a5122-f82c-11ea-1eff-c745d65ab11b
# ╟─dc3af532-f82c-11ea-3212-f1b3c852513b
# ╠═044a19a2-f866-11ea-2b89-a51866d89a50
# ╟─dc476376-f82c-11ea-10ac-97fcb8c78627
# ╠═25318e22-f866-11ea-015b-d736c83ebfaa
# ╟─2814f2fc-f826-11ea-3fbc-0541fe904b97
