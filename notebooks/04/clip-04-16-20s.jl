### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ f6b0972a-f8e8-11ea-335d-09e48cf0ff26
using Pkg, DrWatson

# ╔═╡ f6b0ee5a-f8e8-11ea-3427-c5eda9c27be2
begin
	@quickactivate "StatisticalRethinkingTuring"
	using StanSample
	using StatisticalRethinking
end

# ╔═╡ eb0545fc-f8e7-11ea-0c5f-d93cf6088ae0
md"## Clip-04-16-22s.jl"

# ╔═╡ f6b189fa-f8e8-11ea-3c7e-db6b39ae5110
df = CSV.read(sr_datadir("Howell1.csv"), DataFrame; delim=';');

# ╔═╡ f6c1db52-f8e8-11ea-139b-cb0a529c7a67
md"### snippet 4.8"

# ╔═╡ f6c28070-f8e8-11ea-148a-6364b90cc4fe
md"##### Use only adults."

# ╔═╡ f6cec6c0-f8e8-11ea-3b68-1dbc8afb5206
df2 = filter(row -> row[:age] >= 18, df);

# ╔═╡ f6cfe79c-f8e8-11ea-2c05-c1cbd48b9c72
md"##### Show first 5 rows of DataFrame df."

# ╔═╡ f6db9e66-f8e8-11ea-280a-93d48faed703
first(df2, 5)

# ╔═╡ f6dceb68-f8e8-11ea-3258-f3c3b6f4ef76
md"### Snippet 4.16"

# ╔═╡ f6ea346c-f8e8-11ea-375b-19bd5ff09fb7
md"##### Generate approximate probabilities."

# ╔═╡ f6eae056-f8e8-11ea-27d8-9b4dca29c5ef
function grid_prob(x, y, prior_x, prior_y, obs)

	# Create an x vs. y grid (vector of vectors), e.g.
	# 10000-element Array{Array{Float64,1},1}:
 	#	[150.0, 7.0]
 	#	[150.1010101010101, 7.0]
 	#	[150.2020202020202, 7.0]
 	#   ...

 	df = DataFrame()
	grid = reshape([ [x,y]  for x=x, y=y ], length(x)*length(y))

	# Define the priors

	d2 = Normal(178.0, 20.0)
	d3 = Uniform(0, 50)

	# Compute the log(likelihood * prior)

	the_prod = []
	for i in 1:length(grid)
	    d1 = Normal(grid[i][1], grid[i][2])
	    ll = sum(log.(pdf.(d1, obs)))
	    append!(df, DataFrame(mu=grid[i][1], sigma=grid[i][2],
	    	ll=ll))
		append!(the_prod, ll + log.(pdf.(prior_x, grid[i][1])) + 
			log.(pdf.(prior_y, grid[i][2])))
	end

	# Make it a probability

	df[!, :prob] = exp.(the_prod .- maximum(the_prod))
	df
end

# ╔═╡ f6f2cb18-f8e8-11ea-1e7c-0bcb5a6cf62e
begin
	mu_list = range(150, 160, length=100)
	sigma_list = range(7, 9, length=100)
	prior_mu = Normal(178.0, 20.0)
	prior_sigma = Uniform(0, 50)
end

# ╔═╡ f6fc86b2-f8e8-11ea-3155-b13066e21ab6
md"### snippet 4.17"

# ╔═╡ f704264e-f8e8-11ea-24ad-2f534da03d56
post_df = grid_prob(mu_list, sigma_list, prior_mu, prior_sigma,
	df2[:, :height]);

# ╔═╡ f70b8194-f8e8-11ea-1418-698e93000714
p1 = contour(mu_list, sigma_list, post_df[:, :prob],
	xlim = (153.5, 155.7),
	ylim = (7.0, 8.5),
	xlab="height",
	ylab="sigma",
	title="Contour")

# ╔═╡ f713c0ca-f8e8-11ea-3fb6-95068a226a6d
md"### snippet 4.18"

# ╔═╡ f720c22a-f8e8-11ea-34ca-3ffd4eb34e01
p2 = heatmap(mu_list, sigma_list, transpose(reshape(post_df[:, :prob], 100,100)),
	xlim = (153.5, 155.7),
	ylim = (7.0, 8.5),
	xlab="height",
	ylab="sigma",
	title="Heatmap")

# ╔═╡ f7240cf0-f8e8-11ea-224a-89543e618960
md"### Snippet 4.19"

# ╔═╡ f72d36f4-f8e8-11ea-3659-3f9ce23635d2
md"##### Sample post_df."

# ╔═╡ f7359696-f8e8-11ea-04b4-add413a60b41
samples = post_df[sample(1:size(post_df, 1), Weights(post_df[:, :prob]), 
	10000, replace=true), :];

# ╔═╡ f73e1758-f8e8-11ea-1812-8d4723d20456
md"### Snippet 4.22"

# ╔═╡ f74778c0-f8e8-11ea-31b2-2b4adb0f865b
md"##### Convert to an MCMCChains.Chains object."

# ╔═╡ f751b452-f8e8-11ea-33a1-23bf97c35a18
begin
	a2d = hcat(samples[:, :mu], samples[:, :sigma])
	a3d = reshape(a2d, (size(a2d, 1), size(a2d, 2), 1))
	chn = StanSample.convert_a3d(a3d, ["mu", "sigma"], Val(:mcmcchains); start=1)
end

# ╔═╡ f75b7bf6-f8e8-11ea-3702-7df2f3263677
md"##### Show hpd regions."

# ╔═╡ f7661370-f8e8-11ea-0af3-97ae384db28c
bnds = MCMCChains.hpd(chn)

# ╔═╡ f76fece2-f8e8-11ea-1a99-7148dc026dad
md"### Snippet 4.21"

# ╔═╡ f779f750-f8e8-11ea-2d2a-2192846ac6d8
md"##### Density of mu."

# ╔═╡ f7844f52-f8e8-11ea-3b84-4d57615c237b
begin
	p3 = density(samples[:, :mu],
		xlab="height",
		ylab="density",
		lab="mu",
		title="posterior mu")
	vline!(p3, [bnds[:mu, :upper]], line=:dash, lab="Lower bound")
	vline!(p3, [bnds[:mu, :lower]], line=:dash, lab="Upper bound")
end;

# ╔═╡ f78e0d58-f8e8-11ea-1a60-87fc0b4e2625
md"##### Density of sigma."

# ╔═╡ f7982d94-f8e8-11ea-3b9d-9d9542626cc1
begin
	p4 = density(samples[:, :sigma],
		xlab="sigma",
		ylab="density",
		lab="sigma",
		title="posterior sigma")
	vline!(p4, [bnds[:sigma, :upper]], line=:dash, lab="Lower bound")
	vline!(p4, [bnds[:sigma, :lower]], line=:dash, lab="Upper bound")
end;

# ╔═╡ f7a39178-f8e8-11ea-06f9-13204fc0c56f
plot(p1, p2, p3, p4, layout=(2,2))

# ╔═╡ f7adf87a-f8e8-11ea-36e1-d95fd1feb4c6
md"## End of clip-04-16-20s.jl"

# ╔═╡ Cell order:
# ╟─eb0545fc-f8e7-11ea-0c5f-d93cf6088ae0
# ╠═f6b0972a-f8e8-11ea-335d-09e48cf0ff26
# ╠═f6b0ee5a-f8e8-11ea-3427-c5eda9c27be2
# ╠═f6b189fa-f8e8-11ea-3c7e-db6b39ae5110
# ╟─f6c1db52-f8e8-11ea-139b-cb0a529c7a67
# ╟─f6c28070-f8e8-11ea-148a-6364b90cc4fe
# ╠═f6cec6c0-f8e8-11ea-3b68-1dbc8afb5206
# ╟─f6cfe79c-f8e8-11ea-2c05-c1cbd48b9c72
# ╠═f6db9e66-f8e8-11ea-280a-93d48faed703
# ╟─f6dceb68-f8e8-11ea-3258-f3c3b6f4ef76
# ╟─f6ea346c-f8e8-11ea-375b-19bd5ff09fb7
# ╠═f6eae056-f8e8-11ea-27d8-9b4dca29c5ef
# ╠═f6f2cb18-f8e8-11ea-1e7c-0bcb5a6cf62e
# ╟─f6fc86b2-f8e8-11ea-3155-b13066e21ab6
# ╠═f704264e-f8e8-11ea-24ad-2f534da03d56
# ╠═f70b8194-f8e8-11ea-1418-698e93000714
# ╠═f713c0ca-f8e8-11ea-3fb6-95068a226a6d
# ╠═f720c22a-f8e8-11ea-34ca-3ffd4eb34e01
# ╠═f7240cf0-f8e8-11ea-224a-89543e618960
# ╠═f72d36f4-f8e8-11ea-3659-3f9ce23635d2
# ╠═f7359696-f8e8-11ea-04b4-add413a60b41
# ╟─f73e1758-f8e8-11ea-1812-8d4723d20456
# ╠═f74778c0-f8e8-11ea-31b2-2b4adb0f865b
# ╠═f751b452-f8e8-11ea-33a1-23bf97c35a18
# ╟─f75b7bf6-f8e8-11ea-3702-7df2f3263677
# ╠═f7661370-f8e8-11ea-0af3-97ae384db28c
# ╟─f76fece2-f8e8-11ea-1a99-7148dc026dad
# ╟─f779f750-f8e8-11ea-2d2a-2192846ac6d8
# ╠═f7844f52-f8e8-11ea-3b84-4d57615c237b
# ╟─f78e0d58-f8e8-11ea-1a60-87fc0b4e2625
# ╠═f7982d94-f8e8-11ea-3b9d-9d9542626cc1
# ╠═f7a39178-f8e8-11ea-06f9-13204fc0c56f
# ╟─f7adf87a-f8e8-11ea-36e1-d95fd1feb4c6
