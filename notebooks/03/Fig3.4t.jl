### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ f394ac3a-f787-11ea-314f-bf1bac91a4e9
using Pkg, DrWatson

# ╔═╡ f394e7e2-f787-11ea-1c4e-71f58d4abbea
begin
	@quickactivate "SR2TuringPluto"
	using StatisticalRethinking
end

# ╔═╡ 7c07959e-f787-11ea-14c6-23630e43ea67
md"## Fig3.4t.jl"

# ╔═╡ f3956618-f787-11ea-2032-fbf3b79ef810
p_grid = range(0, stop=1, length=10000)

# ╔═╡ 8171efda-f788-11ea-2c6b-834d66a19f38
prior = ones(length(p_grid));

# ╔═╡ 8172298c-f788-11ea-1cc5-4779c3e53066
likelihood = pdf.(Binomial.(3, p_grid), 3)

# ╔═╡ 8172c05c-f788-11ea-0fe9-d513d34fd0cd
begin
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
end

# ╔═╡ 8182be1c-f788-11ea-3fba-8966e136c7d4
samples = sample(p_grid, Weights(posterior), length(p_grid));

# ╔═╡ f3a1f566-f787-11ea-1901-21b24e86c6ba
begin
	p1 = density(samples;
	  xlab="Proportion water (p)", ylab="Density", lab="density",
	  title="Mean, median & mode", leg=:topleft);
	vline!(p1, [mode(samples)], lab="mode")
	vline!(p1, [median(samples)], lab="median")
	vline!(p1, [mean(samples)], lab="mean")
end;

# ╔═╡ f3a27664-f787-11ea-0a07-0f5212c9a2b7
begin
	loss = map(p -> sum(posterior .* abs.(p .- p_grid)), p_grid)
	m = findmin(loss)
	p2 = plot(loss;
	  xlab="Decision (p)", ylab="Expected proportional loss",
	  title="Loss value", lab="Loss function");
	scatter!([m[2]], [m[1]], lab="p_grid[$(m[2])]=$(round(m[1], digits=3))")
end;

# ╔═╡ f3b12006-f787-11ea-3f0c-27f71c4ec2da
plot(p1, p2, layout=(1, 2))

# ╔═╡ f3b1aed6-f787-11ea-30cf-515030a997c3
md"## End of Fig3.4t.jl"

# ╔═╡ Cell order:
# ╟─7c07959e-f787-11ea-14c6-23630e43ea67
# ╠═f394ac3a-f787-11ea-314f-bf1bac91a4e9
# ╠═f394e7e2-f787-11ea-1c4e-71f58d4abbea
# ╠═f3956618-f787-11ea-2032-fbf3b79ef810
# ╠═8171efda-f788-11ea-2c6b-834d66a19f38
# ╠═8172298c-f788-11ea-1cc5-4779c3e53066
# ╠═8172c05c-f788-11ea-0fe9-d513d34fd0cd
# ╠═8182be1c-f788-11ea-3fba-8966e136c7d4
# ╠═f3a1f566-f787-11ea-1901-21b24e86c6ba
# ╠═f3a27664-f787-11ea-0a07-0f5212c9a2b7
# ╠═f3b12006-f787-11ea-3f0c-27f71c4ec2da
# ╟─f3b1aed6-f787-11ea-30cf-515030a997c3
