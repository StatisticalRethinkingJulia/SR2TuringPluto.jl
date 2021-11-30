### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 175c73f6-f2e5-11ea-3458-eddf1465c5c4
using Pkg, DrWatson

# ╔═╡ 175d464e-f2e5-11ea-318a-dd4ad9573a95
using StatisticalRethinking

# ╔═╡ 0b40c1ca-f2e4-11ea-1658-0d38b13854c0
md"## Clip-03-02-05t.jl"

# ╔═╡ 175ca6a0-f2e5-11ea-066c-49e37107aea2
@quickactivate "SR2TuringPluto"

# ╔═╡ 176c899a-f2e5-11ea-3f17-59b92f0a2483
md"### snippet 3.2"

# ╔═╡ 176d23ca-f2e5-11ea-3a5c-93353273f6eb
begin
	p_grid = range(0, stop=1, length=1000)
	prior = ones(length(p_grid))
	likelihood = pdf.(Binomial.(9, p_grid), 6)
	posterior = likelihood .* prior
	posterior = posterior / sum(posterior)
end

# ╔═╡ 177885bc-f2e5-11ea-0d11-f9ae48c08d8f
md"### snippet 3.3"

# ╔═╡ 17798f98-f2e5-11ea-299e-f3424848225c
md"##### Draw 10000 samples from this posterior distribution."

# ╔═╡ 1784a202-f2e5-11ea-051b-09f8abb69875
begin
	N = 10000
	samples = sample(p_grid, Weights(posterior), N)
end;

# ╔═╡ 178bd8ba-f2e5-11ea-21dd-61c16076a938
md"##### Create an MCMCChains.Chains object."

# ╔═╡ 178d0d36-f2e5-11ea-2593-bf4c46edb720
chn = MCMCChains.Chains(reshape(samples, N, 1, 1), [:p]);

# ╔═╡ df229f44-7dd0-11eb-2ebe-9b12eff28aa1
typeof(chn)

# ╔═╡ 179c95ce-f2e5-11ea-2358-714343a64afa
md"##### Describe the chain."

# ╔═╡ 17a5374c-f2e5-11ea-1bfe-3fb1cce48d48
Text(sprint(show, "text/plain", summarize(chn)))

# ╔═╡ 17a70c16-f2e5-11ea-2e3f-2dc8006c1ee7
md"##### Plot the chain."

# ╔═╡ 17aece92-f2e5-11ea-311a-1b4591b1172a
p1 = plot(chn; seriestype=:traceplot)

# ╔═╡ 17bdf552-f2e5-11ea-28e5-1f5f22b43146
md"### snippet 3.4"

# ╔═╡ 17c6af64-f2e5-11ea-0b28-c55bf548b3ea
p2 = plot(chn; seriestype=:density)

# ╔═╡ 17cf1094-f2e5-11ea-1a80-858534e81348
md"### snippet 3.5"

# ╔═╡ 17d80aa0-f2e5-11ea-1c8f-01770152e142
md"##### Compare with analytical (conjugate) solution."

# ╔═╡ 17e08cca-f2e5-11ea-31a6-8931f98059ba
begin
	w = 6
	n = 9
	x = 0:0.01:1
	plot( x, pdf.(Beta( w+1 , n-w+1 ) , x ), lab="Conjugate solution")
	density!(samples, lab="Sample density")
end

# ╔═╡ 17f4a430-f2e5-11ea-048a-73ae33e192ed
begin
	density(samples, lab="Sample2 density")
	vline!(hpdi(samples), lab="hpdi samples2")
	vline!(quantile(samples, [0.25, 0.75]), lab="quantiles [0.25, 0.75]")
end

# ╔═╡ 17fecdb6-f2e5-11ea-2c5b-7bc56e8b1360
md"## End of clip-03-02-05t.jl"

# ╔═╡ Cell order:
# ╟─0b40c1ca-f2e4-11ea-1658-0d38b13854c0
# ╠═175c73f6-f2e5-11ea-3458-eddf1465c5c4
# ╠═175ca6a0-f2e5-11ea-066c-49e37107aea2
# ╠═175d464e-f2e5-11ea-318a-dd4ad9573a95
# ╟─176c899a-f2e5-11ea-3f17-59b92f0a2483
# ╠═176d23ca-f2e5-11ea-3a5c-93353273f6eb
# ╟─177885bc-f2e5-11ea-0d11-f9ae48c08d8f
# ╟─17798f98-f2e5-11ea-299e-f3424848225c
# ╠═1784a202-f2e5-11ea-051b-09f8abb69875
# ╟─178bd8ba-f2e5-11ea-21dd-61c16076a938
# ╠═178d0d36-f2e5-11ea-2593-bf4c46edb720
# ╠═df229f44-7dd0-11eb-2ebe-9b12eff28aa1
# ╟─179c95ce-f2e5-11ea-2358-714343a64afa
# ╠═17a5374c-f2e5-11ea-1bfe-3fb1cce48d48
# ╟─17a70c16-f2e5-11ea-2e3f-2dc8006c1ee7
# ╠═17aece92-f2e5-11ea-311a-1b4591b1172a
# ╟─17bdf552-f2e5-11ea-28e5-1f5f22b43146
# ╠═17c6af64-f2e5-11ea-0b28-c55bf548b3ea
# ╟─17cf1094-f2e5-11ea-1a80-858534e81348
# ╟─17d80aa0-f2e5-11ea-1c8f-01770152e142
# ╠═17e08cca-f2e5-11ea-31a6-8931f98059ba
# ╠═17f4a430-f2e5-11ea-048a-73ae33e192ed
# ╟─17fecdb6-f2e5-11ea-2c5b-7bc56e8b1360
