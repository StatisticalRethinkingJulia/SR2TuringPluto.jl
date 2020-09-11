### A Pluto.jl notebook ###
# v0.11.10

using Markdown
using InteractiveUtils

# ╔═╡ ea02e53e-e0cc-11ea-2265-3165c74e19b3
using Pkg, DrWatson

# ╔═╡ 02ea1c3e-e0cd-11ea-04d1-f34150f81c89
begin
	@quickactivate "StatisticalRethinkingStan"
	using StatisticalRethinking
	using StanSample
end

# ╔═╡ d8ea5d90-e0cc-11ea-0d2d-25c807c1ae80
md"## # Fig 2.5s"

# ╔═╡ 11533bb6-e0cd-11ea-331e-278de5d6b26f
md"""
### This clip is only intended to generate Fig 2.5. 

It is not intended to show how to use Stan (yet)!"""

# ╔═╡ 3067f2d0-e0cd-11ea-17d4-ab40276a0379
m2_0 = "
// Inferring a rate
data {
  int n;
  int k;
}
parameters {
  real<lower=0,upper=1> theta;
}
model {
  // Prior distribution for θ
  theta ~ uniform(0, 1);

  // Observed Counts
  k ~ binomial(n, theta);
}
";

# ╔═╡ 40b6f53c-e0cd-11ea-298d-457305accab3
md"### 1. Create a SampleModel object:"

# ╔═╡ 5582ffb8-e0cd-11ea-358b-a1e5bac536af
m2_0s = SampleModel("m2.0s", m2_0);

# ╔═╡ 7d380960-e0cd-11ea-1401-736a8ff3f998
md"##### n will go from 1:9"

# ╔═╡ 89ee4714-e0cd-11ea-2378-dffe0e058ea1
begin
	k = [1,0,1,1,1,0,1,0,1]       # Sequence observed
	x = range(0, stop=9, length=10)
end

# ╔═╡ c8d96392-e0cc-11ea-2650-a950d62c37a5
begin
	p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 9)
	dens = Vector{DataFrame}(undef, 10)
	for n in 1:9

		p[n] = plot(xlims=(0.0, 1.0), ylims=(0.0, 3.0), leg=false)
		m2_0_data = Dict("n" => n, "k" => sum(k[1:n]));
		rc = stan_sample(m2_0s, data=m2_0_data);
		dfs = read_samples(m2_0s; output_format=:dataframe)
		if n == 1
			hline!([1.0], line=(:dash))
		else
			density!(dens[n][:, :theta], line=(:dash))
		end
		density!(dfs[:, :theta])
		dens[n+1] = dfs

	end
end

# ╔═╡ dd63b5f0-e0cd-11ea-0063-61ba73f99cac
plot(p..., layout=(3, 3))

# ╔═╡ ee6b3094-e0cd-11ea-1ceb-6f178f55cb23
md"## End of Fig2.5s.jl"

# ╔═╡ Cell order:
# ╟─d8ea5d90-e0cc-11ea-0d2d-25c807c1ae80
# ╠═ea02e53e-e0cc-11ea-2265-3165c74e19b3
# ╠═02ea1c3e-e0cd-11ea-04d1-f34150f81c89
# ╟─11533bb6-e0cd-11ea-331e-278de5d6b26f
# ╠═3067f2d0-e0cd-11ea-17d4-ab40276a0379
# ╟─40b6f53c-e0cd-11ea-298d-457305accab3
# ╠═5582ffb8-e0cd-11ea-358b-a1e5bac536af
# ╟─7d380960-e0cd-11ea-1401-736a8ff3f998
# ╠═89ee4714-e0cd-11ea-2378-dffe0e058ea1
# ╠═c8d96392-e0cc-11ea-2650-a950d62c37a5
# ╠═dd63b5f0-e0cd-11ea-0063-61ba73f99cac
# ╟─ee6b3094-e0cd-11ea-1ceb-6f178f55cb23
