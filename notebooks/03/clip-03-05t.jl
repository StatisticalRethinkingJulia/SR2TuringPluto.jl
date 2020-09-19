### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 77ce1956-f365-11ea-23ec-bffcd0523615
using Pkg, DrWatson

# ╔═╡ 77ce5950-f365-11ea-209d-e946b7e20a26
begin
	@quickactivate "StatisticalRethinkingTuring"
	using StanSample
	using StatisticalRethinking
end

# ╔═╡ 0072bcc0-f362-11ea-18eb-8df094618d42
md"## Clip-03-05t.jl"

# ╔═╡ 77ceda44-f365-11ea-1dd7-e3c4f23098a2
md"##### Define the Stan language model."

# ╔═╡ 77e2a9b6-f365-11ea-3d05-4d23a5ce0eed
model_05s = "
// Inferring a Rate
data {
  int N;
  int<lower=0> k[N];
  int<lower=1> n[N];
}
parameters {
  real<lower=0,upper=1> theta;
  real<lower=0,upper=1> thetaprior;
}
model {
  // Prior Distribution for Rate Theta
  theta ~ beta(1, 1);
  thetaprior ~ beta(1, 1);

  // Observed Counts
  k ~ binomial(n, theta);
}
";

# ╔═╡ 77e32788-f365-11ea-2c71-b1eac7897987
md"##### Define the Stanmodel and set the output format to :mcmcchains."

# ╔═╡ 77ef18f4-f365-11ea-17f7-5b75918ed324
sm = SampleModel("m_05s", model_05s);

# ╔═╡ 77efb390-f365-11ea-16f0-23312c3610e3
md"###### Use 16 observations."

# ╔═╡ 77fb3738-f365-11ea-1887-f7cecce15c32
begin
	N2 = 4^2
	d = Binomial(9, 0.66)
	n2 = Int.(9 * ones(Int, N2))
	k2 = rand(d, N2)
end

# ╔═╡ 77fde430-f365-11ea-3a92-5f60e48ae4ba
md"##### Input data for cmdstan."

# ╔═╡ 780cb418-f365-11ea-3311-a34d10d82a1f
m_05s_data = Dict("N" => length(n2), "n" => n2, "k" => k2);

# ╔═╡ 780f9ebc-f365-11ea-224f-1d86e552d93d
md"##### Sample using cmdstan."

# ╔═╡ 7814420a-f365-11ea-12be-d308c29a7481
rc = stan_sample(sm, data=m_05s_data);

# ╔═╡ 7bbbbb1e-f369-11ea-02a4-8b0f27dc1044
md"##### Retrieve samples as an MCMCChains.Chain object and as a Particles summary.."

# ╔═╡ add0c3ac-f368-11ea-1e7d-650962e47ef9
if success(rc)
  chn = read_samples(sm; output_format=:mcmcchains)
  dict = read_samples(sm; output_format=:particles)
end;

# ╔═╡ 90b9a838-f368-11ea-11d9-61225ef012fc
md"##### Describe the chains."

# ╔═╡ 9e931692-f368-11ea-2359-4b617b0cd65c
chn

# ╔═╡ 781b804c-f365-11ea-21d8-25088dbbbceb
md"##### Plot the chains."

# ╔═╡ bdb7dd78-f368-11ea-2039-2dbfa10cc364
plot(chn)

# ╔═╡ 78275f52-f365-11ea-24fe-132bee9323d0
md"##### Particles summary of the chains,"

# ╔═╡ 0d24edbc-f369-11ea-0ebc-65b4a405820e
dict

# ╔═╡ 0d25985a-f369-11ea-1cea-8137bcb666dc
md"##### Notice in this example that the prior theta (thetaprior), the `unconditioned-on-the data theta`, shows a mean of 0.5 and a std of 0.29."

# ╔═╡ 782867ee-f365-11ea-0f34-6b27f578bb94
md"## End of clip-03-05t.jl"

# ╔═╡ Cell order:
# ╟─0072bcc0-f362-11ea-18eb-8df094618d42
# ╠═77ce1956-f365-11ea-23ec-bffcd0523615
# ╠═77ce5950-f365-11ea-209d-e946b7e20a26
# ╟─77ceda44-f365-11ea-1dd7-e3c4f23098a2
# ╠═77e2a9b6-f365-11ea-3d05-4d23a5ce0eed
# ╟─77e32788-f365-11ea-2c71-b1eac7897987
# ╠═77ef18f4-f365-11ea-17f7-5b75918ed324
# ╟─77efb390-f365-11ea-16f0-23312c3610e3
# ╠═77fb3738-f365-11ea-1887-f7cecce15c32
# ╟─77fde430-f365-11ea-3a92-5f60e48ae4ba
# ╠═780cb418-f365-11ea-3311-a34d10d82a1f
# ╟─780f9ebc-f365-11ea-224f-1d86e552d93d
# ╠═7814420a-f365-11ea-12be-d308c29a7481
# ╟─7bbbbb1e-f369-11ea-02a4-8b0f27dc1044
# ╠═add0c3ac-f368-11ea-1e7d-650962e47ef9
# ╟─90b9a838-f368-11ea-11d9-61225ef012fc
# ╠═9e931692-f368-11ea-2359-4b617b0cd65c
# ╟─781b804c-f365-11ea-21d8-25088dbbbceb
# ╠═bdb7dd78-f368-11ea-2039-2dbfa10cc364
# ╠═78275f52-f365-11ea-24fe-132bee9323d0
# ╠═0d24edbc-f369-11ea-0ebc-65b4a405820e
# ╟─0d25985a-f369-11ea-1cea-8137bcb666dc
# ╟─782867ee-f365-11ea-0f34-6b27f578bb94
