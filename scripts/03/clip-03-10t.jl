
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StanSample
	using StatisticalRethinking
end

md"## Clip-03-10t.jl"

md"##### Define the Stan language model."

m_17s = "
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

md"##### Define the SampleModel."

sm = SampleModel("m_17s", m_17s);

md"##### Use 4 observations."

begin
	N2 = 4
	n2 = Int.(9 * ones(Int, N2))
	k2 = [6, 5, 7, 6]
end

md"##### Input data for stan_sample()."

m_17s_data = Dict("N" => length(n2), "n" => n2, "k" => k2);

md"##### Sample using stan_sample()."

rc = stan_sample(sm, data=m_17s_data)

if success(rc)
  chn = read_samples(sm; output_format=:mcmcchains)
end

chn

md"##### Plot the chains."

begin
	mixeddensity(chn)
	bnds = MCMCChains.hpd(chn)
	vline!([bnds[:theta, :lower]], line=:dash)
	vline!([bnds[:theta, :upper]], line=:dash)
end

md"##### Look at area of hpd."

MCMCChains.hpd(chn)

md"## End of clip-03-10t.jl"

