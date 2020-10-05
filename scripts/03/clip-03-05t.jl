
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using StanSample
	using StatisticalRethinking
end

md"## Clip-03-05t.jl"

md"##### Define the Stan language model."

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

md"##### Define the Stanmodel and set the output format to :mcmcchains."

sm = SampleModel("m_05s", model_05s);

md"###### Use 16 observations."

begin
	N2 = 4^2
	d = Binomial(9, 0.66)
	n2 = Int.(9 * ones(Int, N2))
	k2 = rand(d, N2)
end

md"##### Input data for cmdstan."

m_05s_data = Dict("N" => length(n2), "n" => n2, "k" => k2);

md"##### Sample using cmdstan."

rc = stan_sample(sm, data=m_05s_data);

md"##### Retrieve samples as an MCMCChains.Chain object and as a Particles summary.."

if success(rc)
  chn = read_samples(sm; output_format=:mcmcchains)
  dict = read_samples(sm; output_format=:particles)
end;

md"##### Describe the chains."

chn

md"##### Plot the chains."

plot(chn)

md"##### Particles summary of the chains,"

dict

md"##### Notice in this example that the prior theta (thetaprior), the `unconditioned-on-the data theta`, shows a mean of 0.5 and a std of 0.29."

md"## End of clip-03-05t.jl"

