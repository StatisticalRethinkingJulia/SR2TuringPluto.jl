# # m2.1t.jl

# m0.1t.jl

using DrWatson
@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking

# Define the data

k = 6; n = 9;

# Define the model

@model m2_1t(n, k) = begin
  theta ~ Beta(1, 1) # prior
  k ~ Binomial(n, theta) # model
  return k, theta
end;

# Use Turing mcmc

chns2_1t = sample(m2_1t(n, k), NUTS(0.65), 1000)

# Look at the proper draws (in corrected chn2)

chns2_1t |> display

# Show the hpd region

hpd(chns2_1t, alpha=0.055) |> display

df=DataFrame(chns2_1t[[:theta]])
part0_1t = Particles(df[:, [:theta]])
part0_1t |> display

quap0_1t = quap(m2_1t(n, k))
quap0_1t.coef |> display
quap0_1t.vcov[1] |> sqrt |> display


# End of `02/m2.1t.jl`
