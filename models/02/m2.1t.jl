## m2.1t.jl

using DrWatson
@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking

# Define the data

k = 6; n = 9;

# Define the model

@model function ppl2_1(n, k)
  theta ~ Beta(1, 1) # prior
  k ~ Binomial(n, theta) # model
  return k, theta
end;

# Use Turing mcmc

m2_1t = ppl2_1(n, k)
chns2_1t = sample(m2_1t, NUTS(0.65), 1000)

# Look at the proper draws (in corrected chn2)

chns2_1t |> display

# Show the hpd region

hpd(chns2_1t, alpha=0.055) |> display

Particles(chns2_1t[:theta]) |> display

q2_1t = quap(m2_1t)
q2_1t.coef |> display
q2_1t.vcov .|> sqrt |> display


# End of m2.1t.jl
