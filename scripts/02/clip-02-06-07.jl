# Clip-02-06-07.jl

using DrWatson
@quickactivate "StatisticalRethinkingTuring"
using StatisticalRethinking
using Turing

include(srcdir("quap.jl"))

# snippet 2.3

p_grid = range(0, 1, length = 20)
prior = ones(20)
likelihood = pdf.(Binomial.(9, p_grid), 6)
posterior = likelihood .* prior
posterior ./= sum(posterior)

# %% 2.4
plot(p_grid, posterior, m = 3, legend = false,
    xlabel = "probability of water", ylabel = "posterior probability",
    title = "20 points")

# snippet 2.6

@model globethrowing(W, L) = begin
    p ~ Uniform(0, 1)
    W ~ Binomial(W + L, p)
end
m = globethrowing(6, 3)
r = quap(m)

d = MvNormal([r.coef.p], collect(reshape(r.vcov, 1, 1)))

s = rand(d, 10000)'

#histogram!(collect(s), normalize = :probability)
@show r.vcov[1] |> sqrt

plot!(p_grid, posterior, m = 3)

#histogram!(rand(d, 10000)', normalize = :probability)

# End clip-02-06-07.jl