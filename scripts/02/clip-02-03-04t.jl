# Clip-02-03-06.jl

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

#=
This is another difference between Julia and R. You don't have automatic
broadcasting (making functions work the same whether you run it on a single
number or on an array of numbers). However, fixing this is very easy,
you only have to annotate your function call with a dot:
  `f(single_number)` -> `f.(array_of_numbers)`.

That means the translation of the line
  `likelihood = pdf.(Binomial.(9, p_grid), 6)`
goes something like: first, for every value in `p_grid` make a binomial
distribution with that p value. For every distribution then take the pdf
of 6. The same in the next line, this is elementwise multiplication.

Related advice: when you get an error of the kind "no method matching 
..(::Float64, ::Array{Float64})" or similar, first check if all your
dots are in order.

For more infos read https://julialang.org/blog/2017/01/moredots/
=#

# %% 2.4
plot(p_grid, posterior, m = 3, legend = false,
    xlabel = "probability of water", ylabel = "posterior probability", title = "20 points")

# %% 2.5
prior = ifelse.(p_grid .< 0.5, 0, 1)
prior = @. exp(-5 * abs(p_grid - 0.5))  # @. means broadcast everything that follows

# %% 2.6
@model globethrowing(W, L) = begin
    p ~ Uniform(0, 1)
    W ~ Binomial(W + L, p)
end
m = globethrowing(6, 3)
r = quap(m)

d = MvNormal([r.coef.p], collect(reshape(r.vcov, 1, 1)))

s = rand(d, 10000)'

histogram!(collect(s), normalize = :probability)
@show r.vcov[1] |> sqrt

plot!(p_grid, posterior ./ 2.5, m = 3)

histogram!(rand(d, 10000)', normalize = :probability)

# End clip-02-03-06.jl