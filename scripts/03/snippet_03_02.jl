using DrWatson
@quickactivate "StatReth"

# %%
using StatsBase
using Distributions
using StatsPlots
using Turing  # for MCMCChains

# %% 3.2
n = 1000
p_grid = range(0, 1, length = n)
prob_p = ones(n)
prob_data = @. pdf(Binomial(9, p_grid), 6)
posterior = prob_data .* prob_p
posterior ./= sum(posterior)

plot(p_grid, posterior)

# %% 3.3 - 3.5
weights = pweights(posterior)
samples = sample(p_grid, weights, 10_000)

scatter(samples)
density(samples)
# plot!(p_grid, posterior * n)

# %% 3.6, 3.7
sum(posterior[p_grid .< 0.5])

sum(samples .< 0.5) / 10_000    # either
mean(samples .< 0.5)            # or

# %% 3.8
mean(0.5 .< samples .< 0.75)    # either
mean(samples) do s              # or (this will be faster than above)
    0.5 < s < 0.75
end
mean(samples) do s              # or (can do multi-line; return variable on the last line)
    condition1 = 0.5 < s
    condition2 = s < 0.75
    both_conditions = condition1 & condition2
    both_conditions
end

# %% 3.9, 3.10
quantile(samples, 0.8)
quantile(samples, (0.1, 0.9))

# %% 3.11
p_grid = range(0, 1, length = n)
prior = ones(n)
likelihood = @. pdf(Binomial(3, p_grid), 3)
posterior = likelihood .* prior
posterior ./= sum(posterior)
samples = sample(p_grid, pweights(posterior), 10_000)

# %% 3.12, 3.13
quantile(samples, (0.25, 0.75))

chn = Chains(samples)
hpd(chn; alpha = 0.5)

# %% 3.14, 3.16, TODO: not sure about 3.15
p_grid[argmax(posterior)]

mean(samples)
median(samples)

# %% 3.17
sum(posterior .* abs.(0.5 .- p_grid))

# %% 3.18, 3.19
loss = map(p_grid) do d
    sum(posterior .* abs.(d .- p_grid))
end
plot(loss)

p_grid[argmin(loss)]
