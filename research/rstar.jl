### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 8a4bdc28-f9b7-11ea-068e-eb70e542a539
using Pkg, DrWatson

# ╔═╡ 8a4c20a2-f9b7-11ea-2d67-a1e2531dcbd6
begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
	using StatsPlots
end

# ╔═╡ 2a9eafae-1246-11eb-0073-f7df8df64915
using XGBoost

# ╔═╡ 4eb6e25c-f9b7-11ea-3e43-05d4d6f124f5
md"## m0.2t.jl"

# ╔═╡ 9de23082-04c7-11eb-3b7c-393a11a68cf4
md"##### This notebook takes a look at the [rstar()](https://arxiv.org/pdf/2003.07900.pdf) diagnostic. Currently I can't get MLJ & MLJModels to work (serialization issue?)."

# ╔═╡ 8a4c985e-f9b7-11ea-22cb-717b854951d6
md"##### Define a simple Normal model with unknown mean and variance."

# ╔═╡ 8a53fb6a-f9b7-11ea-1a81-5b84e8ec7e71
@model gdemo(x, y) = begin
  s ~ InverseGamma(2, 3)
  m ~ Normal(0, sqrt(s))
  x ~ Normal(m, sqrt(s))
  y ~ Normal(m, sqrt(s))
end

# ╔═╡ 8a5473b8-f9b7-11ea-3b00-75cd16c17551
md"#####  Run sampler, collect results."

# ╔═╡ 8a5e1b36-f9b7-11ea-36ef-ebe92bc31e40
chns = mapreduce(c -> sample(gdemo(1.5, 2), NUTS(0.65), 2000), chainscat, 1:4)

# ╔═╡ 8a5ea05e-f9b7-11ea-2cbe-e15b867548b3
md"##### Plot the chains."

# ╔═╡ 8a6d9e26-f9b7-11ea-1d0f-f95263880beb
plot(chns; seriestype=:traceplot)

# ╔═╡ 45425d72-f9b8-11ea-1faa-41e93626b8b5
plot(chns; seriestype=:density)

# ╔═╡ 5001e4fc-1246-11eb-38f0-a5435c180d1a
"""
    rstar(chains::Chains; subset = 0.8, niter = 1_000, eta = 0.5, XGBoostParams)
    rstar(chains::Chains, iterations::Int; subset = 0.8, niter = 1_000, eta = 0.5, XGBoostParams)
    rstar(x::AbstractMatrix, y::AbstractVector, nchains::Int, iterations::Int; subset = 0.8, niter = 1_000, eta = 0.5, XGBoostParams)

Compute the R* statistic for convergence diagnostic of MCMC. This implementation is an adaption of Algorithm 1 & 2, described in [Lambert & Vehtari]. Note that the correctness of the statistic depends on the convergence of the classifier used internally in the statistic. You can track if the training of the classifier converged by inspection of the printed RMSE values from the XGBoost backend. To adjust the number of iterations used to train the classifier set `niter` accordingly.

# Usage

```julia-repl
using XGBoost
# You need to load XGBoost before using MCMCChains.Rstar

...

chn = ...

# Compute R⋆ using defaults settings for the  gradient boosting classifier used to compute the statistic.
# This is the recomended use.
R = rstar(chn)

# Compute 100 samples of the R⋆ statistic using sampling from according to the prediction probabilities.
# This approach can be slow and results in a less accurate estimation of the R* statistic.
# See discussion in Section 3.1.3 in the paper.
Rs = rstar(chn, 100)

# estimate Rstar
R = mean(Rs)

# visualize distribution
histogram(Rs)
```

## References:
[Lambert & Vehtari] Ben Lambert and Aki Vehtari. "R∗: A robust MCMC convergence diagnostic with uncertainty using gradient-boostined machines." Arxiv 2020.
"""
function rstar(x::AbstractMatrix, y::AbstractVector, nchains::Int, iterations::Int; subset = 0.8, niter = 1_000, eta = 0.5, xgboostparams...)

    N = length(y)

    # randomly sub-select training and testing set
    Ntrain = round(Int, N*subset)
    Ntest = N - Ntrain

    ids = Random.randperm(N)
    train_ids = view(ids, 1:Ntrain)
    test_ids = view(ids, (Ntrain+1):N)

    @assert length(test_ids) == Ntest

    # use predicted probabilities?
    mode = iterations > 1 ? "multi:softprob" : "multi:softmax"

    @info "Training classifier"
    # train classifier using XGBoost
    classif = XGBoost.xgboost(x[train_ids,:], niter; label = y[train_ids], 
                      objective = mode, num_class = nchains,
                      xgboostparams...)

    @info "Computing R* statistics"
    Rstats = ones(iterations) * Inf
    for i in 1:iterations
        # predict labels for "test" data
        p = XGBoost.predict(classif, x[test_ids,:])

        pred = if length(p) == Ntest*nchains
            probs = reshape(p, Ntest, nchains)
            map(s -> rand(Categorical(s / sum(s))), eachrow(probs))
        else
            p
        end

        # compute statistic
        a = mean(pred .== y[test_ids])
        Rstats[i] = nchains*a
    end

    return Rstats
end

# ╔═╡ 5002183c-1246-11eb-2867-adbcc9adab3c
function rstar(chn::Chains, iterations; kwargs...)
    nchains = size(chn, 3)
    @assert nchains > 1

    # collect data
    x = mapreduce(c -> Array(chn[:,:,c]), vcat, chains(chn))
    y = mapreduce(c -> ones(Int, length(chn))*c, vcat, chains(chn)) .- 1

    return rstar(x, y, nchains, iterations, kwargs...)
end

# ╔═╡ 5002de04-1246-11eb-0cfb-e931786d8c7f
rstar(chn::Chains; kwargs...) = first(rstar(chn, 1; kwargs...))

# ╔═╡ 4d6193a0-1246-11eb-2b25-1d555f611e5c
rstar(chns, 1)

# ╔═╡ 727fa3e8-1246-11eb-3057-6fe955204f01
rs = rstar(chns, 100)

# ╔═╡ 7280a8e2-1246-11eb-3cd5-8b6f9981eada
begin
	density(rs, lab="rs density")
	vline!([mean(rs)], lab="mean(rs)")
end

# ╔═╡ 8a6e37c6-f9b7-11ea-2abd-5db63b061a7d
md"## End m0.2t.jl"

# ╔═╡ Cell order:
# ╟─4eb6e25c-f9b7-11ea-3e43-05d4d6f124f5
# ╟─9de23082-04c7-11eb-3b7c-393a11a68cf4
# ╠═8a4bdc28-f9b7-11ea-068e-eb70e542a539
# ╠═8a4c20a2-f9b7-11ea-2d67-a1e2531dcbd6
# ╟─8a4c985e-f9b7-11ea-22cb-717b854951d6
# ╠═8a53fb6a-f9b7-11ea-1a81-5b84e8ec7e71
# ╟─8a5473b8-f9b7-11ea-3b00-75cd16c17551
# ╠═8a5e1b36-f9b7-11ea-36ef-ebe92bc31e40
# ╟─8a5ea05e-f9b7-11ea-2cbe-e15b867548b3
# ╠═8a6d9e26-f9b7-11ea-1d0f-f95263880beb
# ╠═45425d72-f9b8-11ea-1faa-41e93626b8b5
# ╠═2a9eafae-1246-11eb-0073-f7df8df64915
# ╠═5001e4fc-1246-11eb-38f0-a5435c180d1a
# ╠═5002183c-1246-11eb-2867-adbcc9adab3c
# ╠═5002de04-1246-11eb-0cfb-e931786d8c7f
# ╠═4d6193a0-1246-11eb-2b25-1d555f611e5c
# ╠═727fa3e8-1246-11eb-3057-6fe955204f01
# ╠═7280a8e2-1246-11eb-3cd5-8b6f9981eada
# ╟─8a6e37c6-f9b7-11ea-2abd-5db63b061a7d
