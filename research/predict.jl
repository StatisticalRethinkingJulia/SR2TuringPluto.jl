using Pkg, DrWatson, Test

@quickactivate "StatisticalRethinkingTuring"
using Turing 
using StatisticalRethinking
using Turing; Turing.turnprogress(false);

@model function linear_reg(x, y, σ = 0.1)
    β ~ Normal(0, 1)

    for i ∈ eachindex(y)
        y[i] ~ Normal(β * x[i], σ)
    end
end

σ = 0.1; f(x) = 2 * x + 0.1 * randn()

Δ = 0.1; xs_train = 0:Δ:10; ys_train = f.(xs_train);
xs_test = [10 + Δ, 10 + 2 * Δ]; ys_test = f.(xs_test);

# Infer
m_train = linear_reg(xs_train, ys_train, σ);
chain_lin_reg = sample(m_train, NUTS(100, 0.65), 200);
chain_lin_reg |> display

m_test = linear_reg(
    xs_test, 
    Vector{Union{Missing, Float64}}(undef, length(ys_test)), 
    σ
);

# Use the new predict function!
predictions = predict(m_test, chain_lin_reg)

# Get the mean predicted values.
ys_pred = vec(mean(Array(group(predictions, :y)); dims = 1))

# Get the prediction error:
ys_test - ys_pred |> display

@test sum(abs2, ys_test - ys_pred) ≤ σ
