
using Pkg, DrWatson, Random, Test

Random.seed!(123)

@quickactivate "SR2TuringPluto"
using Turing 
using StatisticalRethinking

#=
using Turing, DataFrames
=#

Turing.setprogress!(false);

@model function linear_reg(x, y)
    β ~ Normal(0, 1)
    α ~ Normal(0, 5)
    σ ~ Exponential(1)

    for i ∈ eachindex(y)
        y[i] ~ Normal(α + β * x[i], σ)
    end
    return (α, β, σ)
end

f(x) = 5 + 2 * x + 0.5 * randn()

Δ = 1.0; xs_train = 0:Δ:20; ys_train = f.(xs_train);
xs_test = [10 + i * Δ for i in 1:4]; ys_test = f.(xs_test);

# Infer
m_train = linear_reg(xs_train, ys_train);
#precis(m_train) |> display  # From StatisticalRethinking.precis

chain_lin_reg = sample(m_train, NUTS(100, 0.65), 200);
chain_lin_reg |> display

#m_test = linear_reg(xs_test, fill(missing, length(ys_test)));
m_test = linear_reg(xs_test, similar(xs_test, Missing));
#m_test = linear_reg(xs_test, missing);

# Use the new predict function!
predictions = predict(m_test, chain_lin_reg)
predictions |> display

# Get the mean predicted values.
ys_mean_pred = vec(mean(Array(group(predictions, :y)); dims = 1))
ys_std_pred = vec(std(Array(group(predictions, :y)); dims = 1))


# Get the prediction error:
ys_test - ys_mean_pred |> display
ys_test - ys_mean_pred2 |> display

#@test sum(abs2, ys_test - ys_pred) ≤ σ
