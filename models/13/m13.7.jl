# ### m13.7t.jl

using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"
using Turing
using StatisticalRethinking
Turing.setprogress!(false)

df = CSV.read(sr_datadir("Kline2.csv"), DataFrame);
Dmat = CSV.read(sr_datadir("islandDistMatrix.csv"), DataFrame)

df.society = 1:10

@model ppl13_7(Dmat, society, logpop, total_tools) = begin
    rhosq ~ truncated(Cauchy(0, 1), 0, Inf)
    etasq ~ truncated(Cauchy(0, 1), 0, Inf)
    bp ~ Normal(0, 1)
    a ~ Normal(0, 10)
    
    # GPL2
    SIGMA_Dmat = etasq * exp.(-rhosq * Dmat.^2)
    SIGMA_Dmat = SIGMA_Dmat + 0.01I
    SIGMA_Dmat = (SIGMA_Dmat' + SIGMA_Dmat) / 2
    g ~ MvNormal(zeros(10), SIGMA_Dmat)
    
    log_lambda = a .+ g[society] .+ bp * logpop
    
    total_tools .~ Poisson.(exp.(log_lambda))
end

m13_7t = ppl13_7(Dmat, df.society, df.logpop, df.total_tools)
nchains = 4; sampler = NUTS(0.65); nsamples=2000
#nchains = 1 # Fails with nchains=4
chns13_7t = mapreduce(c -> sample(m13_7t, sampler, nsamples), chainscat, 1:nchains)

m_13_7s_results = """
            mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
    g[1]   -0.27    0.01  0.44  -1.26  -0.50  -0.24  -0.01   0.53  3278    1
    g[2]   -0.13    0.01  0.43  -1.07  -0.34  -0.10   0.12   0.67  3144    1
    g[3]   -0.17    0.01  0.42  -1.10  -0.37  -0.14   0.07   0.59  3096    1
    g[4]    0.30    0.01  0.37  -0.51   0.12   0.30   0.50   1.02  3145    1
    g[5]    0.03    0.01  0.37  -0.77  -0.15   0.04   0.22   0.72  3055    1
    g[6]   -0.46    0.01  0.38  -1.31  -0.64  -0.42  -0.23   0.19  3306    1
    g[7]    0.10    0.01  0.36  -0.70  -0.07   0.11   0.29   0.78  3175    1
    g[8]   -0.26    0.01  0.37  -1.07  -0.43  -0.24  -0.06   0.40  3246    1
    g[9]    0.23    0.01  0.35  -0.52   0.07   0.25   0.42   0.89  3302    1
    g[10]  -0.12    0.01  0.45  -1.04  -0.36  -0.11   0.13   0.77  5359    1
    a       1.31    0.02  1.15  -0.98   0.61   1.31   2.00   3.69  4389    1
    bp      0.25    0.00  0.11   0.02   0.18   0.24   0.31   0.47  5634    1
    etasq   0.34    0.01  0.53   0.03   0.10   0.20   0.39   1.52  5643    1
    rhosq   1.52    0.09 11.82   0.03   0.16   0.39   0.96   7.96 15955    1
    lp__  925.98    0.03  2.96 919.16 924.20 926.34 928.14 930.67  7296    1    
"""

m_13_7_rethinking = """
       mean   sd  5.5% 94.5% n_eff Rhat4
k[1]  -0.16 0.34 -0.71  0.33   615  1.01
k[2]  -0.02 0.33 -0.55  0.47   550  1.01
k[3]  -0.07 0.32 -0.58  0.40   566  1.01
k[4]   0.35 0.29 -0.05  0.81   557  1.01
k[5]   0.08 0.29 -0.36  0.51   509  1.01
k[6]  -0.38 0.30 -0.88  0.03   585  1.00
k[7]   0.14 0.28 -0.29  0.56   567  1.00
k[8]  -0.21 0.29 -0.66  0.20   526  1.01
k[9]   0.26 0.27 -0.15  0.65   534  1.01
k[10] -0.17 0.36 -0.73  0.36   768  1.00
g      0.61 0.60  0.07  1.70  1321  1.00
b      0.28 0.09  0.14  0.42  1011  1.00
a      1.38 1.05  0.25  3.30  1800  1.00
etasq  0.21 0.22  0.03  0.61   982  1.00
rhosq  1.28 1.56  0.08  4.32  2020  1.00
"""

# End of m13.7t.jl
