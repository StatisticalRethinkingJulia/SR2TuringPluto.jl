### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ d9015225-b635-42ad-ae96-fe7772a031b7
using Pkg, DrWatson

# ╔═╡ 727e041b-fe08-4f7d-ad19-e3c57772c967
begin
	using Optim
	using GLM
	using CSV
	using Random
	using StatsBase
	using DataFrames
	using Dagitty
	using Turing
	using StatsPlots
	using StatisticalRethinking
	using StatisticalRethinkingPlots
	using Logging
end

# ╔═╡ 078fdb65-5bd6-4b0d-b672-617f5912480d
begin
	default(labels=false)
	Logging.disable_logging(Logging.Warn);
end;

# ╔═╡ 5141dfdf-490d-479e-ac58-63c929ae8fd1
md"## 7.1 The problem with parameters."

# ╔═╡ 611c3828-f01f-4721-830b-736b28620e7e
md"### Code 7.1"

# ╔═╡ ad511be2-fd88-4aba-a830-21e242d994f1
begin
	sppnames = ["afarensis", "africanus", "habilis", "boisei",
		"rudolfensis", "ergaster", "sapiens"]
	brainvolcc = [438, 452, 612, 521, 752, 871, 1350]
	masskg = [37.0, 35.5, 34.5, 41.5, 55.5, 61.0, 53.5]
	d = DataFrame(:species => sppnames, :brain => brainvolcc, :mass => masskg)
end;

# ╔═╡ 4ee38b80-aded-40af-93ac-11a4ac2b9378
md"### Code 7.2"

# ╔═╡ 5c7f1f26-02ac-4325-a06c-8c26c75e6cf4
begin
	d[!,:mass_std] = (d.mass .- mean(d.mass))./std(d.mass)
	d[!,:brain_std] = d.brain ./ maximum(d.brain)
	d
end

# ╔═╡ 5a2f79ca-2d3e-4cd4-9de9-7c7e283a5052
md"### Code 7.3"

# ╔═╡ 94b5539b-4374-40b8-9501-8d3df49f20c3
@model function model_m7_1(mass_std, brain_std)
    a ~ Normal(0.5, 1)
    b ~ Normal(0, 10)
    μ = @. a + b*mass_std
    log_σ ~ Normal()
    brain_std ~ MvNormal(μ, exp(log_σ))
end

# ╔═╡ 2aafe78c-b1b1-4963-ad89-dbe09b546d96
begin
	m7_1_ch = sample(model_m7_1(d.mass_std, d.brain_std), NUTS(), 1000)
	m7_1 = DataFrame(m7_1_ch)
	describe(m7_1)
end

# ╔═╡ 45122674-c5d2-4f61-af03-c74482618e90
md"### Code 7.4"

# ╔═╡ 3b7ca48f-078f-4eef-9d0c-fa8724e711d3
X = hcat(ones(length(d.mass_std)), d.mass_std)

# ╔═╡ 40c173f3-3dbb-4907-9188-632d34b3d7cf
m = lm(X, d.brain_std)

# ╔═╡ 39dc811c-2742-4977-a203-c9fdcb5b218d
md"### Code 7.5"

# ╔═╡ 85bb82da-efe7-4b34-8af5-cf60ae6e4757
Random.seed!(12);

# ╔═╡ 53806285-ea68-4f5c-b3fb-a382a894ec56
md"## Do explicit simulation due to log_σ."

# ╔═╡ 59b9a24b-6967-418b-a4f8-e67ceb0a5474
begin
	s = [
		rand(MvNormal((@. r.a + r.b * d.mass_std), exp(r.log_σ)))
		for r ∈ eachrow(m7_1)
	]
	s = vcat(s'...);

	r = mean.(eachcol(s)) .- d.brain_std;
	resid_var = var(r, corrected=false)
	outcome_var = var(d.brain_std, corrected=false)
	1 - resid_var/outcome_var
end

# ╔═╡ ed541621-2245-49f9-9672-31e20b9b7765
md"### Code 7.6"

# ╔═╡ 47d9566c-93ec-4f70-b7b6-5586193bcdef
md"#### Function is implemented in a generic way to support any amount of b[x] coefficients."

# ╔═╡ 1e5a39c0-6b5b-4d29-826e-e68c31c8f34e
function R2_is_bad(df; sigma=missing)
    degree = ncol(df[!,r"b"])
    # build mass_std*degree matrix, with each col exponentiated to col's index
    t = repeat(d.mass_std, 1, degree)
    t = hcat(map(.^, eachcol(t), 1:degree)...)
    s = [
        begin
            # calculate product on coefficient's vector
            b = collect(r[r"b"])
            μ = r.a .+ t * b
            s = ismissing(sigma) ? exp(r.log_σ) : sigma
            rand(MvNormal(μ, s))
        end
        for r ∈ eachrow(df)
    ]
    s = vcat(s'...);

    r = mean.(eachcol(s)) .- d.brain_std;
    v1 = var(r, corrected=false)
    v2 = var(d.brain_std, corrected=false)
    1 - v1 / v2
end

# ╔═╡ e8f6166f-20e5-40a2-a8fc-692ea8704b6d
md"### Code 7.7"

# ╔═╡ f74d5bc3-747b-4b12-9fb2-7e4a19cfe2c1
@model function model_m7_2(mass_std, brain_std)
    a ~ Normal(0.5, 1)
    b ~ MvNormal([0, 0], 10)
    μ = @. a + b[1]*mass_std + b[2]*mass_std^2
    log_σ ~ Normal()
    brain_std ~ MvNormal(μ, exp(log_σ))
end

# ╔═╡ fa9522ce-d91d-4cbd-9a1a-fe63fd2b7d15
begin
	m7_2_ch = sample(model_m7_2(d.mass_std, d.brain_std), NUTS(), 10000)
	m7_2 = DataFrame(m7_2_ch)
	describe(m7_2)
end

# ╔═╡ 41fc983a-9209-4ff5-b9fa-5798c7d387af
md"### Code 7.8"

# ╔═╡ 0f8d2476-7379-4bec-8a6f-e6ef81066ebb
md"#### Implemented the sample in a general way."

# ╔═╡ 4632e385-2417-4683-9e31-72825cf5501c
@model function model_m7_n(mass_std, brain_std; degree::Int)
    a ~ Normal(0.5, 1)
    b ~ MvNormal(zeros(degree), 10)
    # build matrix n*degree
    t = repeat(mass_std, 1, degree)
    # exponent its columns
    t = hcat(map(.^, eachcol(t), 1:degree)...)
    # calculate product on coefficient's vector
    μ = a .+ t * b
    
    log_σ ~ Normal()
    brain_std ~ MvNormal(μ, exp(log_σ))
end

# ╔═╡ bf763ccd-2c72-4529-b6ac-747f39f1b9cf
begin
	m7_3_ch = sample(model_m7_n(d.mass_std, d.brain_std, degree=3), NUTS(), 1000)
	m7_3 = DataFrame(m7_3_ch)
	describe(m7_3)
end

# ╔═╡ ef2e821a-d947-40e1-9855-79819aec5dbb
begin
	m7_4_ch = sample(model_m7_n(d.mass_std, d.brain_std, degree=4), NUTS(), 1000)
	m7_4 = DataFrame(m7_4_ch)
	describe(m7_4)
end

# ╔═╡ e1e48fb9-bb0f-4fe6-a34e-46ea2ced2510
begin
	m7_5_ch = sample(model_m7_n(d.mass_std, d.brain_std, degree=5), NUTS(), 1000)
	m7_5 = DataFrame(m7_5_ch)
	describe(m7_5)
end

# ╔═╡ bff2a56e-e799-4d50-b97b-79983df86565
@model function model_m7_6(mass_std, brain_std)
    a ~ Normal(0.5, 1)
    b ~ MvNormal(zeros(6), 10)
    μ = @. a + b[1]*mass_std + b[2]*mass_std^2 + b[3]*mass_std^3 + 
               b[4]*mass_std^4 + b[5]*mass_std^5 + b[6]*mass_std^6 
    brain_std ~ MvNormal(μ, 0.001)
end

# ╔═╡ 5c26fd02-6695-4f3d-8f17-37c9b2ada2e6
begin
	m7_6_ch = sample(model_m7_6(d.mass_std, d.brain_std), NUTS(), 1000)
	m7_6 = DataFrame(m7_6_ch)
	describe(m7_6)
end

# ╔═╡ 941248e9-c1ec-4d6a-8a2a-dbdd9aee17d8
md"### Code 7.10"

# ╔═╡ 97a58b2a-a200-4249-acf3-c52d28ca6aea
mass_seq = range(extrema(d.mass_std)...; length=100)

# ╔═╡ 7174161a-6536-4ec7-a89e-2c385de7052c
begin
	l = [
    	@. r.a + r.b * mass_seq
    	for r ∈ eachrow(m7_1)
	]
	l = vcat(l'...)
	μ = mean.(eachcol(l))
end

# ╔═╡ f4955154-ea66-4523-aacd-1a4beef6d41b
begin
	ci = PI.(eachcol(l))
	ci = vcat(ci'...)
end

# ╔═╡ caa5e742-7606-43fa-8ae3-a1287ff24951
scatter(d.mass_std, d.brain_std; title="1: R² = $(round(R2_is_bad(m7_1); digits=3))")

# ╔═╡ 2f524064-abcc-40a1-a21e-99a0a4636545
plot!(mass_seq, [μ μ]; fillrange=ci, c=:black, fillalpha=0.3)

# ╔═╡ 0f92b3fa-1b99-45f4-8a8e-5b7ea33a464d
md"### Reimplemented the brain_plot function to check my results."

# ╔═╡ cb2d9bdf-3cf3-48f2-a984-5440dba08ad9
function brain_plot(df; sigma=missing)
    degree = ncol(df[!,r"b"])
    # build mass_seq*degree matrix, with each col exponentiated to col's index
    t = repeat(mass_seq, 1, degree)
    t = hcat(map(.^, eachcol(t), 1:degree)...)
    l = [
        r.a .+ t * collect(r[r"b"])
        for r ∈ eachrow(df)
    ]
    l = vcat(l'...)
    μ = mean.(eachcol(l))
    ci = PI.(eachcol(l))
    ci = vcat(ci'...)

    r2 = round(R2_is_bad(df, sigma=sigma); digits=3)
    scatter(d.mass_std, d.brain_std; title="$degree: R² = $r2")
    plot!(mass_seq, [μ μ]; fillrange=ci, c=:black, fillalpha=0.3)
end;

# ╔═╡ 29f7fc44-1d5c-4912-a0b0-07a7e97cfe26
plot(
    brain_plot(m7_1),
    brain_plot(m7_2),
    brain_plot(m7_3),
    brain_plot(m7_4),
    brain_plot(m7_5),
    brain_plot(m7_6, sigma=0.001);
    size=(1000, 600)
)

# ╔═╡ bb21db38-31cf-4e88-86d2-922d62533769
md"### Code 7.11"

# ╔═╡ 4f7bfda2-297d-432d-83c9-c931f5f62868
i = 3;

# ╔═╡ 90d5c3e8-a045-4d96-aebc-2ec99609e10b
d_minus_i = d[setdiff(1:end,i),:];

# ╔═╡ 3babf620-dd01-45ae-b21c-e46e136057f6
function brain_loo_plot(model, data; title::String)
    (a, b) = extrema(data.brain_std)
    p = scatter(data.mass_std, data.brain_std; title=title, ylim=(a-0.1, b+0.1))
    mass_seq = range(extrema(data.mass_std)...; length=100)
    
    for i ∈ 1:nrow(data)
        d_minus_i = data[setdiff(1:end,i),:]
        df = DataFrame(sample(
				model(d_minus_i.mass_std, d_minus_i.brain_std), 
				NUTS(), 
				1000))

        degree = ncol(df[!,r"b"])
		
        # build mass_seq*degree matrix, with each col exponentiated to col's index
		
        t = repeat(mass_seq, 1, degree)
        t = hcat(map(.^, eachcol(t), 1:degree)...)
        l = [
            r.a .+ t * collect(r[r"b"])
            for r ∈ eachrow(df)
        ]
        l = vcat(l'...)
        μ = mean.(eachcol(l))
        plot!(mass_seq, μ; c=:black)
    end
    p
end

# ╔═╡ 9a5d927d-e816-4f5f-a460-c89a507ce1ae
model_m7_4 = (mass, brain) -> model_m7_n(mass, brain, degree=4)

# ╔═╡ b3778ab5-03b7-4ba5-9501-3e6bb643e4d5
plot(
    brain_loo_plot(model_m7_1, d, title="m7.1"),
    brain_loo_plot(model_m7_4, d, title="m7.4");
    size=(800, 400)
)

# ╔═╡ b55cffd3-d5ce-41db-968d-0b54b1b872a7
md"## 7.2 Entropy and accuracy."

# ╔═╡ e8765238-eec6-41c0-9c8c-148c9e51aad9
md"### Code 7.12"

# ╔═╡ fd1f8ae2-6839-4f7e-9141-66824490df72
p = [0.3, 0.7];

# ╔═╡ 3a3537ec-039a-434d-a7d4-15ab41eb97d8
-sum(p .* log.(p))

# ╔═╡ Cell order:
# ╠═d9015225-b635-42ad-ae96-fe7772a031b7
# ╠═727e041b-fe08-4f7d-ad19-e3c57772c967
# ╠═078fdb65-5bd6-4b0d-b672-617f5912480d
# ╟─5141dfdf-490d-479e-ac58-63c929ae8fd1
# ╟─611c3828-f01f-4721-830b-736b28620e7e
# ╠═ad511be2-fd88-4aba-a830-21e242d994f1
# ╟─4ee38b80-aded-40af-93ac-11a4ac2b9378
# ╠═5c7f1f26-02ac-4325-a06c-8c26c75e6cf4
# ╟─5a2f79ca-2d3e-4cd4-9de9-7c7e283a5052
# ╠═94b5539b-4374-40b8-9501-8d3df49f20c3
# ╠═2aafe78c-b1b1-4963-ad89-dbe09b546d96
# ╟─45122674-c5d2-4f61-af03-c74482618e90
# ╠═3b7ca48f-078f-4eef-9d0c-fa8724e711d3
# ╠═40c173f3-3dbb-4907-9188-632d34b3d7cf
# ╟─39dc811c-2742-4977-a203-c9fdcb5b218d
# ╠═85bb82da-efe7-4b34-8af5-cf60ae6e4757
# ╟─53806285-ea68-4f5c-b3fb-a382a894ec56
# ╠═59b9a24b-6967-418b-a4f8-e67ceb0a5474
# ╟─ed541621-2245-49f9-9672-31e20b9b7765
# ╟─47d9566c-93ec-4f70-b7b6-5586193bcdef
# ╠═1e5a39c0-6b5b-4d29-826e-e68c31c8f34e
# ╟─e8f6166f-20e5-40a2-a8fc-692ea8704b6d
# ╠═f74d5bc3-747b-4b12-9fb2-7e4a19cfe2c1
# ╠═fa9522ce-d91d-4cbd-9a1a-fe63fd2b7d15
# ╟─41fc983a-9209-4ff5-b9fa-5798c7d387af
# ╟─0f8d2476-7379-4bec-8a6f-e6ef81066ebb
# ╠═4632e385-2417-4683-9e31-72825cf5501c
# ╠═bf763ccd-2c72-4529-b6ac-747f39f1b9cf
# ╠═ef2e821a-d947-40e1-9855-79819aec5dbb
# ╠═e1e48fb9-bb0f-4fe6-a34e-46ea2ced2510
# ╠═bff2a56e-e799-4d50-b97b-79983df86565
# ╠═5c26fd02-6695-4f3d-8f17-37c9b2ada2e6
# ╟─941248e9-c1ec-4d6a-8a2a-dbdd9aee17d8
# ╠═97a58b2a-a200-4249-acf3-c52d28ca6aea
# ╠═7174161a-6536-4ec7-a89e-2c385de7052c
# ╠═f4955154-ea66-4523-aacd-1a4beef6d41b
# ╠═caa5e742-7606-43fa-8ae3-a1287ff24951
# ╠═2f524064-abcc-40a1-a21e-99a0a4636545
# ╟─0f92b3fa-1b99-45f4-8a8e-5b7ea33a464d
# ╠═cb2d9bdf-3cf3-48f2-a984-5440dba08ad9
# ╠═29f7fc44-1d5c-4912-a0b0-07a7e97cfe26
# ╟─bb21db38-31cf-4e88-86d2-922d62533769
# ╠═4f7bfda2-297d-432d-83c9-c931f5f62868
# ╠═90d5c3e8-a045-4d96-aebc-2ec99609e10b
# ╠═3babf620-dd01-45ae-b21c-e46e136057f6
# ╠═9a5d927d-e816-4f5f-a460-c89a507ce1ae
# ╠═b3778ab5-03b7-4ba5-9501-3e6bb643e4d5
# ╟─b55cffd3-d5ce-41db-968d-0b54b1b872a7
# ╟─e8765238-eec6-41c0-9c8c-148c9e51aad9
# ╠═fd1f8ae2-6839-4f7e-9141-66824490df72
# ╠═3a3537ec-039a-434d-a7d4-15ab41eb97d8
