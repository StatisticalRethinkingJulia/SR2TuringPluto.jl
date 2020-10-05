
using Markdown
using InteractiveUtils

macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingTuring"
	using Turing
	using StatisticalRethinking
	using PlutoUI
end

md"## Fig 2.5.1t"

md"""

It is not intended to show how to use Stan (yet)!

This notebook demonstrates simple PlutoUI interactivity."""

md"### 1. Create a Turing model:"

@model m2_0t(n, k) = begin
  theta ~ Beta(1, 1) # prior
  k ~ Binomial(n, theta) # model
  return k, theta
end;

md"### 2. Generate observed data."

md"##### n can go from 1:9"

@bind n Slider(1:18, default=9)

md"### 3. Sample for varying n:"

begin
	dfa2_0t = DataFrame()
	k = sum([1,0,1,1,1,0,1,0,1,1,0,1,1,1,0,1,0,1][1:n])
	chns = sample(m2_0t(n, k), NUTS(0.65), 1000)
	dfa2_0t.theta = reshape(chns[:theta].data, length(chns))
	Text(precis(dfa2_0t; io=String))
end

md"### 6. Show the posterior."

begin
  plot(xlims=(0.0, 1.0), ylims=(0.0, 4.0), leg=false)
  hline!([1.0], line=(:dash))
  density!(dfa2_0t.theta, line=(:dash))
 end

md"##### A few ways to summarize these samples."

md"###### Normal distribution."

part2_0t = Particles(dfa2_0t)

md"###### MAP based on samples."

quap2_0 = quap(dfa2_0t)

md"###### Turing quap estimate."

quap2_0t = quap(m2_0t(n, k))

md"## End of Fig2.5.1t.jl"

