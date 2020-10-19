
using Markdown
using InteractiveUtils

using Pkg, DrWatson

begin
	@quickactivate "StatisticalRethinkingStan"
	using StructuralCausalModels
	using StatisticalRethinking
end

md"## Clip-04-07-09s.jl"

md"### snippet 5.7"

DMA_1 = OrderedDict(
  :d => [:a, :m],
  :m => :a
);

DMA_dag_1 = DAG("DMA_dag1", DMA_1)

begin
	fname1 = joinpath(mktempdir(), "DMA_dag_1.dot")
	to_graphviz(DMA_dag_1, fname1)
	Sys.isapple() && run(`open -a GraphViz.app $(fname1)`)
end;

md"### snippet 5.9"

Text(pluto_string(basis_set(DMA_dag_1)))

md"##### `DMA_dag_1` has no implied conditional independencies.  This means all 3 variables are associated (even after conditioning). This is a testable implications."

adjustment_sets(DMA_dag_1, :a, :d)

md"##### `DMA_dag_1` has no adjustment sets when regressing `:d` on `:a`."

md"### snippet 5.8"

DMA_2 = OrderedDict(
  [:d, :m] => :a
)

DMA_dag_2 = DAG("DMA_dag_2", DMA_2)

begin
	fname2 = joinpath(mktempdir(), "DMA_dag_2.dot")
	to_graphviz(DMA_dag_2, fname2)
	Sys.isapple() && run(`open -a GraphViz.app $(fname2)`)
end;

Text(pluto_string(basis_set(DMA_dag_2)))

md"##### `DMA_dag_2` has 1 testable conditional independence (:d and :m are independent after conditioning on :a). Otherwise all pairs of variables are associated"

adjustment_sets(DMA_dag_2, :a, :d)

md"##### `DMA_dag_2` has no adjustment sets."

md"## End of clip-05-07-09s.jl"

