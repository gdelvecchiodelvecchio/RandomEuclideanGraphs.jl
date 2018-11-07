__precompile__()

module RandomEuclideanGraphs

using Erdos, ErdosExtras, Distributions, LinearAlgebra #this only for Julia 0.7.0.+


include("BasicFunctions.jl")
include("frustrationOBC.jl")
include("frustrationPBC.jl")
#=
export pk, euclidean_cost, euclidean_distance
export pcost, nearest_neighbors_bipartite, nearest_neighbors_monopartite
export nearest_neighbors_bipartite_inf, nearest_neighbors_monopartite_inf
export randomVerteMap!, fillEdgeMap!, weight
=#

end #module
