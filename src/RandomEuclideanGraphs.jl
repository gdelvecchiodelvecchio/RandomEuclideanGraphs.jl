__precompile__()

module RandomEuclideanGraphs

using Erdos, ErdosExtras, Distributions


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
