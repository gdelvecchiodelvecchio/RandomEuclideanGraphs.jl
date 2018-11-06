
"""
euclidean_cost(f::Function, z::Array{Float64}, param...) = f(norm(z, 2), param...)
Evaluates a scalar field at a given point z: f must be a function which takes Array{Float64} as first argument.
Arbitrary number of parameters can be given.
cost(x .- y, 1) give the distance between x and y.
"""
euclidean_cost(f::Function, z::Array{Float64}, param...) = f(norm(z, 2), param...)
export euclidean_cost
"""
pcost(z::Float64, p::Float64) = z^p
"""
pcost(z::Float64, p::Float64) = z^p
export pcost
"""
function distancePBC(a::Float64,R::Float64)
	if abs(a)<0.5*R
		return abs(a)
	else
		return R-abs(a)
	end
end
"""
function distancePBC(a::Array{Float64},R::Float64)
	if norm(a,2)<0.5*R
		return norm(a,2)
	else
		return R-norm(a,2)
	end
end
export distancePBC
"""
circle_cost(f::Function, z::Array{Float64}, param...) = f(distancePBC(z, 1.), param...)
"""
circle_cost(f::Function, z::Array{Float64}, param...) = f(distancePBC(z, 1.), param...)
export circle_cost
"""
euclidean_distance(z::Array{Float64}) = norm(z, 2)
"""
euclidean_distance(z::Array{Float64}) = norm(z, 2)
export euclidean_distance
"""
wheight(e::AEdge, vm::VertexMap, f::Function)

Assosiate a weight to an edge e = (v,w) as f(v,w).
"""
function weight(e::AEdge, vm::VertexMap, f::Function)
    return f(vm.data[src(e)] .- vm.data[dst(e)])
end
export weight
"""
randomVertexMap!(vm::VertexMap, f::UnionAll, param...)

Function to randomly fill an existing vertex map according to
a given distribution f with parameters given by param...
"""
function randomVertexMap!(vm::VertexMap, f::UnionAll, param...)
    d = length(vm.data[1])
    vm.data = [rand(f(param...), d) for v in 1:length(vm.data)]
end
export randomVertexMap!
"""
fillEdgeMap!(g::AGraph, em::EdgeMap, f::Function)
Function to fill an existing edge map associated with g according to f.
"""
function fillEdgeMap!(g::AGraph, em::EdgeMap, f::Function)
    em.data = Dict((e,f(e)) for e in edges(g))
end

export fillEdgeMap!

#=
#enhanced src and dst to access edges

src2(g::AGraph, e::AEdge) = has_edge(g, e) ? src(e) : error("e=($(e.src),$(e.dst))∉g")
dst2(g::AGraph, e::AEdge) = has_edge(g, e) ? dst(e) : error("e=($(e.src),$(e.dst))∉g")
src2(g::AGraph, v::Int, w::Int) = has_edge(g, v, w) ? v : error("e=($(v),$(w))∉g")
dst2(g::AGraph, v::Int, w::Int) = has_edge(g, v, w) ? w : error("e=($(v),$(w))∉g")

=#

