"""
nearest_neighbors_bipartite(g::AGraph, vm::VertexMap, em::EdgeMap)

Return an array of nv1 elements (nv1 is the number of (say) black vertices)
where in position i there is the number of white points between it and its
matched white point this could be a measure of the frustration because the black
point does not match the nearest white point but the k-th.
"""
function nearest_neighbors_bipartite(g::AGraph, vm::VertexMap, em::EdgeMap, p::Float64)
    n = Int64(nv(g)/2)
    ris = Array{Int64}(n)
    match_solution = minimum_weight_perfect_matching(g, em).mate
    @inbounds @fastmath for blackpoint in 1:n
        distances_from_black =  [euclidean_cost(pcost, vm.data[blackpoint] .- white ,p) for white in vm.data]     #verlet list
        dist = euclidean_cost(pcost, vm.data[blackpoint] .- vm.data[match_solution[blackpoint]], p)
        ris[blackpoint] = count( x -> (0. < x <= dist) , distances_from_black)
    end
    return ris
end
export nearest_neighbors_bipartite
"""
nearest_neighbors_monopartite(g::AGraph, vm::VertexMap, em::EdgeMap)

Return an array of nv1 elements (nv1 is the number of (say) black vertices)
where in position i there is the number of white points between it and its
matched white point this could be a measure of the frustration because the black
point does not match the nearest white point but the k-th.
"""
function nearest_neighbors_monopartite(g::AGraph, vm::VertexMap, em::EdgeMap, p::Float64)
    n = Int64(nv(g))
    ris = Array{Int64}(n)
    match_solution = minimum_weight_perfect_matching(g, em).mate
    @inbounds @fastmath for point in 1:n
        distances_from_point = [euclidean_cost(pcost,vm.data[point] .- others ,p) for others in vm.data]
        dist = euclidean_cost(pcost, vm.data[point] .- vm.data[match_solution[point]], p)
        ris[point] = count( x -> (0. < x < dist) , distances_from_point)
    end
    return ris
end
export nearest_neighbors_monopartite
"""
nearest_neighbors_bipartite_inf(g::AGraph, vm::VertexMap, em::EdgeMap)

Same as without inf but here consider the mean field model.
"""
function nearest_neighbors_bipartite_inf(g::AGraph, vm::VertexMap, em::EdgeMap)
    n = Int64(nv(g)/2)
    ris = Array{Int64}(n)
    match_solution = minimum_weight_perfect_matching(g, em).mate[1:n]
    @inbounds @fastmath for blackpoint in 1:n
        distances_from_black = [em.data[edge(g, blackpoint, white)] for white in n+1:2*n]
        dist = em.data[edge(g, blackpoint, match_solution[blackpoint])]
        ris[blackpoint] = count(x -> x <= dist, distances_from_black)
    end
    return ris
end

export nearest_neighbors_bipartite_inf
"""
nearest_neighbors_monopartite_inf(g::AGraph, vm::VertexMap, em::EdgeMap)

Same as without inf but here consider the mean field model.
"""
function nearest_neighbors_monopartite_inf(g::AGraph, vm::VertexMap, em::EdgeMap)
    n = nv(g)
    ris = Array{Int64}(n)
    match_solution = minimum_weight_perfect_matching(g, em).mate[1:n]
    ind = 1
    @inbounds for v in 1:n
            val = append!([em.data[edge(g,i,v)] for i in 1:v-1],[em.data[edge(g,v,i)] for i in v+1:n])
            dist = em.data[v < match_solution[v] ? edge(g, v, match_solution[v]) : edge(g, match_solution[v],v)]
            ris[ind] = count(x -> abs(x) <= dist, val)
            ind = ind + 1
    end
    return ris
end
export nearest_neighbors_monopartite_inf

#With this scheme it's easier to parralelize the code. Each function call on
#a single core with a low number of instances. It is also possible to test the monopartite and the bipartite
#case at the same time.
#It is possible to add an argument (a path) to this function in order to produce a file
"""
pk(G::AGraph, D::Real, nI::Int64, P::Float64, f::UnionAll, param...)

Given graph G (bipartite or monopartite) the function
computes the average probability that in the random euclidean
matching (assignment) problem, in the optimal matching, a given
point (black point) is linked with its k-th nearest (white) neighbor.
The function returns a tuple of vectors: the probability and its standard
error computed over the nI instances. The parameters are the following:
G is the graph
D is the dimension of the euclidean space
P is the exponent of the cost function z^P
f is the distribution used to generate points
param represents an indefinite number which f accepts: f(param)
"""
function pk(G::AGraph, D::Int64, nI::Int64, P::Float64, f::UnionAll, param...)
    n = is_bipartite(G) ? Int64(nv(G)/2) : nv(G)    #N for bi 2N for mono
    nInst = nI      #number of instances
    d = D           #dimension
    p = P           #exponent
    g = G           #graph
    ρ = f           #probability density
    par = param     #parameters fo the distribution e.g Unform(0.,1.)
    fr = Matrix{Int64}(n,nInst)



    #to speedup the computation it is better to iterate over the rows
    #(Julia native's Array type is the column vector)

    #initialization of graph variables for the computation

    vm = d < Inf ? VertexMap(g, x->rand(d)) : VertexMap(g,x->rand())
    em = EdgeMap(g, e->rand())
    #possibility to test infinte dimension


    if is_bipartite(g)
        if d < Inf
            @inbounds for i in 1:nInst
                randomVertexMap!(vm, ρ, par...)
                fillEdgeMap!(g, em, e->weight(e, vm, x->euclidean_cost(pcost, x, p)))
                fr[:,i] = nearest_neighbors_bipartite(g, vm, em, p)
            end
        else
            @inbounds for i in 1:nInst
                fillEdgeMap!(g, em, e->rand())
                fr[:,i] = nearest_neighbors_bipartite_inf(g, vm, em)
            end
        end
    else
        if d < Inf
            @inbounds for i in 1:nInst
                randomVertexMap!(vm, ρ, par...)
                fillEdgeMap!(g, em, e->weight(e, vm, x->euclidean_cost(pcost, x, p)))
                fr[:,i] = nearest_neighbors_monopartite(g, vm, em, p)
            end
        else
            @inbounds for i in 1:nInst
                fillEdgeMap!(g, em, e->rand())
                fr[:,i] = nearest_neighbors_monopartite_inf(g, vm, em)
            end
        end
    end



    probmat = Matrix{Float64}(nInst,n)
    prob = Array{Float64}(n)
    errors = Array{Float64}(n)



    #probability that inside a (hyper)sphere around a black matched point there are j-th white points (no frustration correpsonds to prob[:,1])
    @inbounds @fastmath for inst in 1:nInst
        for k in 1:n
            r = fr[:,inst] .== k
            probmat[inst,k] = mean(r)
        end
    end

    @inbounds  @fastmath for i in 1:n
        prob[i] = mean(probmat[:,i])
        errors[i] = std(probmat[:, i]) / sqrt(nInst)
    end

    return prob, errors
end

export pk
