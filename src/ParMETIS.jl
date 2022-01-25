module ParMETIS

using SparseArrays
using LinearAlgebra
import LightGraphs
using PARMETIS_jll: libparmetis

# Metis C API
include("parmetis_h.jl")
const options = fill(Cint(-1), METIS_NOPTIONS)
options[METIS_OPTION_NUMBERING] = 1

# Julia interface
"""
    ParMETIS.Graph
1-based CSR representation of a graph as defined in section 4.2
"Input and Output Formats used by ParMETIS" in the ParMETIS manual.
"""
struct Graph
    nvtxs::idx_t
    xadj::Vector{idx_t}
    adjncy::Vector{idx_t}
    vwgt::Vector{idx_t}
    Graph(nvtxs, xadj, adjncy) = new(nvtxs, xadj, adjncy)
    Graph(nvtxs, xadj, adjncy, vwgt) = new(nvtxs, xadj, adjncy, vwgt)
end

"""
    ParMETIS.graph(G::SparseMatrixCSC; check_hermitian=true)
Construct the 1-based CSR representation of the sparse matrix `G`.
If `check_hermitian` is `false` the matrix is not checked for being hermitian
before constructing the graph.
"""
function graph(G::SparseMatrixCSC; check_hermitian=true)
    if check_hermitian
        ishermitian(G) || throw(ArgumentError("matrix must be Hermitian"))
    end
    N = size(G, 1)
    xadj = Vector{idx_t}(undef, N+1)
    xadj[1] = 1
    adjncy = Vector{idx_t}(undef, nnz(G))
    adjncy_i = 0
    @inbounds for j in 1:N
        n_rows = 0
        for k in G.colptr[j] : (G.colptr[j+1] - 1)
            i = G.rowval[k]
            if i != j # don't include diagonal elements
                n_rows += 1
                adjncy_i += 1
                adjncy[adjncy_i] = i
            end
        end
        xadj[j+1] = xadj[j] + n_rows
    end
    resize!(adjncy, adjncy_i)
    return Graph(idx_t(N), xadj, adjncy)
end

"""
    ParMETIS.graph(G::LightGraphs.AbstractSimpleGraph)
Construct the 1-based CSR representation of the `LightGraphs` graph `G`.
"""
function graph(G::LightGraphs.AbstractSimpleGraph)
    N = LightGraphs.nv(G)
    xadj = Vector{idx_t}(undef, N+1)
    xadj[1] = 1
    adjncy = Vector{idx_t}(undef, 2*LightGraphs.ne(G))
    adjncy_i = 0
    for j in 1:N
        ne = 0
        for i in LightGraphs.outneighbors(G, j)
            ne += 1
            adjncy_i += 1
            adjncy[adjncy_i] = i
        end
        xadj[j+1] = xadj[j] + ne
    end
    resize!(adjncy, adjncy_i)
    return Graph(idx_t(N), xadj, adjncy)
end

"""
    perm, iperm = ParMETIS.permutation(G)
Compute the fill reducing permutation `perm`
and its inverse `iperm` of `G`.
"""
permutation(G) = permutation(graph(G))

function permutation(G::Graph)
    perm = Vector{idx_t}(undef, G.nvtxs)
    iperm = Vector{idx_t}(undef, G.nvtxs)
    vwgt = isdefined(G, :vwgt) ? G.vwgt : C_NULL
    METIS_NodeND(G.nvtxs, G.xadj, G.adjncy, vwgt, options, perm, iperm)
    return perm, iperm
end

"""
    ParMETIS.partition(G, n; alg = :KWAY)
Partition the graph `G` in `n` parts.
The partition algorithm is defined by the `alg` keyword:
 - :KWAY: multilevel k-way partitioning
 - :RECURSIVE: multilevel recursive bisection
"""
partition(G, nparts; alg = :KWAY) = partition(graph(G), nparts, alg = alg)

function partition(G::Graph, nparts::Integer; alg = :KWAY)
    part = Vector{idx_t}(undef, G.nvtxs)
    vwgt = isdefined(G, :vwgt) ? G.vwgt : C_NULL
    edgecut = fill(idx_t(0), 1)
    if alg === :RECURSIVE
        METIS_PartGraphRecursive(G.nvtxs, idx_t(1), G.xadj, G.adjncy, vwgt, C_NULL, C_NULL,
                                 idx_t(nparts), C_NULL, C_NULL, options, edgecut, part)
    elseif alg === :KWAY
        METIS_PartGraphKway(G.nvtxs, idx_t(1), G.xadj, G.adjncy, vwgt, C_NULL, C_NULL,
                            idx_t(nparts), C_NULL, C_NULL, options, edgecut, part)
    else
        throw(ArgumentError("unknown algorithm $(repr(alg))"))
    end
    return part
end

end # module
