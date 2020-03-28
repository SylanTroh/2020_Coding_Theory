using Combinatorics

struct binaryCode
    H::Array{Int1,2}
    G::Array{Int1,2}
    n::Int64
    k::Int64
    d::Int64
    M::Int64
end

function hammingCode(n::Int)
    H = [Int1((j>>i)%2) for i=0:n-1, j=1:2^n-1]; #Generate all nonzero columns
    H = standardGen(binGauss(H))
    G = nullspace(H)
    return H,G

end

function listWords(G::Array{Int1})
    S = combinations([G[i,:] for i=1:size(G)[1]])
    C = zeros(Int1,1,size(G)[2])
    C = vcat(C,hcat(sum.(S)...)')
    return(C)
end

function minWeight(C::Array{Int1})
    return minimum(filter(x->x>0,sum(Int64.(C),dims=2)))
end
