using BinaryIntegers: Int1, BinGauss

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
    H = BinGauss.standardGen(BinGauss.binGauss(H))
    G = BinGauss.nullspace(H)
    return H,G
end
