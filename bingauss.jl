using BinaryIntegers
using LinearAlgebra
using Combinatorics

function swapRows(A,i,j)
    A[i,:],A[j,:] = A[j,:],A[i,:];
end

function swapCols(A,i,j)
    A[:,i],A[:,j] = A[:,j],A[:,i];
end

function binGauss(A::Array{Int1})
    m,n = size(A);
    h = 1; k = 1;
    while h ≤ m && k ≤ n
        i_max = argmax(A[h:m,k])+h-1;
        if A[i_max,k] == 0
            k += 1;
        else
            swapRows(A,h,i_max);
            for i=1:m
                if i == h continue end #Skip current row
                A[i,k:n] = A[i,k:n] + A[i,k]*A[h,k:n]
            end
            h += 1; k += 1;
        end
    end
    return A
end

function binGauss(A::Array{T}) where T<:Integer
    return binGauss(Int1.(A))
end

function standardGen(G)
    m,n = size(G);
    for h=1:m
        p = findfirst(isequal(1),G[h,:]);
        if p != h swapCols(G,h,p); end
    end
    return G #Standardized Generator Matrix
end

"""Computes a basis for the nullspace of a binary matrix"""
function nullspace(H::Array{Int1})
    h,k = size(H);
    d = min(h,k);
    H_ = standardGen(binGauss(H));
    P = H_[:,d+1:end]'
    G = hcat(P,I)
    return G
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

H,G = hammingCode(3)
