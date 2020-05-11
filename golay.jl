using BinaryIntegers: Int1, BinGauss
using Random
using LinearAlgebra
using Nemo
using Primes

function ord(a)
    F = a.parent
    if a == F(1) return 1 end
    p = F.p
    n = degree(F)
    q = p^n
    x = a
    for i=2:q-1
        x = x*a
        if x == F(1) return i end
    end
end

function generator(F)
    p = F.p
    n = degree(F)
    q = p^n
    for i=0:q-1
        a = elem(F,gen(F),i)
        if ord(a) == q-1 return a; end
    end
end

function elem(F,α,a)
    p = F.p
    n = degree(F)
    q = p^n
    a = a%q
    return sum(((a,r) = divrem(a,p);r)*α^i for i=0:q-1)
end

function Jacobsthal(p,n)
    if !Primes.isprime(p)
        print("Error: Not Prime")
        return
    end
    q = p^n
    R,x = PolynomialRing(GF(p),"x")
    g = x^n
    while !isirreducible(g)
        g = x^n
        for i = 0:n-1
            g += Random.rand(0:p-1)*x^i
        end
    end
    F,y = FiniteField(g,"y")
    α = generator(F)

    squares = Set((α^i)^2 for i=0:q)
    character = x -> if x == 0 return 0; else return (x ∈ squares ? 1 : -1) end
    J = [character(elem(F,α,a) - elem(F,α,b)) for a=0:q-1,b=0:q-1]
    return J
end

function isprimepower(q::T) where T<:Integer
    l = log(q)
    for n = 1:floor(Int,l)
        p = q^(1/n)
        if isinteger(p)
            p = Int(p)
            if Primes.isprime(p)
                return p,n
            end
        end
    end
    return (0,0)
end

function payley(q)
    if q%4 != 0 return Nothing end
    p,n = isprimepower(q-1)
    if p > 0 && q > 0
        Q = Jacobsthal(p,n)
        Q = Q-I
        H = vcat(ones(Int,1,q-1),Q)
        H = hcat(ones(Int,q,1),H)
        return H
    else
        if div(q,2)%4 != 2 return Nothing end
        p,n = isprimepower(div(q,2)-1)
        if p > 0 && q > 0
            Q = Jacobsthal(p,n)
            H_ = vcat(ones(Int,1,div(q,2)-1),Q)
            H_ = hcat(ones(Int,div(q,2),1),H_)
            H_[1,1] = 0
            H_ = [b == 0 ? [1 -1; -1 -1] : b*[1 1;1 -1] for b in H_];
            H = zeros(Int,q,q)
            for i=1:q
                for j=1:q
                    H[i,j] = H_[ceil(Int,i/2),ceil(Int,j/2)][-1*(i%2-1)+1,-1*(j%2-1)+1]
                end
            end
            return H
        end
    end
end

function normHadamard!(H)
    for (j,a) in enumerate(H[1,:])
        if a == -1
            H[:,j] = H[:,j]*(-1)
        end
    end
    for (i,a) in enumerate(H[:,1])
        if a == -1
            H[i,:] = H[i,:]*(-1)
        end
    end
end

function hadamard(n)
    if n == 1 return [1] end;
    if n == 2 return [1 1; 1 -1] end;
    k = 0;
    while (H=payley(n)) == Nothing
        if n%4 != 0 return end
        n = div(n,2)
        k += 1
    end
    for i=1:k
        H = kron(H,[1 1; 1 -1])
    end
    normHadamard!(H)
    return H
end

function binaryHadamard(H::Array{Int1})
    H = [a==1 ? 0 : 1 for a in H]
    return Int1.(H)
end

function binaryHadamard(n::Int)
    H = hadamard(n)
    if H == nothing return H end
    H = [a==1 ? 0 : 1 for a in H]
    return Int1.(H)
end

function golay()
    G = hcat([ones(Int1,11,1),I,zeros(Int1,11,1),binaryHadamard(12)[2:end,2:end]]...)
    G = vcat(G,hcat(zeros(Int1,1,12),ones(Int1,1,12)))
    return Int1.(G)
end

function nordstromRobinson()
    G = golay()
    c = Int1.([1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0])
    BinGauss.equivCode!(c,G)
    C = BinGauss.listWords(G)
    function isNordstrom(x)
        if sum(Int.(x[1:8])) == 0 return true end
        if x[8] == 1 && sum(Int.(x[1:8])) == 2 return true end
        return false
    end
    C = [C[i,:]' for i=1:size(C)[1] if isNordstrom(C[i,:])]
    C = vcat(C...)
    return C[:,9:end]
end
