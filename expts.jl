module Expt8001

    empty(A::Diagonal, eltype::Type, sz) = Diagonal(Vector{Float64}(sz[1]))
    empty(A::Bidiagonal, eltype::Type, sz) = Bidiagonal(Vector{Float64}(sz[1]), Vector{Float64}(sz[1]-1), true)

    compat(*, A::Diagonal, B::Diagonal) = empty(A, Float64, (size(A,1), size(B,2)) )
    compat(*, A::Diagonal, B::AbstractMatrix) = empty(B, Float64, (size(A,1), size(B,2)) )
    compat(*, A::AbstractMatrix, B::Diagonal) = empty(A, Float64, (size(A,1), size(B,2)) )
    compat(*, A::Bidiagonal, B::Bidiagonal) = UpperTriangular(Array{Float64}(size(A,1), size(A,1)))

    dot(A::Diagonal,  i, B::AbstractMatrix, j) = A.diag[i] * B[i, j]
    function dot(A::Bidiagonal, i, B::AbstractMatrix, j)
        if i == size(A, 2)
            return A.dv[i] * B[i, j]
        else
            return A.dv[i] * B[i, j]  + A.ev[i] * B[i+1, j]
        end
    end

    indexes(D::Diagonal) = ((i,i) for i in 1:length(D.diag))
    function indexes(B::Bidiagonal)
        n = size(B, 1)
        m = n + 1
        ((mod(i,m) + div(i,m), mod(i,m) + 2div(i,m)) for i in 1:(2n-1))
    end
    function indexes(U::UpperTriangular)
        n = size(U, 1)
        ( (i,j)  for i in 1:n, j in 1:n if i<=j )
    end


    function A_mul_B(A::AbstractMatrix, B::AbstractMatrix)
        T = compat(*, A, B)

        for (i, j) in indexes(T)
            T[i,j] = dot(A, i, B, j)
        end
        T
    end    

end

ops = [(Expt8001.A_mul_B, A_mul_B!)]

mats = [n -> Diagonal(randn(n)),
        n -> Bidiagonal(randn(n), randn(n-1), true)]

for (opnew, opold) in ops, m1f in mats, m2f in mats
    n = 10
    m1, m2 = m1f(n), m2f(n)
    
    fullres = Matrix{Float64}(n, n)
    try
        opold(fullres, full(m1), full(m2))
    catch err
        print("FULL OP FAILED ", typeof(err), "\t: ", typeof(m1), " ", opold, " ", typeof(m2), "\n")
        continue
    end
    try
        if !(full(opnew(m1, m2)) ≈ fullres)
            print("NOT EQUAL TO FULL OP ", typeof(m1), " ", opnew, " ", typeof(m2), " ", "\n")
        end
    catch err
        print(typeof(err), "\t: ", typeof(m1), " ", opnew, " ", typeof(m2), "\n")
    end
end            