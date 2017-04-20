module Expt8001

    # Return an empty matrix of the same type as `A`, with element `eltype`
    # and of size `sz`
    empty(A::Diagonal, eltype::Type, sz) = Diagonal(Vector{Float64}(sz[1]))
    empty(A::Bidiagonal, eltype::Type, sz) = Bidiagonal(Vector{Float64}(sz[1]), Vector{Float64}(sz[1]-1), true)

    # Return an empty matrix which is the tightest type possible to hold
    # the result of performing * on `A` and `B`
    compat(*, A::Diagonal, B::Diagonal) = empty(A, Float64, (size(A,1), size(B,2)) )
    compat(*, A::Diagonal, B::AbstractMatrix) = empty(B, Float64, (size(A,1), size(B,2)) )
    compat(*, A::AbstractMatrix, B::Diagonal) = empty(A, Float64, (size(A,1), size(B,2)) )
    compat(*, A::Bidiagonal, B::Bidiagonal) = UpperTriangular(zeros(Float64, (size(A,1), size(A,1))))

    # Indexing
    struct Index
        i::Int64
        j::Int64
    end

    # FIXME this returns some indices more than once
    banded(n::Int64, m::Int64, above::Int64, below::Int64) =
      ( Index(i, clamp(i+j, 1, m))  for i in 1:n, j in -above:below )

    indexes(B::Bidiagonal) = banded(size(B,1), size(B,2), 0, 1)

    indexes(*, A::Diagonal, B::Diagonal) = (Index(i,i) for i in 1:length(A.diag))
    indexes(*, A::Diagonal, B::AbstractMatrix) = indexes(B)
    indexes(*, A::AbstractMatrix, B::Diagonal) = indexes(A)
    indexes(*, A::Bidiagonal, B::Bidiagonal) = banded(size(A,1), size(B,2), 0, 2)

    # Dot product of the `i`th row of `A` with the `j`th column of `B`
    dot(A::Diagonal,  i::Int64, B::AbstractMatrix, j::Int64) = @inbounds return A.diag[i] * B[i, j]
    function dot(A::Bidiagonal, i::Int64, B::AbstractMatrix, j::Int64)
        @inbounds if i == size(A, 2)
            return A.dv[i] * B[i, j]
        else
            return A.dv[i] * B[i, j]  + A.ev[i] * B[i+1, j]
        end
    end


    function A_mul_B(A::AbstractMatrix, B::AbstractMatrix)
        T = compat(*, A, B)

        for ndx in indexes(*, A, B)
            @inbounds T[ndx.i, ndx.j] = dot(A, ndx.i, B, ndx.j)
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