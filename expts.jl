module Expt8001
    vec(eltype::Type, len::Int64) = Vector{eltype}(len)

    # Return an empty matrix of the same type as `A`, with element `eltype`
    # and of size `sz`
    empty(A::Diagonal, eltype::Type, sz) = Diagonal(vec(eltype, sz[1]))
    empty(A::Bidiagonal, eltype::Type, sz) = Bidiagonal(vec(eltype, sz[1]), vec(eltype, sz[1]-1), A.isupper)
    empty(A::Tridiagonal, eltype::Type, sz) = empty(Tridiagonal, eltype, sz)


    empty(::Type{UpperTriangular}, eltype::Type, sz) = UpperTriangular(zeros(Float64, sz))
    empty(::Type{LowerTriangular}, eltype::Type, sz) = LowerTriangular(zeros(Float64, sz))
    empty(::Type{Tridiagonal}, eltype::Type, sz) = Tridiagonal(vec(eltype, sz[1]-1), vec(eltype, sz[1]), vec(eltype, sz[1]-1))
    empty(::Type{Matrix}, eltype::Type, sz) = zeros(eltype, sz...)

    # Return an empty matrix which is the tightest type possible to hold
    # the result of performing * on `A` and `B`
    compat(*, A::Diagonal, B::Diagonal) = empty(A, Float64, (size(A,1), size(B,2)) )
    compat(*, A::Diagonal, B::AbstractMatrix) = empty(B, Float64, (size(A,1), size(B,2)) )
    compat(*, A::AbstractMatrix, B::Diagonal) = empty(A, Float64, (size(A,1), size(B,2)) )
    compat(*, A::Bidiagonal, B::Bidiagonal) = 
        if A.isupper && B.isupper
            empty(UpperTriangular, Float64, (size(A,1), size(A,1)))
        elseif !A.isupper && !B.isupper
            empty(LowerTriangular, Float64, (size(A,1), size(A,1)))
        else
            empty(Tridiagonal, Float64, (size(A,1), size(A,1)) )
        end
    compat(*, A::Bidiagonal, B::Tridiagonal) = empty(Matrix, Float64, (size(A,1), size(A,1)) )
    compat(*, A::Tridiagonal, B::Bidiagonal) = empty(Matrix, Float64, (size(A,1), size(A,1)) )
    compat(*, A::Tridiagonal, B::Tridiagonal) = empty(Matrix, Float64, (size(A,1), size(A,1)) )

    # Indexing
    struct Index
        i::Int64
        j::Int64
    end

    # Represents the shape/structure/indexes of a banded matrix
    abstract type AbstractShape end
    struct DiagonalShape <: AbstractShape
        m::Int64
        n::Int64
    end
    struct BandedShape <: AbstractShape
        m::Int64
        n::Int64
        below::Int64
        above::Int64
    end

    compose(*, A::DiagonalShape, B::DiagonalShape) = A   
    compose(*, A::DiagonalShape, B::AbstractShape) = B
    compose(*, A::AbstractShape, B::DiagonalShape) = A
    compose(*, A::BandedShape, B::BandedShape) = BandedShape(A.m, B.n, A.below+B.below, A.above+B.above)

    structure(D::Diagonal)::DiagonalShape = DiagonalShape(size(D,1), size(D,2))
    structure(B::Bidiagonal)::BandedShape = 
        if B.isupper
            BandedShape(size(B,1), size(B,2), 0, 1)
        else
            BandedShape(size(B,1), size(B,2), 1, 0)
        end
    structure(T::Tridiagonal)::BandedShape = BandedShape(size(T,1), size(T,2), 1, 1)
    
    iterateshape(D::DiagonalShape) = (Index(i,i) for i in 1:D.m)
    
    # FIXME this returns some indices more than once
    banded(n::Int64, m::Int64, below::Int64, above::Int64) =
        ( Index(i, clamp(i+j, 1, m))  for i in 1:n, j in -below:above )
    iterateshape(B::BandedShape) = banded(B.n, B.m, B.below, B.above)

    indexes(*, A::AbstractMatrix, B::AbstractMatrix) = 
        iterateshape(compose(*, structure(A), structure(B)))
    

    # Dot product of the `i`th row of `A` with the `j`th column of `B`
    dot(A::Diagonal,  i::Int64, B::AbstractMatrix, j::Int64) = @inbounds return A.diag[i] * B[i, j]

    @inline function dot(A::Bidiagonal, i::Int64, B::AbstractMatrix, j::Int64)
        (lim, evoff, Boff) = if A.isupper
            (size(A,2), 0, 1)
        else
            (1, -1, -1)
        end
        @inbounds return if i == lim
            A.dv[i] * B[i, j]
        else 
            A.dv[i] * B[i, j]  + A.ev[i + evoff] * B[i + Boff, j]
        end
    end

    @inline function dot(A::Tridiagonal, i::Int64, B::AbstractMatrix, j::Int64)
        # FIXME Try copying diagonals of Tridiagonal to banded form and working with them instead
        # With zeros in place can remove branching at edges
        # FIXME out of bounds if A is 2x2?
        @inbounds return if i == 1
            A.d[i] * B[i, j] + A.du[i] * B[i+1, j]
        elseif i == size(A, 2)
            A.d[i] * B[i, j] + A.dl[end] * B[i-1, j]
        else 
            A.d[i] * B[i, j] + A.dl[i-1] * B[i-1, j] + A.du[i] * B[i+1, j] 
        end
    end

    # Actual multiplication
    function A_mul_B(A::AbstractMatrix, B::AbstractMatrix)
        T = compat(*, A, B)

        for ndx in indexes(*, A, B)
            @inbounds T[ndx.i, ndx.j] = dot(A, ndx.i, B, ndx.j)
        end
        T
    end  

    # Performances tests to see where this approach looses time
    function a(A::AbstractMatrix, B::AbstractMatrix)
        # Cost of indexes()
        sum = 0
        for ndx in indexes(*, A, B)
            sum += ndx.i
        end
        sum
    end 

    function b(A::AbstractMatrix, B::AbstractMatrix)
        # Cost of dot()
        sum = 0.0
        for i in 1:100, j in 1:4
            sum += dot(A, i, B, j)
        end
        sum
    end 

    function c(A::AbstractMatrix, B::AbstractMatrix)
        # Cost of indexes() and dot()
        sum = 0.0
        for ndx in indexes(*, A, B)
            @inbounds sum += dot(A, ndx.i, B, ndx.j)
        end
        sum
    end 

end

ops = [(Expt8001.A_mul_B, A_mul_B!)]

mats = [n -> Diagonal(randn(n)),
        n -> Bidiagonal(randn(n), randn(n-1), true),
        n -> Bidiagonal(randn(n), randn(n-1), false),
        n -> Tridiagonal(randn(n-1), randn(n), randn(n-1))
        ]

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
        println("\n>>> ", typeof(err), "\t: ", typeof(m1), " ", opnew, " ", typeof(m2))
        showerror(STDOUT, err) #, catch_backtrace())
        println()
    end
end            