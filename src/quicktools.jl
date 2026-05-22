# This file contains miscellanious definitions of useful tools, shorthands, etc.

using LinearAlgebra, StaticArrays

"""
    dim1(x::AbstractArray)

Infer the primaray dimension for processing arrays vector-wise and return its index:
* index of first dimension whose length is more than 1
* 1 if **every** dimension has unit length (there is no direction, but dimension 1 is usable)
* error (or 0 if allowEmpty is on) if **any** dimension has zero length (there is no data in the array)
This is conceptually the direction in which the array elements are laid out in memory.
"""
function dim1(x::AbstractArray; allowEmpty::Bool=false)
    s = size(x)
    any(s .== 0) ? allowEmpty ? 0 : error("Must be non-empty array") : something(findfirst(>(1), s), 1)
end
rsqrt(x::Real) = 1/sqrt(x)  # reciprocal sqrt. on some machines this is faster than the sqrt function itself
plusminus(A, B) = (A+B, A-B)  # plus and minus
const ± = plusminus  # plus and minus
outer(a::AbstractVector, b::AbstractVector) = a * b'  # outer product (quck version for vectors)
outer(A, B) = A .* reshape(B, (ntuple(_ -> 1, ndims(A))..., size(B)...))  # outer product (general)
outer(A) = outer(A, A)  # self outer product
norm2(v::AbstractVector) = dot(v, v)
norm2(A::StaticVector{<:Any,<:Real}) = dot(A, A)
sqr(x) = x*x
cube(x) = x*x*x
Σsq(args...) = sum(sqr,args) # sum of squares
rss(args...) = sqrt(Σsq(args...)) # root of sum of squares
Base.:/(x) = inv(x)  # unary division (reciprocal)
const ⊗ = outer  # outer product
const ⊠ = kron  # Kronecker product

# unit basis vectors (3D, Float64)
const unitX = @SVector [1.0, 0.0, 0.0]
const unitY = @SVector [0.0, 1.0, 0.0]
const unitZ = @SVector [0.0, 0.0, 1.0]

checksize(v::AbstractVector, len::Integer) = length(v) == len || error("Vector length mismatch: Expected length $len, got length $(length(v))")
checksize(M::AbstractMatrix, rows::Integer, cols::Integer) = size(M,1) == rows && size(M,2) == cols || error("Matrix size mismatch: Expected size ($rows, $cols), got size $(size(M))")
checksize(A::AbstractArray, args...) = checksize(A, (args...))
SV2(v::AbstractVector)             = (checksize(v, 2); SVector{2}(v))
SV3(v::AbstractVector)             = (checksize(v, 3); SVector{3}(v))
SV4(v::AbstractVector)             = (checksize(v, 4); SVector{4}(v))
SVn(v::AbstractVector, n::Integer) = (checksize(v, n); SVector{n}(v))  # NOTE for compile time known n, use SVn(v,Val(n)). For runtime-only known n, regular Vector is generally better than SVector; use this function only if you really really need an SVector.


"""
    horner(x::Number, c::Union{Tuple,AbstractVector}; k0::Integer=0)

Sum a series of multiples of consecutive powers of x using Horner's method.
x is a number to raise to consecutive integer powers.
c is a list of factors for the different powers of x.
By default, the first value in c corresponds to x^0 (i.e. it is just a constant added onto the rest of the series).
The power of x corresponding to the first value in c can be specified using k0.

EXAMPLE
# compute   y = 5 + 7x + 6.8x^2 - x^3   with already-known x
c = (5, 7, 6.8, -1)
y = horner(x, c)

# compute   y = 3x + x^4     where x=8.1
y = horner(8.1, (0, 3, 0, 0, 1))
# or
x = 8.1;    y = horner(x, (3, 0, 0, 1); k0=1)

# write the base10 number 8675309 in a fun way
horner(10,(9,0,3,5,7,6,8))
"""
function horner(x::Number, c::Union{Tuple,AbstractVector}; k0::Integer=0)
    n=length(c)
    y=c[n];
    for i = n-1:-1:1
        y = fma(y,x,c[i])
    end
    k0==0 ? y : k0==1 ? x*y : x^k0*y  # manually handle trivial cases without calling power function
end

