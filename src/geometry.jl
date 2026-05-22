using LinearAlgebra, StaticArrays

const rad2deg = 180 / π
const deg2rad = π / 180

"""
    vecmag2(x::AbstractMatrix, dims=nothing)

Squared magnitudes of vectors stored along a specific dimension of an array
"""
function vecmag2(x::AbstractArray{<:Real}, dims=nothing)
    dims === nothing && (dims = dim1(x))
    sum(x.^2, dims=dims)
end

"""
    vecmag(x::AbstractMatrix, dims=nothing)

Magnitudes of vectors stored along a specific dimension of an array
"""
function vecmag(x::AbstractArray{<:Real}, dims=nothing)
    dims === nothing && (dims = dim1(x))
    sqrt.(vecmag2(x, dims))
end

"""
    unit(x::AbstractMatrix, dims=nothing)

Unit vectors from vectors stored along a specific dimension of an array
"""
function unit(x::AbstractArray{<:Real}, dims=nothing)
    dims === nothing && (dims = dim1(x))
    x./vecmag(x, dims)
end

"""
    circle(n::Int=72; r::Real=1, o::AbstractVector{<:Real}=[0,0], includeEnd::Bool=false)

Make a regular n-gon, approximating a circle at increasing values of n.
Optional:
r: radius (default 1)
0: center (default [0,0])
includeEnd: flag to include a copy of the first point (off by default)
"""
circle(n::Int=72; r::Real=1, o::AbstractVector{<:Real}=[0,0], includeEnd::Bool=false) = ((0:n-1+includeEnd).*(360/n) .|> [cosd;;sind]).*r .+ o'

"""
    ellperim(a, b, n=3)

perimeter of an ellipse from semi-major axis a and semi-minor axis b using the first n+1 terms of the Ivory-Kummer series
"""
function ellperim(a, b, N=3)
    h = ((a-b)/(a+b))^2
    n = 1:N
    pi*(a+b)*(1 + sum(h.^n.*(binomial.(2*n,n)./((2*n.-1).*(4 .^n))).^2))
end

