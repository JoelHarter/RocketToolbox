using LinearAlgebra

function rk4(f, y₀, h)
    k₁ = f(y₀)
    k₂ = f(y₀ + h/2 * k₁)
    k₃ = f(y₀ + h/2 * k₂)
    k₄ = f(y₀ + h * k₃)
    y = y₀ + h/6 * (k₁ + 2k₂ + 2k₃ + k₄)
end

# Multi-argument overload: rk4(f, h, x₁, x₂, ...) where f(x₁, x₂, ...) → (ẋ₁, ẋ₂, ...)
# Flattens inputs into a single SVector, runs rk4, splits results back to original types.
function rk4(f, h, args...)
    # Flatten all args into one SVector
    flat = SVector{sum(_nscalars, args)}(vcat(map(_tovec, args)...))

    # Record sizes for splitting
    sizes = map(_nscalars, args)
    types = map(typeof, args)

    # Wrap f to pack/unpack
    function f_flat(y)
        xs = _split(y, sizes, types)
        dx = f(xs...)
        SVector{length(y)}(vcat(map(_tovec, dx)...))
    end

    result = rk4(f_flat, flat, h)
    return _split(result, sizes, types)
end

_nscalars(::Real) = 1
_nscalars(v::StaticVector) = length(v)

_tovec(x::Real) = SVector(Float64(x))
_tovec(v::StaticVector) = SVector{length(v), Float64}(v)

function _split(flat, sizes, types)
    out = []
    i = 1
    for (n, T) in zip(sizes, types)
        chunk = flat[SVector{n}(i:i+n-1)]
        push!(out, T <: Real ? chunk[1] : T(chunk...))
        i += n
    end
    return Tuple(out)
end
