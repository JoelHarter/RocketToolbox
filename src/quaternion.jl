# Quat <: Number — quaternions behave like a scalar (non-commutative) division
# ring. `*` is the Hamilton product (no risk of element-wise traps), and
# generic numeric code (zero/one/promote) just works. Heads-up: this means
# Quat is not a StaticArrays type; downstream code that wants index/iterate
# behavior should go through the `.q` accessor.
using LinearAlgebra, StaticArrays

struct Quat{T <: Real} <: Number
    w::T
    x::T
    y::T
    z::T
end

# constructors
Quat(w::Real, x::Real=false, y::Real=false, z::Real=false) = Quat(promote(w, x, y, z)...)
Quat(q::NTuple{4, Real}) = Quat(q...)  # 4-vec constructor
Quat(q::AbstractArray) = Quat(Tuple(q))
Quat(w::Real, v::NTuple{3, Real}) = Quat(w, v...)  # scalar & 3-vec constructor
Quat(w::Real, v::AbstractArray) = Quat(w, Tuple(v))
Quat(z::Complex{T}) where T = Quat(z.re, z.im)  # imaginary numbers interpreted as i component

# promotions
Base.promote_rule(::Type{Quat{T1}}, ::Type{Quat{T2}}) where {T1,T2}           = Quat{promote_type(T1, T2)}
Base.promote_rule(::Type{Quat{T1}}, ::Type{T2}) where {T1, T2<:Real}          = Quat{promote_type(T1, T2)}
Base.promote_rule(::Type{Quat{T1}}, ::Type{Complex{T2}}) where {T1, T2<:Real} = Quat{promote_type(T1, T2)}
Base.convert(::Type{Quat{T}}, q::Quat) where T = Quat(T(q.w), T(q.x), T(q.y), T(q.z))
Base.convert(::Type{Quat{T}}, w::Real) where T = Quat(T(w))
Base.convert(::Type{Quat{T}}, c::Complex) where T = Quat(Complex{T}(c))

# accessors
function Base.getproperty(q::Quat, S::Symbol)
    S === :q  ? SVector(q.w, q.x, q.y, q.z) :  # all 4 components as SVector
    S === :v  ? SVector(q.x, q.y, q.z) :        # vector part
    S === :Q  ? norm(q) :                       # quaternion norm
    S === :V  ? norm(q.v) :                     # norm of vector part
    S === :v1 ? normalize(q.v) :                # normalized vector part
    S === :q1 ? normalize(q) :                  # normalized quaternion (versor)
    getfield(q, S)
end
Base.real(q::Quat) = q.w
Base.imag(q::Quat) = q.v
Base.isreal(q::Quat) = all(q.v .== 0)
Base.isfinite(q::Quat) = all(isfinite, q.q)
Base.isinf(q::Quat) = any(isinf, q.q)
Base.zero(::Type{Quat{T}}) where T = Quat(zero(T))
Base.zero(q::Quat) = zero(typeof(q))
Base.one(::Type{Quat{T}}) where T = Quat(one(T))
Base.one(q::Quat{Float64})::Quat{Float64} = Quat(1.0)  # transpiler entry point; must precede generic overload
Base.one(q::Quat{T}) where T = Quat(one(T))

# Operators and Functions
Base.:+(q₁::Quat, q₂::Quat) = Quat(q₁.q + q₂.q)                  #   q₁ +  q₂
Base.:-(q₁::Quat, q₂::Quat) = Quat(q₁.q - q₂.q)                  #   q₁ -  q₂
Base.:-(q::Quat) = Quat(-q.q)                                    #  -q
Base.:+(r::Real, q::Quat) = Quat(r + q.w, q.x, q.y, q.z)         #   r  +  q
Base.:+(q::Quat, r::Real) = r  +  q                              #   q  +  r
Base.:-(r::Real, q::Quat) = Quat(r - q.w, -q.x, -q.y, -q.z)      #   r  -  q
Base.:-(q::Quat, r::Real) = Quat(q.w - r,  q.x,  q.y,  q.z)      #   q  -  r
Base.conj(q::Quat) = Quat(q.w, -q.x, -q.y, -q.z)                 #   ̅q             conjugate
Base.adjoint(q::Quat) = conj(q)                                   #   q' ≡ ̅q        use adjoint as conjugate
Base.:*(r::Real, q::Quat) = Quat(r * q.q)                        #   r  *  q       scaling
Base.:*(q::Quat, r::Real) = r * q                                #   q  *  r       scaling
Base.:*(q₁::Quat, q₂::Quat) = Quat(                              #   q  *  q       Hamilton product
        q₁.w*q₂.w -q₁.x*q₂.x -q₁.y*q₂.y -q₁.z*q₂.z,
        q₁.w*q₂.x +q₁.x*q₂.w +q₁.y*q₂.z -q₁.z*q₂.y,
        q₁.w*q₂.y -q₁.x*q₂.z +q₁.y*q₂.w +q₁.z*q₂.x,
        q₁.w*q₂.z +q₁.x*q₂.y -q₁.y*q₂.x +q₁.z*q₂.w)
LinearAlgebra.dot(q₁::Quat, q₂::Quat) = dot(q₁.q, q₂.q)          #   q  ·  q
LinearAlgebra.norm(q::Quat) = norm(q.q)                          #  ‖q‖
Base.abs2(q::Quat) = dot(q, q)                                    #  |q|²
Base.abs(q::Quat) = norm(q)                                       #  |q| ≡ ‖q‖
Base.:/(q::Quat, r::Real) = Quat(q.q / r)                        #   q  /  r
Base.:\(r::Real, q::Quat) = q / r                                #   r  \  q
Base.inv(q::Quat) = conj(q) / abs2(q)                            #   q⁻¹
Base.:/(n::Number, q::Quat) = n * inv(q)                         #   n  /  q
Base.:\(q::Quat, n::Number) = inv(q) * n                         #   q  \  n
function Base.exp(q::Quat)
    V = q.V
    e = exp(q.w)
    V < 1e-10 ? e*q : e*Quat(cos(V), sin(V)/V .* q.v)
end
function Base.log(q::Quat)
    Q = q.Q
    V = q.V
    V < 1e-10 ? Quat(log(Q)) : Quat(log(Q), acos(q.w / Q)/V * q.v)
end
Base.:^(q::Quat, n::Integer) = Base.power_by_squaring(q, n)  # transpiler emits C while-loop via emit_power_by_squaring!; literal n still inlined by emit_int_pow!
Base.:^(q::Quat, n::Real)    = isinteger(n) ? q^Int(n) : exp(n * log(q))
function Base.:sqrt(q::Quat)
    Q = norm(q)
    sqrt(Q) * normalize(q+Q)
end

# Versors (unit quaternions used to describe orientation and rotation)
function slerp(q₁::Quat, q₂::Quat, t::Real)  # spherical linear interpolation
    d = dot(q₁, q₂)
    if d < 0
        q₂ = -q₂
        d = -d
    end
    if d > 0.9995  # use linear interpolation (lerp) and renormalize for tiny angles
        return normalize((1-t)*q₁ + t*q₂)
    end
    θ = acos(d)
    cscθ = rsqrt(1 - d^2)
    sin((1-t)*θ)*cscθ * q₁ + sin(t*θ)*cscθ * q₂
end
squad(q₁::Quat, q₂::Quat, q₃::Quat, q₄::Quat; t::Real) = slerp(slerp(q₁,q₄,t), slerp(q₂,q₃,t), 2*t*(t-1))  # spherical quadrangle interpolation

# ToDo:
# function snerp(q::Quat; t::Real)
# end

# 3-vec (interpreted as rotation vector) constructor
function Quat(v::NTuple{3, Real}; useDeg::Bool=false)
    V = norm(v)
    s,c = useDeg ? sincosd(V/2) : sincos(V/2)
    Quat(c, (V < 1e-10 ? 0.5 : s/V) .* v)
end

# imaginary part global variables
const iq = Quat{Bool}(0,1,0,0)  # quaternion unit i
const jq = Quat{Bool}(0,0,1,0)  # quaternion unit j
const kq = Quat{Bool}(0,0,0,1)  # quaternion unit k

# display
Base.show(io::IO, ::MIME"text/plain", q::Quat) = show(io, q)
function Base.show(io::IO, q::Quat)
    s = get(io, :compact, false)::Bool ? "" : " "
    p=s*"+"*s; m=s*"-"*s
    sep(x) = isfinite(x) ? "" : "*"   # juxtaposition fails for NaN/Inf — emit explicit * so output is paste-able
    print(io, "$(q.w)"*(signbit(q.x) ? m : p)*"$(abs(q.x))$(sep(q.x))iq"*(signbit(q.y) ? m : p)*"$(abs(q.y))$(sep(q.y))jq"*(signbit(q.z) ? m : p)*"$(abs(q.z))$(sep(q.z))kq")
end
function Base.show(io::IO, q::Quat{Bool})
    function QuatBoolStr(q)
        I="1";O="0";c=","
        "Quat{Bool}("*(q.w ? I : O)*c*(q.x ? I : O)*c*(q.y ? I : O)*c*(q.z ? I : O)*")"
    end
    print(io,
        q.x ? !(q.y|q.z|q.w) ? "iq" : QuatBoolStr(q) :
        q.y ? !(    q.z|q.w) ? "jq" : QuatBoolStr(q) :
        q.z ? !(        q.w) ? "kq" : QuatBoolStr(q) :
        q.w ?          "Quat(true)" : QuatBoolStr(q) )
end

# ToDo:
# * implement spherical n-degree interpolation (I'll call it snerp)
# * from rotmat
# * to/from euler angles
