using StaticArrays
using LinearAlgebra

# Define types for position data in various coordinate systems
abstract type AbstractCoord3D{T<:Real} <: FieldVector{3,T} end  # supertype for all 3D position coordinates
struct ECI{T}<:AbstractCoord3D{T}  # Earth-Centered, Inertial
    X::T  # [m]
    Y::T  # [m]
    Z::T  # [m]
end
struct ECEF{T}<:AbstractCoord3D{T}  # Earth-Centered, Earth-Fixed
    X::T  # [m]
    Y::T  # [m]
    Z::T  # [m]
end
struct Geo{T}  # Geodetic — not a Euclidean vector; lat/lon in degrees, alt in meters
    lat::T  # [°] (ϕ) latitude
    lon::T  # [°] (λ) longitude
    alt::T  # [m] (h) altitude, i.e. height above ellipsoid
end
struct ENU{T}<:AbstractCoord3D{T}  # East-North-Up
    E::T  # [m] east
    N::T  # [m] north
    U::T  # [m] up
end
struct NED{T}<:AbstractCoord3D{T}  # North-East-Down
    N::T  # [m] north
    E::T  # [m] east
    D::T  # [m] down
end
struct LVLH{T}<:AbstractCoord3D{T}  # Local Vertical Local Horizontal (orbit-relative frame)
    x::T  # [m] in-track
    y::T  # [m] cross-track (opposite angular momentum)
    z::T  # [m] radial (toward Earth center)
end
struct BodyCoord{T}<:AbstractCoord3D{T}  # Body-fixed coordinate frame
    x::T  # [m] or [rad/s] or whatever — developer knows the units
    y::T
    z::T
end

# ── 6D position + velocity state ──────────────────────────────────────────────
# Both fields transform together: r by the position transform, v by the
# corresponding relative-vector transform for that frame pair.
struct ECI6{T}
    r::ECI{T}
    v::ECI{T}
end
struct ECEF6{T}
    r::ECEF{T}
    v::ECEF{T}
end
struct LVLH6{T}
    r::LVLH{T}
    v::LVLH{T}
end
struct BodyCoord6{T}
    r::BodyCoord{T}
    v::BodyCoord{T}
end

# TODO: LVLH (Local Vertical Local Horizontal) for relative navigation
#       (rendezvous, proximity ops): z toward Earth center, y opposite orbit
#       angular momentum, x in-track. LVLH_rel and LVLH_pv to follow.

function Base.show(io::IO, ::MIME"text/plain", p::AbstractCoord3D)
    show(io, p)
end
function Base.show(io::IO, p::AbstractCoord3D)
    compact = get(io, :compact, false)::Bool
    T = compact ? nameof(typeof(p)) : typeof(p)
    s = compact ? "," : ", "
    compact && (p=round.(p,sigdigits=13))
    print(io,T,"(",p[1],s,p[2],s,p[3],")")
end
function Base.show(io::IO, ::MIME"text/plain", p::Geo)
    show(io, p)
end
function Base.show(io::IO, p::Geo)
    compact = get(io, :compact, false)::Bool
    T = compact ? nameof(typeof(p)) : typeof(p)
    s = compact ? "," : ", "
    lat = compact ? round(p.lat, sigdigits=13) : p.lat
    lon = compact ? round(p.lon, sigdigits=13) : p.lon
    alt = compact ? round(p.alt, sigdigits=13) : p.alt
    print(io, T, "(", lat, s, lon, s, alt, ")")
end

function Base.show(io::IO, ::MIME"text/plain", rv::Union{ECI6,ECEF6,BodyCoord6})
    show(io, rv)
end
function Base.show(io::IO, rv::Union{ECI6,ECEF6,BodyCoord6})
    compact = get(io, :compact, false)::Bool
    T = compact ? nameof(typeof(rv)) : typeof(rv)
    s = compact ? "," : ", "
    r = compact ? round.(rv.r, sigdigits=13) : rv.r
    v = compact ? round.(rv.v, sigdigits=13) : rv.v
    print(io, T, "(", r[1], s, r[2], s, r[3], s, v[1], s, v[2], s, v[3], ")")
end

# Mixed-frame arithmetic guard — prevent silent ECI+ECEF etc.
Base.:+(a::AbstractCoord3D, b::AbstractCoord3D) = error("Frame mismatch: $(nameof(typeof(a))) + $(nameof(typeof(b)))")
Base.:-(a::AbstractCoord3D, b::AbstractCoord3D) = error("Frame mismatch: $(nameof(typeof(a))) - $(nameof(typeof(b)))")
# Same-frame arithmetic — FieldVector provides these but the mixed-frame guard above is more specific,
# so we need explicit same-type methods to win dispatch.
Base.:+(a::P, b::P) where {P<:AbstractCoord3D} = P(a[1]+b[1], a[2]+b[2], a[3]+b[3])
Base.:-(a::P, b::P) where {P<:AbstractCoord3D} = P(a[1]-b[1], a[2]-b[2], a[3]-b[3])

# SVector interop: allows raw computation results to be added to typed coords
Base.:+(a::P, b::SVector{3}) where {P<:AbstractCoord3D} = P(a[1]+b[1], a[2]+b[2], a[3]+b[3])
Base.:+(a::SVector{3}, b::P) where {P<:AbstractCoord3D} = P(a[1]+b[1], a[2]+b[2], a[3]+b[3])
Base.:-(a::P, b::SVector{3}) where {P<:AbstractCoord3D} = P(a[1]-b[1], a[2]-b[2], a[3]-b[3])

θ_E(t::Real) = cis(t*Earth.Ω) # rotation angle of earth based on time as a complex number
function ECEF(r::ECI, θ::Complex)
    ECEF(θ.re*r.X-θ.im*r.Y,θ.im*r.X+θ.re*r.Y,r.Z)
end
ECEF(r::ECI, t::Real) = ECEF(r, θ_E(t))
function ECI(r::ECEF, θ::Complex)
    ECI(θ.re*r.X+θ.im*r.Y,-θ.im*r.X+θ.re*r.Y,r.Z)
end
ECI(r::ECEF, t::Real) = ECI(r, θ_E(t))

Geo(lat::T, lon::T) where T = Geo(lat, lon, zero(T)) # altitude h defaults to zero if you just call Geo(lat,lon)
"""
Convert from ECEF to Geodetic using Vermeille's method
https://www.researchgate.net/publication/225460650_An_analytical_method_to_transform_geocentric_into_geodetic_coordinates
This model will break if the ECEF point is inside the evolute or on the singular disc described in the paper - way down deep inside the earth.
The paper provides a test to tell if points are there and other equations to use if so. We will skip them and assume we only use points in space, or on or just below the earth's surface, and not deep inside the earth.
This method is technically an approximation of the true answer, but it is within computer precision for locations anywhere near the earth, and is much faster than iterative methods.
"""
function Geo(r::ECEF)
    o=Σsq(r.X,r.Y)
    p=o*Earth.a⁻²  # called "p" in the paper (Vermeille)
    q=Earth.a⁻⁴b²*sqr(r.Z)
    rr=(p+q-Earth.e⁴)/6  # called "r" in the paper; renamed to avoid collision
    s=Earth.e⁴*p*q
    t=.5*sqr(cbrt(sqrt(8*cube(rr)+s)+sqrt(s)))
    u=rr+t+sqr(rr)/t
    v=sqrt(sqr(u)+Earth.e⁴*q)
    w=Earth.e²*(u+v-q)/(2*v)
    k=(u+v)/(sqrt(sqr(w)+u+v)+w)
    D=k*sqrt(o)/(k+Earth.e²)
    g=rss(D,r.Z)
    Geo(2*atand(r.Z/(g+D)),atand(r.Y,r.X),(k-Earth.a⁻²b²)/k*g)
end

function ECEF(p::Geo)
    (sϕ,cϕ) = sincosd(p.lat)
    (sλ,cλ) = sincosd(p.lon)
    N=Earth.a*rsqrt(1. - Earth.e²*sqr(sϕ))
    Nhcϕ = cϕ*(N+p.alt)
    ECEF(Nhcϕ*cλ, Nhcϕ*sλ, (Earth.a⁻²b²*N+p.alt)*sϕ)
end

function ENU(r::ECEF, o::ECEF, R::SMatrix{3,3,<:Real})
    ENU(r-o, R)
end
function ENU(r::ECEF, o::ECEF)
    og=Geo(o)
    ENU(r-o, og.lat, og.lon)
end
function ENU(r::ECEF, o::Geo)
    ENU(r-ECEF(o), o.lat, o.lon)
end
function ENU(v::ECEF, R::SMatrix{3,3,<:Real})
    ENU((R*v)...)
end

# Ω is purely in z: cross([0,0,Ω], r) = [-Ω*r[2], Ω*r[1], 0]
# Coriolis-only helpers (no frame rotation) — operate on raw SVectors
rel_eci2ecef(v::AbstractVector, r::AbstractVector) = SVector(v[1] + Earth.Ω*r[2], v[2] - Earth.Ω*r[1], v[3])
rel_ecef2eci(v::AbstractVector, r::AbstractVector) = SVector(v[1] - Earth.Ω*r[2], v[2] + Earth.Ω*r[1], v[3])
# Full relative-vector transforms: frame rotation (θ) + Coriolis correction (r_ecef)
function ECEF(v::ECI, θ::Complex, r_ecef::ECEF; rel::Bool)
    vr = SVector(θ.re*v[1]-θ.im*v[2], θ.im*v[1]+θ.re*v[2], v[3])
    ECEF(rel_eci2ecef(vr, r_ecef)...)
end
ECEF(v::ECI, t::Real, r_ecef::ECEF; rel::Bool) = ECEF(v, θ_E(t), r_ecef; rel)
function ECI(v::ECEF, θ::Complex, r_ecef::ECEF; rel::Bool)
    vr = rel_ecef2eci(v, r_ecef)
    ECI(θ.re*vr[1]+θ.im*vr[2], -θ.im*vr[1]+θ.re*vr[2], vr[3])
end
ECI(v::ECEF, t::Real, r_ecef::ECEF; rel::Bool) = ECI(v, θ_E(t), r_ecef; rel)

# ── BodyCoord↔ECI transforms ──────────────────────────────────────────────────
# All transforms share this parameter convention:
#   q     — BodyCoord→ECI attitude quaternion  (rot(q, v_body) gives v in ECI; active)
#   r₀    — ECI position of the BodyCoord origin [m]
#   v₀    — ECI velocity of the BodyCoord origin [m/s]  (pv transforms only)
#   ω_eci — angular velocity of body, expressed in ECI frame [rad/s] (pv transforms only)

# Position
BodyCoord(r::ECI, q::Quat, r₀::ECI) = BodyCoord(rot(conj(q), r - r₀)...)
ECI(r::BodyCoord, q::Quat, r₀::ECI) = ECI((r₀ + rot(q, r))...)

# Free vector — re-expression only (gravity, thrust, sensor readings, etc.)
# No Coriolis: the vector is simply rotated into the new frame.
BodyCoord(v::ECI, q::Quat; rel::Bool) = BodyCoord(rot(conj(q), v)...)
ECI(v::BodyCoord, q::Quat; rel::Bool) = ECI(rot(q, v)...)

# 6D position + velocity state
# Velocity includes the Coriolis correction for the rotating body frame:
#   v_bc = rot(q', v_eci - v₀) - ω_bc × r_bc
#   v_eci = rot(q,  v_bc + ω_bc × r_bc) + v₀
function BodyCoord6(rv::ECI6, q::Quat, r₀::ECI, v₀::ECI, ω_eci::ECI)
    r_bc = BodyCoord(rv.r, q, r₀)
    ω_bc = BodyCoord(ω_eci, q; rel=true)
    vrot = rot(conj(q), rv.v - v₀)          # rot(q', v_eci - v₀)
    vcor = cross(ω_bc, r_bc)                # ω_bc × r_bc
    BodyCoord6(r_bc, BodyCoord(vrot[1]-vcor[1], vrot[2]-vcor[2], vrot[3]-vcor[3]))
end
function ECI6(rv::BodyCoord6, q::Quat, r₀::ECI, v₀::ECI, ω_eci::ECI)
    r_eci = ECI(rv.r, q, r₀)
    ω_bc  = BodyCoord(ω_eci, q; rel=true)
    vcor  = cross(ω_bc, rv.r)               # ω_bc × r_bc
    vrot  = rot(q, SVector(rv.v[1]+vcor[1], rv.v[2]+vcor[2], rv.v[3]+vcor[3]))
    ECI6(r_eci, ECI(vrot[1]+v₀[1], vrot[2]+v₀[2], vrot[3]+v₀[3]))
end

ecef2eci(t::Real) = let (s,c)=sincos(Earth.Ω*t); @SMatrix[c;-s;0.;;s;c;0.;;0.;0.;1.] end
eci2ecef(t::Real) = let (s,c)=sincos(Earth.Ω*t); @SMatrix[c;s;0.;;-s;c;0.;;0.;0.;1.] end
ecef2enu(ϕ::Real,λ::Real) = let (sϕ,cϕ)=sincosd(ϕ), (sλ,cλ)=sincosd(λ); @SMatrix[-sλ;cλ;0.;;-cλ*sϕ;-sλ*sϕ;cϕ;;cλ*cϕ;sλ*cϕ;sϕ] end
enu2ecef(ϕ::Real,λ::Real) = let (sϕ,cϕ)=sincosd(ϕ), (sλ,cλ)=sincosd(λ); @SMatrix[-sλ;-cλ*sϕ;cλ*cϕ;;cλ;-sλ*sϕ;sλ*cϕ;;0.;cϕ;sϕ] end
ecef2ned(ϕ::Real,λ::Real) = let (sϕ,cϕ)=sincosd(ϕ), (sλ,cλ)=sincosd(λ); @SMatrix[-cλ*sϕ;-sλ*sϕ;cϕ;;-sλ;cλ;0.;;-cλ*cϕ;-sλ*cϕ;-sϕ] end
ned2ecef(ϕ::Real,λ::Real) = let (sϕ,cϕ)=sincosd(ϕ), (sλ,cλ)=sincosd(λ); @SMatrix[-cλ*sϕ;-sλ;-cλ*cϕ;;-sλ*sϕ;cλ;-sλ*cϕ;;cϕ;0.;-sϕ] end
const enu2ned = const ned2enu = @SMatrix[0;1;0;;1;0;0;;0;0;-1]

ENU(v::ECEF, ϕ::Real, λ::Real) = ENU(v, ecef2enu(ϕ,λ))

Geo(r::ECI, t::Union{Real,Complex}) = Geo(ECEF(r,t))
ECI(r::Geo, t::Union{Real,Complex}) = ECI(ECEF(r),t)

NED(p::ENU)=NED(p.N,p.E,-p.U)
ENU(p::NED)=ENU(p.E,p.N,-p.D)
