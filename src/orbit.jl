using StaticArrays
using LinearAlgebra

"""
OrbEl{T<:Real}
6-Element Kepler Orbit Type
"""
struct OrbEl{T<:Real}
    a::T # [m] semi-major axis
    e::T # [1] eccentricity
    i::T # [°] inclination
    Ω::T # [°] longitude of ascending node
    ω::T # [°] argument of perigee
    θ::T # [°] true anomaly
end

"""
    eci2orbitplane_mat(O::OrbEl)

2×3 matrix that projects an ECI vector onto the orbit plane (u, v) coordinates.
Rows are the perifocal basis vectors f̂ (toward perigee) and ĝ (90° ahead).
The transpose (3×2) maps orbit plane back to ECI.
"""
function eci2orbitplane_mat(O::OrbEl)
    (si, ci) = sincosd(O.i)
    (sΩ, cΩ) = sincosd(O.Ω)
    (sω, cω) = sincosd(O.ω)
    @SMatrix [
        cΩ*cω-sΩ*ci*sω   sΩ*cω+cΩ*ci*sω   si*sω ;
       -cΩ*sω-sΩ*ci*cω  -sΩ*sω+cΩ*ci*cω   si*cω
    ]
end

"""
    project_to_orbit(r::ECI, P::SMatrix{2,3}, a, e)

Project position `r` onto the target orbit defined by semi-major axis `a`,
eccentricity `e`, and 2×3 ECI-to-orbit-plane matrix `P`.

Returns `(r_D, v_D, θ)` — the position, velocity, and true anomaly [deg]
on the orbit at the projected angle.

Steps:
Projects an ECI position onto the orbit plane, calculates its true anomaly,
projects that outward from the center of the earth onto the elliptical orbit,
converts back to ECI position, and calculates the velocity at that point on the orbit.
"""
function project_to_orbit(r::ECI, P::SMatrix{2,3}, a::Float64, e::Float64)
    e2_1 = 1.0 - e * e

    # Project to orbit plane (2D)
    rpo = P * SVector{3}(r[1], r[2], r[3])
    x, y = rpo[1], rpo[2]
    x2, y2 = x * x, y * y

    # Scale from angular direction onto the ellipse
    s = a / (x2 + y2 / e2_1) * (sqrt(x2 + y2) - e * x)
    u, v = x * s, y * s

    # Back to ECI
    Pt = P'  # 3×2
    r_D = ECI((Pt * SVector{2}(u, v))...)

    # True anomaly
    θ = atand(v, u)

    # Velocity: vis-viva magnitude, direction from conic geometry
    r_mag = norm(r_D)
    v_mag = sqrt(Earth.μ * (2.0 / r_mag - 1.0 / a))
    vdir = SVector{2}(-v, (u + a * e) * e2_1)
    vdir_n = vdir / norm(vdir)
    v_D = ECI((v_mag * Pt * vdir_n)...)

    return (r_D, v_D, θ)
end

"""
    project_to_orbit(r::ECI, O::OrbEl)

Convenience overload: computes the 2×3 matrix from orbital elements.
"""
function project_to_orbit(r::ECI, O::OrbEl)
    P = eci2orbitplane_mat(O)
    return project_to_orbit(r, P, Float64(O.a), Float64(O.e))
end

"""
    closest_on_orbit(r::ECI, P::SMatrix{2,3}, a, e)

Find the point on the target orbit closest (angularly) to `r`.
Equivalent to `project_to_orbit` — the name clarifies intent at the call site.
Returns `(r_D, v_D, θ)`.
"""
closest_on_orbit(r::ECI, P::SMatrix{2,3}, a::Float64, e::Float64) =
    project_to_orbit(r, P, a, e)

"""
    closest_on_orbit(r::ECI, O::OrbEl)

Convenience overload: computes the 2×3 matrix from orbital elements.
"""
closest_on_orbit(r::ECI, O::OrbEl) = project_to_orbit(r, O)

"""
ECI6(O::OrbEl)
ECI position and velocity.
"""
function ECI6(O::OrbEl)
    (si,ci) = sincosd(O.i)
    (sΩ,cΩ) = sincosd(O.Ω)
    (sω,cω) = sincosd(O.ω)
    (sθ,cθ) = sincosd(O.θ)
    p = O.a * (1 - sqr(O.e))
    c₁ = ECI(si*sΩ, -si*cΩ, ci)
    f₁ = ECI(cΩ*cω-sΩ*ci*sω, sΩ*cω+cΩ*ci*sω, si*sω)
    y₁ = c₁ × f₁
    ECI6(
        p / (1+O.e*cθ) * (cθ*f₁+sθ*y₁),  # r
        sqrt(Earth.μ/p) * (-sθ*f₁+(O.e+cθ) * y₁)  # v
    )
end

function OrbEl(r::ECI, v::ECI)
    ε = 1e-10  # const
    ε² = sqr(ε)  # const
    r̂ = unit(r)
    c = r × v
    f = v × c - Earth.μ * r̂  # Laplace vector
    ĉ = unit(c)
    f¹ = norm(f)
    f̂ = f/f¹
    i = acosd(ĉ[3])
    if i>ε
        Ω = atand(ĉ[1], -ĉ[2])
        ω = acosd(ĉ[1:2] ⋅ f̂[1:2]) * sign(f[3])
    else
        Ω = 0.0
        norm2(f) > ε² ? ω = atand(f[2], f[1]) : ω = 0.0
    end
    p = norm2(c) / Earth.μ
    e = f¹ / Earth.μ
    a = p / (1 - sqr(e))
    θ = acosd(f̂ ⋅ r̂) * sign(r̂ ⋅ v)
    return OrbEl(a, e, i, Ω, ω, θ)
end
OrbEl(s::ECI6) = OrbEl(s.r, s.v)
