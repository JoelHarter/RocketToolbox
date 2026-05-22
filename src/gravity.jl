using LinearAlgebra

# acceleration of gravity at ECI point
function grav(r::ECI, useJ₂::Bool=false)
    r⁻² = /(norm2(r))
    r⁻³ = r⁻² * sqrt(r⁻²)
    g = -Earth.μ * r * r⁻³
    if useJ₂
        g += Earth.K_J₂ * r⁻²*r⁻³ * (5 * sqr(r.Z) * r⁻² .- SVector(1, 1, 3)) .* r
    end
    return ECI(g...)
end

# gravity gradient ∂g/∂r at ECI point.
# Spherical: Γ_sph = (μ/r⁵)(3 r rᵀ − r² I₃).
# J₂ adds the Hessian of the J₂ potential; symmetric and traceless.
function gravgrad(r::ECI, useJ₂::Bool=false)
    r²  = norm2(r)
    r⁻² = /(r²)
    r⁻⁵ = r⁻² * r⁻² * sqrt(r⁻²)
    Γ = Earth.μ * r⁻⁵ * (3 * r⊗r - r² * I(3))
    if useJ₂
        W = 5 * sqr(r.Z) * r⁻² .- SVector(1, 1, 3)
        Γ += Earth.K_J₂ * r⁻⁵ * (
            r ⊗ ((10 * r.Z * r⁻²) * (unitZ - (r.Z * r⁻²) * r)) +
            diagm(W) - 5 * r⁻² * (W .* r) ⊗ r)
    end
    return Γ
end
