using LinearAlgebra, StaticArrays

# Helpers
cisd(θ) = cosd(θ) + 1im*sind(θ)
axind(ax::Char) = UInt8(ax) - 0x77
axind(ax::Integer) = ax
θʜrd(useDeg::Bool) = useDeg ? 180 : π  # half rotation angle
cisrd(θ::Real, useDeg::Bool) = useDeg ? cisd(θ) : cis(θ)
sinrd(θ::Real, useDeg::Bool) = useDeg ? sind(θ) : sin(θ)
cosrd(θ::Real, useDeg::Bool) = useDeg ? cosd(θ) : cos(θ)
tanrd(θ::Real, useDeg::Bool) = useDeg ? tand(θ) : tan(θ)
asinrd(x::Real, useDeg::Bool) = useDeg ? asind(x) : asin(x)
acosrd(x::Real, useDeg::Bool) = useDeg ? acosd(x) : acos(x)
atanrd(x::Real, useDeg::Bool) = useDeg ? atand(x) : atan(x)
atanrd(y::Real, x::Real, useDeg::Bool) = useDeg ? atand(y,x) : atan(y,x)


# Rotation Matrices
rotmat(θ::Complex) = @SMatrix [θ.re -θ.im; θ.im θ.re]
rotmat(θ::Real; useDeg::Bool=false) = rotmat(cisrd(θ,useDeg))  # 2D from angle
function rotmat(ax::Union{Char, Integer}, θ::Complex)  # 3D principal from complex
    ax = axind(ax)
    (c,s) = reim(θ)
    ax == 1 ? (@SMatrix [1 0 0; 0 c -s; 0 s c]) :
    ax == 2 ? (@SMatrix [c 0 s; 0 1 0; -s 0 c]) :
    ax == 3 ? (@SMatrix [c -s 0; s c 0; 0 0 1]) :
    error("Axis must be 1, 2, 3 or 'x', 'y', 'z'")
end
function rotmat(ax::AbstractVector{<:Real}, θ::Complex)
    checksize(ax, 3)
    c, s = reim(θ)
    x, y, z = ax
    t = 1-c
    txx = t*x*x
    tyy = t*y*y
    tzz = t*z*z
    tyz = t*y*z
    tzx = t*z*x
    txy = t*x*y
    sx = s*x
    sy = s*y
    sz = s*z
    @SMatrix [
        txx+c  txy-sz tzx+sy
        txy+sz tyy+c  tyz-sx
        tzx-sy tyz+sx tzz+c
    ]
end
rotmat(ax::Union{Char, Integer, AbstractVector{<:Real}}, θ::Real; useDeg::Bool=false) = rotmat(ax, cisrd(θ,useDeg))  # 3D from axis (principal or unit vector) and angle
rotmat(ax₁, θ₁, ax₂, θ₂; useDeg::Bool=false) = rotmat(ax₁, θ₁; useDeg=useDeg) * rotmat(ax₂, θ₂; useDeg=useDeg)  # 3D from two elementary rotations
rotmat(ax₁, θ₁, ax₂, θ₂, ax₃, θ₃; useDeg::Bool=false) = rotmat(ax₁, θ₁, ax₂, θ₂; useDeg) * rotmat(ax₃, θ₃; useDeg=useDeg)  # 3D from three elementary rotations
rotmat(ax₁, θ₁, ax₂, θ₂, ax₃, θ₃, args... ; useDeg::Bool=false) = rotmat(ax₁, θ₁, ax₂, θ₂, ax₃, θ₃; useDeg=useDeg) * rotmat(args...; useDeg=useDeg)  # 3D from n elementary rotations (Only 3 are required for any general 3D rotation. This function may use heap allocation.)
function rotmat(v::AbstractVector{<:Real}; useDeg::Bool=false)  # 3D from rotation vector
    checksize(v, 3)
    V = norm(v)
    abs(V) < eps() ? SMatrix{3,3}(I) : rotmat(v/V, V; useDeg=useDeg)
end
function rotmat(q::Quat)  # 3D from quaternion
    wx = q.w*q.x
    wy = q.w*q.y
    wz = q.w*q.z
    xx = q.x*q.x
    xy = q.x*q.y
    xz = q.x*q.z
    yy = q.y*q.y
    yz = q.y*q.z
    zz = q.z*q.z
    @SMatrix [
      1-2(yy+zz)   2(xy-wz)   2(xz+wy)
        2(xy+wz) 1-2(xx+zz)   2(yz-wx)
        2(xz-wy)   2(yz+wx) 1-2(xx+yy)
    ]
end
rotmatd(args...) = rotmat(args...; useDeg=true)  # degree version of rotmat function

# Rotation Vectors
function rotvec(q::Quat; useDeg::Bool=false)  # from quaternion (versor)
    V = norm(q.v)
    V < eps() ? zero(q.v) :
    2 * atanrd(V, q.w, useDeg) / V .* q.v
end
function rotvec(R::AbstractMatrix{<:Real}; useDeg::Bool=false)  # from rotation matrix
    checksize(R, 3, 3)
    θʜ = θʜrd(useDeg)  # half cycle angle
    c=clamp((tr(R)-1)/2,-1,1)  # cos(θ)
    θ = acosrd(c, useDeg)
    abs(θ) < eps() ? SVector(0,0,0) :
    abs(θʜ - θ) < eps(Float64(θʜ)) ? let M=R+I; θʜ * sqrt.(max.(diag(M)/2,0)) .* sign.(M[:,findmax(diag(M))[2]]) end :
    θ / (2*sinrd(θ, useDeg)) * SVector(R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2])
end
rotvecd(args...) = rotvec(args...; useDeg=true)  # degree version of rotvec function

Quat(R::AbstractMatrix{<:Real}) = Quat(rotvec(R))  # Quat from rotation matrix

function rot(q::Quat, v::AbstractVector{<:Real})  # use a quaternion to rotate a point; checksize validates at transpile time
    checksize(v, 3)
    r = q * Quat(0, v[1], v[2], v[3]) * q'
    return SVector(r.x, r.y, r.z)
end
const ↺(R, X) = rot(R, X)

# ToDo
# * make a function: rot(rotation, object)
#     * rotation can be a rotmat, rotvec, quaternion, complex, euler angles and axes, etc.
#     * object can be a 2D or 3D vector (column, row, or some weird higher dimension) or a set of vectors or another rotation object.
#     * i would then use the ↺ symbol for this operation
#     * i might also make a rotinv function and use the ↻ symbol for that
#     * adjust rot(q::Quat, v::AbstractVector{<:Real}) overload to use cross product definition for efficiency
# * make a _safe version of functions that check sizes; for StaticArrays this can be handled at compile time by an overload for that specific size, and for other types it can be handled by checksize; both can call the _safe version. don't export the _safe versions; they are just for the nasty guts of those functions.

