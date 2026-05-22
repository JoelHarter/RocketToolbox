# RocketToolbox functional test suite.
# Run from the repo root: julia test/test.jl
# Covers all non-graphical modules (plotez, showearth require Plots/GLMakie).

include(joinpath(@__DIR__, "..", "include.jl"))

using LinearAlgebra, StaticArrays

pass = 0; fail = 0
function check(name, val)
    if val; global pass += 1
    else;   global fail += 1; println("FAIL: $name"); end
end

# ── QuickTools ────────────────────────────────────────────────────────────────
check("dim1 vector",      dim1([1,2,3]) == 1)
check("dim1 matrix",      dim1([1 2; 3 4]) == 1)
check("dim1 row vec",     dim1([1 2 3]) == 2)
check("rsqrt",            rsqrt(4.0) ≈ 0.5)
check("plusminus",        plusminus(3,1) == (4,2))
check("± alias",          (3 ± 1) == (4,2))
check("sqr",              sqr(5) == 25)
check("cube",             cube(3) == 27)
check("norm2",            norm2(SVector(3.0,4.0)) ≈ 25.0)
check("Σsq",              Σsq(3,4) == 25)
check("rss",              rss(3.0,4.0) ≈ 5.0)
check("outer vecs",       outer([1,2],[3,4]) == [3 4; 6 8])
check("SV3",              SV3([1,2,3]) == SVector(1,2,3))
check("horner",           horner(2.0,(1,2,3)) ≈ 17.0)

# ── Geometry ──────────────────────────────────────────────────────────────────
check("rad2deg",          rad2deg * π ≈ 180)
check("deg2rad",          deg2rad * 180 ≈ π)
check("vecmag single",    only(vecmag([3.0,4.0,0.0])) ≈ 5.0)
check("vecmag2 single",   only(vecmag2([3.0,4.0,0.0])) ≈ 25.0)
check("unit single norm", norm(unit([3.0,4.0,0.0])) ≈ 1.0)
M = [3. 4.; 4. 3.; 0. 0.]
check("vecmag matrix",    vecmag(M) ≈ [5. 5.])
check("unit matrix",      all(isapprox.(vecmag(unit(M)), 1.0; atol=1e-14)))
check("circle size",      size(circle(8)) == (8,2))
check("circle radius",    all(norm.(eachrow(circle(36))) .≈ 1.0))
check("ellperim circle",  ellperim(1.0,1.0) ≈ 2π)

# ── Quaternion ────────────────────────────────────────────────────────────────
qi = Quat(0.,1.,0.,0.)
check("Quat identity",    Quat(1.0) == Quat(1,0,0,0))
check("iq constant",      iq == Quat(0,1,0,0))
check("Hamilton qi²=-1",  qi*qi ≈ Quat(-1.0))
check("Quat conj",        conj(qi) ≈ -qi)
check("Quat norm",        norm(Quat(1.,2.,3.,4.)) ≈ sqrt(30))
check("Quat inv",         qi * inv(qi) ≈ Quat(1.0))
check("slerp t=0",        slerp(Quat(1.0), qi, 0.0) ≈ Quat(1.0))
check("slerp t=1",        slerp(Quat(1.0), qi, 1.0) ≈ qi)
check("Quat(SVector{4})", Quat(SVector(0.,1.,0.,0.)) == qi)
check("Quat(SVector{3})", abs(Quat(SVector(π/2, 0., 0.))) ≈ 1.0)
check("Quat(SVector{3})==NTuple", Quat(SVector(π,0.,0.)) ≈ Quat((π,0.,0.)))

# ── Rotation ──────────────────────────────────────────────────────────────────
R1 = rotmat(1, π/2)
check("rotmat x90",       R1 * [0.,1.,0.] ≈ [0.,0.,1.])
check("rotmatd z90",      rotmatd(3, 90.0) * [1.,0.,0.] ≈ [0.,1.,0.])
check("rotvec roundtrip", rotvec(rotmat(2, 1.2)) ≈ [0.,1.2,0.])
check("Quat(SVector) rotvec", rot(Quat(rotvec(R1)), [0.,1.,0.]) ≈ [0.,0.,1.])
check("rot via quat",     rot(Quat((0., π/2, 0.)), [1.,0.,0.]) ≈ [0.,0.,-1.])
check("rotmat from quat", rotmat(Quat((0.,0.,π/4))) * [1.,0.,0.] ≈ [cos(π/4),sin(π/4),0.])
check("rotmat↔Quat roundtrip", rotmat(Quat(rotvec(rotmat(2,0.7)))) ≈ rotmat(2,0.7))

# ── Gravity ───────────────────────────────────────────────────────────────────
r = ECI(7e6, 0.0, 0.0)
g = grav(r)
check("grav direction",   g[1] < 0 && abs(g[2]) < 1e-10)
check("grav magnitude",   norm(g) ≈ Earth.μ / 7e6^2)
check("grav returns ECI", g isa ECI)

# ── Coordinate ────────────────────────────────────────────────────────────────
ecef0 = ECEF(1178979.204878298, 3239218.743409690, 5933489.835830858)
geo   = Geo(ecef0)
check("ecef→geo lat",     abs(geo.lat - 60.0) < 0.001)
check("ecef→geo lon",     abs(geo.lon - 70.0) < 0.001)
check("ecef→geo alt",     abs(geo.alt - 5e5) < 1.0)
check("geo→ecef roundtrip", norm(ECEF(geo) - ecef0) < 1e-4)
ecef_rt = ECEF(ECI(ecef0, 3600.), 3600.)
check("eci↔ecef roundtrip", norm(ecef_rt - ecef0) < 1e-6)
v_eci = SVector(100.,200.,300.); r_pos = SVector(7e6,0.,0.)
check("vel eci↔ecef roundtrip", norm(rel_ecef2eci(rel_eci2ecef(v_eci,r_pos),r_pos) - v_eci) < 1e-10)

# ── New coordinate types ──────────────────────────────────────────────────────
check("BodyCoord",    BodyCoord(1.,2.,3.) isa AbstractCoord3D)
check("ECI",          ECI(0.,7800.,0.) isa AbstractCoord3D)
check("ECEF",         ECEF(0.,7800.,0.) isa AbstractCoord3D)
check("ENU",          ENU(1.,0.,0.) isa AbstractCoord3D)
check("NED",          NED(0.,1.,0.) isa AbstractCoord3D)
r_ecef = ECEF(7e6, 0., 0.)
vi = ECI(0., 7800., 0.)
ve = ECEF(vi, 0.0, r_ecef; rel=true)
check("ECI→ECEF vel", ve isa ECEF)
check("ECEF→ECI vel", ECI(ve, 0.0, r_ecef; rel=true) ≈ vi)
pv = ECI6(ECI(7e6,0.,0.), vi)
check("ECI6 fields",    pv.r == ECI(7e6,0.,0.) && pv.v == vi)
check("ECEF6 fields",   ECEF6(ECEF(7e6,0.,0.), ve).r == ECEF(7e6,0.,0.))

# ── Coordinate type algebra ────────────────────────────────────────────────────
r1 = ECI(7e6, 0., 0.); r2 = ECI(6e6, 1e5, 2e5); dv = ECI(100., 200., 300.)
check("ECI + ECI → ECI",     (r1 + dv) isa ECI)
check("ECI - ECI → ECI",     (r1 - r2) isa ECI)
check("scalar * ECI → ECI",  (2.0 * dv) isa ECI)
check("ECI * scalar → ECI",  (dv * 2.0) isa ECI)
check("ECI / scalar → ECI",  (dv / 2.0) isa ECI)
check("-ECI → ECI",          (-dv) isa ECI)
check("cross ECI → ECI",     cross(dv, dv) isa ECI)
check("ECI + SVector → ECI", (r1 + SVector(1.,2.,3.)) isa ECI)
check("ECI - ECI values",            r1 - r2 ≈ ECI(1e6, -1e5, -2e5))
check("ECI + ECI values",        r1 + dv ≈ ECI(7e6+100., 200., 300.))
check("scalar * ECI values",     2.0 * dv ≈ ECI(200., 400., 600.))
check("ECI + ECEF → error",  try (ECI(1.,0.,0.) + ECEF(1.,0.,0.)); false catch; true end)

# ── Rotation matrix helpers return SMatrix ────────────────────────────────────
check("ecef2eci SMatrix",  ecef2eci(0.) isa SMatrix{3,3})
check("eci2ecef SMatrix",  eci2ecef(0.) isa SMatrix{3,3})
check("ecef2enu SMatrix",  ecef2enu(0.,0.) isa SMatrix{3,3})
check("ecef2ned SMatrix",  ecef2ned(0.,0.) isa SMatrix{3,3})

# ── BodyCoord↔ECI transforms ─────────────────────────────────────────────────────
# Body frame rotated -90° about z (body→ECI active): ECI x → BodyCoord y
q90z = Quat((0., 0., -π/2))
r₀   = ECI(0., 0., 0.)
v₀   = ECI(0., 0., 0.)
ω0   = ECI(0., 0., 0.)
r_eci = ECI(1., 0., 0.)
r_bf  = BodyCoord(r_eci, q90z, r₀)
check("BodyCoord pos transform",  isapprox(r_bf, BodyCoord(0., 1., 0.), atol=1e-14))
check("BodyCoord pos roundtrip",  norm(ECI(r_bf, q90z, r₀) - r_eci) < 1e-14)
v_eci = ECI(1., 0., 0.)
v_bf  = BodyCoord(v_eci, q90z; rel=true)
check("BodyCoord vel transform",  isapprox(v_bf, BodyCoord(0., 1., 0.), atol=1e-14))
check("BodyCoord vel roundtrip",  norm(ECI(v_bf, q90z; rel=true) - v_eci) < 1e-14)
# pv: stationary ECI point at (1,0,0), identity attitude, body spinning at 1 rad/s about z
# → body-frame velocity = (0,-1,0): point drifts in -y as frame rotates counterclockwise
ω_eci  = ECI(0., 0., 1.)
pv_eci = ECI6(ECI(1.,0.,0.), ECI(0.,0.,0.))
pv_bf  = BodyCoord6(pv_eci, Quat(1.0), r₀, v₀, ω_eci)
check("BodyCoord6 vel Coriolis", isapprox(pv_bf.v, BodyCoord(0.,-1.,0.), atol=1e-14))
pv_rt  = ECI6(pv_bf, Quat(1.0), r₀, v₀, ω_eci)
check("ECI6 roundtrip",      norm(pv_rt.r - pv_eci.r) < 1e-14 && norm(pv_rt.v - pv_eci.v) < 1e-14)
# pv roundtrip with nonzero attitude (q90z) and ω
pv_bf2 = BodyCoord6(pv_eci, q90z, r₀, v₀, ω_eci)
pv_rt2 = ECI6(pv_bf2, q90z, r₀, v₀, ω_eci)
check("ECI6 roundtrip q90z", norm(pv_rt2.r - pv_eci.r) < 1e-14 && norm(pv_rt2.v - pv_eci.v) < 1e-14)

# ── Orbit ─────────────────────────────────────────────────────────────────────
o = OrbEl(7e6, 0.01, 51.6, 45.0, 90.0, 0.0)
check("OrbEl fields",     o.a == 7e6 && o.e == 0.01 && o.i == 51.6)

# ── Result ────────────────────────────────────────────────────────────────────
println("\n$pass passed, $fail failed")
fail > 0 && exit(1)
