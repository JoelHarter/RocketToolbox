# UNITS
const Units = let
     # all defined in terms of SI base units and radians unless otherwise stated
    rev = 2π # [rad] full revolution
    deg = rev/360 # [rad] degree of arc
    arcmin = deg/60 # [rad] arc minute
    arcsec = arcmin/60 # [rad] arc second
    minute = 60.0 # [s]
    hour = 60*minute # [s]
    day = 24*hour # [s] Julian day
    year = 365.25*day # [s] Julian year
    cy = 100*year # [s] Julian century
    in = 0.0254 # [m] inch
    Nmi = 1852.0 # [m] nautical mile
    au  = 149_597_870_700.0 # <IAU> [m] astronomical unit
    knot = Nmi/hour # [m/s] knot (nautical mile per hour)
    lb = 0.45359237 # [kg] pound-mass
    t = 1e3 # [kg] tonne (metric ton)
    l = 1e-3 # [m³] liter
    atm = 101325.0 # <CGPM> [Pa] standard atmosphere
    bar = 1e5 # [Pa] bar
    g₀  = 9.80665 # <CGPM> [m/s²] standard gravity
    lbf = lb * g₀ # [N] pound-force
    °C₀ = 273.15 # [K] 0 degrees Celsius
    °R = 5/9 # [K] degree Rankine
    °F₀ = °C₀ - 32°R # [K] 0 degrees Fahrenheit
    (; rev, deg, arcmin, arcsec, minute, hour, day, year, cy, in, Nmi, au, knot, lb, t, l, atm, bar, g₀, lbf, °C₀, °R, °F₀)
end

# PHYSICAL CONSTANTS
const Phys = let
    G  = 6.67428e-11 #      <IAU>                   [N·m²/kg²]  universal gravitational constant
    c  = 299_792_458.0 #    <SI>                    [m/s]       speed of light in vacuum
    h  = 6.62607015e-34 #   <SI>                    [J/Hz]      Planck's constant
    h̄  = h / (2π) #                                 [J·s/rad]   reduced Planck's constant
    k  = 1.380649e-23 #     <SI>                    [J/K]       Boltzmann's constant
    NA = 6.02214076e23 #    <SI>                    [1/mol]     Avogadro's constant
    R  = NA * k #                                   [J/(mol·K)] universal gas constant
    (; G, c, h, h̄, k, NA, R)
end

# EARTH AND ASTRONOMICAL PARAMETERS
const Earth = let
    a     = 6_378_137.0                         # [m] <WGS> semi-major axis (equatorial radius)
    a²    = sqr(a)                              # [m²]
    f⁻¹   = 298.257223563                       # [1] <WGS> inverse flattening
    f     = /(f⁻¹)                              # [1] flattening
    b     = (1.0 - f) * a                       # [m] semi-minor axis (polar radius)
    b²    = sqr(b)                              # [m²]
    a⁻²b² = b² / a²                             # [1] ≡ 1-e²
    e²    = 1.0 - a⁻²b²                         # [1]
    e     = sqrt(e²)                            # [1] eccentricity
    a⁻²   = inv(a²)                             # [1/m²] inverse of a squared
    a⁻⁴b² = a⁻²b² * a⁻²                         # [1/m²]
    e⁴    = sqr(e²)                             # [1] eccentricity to the 4th power
    θ₀    = 0.77905727326403(Units.rev)         # [rad] <IAU> rotation angle at J2000 epoch
    Ω     = 72.92115e-6                         # [rad/s] <WGS> earth's rotational speed
    ε     = 84381.406(Units.arcsec)             # [rad] <IAU> obliquity of ecliptic (axial tilt)
    μ     = 3986004.415e8                       # [m³/s²] <EGM> gravitational constant, including mass of atmosphere
    μS    = 1.32712442099e20                    # [m³/s²] <IAU> gravitational constant of the sun
    g₀    = Units.g₀                            # [m/s²] <CGPM> standard acceleration of gravity set by ISO 80000 (international standard)
    atm   = Units.atm                           # [Pa] <CGPM> standard atmosphere set by ISO 80000 (international standard)
    au    = Units.au                            # [m] semi-major axis of earth's orbit
    p     = 5028.796195(Units.arcsec/Units.cy)  # [rad/s] <IAU> axial precession rate
    R₁    = (2a + b)/3                          # [m] arithmetic mean radius
    R₂    = sqrt((a² + b² * atanh(e) / e) * 0.5) # [m] authalic (equal area) radius
    R₃    = cbrt(a² * b)                        # [m] volumetric (sphere with equal volume) radius
    M     = μ / Phys.G;                         # [kg] earth mass; NOTE: current best estimate of earth's mass is (5.9722±0.0006)*10^24 kg. The value assigned here keeps things consistent with the defined μ and G
    C̅₂₀   = -0.484165143790815e-3               # [1] <EGM> 
    J₂    = C̅₂₀*-sqrt(5)                        # [1] second zonal harmonic coefficient
    K_J₂  = 1.5J₂*μ*a²                         # [m⁵/s²] J₂ perturbation constant
    ωS    = sqrt(au^3\(μS+μ))                   # [rad/s] mean angular velocity of earth around sun
    yearS = 2π/ωS                               # [s] siderial year, the length of time it takes the sun to return to the same relative position with respect to the fixed stars, taking the earth's axial precession into account
    yearT = 2π/(ωS + p)                         # [s] tropical year, the length of time it takes the earth to go from one vernal equinox to the next
    (; a, b, f, f⁻¹, e, e², a², b², a⁻², a⁻²b², a⁻⁴b², e⁴, R₁, R₂, R₃, M, μ, C̅₂₀, J₂, K_J₂, Ω, θ₀, ε, p, μS, au, yearT, yearS, g₀, atm)
end


# REFERENCES
#   <SI> Système international d'unités 2019
#       https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf"
#   <WGS> World Geodetic System 1984
#      https://apps.dtic.mil/sti/pdfs/ADA280358.pdf
#   <EGM> Earth Gravitational Model 2008
#      https://earth-info.nga.mil/php/download.php?file=egm-08spherical
#   <IAU> International Astronomical Union 2012
#      https://aa.usno.navy.mil/downloads/publications/Constants_2021.pdf
#   <CGPM> Conférence générale des poids et mesures 1901
#      https://www.bipm.org/documents/20126/38096616/3rd+CGPM+%281901%29/9993082b-5a57-0b7e-1dfd-ea3245536d94
