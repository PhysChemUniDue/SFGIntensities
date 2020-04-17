using SFGIntensities
using Plots

## SETUP
# Wavelengths
λ_vis = 515
λ_ir = 3350
λ_sf = sfwavelength(λ_vis, λ_ir)

# Refractive Indices
# Air
n1_vis = 1.0 + 0im
n1_sf = 1.0 + 0im
n1_ir = 1.0 + 0im
# BK7 Glass
n2_vis = 1.52 + 0im
n2_sf = 1.52 + 0im
n2_ir = 1.34 + 0im # guessed value
# Interface
ni_vis = ninterface(n2_vis)
ni_sf  = ninterface(n2_sf )
ni_ir  = ninterface(n2_ir )

# Incident Angles
β_ir = deg2rad(60)
β_vis = deg2rad(45)
β_sf = sfangle(λ_vis, λ_ir, β_vis, β_ir, n1_sf, n1_vis, n1_ir)

setup = SFGIntensities.Setup(
    (sf=λ_sf, vis=λ_vis, ir=λ_ir),
    (sf=β_sf, vis=β_vis, ir=β_ir),
    (sf=n1_sf, vis=n1_vis, ir=n1_ir),
    (sf=n2_sf, vis=n2_vis, ir=n2_ir),
    (sf=n2_sf, vis=n2_vis, ir=n2_ir)
)

## Block
ρ = 0.03 #guess
pol = [:ssp, :ppp]
θavg = deg2rad(0) #guess
symmetry = :c2v
mode = [:ss, :as]
ψavg = rad2deg(90)
ψ_dist = :fixed
τ = deg2rad(109.5)
M_c = 12
M_b = 1.0
ω = [2880.0, 2960.0]

χeff = [
    effective_susceptibility(setup, _p, symmetry, mode[j], θavg, ψavg, τ,
                             ρ, M_c, M_b, ω[j])
    for _p in pol, j in eachindex(mode)
]

φ = angle.(χeff)

println("χ effective (SS) | (AS):\nssp\nppp")
display(χeff)
println("Phase Angle: (SS) | (AS):\nssp\nppp")
display(φ)
