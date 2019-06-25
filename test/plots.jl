using PyPlot
using SFGIntensities
using Statistics

plt[:style][:use]("/Users/lackner/styles/mystyle.rc")

pygui(true)


λ_vis = 532.0
λ_ir = 3450
λ_sf = SFGIntensities.sfwavelength(λ_vis, λ_ir)

β_vis = deg2rad(60)
β_ir = deg2rad(55)
β_sf = SFGIntensities.sfangle(λ_vis, λ_ir, β_vis, β_ir)

setup = SFGIntensities.Setup(
    (sf=λ_sf, vis=λ_vis, ir=λ_ir),
    (sf=β_sf, vis=β_vis, ir=β_ir),
    (sf=1.0+0.3im, vis=1.0+0.3im, ir=1.0+0.3im),
    (sf=1.44+0.3im, vis=1.44+0.3im, ir=1.44+0.3im),
    (sf=1.18+0.3im, vis=1.18+0.3im, ir=1.18+0.3im)
)

θ = range(0, stop=π/2, length=100)
θ_deg = rad2deg.(θ)
τ = deg2rad(109)
ρ = 1e-6
ψ = 0.0
ω = 3000.0
M_c = 12.0
M_b = 1.0
polarization_combinations = [:ssp, :ppp, :sps, :pss]

figure()
suptitle("psi = $ψ")

M = fill(zero(typeof(setup.n1[:sf])), length(θ), 4)
for (i, pol) in enumerate(polarization_combinations), (j, t) in enumerate(θ)
    M[j,i] = SFGIntensities.effective_susceptibility(setup, pol, :c3v, :ss, t, ψ, τ, ρ, M_c, M_b, ω; ψ_dist=:iso)
end
subplot(2,2,1)
title("C3V SS")
plot(SFGIntensities.sfintensity(M, setup.β[:sf]))

for (i, pol) in enumerate(polarization_combinations), (j, t) in enumerate(θ)
    M[j,i] = SFGIntensities.effective_susceptibility(setup, pol, :c2v, :ss, t, ψ, τ, ρ, M_c, M_b, ω; ψ_dist=:iso)
end
subplot(2,2,2)
title("C2V SS")
plot(SFGIntensities.sfintensity(M, setup.β[:sf]))

for (i, pol) in enumerate(polarization_combinations), (j, t) in enumerate(θ)
    M[j,i] = SFGIntensities.effective_susceptibility(setup, pol, :c2v, :as, t, ψ, τ, ρ, M_c, M_b, ω; ψ_dist=:iso)
end
subplot(2,2,4)
title("C2V AS")
plot(SFGIntensities.sfintensity(M, setup.β[:sf]))

for (i, pol) in enumerate(polarization_combinations), (j, t) in enumerate(θ)
    M[j,i] = SFGIntensities.effective_susceptibility(setup, pol, :c3v, :as, t, ψ, τ, ρ, M_c, M_b, ω; ψ_dist=:iso)
end
subplot(2,2,3)
title("C3V AS")
plot(SFGIntensities.sfintensity(M, setup.β[:sf]))

# blub
