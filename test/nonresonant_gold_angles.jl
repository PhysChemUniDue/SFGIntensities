using SFGIntensities
using PyPlot

# Wavelengths
λ_vis = 512
λ_ir = 3400
λ_sf = sfwavelength(λ_vis, λ_ir)

# Refractive Indices
# Air
n1_vis = 1.0 + 0im
n1_sf = 1.0 + 0im
n1_ir = 1.0 + 0im
# Gold @ 525 nm
n2_vis = 0.48130 + 2.1355im
n2_sf = 1.2070 + 1.7202im
n2_ir = 6.2925 + 46.619im
# Interface

# Incident Angles
β_vis = deg2rad(45)
β_ir = range(deg2rad(-80), deg2rad(80), length=199)
β_ir = deg2rad.([-45, 45])
β_sf = sfangle.(λ_ir, λ_vis, β_ir, β_vis)

s = SFGIntensities.Setup(
    (sf=λ_sf, vis=λ_vis, ir=λ_ir),
    (sf=β_sf[1], vis=β_vis, ir=β_ir[1]),
    (sf=n1_sf, vis=n1_vis, ir=n1_ir),
    (sf=n2_sf, vis=n2_vis, ir=n2_ir),
    (sf=n2_sf, vis=n2_vis, ir=n2_ir)
)

χ_xxz = 0.0
χ_xzx = 0.0
χ_zxx = 0.0
χ_zzz = 1.0

χ_eff = Array{ComplexF64}(undef, length(β_ir))
I_sf = Array{Float64}(undef, length(β_ir))
for (i, b) in β_ir |> enumerate
    β_sf = sfangle.(λ_ir, λ_vis, b, β_vis)
    s.β = (sf=β_sf, vis=s.β.vis, ir=b)
    χ_eff[i] = effective_susceptibility(s, :ppp, :zzz)
    I_sf[i] = sfintensity(χ_eff[i], s)
end

χ_eff ./= maximum(real(χ_eff))
I_sf ./= maximum(I_sf)

# I1 = I_sf[1:100]
# I2 = I_sf[100:end]

##--PLOT
pygui(true)
figure()
title("Nonresonant Signal Au")
plot(rad2deg.(β_ir), χ_eff, label="Effective Susceptibility")
plot(rad2deg.(β_ir), I_sf, label="SF Signal Intensity")
xlabel("IR Angle (deg)")

figure()
plot(I1 ./ reverse(I2))
