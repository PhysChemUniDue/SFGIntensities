using SFGIntensities
using ProgressMeter
using PyPlot

λ_vis = 532.0
λ_ir = 3450.0
λ_sf = SFGIntensities.sfwavelength(λ_vis, λ_ir)

β_vis = range(0, stop=π/2, length=100) |> collect
β_ir = range(π/3, stop=π/2, length=100) |> collect
β_sf = sfangle(λ_vis, λ_ir, β_vis, β_ir)

n2_sf = 1.4821 + 1.7888im
n2_vis = 0.58076 + 1.8032im
n2_ir = 0.54182 + 19.716im

# NH2
M_c = 14.0
M_b = 1.0
ρ = 0.0005 #???
τ = deg2rad(109.5)

ω = 3444.0

θ_mean = deg2rad(20)
σ = deg2rad(20)
ψ_mean = deg2rad(0)
σ_ψ = deg2rad(20)

χ_eff_ss = Array{Complex{Float64},2}(undef, length(β_vis), length(β_ir))
χ_eff_as = copy(χ_eff_ss)

@showprogress 1 "Dada..." for i in 1:length(β_vis), j in 1:length(β_ir)
    setup = SFGIntensities.Setup(
        (sf=λ_sf, vis=λ_vis, ir=λ_ir),
        (sf=β_sf[i,j], vis=β_vis[i], ir=β_ir[j]),
        (sf=1.0+0.0im, vis=1.0+0.0im, ir=1.0+0.0im),
        (sf=n2_sf, vis=n2_vis, ir=n2_ir),
        (sf=n2_sf, vis=n2_vis, ir=n2_ir)
    )

    χ_eff_ss[i,j] = effective_susceptibility(setup, :ssp, :c2v, :ss, θ_mean, ψ_mean, τ, ρ, M_c, M_b, ω;
                                                            ψ_dist=:iso, θ_dist=:δ, σ_θ=σ, σ_ψ=σ_ψ, n=10)
    χ_eff_as[i,j] = effective_susceptibility(setup, :ssp, :c2v, :as, θ_mean, ψ_mean, τ, ρ, M_c, M_b, ω;
                                                            ψ_dist=:iso, θ_dist=:δ, σ_θ=σ, σ_ψ=σ_ψ, n=10)
end

int_ss = sfintensity(χ_eff_ss, β_sf)
int_as = sfintensity(χ_eff_as, β_sf)

figure()
title("NH2 Symmetric Stretch")
pcolor(rad2deg.(β_ir), rad2deg.(β_vis), int_ss)
xlabel("IR Angle")
ylabel("Vis Angle")
colorbar()

figure()
title("NH2 Asymmetric Stretch")
pcolor(rad2deg.(β_ir), rad2deg.(β_vis), int_as)
xlabel("IR Angle")
ylabel("Vis Angle")
colorbar()
