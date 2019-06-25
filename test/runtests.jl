using Test
using SFGIntensities
using Plots

reference_image_path = joinpath(@__DIR__, "plots")
!isdir(reference_image_path) && mkdir(reference_image_path)

## TEST AGAINST http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=20108300&retmode=ref&cmd=prlinks

# Gold Surface @ 1350 1/cm

λ_vis = 1 / 18797 * 1e7
λ_ir = 1 / 1350.0 * 1e7
λ_sf = sfwavelength(λ_vis, λ_ir)
β_vis = deg2rad(55)
β_ir = deg2rad(65)
β_sf = sfangle(λ_vis, λ_ir, β_vis, β_ir)

@test round(λ_sf, digits=1) == 496.4
@test round(rad2deg(β_sf), digits=2) == 55.59

setup_Au_NO2 = SFGIntensities.Setup(
    (sf=λ_sf, vis=λ_vis, ir=λ_ir),
    (sf=β_sf, vis=β_vis, ir=β_ir),
    (sf=1.0, vis=1.0, ir=1.0),
    (sf=0.910 + 1.850im, vis=0.465 + 2.400im, ir=7.477 + 43.600im),
    (sf=1.0, vis=1.0, ir=1.0)
)

println("\nFresnel Factors Au Surface @ 1350 1/cm")
pol = [:ppp, :ppp, :ppp, :ppp, :ssp, :sps, :pss]
geo = [:zzz, :xxz, :xzx, :zxx, :xxz, :xzx, :zxx]
ref_val = [4.7e-0, 4.9e-1, 1.1e-3, 7.9e-4, 1.2e-1, 1.1e-4, 6.0e-5]
for (p, g, r) in zip(pol, geo, ref_val)
    val = effective_susceptibility(setup_Au_NO2, p, g)
    println("$p - $g - $(round(abs2(val), sigdigits=3)) -- $r")
    # @test round(val, sigdigits=2) == r
end


# Platinum Surface @ 1350 1/cm

setup_Pt_NO2 = SFGIntensities.Setup(
    (sf=λ_sf, vis=λ_vis, ir=λ_ir),
    (sf=β_sf, vis=β_vis, ir=β_ir),
    (sf=1.0, vis=1.0, ir=1.0),
    (sf=1.960 + 3.420im, vis=2.080 + 3.630im, ir=6.961 + 27.700im),
    (sf=1.0, vis=1.0, ir=1.0)
)

println("\nFresnel Factors Pt Surface @ 1350 1/cm")
pol = [:ppp, :ppp, :ppp, :ppp, :ssp, :sps, :pss]
geo = [:zzz, :xxz, :xzx, :zxx, :xxz, :xzx, :zxx]
ref_val = [8.5e-0, 7.2e-2, 1.2e-3, 1.1e-3, 1.4e-2, 1.0e-4, 9.0e-5]
for (p, g, r) in zip(pol, geo, ref_val)
    val = effective_susceptibility(setup_Pt_NO2, p, g)
    println("$p - $g - $(round(abs2(val), sigdigits=3)) -- $r")
    # @test round(val, sigdigits=2) == r
end


# Platinum Surface @ 3000 1/cm

λ_vis = 1 / 18797 * 1e7
λ_ir = 1 / 3000.0 * 1e7
λ_sf = sfwavelength(λ_vis, λ_ir)
β_vis = deg2rad(55)
β_ir = deg2rad(65)
β_sf = sfangle(λ_vis, λ_ir, β_vis, β_ir)

@test round(λ_sf, digits=1) == 458.8
@test round(rad2deg(β_sf), digits=2) == 56.22

setup_Pt_CH3 = SFGIntensities.Setup(
    (sf=λ_sf, vis=λ_vis, ir=λ_ir),
    (sf=β_sf, vis=β_vis, ir=β_ir),
    (sf=1.0, vis=1.0, ir=1.0),
    (sf=1.870 + 3.200im, vis=2.080 + 3.630im, ir=3.10 + 12.800im),
    (sf=1.0, vis=1.0, ir=1.0)
)

println("\nFresnel Factors Pt Surface @ 3000 1/cm")
pol = [:ppp, :ppp, :ppp, :ppp, :ssp, :sps, :pss]
geo = [:zzz, :xxz, :xzx, :zxx, :xxz, :xzx, :zxx]
ref_val = [7.7e0, 7.3e-2, 5.9e-3, 4.7e-3, 1.4e-2, 5.1e-4, 4.1e-4]
for (p, g, r) in zip(pol, geo, ref_val)
    val = effective_susceptibility(setup_Pt_CH3, p, g)
    println("$p - $g - $(round(abs2(val), sigdigits=3)) -- $r")
    # @test round(val, sigdigits=2) == r
end

# NO2 r+ on Au ψ angle dependence
ψ_values = [0.0, π/2]
θ_values = range(0.0, stop=π/2, length=100)
τ = deg2rad(122.0)
# raman depolarization ratio
ρ = 0.17
M_c = 14.0
M_b = 16.0
ω = 1 / 1350 * 1e7

p = plot()
for ψ in ψ_values
    χ_eff = [effective_susceptibility(setup_Au_NO2, :ppp, :c2v, :ss, θ,
                                      ψ, τ, ρ, M_c, M_b, ω)
                                      for θ in θ_values]
    I_sf = sfintensity(χ_eff)
    p = plot!(θ_values, I_sf, label="ψ = $(rad2deg(ψ))̊")
end
savefig(p, joinpath(reference_image_path, "NO2_ppp.png"))


p = plot()
for ψ in ψ_values
    χ_eff = [effective_susceptibility(setup_Au_NO2, :ssp, :c2v, :ss, θ,
                                      ψ, τ, ρ, M_c, M_b, ω)
                                      for θ in θ_values]
    I_sf = sfintensity(χ_eff)
    p = plot!(θ_values, I_sf, label="ψ = $(rad2deg(ψ))̊")
end
savefig(p, joinpath(reference_image_path, "NO2_ssp.png"))


τ = deg2rad(122.0)
# raman depolarization ratio
ρ = 0.0005
M_c = 12.0
M_b = 1.0
ω = 1 / 3000 * 1e7

p = plot()
for ψ in ψ_values
    χ_eff = [effective_susceptibility(setup_Pt_CH3, :ppp, :c3v, :ss, θ,
                                      ψ, τ, ρ, M_c, M_b, ω)
                                      for θ in θ_values]
    I_sf = sfintensity(χ_eff)
    p = plot!(θ_values, I_sf, label="ψ = $(rad2deg(ψ))̊")
end
savefig(p, joinpath(reference_image_path, "Pt_CH3_ppp.png"))

## STUFF

λ_vis = 532.0
λ_ir = 4762.0
λ_sf = SFGIntensities.sfwavelength(λ_vis, λ_ir)

β_vis = deg2rad(62)
β_ir = deg2rad(53)
β_sf = SFGIntensities.sfangle(λ_vis, λ_ir, β_vis, β_ir)



setup = SFGIntensities.Setup(
    (sf=λ_sf, vis=λ_vis, ir=λ_ir),
    (sf=β_sf, vis=β_vis, ir=β_ir),
    (sf=1.0, vis=1.0, ir=1.0),
    (sf=1.30, vis=1.33, ir=1.33),
    (sf=1.15, vis=1.15, ir=1.15)
)

# setup = SFGIntensities.Setup(
#     (vis=deg2rad(55), ir=deg2rad(38)),
#     (sf=1.0, vis=1.0, ir=1.0),
#     (sf=1.3, vis=1.3, ir=1.5),
#     (sf=1.2, vis=1.2, ir=1.2)
# )

# setup = SFGIntensities.Setup(
#     (deg2rad(55), deg2rad(38)),
#     (1.0, 1.0, 1.0),
#     (1.3, 1.3, 1.5),
#     (1.2, 1.2, 1.2)
# )

for ω in [:sf, :vis, :ir]
    SFGIntensities.l_xx(setup, ω)
    SFGIntensities.l_yy(setup, ω)
    SFGIntensities.l_zz(setup, ω)
end

θ = deg2rad(10)
ψ = 0.0
τ = deg2rad(109.5)
ρ = 0.005
M_c = 12.0
M_b = 1.0
ω = 3000.0
for geo in [:xzx, :zxx, :yzy, :zyy, :xxz, :yyz, :zzz], sym in [:c2v]
    SFGIntensities.susceptibility(geo, sym, :ss, θ, ψ, τ, ρ, M_c, M_b, ω)
end

θ = deg2rad(10)

for pol in [:ssp, :pss, :sps, :ppp], sym in [:c2v], mode in [:ss, :as]
    SFGIntensities.effective_susceptibility(setup, pol, sym, mode, θ, ψ, τ, ρ, M_c, M_b, ω)
end

# Water molecule
τ = deg2rad(104.5)
ρ = 0.005
M_O = 16.0
M_H = 1.0
ω = 3000.0
β_aac = hyperpolarizability(:aac, :c2v, τ, ρ, M_O, M_H, ω)
β_bbc = hyperpolarizability(:bbc, :c2v, τ, ρ, M_O, M_H, ω)
β_ccc = hyperpolarizability(:ccc, :c2v, τ, ρ, M_O, M_H, ω)
@test β_aac + β_bbc - 2β_ccc <= (β_aac + β_bbc + β_ccc) / 30

θ = range(0, stop=π/2, length=100)
θ_deg = rad2deg.(θ)
χ_sps = [effective_susceptibility(setup, :sps, :c2v, :ss, θx, ψ, τ, ρ, M_c, M_b, ω) for θx in θ]
χ_ssp = [effective_susceptibility(setup, :ssp, :c2v, :ss, θx, ψ, τ, ρ, M_c, M_b, ω) for θx in θ]
χ_ppp = [effective_susceptibility(setup, :ppp, :c2v, :ss, θx, ψ, τ, ρ, M_c, M_b, ω) for θx in θ]
χ_pss = [effective_susceptibility(setup, :pss, :c2v, :ss, θx, ψ, τ, ρ, M_c, M_b, ω) for θx in θ]


ψ = 0.0
for pol in (:ssp, :ppp, :sps, :pss), mode in (:ss, :as)
    [effective_susceptibility(setup, pol, :c2v, mode, θx, ψ, τ, ρ, M_c, M_b, ω) for θx in θ]
end


p = plot(θ_deg, χ_ssp, label="ssp")
plot!(θ_deg, χ_sps * 100, label="sps × 100")
plot!(θ_deg, χ_ppp, label="ppp")
plot!(θ_deg, χ_pss * 100, label="pss × 100")

savefig(p, joinpath(reference_image_path, "polarization_combinations.png"))
