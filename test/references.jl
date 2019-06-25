using PyPlot
using SFGIntensities

## ----
# First Reference is taken from https://pubs.acs.org/doi/10.1021/jp306838a

# Wavelengths
λ_vis = 800
λ_ir = 3350
λ_sf = sfwavelength(λ_vis, λ_ir)

# Refractive Indices
# Silica
n1_vis = 1.453 + 0im
n1_sf = 1.455 + 0im
n1_ir = 1.414 + 0im
# Water
n2_vis = 1.33 + 0im
n2_sf = 1.33 + 0im
n2_ir = 1.34 + 0.08im
# Interface

# Incident Angles
β_ir = range(0, stop=π/2, length=200)
β_vis = range(0, stop=π/2, length=200)
β_sf = [sfangle(λ_vis, λ_ir, v, i, n1_sf, n1_vis, n1_ir) for
             v in β_vis, i in β_ir]

setup = SFGIntensities.Setup(
    (sf=λ_sf, vis=λ_vis, ir=λ_ir),
    (sf=β_sf[1,1], vis=β_vis[1], ir=β_ir[1]),
    (sf=n1_sf, vis=n1_vis, ir=n1_ir),
    (sf=n2_sf, vis=n2_vis, ir=n2_ir),
    (sf=n2_sf, vis=n2_vis, ir=n2_ir)    )

function makemat(setup, β_vis, β_ir, β_sf, pol, geo=:none)
    χ_eff = Array{ComplexF64,2}(undef, length(β_ir), length(β_vis))
    for v = 1:length(β_vis), i = 1:length(β_ir)
        setup.β = (sf=β_sf[v,i], vis=β_vis[v], ir=β_ir[i])
        χ_eff[i,v] = effective_susceptibility(setup, pol, geo)
    end
    χ_eff
end

## -- ssp sps pss polarizations angle dependence
figure()
for (i, pol) in [:ssp, :sps, :pss] |> enumerate
    χ_eff = makemat(setup, β_vis, β_ir, β_sf, pol)
    intensity = sfintensity(χ_eff, setup) |> real

    subplot(1,3,i)
    title("$pol")
    pcolor(rad2deg.(β_ir), rad2deg.(β_vis), intensity, cmap=ColorMap("magma"),
        norm=plt[:matplotlib][:colors][:LogNorm](vmin=minimum(intensity), vmax=maximum(intensity)),
        vmin=1e-3, vmax=1e2)
    ylabel("IR Angle")
    xlabel("Vis Angle")
end
tight_layout()
colorbar(extend="both", label="Fresnel Enhancement")

## --
figure()
for (i, geo) in [:zzz, :xzx, :zxx, :xxz] |> enumerate
    χ_eff = makemat(setup, β_vis, β_ir, β_sf, :ppp, geo)
    global intensity = sfintensity(χ_eff, setup) |> real
    intensity[intensity .== 0.0] .= 1e-50

    subplot(2,2,i)
    title("$geo")
    pcolor(rad2deg.(β_ir), rad2deg.(β_vis), intensity, cmap=ColorMap("magma"),
             norm=plt[:matplotlib][:colors][:LogNorm](vmin=minimum(intensity), vmax=maximum(intensity)),
             vmin=1e-3, vmax=1e2)
    ylabel("IR Angle")
    xlabel("Vis Angle")
end
colorbar(label="Fresnel Enhancement", extend="both")
tight_layout()
