using PyPlot
using SFGIntensities

## ----
N = 200

angles_old = [53, 42]
angles_new = [42, 46]

# Wavelengths
λ_vis = 512
λ_ir = 7500
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
β_ir = range(0, stop=π/2, length=N)
β_vis = range(0, stop=π/2, length=N)
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
vmin = [1e23, 1e24, 1e24]
vmax = [1e25, 1e28, 1e28]
figure()
for (i, pol) in [:ssp, :sps, :pss] |> enumerate
    χ_eff = makemat(setup, β_vis, β_ir, β_sf, pol)
    global intensity = sfintensity(χ_eff, setup) |> real

    subplot(1,3,i)
    title("$pol")
    pcolor(rad2deg.(β_ir), rad2deg.(β_vis), intensity, cmap=ColorMap("magma"),
        norm=plt[:matplotlib][:colors][:LogNorm](vmin=minimum(intensity), vmax=maximum(intensity)),
        vmin=vmin[i], vmax=vmax[i])

    colorbar(extend="both", label="Fresnel Enhancement")
    xlabel("Vis Angle")

    angles_old_idx = N ÷ 90 * angles_old
    angles_new_idx = N ÷ 90 * angles_new

    @printf "\n"
    @printf "%s" "----- $pol -----\n"
    @printf "%s:\t%2.2g\n" "old" intensity[angles_old[2], angles_old[1]]
    @printf "%s:\t%2.2g\n" "new" intensity[angles_new[2], angles_new[1]]
    @printf "%s:\t%5.2f\n" "ratio" (intensity[angles_new[2], angles_new[1]] ./ intensity[angles_old[2], angles_old[1]])

    scatter(angles_old[1], angles_old[2], marker="o", color="k")
    scatter(angles_new[1], angles_new[2], marker="x", color="k")
end
subplot(131)
ylabel("IR Angle")
# tight_layout()

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
             vmin=1e20, vmax=1e28)
    ylabel("IR Angle")
    xlabel("Vis Angle")
    colorbar(label="Fresnel Enhancement", extend="both")

    angles_old_idx = N ÷ 90 * angles_old
    angles_new_idx = N ÷ 90 * angles_new

    @printf "\n"
    @printf "%s" "----- $geo -----\n"
    @printf "%s:\t%2.2g\n" "old" intensity[angles_old[2], angles_old[1]]
    @printf "%s:\t%2.2g\n" "new" intensity[angles_new[2], angles_new[1]]
    @printf "%s:\t%5.2f\n" "ratio" (intensity[angles_new[2], angles_new[1]] ./ intensity[angles_old[2], angles_old[1]])

    scatter(angles_old[1], angles_old[2], marker="o", color="k")
    scatter(angles_new[1], angles_new[2], marker="x", color="k")
end
tight_layout()
