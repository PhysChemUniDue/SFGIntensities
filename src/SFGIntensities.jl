module SFGIntensities

using Roots
using Distributions
using Statistics

mutable struct Setup
    λ
    β
    n1
    n2
    n′
end

export effective_susceptibility,
       sfwavelength,
       sfintensity,
       hyperpolarizability,
       ninterface,
       sfangle

refracangle(β, n1, n2) = @. asin(n1 / n2 * sin(β))

sfwavelength(λ_vis, λ_ir) = @. 1 / (1 / λ_vis + 1 / λ_ir)

function sfangle(λ_vis, λ_ir, β_vis::T, β_ir::T, n1_sf=1.0, n1_vis=1.0, n1_ir=1.0) where T <: Real
    λ_sf = sfwavelength(λ_vis, λ_ir)
    β_sf = @. asin(λ_sf / λ_vis * n1_vis / n1_sf * sin(β_vis) + λ_sf / λ_ir * n1_ir / n1_sf * sin(β_ir))
    real(β_sf)
end

function sfangle(λ_vis, λ_ir, β_vis::Array{T}, β_ir::Array{T}, n1_sf=1.0, n1_vis=1.0, n1_ir=1.0) where T <: Real
    β_sf = Array{T,2}(undef, length(β_vis), length(β_ir))
    for i ∈ 1:length(β_vis), j ∈ 1:length(β_ir)
        β_sf[i,j] = sfangle(λ_vis, λ_ir, β_vis[i], β_ir[j], n1_sf, n1_vis, n1_ir)
    end
    β_sf
end

function sfintensity(χ_eff, β_sf, ω_sf=1.0, n1_sf=1.0, n1_vis=1.0, n1_ir=1.0, I_vis=1.0, I_ir=1.0)
    c0 = 3e8
    @. 8π^3 * ω_sf^2 * sec(β_sf)^2 / (c0 * n1_sf * n1_vis * n1_ir) * abs(χ_eff)^2 * I_vis * I_ir
end
sfintensity(χ_eff, s::Setup, I_vis=1.0, I_ir=1.0) = sfintensity(χ_eff, s.β.sf, s.λ.sf, s.n1.sf, s.n1.vis, s.n1.ir, I_vis, I_ir)
sfintensity(χ_eff) = abs2.(χ_eff)

function l_xx(β, n1, n2)
    γ = refracangle(β, n1, n2)
    @. 2n1 * cos(γ) / (n1 * cos(γ) + n2 * cos(β))
end
l_xx(s, ω::Symbol) = l_xx(s.β[ω], s.n1[ω], s.n2[ω])

function l_yy(β, n1, n2)
    γ = refracangle(β, n1, n2)
    @. 2n1 * cos(β) / (n1 * cos(β) + n2 * cos(γ))
end
l_yy(s, ω::Symbol) = l_yy(s.β[ω], s.n1[ω], s.n2[ω])

function l_zz(β, n1, n2, n′)
    γ = refracangle(β, n1, n2)
    @. 2n2 * cos(β) / (n1 * cos(γ) + n2 * cos(β)) * (n1 / n′)^2
end
l_zz(s, ω::Symbol) = l_zz(s.β[ω], s.n1[ω], s.n2[ω], s.n′[ω])

"""
r is the single bond polarizability derivate ratio.
τ is the bond angle
"""
function ρ2r(symmetry::Symbol, ρ, τ)
    if symmetry == :c2v
        f = r -> 3 / (4 + 20 * ((1 + 2r)^2 / ((1 - r)^2 * (1 + 3cos(τ)^2)))) - ρ
        return find_zero(f, (0, 1))
    elseif symmetry == :c3v
        f = r -> 3 / (4 + 20 * ((1 + 2r) / ((1 - r) * (1 - 3cos(τ)^2) ))^2) - ρ
        return find_zero(f, (0, 1))
    else
        error("Unknown symmetry $symmetry.")
    end
end


"""
parameter
"""
function hyperpolarizability(geo::Symbol, symmetry::Symbol, τ, ρ, M_c, M_b, ω; β=1.0)

    r = ρ2r(symmetry, ρ, τ)

    if symmetry == :c2v

        G_a = @. (1 + cos(τ)) / M_c + 1 / M_b
        G_b = @. (1 - cos(τ)) / M_c + 1 / M_b

        if geo == :aac
            β = @. G_a * β / ω * ((1 + r) - (1 - r) * cos(τ)) * cos(τ/2)
        elseif geo == :bbc
            β = @. 2 * G_a * β / ω * r * cos(τ/2)
        elseif geo == :ccc
            β = @. G_a * β / ω * ((1 + r) + (1 - r) * cos(τ)) * cos(τ/2)
        elseif any(geo .== (:aca, :caa))
            β = @. G_b * β / ω * ((1 - r) * sin(τ)) * sin(τ/2)
        elseif any(geo .== (:bcb, :cbb))
            β = 0
        end
    end

    if symmetry == :c3v

        G_s = @. (1 + 2cos(τ)) / M_c + 1 / M_b
        G_d = @. (1 - cos(τ)) / M_c + 1 / M_b

        if any(geo .== (:aac, :bbc))
            β = @. -3/2 * G_s * β / ω * ((1 + r) - (1 -r) * cos(τ)^2) * cos(τ)
        elseif geo == :ccc
            β = @. -3 * G_s * β / ω * (r + (1 - r) * cos(τ)^2) * cos(τ)
        elseif any(geo .== (:bcb, :cbb, :caa, :aca))
            β = @. -3/2 * G_d * β / ω * (1 - r) * sin(τ)^2 * cos(τ)
        elseif any(geo .== (:bba, :abb, :bab))
            β = @. -3/4 * G_d * β / ω * (1 - r) * sin(τ)^3
        elseif geo == :aaa
            β = @. 3/4 * G_d * β / ω * (1 - r) * sin(τ)^3
        end
    end

    β
end

# r -> R in c3v
function susceptibility(geo::Symbol, symmetry::Symbol, mode::Symbol, θ, ψ, τ, ρ, M_c, M_b, ω; N=1.0)

    supported_symmetries = [:c2v, :c3v]
    supported_modes = [:ss, :as]

    @assert any(symmetry .== supported_symmetries) "Unknown symmety $symmetry. Available symmetries are: $supported_symmetries."
    @assert any(mode .== [:ss, :as]) "Unknown mode $mode. Available modes are: $supported_modes."

    β_aac = hyperpolarizability(:aac, symmetry, τ, ρ, M_c, M_b, ω)
    β_bbc = hyperpolarizability(:bbc, symmetry, τ, ρ, M_c, M_b, ω)
    β_ccc = hyperpolarizability(:ccc, symmetry, τ, ρ, M_c, M_b, ω)
    β_aca = hyperpolarizability(:aca, symmetry, τ, ρ, M_c, M_b, ω)

    # C3V symmetric stretch
    if symmetry == :c3v && mode == :ss
        R = β_aac ./ β_ccc
        if any(geo .== (:xzx, :zxx, :yzy, :zyy))
            χ = @. 1/2 * N * β_ccc * (cos(θ) - cos(θ)^3) * (1 - R)
        elseif any(geo .== (:xxz, :yyz))
            χ = @. 1/2 * N * β_ccc * (cos(θ) * (1 + R) - cos(θ)^3 * (1 - R))
        elseif geo == :zzz
            χ = @.       N * β_ccc * (R * cos(θ) + cos(θ)^3 * (1 - R))
        else
            error("$geo is not a valid combination.")
        end
    end

    # C3V asymmetric stretch
    if symmetry == :c3v && mode == :as
        if any(geo .== (:xxz, :yyz))
            χ = @.    -N * β_aca * (cos(θ) - cos(θ)^3)
        elseif any(geo .== (:xzx, :zxx, :yzy, :zyy))
            χ = @.     N * β_aca * cos(θ)^3
        elseif geo == :zzz
            χ = @. 2 * N * β_aca * (cos(θ) - cos(θ)^3)
        end
    end

    # C2V symmetric stretch
    if symmetry == :c2v && mode == :ss
        if any(geo .== (:xxz, :yyz))
            χ = @. (1/2 * N * (cos(ψ)^2 * β_aac + sin(ψ)^2 * β_bbc + β_ccc) * cos(θ) +
                    1/2 * N * (sin(ψ)^2 * β_aac + cos(ψ)^2 * β_bbc - β_ccc) * cos(θ)^3)
        elseif any(geo .== (:xzx, :zxx, :yzy, :zyy))
            χ = @. -1/2 * N * (sin(ψ)^2 * β_aac + cos(ψ)^2 * β_bbc - β_ccc) * (cos(θ) - cos(θ)^3)
        elseif geo == :zzz
            χ = @.       (N * (sin(ψ)^2 * β_aac + cos(ψ)^2 * β_bbc) * cos(θ) -
                          N * (sin(ψ)^2 * β_aac + cos(ψ)^2 * β_bbc - β_ccc) * cos(θ)^3)
        else
            error("$geo is not a valid tensor element for :$symmetry with :$mode mode.")
        end
    end

    # C2V asymmetric stretch
    if symmetry == :c2v && mode == :as
        if any(geo .== (:xxz, :yyz))
            χ = @.      -N * β_aca * sin(ψ)^2 * (cos(θ) - cos(θ)^3)
        elseif any(geo .== (:xzx, :zxx, :yzy, :zyy))
            # Different Formulae @Wang and @Tyrode
            χ = @. (1/2 * N * β_aca * (cos(ψ)^2 - sin(ψ)^2) * cos(θ) +
                          N * β_aca * sin(ψ)^2 * cos(θ)^3)
        elseif geo == :zzz
            χ = @.   2 * N * β_aca * sin(ψ)^2 * (cos(θ) - cos(θ)^3)
        else
            error("$geo is not a valid tensor element for $symmetry with $mode mode.")
        end
    end

    χ
end

function effective_susceptibility(s::Setup, pol::Symbol, geo=:none)

    β_sf = s.β[:sf]
    β_vis = s.β[:vis]
    β_ir = s.β[:ir]

    if pol == :ssp
        χ_eff = l_yy(s, :sf)  .* l_yy(s, :vis) .* l_zz(s, :ir) .* sin(β_ir)
    elseif pol == :sps
        χ_eff = l_yy(s, :sf)  .* l_zz(s, :vis) .* l_yy(s, :ir) .* sin(β_vis)
    elseif pol == :pss
        χ_eff = l_zz(s, :sf)  .* l_yy(s, :vis) .* l_yy(s, :ir) .* sin(β_sf)
    elseif pol == :ppp
        if geo == :xxz
            χ_eff = -l_xx(s, :sf) .* l_xx(s, :vis) .* l_zz(s, :ir) .* cos(β_sf) .* cos(β_vis) .* sin(β_ir)
        elseif geo == :xzx
            χ_eff = -l_xx(s, :sf) .* l_zz(s, :vis) .* l_xx(s, :ir) .* cos(β_sf) .* sin(β_vis) .* cos(β_ir)
        elseif geo == :zxx
            χ_eff = l_zz(s, :sf) .* l_xx(s, :vis) .* l_xx(s, :ir) .* sin(β_sf) .* cos(β_vis) .* cos(β_ir)
        elseif geo == :zzz
            χ_eff = l_zz(s, :sf) .* l_zz(s, :vis) .* l_zz(s, :ir) .* sin(β_sf) .* sin(β_vis) .* sin(β_ir)
        else
            error("$geo is not a valid tensor element for ppp polarization")
        end
    else
        error("$pol is not a valid polarization combination")
    end
    χ_eff
end

"""
pol: polarization :ssp, :ppp
symmetry: :c3v, :c2v
mode: :ss (symmetric stretch), :as (antisymmetric stretch)
θavg: average tilt angle of the moiety
ψavg: average ratational angle of the moiety
τ: bond angle
ρ: Raman depolarization ratio
M_c: Mass of the center atom
M_b: Mass of the bonded atoms
ω: wavenumber of the mode
θ_dist: distribution of tilt angles - :δ, :fixed, :Normal
σ_θ: standard deviation for :Normal distribution
ψ_dist: distribution of rotational angles - :iso, :Normal
σ_ψ: standard deviation for :Normal distribution
n: number of angles to calculate
"""
function effective_susceptibility(s::Setup, pol::Symbol, symmetry::Symbol,
                                  mode::Symbol, θavg::Real, ψavg::Real, τ::Real, ρ::Real,
                                  M_c::Real, M_b::Real, ω::Real;
                                  θ_dist=:δ,
                                  σ_θ=π/8,
                                  ψ_dist=:iso,
                                  σ_ψ=π/8,
                                  n=100,
                                  kwargs...)

        β_sf = s.β[:sf]
        β_vis = s.β[:vis]
        β_ir = s.β[:ir]

        if θ_dist == :δ || θ_dist==:fixed
            θ = [θavg]
            θn = [one(eltype(θ))]
        elseif θ_dist == :Normal
            θ = range(-π, stop=π, length=n)
            d = Ref(Normal(θavg, σ_θ))
            θn = pdf.(d, θ)
        else
            error("Distribution $θ_dist not defined.")
        end

        if symmetry == :c3v || ψ_dist == :fixed
            # For C3v there is no dependence on ψ
            ψ = [ψavg]
            ψn = [one(eltype(ψ))]
        elseif ψ_dist == :iso
            ψ = collect(range(-π, stop=π, length=n))
            ψn = fill(one(eltype(ψ)), n)
        elseif ψ_dist == :Normal
            ψ = range(-π, stop=π, length=n)
            d = Ref(Normal(ψavg, σ_ψ))
            ψn = pdf.(d, ψ)
        end

        χ_eff_mat = Array{ComplexF64,2}(undef, length(θ), length(ψ))
        pdf_mat = pdfgrid(θn, ψn)

        for i = 1:length(θn), j = 1:length(ψ)
            if pol == :ssp
                χ_yyz = susceptibility(:yyz, symmetry, mode, θ[i], ψ[j], τ, ρ, M_c, M_b, ω; kwargs...)
                χ_eff_mat[i,j] = l_yy(s, :sf)  .* l_yy(s, :vis) .* l_zz(s, :ir) .* sin(β_ir)  .* χ_yyz
            elseif pol == :sps
                χ_yzy = susceptibility(:yzy, symmetry, mode, θ[i], ψ[j], τ, ρ, M_c, M_b, ω; kwargs...)
                χ_eff_mat[i,j] = l_yy(s, :sf)  .* l_zz(s, :vis) .* l_yy(s, :ir) .* sin(β_vis) .* χ_yzy
            elseif pol == :pss
                χ_zyy = susceptibility(:zyy, symmetry, mode, θ[i], ψ[j], τ, ρ, M_c, M_b, ω; kwargs...)
                χ_eff_mat[i,j] = l_zz(s, :sf)  .* l_yy(s, :vis) .* l_yy(s, :ir) .* sin(β_sf)  .* χ_zyy
            elseif pol == :ppp
                χ_xxz = susceptibility(:xxz, symmetry, mode, θ[i], ψ[j], τ, ρ, M_c, M_b, ω; kwargs...)
                χ_xzx = susceptibility(:xzx, symmetry, mode, θ[i], ψ[j], τ, ρ, M_c, M_b, ω; kwargs...)
                χ_zxx = susceptibility(:zxx, symmetry, mode, θ[i], ψ[j], τ, ρ, M_c, M_b, ω; kwargs...)
                χ_zzz = susceptibility(:zzz, symmetry, mode, θ[i], ψ[j], τ, ρ, M_c, M_b, ω; kwargs...)
                χ_eff_mat[i,j] = -l_xx(s, :sf) .* l_xx(s, :vis) .* l_zz(s, :ir) .* cos(β_sf) .* cos(β_vis) .* sin(β_ir) .* χ_xxz .-
                                  l_xx(s, :sf) .* l_zz(s, :vis) .* l_xx(s, :ir) .* cos(β_sf) .* sin(β_vis) .* cos(β_ir) .* χ_xzx .+
                                  l_zz(s, :sf) .* l_xx(s, :vis) .* l_xx(s, :ir) .* sin(β_sf) .* cos(β_vis) .* cos(β_ir) .* χ_zxx .+
                                  l_zz(s, :sf) .* l_zz(s, :vis) .* l_zz(s, :ir) .* sin(β_sf) .* sin(β_vis) .* sin(β_ir) .* χ_zzz
            end
        end

        χ_eff = mean(χ_eff_mat .* pdf_mat)

end

function pdfgrid(x, y)
    # Sqrt ????
    sqrt.(x * y')
end

function ninterface(n2)
    @. 1 / √((4n2^2 + 2) / (n2^2 * (n2^2 + 5)))
end

end # module
