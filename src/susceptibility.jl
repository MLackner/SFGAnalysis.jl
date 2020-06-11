"""
Susceptibility for tensor elements xxz and yyz.
"""
function susceptibility_c3v_ss_xxz(θ, R; β_ccc=1.0, N=1.0)
    1/2 * N * β_ccc * ( (1+R) * cos(θ) - (1-R) * cos(θ)^3 )
end

"""
Susceptibility for tensor elements xzx, zxx, yzy and zyy.
"""
function susceptibility_c3v_ss_xzx(θ, R; β_ccc=1.0, N=1.0)
    1/2 * N * β_ccc * (1-R) * ( cos(θ) - cos(θ)^3 )
end

"""
Susceptibility for tensor element zzz
"""
function susceptibility_c3v_ss_zzz(θ, R; β_ccc=1.0, N=1.0)
    N * β_ccc * ( R * cos(θ) + (1-R) * cos(θ)^3 )
end

"""
Susceptibility for tensor elements xxz and yyz.
"""
function susceptibility_c3v_as_xxz(θ; β_aca=1.0, N=1.0)
    -N * β_aca * (cos(θ) - cos(θ)^3)
end

"""
Susceptibility for tensor elements xzx, zxx, yzy and zyy.
"""
function susceptibility_c3v_as_xzx(θ; β_aca=1.0, N=1.0)
    N * β_aca * cos(θ)^3
end

"""
Susceptibility for tensor element zzz
"""
function susceptibility_c3v_as_zzz(θ; β_aca=1.0, N=1.0)
    2N * β_aca * (cos(θ) - cos(θ)^3)
end

"""
    susceptibility(θ::Real, t::Symbol; pointgroup, mode, R, beta, N)

Compute the susceptibility for angle `θ` and tensor element `t`.
Tensor elements important for specific polarization combinations:
* ssp: `:yyz`
* sps: `:yzy`
* pss: `:zyy`
* ppp: `:xxz`, `:xzx`, `:zxx` and `:zzz`

## Keyword Arguments
* `pointgroup::Symbol`: point group of the moiety (possible values: `:c3v`)
* `mode::Symbol`: vibrational mode (possible values: `:ss` (symmetric stretch), `:as` (asymmetric stretch))
* `R::Real`: the hyperpolarizability ratio \$R = \\beta_{aac} / \\beta_{ccc} = \\beta_{bbc} / \\beta_{ccc}\$. Has only to be provided for point group `:c3v` and mode `:ss`.
* `beta`: hyperpolarizability
    * \$\\beta_{ccc}\$ for `:c3v` and `:ss`
    * \$\\beta_{aca} = \\beta_{bcb}\$ for `:c3v` and `:as`

| Molecule                  	| Group  	| R   	|
|:--------------------------	|:-------	|:----	|
| methanol                  	| methyl 	| 1.7 	|
| ethanol and longer X-1-ol 	| methyl 	| 3.4 	|
| acetone                   	| methyl 	| 1.9 	|
| DMSO                      	| methyl 	| 2.3 	|
"""
function susceptibility(θ::Real, t::Symbol;
        pointgroup::Union{Nothing,Symbol}=nothing, mode::Union{Nothing,Symbol}=nothing, 
        R::Union{Nothing,Real}=nothing, beta::Real=1.0, N::Real=1.0)
    # Run checks
    isnothing(pointgroup) && error("no pointgroup provided")
    isnothing(mode) && error("no mode provided")
    pointgroup in [:c3v] || error("unknown point group $pointgroup")
    mode in [:ss, :as] || error("unknown mode $mode")
    t in [:xxz, :yyz, :xzx, :zxx, :yzy, :zyy, :zzz] || error("unknown tensor element $t")
    # For the point group of c3v we have to provide an R value
    pointgroup in [:c3v] && mode == :ss && isnothing(R) && error("You have to provide the keyword argument R.")

    if pointgroup == :c3v
        if mode == :ss
            if t in [:xxz, :yyz]
                χ = susceptibility_c3v_ss_xxz(θ, R; β_ccc=beta, N=N)
            elseif t in [:xzx, :zxx, :yzy, :zyy]
                χ = susceptibility_c3v_ss_xzx(θ, R; β_ccc=beta, N=N)
            elseif t in [:zzz]
                χ = susceptibility_c3v_ss_zzz(θ, R; β_ccc=beta, N=N)
            end
        elseif mode == :as
            if t in [:xxz, :yyz]
                χ = susceptibility_c3v_as_xxz(θ; β_aca=beta, N=N)
            elseif t in [:xzx, :zxx, :yzy, :zyy]
                χ = susceptibility_c3v_as_xzx(θ; β_aca=beta, N=N)
            elseif t in [:zzz]
                χ = susceptibility_c3v_as_zzz(θ; β_aca=beta, N=N)
            end
        end
    end
    return χ
end

# """
#     susceptibility(θ::AbstractArray{<:Real,1}, t::Array{Symbol}; pointgroup, mode, R, beta, N)

# Compute the susceptibility for angles `θ` and tensor elements `t`. Returns an array with
# length `length(θ)` where each entry represents the susceptibility for all tensor elements
# at the given angle.
# """
# function susceptibility(θ::Real, t::Array{Symbol,1};
#     pointgroup::Union{Nothing,Symbol}=nothing, mode::Union{Nothing,Symbol}=nothing, 
#     R::Union{Nothing,Real}=nothing, beta::Real=1.0, N::Real=1.0)
#     [susceptibility(θ, _t; pointgroup=pointgroup, mode=mode, R=R) for _t in t]
# end

"""
    effective_susceptibility_ssp(n1, n2, n′, β, χ_yyz)

Compute the effective susceptibility in ppp polarization.

# Arguments
* `n1`: refractive indices of the medium in which the waves are propagating towards the interface
* `n2`: refractive indices of the other medium
* `n′`: refractive index of the interface for \$\\omega_2\$
* `β`:  incident or reflected angles
* `χ_yyz`: susceptibility

The arguments `n1`, `n2`, and `β` have to be arrays where the first element
corresponds to the property for the sum frequency, the second to \$\\omega_1\$
and the third to \$\\omega_2\$.

https://doi.org/10.1080/01442350500225894
"""
function effective_susceptibility_ssp(
        n1::Array{<:Number},
        n2::Array{<:Number}, 
        n′::Array{<:Number}, 
        β::Array{<:Number}, 
        χ_yyz::Number
    )
    χ2_eff = fresnel_y(n1[1], n2[1], β[1]) * 
             fresnel_y(n1[2], n2[2], β[2]) * 
             fresnel_z(n1[3], n2[3], n′[3], β[3]) *
             sin(β[3]) * χ_yyz
end

"""
    effective_susceptibility_ssp(s::SFGAnalysis.Setup, χ_yyz)

# Arguments
* `s`: setup object
* `χ_yyz`: susceptibility
"""
function effective_susceptibility_ssp(s::Setup, χ_yyz::Number)
    n1 = [s.Ω.n1, s.Ω1.n1, s.Ω2.n1]
    n2 = [s.Ω.n2, s.Ω1.n2, s.Ω2.n2]
    n′ = [s.Ω.n′, s.Ω1.n′, s.Ω2.n′]
    β = [s.Ω.β, s.Ω1.β, s.Ω2.β]
    effective_susceptibility_ssp(n1, n2, n′, β, χ_yyz)
end

"""
    effective_susceptibility_ppp(n1, n2, n′, β, χ)

Compute the effective susceptibility in ppp polarization.

# Arguments
* `n1`: refractive indices of the medium in which the waves are propagating towards the interface
* `n2`: refractive indices of the other medium
* `n′`: refractive indices of the interface
* `β`:  incident or reflected angles
* `χ`: susceptibilities as an array `[χ_xxz, χ_xzx, χ_zxx, χ_zzz]`

The arguments `n1`, `n2`, `n′` and `β` have to be arrays where the first element
corresponds to the property for the sum frequency, the second to \$\\omega_1\$
and the third to \$\\omega_2\$.

https://doi.org/10.1080/01442350500225894
"""
function effective_susceptibility_ppp(n1::Array{<:Number}, 
                n2::Array{<:Number}, n′::Array{<:Number}, 
                β::Array{<:Number}, χ::Array{<:Number})
    χ_xxz, χ_xzx, χ_zxx, χ_zzz = χ
    χ2_eff = -fresnel_x(n1[1], n2[1], β[1]) * 
             fresnel_x(n1[2], n2[2], β[2]) * 
             fresnel_z(n1[3], n2[3], n′[3], β[3]) *
             cos(β[1]) * cos(β[2]) * sin(β[3]) * χ_xxz -
             fresnel_x(n1[1], n2[1], β[1]) * 
             fresnel_z(n1[2], n2[2], n′[2], β[2]) *
             fresnel_x(n1[3], n2[3], β[3]) * 
             cos(β[1]) * sin(β[2]) * cos(β[3]) * χ_xzx +
             fresnel_z(n1[1], n2[1], n′[1], β[1]) *
             fresnel_x(n1[2], n2[2], β[2]) * 
             fresnel_x(n1[3], n2[3], β[3]) * 
             sin(β[1]) * cos(β[2]) * cos(β[3]) * χ_zxx +
             fresnel_z(n1[1], n2[1], n′[1], β[1]) *
             fresnel_z(n1[2], n2[2], n′[2], β[2]) *
             fresnel_z(n1[3], n2[3], n′[3], β[3]) *
             sin(β[1]) * sin(β[2]) * sin(β[3]) * χ_zzz
end

"""
    effective_susceptibility_ppp(s::SFGAnalysis.Setup, χ)

# Arguments
* `s`: setup object
* `χ`: susceptibilities as an array `[χ_xxz, χ_xzx, χ_zxx, χ_zzz]`
"""
function effective_susceptibility_ppp(s::Setup, χ::Array{<:Number})
    n1 = [s.Ω.n1, s.Ω1.n1, s.Ω2.n1]
    n2 = [s.Ω.n2, s.Ω1.n2, s.Ω2.n2]
    n′ = [s.Ω.n′, s.Ω1.n′, s.Ω2.n′]
    β = [s.Ω.β, s.Ω1.β, s.Ω2.β]
    effective_susceptibility_ppp(n1, n2, n′, β, χ)
end

"""
"""
function effective_susceptibility(θ::Real, s::Setup, p::Symbol;
    pointgroup::Union{Nothing,Symbol}=nothing, 
    mode::Union{Nothing,Symbol}=nothing, 
    R::Union{Nothing,Real}=nothing, 
    beta::Real=1.0, N::Real=1.0)

    p in [:ssp, :ppp] || error("unknown polarization combination $p")

    if p == :ppp
        χ = susceptibility.(θ, [:xxz, :xzx, :zxx, :zzz]; 
            pointgroup=pointgroup, 
            mode=mode, R=R,
            beta=beta, N=N
        )
        χ_eff = effective_susceptibility_ppp(s, χ)
    elseif p == :ssp
        χ = susceptibility(θ, :yyz; 
            pointgroup=pointgroup,
            mode=mode, R=R,
            beta=beta, N=N
        )
        χ_eff = effective_susceptibility_ssp(s, χ)
    end
    χ_eff
end