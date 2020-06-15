using SFGAnalysis
using Test
import Distributions: Normal

@testset "SFGAnalysis.jl" begin
    λ1 = 515
    λ2 = 3450
    ω1 = wl2freq(λ1)
    ω2 = wl2freq(λ2)
    ω = ω1 + ω2
    λ = freq2wl(ω)
    β1 = π/4
    β2 = π/6
    β = sfangle(ω, ω1, ω2, β1, β2)

    #TODO: are these really frequencies or angular frequencies or wavenumbers?
    @test wl2freq(λ1) ≈ 5.821212776699029e14 # source?
    @test λ ≈ 448.1084489281211  # source? but looks correct

    # Test Setup Type
    n1  = 1.0
    n11 = 1.0
    n12 = 1.0
    n2  = 1.3
    n21 = 1.3
    n22 = 1.4
    s = Setup((ω1, ω2), (β1, β2), (n1, n11, n12), (n2, n21, n22))

    # test sfgintensity
    χ2_eff = 1e-5 + 1e-6im   # completely arbitrary. change for calculated value
    @test sfgintensity(ω, β, χ2_eff) == sfgintensity(s, χ2_eff)

    # Test ninterface
    # https://doi.org/10.1103/PhysRevB.59.1263 for verification
    @test round(ninterface(1.5), digits=2) == 1.22
    @test round(ninterface(1.34), digits=2) == 1.15

    # Test susceptibility
    for t in [:yyz, :yzy, :zyy, :xxz, :xzx, :zxx, :zzz]
        susceptibility(π/4, t; pointgroup=:c3v, mode=:ss, R=1.7)
    end
    for t in [:yyz, :yzy, :zyy, :xxz, :xzx, :zxx, :zzz]
        susceptibility(π/4, t; pointgroup=:c3v, mode=:as)
    end
    @test_throws ErrorException susceptibility(π/4, :yyz)
    @test_throws ErrorException susceptibility(π/4, :yyz; pointgroup=:c3v)
    @test_throws ErrorException susceptibility(π/4, :yyz; mode=:c3v)
    @test_throws ErrorException susceptibility(π/4, :yyz; pointgroup=:c3v, mode=:c3v)

    # Test effective_susceptibility
    χ_ssp = susceptibility(π/8, :yyz; pointgroup=:c3v, mode=:as)
    χ_ppp = [susceptibility(π/8, t; pointgroup=:c3v, mode=:as) for t in [:xxz, :xzx, :zxx, :zzz]]
    χ2_eff_ssp = effective_susceptibility_ssp(s, χ_ssp)
    χ2_eff_ppp = effective_susceptibility_ppp(s, χ_ppp)

end

@testset "J. Phys. Chem. C 2007, 111, 8716-8725" begin
    # Test https://doi.org/10.1021/jp067062h
    ω1 = wl2freq(532)
    ω2 = wl2freq(3333)
    β1 = [deg2rad(62), deg2rad(37)]
    β2 = [deg2rad(53), deg2rad(51)]
    n1  = 1
    n11 = 1
    n12 = 1
    n2  = 1.33
    n21 = 1.33
    n22 = 1.33
    R = 1.7
    s = Setup((ω1, ω2), (β1[1], β2[1]), (n1, n11, n12), (n2, n21, n22))

    θ = range(0, π/2, length=91)

    # a factor is appended to the effective susceptibilities
    χ_eff_ssp_ss = effective_susceptibility.(θ, s, :ssp; pointgroup=:c3v, mode=:ss, R=R) * sec(s.Ω.β)
    χ_eff_ssp_as = effective_susceptibility.(θ, s, :ssp; pointgroup=:c3v, mode=:as) * sec(s.Ω.β)
    χ_eff_ppp_ss = effective_susceptibility.(θ, s, :ppp; pointgroup=:c3v, mode=:ss, R=R) * sec(s.Ω.β)
    χ_eff_ppp_as = effective_susceptibility.(θ, s, :ppp; pointgroup=:c3v, mode=:as) * sec(s.Ω.β)

    @test round(χ_eff_ssp_ss[1], digits=1) == 0.9
    @test round(χ_eff_ssp_ss[45], digits=1) == 0.6
    @test round(χ_eff_ssp_ss[90], digits=1) == 0.0

    @test round(χ_eff_ppp_ss[1], digits=1) == -0.1
    @test round(χ_eff_ppp_ss[45], digits=1) == 0.1
    @test round(χ_eff_ppp_ss[90], digits=1) == 0.0

    @test round(χ_eff_ssp_as[1], digits=1) == 0.0
    @test round(χ_eff_ssp_as[45], digits=1) == -0.2
    @test round(χ_eff_ssp_as[90], digits=1) == 0.0

    @test round(χ_eff_ppp_as[1], digits=1) == 0.0
    @test round(χ_eff_ppp_as[45], digits=1) == 0.4
    @test round(χ_eff_ppp_as[90], digits=1) == 0.0
end

@testset "Distributions" begin
    # Test the distribution functionality
    # Do everything like in the testse "J. Phys. Chem. C 2007, 111, 8716-8725"
    # This time we just assume a very narrow distribution. This way the results
    # should be the same.
    ω1 = wl2freq(532)
    ω2 = wl2freq(3333)
    β1 = [deg2rad(62), deg2rad(37)]
    β2 = [deg2rad(53), deg2rad(51)]
    n1  = 1
    n11 = 1
    n12 = 1
    n2  = 1.33
    n21 = 1.33
    n22 = 1.33
    R = 1.7
    s = Setup((ω1, ω2), (β1[1], β2[1]), (n1, n11, n12), (n2, n21, n22))
    display(s)

    dist = [Normal(θ, π/100) for θ in range(0, π/2, length=91)]

    # a factor is appended to the effective susceptibilities
    χ_eff_ssp_ss = effective_susceptibility.(dist, s, :ssp; pointgroup=:c3v, mode=:ss, R=R) * sec(s.Ω.β)
    χ_eff_ssp_as = effective_susceptibility.(dist, s, :ssp; pointgroup=:c3v, mode=:as) * sec(s.Ω.β)
    χ_eff_ppp_ss = effective_susceptibility.(dist, s, :ppp; pointgroup=:c3v, mode=:ss, R=R) * sec(s.Ω.β)
    χ_eff_ppp_as = effective_susceptibility.(dist, s, :ppp; pointgroup=:c3v, mode=:as) * sec(s.Ω.β)

    @test round(χ_eff_ssp_ss[1], digits=1) == 0.9
    @test round(χ_eff_ssp_ss[45], digits=1) == 0.6
    @test round(χ_eff_ssp_ss[90], digits=1) == 0.0

    @test round(χ_eff_ppp_ss[1], digits=1) == -0.1
    @test round(χ_eff_ppp_ss[45], digits=1) == 0.1
    @test round(χ_eff_ppp_ss[90], digits=1) == 0.0

    @test round(χ_eff_ssp_as[1], digits=1) == 0.0
    @test round(χ_eff_ssp_as[45], digits=1) == -0.2
    @test round(χ_eff_ssp_as[90], digits=1) == 0.0

    @test round(χ_eff_ppp_as[1], digits=1) == 0.0
    @test round(χ_eff_ppp_as[45], digits=1) == 0.4
    @test round(χ_eff_ppp_as[90], digits=1) == 0.0
end
