var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = SFGAnalysis","category":"page"},{"location":"#SFGAnalysis","page":"Home","title":"SFGAnalysis","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SFGAnalysis]","category":"page"},{"location":"#SFGAnalysis.Setup-Tuple{Tuple{Real,Real},Tuple{Real,Real},Tuple{Number,Number,Number},Tuple{Number,Number,Number}}","page":"Home","title":"SFGAnalysis.Setup","text":"Setup(ω, β, n1, n2)\n\nCreate an object representing the setup.\n\nArguments\n\nω: Tuple{Real,Real} containing omega_1 and omega_2 so that omega = omega_1 + omega_2.\nβ: Tuple{Real,Real} containing the incident angles corresponing to omega_1 and omega_2.\nn1: Tuple{Number,Number,Number} containing the refractive indices of medium 1 corresponding to omega, omega_1 and omega_2.\nn2: Tuple{Number,Number,Number} containing the refractive indices of medium 2 corresponding to omega, omega_1 and omega_2.\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.angleofrefraction-Tuple{Any,Any,Any}","page":"Home","title":"SFGAnalysis.angleofrefraction","text":"angleofrefraction(n1, n2, β)\n\nCompute the angle of refraction.\n\nArguments\n\nn1: refractive index of the medium in which the wave is propagating towards the interface\nn2: refractive index of the other medium\nβ:  incident or reflected angle\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.effective_susceptibility-Tuple{Distributions.Distribution,Setup,Symbol}","page":"Home","title":"SFGAnalysis.effective_susceptibility","text":"effective_susceptibility(d::Distribution, s::Setup, p::Symbol; <kwargs>)\n\nSame as above. In this case the tilt angle θ is expressed by a Distribution d.\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.effective_susceptibility-Tuple{Real,Setup,Symbol}","page":"Home","title":"SFGAnalysis.effective_susceptibility","text":"effective_susceptibility(θ::Real, s::Setup, p::Symbol; <kwargs>)\n\nCompute the effective susceptibility for a moiety with tilt angle θ for a Setup s in polarization p.\n\nArguments\n\nθ: tilt angle of the moiety\ns: Setup object\np: polarization combination (:ppp or :ssp)\n\nKeyword Arguments\n\npointgroup::Symbol: point group of the moiety (possible values: :c3v)\nmode::Symbol: vibrational mode (possible values: :ss (symmetric stretch), :as (asymmetric stretch))\nR::Real: the hyperpolarizability ratio R = beta_aac  beta_ccc = beta_bbc  beta_ccc. Has only to be provided for point group :c3v and mode :ss.\nbeta::Real: hyperpolarizability\nbeta_ccc for :c3v and :ss\nbeta_aca = beta_bcb for :c3v and :as\nN::Real: number of oscillators\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.effective_susceptibility_ppp-Tuple{Array{#s14,N} where N where #s14<:Number,Array{#s13,N} where N where #s13<:Number,Array{#s12,N} where N where #s12<:Number,Array{#s21,N} where N where #s21<:Number,Array{#s22,N} where N where #s22<:Number}","page":"Home","title":"SFGAnalysis.effective_susceptibility_ppp","text":"effective_susceptibility_ppp(n1, n2, n′, β, χ)\n\nCompute the effective susceptibility in ppp polarization.\n\nArguments\n\nn1: refractive indices of the medium in which the waves are propagating towards the interface\nn2: refractive indices of the other medium\nn′: refractive indices of the interface\nβ:  incident or reflected angles\nχ: susceptibilities as an array [χ_xxz, χ_xzx, χ_zxx, χ_zzz]\n\nThe arguments n1, n2, n′ and β have to be arrays where the first element corresponds to the property for the sum frequency, the second to omega_1 and the third to omega_2.\n\nhttps://doi.org/10.1080/01442350500225894\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.effective_susceptibility_ppp-Tuple{Setup,Array{#s21,N} where N where #s21<:Number}","page":"Home","title":"SFGAnalysis.effective_susceptibility_ppp","text":"effective_susceptibility_ppp(s::SFGAnalysis.Setup, χ)\n\nArguments\n\ns: setup object\nχ: susceptibilities as an array [χ_xxz, χ_xzx, χ_zxx, χ_zzz]\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.effective_susceptibility_pss-Tuple{Array{#s16,N} where N where #s16<:Number,Array{#s15,N} where N where #s15<:Number,Array{#s14,N} where N where #s14<:Number,Array{#s13,N} where N where #s13<:Number,Number}","page":"Home","title":"SFGAnalysis.effective_susceptibility_pss","text":"effective_susceptibility_pss(n1, n2, n′, β, χ_yyz)\n\nCompute the effective susceptibility in pss polarization.\n\nArguments\n\nn1: refractive indices of the medium in which the waves are propagating towards the interface\nn2: refractive indices of the other medium\nn′: refractive index of the interface for omega_2\nβ:  incident or reflected angles\nχ_yyz: susceptibility\n\nThe arguments n1, n2, and β have to be arrays where the first element corresponds to the property for the sum frequency, the second to omega_1 and the third to omega_2.\n\nhttps://doi.org/10.1080/01442350500225894\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.effective_susceptibility_pss-Tuple{Setup,Number}","page":"Home","title":"SFGAnalysis.effective_susceptibility_pss","text":"effective_susceptibility_pss(s::SFGAnalysis.Setup, χ_zyy)\n\nCompute the effective susceptibility in pss polarization.\n\nArguments\n\ns: setup object\nχ_zyy: susceptibility\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.effective_susceptibility_sps-Tuple{Array{#s16,N} where N where #s16<:Number,Array{#s15,N} where N where #s15<:Number,Array{#s14,N} where N where #s14<:Number,Array{#s13,N} where N where #s13<:Number,Number}","page":"Home","title":"SFGAnalysis.effective_susceptibility_sps","text":"effective_susceptibility_sps(n1, n2, n′, β, χ_yyz)\n\nCompute the effective susceptibility in sps polarization.\n\nArguments\n\nn1: refractive indices of the medium in which the waves are propagating towards the interface\nn2: refractive indices of the other medium\nn′: refractive index of the interface for omega_2\nβ:  incident or reflected angles\nχ_yyz: susceptibility\n\nThe arguments n1, n2, and β have to be arrays where the first element corresponds to the property for the sum frequency, the second to omega_1 and the third to omega_2.\n\nhttps://doi.org/10.1080/01442350500225894\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.effective_susceptibility_sps-Tuple{Setup,Number}","page":"Home","title":"SFGAnalysis.effective_susceptibility_sps","text":"effective_susceptibility_sps(s::SFGAnalysis.Setup, χ_yzy)\n\nCompute the effective susceptibility in sps polarization.\n\nArguments\n\ns: setup object\nχ_yzy: susceptibility\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.effective_susceptibility_ssp-Tuple{Array{#s16,N} where N where #s16<:Number,Array{#s15,N} where N where #s15<:Number,Array{#s14,N} where N where #s14<:Number,Array{#s13,N} where N where #s13<:Number,Number}","page":"Home","title":"SFGAnalysis.effective_susceptibility_ssp","text":"effective_susceptibility_ssp(n1, n2, n′, β, χ_yyz)\n\nCompute the effective susceptibility in ssp polarization.\n\nArguments\n\nn1: refractive indices of the medium in which the waves are propagating towards the interface\nn2: refractive indices of the other medium\nn′: refractive index of the interface for omega_2\nβ:  incident or reflected angles\nχ_yyz: susceptibility\n\nThe arguments n1, n2, and β have to be arrays where the first element corresponds to the property for the sum frequency, the second to omega_1 and the third to omega_2.\n\nhttps://doi.org/10.1080/01442350500225894\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.effective_susceptibility_ssp-Tuple{Setup,Number}","page":"Home","title":"SFGAnalysis.effective_susceptibility_ssp","text":"effective_susceptibility_ssp(s::SFGAnalysis.Setup, χ_yyz)\n\nCompute the effective susceptibility in ssp polarization.\n\nArguments\n\ns: setup object\nχ_yyz: susceptibility\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.freq2wl-Tuple{Any}","page":"Home","title":"SFGAnalysis.freq2wl","text":"freq2wl(ω)\n\nCompute the wavelength of a field with frequency ω in 1/s. The wavelength will be given in nm.\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.fresnel_x-Tuple{Any,Any,Any}","page":"Home","title":"SFGAnalysis.fresnel_x","text":"fresnel_x(n1, n2, β)\n\nCompute the diagonal elements L_xx of the fresnel factor.\n\nArguments\n\nn1: refractive index of the medium in which the wave is propagating towards the interface\nn2: refractive index of the other medium\nβ:  incident or reflected angle\n\nhttps://doi.org/10.1103/PhysRevB.59.1263\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.fresnel_y-Tuple{Any,Any,Any}","page":"Home","title":"SFGAnalysis.fresnel_y","text":"fresnel_y(n1, n2, β)\n\nCompute the diagonal elements L_yy of the fresnel factor.\n\nArguments\n\nn1: refractive index of the medium in which the wave is propagating towards the interface\nn2: refractive index of the other medium\nβ:  incident or reflected angle\n\nhttps://doi.org/10.1103/PhysRevB.59.1263\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.fresnel_z-NTuple{4,Any}","page":"Home","title":"SFGAnalysis.fresnel_z","text":"fresnel_z(n1, n2, n′, β)\n\nCompute the diagonal elements L_yy of the fresnel factor.\n\nArguments\n\nn1: refractive index of the medium in which the wave is propagating towards the interface\nn2: refractive index of the other medium\nn′: refractive index of the interface\nβ:  incident or reflected angle\n\nhttps://doi.org/10.1103/PhysRevB.59.1263\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.gaussian-Tuple{Any,Any,Any}","page":"Home","title":"SFGAnalysis.gaussian","text":"\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.lorentzian-NTuple{4,Any}","page":"Home","title":"SFGAnalysis.lorentzian","text":"Documentation\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.ninterface-Tuple{Any}","page":"Home","title":"SFGAnalysis.ninterface","text":"ninterface(n2)\n\nCalculate the index of refraction of the interfacial layer according to\n\nleft(frac1n^primeright)^2 = frac4n_2^2 + 2n_2^2(n_2 ^2+5).\n\nArguments\n\nn2: refractive index of the medium in which the wave is not propagating towards the interface\n\nhttps://doi.org/10.1103/PhysRevB.59.1263\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.sfangle-NTuple{5,Any}","page":"Home","title":"SFGAnalysis.sfangle","text":"sfangle(ω, ω1, ω2, β1, β2)\n\nCompute the angle of the reflected sum frequency field. In a counter-propagating geometry one of the βs gets negative.\n\nArguments\n\nω:  sum frequency\nω1: frequency of wave 1\nω2: frequency of wave 2\nβ1:  incident angle for field with frequency ω1\nβ2:  incident angle for field with frequency ω2\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.sfgintensity-Tuple{Any,Any,Any}","page":"Home","title":"SFGAnalysis.sfgintensity","text":"sfgintensity(ω, β, χ2_eff; <keyword arguments>)\n\nCompute the intentensity of the sum frequency in the reflected direction according to\n\nI_mathrmSF = frac8pi^3omega^2sec^2betac^3 n n_1 n_2 left chi_mathrmeff^(2) right^2 I_1 I_2.\n\nArguments\n\nω:  sum frequency\nβ:  reflection angle of the sum frequency field\nχ2_eff: effective susceptibility\n\nKeyword Arguments\n\nn:  refractive index of the bulk medium for ω\nn1: refractive index of the bulk medium for ω1\nn2: refractive index of the bulk medium for ω2\nI1: intensity of ω1\nI2: intensity of ω2\n\nhttps://doi.org/10.1103/PhysRevB.59.12632\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.sfgintensity-Tuple{Any,Any}","page":"Home","title":"SFGAnalysis.sfgintensity","text":"sfgintensity(s::Setup, χ2_eff::Number; <keyword arguments>)\n\nSame as sfgintensity(ω, β, χ2_eff; <keyword arguments>).\n\nArguments\n\ns:  Setup object\nχ2_eff: effective susceptibility\n\nKeyword Arguments\n\nI1: intensity of ω1\nI2: intensity of ω2\n\nhttps://doi.org/10.1103/PhysRevB.59.12632\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.sfspectrum-NTuple{4,Any}","page":"Home","title":"SFGAnalysis.sfspectrum","text":"If no phase information and non-resonant background is passed to the function, we set it to 0.\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.sfspectrum-NTuple{5,Any}","page":"Home","title":"SFGAnalysis.sfspectrum","text":"If no phase information and but the non-resonant background is passed it's assumed that the relative phase of the resonances is +/- π, which can be controlled via the sign of A.\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.sfspectrum-Tuple{Any,Any,Any,Any,Number,Number}","page":"Home","title":"SFGAnalysis.sfspectrum","text":"If just one phase is passed, we assume that this belongs to the non-resonant background.\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.sfspectrum-Union{Tuple{T}, Tuple{Number,T,T,T,T,Number}} where T<:AbstractArray","page":"Home","title":"SFGAnalysis.sfspectrum","text":"sfspectrum(x::Number, A::T, ω::T, Γ::T[, φ::AbstractArray, χnr]) where {T<:AbstractArray}\n\nCalculate the sum frequency spectrum according to\n\nI = left chi_NR + sum_q fracA_qx - ω_q - iGamma_q right^2\n\nArguments\n\nx the frequency/wavenumber at which to evaluate the sfspectrum\nA oscillator strength\nω resonances\nΓ damping coefficients\n\nKeyword Arguments\n\nExamples\n\nx = range(2800, 3000, length=201)\nA = [1, 1]\nω = [2880, 2930]\nΓ = [8, 7]\n\nIn order to use the broadcast functionality we have to wrap the parameters in Refs:\n\ny = sfspectrum.(x, Ref(A), Ref(ω), Ref(Γ))\n\nAlternatively use array comprehensions:\n\ny = [sfspectrum(_x, A, ω, Γ) for _x in x]\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.sfspectrum_voigt-NTuple{5,Any}","page":"Home","title":"SFGAnalysis.sfspectrum_voigt","text":"\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.sfspectrum_voigt-Union{Tuple{T}, Tuple{Number,T,T,T,T,T,Number}} where T<:AbstractArray","page":"Home","title":"SFGAnalysis.sfspectrum_voigt","text":"\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.susceptibility-Tuple{Real,Symbol}","page":"Home","title":"SFGAnalysis.susceptibility","text":"susceptibility(θ::Real, t::Symbol; pointgroup, mode, R, beta, N; <kwargs>)\n\nCompute the susceptibility for angle θ and tensor element t. Tensor elements important for specific polarization combinations:\n\nssp: :yyz\nsps: :yzy\npss: :zyy\nppp: :xxz, :xzx, :zxx and :zzz\n\nArguments\n\nθ: title angle of the moiety\nt: tensor element (see above)\n\nKeyword Arguments\n\npointgroup::Symbol: point group of the moiety (possible values: :c3v)\nmode::Symbol: vibrational mode (possible values: :ss (symmetric stretch), :as (asymmetric stretch))\nR::Real: the hyperpolarizability ratio R = beta_aac  beta_ccc = beta_bbc  beta_ccc. Has only to be provided for point group :c3v and mode :ss.\nbeta::Real: hyperpolarizability\nbeta_ccc for :c3v and :ss\nbeta_aca = beta_bcb for :c3v and :as\nN::Real: number of oscillators\n\nMolecule Group R\nmethanol methyl 1.7\nethanol and longer X-1-ol methyl 3.4\nacetone methyl 1.9\nDMSO methyl 2.3\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.susceptibility_c3v_as_xxz-Tuple{Any}","page":"Home","title":"SFGAnalysis.susceptibility_c3v_as_xxz","text":"Susceptibility for tensor elements xxz and yyz.\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.susceptibility_c3v_as_xzx-Tuple{Any}","page":"Home","title":"SFGAnalysis.susceptibility_c3v_as_xzx","text":"Susceptibility for tensor elements xzx, zxx, yzy and zyy.\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.susceptibility_c3v_as_zzz-Tuple{Any}","page":"Home","title":"SFGAnalysis.susceptibility_c3v_as_zzz","text":"Susceptibility for tensor element zzz\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.susceptibility_c3v_ss_xxz-Tuple{Any,Any}","page":"Home","title":"SFGAnalysis.susceptibility_c3v_ss_xxz","text":"Susceptibility for tensor elements xxz and yyz.\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.susceptibility_c3v_ss_xzx-Tuple{Any,Any}","page":"Home","title":"SFGAnalysis.susceptibility_c3v_ss_xzx","text":"Susceptibility for tensor elements xzx, zxx, yzy and zyy.\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.susceptibility_c3v_ss_zzz-Tuple{Any,Any}","page":"Home","title":"SFGAnalysis.susceptibility_c3v_ss_zzz","text":"Susceptibility for tensor element zzz\n\n\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.voigt-NTuple{5,Any}","page":"Home","title":"SFGAnalysis.voigt","text":"\n\n\n\n","category":"method"},{"location":"#SFGAnalysis.wl2freq-Tuple{Any}","page":"Home","title":"SFGAnalysis.wl2freq","text":"wl2freq(λ)\n\nCompute the frequency of a field with wavelength λ in nm.\n\n\n\n\n\n","category":"method"}]
}
