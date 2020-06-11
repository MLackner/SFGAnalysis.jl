
"""
    sfgintensity(ω, ω1, ω2, β, χ2_eff; <keyword arguments>)

Compute the intentensity of the sum frequency of `ω1` and `ω2`
in the reflected direction according to

``I_\\mathrm{SF} = \\frac{8\\pi^3\\omega^2\\sec^2\\beta}{c^3 n n_1 n_2}
\\left| \\chi_\\mathrm{eff}^{(2)} \\right|^2 I_1 I_2``.

# Arguments
* `ω`:  sum frequency
* `ω1`: frequency of wave 1
* `ω2`: frequency of wave 2
* `β`:  reflection angle of the sum frequency field

# Keyword Arguments
* `n`:  refractive index of the bulk medium for `ω`
* `n1`: refractive index of the bulk medium for `ω1`
* `n2`: refractive index of the bulk medium for `ω2`
* `I1`: intensity of `ω1`
* `I2`: intensity of `ω2`

https://doi.org/10.1103/PhysRevB.59.12632
"""
function sfgintensity(
        ω, ω1, ω2, β, χ2_eff;
        n=1.0, n1=1.0, n2=1.0, I1=1.0, I2=1.0
    )
    c = 299792458.0 # speed of light in vacuum
    8π^3 * ω^2 * sec(β)^2 / (c^3 * n * n1 * n2) * abs2(χ2_eff) * I1 * I2
end


"""
    sfangle(ω, ω1, ω2, β1, β2)

Compute the angle of the reflected sum frequency field.
In a counter-propagating geometry one of the `β`s gets negative.

# Arguments
* `ω`:  sum frequency
* `ω1`: frequency of wave 1
* `ω2`: frequency of wave 2
* `β1`:  incident angle for field with frequency `ω1`
* `β2`:  incident angle for field with frequency `ω2`
"""
function sfangle(ω, ω1, ω2, β1, β2)
    asin(ω1 / ω * sin(β1) + ω2 / ω * sin(β2))
end

"""
    wl2freq(λ)

Compute the frequency of a field with wavelength `λ` in nm.
"""
wl2freq(λ) = 299792458.0 / (λ * 1e-9)

"""
    freq2wl(ω)

Compute the wavelength of a field with frequency `ω` in 1/s.
The wavelength will be given in nm.
"""
freq2wl(ω) = 299792458.0 / ω * 1e9