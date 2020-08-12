"""
    sfspectrum(x::Number, A::T, ω::T, Γ::T[, φ::AbstractArray, χnr]) where {T<:AbstractArray}

Calculate the sum frequency spectrum according to

``I = \\left| \\chi_{NR} + \\sum_q \\frac{A_q}{x - ω_q - i\\Gamma_q} \\right|^2``

# Arguments
* `x` the frequency/wavenumber at which to evaluate the sfspectrum
* `A` oscillator strength
* `ω` resonances
* `Γ` damping coefficients

# Keyword Arguments

# Examples
```
x = range(2800, 3000, length=201)
A = [1, 1]
ω = [2880, 2930]
Γ = [8, 7]
```
In order to use the broadcast functionality we have to wrap the parameters in Refs:
```
y = sfspectrum.(x, Ref(A), Ref(ω), Ref(Γ))
```
Alternatively use array comprehensions:
```
y = [sfspectrum(_x, A, ω, Γ) for _x in x]
```
"""
function sfspectrum(x::Number, A::T, ω::T, Γ::T, φ::T, χnr::Number) where {T<:AbstractArray}
    y = χnr
    for i in eachindex(A)
        y += lorentzian(x, A[i], ω[i], Γ[i]) * exp(1im * φ[i])
    end
    abs2(y)
end

"""
If no phase information and non-resonant background is passed to the function, we set it to 0.
"""
function sfspectrum(x, A, ω, Γ)
    T = eltype(A)
    N = length(A)
    sfspectrum(x, A, ω, Γ, zeros(T,N), zero(T))
end

"""
If just one phase is passed, we assume that this belongs to the non-resonant background.
"""
function sfspectrum(x, A, ω, Γ, φ::Number, χnr::Number)
    T = eltype(A)
    N = length(A)
    φ_array = zeros(T,N)
    χnr_complex = χnr * exp(1im * φ)
    sfspectrum(x, A, ω, Γ, φ_array, χnr_complex)
end

"""
If no phase information and but the non-resonant background is passed it's assumed that the relative phase of the resonances is +/- π, which can be controlled via the sign of `A`.
"""
function sfspectrum(x, A, ω, Γ, χnr)
    T = eltype(A)
    N = length(A)
    φ = zeros(T,N)
    sfspectrum(x, A, ω, Γ, φ, χnr)
end

"""
Documentation
"""
function lorentzian(x, A, ω, Γ)
    A / (x - ω + 1im * Γ)
end