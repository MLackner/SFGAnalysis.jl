"""
    angleofrefraction(n1, n2, β)

Compute the angle of refraction.

# Arguments
* `n1`: refractive index of the medium in which the wave is propagating towards the interface
* `n2`: refractive index of the other medium
* `β`:  incident or reflected angle
"""
function angleofrefraction(n1, n2, β)
    asin(n1 / n2 * sin(β))
end

"""
    fresnel_x(n1, n2, β)

Compute the diagonal elements \$L_{xx}\$ of the fresnel factor.

# Arguments
* `n1`: refractive index of the medium in which the wave is propagating towards the interface
* `n2`: refractive index of the other medium
* `β`:  incident or reflected angle

https://doi.org/10.1103/PhysRevB.59.1263
"""
function fresnel_x(n1, n2, β)
    γ = angleofrefraction(n1, n2, β)
    2n1 * cos(γ) / (n1 * cos(γ) + n2 * cos(β))
end

"""
    fresnel_y(n1, n2, β)

Compute the diagonal elements \$L_{yy}\$ of the fresnel factor.

# Arguments
* `n1`: refractive index of the medium in which the wave is propagating towards the interface
* `n2`: refractive index of the other medium
* `β`:  incident or reflected angle

https://doi.org/10.1103/PhysRevB.59.1263
"""
function fresnel_y(n1, n2, β)
    γ = angleofrefraction(n1, n2, β)
    2n1 * cos(β) / (n1 * cos(β) + n2 * cos(γ))
end

"""
    fresnel_z(n1, n2, n′, β)

Compute the diagonal elements \$L_{yy}\$ of the fresnel factor.

# Arguments
* `n1`: refractive index of the medium in which the wave is propagating towards the interface
* `n2`: refractive index of the other medium
* `n′`: refractive index of the interface
* `β`:  incident or reflected angle

https://doi.org/10.1103/PhysRevB.59.1263
"""
function fresnel_z(n1, n2, n′, β)
    γ = angleofrefraction(n1, n2, β)
    2n2 * cos(β) / (n1 * cos(γ) + n2 * cos(β)) * (n1 / n′)^2
end

"""
    ninterface(n2)

Calculate the index of refraction of the interfacial layer according to

``\\left(\\frac{1}{n^\\prime}\\right)^2 = \\frac{4n_2^2 + 2}{n_2^2(n_2 ^2+5)}``.

# Arguments
* `n2`: refractive index of the medium in which the wave is **not** propagating towards the interface

https://doi.org/10.1103/PhysRevB.59.1263
"""
function ninterface(n2)
    1 / sqrt((4 * n2^2 + 2) / (n2^2 * (n2^2 + 5)))
end
