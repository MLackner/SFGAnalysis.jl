import PrettyTables: pretty_table

mutable struct Wave
    ω::Real
    β::Real
    n1::Number
    n2::Number
    n′::Number
end

mutable struct Setup
    Ω::Wave
    Ω1::Wave
    Ω2::Wave
end

"""
    Setup(ω, β, n1, n2)

Create an object representing the setup.

# Arguments
* `ω`: `Tuple{Real,Real}` containing \$\\omega_1\$ and \$\\omega_2\$ so that \$\\omega = \\omega_1 + \\omega_2\$.
* `β`: `Tuple{Real,Real}` containing the incident angles corresponing to \$\\omega_1\$ and \$\\omega_2\$.
* `n1`: `Tuple{Number,Number,Number}` containing the refractive indices of medium 1 corresponding to \$\\omega\$, \$\\omega_1\$ and \$\\omega_2\$.
* `n2`: `Tuple{Number,Number,Number}` containing the refractive indices of medium 2 corresponding to \$\\omega\$, \$\\omega_1\$ and \$\\omega_2\$.
"""
function Setup(
        (ω1, ω2)::Tuple{Real,Real}, 
        (β1, β2)::Tuple{Real,Real},
        (n1, n11, n12)::Tuple{Number,Number,Number}, 
        (n2, n21, n22)::Tuple{Number,Number,Number},
    )
    ω = ω1 + ω2
    β = sfangle(ω, ω1, ω2, β1, β2)
    n′, n′1, n′2 = ninterface.((n2, n21, n22))

    wave  = Wave(ω, β, n1, n2, n′)
    wave1 = Wave(ω1, β1, n11, n21, n′1)
    wave2 = Wave(ω2, β2, n12, n22, n′2)

    Setup(wave, wave1, wave2)
end

function Base.show(io::IO, s::Setup)
    ω = s.Ω.ω
    ω1= s.Ω1.ω
    ω2= s.Ω2.ω
    λ = freq2wl(ω)
    λ1= freq2wl(ω1)
    λ2= freq2wl(ω2)

    data = [
        "Ω"  λ  ω  rad2deg(s.Ω.β)  s.Ω.n1  s.Ω.n2  s.Ω.n′;
        "Ω1" λ1 ω1 rad2deg(s.Ω1.β) s.Ω1.n1 s.Ω1.n2 s.Ω1.n′;
        "Ω2" λ2 ω2 rad2deg(s.Ω2.β) s.Ω2.n1 s.Ω2.n2 s.Ω2.n′;
    ]
    formatter = (v,i,j) -> typeof(v) <: Number ? round(v, digits=3) : v
    header = ["", "λ (nm)", "ω (s⁻¹)", "β (deg)", "n₁", "n₂", "n′"]

    println("SFGAnalysis.Setup")
    pretty_table(data; formatters=formatter, header=header)
end

# In order to not broadcast Setup and treat it as a scalar
Broadcast.broadcastable(s::Setup) = (s,)