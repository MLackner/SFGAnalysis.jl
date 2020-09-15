module SFGAnalysis

using Distributions
using QuadGK

include("types.jl")
include("general.jl")
include("fresnel.jl")
include("susceptibility.jl")
include("spectra.jl")

export # types.jl
       Setup,
       # general.jl
       sfgintensity,
       sfangle,
       wl2freq,
       freq2wl,
       # fresnel.jl
       angleofrefraction,
       fresnel_x,
       fresnel_y,
       fresnel_z,
       ninterface,
       # susceptibility.jl
       susceptibility,
       effective_susceptibility_ppp,
       effective_susceptibility_ssp,
       effective_susceptibility,
       # spectra.jl
       sfspectrum,
       sfspectrum_voigt

end
