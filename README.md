# SFGAnalysis 
[![Build Status](https://travis-ci.com/MLackner/SFGAnalysis.jl.svg?branch=master)](https://travis-ci.com/MLackner/SFGAnalysis.jl) [![Coverage](https://codecov.io/gh/MLackner/SFGAnalysis.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MLackner/SFGAnalysis.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://MLackner.github.io/SFGAnalysis.jl/dev)

## Effective Susceptibilities

Let us replicate the data in Figure 4 from https://doi.org/10.1021/jp067062h (*J.
Phys. Chem. C 2007, 111, 8716-8725*).

First we have to define the correct experimental settings. Internally the
package works with frequencies `ω`. We can convert wavelengths to frequencies
easily with `wl2freq`. Converting the other way around works with `freq2wl`. 

```julia
using SFGAnalysis

ω1 = wl2freq(532)
ω2 = wl2freq(3333)
β1 = deg2rad(62)
β2 = deg2rad(53)
n1  = 1
n11 = 1
n12 = 1
n2  = 1.33
n21 = 1.33
n22 = 1.33
R = 1.7
```

In the next step we put these parameter into a `Setup` object so that we can
handle these parameters more easily. The wavelength of the sum frequency as well
as the reflected angle and the refractive index of the boundary layer are
calculated automatically via `ω = ω1 + ω2`, the `sfangle` function and the
`ninterface` function respectively. 

```julia
s = Setup((ω1, ω2), (β1[1], β2[1]), (n1, n11, n12), (n2, n21, n22))

SFGAnalysis.Setup
===== ========= ====================== ========= ===== ====== ========
        λ (nm)                ω (1/s)   β (deg)    n₁     n₂      n′  
===== ========= ====================== ========= ===== ====== ========
   Ω   458.773   6.534663899679439e14    60.615   1.0   1.33   1.149  
  Ω1     532.0   5.635196578947368e14      62.0   1.0   1.33   1.149  
  Ω2    3333.0   8.994673207320731e13      53.0   1.0   1.33   1.149  
===== ========= ====================== ========= ===== ====== ========
```

In the next step we are calculating the effective susceptibilities of the methyl
group for tilt angles between 0 and 90 degrees. Note that we append a factor of
`sec(β)` to the effective susceptibilities just like in the reference paper.

```julia
θ = range(0, π/2, length=91)

χ_eff_ssp_ss = effective_susceptibility.(θ, s, :ssp; pointgroup=:c3v, mode=:ss, R=R) * sec(s.Ω.β)
χ_eff_ssp_as = effective_susceptibility.(θ, s, :ssp; pointgroup=:c3v, mode=:as) * sec(s.Ω.β)
χ_eff_ppp_ss = effective_susceptibility.(θ, s, :ppp; pointgroup=:c3v, mode=:ss, R=R) * sec(s.Ω.β)
χ_eff_ppp_as = effective_susceptibility.(θ, s, :ppp; pointgroup=:c3v, mode=:as) * sec(s.Ω.β)
```

Finally, we can plot the effective susceptibilities and we should be able to see the expected result.

```julia
using Plots

θdeg = rad2deg.(θ)
plt = plot(xlabel="Tilt Angle (deg)", ylabel="Effective Susceptibility")
plot!(θdeg, χ_eff_ssp_ss, label="ssp ss")
plot!(θdeg, χ_eff_ssp_as, label="ssp as")
plot!(θdeg, χ_eff_ppp_ss, label="ppp ss")
plot!(θdeg, χ_eff_ppp_as, label="ppp as")
```

![example](https://user-images.githubusercontent.com/8495596/84396387-9e0e4680-abfe-11ea-8d54-cd1fa7a246a8.png)

## Angle Distributions

To input a distribution of tilt angles instead of a fixed angle we can use the
`Distributions.jl` package. In this example we use a distribution of the tilt
angle with the central angle being `0` radians and a standard deviation of `π/4`.

```julia
using Distributions

dist = Normal(0, π/4)
χ_eff = effective_susceptibility(dist, s, :ssp; pointgroup=:c3v, mode=:as)
```
