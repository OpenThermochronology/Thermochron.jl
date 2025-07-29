# Thermochron.jl

[![DOI](osf_io_WQ2U5.svg)](https://doi.org/10.17605/OSF.IO/WQ2U5)
[![Dev][docs-dev-img]][docs-dev-url]
[![CI][ci-img]][ci-url]
[![CI (Julia nightly)][ci-nightly-img]][ci-nightly-url]
[![codecov.io][codecov-img]][codecov-url]

Open-source time-Temperature inversion of thermochronometric data. 

Implements a transdimensional Bayesian Markov chain Monte Carlo (MCMC) time-Temperature inversion with optional Simulated Annealing (e.g., [1](https://en.wikipedia.org/wiki/Simulated_annealing), [2](https://doi.org/10.1007/978-94-015-7744-1_2)) and kinetic uncertainty propagation/inversion.

Currently, this package supports the inversion of mineral helium and argon ages using a spherical Crank-Nicolson forward diffusion model following the equations of [Gallagher 1995](https://doi.org/10.1016/0012-821X(95)00197-K) and [Ketcham 2005](https://doi.org/10.2138/rmg.2005.58.11) along with the damage and annealing model of [Guenthner et al. 2013](https://doi.org/10.2475/03.2013.01) for zircon U-Th/He (ZRDAAM), the damage and annealing model of [Flowers et al. 2009](https://doi.org/10.1016/j.gca.2009.01.015) for apatite U-Th/He (RDAAM), or any constant user-specified $D_0$ and $E_a$ for any other generic He or Ar chronometer. The K-feldspar multi-diffusion domain (MDD) method of [Lovera et al. 1989](https://doi.org/10.1029/JB094iB12p17917) is also included. 

Apatite fission track age and fission track length data are supported with the annealing models of [Ketcham et al. 1999](https://doi.org/10.2138/am-1999-0903) and [Ketcham et al. 2007](https://doi.org/10.2138/am.2007.2281), while zircon fission track data are supported with the annealing model of [Yamada et al. 2007](https://doi.org/10.1016/j.chemgeo.2006.09.002) and the simultaneous-fit fanning curvilinear ZFT annealing model from [Guenthner et al. 2013](https://doi.org/10.2475/03.2013.01), discussed further by [Ketcham 2019](https://doi.org/10.1007/978-3-319-89421-8_3) (sec. 3.8, pgs. 65-67). Monazite fission track is also included with the re-fit annealing model of [Jones et al. 2021](https://doi.org/10.5194/gchron-3-89-2021). 

Diffusion-based chronometers:
| Chronometer                  | Mineral  | System         | Geometry   | Diffusion and/or annealing model(s)   |
| :---                         | :---     | :---           | :---       | :---                                  |
| [`ZirconHe`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.ZirconHe) | zircon | helium | spherical  | [`ZRDAAM`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.ZRDAAM) |
| [`ApatiteHe`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.ApatiteHe) | apatite  | helium  | spherical | [`RDAAM`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.RDAAM) |
| [`SphericalHe`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.SphericalHe) | any | helium | spherical | [`Diffusivity`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Diffusivity) (user-specified D₀, Eₐ) |
| [`PlanarHe`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.PlanarHe) | any | helium | slab | [`Diffusivity`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Diffusivity) (user-specified D₀, Eₐ) |
| [`SphericalAr`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.SphericalAr) | any | argon | spherical | [`Diffusivity`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Diffusivity) (user-specified D₀, Eₐ) |
| [`PlanarAr`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.PlanarAr) | any | argon | slab | [`Diffusivity`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Diffusivity) (user-specified D₀, Eₐ) |
| [`MultipleDomain`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.MultipleDomain) | feldspar | argon | slab/sphere| [`MDDiffusivity`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.MDDiffusivity) |

Fission track chronometers:
| Chronometer                  | Mineral  | System            | Annealing model(s)                  |
| :---                         | :---     | :---              | :---                                |
| [`ZirconFT`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.ZirconFT) | zircon | fission track age | [`Yamada2007PC`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Yamada2007PC),[`Guenthner2013FC`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Guenthner2013FC) |
| [`MonaziteFT`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.MonaziteFT) | monazite | fission track age | [`Jones2021FA`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Jones2021FA) |
| [`ApatiteFT`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.ApatiteFT) | apatite  | fission track age | [`Ketcham1999FC`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Ketcham1999FC) [`Ketcham2007FC`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Ketcham2007FC) |
| [`ZirconTrackLength`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.ZirconTrackLength) | zircon | track length | [`Yamada2007PC`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Yamada2007PC), [`Guenthner2013FC`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Guenthner2013FC) |
| [`MonaziteTrackLength`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.MonaziteTrackLength)  | monazite | track length | [`Jones2021FA`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Jones2021FA) |
| [`ApatiteTrackLength`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.ApatiteTrackLength)  | apatite | track length | [`Ketcham1999FC(:unoriented)`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Ketcham1999FC) |
| [`ApatiteTrackLengthOriented`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.ApatiteTrackLengthOriented) | apatite  | track length | [`Ketcham1999FC`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Ketcham1999FC), [`Ketcham2007FC`](https://openthermochronology.github.io/Thermochron.jl/dev/#Thermochron.Ketcham2007FC) |

Additional systems and models are expected to be added in future releases.

## Installation
Thermochron.jl is written in the open-source programming language Julia, for which installation instructions are available at [julialang.org/install](https://julialang.org/install/).
As a registered Julia package, Thermochron.jl can by installed by typing
```julia
] add Thermochron
```
at the Julia REPL, or alternatively
```julia
using Pkg
Pkg.add("Thermochron")
```

## Usage
Download an example script such as [tTinversion.jl](examples/tTinversion.jl) from the [examples](examples) folder, along with any relevant data files, and run it in your favorite Julia-connected editor or IDE (e.g. [vscode/vscodium](https://github.com/julia-vscode/julia-vscode#installing-juliavs-codevs-code-julia-extension)). A Manifest.toml is also provided in the examples folder, which you may `Pkg.instantiate` to ensure you have the same versions of all relevant packages that the example was built for.

See also the test suite for more information and examples.

## Citation
Thermochron.jl may be cited as:
> Keller, C.B., McDannell, K.T., Guenthner, W.R., and Shuster, D.L. (2022). *Thermochron.jl: Open-source time-Temperature inversion of thermochronometric data.* [10.17605/osf.io/wq2U5](https://doi.org/10.17605/osf.io/wq2U5)

## Acknowledgments
The development of Thermochron.jl has been supported in part by National Science Foundation under grant EAR-2044800.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://OpenThermochronology.github.io/Thermochron.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://OpenThermochronology.github.io/Thermochron.jl/dev/
[ci-img]: https://github.com/OpenThermochronology/Thermochron.jl/actions/workflows/CI.yml/badge.svg?branch=main
[ci-url]: https://github.com/OpenThermochronology/Thermochron.jl/actions/workflows/CI.yml
[ci-nightly-img]: https://github.com/OpenThermochronology/Thermochron.jl/workflows/CI%20(Julia%20nightly)/badge.svg
[ci-nightly-url]: https://github.com/OpenThermochronology/Thermochron.jl/actions/workflows/CI-julia-nightly.yml
[codecov-img]: https://codecov.io/gh/OpenThermochronology/Thermochron.jl/branch/main/graph/badge.svg
[codecov-url]: http://codecov.io/github/OpenThermochronology/Thermochron.jl?branch=main
