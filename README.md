# Thermochron.jl

[![DOI](osf_io_WQ2U5.svg)](https://doi.org/10.17605/OSF.IO/WQ2U5)
[![Dev][docs-dev-img]][docs-dev-url]
[![CI][ci-img]][ci-url]
[![CI (Julia nightly)][ci-nightly-img]][ci-nightly-url]
[![codecov.io][codecov-img]][codecov-url]

Open-source time-Temperature inversion of thermochronometric data.

Implements a transdimensional Bayesian Markov chain Monte Carlo (MCMC) time-Temperature inversion with optional Simulated Annealing (e.g., [1](https://en.wikipedia.org/wiki/Simulated_annealing), [2](https://doi.org/10.1007/978-94-015-7744-1_2)) and kinetic uncertainty propagation/inversion.

Currently, this package supports the inversion of mineral helium and argon ages using a spherical Crank-Nicholson forward diffusion model following the equations of [Gallagher, 1995](https://doi.org/10.1016/0012-821X(95)00197-K) and [Ketcham, 2005](https://doi.org/10.2138/rmg.2005.58.11) along with the damage and annealing model of [Guenthner et al., 2013](https://doi.org/10.2475/03.2013.01) for zircon helium (ZRDAAM), the damage and annealing model of [Flowers et al. 2009](https://doi.org/10.1016/j.gca.2009.01.015) for apatite helium (RDAAM), or any constant user-specified $D_0$ and $E_a$ for any other generic He or Ar chronometer.

Apatite fission track age and fission track length data are supported with the annealing models of [Ketcham et al. 1999](https://doi.org/10.2138/am-1999-0903) and [Ketcham et al. 2007](https://doi.org/10.2138/am.2007.2281), while zircon fission track data are supported with the annealing model of [Yamada et al. 2007](https://doi.org/10.1016/j.chemgeo.2006.09.002) and the simultaneous-fit fanning curvilinear ZFT annealing model from [Guenthner et al., 2013](https://doi.org/10.2475/03.2013.01), discussed further by [Ketcham, 2019](https://doi.org/10.1007/978-3-319-89421-8_3) (sec. 3.8, pgs. 65-67). Monazite fission track is also included with the re-fit annealing model of [Jones et al. 2021](https://doi.org/10.5194/gchron-3-89-2021). 

Additional systems and models are expected to be added in future releases.

## Installation
Thermochron.jl can be installed in the same ways you would install any other registered Julia package, i.e.
```julia
] add Thermochron
```
at the Julia REPL, or alternatively
```julia
using Pkg
Pkg.add("Thermochron")
```

## Usage
Download an example script such as [tTinversion.jl](examples/tTinversion.jl) from the [examples](examples) folder, along with any relevant data files, and run it in your favorite Julia-connected editor or IDE. A Manifest.toml is also provided in the examples folder, which you may `Pkg.instantiate` to ensure you have the same versions of all relevant packages that the example was built for.

See also the test suite for more information and examples.

## Citation
When using model results derived from Thermochron.jl, you may cite this package as:
> Keller, C.B., McDannell, K.T., Guenthner, W.R., and Shuster, D.L. (2022). *Thermochron.jl: Open-source time-Temperature inversion of thermochronometric data.* https://doi.org/10.17605/osf.io/wq2U5

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
