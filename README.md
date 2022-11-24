# Thermochron.jl

[![DOI](osf_io_WQ2U5.svg)](https://doi.org/10.17605/OSF.IO/WQ2U5)
[![Dev][docs-dev-img]][docs-dev-url]
[![CI][ci-img]][ci-url]
[![CI (Julia nightly)][ci-nightly-img]][ci-nightly-url]
[![codecov.io][codecov-img]][codecov-url]

Open-source time-Temperature inversion of thermochronometric data.

Implements a transdimensional Markov chain Monte Carlo (MCMC) time-Temperature inversion with optional Simulated Annealing (e.g., [1](https://en.wikipedia.org/wiki/Simulated_annealing), [2](https://doi.org/10.1007/978-94-015-7744-1_2)).
Currently supports zircon helium data with the damage and annealing model of [Guenthner et al., 2013](https://doi.org/10.2475/03.2013.01), with additional systems and models expected to be added in future releases.

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
Download an example script such as [ZrnHeInversionVartCryst.jl](examples/ZrnHeInversionVartCryst.jl) from the [examples](examples) folder, along with any relevant data files, and run it in your favorite Julia-connected editor or IDE. A Manifest.toml is also provided in the examples folder, which you may `Pkg.instantiate` to ensure you have the same versions of all relevant packages that the example was built for.

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
[codecov-img]: http://codecov.io/github/OpenThermochronology/Thermochron.jl/coverage.svg?branch=main
[codecov-url]: http://codecov.io/github/OpenThermochronology/Thermochron.jl?branch=main
