# Thermochron.jl

[![Dev][docs-dev-img]][docs-dev-url]
[![CI][ci-img]][ci-url]
[![CI (Julia nightly)][ci-nightly-img]][ci-nightly-url]
[![codecov.io][codecov-img]][codecov-url]

Open-source time-temperature inversion of thermochronometric data.

Implements a Markov chain Monte Carlo (MCMC) time-temperature inversion with optional simulated annealing.
Currently supports zircon helium data with the damage and annealing model of [Guenthner et al., 2013](https://doi.org/10.2475/03.2013.01), with additional systems and models expected to be added in future releases.

## Installation
Thermochron.jl can be installed in the same ways you would install any other registered Julia package, i.e.
```julia
] add Thermochron
```
at the Julia REPL, or alternatively
```julia
using Pkg
Pkg.add("StaticCompiler")
```

## Usage
Download an example script such as [ZrnHeInversionVartCryst.jl](examples/ZrnHeInversionVartCryst.jl) from the [examples](examples) folder, along with any relevant data files, and run it in your favorite Julia-connected editor or IDE. A Manifest.toml is also provided in the examples folder, which you may `Pkg.instantiate` to ensure you have the same versions of all relevant packages that the example was built for.

See also the test suite for more information and examples.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://brenhinkeller.github.io/Thermochron.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://brenhinkeller.github.io/Thermochron.jl/dev/
[ci-img]: https://github.com/brenhinkeller/Thermochron.jl/actions/workflows/CI.yml/badge.svg?branch=main
[ci-url]: https://github.com/brenhinkeller/Thermochron.jl/actions/workflows/CI.yml
[ci-nightly-img]: https://github.com/brenhinkeller/Thermochron.jl/workflows/CI%20(Julia%20nightly)/badge.svg
[ci-nightly-url]: https://github.com/brenhinkeller/Thermochron.jl/actions/workflows/CI-julia-nightly.yml
[codecov-img]: http://codecov.io/github/brenhinkeller/Thermochron.jl/coverage.svg?branch=main
[codecov-url]: http://codecov.io/github/brenhinkeller/Thermochron.jl?branch=main
