var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Thermochron","category":"page"},{"location":"#Thermochron","page":"Home","title":"Thermochron","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for the Thermochron.jl package.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Thermochron]","category":"page"},{"location":"#Thermochron.ApatiteHe","page":"Home","title":"Thermochron.ApatiteHe","text":"ApatiteHe(T=Float64;\n    age::Number = T(NaN),\n    age_sigma::Number = T(NaN),\n    offset::Number = zero(T),\n    r::Number, \n    dr::Number = one(T), \n    U238::Number, \n    Th232::Number, \n    Sm147::Number = zero(T), \n    U238_matrix::Number = zero(T), \n    Th232_matrix::Number = zero(T), \n    Sm147_matrix::Number = zero(T), \n    volumeweighting::Symbol = :spherical,\n    agesteps::AbstractRange,\n)\n\nConstruct an ApatiteHe chronometer representing an apatite with a raw  helium age of age ± age_sigma [Ma], a  radius of r [μm], and uniform  U, Th and Sm concentrations specified by U238, Th232, and Sm147 [PPMW].  A present day U-235/U-238 ratio of 1/137.818 is assumed.\n\nSpatial discretization follows a radius step of dr [μm], and temporal discretization follows the age steps specified by the agesteps range, in Ma.\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.Diffusivity","page":"Home","title":"Thermochron.Diffusivity","text":"Diffusivity(\n    D0::T = 59.98               # [cm^2/sec] Maximum diffusion coefficient\n    D0_logsigma::T = log(2)/2   # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)\n    Ea::T = 205.94              # [kJ/mol] Activation energy\n    Ea_logsigma::T = log(2)/2   # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)\n)\n\nA generic diffusivity model, with user-specified D0 and Ea. Default values are appropriate for argon in k-feldspar.\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.Guenthner2013FC","page":"Home","title":"Thermochron.Guenthner2013FC","text":"Guenthner2013FC(\n    C0 = 6.24534         # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon\n    C1 = -0.11977        # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon\n    C2 = -314.93688      # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon\n    C3 = -14.2868        # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon\n    alpha = -0.057206897 # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon\n    l0 = 11.17           # [um] Initial track length\n    l0_sigma = 0.051     # [um] Initial track length uncertainty\n)\n\nFanning Curvilinear zircon annealing model with simplified Box-Cox transform  from Guenthner et al. 2013 (doi: 10.2475/03.2013.01)\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.Jones2021FA","page":"Home","title":"Thermochron.Jones2021FA","text":"Jones2021FA(\n    C0 = 1.374           # Annealing parameter\n    C1 = -4.192812e-5    # Annealing parameter\n    C2 = -22.70885029    # Annealing parameter\n    C3 = 0.0             # Annealing parameter\n    l0 = 10.60           # [um] Initial track length\n    l0_sigma = 0.19      # [um] Initial track length uncertainty\n)\n\nParallel Arrhenius monazite annealing model modified from Jones et al., 2021 (doi: 10.5194/gchron-3-89-2021)\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.Ketcham1999FC","page":"Home","title":"Thermochron.Ketcham1999FC","text":"Ketcham1999FC(\n    C0::T = -19.844     # \"Simultaneous fit\" from Ketcham et al. 1999 apatite\n    C1::T = 0.38951     # \"Simultaneous fit\" from Ketcham et al. 1999 apatite\n    C2::T = -51.253     # \"Simultaneous fit\" from Ketcham et al. 1999 apatite\n    C3::T = -7.6423     # \"Simultaneous fit\" from Ketcham et al. 1999 apatite\n    alpha::T = -0.12327 # Box-Cox transform parameter\n    beta::T = -11.988   # Box-Cox transform parameter\n    l0::T = 16.38       # [um] Initial track length\n    l0_sigma::T = 0.09  # [um] Initial track length unertainty\n)\n\nFanning Curvilinear apatite annealing model from Ketcham, 1999 (doi: 10.2138/am-1999-0903)\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.Ketcham2007FC","page":"Home","title":"Thermochron.Ketcham2007FC","text":"Ketcham2007FC(\n    C0 = 0.39528     # \"Simultaneous fit\" from Ketcham et al. 2007 apatite\n    C1 = 0.01073     # \"Simultaneous fit\" from Ketcham et al. 2007 apatite\n    C2 = -65.12969   # \"Simultaneous fit\" from Ketcham et al. 2007 apatite\n    C3 = -7.91715    # \"Simultaneous fit\" from Ketcham et al. 2007 apatite\n    alpha = 0.04672  # \"Simultaneous fit\" from Ketcham et al. 2007 apatite\n    l0 = 16.38       # [um] Initial track length\n    l0_sigma = 0.09  # [um] Initial track length unertainty\n)\n\nFanning Curvilinear apatite annealing model with simplified Box-Cox transform  from Ketcham, 2007 (doi: 10.2138/am.2007.2281)\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.MDDiffusivity","page":"Home","title":"Thermochron.MDDiffusivity","text":"MDDiffusivity(\n    D0::NTuple{N,T}             # [cm^2/sec] Maximum diffusivity\n    D0_logsigma::NTuple{N,T}    # [unitless] log uncertainty (default = 1/2 = a factor of ℯ two-sigma)\n    Ea::T                       # [kJ/mol] Activation energy\n    Ea_logsigma::T              # [unitless] log uncertainty (default = 1/2 = a factor of ℯ two-sigma)\n)\n\nMultiple diffusivities for multiple domains\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.MultipleDomain","page":"Home","title":"Thermochron.MultipleDomain","text":"MultipleDomain(T=Float64, C=PlanarAr;\n    age::AbstractVector,\n    age_sigma::AbstractVector,\n    fraction_experimental::AbstractVector,\n    fraction_experimental_sigma::Number=T(0.01),\n    tsteps_experimental::AbstractVector,\n    Tsteps_experimental::AbstractVector,\n    fit::AbstractVector,\n    offset::Number = zero(T),\n    fuse::Bool = true,\n    volume_fraction::AbstractVector,\n    r::Number = 100,\n    dr::Number = one(T),\n    K40::Number = 16.34, \n    agesteps::AbstractRange,\n    tsteps::AbstractRange=reverse(agesteps),\n)\n\nConstruct a MultipleDomain diffusion chronometer given an observed argon release spectrum, degassing schedule, where domain is represented by a  PlanarAr or SphericalAr chronometer.\n\nDomain diffusivity and volume parameters must be supplied as vectors Ea [kJ/mol], lnD0a2 [log(1/s)], and volume_fraction [unitless] obtained by separately fitting the release spectrum (the former two as an MDDiffusivity object).\n\nSee also: MDDiffusivity, PlanarAr, SphericalAr, degas!\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.PlanarAr","page":"Home","title":"Thermochron.PlanarAr","text":"PlanarAr(T=Float64;\n    age::Number = T(NaN),\n    age_sigma::Number = T(NaN),\n    offset::Number = zero(T),\n    r::Number, \n    dr::Number = one(T), \n    K40::Number=16.34, \n    agesteps::AbstractRange,\n)\n\nConstruct an PlanarAr chronometer representing a mineral with a raw  argon age of age ± age_sigma [Ma], a uniform diffusivity, a radius of r [μm], and uniform K-40 concentrations specified by K40 [PPM].\n\nSpatial discretization follows a halfwidth step of dr [μm], and temporal discretization follows the age steps specified by the agesteps range, in Ma.\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.PlanarHe","page":"Home","title":"Thermochron.PlanarHe","text":"PlanarHe(T=Float64;\n    age::Number = T(NaN),\n    age_sigma::Number = T(NaN),\n    offset::Number = zero(T),\n    stoppingpower::Number = T(1.189),\n    r::Number, \n    dr::Number = one(T), \n    U238::Number, \n    Th232::Number, \n    Sm147::Number = zero(T), \n    U238_matrix::Number = zero(T), \n    Th232_matrix::Number = zero(T), \n    Sm147_matrix::Number = zero(T), \n    agesteps::AbstractRange,\n)\n\nConstruct an PlanarHe chronometer representing a mineral with a raw  helium age of age ± age_sigma [Ma], uniform diffusivity,  a halfwidth of r [μm], and uniform U, Th and Sm concentrations specified by U238, Th232, and Sm147 [PPM]. (A present day U-235/U-238  ratio of 1/137.818 is assumed)\n\nSpatial discretization follows a halfwidth step of dr [μm], and temporal discretization follows the age steps specified by the agesteps range, in Ma.\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.RDAAM","page":"Home","title":"Thermochron.RDAAM","text":"RDAAM(\n    D0L::T=0.6071               # [cm^2/sec] Maximum diffusivity\n    D0L_logsigma::T=log(2)/2    # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)\n    EaL::T=122.3                # [kJ/mol] Activation energy\n    EaL_logsigma::T=log(2)/4    # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)\n    EaTrap::T=34.0              # [kJ/mol] Activation energy\n    EaTrap_logsigma::T=log(2)/4 # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)\n    psi::T=1e-13                # empirical polynomial coefficient\n    omega::T=1e-22              # empirical polynomial coefficient\n    etaq::T=0.91                # Durango ηq\n    rhoap::T=3.19               # Density of apatite [g/cm3]\n    L::T=0.000815               # Etchable fission track half-length, cm\n    lambdaf::T=8.46e-17         # \n    lambdaD::T=1.55125e-10      # \n    beta::T=0.04672             # Apatite annealing parameter. Also caled alpha, but equivalent to beta in ZRDAAM\n    C0::T=0.39528               # Apatite annealing parameter\n    C1::T=0.01073               # Apatite annealing parameter\n    C2::T=-65.12969 - LOG_SEC_MYR # Apatite annealing parameter. Includes conversion factor from seconds to Myr for dt, in addition to traditional C2 value\n    C3::T=-7.91715              # Apatite annealing parameter\n    rmr0::T=0.83                # Damage conversion parameter\n    rmr0_sigma::T=0.15          # Damage conversion parameter uncertainty\n    kappa::T=1.04-0.83          # Damage conversion parameter\n    kappa_rmr0::T=1.04          # Damage conversion parameter (the sum of kappa and rmr0)\n)\n\nApatite Radiation Damage Accumulation and Annealing Model (RDAAM) of Flowers et al. 2009 (doi: 10.1016/j.gca.2009.01.015)\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.SphericalAr","page":"Home","title":"Thermochron.SphericalAr","text":"SphericalAr(T=Float64;\n    age::Number = T(NaN),\n    age_sigma::Number = T(NaN),\n    offset::Number = zero(T),\n    r::Number, \n    dr::Number = one(T), \n    K40::Number=16.34, \n    agesteps::AbstractRange,\n)\n\nConstruct an SphericalAr chronometer representing a mineral with a raw  argon age of age ± age_sigma [Ma], a uniform diffusivity, a radius of r [μm], and uniform K-40 concentrations specified by K40 [PPM].\n\nSpatial discretization follows a radius step of dr [μm], and temporal discretization follows the age steps specified by the agesteps range, in Ma.\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.SphericalHe","page":"Home","title":"Thermochron.SphericalHe","text":"SphericalHe(T=Float64;\n    age::Number = T(NaN),\n    age_sigma::Number = T(NaN),\n    offset::Number = zero(T),\n    stoppingpower::Number = T(1.189),\n    r::Number, \n    dr::Number = one(T), \n    U238::Number, \n    Th232::Number, \n    Sm147::Number = zero(T), \n    U238_matrix::Number = zero(T), \n    Th232_matrix::Number = zero(T), \n    Sm147_matrix::Number = zero(T), \n    agesteps::AbstractRange,\n)\n\nConstruct a SphericalHe chronometer representing a mineral with a raw  helium age of age ± age_sigma [Ma], uniform diffusivity, a radius of r [μm], and uniform U, Th and Sm concentrations specified by U238, Th232, and Sm147 [PPM]. (A present day U-235/U-238  ratio of 1/137.818 is assumed)\n\nSpatial discretization follows a radius step of dr [μm], and temporal discretization follows the age steps specified by the agesteps range, in Ma.\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.Yamada2007PC","page":"Home","title":"Thermochron.Yamada2007PC","text":"Yamada2007PC(\n    c0p = -63.37     # Yamada et al. 2007 zircon\n    c1p = 0.212      # Yamada et al. 2007 zircon\n    bp = 43.00       # Yamada et al. 2007 zircon\n    l0 = 11.17       # [um] effective initial track length (μmax)\n    l0_sigma = 0.051 # [um] effective initial track length uncertainty (σ)\n)\n\nParallel Curvilinear zircon annealing model of Yamada, 2007  (doi: 10.1016/j.chemgeo.2006.09.002)\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.ZRDAAM","page":"Home","title":"Thermochron.ZRDAAM","text":"ZRDAAM(\n    DzD0::T = 193188.0          # [cm^2/sec] Maximum diffusivity, crystalline endmember\n    DzD0_logsigma::T=log(2)/2   # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)\n    DzEa::T=165.0               # [kJ/mol] Activation energy, crystalline endmember\n    DzEa_logsigma::T=log(2)/4   # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)\n    DN17D0::T = 6.367E-3        # [cm^2/sec] Maximum diffusivity, amorphous endmember\n    DN17D0_logsigma::T=log(2)/2 # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)\n    DN17Ea::T=71.0              # [kJ/mol] Activation energy, amorphous endmember\n    DN17Ea_logsigma::T=log(2)/4 # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)\n    lint0::T=45920.0            # [nm]\n    SV::T=1.669                 # [1/nm]\n    Bα::T=5.48E-19              # Amorphous material produced per alpha decay [g/alpha]\n    Phi::T=3.0                  # unitless\n    beta::T=-0.05721            # Zircon anealing parameter\n    C0::T=6.24534               # Zircon anealing parameter\n    C1::T=-0.11977              # Zircon anealing parameter\n    C2::T=-314.937 - LOG_SEC_MYR # Zircon anealing parameter. Includes conversion factor from seconds to Myr for dt (for performance), in addition to traditional C2 value\n    C3::T=-14.2868              # Zircon anealing parameter\n    rmin::T=0.2                 # Damage conversion parameter\n    rmin_sigma::T=0.15          # Damage conversion parameter uncertainty\n)\n\nZircon Radiation Damage Accumulation and Annealing Model (ZRDAAM) of Guenthner et al. 2013 (doi: 10.2475/03.2013.01)\n\n\n\n\n\n","category":"type"},{"location":"#Thermochron.ZirconHe","page":"Home","title":"Thermochron.ZirconHe","text":"ZirconHe(T=Float64;\n    age::Number = T(NaN),\n    age_sigma::Number = T(NaN),\n    T(offset),\n    r::Number = one(T),\n    dr::Number, \n    U238::Number,\n    Th232::Number,\n    Sm147::Number = zero(T),\n    U238_matrix::Number = zero(T), \n    Th232_matrix::Number = zero(T), \n    Sm147_matrix::Number = zero(T), \n    volumeweighting::Symbol=:cylindrical,\n    agesteps::AbstractRange,\n)\n\nConstruct a ZirconHe chronometer representing a zircon with a raw  helium age of age ± age_sigma [Ma], a  radius of r [μm], and uniform  U, Th and Sm concentrations specified by U238, Th232, and Sm147 [PPMw].  A present day U-235/U-238 ratio of 1/137.818 is assumed.\n\nSpatial discretization follows a radius step of dr μm, and temporal discretization follows the age steps specified by the agesteps range, in Ma.\n\n\n\n\n\n","category":"type"},{"location":"#StatGeochemBase.image_from_paths-Tuple{Thermochron.TTResult}","page":"Home","title":"StatGeochemBase.image_from_paths","text":"image_from_paths(tT::TTResult; \n    xresolution::Int=1800, \n    yresolution::Int=1200, \n    xrange = nanextrema(tT.tpointdist), \n    yrange = nanextrema(tT.Tpointdist), \n)\n\nProduce a 2d image (histogram) of path densities given a TTResult\n\njulia> imgcounts, xq, yq = image_from_paths(tT; xrange=boundary.agepoints, yrange=boundary.T₀)\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.MCMC-Union{Tuple{T}, Tuple{NamedTuple, NamedTuple, Boundary{T}}, Tuple{NamedTuple, NamedTuple, Boundary{T}, Constraint{T}}, Tuple{NamedTuple, NamedTuple, Boundary{T}, Constraint{T}, DetailInterval{T}}} where T<:AbstractFloat","page":"Home","title":"Thermochron.MCMC","text":"MCMC(chrons::Vector{<:Chronometer}, model::NamedTuple, npoints::Int, path.agepoints::Vector, path.Tpoints::Vector, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])\n\nMarkov chain Monte Carlo time-Temperature inversion of the thermochronometric chrons  specified as a vector chrons of Chronometer objects (ZirconHe, ApatiteHe,  ApatiteFT, etc.) and model parameters specified by the named tuple model,  with variable diffusion kinetics.\n\nReturns a TTResult object containing posterior time-temperature paths,\n\nSee also MCMC_varkinetics for a variant with variable diffusion kinetics.\n\nExamples\n\ntT = MCMC(chrons::NamedTuple, model::NamedTuple, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.MCMC_varkinetics-Union{Tuple{T}, Tuple{NamedTuple, NamedTuple, Boundary{T}}, Tuple{NamedTuple, NamedTuple, Boundary{T}, Constraint{T}}, Tuple{NamedTuple, NamedTuple, Boundary{T}, Constraint{T}, DetailInterval{T}}} where T<:AbstractFloat","page":"Home","title":"Thermochron.MCMC_varkinetics","text":"MCMC_varkinetics(chrons::Vector{<:Chronometer}, model::NamedTuple, npoints::Int, path.agepoints::Vector, path.Tpoints::Vector, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])\n\nMarkov chain Monte Carlo time-Temperature inversion of the thermochronometric chrons  specified as a vector chrons of Chronometer objects (ZirconHe, ApatiteHe,  ApatiteFT, etc.) and model parameters specified by the named tuple model,  with variable diffusion kinetics.\n\nReturns a TTResult object containing posterior time-temperature paths, and a KineticResult object containing the posterior kinetic parameters.\n\nSee also MCMC for a variant with constant diffusion kinetics.\n\nExamples\n\ntT, kinetics = MCMC_varkinetics(chrons::NamedTuple, model::NamedTuple, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.anneal","page":"Home","title":"Thermochron.anneal","text":"ρᵣ = anneal(dt::Number, tsteps::Vector, Tsteps::Matrix, [model::DiffusivityModel=ZRDAAM()])\n\nZirconHe damage annealing model as in Guenthner et al. 2013 (AJS)\n\n\n\n\n\n","category":"function"},{"location":"#Thermochron.anneal!-Union{Tuple{C}, Tuple{T}, Tuple{Vector{<:Union{ApatiteFT{T}, ApatiteHe{T}, ApatiteTrackLength{T}, MonaziteFT{T}, MonaziteTrackLength{T}, MultipleDomain{T, PlanarAr{T}}, MultipleDomain{T, SphericalAr{T}}, PlanarAr{T}, PlanarHe{T}, SphericalAr{T}, SphericalHe{T}, ZirconFT{T}, ZirconHe{T}, ZirconTrackLength{T}}}, Type{C}, AbstractRange{T}, AbstractVector{T}, Thermochron.DiffusivityModel{T}}} where {T<:AbstractFloat, C<:Thermochron.HeliumSample}","page":"Home","title":"Thermochron.anneal!","text":"anneal!(data::Vector{<:Chronometer}, ::Type{<:HeliumSample}, tsteps, Tsteps, dm::DiffusivityModel)\nanneal!(mineral::ZirconHe, Tsteps, dm::ZRDAAM)\nanneal!(mineral::ApatiteHe, Tsteps, dm::RDAAM)\nanneal!(ρᵣ::Matrix, dt::Number, tsteps, Tsteps, [dm::DiffusivityModel=ZRDAAM()])\n\nIn-place version of anneal\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.chronometers-Tuple{Any, Any}","page":"Home","title":"Thermochron.chronometers","text":"chronometers([T=Float64], data, model)\n\nConstruct a vector of Chronometer objects given a dataset data and model parameters model\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.lcmod-Tuple{ApatiteTrackLength}","page":"Home","title":"Thermochron.lcmod","text":"lcmod(l, θ)\n\nCalculate the model c-axis equivalent length (\"lc,mod\") given a measured \"confined\" fission track length l [microns] and angle from the c-axis  θ [degrees] following the approach of Donelick et al. 1999  (doi: 10.2138/am-1999-0902) \n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.modelage-Tuple{ZirconHe, AbstractVector, AbstractMatrix, ZirconHeliumModel}","page":"Home","title":"Thermochron.modelage","text":"modelage(mineral::ZirconHe, Tsteps, [ρᵣ], dm::ZRDAAM)\nmodelage(mineral::ApatiteHe, Tsteps, [ρᵣ], dm::RDAAM)\nmodelage(mineral::SphericalHe, Tsteps, dm::Diffusivity)\nmodelage(mineral::PlanarHe, Tsteps, dm::Diffusivity)\n\nCalculate the predicted bulk U-Th/He age of a zircon, apatite, or other mineral that has experienced a given t-T path (specified by mineral.tsteps for time and Tsteps for temperature, at a time resolution of step(mineral.tsteps))  using a Crank-Nicolson diffusion solution for a spherical (or planar slab) grain of radius (or halfwidth) mineral.r at spatial resolution mineral.dr.\n\nSpherical implementation based on the the Crank-Nicolson solution for diffusion out of a spherical mineral crystal in Ketcham, 2005 (doi: 10.2138/rmg.2005.58.11).\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.modelage-Union{Tuple{T}, Tuple{SphericalAr{T}, AbstractVector{T}, Diffusivity{T}}} where T<:AbstractFloat","page":"Home","title":"Thermochron.modelage","text":"modelage(mineral::SphericalAr, Tsteps, dm::Diffusivity)\nmodelage(mineral::PlanarAr, Tsteps, dm::Diffusivity)\n\nCalculate the precdicted bulk K/Ar age of a mineral that has experienced a given  t-T path (specified by mineral.tsteps for time and Tsteps for temperature,  at a time resolution of step(mineral.tsteps)) using a Crank-Nicholson diffusion  solution for a spherical (or planar slab) grain of radius (or halfwidth ) mineral.r  at spatial resolution mineral.dr.\n\nSpherical implementation based on the the Crank-Nicolson solution for diffusion out of a spherical mineral crystal in Ketcham, 2005 (doi: 10.2138/rmg.2005.58.11).\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.modelage-Union{Tuple{T}, Tuple{ZirconFT{T}, AbstractVector, Thermochron.ZirconAnnealingModel{T}}} where T<:AbstractFloat","page":"Home","title":"Thermochron.modelage","text":"modelage(mineral::ZirconFT, Tsteps, am::ZirconAnnealingModel)\nmodelage(mineral::MonaziteFT, Tsteps, am::MonaziteAnnealingModel)\nmodelage(mineral::ApatiteFT, Tsteps, am::ApatiteAnnealingModel)\n\nCalculate the precdicted fission track age of an apatite that has experienced a given  t-T path (specified by mineral.tsteps for time and Tsteps for temperature, at a time resolution of step(mineral.tsteps)) and given annealing model parameters am.\n\nPossible annealing model types and the references for the equations  which they respetively implement include    Ketcham1999FC       Fanning Curvilinear apatite model of Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)   Ketcham2007FC       Fanning Curvilinear apatite model of Ketcham et al. 2007 (doi: 10.2138/am.2007.2281)   Yamada2007PC        Parallel Curvilinear zircon model of Yamada et al. 2007 (doi: 10.1016/j.chemgeo.2006.09.002)   Guenthner2013FC     Fanning Curvilinear zircon model of Guenthner et al. 2013 (doi: 10.2475/03.2013.01)   Jones2021FA         Fanning Arrhenius (Fanning Linear) model adapted from Jones et al. 2021 (doi: 10.5194/gchron-3-89-2021)\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.modellength-Union{Tuple{T}, Tuple{ApatiteTrackLength{T}, AbstractVector, Thermochron.ApatiteAnnealingModel{T}}} where T<:AbstractFloat","page":"Home","title":"Thermochron.modellength","text":"modellength(track::ApatiteTrackLength, Tsteps, am::ApatiteAnnealingModel)\n\nCalculate the predicted mean and standard deviation of the distribution of fission   track lengths of an apatite that has experienced a given t-T path (specified by  track.tsteps for time and Tsteps for temperature, at a time resolution of  step(mineral.tsteps)) and given annealing model parameters am.\n\nPossible annealing model types and the references for the equations  which they respetively implement include    Ketcham1999FC       Fanning Curvilinear apatite model of Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)   Ketcham2007FC       Fanning Curvilinear apatite model of Ketcham et al. 2007 (doi: 10.2138/am.2007.2281)\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.reltrackdensity-Tuple{Number, Number, Thermochron.ZirconAnnealingModel}","page":"Home","title":"Thermochron.reltrackdensity","text":"reltrackdensity(t, T, am::AnnealingModel)\n\nCalculate the relative track density ρ corresponding to a given  relative track length r \n\nFollows the relations of Ketcham et al. (2000), equations 7a and 7b  for apatite and Tagami et al. (1999) for zircon\n\nSee also: reltracklength.\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.rmr0fromcl-Tuple{Any}","page":"Home","title":"Thermochron.rmr0fromcl","text":"rmr0fromcl(Cl)\n\nCalculate rmr0 as a function of chlorine content Cl [APFU] for  \"multikinetic\" apatite fission track following the relation (Fig. 7a)  of Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)\n\nrmr0 = 1 - exp(2.107(1 - abs(Cl - 1)) - 1.834)\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.rmr0fromdpar-Tuple{Any}","page":"Home","title":"Thermochron.rmr0fromdpar","text":"rmr0fromdpar(dpar)\n\nCalculate rmr0 as a function of dpar for \"multikinetic\" apatite  fission track following the relation (Fig. 7b) of Ketcham et al. 1999 (doi: 10.2138/am-1999-0903)\n\nrmr0 = 1 - exp(0.647(dpar-1.75) - 1.834)\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.rmr0model","page":"Home","title":"Thermochron.rmr0model","text":"rmr0model(F, Cl, OH, Mn=0, Fe=0, others=0)\n\nCalculate rmr0 as a function of composition (specified in terms of atoms per fomula unit, or APFU) for \"multikinetic\" apatite fission  track thermochronology.\n\nImplements equation 11 from Ketcham et al. 2007  (doi: 10.2138/am.2007.2281)\n\nrmr0 = (-0.0495 -0.0348F +0.3528|Cl - 1| +0.0701|OH - 1| \n        -0.8592Mn -1.2252Fe -0.1721Others)^0.1433\n\n\n\n\n\n","category":"function"},{"location":"#Thermochron.simannealT-Tuple{Integer, Number, Number}","page":"Home","title":"Thermochron.simannealT","text":"simannealT(n::Integer, Tₐ::Number, λₐ::Number)\n\nTo avoid getting stuck in local optima, increase the probability  of accepting new proposals at higher annealing \"temperature\"\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.slabsphereintersectiondensity-Union{Tuple{T}, Tuple{AbstractVector{T}, T, T}} where T<:AbstractFloat","page":"Home","title":"Thermochron.slabsphereintersectiondensity","text":"slabsphereintersectiondensity(redges::Vector, ralpha, d)\n\nCalculate the fractional intersection density of an alpha stopping sphere of radius ralpha with each concentric slab-shell (with edges redges and relative volumes rvolumes) of a planar slab crystal where the two are separated by distance d\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.slabsphereintersectionfraction-Union{Tuple{T}, Tuple{T, T, T}} where T<:AbstractFloat","page":"Home","title":"Thermochron.slabsphereintersectionfraction","text":"slabsphereintersectionfraction(rₚ, rₛ, d)\n\nCalculate the fraction of the surface area of a sphere s with radius rₛ that intersects the interior of a planar slab p of halfwidth rₚ if the two are separated by distance d.\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.sphereintersectiondensity-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{T}, T, T}} where T<:AbstractFloat","page":"Home","title":"Thermochron.sphereintersectiondensity","text":"sphereintersectiondensity(redges::Vector, rvolumes::Vector, ralpha, d)\n\nCalculate the volume-nomalized fractional intersection density of an alpha stopping sphere of radius ralpha with each concentric shell (with shell edges redges and relative volumes rvolumes) of a spherical crystal where the two are separated by distance d\n\n\n\n\n\n","category":"method"},{"location":"#Thermochron.sphereintersectionfraction-Union{Tuple{T}, Tuple{T, T, T}} where T<:AbstractFloat","page":"Home","title":"Thermochron.sphereintersectionfraction","text":"sphereintersectionfraction(r₁, r₂, d)\n\nCalculate the fraction of the surface area of a sphere s₂ with radius r₂ that intersects the interior of a sphere s₁ of radius r₁ if the two are separated by distance d.\n\n\n\n\n\n","category":"method"}]
}
