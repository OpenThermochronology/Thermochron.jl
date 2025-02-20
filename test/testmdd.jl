datapath = joinpath("..", "examples", "ol13.csv")
ds = importdataset(datapath, importas=:Tuple)

agesteps = 995:-10.:5
tsteps = reverse(agesteps)
Tsteps = range(250, 0., length(tsteps))

mdd = MultipleDomain(;
    age = ds.age_Ma,
    age_sigma = ds.age_sigma_Ma,
    fraction_released = ds.fraction_degassed,
    tsteps_degassing = cumsum(ds.time_s),
    Tsteps_degassing = ds.temperature_C,
    fit = ds.fit,
    volume_fraction = ds.volume_fraction[.!isnan.(ds.volume_fraction)],
    lnD0_a_2 = ds.lnD0_a_2[.!isnan.(ds.lnD0_a_2)],
    Ea = ds.Ea_kJ_mol[.!isnan.(ds.Ea_kJ_mol)],
    agesteps,
)
@test mdd isa MultipleDomain{Float64, PlanarAr{Float64}}

age, fraction = modelage(mdd, Tsteps)
@test round.(age, sigdigits=5) == [511.05, 513.28, 518.53, 525.12, 532.65, 540.92, 548.7, 556.79, 564.49, 572.86, 581.1, 590.0, 599.88, 611.49, 621.25, 631.23, 641.34, 652.68, 663.1, 674.4, 684.36, 697.53, 711.61, 725.0, 739.63, 758.83, 778.05, 789.42, 799.38, 816.14, 837.41, 855.47, 869.29, 883.52, 895.81, 906.32, 913.78, 919.78, 923.7, 927.11, 929.4, 931.89, 934.39, 937.39, 940.74, 944.76, 949.37, 954.93, 961.16, 967.47, 973.49, 978.48, 982.35, 986.03]
# println(round.(age, sigdigits=5))
@test round.(fraction, sigdigits=5) == [1.8452e-5, 8.6648e-5, 0.00023124, 0.00036146, 0.00061729, 0.00084201, 0.0012625, 0.0016299, 0.0022882, 0.0028575, 0.003833, 0.0046689, 0.0065147, 0.0080275, 0.010368, 0.012338, 0.016248, 0.019464, 0.025456, 0.030366, 0.039093, 0.052982, 0.073361, 0.10128, 0.13616, 0.17558, 0.21951, 0.27052, 0.32554, 0.37653, 0.41564, 0.44896, 0.47947, 0.50148, 0.51861, 0.53272, 0.54495, 0.55588, 0.56589, 0.57518, 0.5839, 0.59462, 0.60757, 0.62287, 0.64053, 0.66043, 0.68235, 0.70857, 0.73591, 0.76723, 0.80596, 0.85255, 0.92177, 1.0]
# println(round.(fraction, sigdigits=5))