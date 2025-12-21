using Thermochron, Measurements
ds = importdataset("uthsmhe.csv", importas=:Tuple)

raw_He_age = fill(NaN, size(ds.U_ng))
raw_He_age_sigma = fill(NaN, size(ds.U_ng))
for i in eachindex(raw_He_age)
    U = ds.U_ng[i] ± ds.U_sigma_ng[i]
    Th = ds.Th_ng[i] ± ds.Th_sigma_ng[i]
    Sm = haskey(ds, :Sm_ng) ? ds.Sm_ng[i] ± ds.Sm_sigma_ng[i] : 0.0
    He = ds.He_pmol[i] ± ds.He_sigma_pmol[i]
    a = Thermochron.newton_he_age(He, U*137.818/138.818*1000/238, U*1/138.818*1000/235, Th*1000/232, Sm*0.15*1000/147)
    raw_He_age[i] = a.val
    raw_He_age_sigma[i] = a.err
end

ds = (;ds..., raw_He_age=raw_He_age, raw_He_age_sigma=raw_He_age_sigma)
exportdataset(ds, "uthsmhe_calc.csv")