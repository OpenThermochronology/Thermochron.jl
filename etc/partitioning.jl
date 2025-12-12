cd(@__DIR__)

R = 0.008314472         # [kJ/(K*mol)]
r_boundary = 2e-9       # [m] Grain boundary region thickness

using StatGeochem, Plots
ds = importdataset("baxter 2007 partitioning.csv", importas=:tuple)

## --- Arrhenius (D0 and Ea) fit to Ar partitioning data (i.e., D = D0*exp(-Ea/(R*T)))

Ar_boundary = zeros(size(ds.T_C))
nanadd!(Ar_boundary, ds.Ar_at_Pierce)
nanadd!(Ar_boundary, ds.Ar_at_Crush)
nanadd!(Ar_boundary, ds.Ar_at_150)
nanadd!(Ar_boundary, ds.Ar_at_400)
nanadd!(Ar_boundary, ds.Ar_at_750)
nanadd!(Ar_boundary, ds.Ar_at_1050)
nanadd!(Ar_boundary, ds.Ar_at_1250)
Ar_internal = zeros(size(ds.T_C))
nanadd!(Ar_internal, ds.Ar_at_Fuse)
f_gb = Ar_boundary ./(Ar_boundary .+ Ar_internal)
K_Ar_itm = Ar_internal ./ Ar_boundary .* ds.s_m2_m3 .* r_boundary

t = (ds.duration_h .> 20) .& .!isnan.(K_Ar_itm)
h = plot(framestyle=:box, 
    xlabel="T [C]", 
    ylabel="Kⁱₑ Ar [internal Ar/external Ar]", 
    yscale=:log10, 
    ylims=(10^-6, 10^-4),
)
scatter!(h, ds.T_C[t], K_Ar_itm[t], label="Ar", zcolor=ds.duration_h[t])

using LsqFit
f(x, c) = @. c[1] * exp(-c[2]/(R*(x+273.15)))
c0 = [1., 1.]
fobj = curve_fit(f, ds.T_C[t], K_Ar_itm[t], c0, lower=[1e-9, 0], upper=[1., 200])
K0_Ar, Ea_Ar = fobj.param
@info "K0_Ar = $K0_Ar"
@info "Ea_Ar = $Ea_Ar"
x = range(xlims()..., length=100)
plot!(h, x, @.(K0_Ar*exp(-Ea_Ar/(R*(x+273.15)))), label="$(round(K0_Ar, sigdigits=4)) exp(-$(round(Ea_Ar, sigdigits=4)) / RT)")

## --- Arrhenius (D0 and Ea) fit to He partitioning data (i.e., D = D0*exp(-Ea/(R*T)))

He_boundary = zeros(size(ds.T_C))
nanadd!(He_boundary, ds.He_at_Pierce)
nanadd!(He_boundary, ds.He_at_Crush)
nanadd!(He_boundary, ds.He_at_150)
nanadd!(He_boundary, ds.He_at_400)
He_internal = zeros(size(ds.T_C))
nanadd!(He_internal, ds.He_at_750)
nanadd!(He_internal, ds.He_at_1050)
nanadd!(He_internal, ds.He_at_1250)
nanadd!(He_internal, ds.He_at_Fuse)
f_gb = He_boundary ./(He_boundary .+ He_internal)
K_He_itm = He_internal ./ He_boundary .* ds.s_m2_m3 .* r_boundary

t = (ds.duration_h .> 20) .& .!isnan.(K_He_itm)
h = plot(framestyle=:box, 
    xlabel="T [C]", 
    ylabel="Kⁱₑ He [internal He/external He]", 
    yscale=:log10, 
    ylims=(10^-6, 10^-4),
)
scatter!(ds.T_C[t], K_He_itm[t], label="He", zcolor=ds.duration_h[t])

using LsqFit
f(x, c) = @. c[1] * exp(-c[2]/(R*(x+273.15)))
c0 = [1., 1.]
fobj = curve_fit(f, ds.T_C[t], K_He_itm[t], c0, lower=[1e-9, 0], upper=[1., 200])
K0_He, Ea_He = fobj.param
@info "K0_He = $K0_He"
@info "Ea_He = $Ea_He"
x = range(xlims()..., length=100)
plot!(x, @.(K0_He*exp(-Ea_He/(R*(x+273.15)))), label="$(round(K0_He, sigdigits=4)) exp(-$(round(Ea_He, sigdigits=4)) / RT)")

## --- Resulting equations

function phi_boundary(grainsize_mm::Number)
    r = grainsize_mm/1000
    v = r^3
    vb = (r+r_boundary)^3
    return (vb-v)/vb
end

function fraction_internal_Ar(T, grainsize_mm)
    ϕ = phi_boundary(grainsize_mm)
    Kd = K0_Ar * exp(-Ea_Ar/(R*(T+273.15)))
    Ar_boundary = ϕ
    Ar_internal = (1 - ϕ) * Kd
    return Ar_internal/(Ar_internal+Ar_boundary)
end

function fraction_internal_He(T, grainsize_mm)
    ϕ = phi_boundary(grainsize_mm)
    Kd = K0_He * exp(-Ea_He/(R*(T+273.15)))
    He_boundary = ϕ
    He_internal = (1 - ϕ) * Kd
    return He_internal/(He_internal+He_boundary)
end

## --- End of File