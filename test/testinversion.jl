
# Test simulated annealing
simannealmodel = (
    σModel = 10., # Model uncertainty [Ma]
    σAnnealing = 100., # Initial uncertainty [Ma]
    λAnnealing = 10 / 10^5, # lambda [1/n]
)

@test simannealsigma(1, 10.; simannealmodel) ≈ 110.44365174144839
@test simannealsigma(10^5, 10.; simannealmodel) ≈ 14.145346247832224

# Test collectto!
buffer = rand(20)
a = rand(3)
b = rand(4)
c = rand(2)
@test Thermochron.collectto!(buffer, a, b, c) == vcat(a,b,c)

# Test maxdiff
for i=1:10
    local a = rand(100)
    @test Thermochron.maxdiff(a) === maximum(diff(a))
    @test Thermochron.mindiff(a) === minimum(diff(a))
    @test Thermochron.maxabsdiff(a) === maximum(abs.(diff(a)))
end

# Test other utilities
@test Thermochron.isdistinct(collect(1:10), 10, 5, 0.5)
@test !Thermochron.isdistinct(collect(1:10), 10, 5, 1.5)
@test Thermochron.pointsininterval(collect(1:10), 10, 0.5, 5.5) == 5
