
# Test simulated annealing
simannealmodel = (
    σmodel = 10., # Model uncertainty [Ma]
    σannealing = 100., # Initial uncertainty [Ma]
    λannealing = 10 / 10^5, # lambda [1/n]
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

# Test diff_ll
Thermochron.diff_ll(1:10, 1, 1) === 0.0
Thermochron.diff_ll(2*(1:10), 1, 1) ≈ -4.5
Thermochron.diff_ll(2*(1:10), 1, 0.5) ≈ -18.0

# Test other utilities
@test Thermochron.isdistinct(collect(1:10), 10, 5, 0.5)
@test !Thermochron.isdistinct(collect(1:10), 10, 5, 1.5)
@test Thermochron.pointsininterval(collect(1:10), 10, 0.5, 5.5) == 5


# Test boundary condition functions
@test Thermochron.reflecting(101, 0, 100) == 99
@test Thermochron.reflecting(101, 1, 100) == 99
@test Thermochron.reflecting(-1, 0, 100) == 1
@test Thermochron.reflecting(-1, 1, 100) == 3
@test Thermochron.reflecting(150, 0, 100) == 50
@test Thermochron.reflecting(200, 0, 100) == 0
@test Thermochron.reflecting(201, 0, 100) == 1
@test Thermochron.reflecting(-101, 0, 100) == 99

@test Thermochron.periodic(101, 0, 100) == 1
@test Thermochron.periodic(101, 1, 100) == 2
@test Thermochron.periodic(-1, 0, 100) == 99
@test Thermochron.periodic(-1, 1, 100) == 98
@test Thermochron.periodic(150, 0, 100) == 50
@test Thermochron.periodic(200, 0, 100) == 0
@test Thermochron.periodic(201, 0, 100) == 1
@test Thermochron.periodic(-101, 0, 100) == 99

@test Thermochron.hard(101, 0, 100) == 100
@test Thermochron.hard(101, 1, 100) == 100
@test Thermochron.hard(-1, 0, 100) == 0
@test Thermochron.hard(-1, 1, 100) == 1
@test Thermochron.hard(150, 0, 100) == 100
@test Thermochron.hard(200, 0, 100) == 100
@test Thermochron.hard(201, 0, 100) == 100
@test Thermochron.hard(-101, 0, 100) == 0
