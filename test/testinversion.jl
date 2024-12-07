
# Test simulated annealing
@test Thermochron.simannealsigma(1, 10., 100., 10 / 10^5) ≈ 100.48880634173474
@test Thermochron.simannealsigma(10^5, 10., 100., 10 / 10^5) ≈ 10.000001030576758

# Test collectto!
buffer = rand(20)
a = rand(3)
b = rand(4)
c = rand(2)
@test Thermochron.collectto!(buffer, a, b, c) == vcat(a,b,c)

# Test pointsininterval
@test Thermochron.pointsininterval(collect(1:10), 10, 0.5, 5.5) == 5

# Test diff_ll
@test Thermochron.diff_ll(1:10, 1, 1) == 0.0
@test Thermochron.diff_ll(2*(1:10), 1, 1) ≈ -4.5
@test Thermochron.diff_ll(2*(1:10), 1, 0.5) ≈ -18.0
x = rand(100)
d = diff(x)
@test Thermochron.diff_ll(x, 0, 1) ≈ normpdf_ll(0, 1, d[d .> 0])

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

## ---Test generation of Chronometer objects
