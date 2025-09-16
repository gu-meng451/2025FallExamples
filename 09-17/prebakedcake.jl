using Pkg
Pkg.activate(".")

# first time:
# Pkg.instantiate()

using DifferentialEquations
using Plots
using LinearAlgebra

# 
A = [0 1;
     -1 0.]
x0 = [1; 0.]
tspan = (0.0, 10.0)

# Matrix Exponential
x(t) = exp(A*t)*x0
plot(t->x(t)[1],tspan[1],tspan[end])

# Numerical solution
function ode!(du, u, p, t)
    A = p
    du[:] = A*u
end
prob = ODEProblem(ode!, x0, tspan, A)
sol = solve(prob)

plot!(sol, idxs=(1), linestyle=:dash)