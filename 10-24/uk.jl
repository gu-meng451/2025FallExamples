## UK pendulum example
using Pkg
Pkg.activate(".")

using ForwardDiff
using LinearAlgebra
using Plots

include("ForwardEuler.jl")
using .ForwardEuler

include("BackwardEuler.jl")
using .BackwardEuler

## Setup the system
m = 1.
g = 9.81
l = 0.5

M = diagm(0 => [m, m])
Qex = [0, -m * g]
ϕ(q) = -l^2 + q[1]^2 + q[2]^2

## Built UK equations
A(q) = ForwardDiff.gradient(ϕ, q)'
b(q, qd) = -qd' * ForwardDiff.hessian(ϕ, q) * qd

p = (M, Qex, A, b)
function sysode(z, p, t)
    M, Qex, A, b = p
    n = size(M, 1)

    # z = [q; qd]
    q = z[1:n]
    qd = z[n+1:2*n]

    Ms = sqrt(M)
    a = M \ Qex
    An = A(q)
    bn = b(q, qd)

    qddot = a + Ms \ ((An / Ms) \ (bn - An * a))

    dz = vcat(qd, qddot)

    return dz

end

## check to make sure that the syntax is working
sysode([1, 2, 3, 4], p, 0)

## z = [x, y, xdot, ydot]
tspan = (0, 30.)
θ0 = 30 * pi / 180
x0 = l * sin(θ0)
y0 = -l * cos(θ0)
z0 = [x0, y0, 0, 0]

## some function to compute quantities from the solution
length_error(x,y) = @. (sqrt(x^2 + y^2) - l) / l
energy(x, y, ẋ, ẏ, g, m) = @. 1 / 2 * m * (ẋ^2 + ẏ^2) + m * g * y
relative_energy(X) = energy(X[:, 1], X[:, 2], X[:, 3], X[:, 4], g, m)/energy(z0[1], z0[2], z0[3], z0[4], g, m)


## Using forward Euler:
dt = 0.1
X, time = feuler((x, t) -> sysode(x, p, t), tspan[2], dt, z0)
plt = plot(time, X, xlabel="Time [s]", ylabel="States", title="ForwardEuler")
# plot(time, length_error(X[:,1], X[:,2]), xlabel="Time [s]", ylabel="Normalized error of the length")
# plot(time, relative_energy(X), xlabel="Time [s]", ylabel="Total Energy H/H(0)")


## using BackwardEuler
# dt = 0.01
# X, time = beuler((x, t) -> sysode(x, p, t), tspan[2], dt, z0)
# plot(time, X, line_style=:dash, xlabel="Time [s]", ylabel="States", title="BackwardEuler")
# plot(time, length_error(X[:,1], X[:,2]), xlabel="Time [s]", ylabel="Normalized error of the length")
# plot(time, relative_energy(X), xlabel="Time [s]", ylabel="Total Energy H/H(0)")
