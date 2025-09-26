using Pkg
Pkg.activate(".")
using Plots
using DifferentialEquations
using LinearAlgebra
using Polynomials: fromroots


function system!(dx, x, p, t)

    # unpack parameters
    A, B, u = p

    # make an in-place change to du.  Do no write du = [1,2,3] as that will not work.
    dx[:] = A*x + B*u(x,t)
    
end

M = diagm(0=>[10., 1])
C = [10. 0; 0 0]
K = [100 -10; -10 10.]

A = [zeros(2,2) I;
     -M\K -M\C]
B = [ zeros(2,1); M\[0, 1]]
C = [1 0 0 0]

u(x,t) = 0.

# solve the system
p = A,B,u
x0 = [0, 1., 0., 0.]
tspan = (0., 10.)
prob = ODEProblem(system!, x0, tspan, p)
sol = solve(prob)

# Now let's control the system
poles = [-3+3im, -3-3im, -10, -15]
desired_char = fromroots(poles)
d = desired_char.coeffs

#
Cmx = [ B A*B A^2*B A^3*B ]
det(Cmx)

# build Abar
λ = eigvals(A)
a = fromroots(λ).coeffs

Abar = diagm(1=>ones(3))
Abar[end,:] = -a[1:end-1]
Abar
Bbar = zeros(4,1)
Bbar[end] = 1

Cmz = [Bbar Abar*Bbar Abar^2*Bbar Abar^3*Bbar]

P = Cmx/Cmz

Kz = (d[1:end-1] - a[1:end-1])'
Kx = Kz/P

## 
u(x,t) = -Kx*x
p = A, B, u
x0 = [0, 1., 0., 0.]
tspan = (0., 10.)
prob = ODEProblem(system!, x0, tspan, p)
sol = solve(prob)
plot(sol)

eigvals(A-B*Kx)