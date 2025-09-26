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

