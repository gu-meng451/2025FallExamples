using Pkg
Pkg.activate(".")

using NonlinearSolve

# Residual function
# 0 = [ u[1]^2 - p ;
#       u[2]^2 - p]
#
f(u, p) = u .* u .- p

# initial guess
u0 = [1.0, 1.0]

# parameter
p = 2.0

prob = NonlinearProblem(f, u0, p)
sol = solve(prob, NewtonRaphson())

# sol is a struc
sol.stats

# did it work?
sol.retcode

# get the solution
sol.u

# check the result
f(sol.u,p)