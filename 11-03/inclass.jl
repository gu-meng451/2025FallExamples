using Pkg
Pkg.activate(".")

using Plots

## ## demonstration
function odesys(u, p, t)
    A = p
    return A * u
end

A = [0 1; -1 -0.1]
u0 = [1, 1.]
tf = 10.
u_true(t) = exp(A * t) * u0
plot(t -> u_true(t)[1], 0, tf, label="true y1", line=:dash)
plot!(t -> u_true(t)[2], 0, tf, label="true y2", line=:dash)


## Explicit RK step


## constant time-step integration


## Adaptive time-stepping


## Plots

