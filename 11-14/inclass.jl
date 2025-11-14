## Inclass demo of using Euler Angles
using Pkg
Pkg.activate(".")

using LinearAlgebra
using GLMakie

l0 = 0.5
l1 = 1.0
l2 = 0.6

R1(θ) = [1 0 0;
    0 cos(θ) sin(θ);
    0 -sin(θ) cos(θ)]

P01 = [0 0 1;
    1 0 0;
    0 1 0]

N_rA = [0; 0; l0]

N_rB(θ1) = N_rA + P01' * R1(θ1)' * [0; 0; l1]

P12 = [0 0 -1;
    1 0 0;
    0 -1 0]
N_rC(θ1, θ2) = N_rB(θ1) + P01' * R1(θ1)' * P12' * R1(θ2)' * [0; l2; 0]

fig = Figure(size=(800, 600))
ax = Axis3(fig[1, 1], aspect=:equal)

θ1 = 0
θ2 = -1
lines!(ax, [0, N_rA[1]],
    [0, N_rA[2]],
    [0, N_rA[3]],
    linewidth=4)
lines!(ax, [N_rA[1], N_rB(θ1)[1]],
    [N_rA[2], N_rB(θ1)[2]],
    [N_rA[3], N_rB(θ1)[3]],
    linewidth=4)
lines!(ax, [N_rB(θ1)[1], N_rC(θ1, θ2)[1]],
    [N_rB(θ1)[2], N_rC(θ1, θ2)[2]],
    [N_rB(θ1)[3], N_rC(θ1, θ2)[3]],
    linewidth=4)


display(fig)