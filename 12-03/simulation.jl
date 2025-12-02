using Pkg
Pkg.activate(".")
using LinearAlgebra
using ForwardDiff
using GLMakie

include("DAE.jl")
include("Rotations.jl")

# system parameters
m1 = 1.
grav = 9.81
# location of the pivot in the body frame {B}
rp = [-1., 0, 0]

# location of the pivot in the inertial frame {N}
N_rp = [0., 0, 0]

# The mass moment of inertia in the body frame
Ig = [1. 0 0;
    0 2. 0;
    0 0 3.];

# The order of the coordinates for each body are
# [x,y,z,β{4}] ∈ R^7

function unpack(q)
    # for a single body
    X = q[1:3]
    β = q[4:7]
    return X, β
end

# build the mass matrix
function M(q)

    Mass = zeros(eltype(q),7, 7)

    Mass[1, 1] = m1
    Mass[2, 2] = m1
    Mass[3, 3] = m1

    _, β = unpack(q)
    Hβ = Rotations.Hᵦ(β)
    Mass[4:7, 4:7] = 4 * Hβ' * Ig * Hβ

    return Mass

end

# build f(q,v)
function f(q, v)

    Q = zeros(7)

    Q[1:3] = [0, 0, -m1 * grav]

    _, β = unpack(q)
    _, β̇ = unpack(v)
    Hβ = Rotations.Hᵦ(β)
    Ḣβ = Rotations.Hᵦ(β̇)
    Q[4:7] = -4 * Ḣβ' * Ig * Hβ * β̇

    return Q

end

# build the holonomic constraints
function g(q)
    N_rg, β = unpack(q)

    # make output container
    # add eltype to ensure compatibility with ForwardDiff
    constraint = zeros(eltype(q), 4) 

    # fix the pin location
    C_NB = Rotations.Cᵦ(β)
    constraint[1:3] = N_rp - N_rg - C_NB' * rp

    # unit quaternion
    constraint[4] = β'β - 1

    return constraint
end
G(q) = ForwardDiff.jacobian(g, q)

# number of coordinates (6 per 3d rigid body)
n = 7

# number of constraints
nc = 4

p = M, f, g, G, n, nc

## Build GGL system functions
function F(t, x, z, p)

    M, f, g, G, n, nc = p

    q = x[1:n]
    v = x[n+1:2n]
    λ = z[1:nc]
    μ = z[nc+1:2nc]

    return [v - transpose(G(q)) * μ;
        M(q) \ (f(q, v) - transpose(G(q)) * λ)]
end

function H(t, x, z, p)
    M, f, g, G, n, nc = p
    q = x[1:n]
    v = x[n+1:2n]

    return [g(q); G(q) * v]

end

## make consistent initial conditions
C_0 = Rotations.C_321(0.1, 0.1, 0.1)
X_0 = C_0'*(-rp)
β_0 = Rotations.CtoEulerParameters(C_0)

# velocity initial conditions
xd_0 = 0.
yd_0 = 0.
zd_0 = 0.
β0d_0 = 0.
β1d_0 = 0.
β2d_0 = 0.
β3d_0 = 0.

# pack the initial condition
q0 = vcat(X_0, β_0)
v0 = [xd_0, yd_0, zd_0, β0d_0, β1d_0, β2d_0, β3d_0]


# verify they work
g(q0) |> norm
G(q0)*v0 |> norm

## integrate the system
h = 0.01
nsteps = 1000
X, Z, t = DAE.irk_integrator(F, H, g, p, vcat(q0,v0), h, nsteps, irk=DAE.radauIIA3)

## Plot using GLMakie
fig = Figure(size=(800, 600))
ax = Axis3(fig[1, 1];
    aspect=DataAspect(),
    xlabel="x (m)", ylabel="y (m)", zlabel="z (m)")
lines!(ax, X[1, :], X[2, :], X[3, :], label="COM")

display(fig)

