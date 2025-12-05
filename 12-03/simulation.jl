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


# integration parameters
h = 0.1
nsteps = 1000

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

    Q = zeros(eltype(q),7)

    Q[1:3] = [0, 0, -m1 * grav]

    _, β = unpack(q)
    _, β̇ = unpack(v)
    Hβ = Rotations.Hᵦ(β)
    Ḣβ = Rotations.Hᵦ(β̇)
    ω = 2*Hβ*β̇
    f_gyro = cross( ω, Ig*ω )
    Q[4:7] = -4 * Ḣβ' * Ig * Hβ * β̇ - 2 * Hβ' * f_gyro

    return Q

end

# build the holonomic constraints
function g(q)
    N_rg, β = unpack(q)

    # make output container
    # add eltype to ensure compatibility with ForwardDiff
    constraint = zeros(eltype(q), 4) 

    # fix the pin location
    C_BN = Rotations.Cᵦ(β)
    constraint[1:3] = N_rp - N_rg - C_BN' * rp

    # unit quaternion constraint
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

    _, f, _, G, n, nc = p

    q = x[1:n]
    v = x[n+1:2n]
    λ = z[1:nc]
    μ = z[nc+1:2nc]

    # since M can be singular it is handled in the solver directly
    return [v - transpose(G(q)) * μ;
        f(q, v) - transpose(G(q)) * λ]
end

function H(t, x, z, p)
    M, f, g, G, n, nc = p
    q = x[1:n]
    v = x[n+1:2n]

    return [g(q);
        G(q) * v]
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

# verify the IC's work
g(q0) |> norm
G(q0)*v0 |> norm

## integrate the system
X, Z, t = DAE.irk_integrator(F, H, p, vcat(q0,v0), h, nsteps, irk=DAE.radauIIA3)

## Static plot using GLMakie
fig = Figure(size=(800, 600))
ax = Axis3(fig[1, 1];
    aspect=:data,
    xlabel="x (m)", ylabel="y (m)", zlabel="z (m)")
lines!(ax, X[1, :], X[2, :], X[3, :], label="COM")

display(fig)

## Animation using GLMakie
fig_anim = Figure(size=(800, 600))
ax_anim = Axis3(fig_anim[1, 1];
    aspect=:data,
    xlabel="x (m)", ylabel="y (m)", zlabel="z (m)",
    limits=((-2, 2), (-2, 2), (-2, 2)))

# Unpack solution to get position and orientation over time
positions = [X[1:3, i] for i in 1:size(X, 2)]
quaternions = [X[4:7, i] for i in 1:size(X, 2)]

# Precompute rotation matrices from quaternions (for original sample times)
C_matrices = [Rotations.Cᵦ(β) for β in quaternions]

# Resample solution to a fixed playback frame rate (30 fps)
fps = 30
simulation_time = t[end]
num_frames = Int(round(simulation_time * fps))
resample_times = collect(range(t[1], stop=t[end], length=num_frames))

function resample_rows(X_rows::AbstractMatrix, t_old::AbstractVector, t_new::AbstractVector)
    nrows = size(X_rows, 1)
    Y = similar(X_rows, nrows, length(t_new))
    for i in 1:nrows
        yi = X_rows[i, :]
        for (j, τ) in enumerate(t_new)
            if τ <= t_old[1]
                Y[i, j] = yi[1]
            elseif τ >= t_old[end]
                Y[i, j] = yi[end]
            else
                k = searchsortedlast(t_old, τ)
                t0 = t_old[k]; t1 = t_old[k+1]
                α = (τ - t0) / (t1 - t0)
                Y[i, j] = (1 - α) * yi[k] + α * yi[k+1]
            end
        end
    end
    return Y
end

# Resample positions and quaternions
pos_res = resample_rows(X[1:3, :], t, resample_times)
quat_res = resample_rows(X[4:7, :], t, resample_times)

# Normalize quaternions and compute rotation matrices at resampled times
quat_normed = similar(quat_res)
C_matrices_res = Vector{Matrix{eltype(X)}}(undef, num_frames)
for i in 1:num_frames
    q = quat_res[:, i]
    q = q / norm(q)
    quat_normed[:, i] = q
    C_matrices_res[i] = Rotations.Cᵦ(q)
end


# Convert resampled positions to Point3f for plotting
traj_points = [Point3f(pos_res[:, i]...) for i in 1:num_frames]

# Observable for animation (indexed over resampled frames)
frame = Observable(1)

# Draw the full trajectory as a solid blue line (previous look)
# Translucent trajectory (blue with 50% opacity)
lines!(ax_anim, traj_points, label="COM trajectory", color=(0.0, 0.45, 0.8, 0.3), linewidth=2)

# Animated body line
body_line = @lift begin
    idx = $frame
    pos = Point3f(pos_res[:, idx]...)
    C = C_matrices_res[idx]
    p_pivot = Point3f((pos .+ C' * N_rp)...) 
    p_com = Point3f((pos .+ C' * rp)...) 
    [p_pivot, p_com]
end

# Body frame axes
body_axes = @lift begin
    idx = $frame
    pos = pos_res[:, idx]
    C = C_matrices_res[idx]
    p_pivot = pos + C' * N_rp
    scale = 0.3
    x_axis = [Point3f(p_pivot...), Point3f((p_pivot + scale * C' * [1, 0, 0])...)]
    y_axis = [Point3f(p_pivot...), Point3f((p_pivot + scale * C' * [0, 1, 0])...)]
    z_axis = [Point3f(p_pivot...), Point3f((p_pivot + scale * C' * [0, 0, 1])...)]
    (x_axis, y_axis, z_axis)
end

lines!(ax_anim, body_line, color=:darkorange, linewidth=4, label="Body")
lines!(ax_anim, @lift(($body_axes)[1]), color=:red, linewidth=2, label="X-axis")
lines!(ax_anim, @lift(($body_axes)[2]), color=:forestgreen, linewidth=2, label="Y-axis")
lines!(ax_anim, @lift(($body_axes)[3]), color=:royalblue, linewidth=2, label="Z-axis")

# Plot markers for pivot and COM
body_marker = @lift begin
    idx = $frame
    pos = pos_res[:, idx]
    C = C_matrices_res[idx]
    p_pivot = Point3f((pos + C' * N_rp)...) 
    p_com = Point3f((pos + C' * rp)...) 
    [p_pivot, p_com]
end
scatter!(ax_anim, body_marker, color=[:black, :red], markersize=8)

# Time display (use resampled times)
time_text = @lift("t = $(round(resample_times[$frame], digits=2)) s")
text!(ax_anim, -1.8, 1.8, 1.8, text=time_text, fontsize=16)

axislegend(ax_anim)

# Record animation at fixed 30 fps using resampled frames
println("Recording $(num_frames) frames at $(fps) fps (duration $(round(simulation_time, digits=2)) s) ...")
record(fig_anim, "body_animation.mp4", 1:num_frames; framerate=fps) do i
    frame[] = i
end

println("Animation saved as body_animation.mp4")
println("Simulation time: $(round(simulation_time, digits=2)) s")
println("Frames: $(num_frames), fps: $(fps)")

