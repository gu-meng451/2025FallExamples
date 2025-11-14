using Pkg
Pkg.activate(".")

using LinearAlgebra
using GLMakie
using Printf

R1(θ) = [1 0 0;
        0 cos(θ) sin(θ);
        0 -sin(θ) cos(θ)]

l0 = 0.5
l1 = 1.0
l2 = 0.6

N_rA = [0; 0; l0]
P01 = [0 0 1; 
       1 0 0;
       0 1 0]
C_NB(θ1) = R1(θ1)*P01
N_rBA(θ1) = C_NB(θ1)'*[0; 0; l1]
P12 = [ 0 0 -1; 
      1 0 0;
      0 -1 0]

C_NC(θ1, θ2) = R1(θ2)*P12*C_NB(θ1)
N_rCB(θ1,θ2) = C_NC(θ1,θ2)'*[0; l2; 0]

N_rB(θ1) = N_rA + N_rBA(θ1)
N_rC(θ1,θ2) = N_rB(θ1) + N_rCB(θ1,θ2)

# Create a Figure with an Axis3 and two sliders (θ1, θ2)
fig = Figure(size = (1000, 1000))

# Create sliders at the top first
slider1 = Slider(fig[1, 1], range = -2π:0.01:2π, startvalue = 0.0, width = 300, height = 20)
Label(fig[1, 0], "θ1", padding = (4, 4, 4, 4))

slider2 = Slider(fig[2, 1], range = -2π:0.01:2π, startvalue = 0.0, width = 300, height = 20)
Label(fig[2, 0], "θ2", padding = (4, 4, 4, 4))

# Get the observable values from the sliders
θ1 = slider1.value
θ2 = slider2.value

# Create observable text labels with fixed-width formatting to prevent jumping
θ1_deg_text = @lift @sprintf "θ1: %7.1f°" rad2deg($θ1)
θ2_deg_text = @lift @sprintf "θ2: %7.1f°" rad2deg($θ2)

Label(fig[1, 2], θ1_deg_text, padding = (4, 4, 4, 4))
Label(fig[2, 2], θ2_deg_text, padding = (4, 4, 4, 4))

# Set fixed width for column 2 so text doesn't cause plot jumping
colsize!(fig.layout, 2, Fixed(120))

# Plot spans remaining rows, taking up most of the space
ax = Axis3(fig[3, 1:2], xlabel="X", ylabel="Y", zlabel="Z", title="Robot Arm Kinematics")

# Set axis limits
xlims!(ax, -1.5, 1.5)
ylims!(ax, -1.5, 1.5)
zlims!(ax, -0.2, 1.1)

# Make rows 1-2 small (sliders) and row 3 flexible (plot)
rowsize!(fig.layout, 1, Fixed(35))
rowsize!(fig.layout, 2, Fixed(35))

# Reduce padding and margins to minimum
colsize!(fig.layout, 1, Relative(1.0))
colgap!(fig.layout, 0)
rowgap!(fig.layout, 0)

# Use @lift to create reactive observables for joint positions
# These automatically update when θ1 or θ2 change
x_AB = @lift [N_rA[1], N_rB($θ1)[1]]
y_AB = @lift [N_rA[2], N_rB($θ1)[2]]
z_AB = @lift [N_rA[3], N_rB($θ1)[3]]

x_BC = @lift [N_rB($θ1)[1], N_rC($θ1, $θ2)[1]]
y_BC = @lift [N_rB($θ1)[2], N_rC($θ1, $θ2)[2]]
z_BC = @lift [N_rB($θ1)[3], N_rC($θ1, $θ2)[3]]

x_joints = @lift [N_rA[1], N_rB($θ1)[1], N_rC($θ1, $θ2)[1]]
y_joints = @lift [N_rA[2], N_rB($θ1)[2], N_rC($θ1, $θ2)[2]]
z_joints = @lift [N_rA[3], N_rB($θ1)[3], N_rC($θ1, $θ2)[3]]

# Plot the robot arm segments with reactive observables
lines!(ax, [0, N_rA[1]], [0, N_rA[2]], [0, N_rA[3]], color="#1f77b4", linewidth=6)
lines!(ax, x_AB, y_AB, z_AB, color="#ff7f0e", linewidth=6)
lines!(ax, x_BC, y_BC, z_BC, color="#2ca02c", linewidth=6)
scatter!(ax, x_joints, y_joints, z_joints, color="#d62728", markersize=20)

display(fig)
