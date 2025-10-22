## Continuation Example
```
Based on the problem from Fig 2 in 
https://doi.org/10.1016/j.ijsolstr.2020.11.037
```

##
using Pkg
Pkg.activate(".")
using LinearAlgebra
using Plots

# define the system residual function
function F(u, λ, p)

    # unpack the parameters
    α, κ = p
    uₐ, uₚ = u
    fₐ = λ

    R = [(-2 * α * uₚ + 2 * α + sqrt(α^2 + uₚ^2 - 2 * uₚ) * (-κ * uₐ + κ * uₚ + 2 * uₚ - 2)) / sqrt(α^2 + uₚ^2 - 2 * uₚ),
        -fₐ + κ * (uₐ - uₚ)]

    return R

end

# build the jacobians
function F_u(u, λ, p)
    # unpack the parameters
    α, κ = p
    uₐ, uₚ = u
    fₐ = λ

    J = [-κ ((uₚ-1)*(α^2+uₚ^2-2*uₚ)*(2*α*uₚ-2*α+sqrt(α^2 + uₚ^2 - 2 * uₚ)*(κ*uₐ-κ*uₚ-2*uₚ+2))+((-2*α+(κ+2)*sqrt(α^2 + uₚ^2 - 2 * uₚ))*sqrt(α^2 + uₚ^2 - 2 * uₚ)-(uₚ-1)*(κ*uₐ-κ*uₚ-2*uₚ+2))*(α^2+uₚ^2-2*uₚ)^(3//2))/(α^2+uₚ^2-2*uₚ)^(5//2);
        κ -κ]

    return J
end

function F_λ(u, λ, p)
    return [0, -1]
end


# build the tangent
function compute_tangent(J_u, J_λ, u̇_prev, λ̇_prev)
    # Solve the (m+1) × (m+1) system:
    # [ J_u    J_λ  ] [ u̇ ]   [ 0 ]
    # [ u̇ᵀ_p  λ̇_p ] [ λ̇ ] = [ 1 ]
    #
    # where (u̇_p, λ̇_p) is the previous tangent

    m = length(J_λ)

    A = [J_u J_λ;
        u̇_prev' λ̇_prev]

    b = [zeros(m); 1]

    tangent = A \ b

    u̇ = tangent[1:m]
    λ̇ = tangent[m+1]

    # Normalize
    norm_factor = sqrt(dot(u̇, u̇) + λ̇^2)
    u̇ = u̇ / norm_factor
    λ̇ = λ̇ / norm_factor

    return u̇, λ̇
end

function compute_initial_tangent(J_u, J_λ, direction=1)
    # Set λ̇ = direction (typically +1 to go forward)
    # Solve: J_u * u̇ = -J_λ
    u̇ = -J_u \ J_λ

    # Normalize
    norm_factor = sqrt(dot(u̇, u̇) + direction^2)
    u̇ = u̇ / norm_factor
    λ̇ = direction / norm_factor

    return u̇, λ̇
end

function pseudo_arclength_continuation(p, u0, λ0, target_u_a, Δs_initial; 
    direction=1, max_iterations=1000, tolerance=1e-5)
    
    # Initialize
    u = u0
    λ = λ0
    Δs = Δs_initial
    solutions = [(u, λ)]

    # ===== COMPUTE INITIAL TANGENT =====
    J_u = F_u(u,λ,p)
    J_λ = F_λ(u,λ,p)

    (u̇, λ̇) = compute_initial_tangent(J_u, J_λ, direction)

    step = 1

    while u[1] < target_u_a
        
        # PREDICTOR
        u_pred = u + Δs * u̇
        λ_pred = λ + Δs * λ̇

        # CORRECTOR (Newton iteration)
        u_new = u_pred
        λ_new = λ_pred

        for iter = 1:max_iterations
            J_u = F_u(u_new, λ_new, p)
            J_λ = F_λ(u_new, λ_new, p)
            
            r_f = F(u_new, λ_new, p)
            r_N = dot(u_new - u, u̇) + (λ_new - λ) * λ̇ - Δs
            
            if norm([r_f; r_N]) < tolerance
                break
            end
            
            J_aug = [J_u      J_λ;
                     u̇'      λ̇]
            
            r_aug = [r_f; r_N]
            
            Δx = -J_aug \ r_aug
            
            u_new = u_new + Δx[1:end-1]
            λ_new = λ_new + Δx[end]

        end

        # COMPUTE NEW TANGENT
        J_u = F_u(u_new, λ_new, p)
        J_λ = F_λ(u_new, λ_new, p)

        if step == 1
            # Still use Method 1 for second tangent if you prefer
            (u̇_new, λ̇_new) = compute_initial_tangent(J_u, J_λ, direction)
        else
            # Use Method 2 (bordered system with previous tangent)
            (u̇_new, λ̇_new) = compute_tangent(J_u, J_λ, u̇, λ̇)
        end

        # Ensure consistent orientation
        if dot(u̇_new, u̇) + λ̇_new * λ̇ < 0
            u̇_new = -u̇_new
            λ̇_new = -λ̇_new
        end

        # UPDATE
        u = u_new
        λ = λ_new
        u̇ = u̇_new
        λ̇ = λ̇_new

        push!(solutions, (u, λ))
        step += 1

    end

    return solutions

end

## Setup and solve the particular problem
p = 5 / 4, 1 / 2

sol = pseudo_arclength_continuation(p,[0,0.],0.,2,0.05)

# extract the solution for plotting
ua = [ x[1][1] for x in sol]
fa = [ x[2] for x in sol]
plot(ua,fa, xlabel="u_a", ylabel="f_a")
