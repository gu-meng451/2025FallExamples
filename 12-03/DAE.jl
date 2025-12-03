module DAE

using LinearAlgebra
using NonlinearSolve

# ==============================================================================
# 1. MODIFIED SYSTEM FUNCTION INSTRUCTIONS
# ==============================================================================
# The function F must now return the RHS forces, NOT the acceleration.
# DO NOT compute M \ ... inside this function.
#
# Expected return form:
# [ v - G' * μ ;        <-- Kinematic RHS (usually just v if constraints handled via μ)
#   f - G' * λ ]        <-- Dynamic RHS (Forces)
#
# The solver handles the Mass Matrix multiplication internally.
# ==============================================================================


# ==============================================================================
# 2. UPDATED IRK STRUCT (Added invA)
# ==============================================================================
struct irk
    s::Int
    A::Matrix
    invA::Matrix  # <--- NEW: Inverse of A is required for implicit formulation
    b::Vector
    c::Vector
end

function radauIIA2_constructor()
    A = [5/12 -1/12; 3/4 1/4]
    return irk(2, A, inv(A), [3 / 4, 1 / 4], [1 / 3, 1])
end
const radauIIA2 = radauIIA2_constructor()

function radauIIA3_constructor()
    A = [(11/45-7*sqrt(6)/360) (37/225-169*sqrt(6)/1800) (-2/225+sqrt(6)/75);
        (37/225+169*sqrt(6)/1800) (11/45+7*sqrt(6)/360) (-2/225-sqrt(6)/75);
        (4/9-sqrt(6)/36) (4/9+sqrt(6)/36) (1/9)]
    return irk(3, A, inv(A),
        [4 / 9 - sqrt(6) / 36, 4 / 9 + sqrt(6) / 36, 1 / 9],
        [2 / 5 - sqrt(6) / 10, 2 / 5 + sqrt(6) / 10, 1])
end
const radauIIA3 = radauIIA3_constructor()

# ==============================================================================
# 3. REWRITTEN RESIDUAL FUNCTION (Handles Singular M)
# ==============================================================================
function residual_irk(U, F, H, p, xn, tn, h, irk)
    s = irk.s
    # invA is used to recover derivatives from stages
    invA = irk.invA
    M_func, _, _, _, n, nc = p

    # View wrappers
    X = @view U[1:2n, 1:s]     # Differential variables (q, v)
    Z = @view U[2n+1:end, 1:s] # Algebraic variables (λ, μ)

    # -----------------------------------------------------------------------
    # Step A: Recover the Derivatives (K) from the Stage Values (X)
    # Relation: X_stage = x_n + h * sum(A_ij * K_j)
    # Inverted: K_stage = (1/h) * sum(invA_ij * (X_j - x_n))
    # -----------------------------------------------------------------------

    # Calculate (X - xn) for all stages. 
    # X is 2n x s. xn is vector 2n.
    X_diff = X .- xn

    # Compute K (derivatives) via matrix multiplication
    # (2n x s) * (s x s) -> (2n x s)
    K = (X_diff * transpose(invA)) ./ h

    # -----------------------------------------------------------------------
    # Step B: Build Residuals enforcing M*K = Forces
    # -----------------------------------------------------------------------

    # We will accumulate the residuals for the differential part here
    T = eltype(U)
    R_diff = Matrix{T}(undef, 2n, s)

    for j in 1:s
        # Current stage values
        xj = X[:, j] # [q; v]
        zj = Z[:, j] # [λ; μ]
        tj = tn + irk.c[j] * h

        q = xj[1:n]

        # 1. Get RHS Forces from user function F
        # F returns: [ v - G'μ ; f - G'λ ]
        rhs = F(tj, xj, zj, p)

        # 2. Split derivative K into q_dot and v_dot
        k_q = K[1:n, j]      # q_dot
        k_v = K[n+1:2n, j]   # v_dot

        # 3. Enforce Equations
        # Eq 1: Kinematic: q_dot = (v - G'μ) 
        #       => residual = k_q - rhs[1:n]
        R_diff[1:n, j] = k_q - rhs[1:n]

        # Eq 2: Dynamic: M(q) * v_dot = (f - G'λ)
        #       => residual = M(q)*k_v - rhs[n+1:end]

        # Handle singular mass matrix here:
        mass_mat = M_func(q) # Evaluate M(q)
        R_diff[n+1:end, j] = mass_mat * k_v - rhs[n+1:end]
    end

    # -----------------------------------------------------------------------
    # Step C: Algebraic Constraints (H)
    # -----------------------------------------------------------------------
    R_alg = [H(tn + irk.c[j] * h, X[:, j], Z[:, j], p) for j in 1:s]

    # Flatten and return
    return vcat(reduce(hcat, eachcol(R_diff)), reduce(hcat, R_alg))
end

function irk_step(F, H, p, xn, zn, tn, h, irk; verbose=false)
    _, _, _, _, n, nc = p

    R(U, p) = residual_irk(U, F, H, p, xn, tn, h, irk)

    # Initial Guess:
    # Note: With singular mass matrices, using the previous step [xn; zn] 
    # as a guess is sometimes risky if the constraints drift, but it is 
    # the standard simple starting point.
    U0 = zeros(Float64, 2n + 2nc, irk.s)
    U0 .= [xn; zn]

    prob = NonlinearProblem(R, U0, p)
    sol = solve(prob, NewtonRaphson())

    if sol.retcode != ReturnCode.Success
        println("Warning: Solver retcode: $(sol.retcode)")
        # In a robust code, you would reduce h and retry here
    end

    X = @view sol.u[1:2n, :]
    Z = @view sol.u[2n+1:end, :]

    # For Radau IIA (stiffly accurate), the last stage is the solution
    xn1 = X[:, end]
    zn1 = Z[:, end]

    if verbose
        return xn1, zn1, sol
    else
        return xn1, zn1
    end
end

function irk_integrator(F, H, p, x0, h, nsteps; irk=radauIIA2)
    _, _, _, _, n, nc = p
    X = zeros(2 * n, nsteps + 1)
    Z = zeros(2 * nc, nsteps + 1)
    t = 0:h:h*nsteps

    X[:, 1] = x0
    # Initialize Z with zeros or a consistent solve if available
    # Ideally, one should solve for consistent λ, μ at t=0

    for i = 1:nsteps
        X[:, i+1], Z[:, i+1] = irk_step(F, H, p, X[:, i], Z[:, i], t[i], h, irk)
    end

    return X, Z, t
end

end