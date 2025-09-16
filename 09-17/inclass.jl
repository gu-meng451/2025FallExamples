using LinearAlgebra

A = [-5.5 -1.1 -1.8 -2.6;
    1 -3.6 0.70 0.9;
    0.5 0.9 -2.8 0.4;
    0 0.4 0.7 -2.1]
B = [-0.2; 0.3; -0.2; 0.3]
C = [3 5 7 5]

# compute the eigenvalue and eigenvectors
F = eigen(A)

# write a function to replace values near zero with zero
function clip(x; tol=100 * eps(eltype(x)))
    idx = abs.(x) .< tol
    z = copy(x)
    z[idx] .= zero(eltype(x))
    return z
end

P = F.vectors

Abar = P \ A * P |> clip;
display(Abar)

Bbar = P \ B |> clip;
display(Bbar)

Cbar = C * P |> clip;
display(Cbar)
