module Rotations

using LinearAlgebra

R1(θ) = [1 0 0;
    0 cos(θ) sin(θ);
    0 -sin(θ) cos(θ)]
R2(θ) = [cos(θ) 0 -sin(θ);
    0 1 0;
    sin(θ) 0 cos(θ)]
R3(θ) = [cos(θ) sin(θ) 0;
    -sin(θ) cos(θ) 0;
    0 0 1]
C_321(ψ, θ, ϕ) = R1(ϕ) * R2(θ) * R3(ψ)

tilde(x) = [0 -x[3] x[2]; x[3] 0 -x[1]; -x[2] x[1] 0]

function Cᵦ(β::Vector)
    β₀ = β[1]
    v = β[2:4]
    C = v * transpose(v) + β₀^2 * I - 2 * β₀ * tilde(v) + tilde(v)^2
    return C
end

function EulerParams(ϕ, ê)

    return [
        cos(ϕ / 2),
        ê[1] * sin(ϕ / 2),
        ê[2] * sin(ϕ / 2),
        ê[3] * sin(ϕ / 2)
    ]

end

function CtoEulerParameters(C)

    trC = tr(C)

    bsq = 1 / 4 * [1 + trC,
        1 + 2 * C[1, 1] - trC,
        1 + 2 * C[2, 2] - trC,
        1 + 2 * C[3, 3] - trC]
    i = argmax(bsq)

    if i == 1
        β0 = sqrt(bsq[1])

        β1 = 1 / 4 * (C[2, 3] - C[3, 2]) / β0
        β2 = 1 / 4 * (C[3, 1] - C[1, 3]) / β0
        β3 = 1 / 4 * (C[1, 2] - C[2, 1]) / β0

    elseif i == 2
        β1 = sqrt(bsq[2])

        β0 = 1 / 4 * (C[2, 3] - C[3, 2]) / β1
        β2 = 1 / 4 * (C[1, 2] + C[2, 1]) / β1
        β3 = 1 / 4 * (C[3, 1] + C[1, 3]) / β1

    elseif i == 3
        β2 = sqrt(bsq[3])

        β0 = 1 / 4 * (C[3, 1] - C[1, 3]) / β2
        β1 = 1 / 4 * (C[1, 2] + C[2, 1]) / β2
        β3 = 1 / 4 * (C[2, 3] + C[3, 2]) / β2

    elseif i == 4
        β3 = sqrt(bsq[4])

        β0 = 1 / 4 * (C[1, 2] - C[2, 1]) / β3
        β1 = 1 / 4 * (C[3, 1] + C[1, 3]) / β3
        β2 = 1 / 4 * (C[2, 3] + C[3, 2]) / β3

    end

    return [β0, β1, β2, β3]

end


function Hᵦ(β::Vector)
    # in the body frame:
    # ω = 2 H(β) β̇

    β₀ = β[0+1]
    β₁ = β[1+1]
    β₂ = β[2+1]
    β₃ = β[3+1]

    H = [
        -β₁ +β₀ +β₃ -β₂;
        -β₂ -β₃ +β₀ +β₁;
        -β₃ +β₂ -β₁ +β₀
    ]

    return H
end


end