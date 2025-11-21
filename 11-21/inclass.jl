using Pkg
Pkg.activate(".")

using LinearAlgebra

m = 1.
grav = 9.81
l = 1.
M = diagm(0=>[m, m])
f(q,v) = [0; -m*grav]
g(q) = -l^2 + q[1]^2 + q[2]^2
G(q) = [ q[1] q[2] ]

n = 2
nc = 1

function F(t,x,z)

    q = x[1:n]
    v = x[n+1:2n]

    λ = z[1:nc]
    μ = z[nc+1:2nc]

    return [v - G(q)'*μ;
    M\( f(q,v) - G(q)'*λ)]

end

F(0, [1,2,3,4], [5,6] )


function H(t,x,z)

    q = x[1:n]
    v = x[n+1:2n]

    λ = z[1:nc]
    μ = z[nc+1:2nc]

    return [ g(q);
            G(q)*v]
end

H(0,[1,2,3,4],[5,6])