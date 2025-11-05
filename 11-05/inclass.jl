using Pkg
Pkg.activate(".")

using LinearAlgebra
using NonlinearSolve
using Plots

## demonstration
function odesys(u, p, t)
    A = p
    return A * u
end

A = [0 1; -1 -0.1]
p = A
u0 = [1, 1.]
tf = 10.
u_true(t) = exp(A * t) * u0

## Adaptive time-stepping
function rk_adapt_step(f, xn, tn, p, h)

    a21 = 1
    c2 = 1
    b1 = 1 / 2
    b2 = 1 / 2

    k1 = f(xn, p, tn)
    k2 = f(xn + h * a21 * k1, p, tn + c2 * h)

    # 2nd order step
    xn1 = xn + h * (b1 * k1 + b2 * k2)

    # 1st order step
    b1s = 1
    b2s = 0
    z = xn + h * (b1s * k1 + b2s * k2)

    return xn1, z

end

function integ_adapt_h(f, x0, p, h0, tf; tol=1e-4)

    n = length(x0)
    # X = zeros(nsteps + 1, n)
    # t = 0:h:h*nsteps

    X = x0'
    t = [0.]

    i = 1
    h = h0
    while t[end] < tf

        xn = X[i, :]
        tn = t[i]
        xn1, z = rk_adapt_step(f, xn, tn, p, h)

        # TODO replace with better method
        err = norm(xn1 - z)

        if err < tol
            X = vcat(X, xn1')
            append!(t, tn + h)
            i += 1
            h = 1.1 * h

        else
            h = h / 2
        end

    end


    return X, t

end


##
h0 = 0.1
X, t = integ_adapt_h(odesys, u0, p, h0, tf, tol=1e-1)

## Plots
plot(t, X, xlabel="Time t", ylabel="State")
plot!(t -> u_true(t)[1], 0, tf, label="true y1", line=:dash)
plot!(t -> u_true(t)[2], 0, tf, label="true y2", line=:dash)
