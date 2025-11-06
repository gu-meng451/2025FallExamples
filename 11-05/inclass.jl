using Pkg
Pkg.activate(".")

using LinearAlgebra
# using NonlinearSolve
using Plots
using Printf

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

function stepsizefactor(xn1, z, xn; order=1, fac=0.9, facmin=0.1, facmax=4., Atol=1e-3, Rtol=1e-3)
    # order: is based on the lower-order method

    errn = norm(xn1 - z)
    toln = Atol + Rtol * max(norm(xn1), norm(xn))
    En = errn / toln
    return En, min(facmax, max(facmin, fac * (1 / En)^(1 / (order + 1))))
end

function stepsizePI(xn1, z, xn, Enm1; order=1, Atol=1e-3, Rtol=1e-3)
    tol = Atol .+ Rtol * maximum(abs, hcat(xn1, xn), dims=2)
    En = 1 / sqrt(length(xn1)) * norm((xn1 - z) ./ tol)

    k = order + 1
    Ki = 0.3 / k
    Kp = 0.4 / k
    safetyfac = 0.9

    fac = safetyfac * En^(-(Ki + Kp)) * Enm1^Kp

    return En, fac

end

function integ_adapt_h(f, x0, p, h0, tf; tol=1e-4, maxTries=100)

    X = x0'
    t = [0.]

    i = 1
    h = h0
    errn = 1.
    failcount = 0
    iTries = 0
    while t[end] < tf

        xn = X[i, :]
        tn = t[i]
        xn1, z = rk_adapt_step(f, xn, tn, p, h)

        # keep track of number of tries to compute the next step
        iTries += 1

        # compute error estimate and next step size
        errn, fac = stepsizefactor(xn1, z, xn, Atol=tol, Rtol=0) 

        # use PI step control
        # errn, fac = stepsizePI(xn1, z, xn, errn, Atol=tol, Rtol=0)
        
        # update step size
        h *= fac

        # compare normalized error to 1
        if errn < 1
            X = vcat(X, xn1')
            append!(t, tn + h)
            iTries = 0 # reset the local 'try' counter
            i += 1
        else
            failcount += 1
        end

        # kill the loop if we aren't making progress
        if iTries >= maxTries
            error(@sprintf("Failed to make progress\n  i=%d\n  t=%f\n  h=%f\n ", i, t[i], h))
        end

    end

    return X, t, failcount

end


##
h0 = 0.1
X, t, failcount = integ_adapt_h(odesys, u0, p, h0, tf, tol=1e-1)

## Plots
plot(t, X, xlabel="Time t", ylabel="State")
plot!(t -> u_true(t)[1], 0, tf, label="true y1", line=:dash)
plot!(t -> u_true(t)[2], 0, tf, label="true y2", line=:dash)
