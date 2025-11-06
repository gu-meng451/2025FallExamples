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

function integ_adapt_h(f, x0, p, h0, tf; tol=1e-4, maxTries=100)

    X = x0'
    t = [0.]

    i = 1
    h = h0
    # errnm1 = NaN
    failcount = 0
    iTries = 0
    while t[end] < tf

        xn = X[i, :]
        tn = t[i]
        xn1, z = rk_adapt_step(f, xn, tn, p, h)

        ## TODO replace with better method
        errn = norm(xn1 - z)

        iTries += 1
        order = 1
        fac = 0.9
        facmin = 0.1
        facmax = 2.
        Atol = 1e-2
        Rtol = 1e-2
        toln = Atol + Rtol*max(norm(xn1), norm(xn))
        En = errn/toln
        println( (1/En)^(1/(order+1)) )
        # println(min(facmax, max(facmin, fac * (1 / errn)^(1 / (order + 1)))))
        h *= min(facmax, max(facmin, fac * (1/En)^(1/(order+1)) ) )

        if En < 1 
            X = vcat(X, xn1')
            append!(t, tn + h)
            iTries = 0
            
            ## Attempt 1: PI
            # if i == 1
            #     h = h
            #     errnm1 = errn
            # else
            #     p = 2
            #     α = 0.7/p
            #     β = 0.4/p
            #     h *= (tol/errn)^α*(errnm1/tol)^β
            #     errnm1 = errn
            # end

            ## attempt 2
            # fac = 0.9
            # facmin = 0.1
            # facmax = 3.
            # h *= min(facmax, max(facmin, fac * (1 / errn)^(1 / (2 + 1))))

            i += 1
        else
            failcount += 1
        end

        if iTries >= maxTries
            error(@sprintf("Failed to make progress\n  i=%d\n  t=%f\n  h=%f\n ", i, t[i],h))
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
