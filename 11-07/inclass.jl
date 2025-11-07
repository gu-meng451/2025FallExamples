using Pkg
Pkg.activate(".")

using Plots
using NonlinearSolve
using LinearAlgebra
using Printf

## Linear Test Problem
function odesys(u, p, t)
    A = p
    return A * u
end

A = [0 1; -1 -0.1]
p = A
u0 = [1, 1.]
tf = 10.
u_true(t) = exp(A * t) * u0


## implicit RK step
function irk_step(f, xn, tn, p, h)

    c1 = 1/2 - sqrt(3)/6
    c2 = 1/2 + sqrt(3)/6

    a11 = 1/4
    a12 = 1/4-sqrt(3)/6
    a21 = 1/4+sqrt(3)/6
    a22 = 1/4

    b1 = 1/2
    b2 = 1/2

    # 0 = -X1 + xn + h*a11*f(t+c1*h,X1) + h*a12*f(t+c2*h,X2)
    # 0 = -X2 + xn + h*a21*f(t+c1*h,X1) + h*a22*f(t+c2*h,X2)
    #
    X0 = hcat(xn,xn)
    R(X,p) = hcat( -X[:,1] + xn+ h*a11*f(X[:,1],p,tn+c1*h) + h*a12*f(X[:,2],p,tn+c2*h),
    -X[:,2] + xn+ h*a21*f(X[:,1],p,tn+c1*h) + h*a22*f(X[:,2],p,tn+c2*h))

    prob = NonlinearProblem(R, X0, p)
    sol = solve(prob, NewtonRaphson())

    if sol.retcode != ReturnCode.Success
        error("Failed to converge")
    end

    # update next time step
    xn1 = xn + h * b1 * f(sol.u[:, 1], p, tn + c1 * h) + h * b2 * f(sol.u[:, 2], p, tn + c2 * h)

    return xn1

end

# irk_step(odesys,[1,2],0,p,0.1)


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

        # h step
        z = irk_step(f, xn, tn, p, h)

        # h/2 steps
        xn1_2 = irk_step(f, xn, tn, p, h/2)
        xn1 = irk_step(f, xn1_2, tn+h/2, p, h/2)

        # keep track of number of tries to compute the next step
        iTries += 1

        # use PI step control
        errn, fac = stepsizePI(xn1, z, xn, errn, Atol=tol, Rtol=0)

        # update step size
        fac
        # h *= fac
        h *= 1

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

