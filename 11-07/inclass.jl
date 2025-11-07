using Pkg
Pkg.activate(".")

using Plots
using NonlinearSolve
using LinearAlgebra

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
function irk_step()
    # TODO make implicit RK method
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
        z = irk_step() # TODO fix function inputs

        # h/2 steps
        xn1_2 = irk_step() # TODO fix function inputs
        xn1 = irk_step() # TODO fix function inputs

        # keep track of number of tries to compute the next step
        iTries += 1

        # use PI step control
        errn, fac = stepsizePI(xn1, z, xn, errn, Atol=tol, Rtol=0)

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

