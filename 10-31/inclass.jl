using Pkg
Pkg.activate(".")
using Plots
using LinearAlgebra

## Let's make a multistep predictor/corrector method
# xdot = f(x,p,t)
#
function forwardEuler_step(f, xn, tn, p, h)
    Fn = f(xn, p, tn)
    xn1 = xn + h*Fn
    return xn1, Fn
end

function ab2_step(f, xn, tn, p, h, f0)
    Fn = f(xn,p,tn)
    xn1 = xn + h/2*(3*Fn - f0)
    return xn1,Fn
end

function am2_step(f, xn, tn, p, h, fn, xs; tol=1e-5, maxIter=100)
    # this is also known as the trap rule.
    # xs: the guess of xn1 from the predictor

    tn1 = tn + h
    xn1 = xs
    Fn1 = copy(fn)
    for i in 1:maxIter
        Fn1 = f(xs,p,tn1)
        xn1 .= xn + h/2*( Fn1 + fn )
        err = norm(xs-xn1)
        xs .= xn1 # xs[:] = xn1[:]
        if err < tol
            break
        end
    end

    return xn1, Fn1

end

function integrator(f,x0,p,h,nsteps)

    n = length(x0)
    X = zeros(nsteps+1, n)
    t = 0:h:h*nsteps

    X[1,:] = x0

    # 1st step: 
    # Predictor: Forward Euler
    xs, F0 = forwardEuler_step(f, x0, t[1], p, h)

    # Corrector: trap
    xn1, F1 = am2_step(f, x0, t[1], p, h, F0, xs)

    # load the result into the storage container
    X[2,:] .= xn1

    # increment over steps 2:n
    for i = 2:nsteps

        xn = X[i,:]
        tn = t[i]

        # predictor
        xs, F1 = ab2_step(f, xn, tn, p, h, F0)

        # corrector
        xn2, F2 = am2_step(f,xn,tn,p,h,F1,xs)

        # save
        X[i+1,:] .= xn2

        # swap
        F0 .= F1
        F1 .= F2

    end

    return X, t
end


## demonstration
function odesys(u,p,t)
    A = p
    return A*u
end

A = [0 1; -1 -0.1]
u0 = [1, 1.]
h = 0.1
nsteps = 100

X, t = integrator(odesys,u0,A,h,nsteps)

plot( t, X, xlabel="Time t", ylabel="States")

u_true(t) = exp(A*t)*u0

plot!( t->u_true(t)[1], t[1], t[end], label="true y1", line=:dash)
plot!(t -> u_true(t)[2], t[1], t[end], label="true y2", line=:dash)
