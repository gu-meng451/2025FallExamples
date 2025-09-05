function ln(x; tol=1e-8, maxiter=1000)

    # (x-1) - (x-1)^2/2 + (x-1)^3/3 - ...

    flag = 0
    z = 0.
    i = 0
    term = -1.
    while flag == 0
        i += 1

        term *= -(x-1)
        
        z += term/i

        err = abs(term/i)

        if err < tol
            flag = 1
            return z
        elseif i > maxiter
            error("Max iterations reached")
        end

    end

end