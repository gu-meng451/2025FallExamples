
function myexp(x; tol= 1e-10, maxiter=100)
    term = 1.0
    y = term
    n = 1
    flag = 0
    while flag == 0
        term *= x / n
        y += term
        n += 1
        if abs(term) < tol
            flag = 1
            return y,n-1
        elseif n > maxiter
            flag = -1
            error("Max iterations reached")
        end
    end

end

(exp(100.0) - myexp(100.0, maxiter=2000)[1])/exp(100.0)
