# make an example pyton implementation of exp(x)

# optional keyword arguments tolerance=1e-10, max_iter=1000
def exp(x, tolerance=1e-10, max_iter=1000):
    y = 1.0
    term = 1.0
    n = 1
    while term > tolerance and n < max_iter:
        term *= x / n
        y += term
        n += 1
    if n == max_iter:
        print("Warning: max_iter reached")
    return y

# test the function
if __name__ == "__main__":  
    import math
    x = 100.0
    print(f"exp({x}) = {exp(x)}")
    print(f"math.exp({x}) = {math.exp(x)}")
    print(f"Difference = {abs(exp(x) - math.exp(x))}")
    #print relative difference
    # print(f"Relative Difference = {abs(exp(x) - math.exp(x)) / math.exp(x)}")