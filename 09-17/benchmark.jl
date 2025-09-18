using Pkg
Pkg.activate(".")
using BenchmarkTools

t = @benchmark sin(t) setup=(t=rand())