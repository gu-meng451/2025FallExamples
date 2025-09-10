using LinearAlgebra

A = [1 2 3; 4 6 5; -1 0 1.]

L = [1 0 0; 4 1 0; -1 -1 1.]
U = [1 2 3; 0 -2 -7; 0 0 -3.]

L*U - A

