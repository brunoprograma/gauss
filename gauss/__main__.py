from .metodos import *

if __name__ == "__main__":
    M = [
        [10,2,1,7],
        [1,5,1,-8],
        [2,3,10,6]
    ]

    print(M)
    jacobi(M, E=0.000000000000001)
    seidel(M, E=0.000000000000001)