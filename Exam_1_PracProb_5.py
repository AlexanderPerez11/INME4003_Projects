import numpy as np

C = 0.8417
e = 0.648


def f(x):
    a = 1 - np.exp(x ** (0.22) / C * (np.exp(-C * x ** (0.78)) - 1)) - e

    return a


def dfdx(x):
    return -(0.261375787097541 * (
                -1 + np.exp(-0.8417 * x ** 0.78)) / x ** 0.78 - 0.78 * x ** 2.77555756156289e-17 * np.exp(
        -0.8417 * x ** 0.78)) * np.exp(1.18807175953428 * x ** 0.22 * (-1 + np.exp(-0.8417 * x ** 0.78)))


def newton_taphson(f, dfdx, x0, tol=1e-7, max_iter=10000):
    x = x0
    iter = 0
    error = 1
    while error > tol:
        fx = f(x)
        dfx = dfdx(x)
        x -= fx / dfx
        print(x)
        error = abs(fx)
        print(error)

        iter += 1
        if iter >= max_iter:
            raise ValueError(f"Failed to converge after {max_iter} iterations.")
    print(f"Reached solution after {iter} iterations")
    print(f'Error: {error}')
    return x


NTU = newton_taphson(f, dfdx, x0=1)
print(f"Solution: {round(NTU, 3)}")
