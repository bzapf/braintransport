import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--test", choices=["fenics", "numpy"], type=str)
parser.add_argument("-n", "--num", type=int, default=64)
args = vars(parser.parse_args())

import os

def print(x):
    os.system("echo " + str(x))

if args["test"] == "numpy":

    print("Starting scripts")

    import numpy as np

    k = 0
    N = 256
    arrays = []

    print("Starting to create arrays")

    while True:
        k += 1
        arrays.append(np.random.rand(N, N, N).astype(float))
        
        print("Created array " + str(k))

elif args["test"] == "fenics":

    from fenics import *

    print("Imported Fenics")

    try:
        import dolfin as df
        print(df.__version__)
    except:
        print("Could not print version")


    N = args["num"]

    print("N=" + str(N))

    mesh = UnitCubeMesh(N,N,N)

    print("Created mesh")

    V = FunctionSpace(mesh, "CG", 4)

    print("Created function space")

    u = TrialFunction(V)
    uu = Function(V)
    v = TestFunction(V)

    f = Function(V)

    print("Set up functions, defining forms")

    a = inner(grad(u), grad(v)) * dx
    L = inner(f, v) * dx

    print("Trying to solve")

    solve(a == L, uu)

    print("solved, script done")