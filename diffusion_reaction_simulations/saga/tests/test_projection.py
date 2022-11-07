from dolfin import *
from dolfin_adjoint import *

testmesh = UnitSquareMesh(4,4)

V1 = FunctionSpace(testmesh, "CG", 1)
V2 = FunctionSpace(testmesh, "DG", 0)

u = Function(V1)

u2 = project(u, V2)

print("success without specifying solver")

u2 = project(u, V2, solver_type="lu")

print("success with specifying solver")