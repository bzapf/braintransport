from dolfin import *
from dolfin_adjoint import *

from pyadjoint import Block
from pyadjoint.overloaded_function import overload_function

import numpy as np
from preprocessing import preprocessing, dof_mapping_scalar

backend_preprocessing = preprocessing


class PreprocessingBlock(Block):
    def __init__(self, func, V, **kwargs):
        super(PreprocessingBlock, self).__init__()
        self.kwargs = kwargs
        self.add_dependency(func)
        self.B = func.function_space()
        self.V = V
        self.dof_map = dof_mapping_scalar(V.mesh(), func.function_space().mesh())

    def __str__(self):
        return 'PreprocessingBlock'

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):
        # dof map

        adj_input = adj_inputs[0]
        adj = adj_input[self.dof_map]
        adjv = Function(self.B)
        adjv.vector()[:] = adj
        return adjv.vector()

    def recompute_component(self, inputs, block_variable, idx, prepared):
        return backend_preprocessing(inputs[0], self.V)


preprocessing = overload_function(preprocessing, PreprocessingBlock)

if __name__ == "__main__":
    testresolution = 10
    testmesh = UnitSquareMesh(testresolution, testresolution)

    bmesh = BoundaryMesh(testmesh, "exterior")

    V = FunctionSpace(bmesh, "CG", 1)
    Vm = FunctionSpace(testmesh, "CG", 1)

    v = Function(V)
    v = interpolate(Constant(1.0), V)

    vm = preprocessing(v, Vm)

    J = assemble(vm * vm * dx)

    Jhat = ReducedFunctional(J, Control(v))

    rate = taylor_test(Jhat, v, v)
    print(rate)