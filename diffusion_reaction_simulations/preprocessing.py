from dolfin import *
import numpy as np


def dof_mapping_scalar(mesh, bmesh):
    Vb = FunctionSpace(bmesh, "CG", 1)
    W = FunctionSpace(mesh, "CG", 1)

    # from boundary mesh to mesh
    mapb = bmesh.entity_map(0)
    d2v = dof_to_vertex_map(Vb)
    v2d = vertex_to_dof_map(W)

    f = Function(Vb)

    Vb_to_W_map = np.zeros(np.size(f.vector().get_local()))
    for i in range(np.size(f.vector().get_local())):
        GVertID = Vertex(bmesh, d2v[i]).index()  # Local Vertex ID for given dof on boundary mesh
        PVertID = mapb[GVertID]  # Local Vertex ID of parent mesh
        PDof = v2d[PVertID]
        Vb_to_W_map[i] = int(PDof)

    return Vb_to_W_map

def preprocessing(vb, V):
    """
    :param vb: Function on vb, scalar CG1 function on boundary
    :param V: FunctionSpace CG1 on mesh
    :param dofmap: dof-mapping CG1-boundary --> CG1-mesh
    :return: v in V, v|_boundary = vb, 0 on other nodes
    """
    dofmap = dof_mapping_scalar(V.mesh(), vb.function_space().mesh())
    dofs = vb.vector()[:]
    v = Function(V)
    v.vector()[dofmap] = dofs
    return v

def preprocessing_list(vlist, V):
    return [preprocessing(v, V) for v in vlist]

if __name__ == "__main__":
    # mesh
    N = 50
    mesh = UnitSquareMesh(N, N)

    bmesh = BoundaryMesh(mesh, "exterior")

    Vb = FunctionSpace(bmesh, "CG", 1)
    V = FunctionSpace(mesh, "CG", 1)

    vb = interpolate(Expression("x[0]", degree = 1), Vb)

    v = preprocessing(vb, V)

    file = File('../Output/preproc_test.pvd')
    file << vb
    file << v

