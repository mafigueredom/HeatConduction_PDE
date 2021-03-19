"""
Created on Thu Jul 23 21:40:06 2020

@author: santiagoortiz
"""
from dolfin import Mesh, XDMFFile, MeshValueCollection, cpp, FunctionSpace, \
                   Function, TestFunction, Constant, DirichletBC, \
                   TrialFunction, dot, grad, dx, derivative, solve
import os
from PhyProperties import K_UO2


def FuelR_solver(FR_mesh, FR_facets, FR_f, FRod_info, D_post, R_path):

    Qv, Tfr_outer, FR_Tave = FRod_info

    # Reading mesh data stored in .xdmf files.
    mesh = Mesh()
    with XDMFFile(FR_mesh) as infile:
        infile.read(mesh)

    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile(FR_facets) as infile:
        infile.read(mvc, "name_to_read")
    mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
    #File("Circle_facet.pvd").write(mf)

    # Mesh data
    print('FuelRod_mesh data\n',
          'Number of cells: ', mesh.num_cells(),
          '\n Number of nodes: ', mesh.num_vertices())

    # Define variational problem
    V = FunctionSpace(mesh, 'P', 1)
    u = Function(V)
    v = TestFunction(V)

    f = Constant(Qv)

    # Define boundary conditions base on pygmsh mesh mark
    bc = DirichletBC(V, Constant(Tfr_outer), mf, FR_f)

    # Variational formulation
    # Klbda_UO2 = Constant(K_UO2(FR_Tave))
    F = K_UO2(u)*dot(grad(u), grad(v))*dx - f*v*dx

    # Compute solution
    du = TrialFunction(V)
    Gain = derivative(F, u, du)
    solve(F==0, u, bc, J=Gain, \
          solver_parameters={"newton_solver": {"linear_solver": "lu",
                                               "relative_tolerance": 1e-9}}, \
          form_compiler_parameters={"cpp_optimize": True,
                                    "representation": "uflacs",
                                    "quadrature_degree" : 2}
          )  # mumps

    # Save solution
    if D_post:
        fU_out = XDMFFile(os.path.join(R_path, 'FuelRod', 'u.xdmf'))
        fU_out.write_checkpoint(u, "T", 0, XDMFFile.Encoding.HDF5, False)
        fU_out.close()

# In-line visualization with "vedo" module
#vedoPlotter(FuelRod_mesh)

