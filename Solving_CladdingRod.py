"""
Created on Thu Jul 23 21:40:06 2020

@author: santiagoortiz
"""
from dolfin import Mesh, XDMFFile, MeshValueCollection, cpp, FunctionSpace, \
                   Function, TestFunction, Constant, DirichletBC, Measure, \
                   TrialFunction, dot, grad, dx, derivative, solve
import os
from PhyProperties import K_ZrO2


def CladdingR_solver(FC_mesh, FC_facets, FC_fs, FC_info, D_post, R_path):

    FC_fi, FC_fo = FC_fs
    h, Tfc_inner, Tb, FC_Tave = FC_info

    # Reading mesh data stored in .xdmf files.
    mesh = Mesh()
    with XDMFFile(FC_mesh) as infile:
        infile.read(mesh)

    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile(FC_facets) as infile:
        infile.read(mvc, "name_to_read")
    mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
    #File("Circle_facet.pvd").write(mf)

    # Mesh data
    print('CladdingRod_mesh data\n',
          'Number of cells: ', mesh.num_cells(),
          '\n Number of nodes: ', mesh.num_vertices())

    # Define variational problem
    V = FunctionSpace(mesh, 'P', 1)
    u = Function(V)
    v = TestFunction(V)

    # Define boundary conditions base on pygmsh mesh mark
    bc = DirichletBC(V, Constant(Tfc_inner), mf, FC_fi)
    ds = Measure("ds", domain=mesh, subdomain_data=mf, subdomain_id=FC_fo)

    # Variational formulation
    # Klbda_ZrO2 = 18/1000  # W/(m K) --> Nuclear reactor original source
    # Klbda_ZrO2 = Constant(K_ZrO2(FC_Tave))
    F = K_ZrO2(u)*dot(grad(u), grad(v))*dx + h*(u - Tb)*v*ds

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
        fU_out = XDMFFile(os.path.join(R_path, 'FuelCladding', 'u.xdmf'))
        fU_out.write_checkpoint(u, "T", 0, XDMFFile.Encoding.HDF5, False)
        fU_out.close()

# In-line visualization with "vedo" module
#vedoPlotter(FuelCladding_mesh)

