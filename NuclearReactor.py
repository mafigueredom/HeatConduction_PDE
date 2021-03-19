"""
Created on Thu Jul 16 19:25:12 2020

@author: santiagoortiz
"""
from pygmeshio import Generate_Rod_pygmsh
from Solving_CladdingRod import CladdingR_solver
from Solving_FuelRod import FuelR_solver

from Vedo_u import vedoPlotter
from math import pi

R_path = "Output_data"
Data_postprocessing = True
Execute_solvers = True

rU = 4  # Fuel pellet radius (mm)
rZi = 4.08  # Cladding outer radius (mm)
rZo = 4.65  # Cladding inner radius (mm)
meszf = (0.05, 0.03)  # mesh element size factor
r_s = rU, rZi, rZo

# Boundary ids
FRod_f = 10  # Bnd_id: Fuel Rod facet
FRod_s = 100  # Bnd_id: Fuel Rod surface
FClad_fi = 20  # Bnd_id: Fuel Cladding inner facet
FClad_fo = 30  # Bnd_id: Fuel Cladding outer facet
FClad_s = 200  # Dmn_id: Fuel Cladding surface

Physical_ids = ((FRod_f, FRod_s), (FClad_fi, FClad_fo, FClad_s))
pygmsh_msh, pygmsh_fcs = Generate_Rod_pygmsh(r_s, meszf, Physical_ids,
                                           MainF='pygmeshio_data')

FuelRod_mesh, FuelCladding_mesh = pygmsh_msh
FuelRod_facets, FuelCladding_facets = pygmsh_fcs

# Physical properties - UO2
Ql = 300/10  # Linear heat rate, W/cm - W/mm
Qv = Ql/(pi*pow(rU, 2))  # Heat generation per unit volume, (W/mm³)
Tfr_outer = 420.+273.15  # Inner boundary fuel rod temperature, K
FR_Tave = 1120  # Average temperature for fuel rod, K

# Physical properties - ZrO2
h = 41*1000*(1/1000**2)     # convection transfer, kW/(m^2 K) - W/(mm^2 K)
Tfc_inner = 360.+273.15     # Inner boundary fuel cladding temperature, K
Tb = 300+273.15             # bulk temperature, K
FC_Tave = (Tfc_inner+Tb)/2  # Average temperature for fuel cladding, K

if Execute_solvers:
    # Solving fuel rod model
    FRod_info = Qv, Tfr_outer, FR_Tave
    FuelR_solver(FuelRod_mesh, FuelRod_facets, FRod_f,
                 FRod_info, Data_postprocessing, R_path)

    # Solving fuel cladding model
    FClad_fs = FClad_fi, FClad_fo
    FClad_info = h, Tfc_inner, Tb, FC_Tave
    CladdingR_solver(FuelCladding_mesh, FuelCladding_facets,
                 FClad_fs, FClad_info, Data_postprocessing, R_path)


