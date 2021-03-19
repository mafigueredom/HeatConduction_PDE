"""
Created on Sun Jul 19 18:48:18 2020

@author: santiagoortiz
"""

from dolfin import Mesh, XDMFFile, FunctionSpace, Function
from vedo.dolfin import plot, screenshot
import os

def vedoPlotter(pygmsh_ms):

    # Reading mesh data stored in .xdmf files.
    mesh = Mesh()
    with XDMFFile(pygmsh_ms) as infile:
        infile.read(mesh)

    # Define variational problem
    V = FunctionSpace(mesh, 'P', 1)
    u = Function(V)

    R_path = "Output_data"
    fU_in = XDMFFile(os.path.join(R_path, 'FuelRod_m.xdmf'))
    fU_in.read_checkpoint(u, "T", 0)

    axes_opts = dict(
        xyGrid = True,
        axesLineWidth = 1,
        xTickColor = 'black',
        yTickColor = 'black',
        xMinorTicks = 1,       # number of minor ticks btw two major ticks
        yMinorTicks = 1,       # number of minor ticks btw two major ticks
        xLabelSize = 0.02,     # size of the numeric labels along axis
        yLabelSize = 0.02, # offset of numeric labels
                    )
    # nipy_spectral, gnuplot
    plot(u, interactive=True, cmap='hot', axes=0, lw=2,
         scalarbar='vertical', wireframe=True, alpha=10.,
         warpZfactor=0.)  # warpZfactor=0.01
    plot()


if __name__ == "__main__":
    pygmsh_ms = 'pygmeshio_data/FuelRod/FuelRod_m.xdmf'
    vedoPlotter(pygmsh_ms)
    screenshot('FuelRod.png')
