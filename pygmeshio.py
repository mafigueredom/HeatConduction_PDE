"""
Created on Sun Jul 12 22:12:59 2020
@author: santiagoortiz
"""
import pygmsh
import meshio
import numpy as np
import os


def Generate_Rod_pygmsh(r_s, meszf, Physical_ids, MainF='pygmeshio_data'):

    # Pygmeshio files
    Domains = ['FuelRod', 'FuelCladding']
    pygmsh_geos = [os.path.join(MainF, DMn, "{}.geo".format(DMn)) for DMn in Domains]
    pygmsh_mshs = [os.path.join(MainF, DMn, "{}.msh".format(DMn)) for DMn in Domains]
    pygmsh_ms = [os.path.join(MainF, DMn, "{}_m.xdmf".format(DMn)) for DMn in Domains]
    pygmsh_fs = [os.path.join(MainF, DMn, "{}_f.xdmf".format(DMn)) for DMn in Domains]

    xy0 = [0.0, 0.0, 0.0]  # Fuel rod center

    FR_meszf, FC_meszf = meszf
    rU, rZi, rZo = r_s

    FR_f, FR_s = tuple([str(FR_id) for FR_id in Physical_ids[0]])
    FC_fi, FC_fo, FC_s = tuple([str(FC_id) for FC_id in Physical_ids[1]])

    geomFR = pygmsh.opencascade.Geometry(FR_meszf, FR_meszf)
    geomFC = pygmsh.opencascade.Geometry(FC_meszf, FC_meszf)

    FuelRod = geomFR.add_disk(xy0, rU)
    InnerCrl = geomFC.add_disk(xy0, rZi)
    OuterCrl = geomFC.add_disk(xy0, rZo)
    FuelCladding = geomFC.boolean_difference([OuterCrl], [InnerCrl])

    geomFR.add_raw_code('Physical Curve('+FR_f+') = {' + FuelRod.id + '};')
    geomFR.add_raw_code('Physical Surface('+FR_s+') = {' + FuelRod.id + '};')
    geomFC.add_raw_code('Physical Curve('+FC_fi+') = {' + InnerCrl.id + '};')
    geomFC.add_raw_code('Physical Curve('+FC_fo+') = {' + OuterCrl.id + '};')
    geomFC.add_raw_code('Physical Surface('+FC_s+') = {' + FuelCladding.id + '};')

    geom_DMs = [geomFR, geomFC]


    formats_zip = zip(pygmsh_geos, pygmsh_mshs, pygmsh_ms, pygmsh_fs, geom_DMs)
    for pygmsh_geo, pygmsh_msh, pygmsh_m, pygmsh_f, geom_DM in formats_zip:

        mesh = pygmsh.generate_mesh(geom_DM, dim=2, prune_z_0=True,
                                    geo_filename=pygmsh_geo,
                                    msh_filename=pygmsh_msh)
        # Meshio conversion
        cells = np.vstack(np.array([cells.data for cells in mesh.cells
                                    if cells.type == "triangle"]))
        triangle_mesh = meshio.Mesh(points=mesh.points, cells=[("triangle", cells)])

        facet_cells = np.vstack(np.array([cells.data for cells in mesh.cells
                                          if cells.type == "line"]))
        facet_data = mesh.cell_data_dict["gmsh:physical"]["line"]
        facet_mesh = meshio.Mesh(points=mesh.points,
                                 cells=[("line", facet_cells)],
                                 cell_data={"name_to_read": [facet_data]})
        # Write mesh
        meshio.xdmf.write(pygmsh_m, triangle_mesh)
        meshio.xdmf.write(pygmsh_f, facet_mesh)

    return pygmsh_ms, pygmsh_fs
