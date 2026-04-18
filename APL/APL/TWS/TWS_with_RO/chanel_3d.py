import os
import openmc
import numpy as np
import math


def get_plane(H, boundary="reflective/transmissive" , orientation="x"):
    """H размер под ключ"""
    edge_length = H / math.sqrt(3)
    R = edge_length
    bord={}
    # установкавнешних границ для чехла ТВС
    if orientation == "x":
        bord = {
            1: openmc.Plane(a=-1 / (H / 2), b=0,      c=0, d=1, boundary_type=boundary), # || OX +
            2: openmc.Plane(a=-1 / H,       b=1 / R,  c=0, d=1, boundary_type=boundary),
            3: openmc.Plane(a=1 / H,        b=1 / R,  c=0, d=1, boundary_type=boundary),
            4: openmc.Plane(a=1 / (H / 2),  b=0,      c=0, d=1, boundary_type=boundary), # || OX -
            5: openmc.Plane(a=1 / H,        b=-1 / R, c=0, d=1, boundary_type=boundary),
            6: openmc.Plane(a=-1 / H,       b=-1 / R, c=0, d=1, boundary_type=boundary), 
              }

    if orientation == "y":
        bord = {
            1 : openmc.Plane(a=-1 / H, b=1 / (3 * R), c=0, d=1, boundary_type=boundary),  # 1
            2 : openmc.Plane(a=0, b=1 / (R + R / 2), c=0, d=1, boundary_type=boundary), # 2
            3 : openmc.Plane(a=1 / H, b=1 / (3 * R), c=0, d=1, boundary_type=boundary),  # 3
            4 : openmc.Plane(a=1 / H, b=-1 / (3 * R), c=0, d=1, boundary_type=boundary),  # 4
            5 : openmc.Plane(a=0, b=-1 / (R + R / 2), c=0, d=1, boundary_type=boundary),  # 5
            6 : openmc.Plane(a=-1 / H, b=-1 / (3 * R), c=0, d=1, boundary_type=boundary),  # 6
              }

    return bord

def set_chanel_3d(bord_outer, bord_inner, p1,p2,p3, z_up, z_down, material, name=""):
    """Задание чехла
    bord_outer внешняя граница
    bord_inner внутренняя граница
    z_up верхняя граница чехла
    z_down нижняя граница чехла
    p1,p2,p3 плоскости для отсечения углов у чехла
    material материал
    name имя для идентификации

    возвращает
    cells  набор из 6 ячеек
    """


    cells = {}

    name_cell1 = "{} sector 1".format(name)
    coolant_cell1 = openmc.Cell(name=name_cell1, fill=material)
    coolant_cell1.region = -bord_outer[1]  & +bord_inner[1] & -z_up & +z_down   & -p1 & +p3
    cells[1] = coolant_cell1

    name_cell2 = "{} sector 2".format(name)
    coolant_cell2 = openmc.Cell(name=name_cell2, fill=material)
    coolant_cell2.region = -bord_outer[2] & +bord_inner[2] & -z_up & +z_down   & +p1 & -p2
    cells[2] = coolant_cell2

    name_cell3 = "{} sector 3".format(name)
    coolant_cell3 = openmc.Cell(name=name_cell3, fill=material)
    coolant_cell3.region = -bord_outer[3] & +bord_inner[3] & -z_up & +z_down   & +p2 & +p3
    cells[3] = coolant_cell3

    name_cell4 = "{} sector 4".format(name)
    coolant_cell4 = openmc.Cell(name=name_cell4, fill=material)
    coolant_cell4.region = -bord_outer[4]  & +bord_inner[4] & -z_up & +z_down  & -p3 & +p1
    cells[4] = coolant_cell4

    name_cell5 = "{} sector 5".format(name)
    coolant_cell5 = openmc.Cell(name=name_cell5, fill=material)
    coolant_cell5.region = -bord_outer[5]  & +bord_inner[5] & -z_up & +z_down  & -p1 & +p2
    cells[5] = coolant_cell5

    name_cell6 = "{} sector 6".format(name)
    coolant_cell6 = openmc.Cell(name=name_cell6, fill=material)
    coolant_cell6.region = -bord_outer[6] & +bord_inner[6] & -z_up & +z_down  & -p2 & -p3
    cells[6] = coolant_cell6

    return cells




if __name__ == "__main__":
    os.environ["OPENMC_CROSS_SECTIONS"] = "/home/sparrow/APL/sections/endfb-viii.0-hdf5/cross_sections.xml"

    #import h5py

    #filename = "/home/adminsrv/projects/sections/endfb71_hdf5_burn/Mo092.h5"
    #file = h5py.File(filename, "r")

    fuel_temperature = 534.0 + 273.15
    coolant_temperature = 432.5 + 273.15
    structure_temperature = 432.5 + 273.15

    # Add nuclides concentration from outer boundary distance 85.82 cm to uo2
    uo2 = openmc.Material(1, "uo2")
    uo2.add_nuclide('U234', 1.6552E-06)
    uo2.add_nuclide('U235', 3.2250E-05)
    uo2.add_nuclide('U236', 1.4710E-06)
    uo2.add_nuclide('U238', 1.8359E-02)

    uo2.add_nuclide('Np237', 1.0175E-04)

    uo2.add_nuclide('Pu236', 6.8053E-10)
    uo2.add_nuclide('Pu238', 1.6416E-04)
    uo2.add_nuclide('Pu239', 2.8416E-03)
    uo2.add_nuclide('Pu240', 1.7508E-03)
    uo2.add_nuclide('Pu241', 2.8697E-04)
    uo2.add_nuclide('Pu242', 4.1028E-04)

    uo2.add_nuclide('Am241', 1.9339E-04)
    uo2.add_nuclide('Am242_m1', 1.2064E-05)
    uo2.add_nuclide('Am243', 1.3206E-04)

    uo2.add_nuclide('Cm242', 5.1976E-06)
    uo2.add_nuclide('Cm243', 5.9372E-07)
    uo2.add_nuclide('Cm244', 7.8359E-05)
    uo2.add_nuclide('Cm245', 1.9913E-05)
    uo2.add_nuclide('Cm246', 1.0410E-05)

    uo2.add_nuclide('Zr90', 3.745663E-03)
    uo2.add_nuclide('Zr91', 8.168384E-04)
    uo2.add_nuclide('Zr92', 1.248554E-03)
    uo2.add_nuclide('Zr94', 1.265299E-03)
    uo2.add_nuclide('Zr96', 2.038456E-04)

    uo2.add_nuclide('Mo92', 9.819924E-05)
    uo2.add_nuclide('Mo94', 6.120910E-05)
    uo2.add_nuclide('Mo95', 1.053458E-04)
    uo2.add_nuclide('Mo96', 1.103749E-04)
    uo2.add_nuclide('Mo97', 6.319426E-05)
    uo2.add_nuclide('Mo98', 1.596730E-04)
    uo2.add_nuclide('Mo100', 6.372364E-05)

    uo2.temperature = fuel_temperature

    # Na
    coolant = openmc.Material(name='Sodium')
    coolant.add_nuclide('Na23', 2.227200E-02)  # !!!!!!!!
    coolant.temperature = coolant_temperature

    # Na
    coolant2 = openmc.Material(name='Sodium2')
    coolant2.add_nuclide('Na23', 2.227200E-02)  # !!!!!!!!
    coolant2.temperature = coolant_temperature

    # ��������
    clading = openmc.Material(name='Clading')

    clading.add_nuclide('Fe54', 4.074842E-03)
    clading.add_nuclide('Fe56', 6.396630E-02)
    clading.add_nuclide('Fe57', 1.477261E-03)
    clading.add_nuclide('Fe58', 1.965963E-04)

    clading.add_nuclide('Ni58', 2.926222E-04)
    clading.add_nuclide('Ni60', 1.127169E-04)
    clading.add_nuclide('Ni61', 4.900176E-06)
    clading.add_nuclide('Ni62', 1.562038E-05)
    clading.add_nuclide('Ni64', 3.980319E-06)

    clading.add_nuclide('Cr50', 4.504027E-04)
    clading.add_nuclide('Cr52', 8.685568E-03)
    clading.add_nuclide('Cr53', 9.848736E-04)
    clading.add_nuclide('Cr54', 2.451559E-04)

    clading.add_nuclide('Mn55', 4.5921E-04)

    clading.add_nuclide('Mo92', 1.099125E-04)
    clading.add_nuclide('Mo94', 6.851013E-05)
    clading.add_nuclide('Mo95', 1.179115E-04)
    clading.add_nuclide('Mo96', 1.235404E-04)
    clading.add_nuclide('Mo97', 7.073207E-05)
    clading.add_nuclide('Mo98', 1.787188E-04)
    clading.add_nuclide('Mo100', 7.132459E-05)

    clading.temperature = structure_temperature

    mats = openmc.Materials((uo2, coolant, clading, coolant2))
    mats.export_to_xml()

    # z boundaries
    boundary = "reflective"
    z_up = openmc.ZPlane(z0=17.16, boundary_type=boundary)
    z_down = openmc.ZPlane(z0=-17.16, boundary_type=boundary)

    ########################################
    # tvel geome
    r_pin = openmc.ZCylinder(r=0.5 / 2, )
    fuel_cell = openmc.Cell()
    fuel_cell.fill = uo2
    fuel_cell.region = -r_pin & -z_up & +z_down

    r_clad = openmc.ZCylinder(r=0.53 / 2)
    clad_cell = openmc.Cell(fill=clading, region=+r_pin & -r_clad & -z_up & +z_down)

    pin_to_pin_dist = 0.80  # cm 0.12526, то есть 0.77134 + 0.12526

    H = pin_to_pin_dist
    edge_length = H / math.sqrt(3)
    R = edge_length

    p1 = openmc.Plane(a=1 / H, b=1 / R, c=0, d=0) # 120
    p2 = openmc.XPlane(x0=0) # 0 
    p3 = openmc.Plane(a=-1 / H, b=1 / R, c=0, d=0) # 60

    bord_inner = get_plane(pin_to_pin_dist - 0.1, boundary="transmission")

    coolant_cell_inner = openmc.Cell(name='coolant_cell_inner', fill=coolant2)
    coolant_cell_inner.region = +r_clad & -bord_inner[1] & -bord_inner[2] & -bord_inner[3] & -bord_inner[4] & - \
    bord_inner[5] & -bord_inner[6] & -z_up & +z_down


    bord_outer = get_plane(H=pin_to_pin_dist, boundary="reflective")

    cells = set_chanel_3d(bord_outer, bord_inner, p1, p2, p3, z_up, z_down, material=coolant, name="chanel")

    # def set_chanel_3d(bord_outer, bord_inner, p1,p2,p3, z_up, z_down, material, name=""):

    pin_universe = openmc.Universe()
    pin_universe.add_cell(fuel_cell)
    pin_universe.add_cell(clad_cell)
    pin_universe.add_cell(coolant_cell_inner)
    for one_cell in cells.values():
        pin_universe.add_cell(one_cell)

    geom = openmc.Geometry(pin_universe)
    geom.export_to_xml()


    # Plot
    p = openmc.Plot()
    p.filename = 'cell_plot_example'
    p.width = (2, 2)
    p.pixels = (1000, 1000)
    p.color_by = 'material'
    p.colors = {uo2: 'red', coolant: 'blue', clading: 'yellow', coolant2: 'green'}

    plots = openmc.Plots([p])
    plots.export_to_xml()
    openmc.plot_geometry()

    # #Computing settings
    # batches = 100
    # inactive =50
    # particles =50000

    # Computing settings
    batches = 50
    inactive = 10
    particles = 5000

    set = openmc.Settings()
    set.batches = batches
    set.inactive = inactive
    set.particles = particles
    set.output = {'tallies': True}
    set.temperature = {"method": "interpolation"}
    set.export_to_xml()

    # Creating an MGXS-library
    group_edges = np.array([0.0, 0.46500E+01, 0.10000E+02, 0.21500E+02, 0.46500E+02, 0.10000E+03,
                                   0.21500E+03, 0.46500E+03, 0.10000E+04, 0.21500E+04, 0.46500E+04,
                                   0.10000E+05, 0.21500E+05, 0.46500E+05, 0.10000E+06, 0.20000E+06,
                                   0.40000E+06, 0.80000E+06, 0.14000E+07, 0.25000E+07, 0.40000E+07, 0.65000E+07,
                                   20.0e6])
    groups = openmc.mgxs.EnergyGroups(group_edges)

    mgxs_lib = openmc.mgxs.Library(geom)
    mgxs_lib.energy_groups = groups

    mgxs_lib.mgxs_types = ['transport', 'nu-fission', 'fission', 'nu-scatter matrix', 'chi', 'absorption']

    mgxs_lib.domain_type = 'material'
    mgxs_lib.domains = geom.get_all_materials().values()

    mgxs_lib.by_nuclide = False
    mgxs_lib.build_library()

    tallies = openmc.Tallies()

    mgxs_lib.add_to_tallies_file(tallies, merge=True)

    tallies.export_to_xml()

    openmc.run(threads=4)

    sp = openmc.StatePoint('statepoint.50.h5')
    mgxs_lib.load_from_statepoint(sp)

    fuel_mgxs = mgxs_lib.get_mgxs(uo2, 'nu-fission')
    fuel_mgxs_1 = mgxs_lib.get_mgxs(uo2, 'absorption')

    fuel_mgxs_2 = mgxs_lib.get_mgxs(uo2, 'fission')
    fuel_mgxs_3 = mgxs_lib.get_mgxs(uo2, 'transport')
    fuel_mgxs_4 = mgxs_lib.get_mgxs(uo2, 'nu-scatter matrix')
    fuel_mgxs_5 = mgxs_lib.get_mgxs(uo2, 'chi')

    fuel_mgxs.print_xs()
    fuel_mgxs_1.print_xs()
    fuel_mgxs_2.print_xs()
    fuel_mgxs_3.print_xs()
    fuel_mgxs_4.print_xs()
    fuel_mgxs_5.print_xs()

    mgxs_lib.build_hdf5_store(filename='mgxs.h5')
