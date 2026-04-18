import os
import openmc
import numpy as np
import math
import sys

from openmc import model
import math

sys.path.append('../')
from chanel_3d import get_plane, set_chanel_3d
from materials_met_1000 import fuel_type1inner
from materials_met_1000 import fuel_type2inner
from materials_met_1000 import fuel_type3inner
from materials_met_1000 import fuel_type4inner
from materials_met_1000 import fuel_type5inner
from materials_met_1000 import coolant, clading, helium, absorber_enriched, absorber_burnup


def assembly_layer_construct(id=100,
                             cell_name="universe",
                             fuel=None,
                             clading=None,
                             coolant=None,
                             absorber_enriched=absorber_enriched,
                             absorber_burnup=absorber_burnup, 
                             zup=None, zdn=None,
                             boundary="transmission/reflective"):
    r_pin = openmc.ZCylinder(r=0.58 / 2)
    r_clad = openmc.ZCylinder(r=0.68 / 2)
    pin_to_pin_dist = 0.995
    edge_length = pin_to_pin_dist / np.sqrt(3)

    fuel_cell = openmc.Cell(fill=fuel, region=-r_pin & -zup & +zdn)
    clad_cell = openmc.Cell(fill=clading, region=+r_pin & -r_clad & -zup & +zdn)
    coolant_cell = openmc.Cell(fill=coolant, region=+r_clad & -zup & +zdn) 
    root_pins = openmc.Universe(cells=(fuel_cell, clad_cell, coolant_cell))

    clad_cell_cool = openmc.Cell(fill=coolant, region=-r_clad & -zup & +zdn)
    root_pins_cool = openmc.Universe(cells=[clad_cell_cool])

    r_pin_abs = openmc.ZCylinder(r=0.35 / 2)
    r_clad_abs = openmc.ZCylinder(r=0.45 / 2)

    abs_cell = openmc.Cell(fill=absorber_burnup, region=-r_pin_abs & -zup & +zdn)
    abs_clad_cell = openmc.Cell(fill=clading, region=+r_pin_abs & -r_clad_abs & -zup & +zdn)
    abs_coolant_cell = openmc.Cell(fill=coolant, region=+r_clad_abs & -zup & +zdn)
    root_pins_abs = openmc.Universe(cells=(abs_cell, abs_clad_cell, abs_coolant_cell))

    pel_cell = openmc.Cell(fill=absorber_enriched, region=-r_pin & -zup & +zdn)
    pel_clad_cell = openmc.Cell(fill=clading, region=+r_pin & -r_clad & -zup & +zdn)
    pel_coolant_cell = openmc.Cell(fill=coolant, region=+r_clad & -zup & +zdn)
    root_pins_pel = openmc.Universe(cells=(pel_cell, pel_clad_cell, pel_coolant_cell))

    assembly_in = openmc.HexLattice(name='Fuel Assembly {}'.format(cell_name))
    assembly_in.pitch = [0.75]
    assembly_in.center = (0.0, 0.0)
    outer_ring = [root_pins_pel] * 6
    inner_ring = [root_pins_pel]
    assembly_in.universes = [outer_ring, inner_ring]

    assembly_out = openmc.HexLattice(name='Fuel Assembly {}'.format(cell_name))
    assembly_out.pitch = [pin_to_pin_dist]
    assembly_out.center = (0.0, 0.0)
    outer_ring = [root_pins] * 30
    ring_1 = [root_pins] * 24
    ring_2 = [root_pins] * 18
    ring_3 = [root_pins] + [root_pins_abs] + [root_pins] + [root_pins_abs] + \
         [root_pins] + [root_pins_abs] + [root_pins] + [root_pins_abs] + \
         [root_pins] + [root_pins_abs] + [root_pins] + [root_pins_abs]
    ring_4 = [root_pins_cool] * 6
    inner_ring = [root_pins_cool]
    assembly_out.universes = [outer_ring, ring_1, ring_2, ring_3, ring_4, inner_ring]

    coolant_between_lattice_and_chanel = openmc.Universe()
    collant_tvs = openmc.Cell(name="Coolant between lattice and chanel", fill=coolant)
    coolant_between_lattice_and_chanel.add_cell(collant_tvs)
    assembly_in.outer = coolant_between_lattice_and_chanel
    assembly_out.outer = coolant_between_lattice_and_chanel

    subassembly_duct = 9.85
    H_out = subassembly_duct
    R_out = H_out / math.sqrt(3)
    p1 = openmc.Plane(a=1 / H_out, b=1 / R_out, c=0, d=0)
    p2 = openmc.XPlane(x0=0)
    p3 = openmc.Plane(a=-1 / H_out, b=1 / R_out, c=0, d=0)

    thickness_out = 0.165
    bord_duct_inner = get_plane(H=H_out - 2 * thickness_out, boundary="transmission")
    bord_duct_outer = get_plane(H=H_out, boundary="transmission")
    duct_cells = set_chanel_3d(bord_outer=bord_duct_outer, bord_inner=bord_duct_inner, p1=p1, p2=p2, p3=p3,
                              z_up=zup, z_down=zdn, material=clading, name="chanel")

    H_fuel_assembly = 9.685
    bord_inner = get_plane(H=H_out, boundary="transmission", )
    bord_outer = get_plane(H=H_fuel_assembly, boundary=boundary)
    outer_coolant_cells = set_chanel_3d(bord_outer=bord_outer, bord_inner=bord_inner, p1=p1, p2=p2, p3=p3,
                                        z_up=zup, z_down=zdn, material=coolant, name="outer coolant")

    H_hole = 2.95
    H_cluster = 3.0 * assembly_in.pitch[0]
    H_duct = 2.32
    svp_ring_radius = 1.34

    edge_outer = H_fuel_assembly / math.sqrt(3.0)
    edge_hole = H_hole / math.sqrt(3.0)
    edge_cluster = H_cluster / math.sqrt(3.0)
    edge_duct = H_duct / math.sqrt(3.0)

    outer_hex = openmc.model.HexagonalPrism(edge_length=edge_outer, orientation='y')
    hole_hex = openmc.model.HexagonalPrism(edge_length=edge_hole, orientation='y')
    cluster_hex = openmc.model.HexagonalPrism(edge_length=edge_cluster, orientation='y')
    duct_hex = openmc.model.HexagonalPrism(edge_length=edge_duct, orientation='y')

    inner_cluster_cell = openmc.Cell(
        name=f"inner_cluster_cell {cell_name}",
        fill=assembly_in,
        region=-cluster_hex & -zup & +zdn
    )

    inner_duct_cell = openmc.Cell(
        name="inner_duct_cell",
        fill=clading,
        region=-duct_hex & ~(-cluster_hex) & -zup & +zdn
    )

    central_coolant_region = -hole_hex & ~(-duct_hex) & -zup & +zdn

    central_coolant_cell = openmc.Cell(
        name="central_coolant_cell",
        fill=coolant,
        region=central_coolant_region
    )

    assembly_cell_region = -outer_hex & ~(-hole_hex) & -zup & +zdn
    assembly_cell = openmc.Cell(
        name="assembly_cell_with_hole",
        fill=assembly_out,
        region=assembly_cell_region
    )

    fuel_assembly_universe = openmc.Universe(name=f"fuel_assembly_universe {cell_name}")
    fuel_assembly_universe.add_cell(inner_cluster_cell)
    fuel_assembly_universe.add_cell(inner_duct_cell)
    fuel_assembly_universe.add_cell(central_coolant_cell)
    fuel_assembly_universe.add_cell(assembly_cell)

    for one_cell in duct_cells.values():
        fuel_assembly_universe.add_cell(one_cell)
    for one_cell in outer_coolant_cells.values():
        fuel_assembly_universe.add_cell(one_cell)

    root_cell = openmc.Cell(name=cell_name, fill=fuel_assembly_universe)
    root_cell.region = -bord_outer[1] & -bord_outer[2] & -bord_outer[3] & -bord_outer[4] & -bord_outer[5] & -bord_outer[6] & -zup & +zdn

    try:
        print(assembly_out.show_indices(num_rings=6))
    except Exception:
        pass

    return root_cell


if __name__ == "__main__":
    z_down = openmc.ZPlane(z0=0.0, boundary_type="reflective")
    z_up = openmc.ZPlane(z0=10.0, boundary_type="reflective")

    fuel_assembly_cell = assembly_layer_construct(
        id=100,
        cell_name="FuelAssembly_Test",
        fuel=fuel_type1inner,
        clading=clading,
        coolant=coolant,
        zup=z_up,
        zdn=z_down,
        boundary="reflective"
    )

    fuel_type1inner.color = "red"
    absorber_enriched.color = "green"
    clading.color = "gray"
    coolant.color = "blue"
    materials = openmc.Materials([
        fuel_type1inner,
        clading,
        coolant,
        absorber_burnup,
        absorber_enriched,
    ])
    materials.export_to_xml()

    root_universe = openmc.Universe(cells=[fuel_assembly_cell])
    geometry = openmc.Geometry(root_universe)
    geometry.export_to_xml()

    mid_z = (z_up.z0 + z_down.z0) / 2.0
    plot = openmc.Plot()
    plot.file_format = "png"
    plot.filename = "debug_cells"
    plot.width = (20.0, 20.0)
    plot.pixels = (800, 800)
    plot.basis = "xy"
    plot.origin = (0.0, 0.0, mid_z)
    plot.color_by = "material"
    plot.colors = {
        fuel_type1inner: "red",
        clading: "gray",
        coolant: "blue",
        absorber_enriched: "green",
    }

    plots = openmc.Plots([plot])

    plots.export_to_xml()
    openmc.plot_geometry()

    print("Geometry plot saved as debug_cells.png")
