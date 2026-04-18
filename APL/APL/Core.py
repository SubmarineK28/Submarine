import openmc
import numpy as np
import sys
from pathlib import Path
import matplotlib.pyplot as plt

# Подключаем текущую папку TWS/TWS_with_RO
PROJECT_DIR = Path(__file__).resolve().parent
TWS_DIR = PROJECT_DIR / "TWS" / "TWS_with_RO"
if str(TWS_DIR) not in sys.path:
    sys.path.insert(0, str(TWS_DIR))

from fa_calculation import fa_3d
from materials_met_1000 import (
    fuel_type1inner,
    fuel_type2inner,
    fuel_type3inner,
    fuel_type4inner,
    fuel_type5inner,
    fuel_type6inner,
    coolant,
    clading,
    steel,
    helium,
    absorber_enriched,
    absorber_burnup,
)

fa_step = 10.1
tvs_gap = 0.25
lat_step = fa_step - tvs_gap / 2.0

z_down = 0.0
z_up = 168.0

RUN_OPENMC = True

# -------------------------------
# MATERIALS
# -------------------------------
materials_lst = [
    fuel_type1inner,
    fuel_type2inner,
    fuel_type3inner,
    fuel_type4inner,
    fuel_type5inner,
    fuel_type6inner,
    coolant,
    clading,
    steel,
    helium,
    absorber_enriched,
    absorber_burnup,
]
materials = openmc.Materials(materials_lst)
materials.export_to_xml()

openmc.config['cross_sections'] = "/home/sparrow/APL/sections/endfb-viii.0-hdf5/cross_sections.xml"

# -------------------------------
# OUTER / LATTICE
# -------------------------------
all_water_cell = openmc.Cell(fill=coolant)
outer_universe = openmc.Universe(cells=(all_water_cell,))

lat = openmc.HexLattice()
lat.center = (0.0, 0.0)
lat.pitch = (lat_step,)
lat.outer = outer_universe
lat.orientation = "x"

# -------------------------------
# TVS OBJECTS (аналог FA0..FA3)
# -------------------------------
def make_fa(fuel_material, idx, level = 0, absorber_enriched = absorber_enriched, coolant = coolant):
    return fa_3d(
        tvs_ind=idx,
        fuels_lst=[fuel_material] * 6,
        pel_list=[coolant] * (6 - level) + [absorber_enriched] * level,
        boundary="transmission",
        is_void=False,
    )

tvs0 = [make_fa(fuel_type1inner, i) for i in range(1, 103)] # 15
tvs1 = [make_fa(fuel_type2inner, i + 1000) for i in range(1, 20)] # 13
tvs2 = [make_fa(fuel_type3inner, i + 2000) for i in range(1, 18 + 24 + 30 + 36 + 36 + 1)]
tvs3 = [make_fa(fuel_type4inner, i + 3000) for i in range(7 + 12 + 1 + 2)]

tvs1_ro1 = [make_fa(fuel_type2inner, i + 1000, 6) for i in range(1, 2)] # 13
tvs0_ro1 = [make_fa(fuel_type1inner, i, 6) for i in range(1, 7)] # 15

tvs1_ro2 = [make_fa(fuel_type2inner, i + 1000, 5) for i in range(1, 7)] # 13
tvs0_ro3 = [make_fa(fuel_type1inner, i, 6) for i in range(1, 7)] # 15
tvs0_ro4 = [make_fa(fuel_type1inner, i, 6) for i in range(1, 7)] # 15
tvs0_ro5 = [make_fa(fuel_type1inner, i, 6) for i in range(1, 7)] # 15
tvs0_ro6 = [make_fa(fuel_type1inner, i, 6) for i in range(1, 7)] # 15
tvs0_ro7 = [make_fa(fuel_type1inner, i, 6) for i in range(1, 7)] # 15

tvs_water = []
for i in range(7):
    wu = openmc.Universe(name=f"water_fa_{i + 4000}")
    wu.add_cell(openmc.Cell(fill=coolant))
    tvs_water.append(wu)

# -------------------------------
# RINGS 
# -------------------------------
i = 0
j = 0
ring6 = []
ring6 = ring6 + [tvs_water[0]]
ring6 = ring6 + tvs0[j:j + 5]; j += 5                                                                                                                                                                                                                                                                                                                                                                                   
ring6 = ring6 + [tvs_water[1]]
ring6 = ring6 + tvs0[j:j + 5]; j += 5
ring6 = ring6 + [tvs_water[2]]
ring6 = ring6 + tvs0[j:j + 5]; j += 5
ring6 = ring6 + [tvs_water[3]]
ring6 = ring6 + tvs0[j:j + 5]; j += 5
ring6 = ring6 + [tvs_water[4]]
ring6 = ring6 + tvs0[j:j + 5]; j += 5
ring6 = ring6 + [tvs_water[5]]
ring6 = ring6 + tvs0[j:j + 5]; j += 5

i = 0
j = 0
ring5 = []
ring5 = ring5 + tvs0[j:j + 30]; j += 30

i = 0
j = 0
ring4 = []
ring4 = ring4 + tvs0[j:j + 24]; j += 24

i = 0
j = 0
ring3 = [] # ro_1
ring3 = ring3 + tvs0_ro1[i: i + 1]; i += 1
ring3 = ring3 + tvs0[j:j + 2]; j += 2
ring3 = ring3 + tvs0_ro1[i: i + 1]; i += 1
ring3 = ring3 + tvs0[j:j + 2]; j += 2
ring3 = ring3 + tvs0_ro1[i: i + 1]; i += 1
ring3 = ring3 + tvs0[j:j + 2]; j += 2
ring3 = ring3 + tvs0_ro1[i: i + 1]; i += 1
ring3 = ring3 + tvs0[j:j + 2]; j += 2
ring3 = ring3 + tvs0_ro1[i: i + 1]; i += 1
ring3 = ring3 + tvs0[j:j + 2]; j += 2
ring3 = ring3 + tvs0_ro1[i: i + 1]; i += 1
ring3 = ring3 + tvs0[j:j + 2]; j += 2

i = 0
j = 0
second_ring = [] # ro_2
second_ring = second_ring + tvs1[i:i + 1];  i += 1
second_ring = second_ring + tvs1_ro2[j:j + 1]; j +=1
second_ring = second_ring + tvs1[i:i + 1];  i += 1
second_ring = second_ring + tvs1_ro2[j:j + 1]; j +=1
second_ring = second_ring + tvs1[i:i + 1];  i += 1
second_ring = second_ring + tvs1_ro2[j:j + 1]; j +=1
second_ring = second_ring + tvs1[i:i + 1];  i += 1
second_ring = second_ring + tvs1_ro2[j:j + 1]; j +=1
second_ring = second_ring + tvs1[i:i + 1];  i += 1
second_ring = second_ring + tvs1_ro2[j:j + 1]; j +=1
second_ring = second_ring + tvs1[i:i + 1];  i += 1
second_ring = second_ring + tvs1_ro2[j:j + 1]; j +=1


first_ring = []
first_ring = first_ring + tvs1[i:i + 6]; i += 6

inner_ring = [] # ro_1
inner_ring = [tvs1_ro1[0]]

lat.universes = [ring6, ring5, ring4, ring3, second_ring, first_ring, inner_ring]

print("Количество ТВС в решетке:", len(ring6) + len(ring5) + len(ring4) + len(ring3) + len(second_ring) + len(first_ring) + len(inner_ring))


# -------------------------------
# GEOMETRY
# -------------------------------
h_cluster = 200.0 # thick = (150 - 131.3)/2 = 9.35
steel_thickness = 5.0
h_duct = h_cluster + 2.0 * steel_thickness  # 210 см

cluster_border = openmc.model.HexagonalPrism(
    edge_length=h_cluster / np.sqrt(3.0),
    orientation="x",
    boundary_type="transmission",
)

duct_border = openmc.model.HexagonalPrism(
    edge_length=h_duct / np.sqrt(3.0),
    orientation="x",
    boundary_type="vacuum",
)

z_plane_down = openmc.ZPlane(z0=z_down, boundary_type="vacuum")
z_plane_up = openmc.ZPlane(z0=z_up, boundary_type="vacuum")

core_cell = openmc.Cell(fill=lat, region=-cluster_border & +z_plane_down & -z_plane_up)
steel_shell_cell = openmc.Cell(fill=steel, region=+cluster_border & -duct_border & +z_plane_down & -z_plane_up)

main_univ = openmc.Universe(cells=[core_cell, steel_shell_cell])
geomi = openmc.Geometry(main_univ)
geomi.export_to_xml()


# -------------------------------
# PLOTS
# -------------------------------
p = openmc.Plot()
p.basis = "xy"
p.origin = (0.0, 0.0, 25.0)
p.filename = "core_xy"
p.width = (300.0, 300.0)
p.pixels = (5000, 5000)
p.color_by = "material"
p.colors = {
    fuel_type1inner: "coral",
    fuel_type2inner: "darkgray",
    fuel_type3inner: "darkgray",
    fuel_type4inner: "dimgray",
    fuel_type5inner: "gray",
    fuel_type6inner: "black",
    coolant: "blue",
    clading: "silver",
    steel: "slategray",
    helium: "whitesmoke",
    absorber_enriched: "green",
    absorber_burnup: "yellow",
}

plots = openmc.Plots([p])
plots.export_to_xml()
openmc.plot_geometry()

p.origin = (0.0, 0.0, (z_down + z_up) / 2.0)
p.filename = "core_yz"
p.basis = "yz"
p.width = (300.0, z_up - z_down)
p.pixels = (5000, 5000)
p.color_by = "material"
p.colors = {
    fuel_type1inner: "white",
    fuel_type2inner: "lightgray",
    fuel_type3inner: "darkgray",
    fuel_type4inner: "dimgray",
    fuel_type5inner: "gray",
    fuel_type6inner: "black",
    coolant: "blue",
    clading: "silver",
    steel: "slategray",
    helium: "whitesmoke",
    absorber_enriched: "green",
    absorber_burnup: "yellow",
}

plots = openmc.Plots([p])
plots.export_to_xml()
openmc.plot_geometry()

# -------------------------------
# SETTINGS
# -------------------------------
point = openmc.stats.Point((0, 0, 80))
src = openmc.Source(space=point)

settings = openmc.Settings()
settings.source = src
settings.batches = 10
settings.inactive = 1
settings.particles = 20
settings.temperature = {'multipole': True, 'method': 'interpolation', 'range': [290, 2500]}
settings.export_to_xml()

tallies_file = openmc.Tallies()

energy_filter = openmc.EnergyFilter([0., 0.625, 20.0e6])

mesh = openmc.RegularMesh()
mesh.dimension = (360, 210, 60)          
mesh.lower_left = (-180, -105, -30)
mesh.upper_right = (180.0, 105, 30)

mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name='flux')
tally.filters = [mesh_filter]
tally.filters.append(energy_filter)
tally.scores = ['flux']
tallies_file.append(tally)
tallies_file.export_to_xml()


if RUN_OPENMC:
    openmc.run()

sp = openmc.StatePoint('statepoint.10.h5')   
t = sp.get_tally(name='flux')

data = t.get_reshaped_data(value='mean', expand_dims=True).squeeze()
thermal = data[..., 0]  
fast = data[..., 1]      

nx, ny, nz = 360, 210, 60
x0, y0, z0 = -180, -105, -30
x1, y1, z1 = 180, 105, 30

# 3) 2D срез (например, по центру z)
k = nz // 2
plt.figure()
plt.imshow(thermal[:, :, k].T, origin='lower', extent=[x0, x1, y0, y1], aspect='auto')
plt.colorbar(label='Flux')
plt.xlabel('x, cm')
plt.ylabel('y, cm')
plt.title('Thermal flux, z=0 slice')
plt.show()
