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

openmc.config['cross_sections'] = "/home/sparrowmsu/apl/sections/endfb-viii.0-hdf5/cross_sections.xml"

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

tvs0 = [make_fa(fuel_type1inner, i) for i in range(1, 97)] # 15
tvs1 = [make_fa(fuel_type2inner, i + 96) for i in range(1, 13)] # 13
tvs2 = [] # [make_fa(fuel_type3inner, i + ) for i in range(1, )] 
tvs3 = [] # [make_fa(fuel_type4inner, i + ) for i in range(1, )]

tvs1_ro1 = [make_fa(fuel_type2inner, i + 108, 6) for i in range(1, 2)] # 13
tvs0_ro1 = [make_fa(fuel_type1inner, i + 109, 6) for i in range(1, 7)] # 15

tvs1_ro2 = [make_fa(fuel_type2inner, i + 115, 5) for i in range(1, 7)] # 13
tvs0_ro3 = [] # [make_fa(fuel_type1inner, i, 6) for i in range(1, 7)] # 15
tvs0_ro4 = []
tvs0_ro5 = []
tvs0_ro6 = []
tvs0_ro7 = []

all_tvs = (
    tvs0 + tvs1 + tvs2 + tvs3 +
    tvs1_ro1 + tvs0_ro1 + tvs1_ro2 +
    tvs0_ro3 + tvs0_ro4 + tvs0_ro5 + tvs0_ro6 + tvs0_ro7
)

tvs_water = []
for i in range(7):
    wu = openmc.Universe(name=f"water_fa_{i + 4000}")
    wu.add_cell(openmc.Cell(fill=coolant))
    tvs_water.append(wu)

# -------------------------------
# RINGS
# -------------------------------

j0 = 0          
i0_ro1 = 0      
i1 = 0         
j1_ro2 = 0     

ring6 = []
for k in range(6):
    ring6.append(tvs_water[k])
    ring6 += tvs0[j0:j0 + 5]
    j0 += 5

ring5 = tvs0[j0:j0 + 30]
j0 += 30

ring4 = tvs0[j0:j0 + 24]
j0 += 24

ring3 = []
for _ in range(6):
    ring3 += tvs0_ro1[i0_ro1:i0_ro1 + 1]
    i0_ro1 += 1

    ring3 += tvs0[j0:j0 + 2]
    j0 += 2

second_ring = []
for _ in range(6):
    second_ring += tvs1[i1:i1 + 1]
    i1 += 1

    second_ring += tvs1_ro2[j1_ro2:j1_ro2 + 1]
    j1_ro2 += 1

first_ring = tvs1[i1:i1 + 6]
i1 += 6

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
# MATS & CROSS_SECTION
# -------------------------------
openmc.config['cross_sections'] = "/home/sparrowmsu/apl/sections/endfb-viii.0-hdf5/cross_sections.xml"
mats = openmc.Materials(materials_lst)
mats.export_to_xml()

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

fa_tally = openmc.Tally(name='flux_heat')
fa_tally.filters = [openmc.UniverseFilter(all_tvs)]
fa_tally.filters.append(energy_filter)
fa_tally.scores = ['flux', 'heating']

tallies_file.append(fa_tally)
tallies_file.export_to_xml()

if RUN_OPENMC:
    openmc.run()

sp = openmc.StatePoint('statepoint.10.h5')   
t = sp.get_tally(name='flux_heat')

heat_idx = t.get_score_index('heating')
data = t.get_reshaped_data(value='mean')
print("shape =", data.shape)

for i, u in enumerate(all_tvs):

    uid = u.id

    heating_total = np.sum(data[i, :, 0, heat_idx])

    print(f"TVS universe_id = {uid:6d}   mean heating = {heating_total:.6e}")

#mesh = openmc.RegularMesh()
## задать mesh только на область нужной ТВС
#
#t_flux = openmc.Tally(name='fa1007_flux_map')
#t_flux.filters = [
#    openmc.UniverseFilter([fa_id]),
#    openmc.MeshFilter(mesh)
#]
#t_flux.scores = ['flux']