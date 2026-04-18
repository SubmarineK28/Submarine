import openmc

m1 = openmc.Material(name='zone_1')
m1.add_nuclide('U235', 1.0)
m1.set_density('g/cm3', 10.0)

m2 = openmc.Material(name='zone_2')
m2.add_nuclide('U238', 1.0)
m2.set_density('g/cm3', 10.0)

m3 = openmc.Material(name='zone_3')
m3.add_element('Zr', 1.0)
m3.set_density('g/cm3', 6.5)

m4 = openmc.Material(name='zone_4')
m4.add_element('Fe', 1.0)
m4.set_density('g/cm3', 7.8)

m5 = openmc.Material(name='zone_5')
m5.add_nuclide('H1', 2.0)
m5.add_nuclide('O16', 1.0)
m5.set_density('g/cm3', 1.0)

materials = openmc.Materials([m1, m2, m3, m4, m5])
materials.export_to_xml()

def make_universe(material, name):
    cell = openmc.Cell(fill=material)
    return openmc.Universe(name=name, cells=[cell])

u1 = make_universe(m1, 'u1')
u2 = make_universe(m2, 'u2')
u3 = make_universe(m3, 'u3')
u4 = make_universe(m4, 'u4')
u5 = make_universe(m5, 'u5')

outer_mat = openmc.Material(name='outer')
outer_mat.add_nuclide('H1', 2.0)
outer_mat.add_nuclide('O16', 1.0)
outer_mat.set_density('g/cm3', 0.001)

materials.append(outer_mat)
materials.export_to_xml()

outer_cell = openmc.Cell(fill=outer_mat)
outer_univ = openmc.Universe(cells=[outer_cell])

lat = openmc.HexLattice(name='core_121')

lat.center = (0.0, 0.0)
lat.pitch = (9.85,)     
lat.outer = outer_univ

print(openmc.HexLattice.show_indices(num_rings=7, orientation='y'))

void_cell = openmc.Cell(fill=outer_mat)
void_univ = openmc.Universe(cells=[void_cell])

outer_ring = []
for i in range(36):
    if i in (0, 6, 12, 18, 24, 30):
        outer_ring.append(void_univ)
    else:
        outer_ring.append(u5)

ring_5 = [u4] * 30
ring_4 = [u3] * 24
ring_3 = [u2] * 18
ring_2 = [u2] * 12
ring_1 = [u1] * 6
center  = [u1]

lat.universes = [
    outer_ring,
    ring_5,
    ring_4,
    ring_3,
    ring_2,
    ring_1,
    center
]

lat.orientation = 'y'

prism = openmc.model.HexagonalPrism(
    edge_length=7 * lat.pitch[0],
    orientation=lat.orientation,
    boundary_type='vacuum'
)

main_cell = openmc.Cell(
    fill=lat,
    region=-prism
)

root = openmc.Universe(cells=[main_cell])
geom = openmc.Geometry(root)
geom.export_to_xml()

plot = openmc.Plot()
plot.filename = 'core_121_hex'
plot.basis = 'xy'
plot.origin = (0.0, 0.0, 0.0)
plot.width = (160.0, 160.0)  
plot.pixels = (1200, 1200)
plot.color_by = 'material'
plot.colors = {
    m1: 'lightgray',
    m2: 'white', 
    m3: 'silver',
    m4: 'gray',
    m5: 'gainsboro',
    outer_mat: 'white'
}

plots = openmc.Plots([plot])
plots.export_to_xml()

# -----------------------------
# 7. Запуск только в plotting mode
# -----------------------------
openmc.plot_geometry()