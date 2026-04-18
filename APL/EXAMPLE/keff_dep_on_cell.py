import openmc
import numpy as np
import sys
import os
import subprocess
import glob

# === Получаем номер запуска ===
if len(sys.argv) > 1:
    N = int(sys.argv[1])
else:
    N = 0

print(f"▶ Запуск расчёта для N={N}")

# === Геометрические поправки ===
dx = 0.025 * N

x = 0.036
#x = 0.044

# Density g/cm3
ro_UO2 = 10.0
ro_H2O = 0.7
ro_He = 0.15
ro_Zr = 6.6

# Materials
den_U235 = x
den_U238 = 1-x

fuel = openmc.Material(1, "fuel")
fuel.add_nuclide('U235', den_U235, 'ao')
fuel.add_nuclide('U238', den_U238, 'ao')
fuel.add_nuclide('O16', 2.0, 'ao')
fuel.set_density('g/cm3', ro_UO2)
fuel.temperature = 1200.0

water = openmc.Material(2, "water")
water.add_nuclide('H1', 2.0, 'ao')
water.add_nuclide('O16', 1.0, 'ao')
water.add_s_alpha_beta('c_H_in_H2O')
water.set_density('g/cm3', ro_H2O)
water.temperature = 600.0

helium = openmc.Material(3, "helium")
helium.add_element('He', 1.0, 'ao')
helium.set_density('g/cm3', ro_He)
helium.temperature = 900.0

zirconium = openmc.Material(4, "zirconium")
zirconium.add_element('Zr', 1.0, 'ao')
zirconium.set_density('g/cm3', ro_Zr)
zirconium.temperature = 900.0

mats = openmc.Materials()
mats = openmc.Materials([fuel, water, helium, zirconium])
mats.cross_sections = "/home/adminsrv/projects/sections/endfb-viii.0-hdf5/cross_sections.xml"

mats.export_to_xml()

# Geometry

d = 0.60
delta = 0.075
gap = 0.01

fuel_or = openmc.ZCylinder(r=d/2-delta-gap)
clad_ir = openmc.ZCylinder(r=d/2-delta)
clad_or = openmc.ZCylinder(r=d/2)

fuel_region = -fuel_or
gap_region = +fuel_or & -clad_ir
clad_region = +clad_ir & -clad_or

fuel_c = openmc.Cell(1, 'fuel')
fuel_c.fill = fuel
fuel_c.region = fuel_region

gap = openmc.Cell(2, 'gap')
gap.fill = helium
gap.region = gap_region

clad = openmc.Cell(3, 'clad')
clad.fill = zirconium
clad.region = clad_region

# Rectangular_prism

pitch = 0.7 + dx
left = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
right = openmc.XPlane(x0=pitch/2, boundary_type='reflective')
bottom = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
top = openmc.YPlane(y0=pitch/2, boundary_type='reflective')
water_region = +left & -right & +bottom & -top & +clad_or
moderator = openmc.Cell(4, 'moderator')
moderator.fill = water
moderator.region = water_region

# Rectangular_prism

# Hexagonal_prism

#pitch = 1.4
#side = pitch/(3**(1/2))
#hex_region = openmc.model.hexagonal_prism(edge_length=side, orientation='x', boundary_type='reflective')
#out_region = +clad_or
#water_region = hex_region & out_region
#moderator = openmc.Cell(4, 'moderator')
#moderator.fill = water
#moderator.region = water_region

# Hexagonal_prism

root = openmc.Universe(cells=(fuel_c, gap, clad, moderator))

geom = openmc.Geometry()
geom.root_universe = root
geom.export_to_xml()

point = openmc.stats.Point((0, 0, 0))
src = openmc.Source(space=point)

settings = openmc.Settings()
settings.source = src
settings.batches = 1000 # 1000
settings.inactive = 150
settings.particles = 4000 # 5000 - 10000
settings.temperature = {'multipole': True, 'method': 'interpolation', 'range': [290, 2500]}

settings.export_to_xml()

plot = openmc.Plot()
plot.filename = 'pinplot'
plot.basis = 'xy'
plot.origin = (0.0, 0.0, 0.0)
plot.width = (pitch+1, pitch+1)
plot.pixels = (600, 600)
#plot.color_by = 'cell'
plot.color_by = 'material'
plot.colors = {fuel: 'black', zirconium: 'grey', water: 'blue'}
plots = openmc.Plots([plot])

plots.export_to_xml()

openmc.plot_geometry()

# Instantiate an empty Tallies object
tallies_file = openmc.Tallies()

# Instantiate energy filter for multi-group cross-section Tallies
energy_filter = openmc.EnergyFilter([0., 0.625, 20.0e6])

# Instantiate flux Tally in moderator and fuel
tally = openmc.Tally(name='flux')
tally.filters = [openmc.CellFilter([fuel_c, moderator])]
tally.filters.append(energy_filter)
tally.scores = ['flux']
tallies_file.append(tally)

# Instantiate reaction rate Tally in fuel
tally = openmc.Tally(name='fuel rxn rates')
tally.filters = [openmc.CellFilter(fuel_c)]
tally.filters.append(energy_filter)
tally.scores = ['nu-fission', 'scatter']
tally.nuclides = ['U238', 'U235']
tallies_file.append(tally)

# Instantiate reaction rate Tally in moderator
tally = openmc.Tally(name='moderator rxn rates')
tally.filters = [openmc.CellFilter(moderator)]
tally.filters.append(energy_filter)
tally.scores = ['absorption', 'total']
tally.nuclides = ['O16', 'H1']
tallies_file.append(tally)

# Instantiate a tally mesh
#mesh = openmc.RegularMesh(mesh_id=1)
#mesh.dimension = [1, 1, 1]
#mesh.lower_left = [-0.63, -0.63, -100.]
#mesh.width = [1.26, 1.26, 200.]
#meshsurface_filter = openmc.MeshSurfaceFilter(mesh)
#
## Instantiate thermal, fast, and total leakage tallies
#leak = openmc.Tally(name='leakage')
#leak.filters = [meshsurface_filter]
#leak.scores = ['current']
#tallies_file.append(leak)

#thermal_leak = openmc.Tally(name='thermal leakage')
#thermal_leak.filters = [meshsurface_filter, openmc.EnergyFilter([0., 0.625])]
#thermal_leak.scores = ['current']
#tallies_file.append(thermal_leak)
#
#fast_leak = openmc.Tally(name='fast leakage')
#fast_leak.filters = [meshsurface_filter, openmc.EnergyFilter([0.625, 20.0e6])]
#fast_leak.scores = ['current']
#tallies_file.append(fast_leak)

# K-Eigenvalue (infinity) tallies
fiss_rate = openmc.Tally(name='fiss. rate')
abs_rate = openmc.Tally(name='abs. rate')
fiss_rate.scores = ['nu-fission']
abs_rate.scores = ['absorption']
tallies_file += (fiss_rate, abs_rate)

# Resonance Escape Probability tallies
therm_abs_rate = openmc.Tally(name='therm. abs. rate')
therm_abs_rate.scores = ['absorption']
therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625])]
tallies_file.append(therm_abs_rate)

# Thermal Flux Utilization tallies
fuel_therm_abs_rate = openmc.Tally(name='fuel therm. abs. rate')
fuel_therm_abs_rate.scores = ['absorption']
fuel_therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625]),
                               openmc.CellFilter([fuel_c])]
tallies_file.append(fuel_therm_abs_rate)

# Fast Fission Factor tallies
therm_fiss_rate = openmc.Tally(name='therm. fiss. rate')
therm_fiss_rate.scores = ['nu-fission']
therm_fiss_rate.filters = [openmc.EnergyFilter([0., 0.625])]
tallies_file.append(therm_fiss_rate)

# Instantiate energy filter to illustrate Tally slicing
fine_energy_filter = openmc.EnergyFilter(np.logspace(np.log10(1e-2), np.log10(20.0e6), 10))

# Instantiate flux Tally in moderator and fuel
tally = openmc.Tally(name='need-to-slice')
tally.filters = [openmc.CellFilter([fuel_c, moderator])]
tally.filters.append(fine_energy_filter)
tally.scores = ['nu-fission', 'scatter']
tally.nuclides = ['H1', 'U238']
tallies_file.append(tally)


# Export to "tallies.xml"
tallies_file.export_to_xml()

openmc.run()

# === Запуск анализа ===
analyze_script = os.path.join(os.path.dirname(__file__), "Tally_Data_Processing.py")

if os.path.exists(analyze_script):
    print(f" Запуск анализа результатов для N={N}")
    subprocess.run([sys.executable, analyze_script, str(N)])
else:
    print(f" Не найден {analyze_script}")