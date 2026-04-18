import openmc
import openmc.deplete
import numpy as np
import math
from math import pi
import matplotlib
from matplotlib import pyplot

x = 0.03

# Density g/cm3
ro_UO2 = 10.0
ro_H2O = 0.7
ro_Zr = 6.6
ro_Gd2O3 = 7.4
ro_Er2O3 = 8.6

# Materials
den_U235 = x
den_U238 = 1-x

UO2 = openmc.Material(name='UO2')
UO2.add_nuclide('U235', den_U235, 'ao')
UO2.add_nuclide('U238', den_U238, 'ao')
UO2.add_nuclide('O16', 2.0, 'ao')
UO2.set_density('g/cm3', ro_UO2)
UO2.temperature = 1200.0

water = openmc.Material(name='water')
water.add_nuclide('H1', 2.0, 'ao')
water.add_nuclide('O16', 1.0, 'ao')
water.add_s_alpha_beta('c_H_in_H2O')
water.set_density('g/cm3', ro_H2O)
water.temperature = 600.0

zirconium = openmc.Material(name='zirconium')
zirconium.add_element('Zr', 1.0, 'ao')
zirconium.set_density('g/cm3', ro_Zr)
zirconium.temperature = 600.0

Gd2O3=openmc.Material(name='Gd2O3')
Gd2O3.add_element('Gd', 2.0, 'ao')
Gd2O3.add_element('O', 3.0, 'ao')
Gd2O3.set_density('g/cm3', ro_Gd2O3)
Gd2O3.temperature = 1200.0 

Er2O3=openmc.Material(name='Er2O3')
Er2O3.add_element('Er', 2.0, 'ao')
Er2O3.add_element('O', 3.0, 'ao')
Er2O3.set_density('g/cm3', ro_Er2O3)
Er2O3.temperature = 1200.0

fuel = openmc.Material.mix_materials([Gd2O3, UO2], [0.00025, 0.99975], 'wo') # !!!!!!!!!!!!!!!!!!!!
#fuel = openmc.Material.mix_materials([Er2O3, UO2], [0.02, 0.98], 'wo')
fuel.temperature = 1200.0

mats = openmc.Materials()
mats = openmc.Materials([fuel, UO2, Gd2O3, Er2O3, water, zirconium])
mats.cross_sections = "/home/adminsrv/projects/sections/endfb-viii.0-hdf5/cross_sections.xml"

mats.export_to_xml()

# Geometry

d = 0.60
delta = 0.04

fuel_or = openmc.ZCylinder(r=d/2-delta)
clad_or = openmc.ZCylinder(r=d/2)

fuel_region = -fuel_or
clad_region = +fuel_or & -clad_or

fuel_c = openmc.Cell(1, 'fuel')
fuel_c.fill = fuel
fuel_c.region = fuel_region

clad = openmc.Cell(3, 'clad')
clad.fill = zirconium
clad.region = clad_region

pitch = 1.26
side = pitch/(3**(1/2)) # cm
fuel.volume = pi*(d-2*delta)**2/4 # cm**2
UO2.volume = pi*(d-2*delta)**2/4 # cm**2
w_by_vol = 110.0 # w/cm**3
s = 3*3**(1/2)/2*side**2 # cm**2
w = w_by_vol*s # w/cm
hex_region = openmc.model.hexagonal_prism(edge_length=side, orientation='x', boundary_type='reflective')
out_region = +clad_or
water_region = hex_region & out_region
moderator = openmc.Cell(4, 'moderator')
moderator.fill = water
moderator.region = water_region

# Hexagonal_prism

root = openmc.Universe(cells=(fuel_c, clad, moderator))

geom = openmc.Geometry()
geom.root_universe = root
geom.export_to_xml()

point = openmc.stats.Point((0, 0, 0))
src = openmc.Source(space=point)

settings = openmc.Settings()
settings.source = src
settings.batches = 1000
settings.inactive = 10
settings.particles = 1000

settings.export_to_xml()

openmc.run()