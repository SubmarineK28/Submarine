import openmc
import openmc.deplete
import numpy as np
import math
from math import pi
import matplotlib
from matplotlib import pyplot
import matplotlib.pyplot as plt
plt.rcParams['lines.markersize'] = 3

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

fuel = openmc.Material.mix_materials([Gd2O3, UO2], [0.00025, 0.99975], 'wo')
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

# Rectangular_prism

#pitch = 1.26
#fuel.volume = pi*(d-2*delta)**2/4 # cm**2
#UO2.volume = pi*(d-2*delta)**2/4 # cm**2
#w_by_vol = 110.0 # w/cm**3
#s = pitch**2 # cm**2
#w = w_by_vol*s # w/cm
#left = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
#right = openmc.XPlane(x0=pitch/2, boundary_type='reflective')
#bottom = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
#top = openmc.YPlane(y0=pitch/2, boundary_type='reflective')
#water_region = +left & -right & +bottom & -top & +clad_or
#moderator = openmc.Cell(4, 'moderator')
#moderator.fill = water
#moderator.region = water_region

# Rectangular_prism

# Hexagonal_prism

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

#openmc.plot_geometry()

#openmc.run()

model = openmc.Model(geometry=geom, settings=settings, materials=mats)
#chain_file = 'chain_simple.xml'
chain_file = '/home/adminsrv/projects/Vorobiova_K/Lab_8/VVER/gd/chain_casl_pwr.xml'
op = openmc.deplete.CoupledOperator(model, chain_file)
time_steps = [0.01] + [1.0] + [5.0] + [50.0]*20
power = w
#time_steps = [0.1] + [9.9] + [10.0]*7
#time_steps = [0.25]*20 + [5.0] + [10.0]*59
#power = [1.0*w]*80
integrator = openmc.deplete.PredictorIntegrator(op, time_steps, power, timestep_units='d')
integrator.integrate()

results = openmc.deplete.Results("./depletion_results.h5")
time, k = results.get_keff()
time_days = time / 86400
power_MW = w * 1e-6
fuel_mass = fuel.volume * ro_UO2 * 1e-3
burnup = power_MW * time_days / fuel_mass

# ============================== Линейная аппроксимация ==============================
time_steps = np.array(time_steps)
start = np.where(np.isclose(time_steps, 50.0))[0][0]
burnup_work = burnup[start:]
k_work = k[start:]
weights = 1.0 / k_work[:,1]
coeffs = np.polyfit(burnup_work, k_work[:,0], deg=1, w=weights)
a, b = coeffs
print(f" ================= k = {a:.5f} * burnup + {b:.5f} =================")
# =========================== Квадратичная аппроксимация =============================
coeffs = np.polyfit(burnup_work, k_work[:,0], deg=2, w=weights)
A = coeffs[0]
B = coeffs[2]
print(f" ================= k = {A:.5f} * burnup^2 + {B:.5f} =================")
# ============================== Численное решение линейное ==========================
x0 = 1.0
def k_eff(x):
    return a*x+b
def eq_1(x):
    return k_eff(x) - 1
from scipy.optimize import fsolve
x_1 = fsolve(eq_1, x0)

def eq_3(x):
    return (k_eff(x/3) + k_eff(x/2) + k_eff(x))/3 - 1
from scipy.optimize import fsolve
x_3 = fsolve(eq_3, x0)
# ============================== Численное решение квадратичное ====================
x0 = 1.0
def k_eff(x):
    return A*(x**2)+B
def eq_1(x):
    return k_eff(x) - 1
from scipy.optimize import fsolve
X_1 = fsolve(eq_1, x0)

def eq_3(x):
    return (k_eff(x/3) + k_eff(x/2) + k_eff(x))/3 - 1
from scipy.optimize import fsolve
X_3 = fsolve(eq_3, x0)

# ============================== Теоретическая кривая ==============================
burnup_fit = np.linspace(burnup.min(), burnup.max(), 200)
k_fit = a * burnup_fit + b
k_fit_2 = a * (burnup_fit**2) + b
# Значения в найденных точках
k_x1 = k_eff(x_1[0])
k_x3 = k_eff(x_3[0])

K2_x1 = k_eff(X_1[0])
K2_x3 = k_eff(X_3[0]) # - this

print(f"k(x_1) = {k_x1:.5f} при burnup = {x_1[0]:.5f}")
#print(f"k(x_3) = {k_x3:.5f} при burnup = {x_3[0]:.5f}")
pyplot.errorbar(burnup, k[:, 0], yerr=k[:, 1], fmt='o', linestyle='--', label='Расчётные точки')
pyplot.plot(burnup_fit, k_fit, '-', label='Линейная аппроксимация')
pyplot.plot(burnup_fit, k_fit_2, '-', label='Квадратичная аппроксимация')

pyplot.plot(x_1, k_x1, 'ro', label=f'k(x₁)=1 при BU={x_1[0]:.2f}')
pyplot.plot(X_1, K2_x1, 'ro', label=f'k(x₁)=1 при BU_2={X_1[0]:.2f}')
pyplot.plot(x_3, k_x3, 'go', label=f'Среднее k(x/3,x/2,x)=1 при BU={X_3[0]:.2f}')
pyplot.plot(X_3, K2_x3, 'go', label=f'Среднее k(x/3,x/2,x)=1 при BU_2={X_3[0]:.2f}')

pyplot.xlabel("Burnup [MWd/kg]")
pyplot.ylabel("$k_{eff}$")
pyplot.legend()
pyplot.grid(True)
pyplot.savefig("k_vs_burnup.png", dpi=300, bbox_inches='tight')
pyplot.show()

_, u235 = results.get_atoms("1", "U235")
_, pu239 = results.get_atoms("1", "Pu239")
_, pu240 = results.get_atoms("1", "Pu240")
_, xe135 = results.get_atoms("1", "Xe135")
_, gd157 = results.get_atoms("6", "Gd157")

u235_norm = u235 / u235[0]
pu239_norm = pu239 / pu239[-1]
pu240_norm = pu240 / pu240[-1]
gd157_norm = gd157 / gd157[0]

pyplot.plot(burnup, u235, marker='o', linestyle='-', label="U-235")
pyplot.plot(burnup, pu239, marker='s', linestyle='-', label="Pu-239")
pyplot.plot(burnup, pu240, marker='^', linestyle='-', label="Pu-240")
pyplot.plot(burnup, gd157, marker='o', linestyle='-', label="Gd-157")

with open("depletion_printout.txt", "w") as f:
    f.write(f"k(x_1) = {k_x1:.5f} при burnup = {x_1[0]:.5f}\n")
    f.write(f"k(x_3) = {k_x3:.5f} при burnup = {x_3[0]:.5f}\n\n")

    f.write("k (mean, sigma) = \n" + np.array2string(k, precision=6) + "\n\n")
    f.write("U235 = \n" + np.array2string(u235, precision=6) + "\n\n")
    f.write("Pu239 = \n" + np.array2string(pu239, precision=6) + "\n\n")
    f.write("Pu240 = \n" + np.array2string(pu240, precision=6) + "\n\n")
    f.write("Xe135 = \n" + np.array2string(xe135, precision=6) + "\n")
    f.write("Gd157 = \n" + np.array2string(gd157, precision=6) + "\n")
    f.write("burnup = \n" + np.array2string(burnup, precision=3) + "\n")


print(k)
print(u235)
print(pu239)
print(pu240)
print(xe135)
print(gd157)