import openmc
import numpy as np

# temperatures for MET-1000
fuel_temperature = 500 + 273.15
coolant_temperature = 300 + 273.15
structure_temperature = 330 + 273.15
r_pin=0.64720/2 # cm what is that ??? 

ro_UO2_Al = 8.87 # = 0.75*10.96 + 0.25*2.6
ro_B4C = 2.5
ro_Gd_Al = 7.618
ro_H2O = 0.7
ro_Zr = 6.6
ro_He = 0.15
ro_Steel = 7.85 

x = 0.2
den_U235 = x
den_U238 = 1-x

#--------------------------------------------------------------------------------
#Add nuclides concentration from inner core boundary distance 68.66 cm to fuel_type1outer  ALL FUEL ASSEMBLY
fuel_type1inner = openmc.Material( name="fuel_type1")
fuel_type1inner.add_nuclide('U235',0.0668, 'ao') # 20 % # 15 %
fuel_type1inner.add_nuclide('U238',0.3783, 'ao') 
fuel_type1inner.add_nuclide('Al27',0.0935, 'ao')
fuel_type1inner.add_nuclide('Si28',0.0165, 'ao')
fuel_type1inner.add_nuclide('O16',0.44500, 'ao')
fuel_type1inner.set_density('g/cm3', ro_UO2_Al)

fuel_type1inner.volume=     np.pi * r_pin ** 2 # Задаем объем единичной высоты. Нужно для моделирования выгорания
fuel_type1inner.depletable = True
fuel_type1inner.temperature=fuel_temperature

#--------------------------------------------------------------------------------
#Add nuclides concentration from inner core boundary distance 68.66 cm to fuel_type1outer  ALL FUEL ASSEMBLY
fuel_type2inner = openmc.Material( name="fuel_type2")
fuel_type2inner.add_nuclide('U235',0.0579, 'ao') # 25 % # 13 %
fuel_type2inner.add_nuclide('U238',0.3783, 'ao')
fuel_type2inner.add_nuclide('Al27',0.0935, 'ao')
fuel_type1inner.add_nuclide('Si28',0.0165, 'ao')
fuel_type2inner.add_nuclide('O16',0.44500, 'ao')
fuel_type2inner.set_density('g/cm3', ro_UO2_Al)

fuel_type2inner.volume=     np.pi * r_pin ** 2 # Задаем объем единичной высоты. Нужно для моделирования выгорания
fuel_type2inner.depletable = True
fuel_type2inner.temperature=fuel_temperature
#--------------------------------------------------------------------------------
#Add nuclides concentration from inner core boundary distance 51.49 cm to fuel_type1outer  ALL FUEL ASSEMBLY
fuel_type3inner = openmc.Material( name="fuel_type3")
fuel_type3inner.add_nuclide('U235',0.1558, 'ao') # 35 %
fuel_type3inner.add_nuclide('U238',0.2893, 'ao')
fuel_type3inner.add_nuclide('Al27',0.0935, 'ao')
fuel_type1inner.add_nuclide('Si28',0.0165, 'ao')
fuel_type3inner.add_nuclide('O16',0.44500, 'ao')
fuel_type3inner.set_density('g/cm3', ro_UO2_Al)

fuel_type3inner.volume=     np.pi * r_pin ** 2 # Задаем объем единичной высоты. Нужно для моделирования выгорания
fuel_type3inner.depletable = True
fuel_type3inner.temperature=fuel_temperature

#--------------------------------------------------------------------------------
#Add nuclides concentration from inner core boundary distance 68.66 cm to fuel_type4inner  ALL FUEL ASSEMBLY
fuel_type4inner = openmc.Material( name="fuel_type4")
fuel_type4inner.add_nuclide('U235',0.1780, 'ao') # 40 %
fuel_type4inner.add_nuclide('U238',0.2670, 'ao')
fuel_type4inner.add_nuclide('Al27',0.0935, 'ao')
fuel_type1inner.add_nuclide('Si28',0.0165, 'ao')
fuel_type4inner.add_nuclide('O16',0.44500, 'ao')
fuel_type4inner.set_density('g/cm3', ro_UO2_Al)

fuel_type4inner.volume=     np.pi * r_pin ** 2 # Задаем объем единичной высоты. Нужно для моделирования выгорания
fuel_type4inner.depletable = True
fuel_type4inner.temperature=fuel_temperature

#--------------------------------------------------------------------------------
#Add nuclides concentration from inner core boundary distance 68.66 cm to fuel_type5inner  ALL FUEL ASSEMBLY
fuel_type5inner = openmc.Material( name="fuel_type5")
fuel_type5inner.add_nuclide('U235',0.0890, 'ao')
fuel_type5inner.add_nuclide('U238',0.3560, 'ao')
fuel_type5inner.add_nuclide('Al27',0.0935, 'ao')
fuel_type1inner.add_nuclide('Si28',0.0165, 'ao')
fuel_type5inner.add_nuclide('O16',0.44500, 'ao')
fuel_type5inner.set_density('g/cm3', ro_UO2_Al)

fuel_type5inner.volume=     np.pi * r_pin ** 2 # Задаем объем единичной высоты. Нужно для моделирования выгорания
fuel_type5inner.depletable = True
fuel_type5inner.temperature=fuel_temperature

#--------------------------------------------------------------------------------
#Add nuclides concentration from inner core boundary distance 68.66 cm to fuel_type5inner  ALL FUEL ASSEMBLY
fuel_type6inner = openmc.Material( name="fuel_type6")
fuel_type6inner.add_nuclide('U235',0.0890, 'ao')
fuel_type6inner.add_nuclide('U238',0.3560, 'ao')
fuel_type6inner.add_nuclide('Al27',0.0935, 'ao')
fuel_type1inner.add_nuclide('Si28',0.0165, 'ao')
fuel_type6inner.add_nuclide('O16',0.44500, 'ao')
fuel_type6inner.set_density('g/cm3', ro_UO2_Al)

fuel_type6inner.volume=     np.pi * r_pin ** 2 # Задаем объем единичной высоты. Нужно для моделирования выгорания
fuel_type6inner.depletable = True
fuel_type6inner.temperature=fuel_temperature

#--------------------------------------------------------------------------------
#absorber enriched
absorber_enriched = openmc.Material( name="absorber_enriched")
#absorber_enriched.add_nuclide('C12',2.040298E-02)
#absorber_enriched.add_nuclide('C13',2.290152E-04)
absorber_enriched.add_nuclide('C12',0.2 , 'ao')
absorber_enriched.add_nuclide('B10',0.4, 'ao')
absorber_enriched.add_nuclide('B11',0.4, 'ao')
absorber_enriched.set_density('g/cm3', ro_B4C)

absorber_enriched.temperature=structure_temperature

#--------------------------------------------------------------------------------
#absorber burnup
absorber_burnup = openmc.Material( name="absorber_burnup")
#absorber_enriched.add_nuclide('C12',2.040298E-02)
#absorber_enriched.add_nuclide('C13',2.290152E-04)
absorber_burnup.add_nuclide('Gd155',0.4 , 'ao')
absorber_burnup.add_nuclide('Al27',0 ,'ao')
absorber_burnup.add_nuclide('O16',0.6, 'ao')
absorber_burnup.set_density('g/cm3', ro_Gd_Al)

absorber_burnup.temperature=structure_temperature

#--------------------------------------------------------------------------------
# Water
coolant = openmc.Material(name='water')
coolant.add_nuclide('H1', 2.0, 'ao')
coolant.add_nuclide('O16', 1.0, 'ao')
coolant.set_density('g/cm3', ro_H2O)
coolant.temperature=coolant_temperature

#--------------------------------------------------------------------------------
# cladding
clading = openmc.Material(name='Clading')
clading.add_element('Zr', 1.0, 'ao')
clading.set_density('g/cm3', ro_Zr)
clading.temperature=structure_temperature

#--------------------------------------------------------------------------------
# helium
helium = openmc.Material(name = "Helium")
helium.add_element('He', 1.0, 'ao')
helium.set_density('g/cm3', ro_He)
helium.temperature = 900.0

#--------------------------------------------------------------------------------
# steel
steel = openmc.Material(name = 'steel')
steel.add_nuclide('Fe56', 0.98, 'ao')
steel.add_nuclide('Cr52', 0.02, 'ao')
steel.set_density('g/cm3', ro_Steel)
steel.temperature = coolant_temperature





