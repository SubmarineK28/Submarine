import os
import openmc
import numpy as np
import math

import sys

sys.path.append('../') # импорт путей из папки сверху
from chanel_3d import get_plane, set_chanel_3d
from fuel_assembly import assembly_layer_construct
from materials_met_1000 import fuel_type1inner
from materials_met_1000 import fuel_type2inner
from materials_met_1000 import fuel_type3inner
from materials_met_1000 import fuel_type4inner
from materials_met_1000 import fuel_type5inner
from materials_met_1000 import fuel_type6inner


from materials_met_1000 import (
    coolant,
    steel,
    clading,
    helium,
    absorber_enriched,
    absorber_burnup,
)


def fa_3d(tvs_ind=1000, fuels_lst = [fuel_type1inner]*6, pel_list = [coolant]*0 + [absorber_enriched]*6, boundary="transmission/reflective", is_void= False):
    # z boundaries - определяем условия отражения для высотных слоев
    # следует учесть, что источник устанавливается в точку (0,0) и,
    # поэтому слои должны охватывать источник сверху и снизу
    z0 = 0
    z1 = 12 # steel
    z2 = z1 + 12 # coolent
    ztype1 = z2 + 20 # fuel_type1inner
    ztype2 = z2 + 40 #fuel_type2inner
    ztype3 = z2 + 60 #fuel_type3inner
    ztype4 = z2 + 80 #fuel_type4inner
    ztype5 = z2 + 100 #fuel_type5inner
    ztype6 = z2 + 120 #fuel_type6inner
    z3 =  ztype6 + 12 # coolent
    z4 = z3 +  12 # steel

    core_down = openmc.ZPlane(z0= z0 , boundary_type="vacuum")
    z_1 = openmc.ZPlane(z0= z1 , boundary_type="transmission")

    z_2 = openmc.ZPlane(z0= z2 , boundary_type="transmission")
    z_ft1 = openmc.ZPlane(z0= ztype1 , boundary_type="transmission")

    z_ft2 = openmc.ZPlane(z0= ztype2 , boundary_type="transmission")
    z_ft3 = openmc.ZPlane(z0= ztype3 , boundary_type="transmission")

    z_ft4 = openmc.ZPlane(z0= ztype4 , boundary_type="transmission")
    z_ft5 = openmc.ZPlane(z0= ztype5 , boundary_type="transmission")
    z_ft6 = openmc.ZPlane(z0= ztype6 , boundary_type="transmission")

    z_bs = openmc.ZPlane(z0= z3 , boundary_type="transmission")
    core_up = openmc.ZPlane(z0= z4 , boundary_type="vacuum")

    coolant_inside = coolant


    #  создадим послойно ячейки с разными материалами
    # !!! приходиться создать дополнительные ячейки и Universe, чтобы можно было получить
    # !!! не только срдение значения по всей ТВС, но и средние значения по слоям

    cell_name = "fuel assembly low steel"
    root_cell_1 = assembly_layer_construct( cell_name=cell_name, fuel=steel, clading=steel, coolant=steel, absorber_enriched=steel, absorber_burnup=steel, zup=z_1, zdn=core_down, boundary=boundary)
    universe1 = openmc.Universe( name=f'universe_{tvs_ind}_{cell_name}')
    universe1.add_cell(root_cell_1)
    view_cell_1 = openmc.Cell(name="view_low_steel", fill=universe1)
    view_cell_1.region = root_cell_1.region


    cell_name = "fuel assembly lower coolant"
    root_cell_2 = assembly_layer_construct(cell_name=cell_name,fuel=clading, clading=clading, coolant=coolant_inside, absorber_enriched=clading, absorber_burnup=clading, zup=z_2, zdn=z_1, boundary=boundary)
    universe2 = openmc.Universe( name=f'universe_{tvs_ind}_{cell_name}')
    universe2.add_cell(root_cell_2)
    view_cell_2 = openmc.Cell( name="view_lower_coolant", fill=universe2)
    view_cell_2.region = root_cell_2.region


    cell_name = "fuel assembly inner fuel 1"
    root_cell_3 = assembly_layer_construct( cell_name=cell_name,fuel=fuels_lst[0], clading=clading, coolant=coolant_inside, absorber_enriched = pel_list[0], zup=z_ft1, zdn=z_2, boundary=boundary)
    universe3 = openmc.Universe(name=f'universe_{tvs_ind}_{cell_name}')
    universe3.add_cell(root_cell_3)
    view_cell_3 = openmc.Cell( name="view_fuel_1", fill=universe3)
    view_cell_3.region = root_cell_3.region


    cell_name = "fuel assembly inner fuel 2"
    root_cell_4 = assembly_layer_construct(cell_name=cell_name,fuel=fuels_lst[1], clading=clading, coolant=coolant_inside, absorber_enriched = pel_list[1], zup=z_ft2, zdn=z_ft1, boundary=boundary)
    universe4 = openmc.Universe(name=f'universe_{tvs_ind}_{cell_name}')
    universe4.add_cell(root_cell_4)
    view_cell_4 = openmc.Cell( name="view_fuel_2", fill=universe4)
    view_cell_4.region = root_cell_4.region


    cell_name = "fuel assembly inner fuel 3"
    root_cell_5 = assembly_layer_construct( cell_name=cell_name,fuel=fuels_lst[2], clading=clading, coolant=coolant_inside, absorber_enriched = pel_list[2], zup=z_ft3, zdn=z_ft2, boundary=boundary)
    universe5 = openmc.Universe(name=f'universe_{tvs_ind}_{cell_name}')
    universe5.add_cell(root_cell_5)
    view_cell_5 = openmc.Cell( name="view_fuel_3", fill=universe5)
    view_cell_5.region = root_cell_5.region


    cell_name = "fuel assembly inner fuel 4"
    root_cell_6 = assembly_layer_construct( cell_name=cell_name,fuel=fuels_lst[3], clading=clading, coolant=coolant_inside, absorber_enriched = pel_list[3], zup=z_ft4, zdn=z_ft3, boundary=boundary)
    universe6 = openmc.Universe(name=f'universe_{tvs_ind}_{cell_name}')
    universe6.add_cell(root_cell_6)
    view_cell_6 = openmc.Cell(name="view_fuel_4", fill=universe6)
    view_cell_6.region = root_cell_6.region


    cell_name = "fuel assembly inner fuel 5"
    root_cell_7 = assembly_layer_construct( cell_name=cell_name,fuel=fuels_lst[4], clading=clading, coolant=coolant_inside, absorber_enriched = pel_list[4], zup=z_ft5, zdn=z_ft4, boundary=boundary)
    universe7 = openmc.Universe(name=f'universe_{tvs_ind}_{cell_name}')
    universe7.add_cell(root_cell_7)
    view_cell_7 = openmc.Cell( name="view_fuel_5", fill=universe7)
    view_cell_7.region = root_cell_7.region


    cell_name = "fuel assembly inner fuel 6"
    root_cell_8 = assembly_layer_construct( cell_name=cell_name,fuel=fuels_lst[5], clading=clading, coolant=coolant_inside, absorber_enriched = pel_list[5], zup=z_ft6, zdn=z_ft5, boundary=boundary)
    universe8 = openmc.Universe(name=f'universe_{tvs_ind}_{cell_name}')
    universe8.add_cell(root_cell_8)
    view_cell_8 = openmc.Cell( name="view_fuel_6", fill=universe8)
    view_cell_8.region = root_cell_8.region


    cell_name = "fuel assembly upper coolant"
    root_cell_9 = assembly_layer_construct( cell_name=cell_name,fuel=clading, clading=clading, coolant=coolant_inside, absorber_enriched=clading, absorber_burnup=clading, zup=z_bs, zdn=z_ft6, boundary=boundary)
    universe9 = openmc.Universe(name=f'universe_{tvs_ind}_{cell_name}')
    universe9.add_cell(root_cell_9)
    view_cell_9 = openmc.Cell( name="view_upper_coolant", fill=universe9)
    view_cell_9.region = root_cell_9.region


    cell_name = "fuel assembly upper steel"
    root_cell_10 = assembly_layer_construct( cell_name=cell_name,fuel=steel, clading=steel, coolant=steel, absorber_enriched=steel, absorber_burnup=steel, zup=core_up , zdn=z_bs, boundary=boundary)
    universe10 = openmc.Universe(name=f'universe_{tvs_ind}_{cell_name}')
    universe10.add_cell(root_cell_10)
    view_cell_10 = openmc.Cell( name="view_upper_steel", fill=universe10 )
    view_cell_10.region = root_cell_10.region



    # здадим объединенную Universe, которую и будем рассчитывать
    fa_inner_universe_return = openmc.Universe( name='fa_inner_universe')
    fa_inner_universe_return.add_cell(view_cell_1)
    fa_inner_universe_return.add_cell(view_cell_2)
    fa_inner_universe_return.add_cell(view_cell_3)
    fa_inner_universe_return.add_cell(view_cell_4)
    fa_inner_universe_return.add_cell(view_cell_5)
    fa_inner_universe_return.add_cell(view_cell_6)
    fa_inner_universe_return.add_cell(view_cell_7)
    fa_inner_universe_return.add_cell(view_cell_8)
    fa_inner_universe_return.add_cell(view_cell_9)
    fa_inner_universe_return.add_cell(view_cell_10)

    return fa_inner_universe_return
if __name__ == '__main__':
    openmc.config['cross_sections'] = "/home/sparrow/APL/sections/endfb-viii.0-hdf5/cross_sections.xml"

    mats = openmc.Materials((fuel_type1inner, fuel_type2inner, fuel_type3inner, fuel_type4inner, fuel_type5inner, fuel_type6inner,
                             clading, coolant, steel, helium, absorber_enriched, absorber_burnup))
    mats.export_to_xml()

    fuels_inner_lst = [fuel_type1inner, fuel_type2inner, fuel_type3inner, fuel_type4inner, fuel_type5inner, fuel_type6inner]

    is_void = False
    fa_inner_universe = fa_3d(tvs_ind=1000, fuels_lst=fuels_inner_lst, boundary="reflective", is_void=is_void)

    # вводим геометрию корневого объекта
    geom = openmc.Geometry(fa_inner_universe)
    geom.export_to_xml()

    z0 = 0
    z1 = 12 # steel
    z2 = z1 + 12 # coolent
    ztype1 = z2 + 20 # fuel_type1inner
    ztype2 = z2 + 40 #fuel_type2inner
    ztype3 = z2 + 60 #fuel_type3inner
    ztype4 = z2 + 80 #fuel_type4inner
    ztype5 = z2 + 100 #fuel_type5inner
    ztype6 = z2 + 120 #fuel_type6inner
    z3 =  ztype6 + 12 # coolent
    z4 = z3 +  12 # steel

    # печать картинки
    p = openmc.Plot()
    # возможность задания срезов для отображения
    p.origin = (0, 0, z2 + 1)
    p.filename = f'fuel_assembly_xy_is_void_{is_void}'
    p.basis = "xy"
    p.width = (26, 26)
    p.pixels = (1000, 1000)
    p.color_by = 'material'
    p.colors = {
        fuel_type1inner: 'red',
        fuel_type2inner: 'orange',
        fuel_type3inner: 'magenta',
        fuel_type4inner: 'coral',
        fuel_type5inner: 'pink',
        fuel_type6inner: 'salmon',
        coolant: 'blue',
        clading: 'gray',
        steel: 'darkgray',
        helium: 'lime',
        absorber_enriched: 'green',
        absorber_burnup: 'yellow',
    }

    plots = openmc.Plots([p])
    plots.export_to_xml()
    openmc.plot_geometry()

    # возможность задания срезов для отображения
    p.origin = (0.0, 0.0, (z0 + z4) / 2.0)
    p.filename = f'fuel_assembly_yz_is_void_{is_void}'
    p.basis = "yz"
    p.width = (26, z4 - z0)
    p.pixels = (1000, 1000)
    p.color_by = 'material'
    p.colors = {
        fuel_type1inner: 'red',
        fuel_type2inner: 'orange',
        fuel_type3inner: 'magenta',
        fuel_type4inner: 'coral',
        fuel_type5inner: 'pink',
        fuel_type6inner: 'salmon',
        coolant: 'blue',
        clading: 'gray',
        steel: 'darkgray',
        helium: 'lime',
        absorber_enriched: 'green',
        absorber_burnup: 'yellow',
    }

    plots = openmc.Plots([p])
    plots.export_to_xml()
    openmc.plot_geometry()

    # Computing settings
    batches = 5
    inactive = 1
    particles = 10

    settings = openmc.Settings()
    settings.batches = batches
    settings.inactive = inactive
    settings.particles = particles
    settings.output = {'tallies': True}
    # установка точечного источника (могут быть разные)
    source_point = openmc.stats.Point(xyz=(0, 0, ztype3 + 1))
    settings.source = openmc.Source(space=source_point)
    settings.temperature = {"method": "interpolation"}

    settings.export_to_xml()

    openmc.run()
