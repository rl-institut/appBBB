#!/usr/bin/python
# -*- coding: utf-8 -*-

import pandas as pd
import os


def calc_dh_supply_temp(temp, **kwargs):
    '''
    Generates an hourly supply temperature profile depending on the ambient
    temperature.
    For ambient temperatures above T_heat_period the load for tap water
    preparation dominates the heat laod. The district heating system is then
    mass flow controlled and the supply temperature kept at a constant
    temperature of T_supply_min.
    For ambient temperatures below T_heat_period the supply temperature
    increases linearly to T_supply_max with decreasing ambient temperature.
    '''
    T_supply_max = kwargs.get('T_supply_max', 135)  # max supply temperature
    T_supply_min = kwargs.get('T_supply_min', 70)  # min supply temperature
    T_heat_period = kwargs.get('T_heat_period', 15)  # amb. temp. at which heating
                                                     # systems are turned on
    T_amb_min = kwargs.get('T_amb_design', -15)  # amb. temp. where max. supply
                                                 # temp. is reached

    # linear correlation between Q and T_sup
    T_supply = pd.Series(0, index=temp.index)
    slope = (T_supply_min - T_supply_max) / (T_heat_period - T_amb_min)
    y_intercept = T_supply_max - slope * T_amb_min

    T_supply = slope * temp + y_intercept
    T_supply[T_supply < T_supply_min] = T_supply_min
    T_supply[T_supply > T_supply_max] = T_supply_max

    return T_supply


#def add_constraint_dh_heating_storage(om, temp, **kwargs):
#    '''
#    Make constraint for post heating of district heating thermal storage in 
#    case the supply temperatures is higher than the storage temperature.
#    '''
#    T_dh_storage = kwargs.get('T_dh_storage', 99)  # storage temperature
#    T_dh_return = kwargs.get('T_dh_return', 50)  # return temperature of dh
#    # get dh supply temperature
#    T_supply = calc_dh_supply_temp(temp)
#    # returns all district heating storages
#    storages = [obj for obj in om.energysystem.entities
#        if 'dh_thermal_storage' in obj.uid]
#    # write list to hand over to constraint
#    transports_ex = []
#    for export in exports:
#        transports_ex += [(export.uid, export.outputs[0].uid)]
#    # write list to hand over to constraint
#    transports_im = []
#    for imp in imports:
#        transports_im += [(imp.uid, imp.outputs[0].uid)]
#    # add new constraint
#    om.export_minimum_constraint = po.Constraint(expr=(
#        sum(om.w[i, o, t] for i, o in transports_ex for t in om.timesteps) -
#        sum(om.w[i, o, t] for i, o in transports_im for t in om.timesteps)
#        >= float(constraints.query('constr=="export_min"')['val'])))
#        
#    # iterate through thermal storages
#    for storage in storages:
#        for hour in list(range(len(T_supply.index))):
#            if T_supply[hour] > T_dh_storage:
#                prob += (lp_variables['DH Thermal Storage Boiler'][hour]
#                    * T_dh_storage - T_dh_return)
#                    - lp_variables['DH Storage Thermal Discharge'][hour]
#                    * (T_sup[hour] - T_dh_storage)
#                    == 0,
#                    'Provide supply temp ' + str(hour))
#
#            else:
#                prob += (lp_variables['DH Thermal Storage Boiler'][hour]
#                    == 0,
#                    'Provide supply temp ' + str(hour)
#    return
#    
#    
#
#
#
#
#filename = os.path.abspath(os.path.join(os.path.dirname(__file__),
#                                        'temp'))
#temp = pd.read_pickle(filename)
#T_supply = calc_dh_supply_temp(temp)
#for i in list(range(len(T_supply.index))):
#    print(i)
##from matplotlib import pyplot as plt
##T_supply.plot()
##plt.show()