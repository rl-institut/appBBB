#!/usr/bin/python
# -*- coding: utf-8
import numpy as np
import pandas as pd

import helper_SmEnOs as hls
from oemof.core.network.entities import Bus
from oemof.core.network.entities.components import sinks as sink
from oemof.core.network.entities.components import transformers as transformer


def calc_brine_supply_temp():
    """
    Returns the supply temperature for a brine heat pump.
    """
    t_brine = np.zeros((8760,))
    for i in range(8760):
        t_brine[i] = (3.42980994954735 * 10 ** -18) * i ** 5 \
                    - (6.28828944308818 * 10 ** -14) * i ** 4 \
                    + (2.44607151047512 * 10 ** -10) * i ** 3 \
                    + (6.25819661168072 * 10 ** -7) * i ** 2 \
                    - 0.0020535109 * i \
                    + 4.1855152734
    return t_brine


def calc_hp_heating_supply_temp(temp, heating_system, **kwargs):
    """
    Generates an hourly supply temperature profile depending on the ambient
    temperature.
    For ambient temperatures below t_heat_period the supply temperature
    increases linearly to t_sup_max with decreasing ambient temperature.

    Parameters
    temp -- pandas Series with ambient temperature
    heating_system -- string specifying the heating system (floor heating or
        radiator)
    """

    t_heat_period = kwargs.get('t_heat_period', 20)  # outdoor temp upto
                                                     # which is heated
    t_amb_min = kwargs.get('t_amb_min', -14)  # design outdoor temp

    # minimum and maximum supply temperature of heating system
    if heating_system == 'floor heating':
        t_sup_max = kwargs.get('t_sup_max', 35)
        t_sup_min = kwargs.get('t_sup_min', 20)
    elif heating_system == 'radiator':
        t_sup_max = kwargs.get('t_sup_max', 55)
        t_sup_min = kwargs.get('t_sup_min', 20)
    else:
        # TODO Raise Warning
        t_sup_max = kwargs.get('t_sup_max', 55)
        t_sup_min = kwargs.get('t_sup_min', 20)

    # calculate parameters for linear correlation for supply temp and
    # ambient temp
    slope = (t_sup_min - t_sup_max) / (t_heat_period - t_amb_min)
    y_intercept = t_sup_max - slope * t_amb_min

    # calculate supply temperature
    t_sup_heating = slope * temp + y_intercept
    t_sup_heating[t_sup_heating < t_sup_min] = t_sup_min
    t_sup_heating[t_sup_heating > t_sup_max] = t_sup_max

    return t_sup_heating


def cop_heating(temp, type_hp, **kwargs):
    """
    Returns the COP of a heat pump for heating.

    Parameters
    temp -- pandas Series with ambient or brine temperature
    type_hp -- string specifying the heat pump type (air or brine)
    """
    # share of heat pumps in new buildings
    share_hp_new_building = kwargs.get('share_hp_new_building', 0.5)
    # share of floor heating in old buildings
    share_fh_old_building = kwargs.get('share_fbh_old_building', 0.25)
    cop_max = kwargs.get('cop_max', 7)
    if type_hp == 'air':
        eta_g = kwargs.get('eta_g', 0.3)  # COP/COP_max
    elif type_hp == 'brine':
        eta_g = kwargs.get('eta_g', 0.4)  # COP/COP_max
    else:
        # TODO Raise Warning
        eta_g = kwargs.get('eta_g', 0.4)  # COP/COP_max
    # get supply temperatures
    t_sup_fh = calc_hp_heating_supply_temp(temp, 'floor heating')
    t_sup_radiator = calc_hp_heating_supply_temp(temp, 'radiator')
    # share of floor heating systems and radiators
    share_fh = (share_hp_new_building + (1 - share_hp_new_building) *
        share_fh_old_building)
    share_rad = (1 - share_hp_new_building) * (1.0 - share_fh_old_building)

    # calculate COP for floor heating and radiators
    cop_hp_heating_fh = eta_g * ((273.15 + t_sup_fh) / (t_sup_fh - temp))
    cop_hp_heating_fh[cop_hp_heating_fh > cop_max] = cop_max
    cop_hp_heating_rad = eta_g * ((273.15 + t_sup_radiator) /
        (t_sup_radiator - temp))
    cop_hp_heating_rad[cop_hp_heating_rad > cop_max] = cop_max
    cop_hp_heating = (share_fh * cop_hp_heating_fh + share_rad *
        cop_hp_heating_rad)

    return cop_hp_heating


def cop_ww(temp, ww_profile_sfh, ww_profile_mfh, **kwargs):
    """
    Returns the COP of a heat pump for warm water

    Parameters
    temp -- pandas Series with temperature profile (ambient or brine temp.)
    ww_profile_sfh -- pandas Dataframe with warm water profile for
        single family houses
    ww_profile_mfh -- pandas Dataframe with warm water profile for
        multi family houses
    """

    t_ww_sfh = kwargs.get('t_ww_sfh', 50)  # warm water temp. in SFH
    t_ww_mfh = kwargs.get('t_ww_mfh', 60)  # warm water temp. in MFH
    cop_max = kwargs.get('cop_max', 7)
    type_hp = kwargs.get('type_hp', 'air')
    if type_hp == 'air':
        eta_g = kwargs.get('eta_g', 0.3)  # COP/COP_max
    elif type_hp == 'brine':
        eta_g = kwargs.get('eta_g', 0.4)  # COP/COP_max
    else:
        # TODO Raise Warning
        eta_g = kwargs.get('eta_g', 0.4)  # COP/COP_max

    # calculate the share of the warm water demand of sfh and mfh for each hour
    share_sfh = ww_profile_sfh.values / (ww_profile_sfh.values +
        ww_profile_mfh.values)

    # calculates mixed WW supply temperature for single and multi family houses
    t_sup_ww = share_sfh * t_ww_sfh + (1 - share_sfh) * t_ww_mfh

    # calculate COP
    cop = eta_g * ((t_sup_ww + 273.15) / (t_sup_ww - temp))
    cop[cop > cop_max] = cop_max

    return cop


def hp_load_profiles(region, year, demand, share_sfh_hp, share_ww, holidays):
    # splitting of heat load profiles in sfh and mfh as well as ww and heating
    time_index = pd.date_range('1/1/{0}'.format(year), periods=8760, freq='H')
    profile_sfh_heating = pd.Series(0, index=time_index)
    profile_mfh_heating = pd.Series(0, index=time_index)
    profile_sfh_ww = pd.Series(0, index=time_index)
    profile_mfh_ww = pd.Series(0, index=time_index)
    # heat pump demand is derived from value for "Sonstige EE" from the energy
    # balance under the assumption that this is mostly "Umgebungswärme" used in
    # heat pumps with 2/3 environmental heat and 1/3 electricity
    demand_hp = demand * 1.5
    if share_sfh_hp != 0:
        profile_sfh_heating = hls.call_heat_demandlib(
            region, time_index, holidays=holidays,
            annual_heat_demand=(
                demand_hp * share_sfh_hp * (1 - share_ww)),
            shlp_type='EFH', ww_incl=False)
        profile_sfh_heating_ww = hls.call_heat_demandlib(
            region, time_index, holidays=holidays,
            annual_heat_demand=demand_hp * share_sfh_hp,
            shlp_type='EFH', ww_incl=True)
        profile_sfh_ww = profile_sfh_heating_ww - profile_sfh_heating

    if share_sfh_hp != 1:
        profile_mfh_heating = hls.call_heat_demandlib(
            region, time_index, holidays=holidays,
            annual_heat_demand=(
                demand_hp * (1 - share_sfh_hp) * (1 - share_ww)),
            shlp_type='MFH', ww_incl=False)
        profile_mfh_heating_ww = hls.call_heat_demandlib(
            region, time_index, holidays=holidays,
            annual_heat_demand=demand_hp * (1 - share_sfh_hp),
            shlp_type='MFH', ww_incl=True)
        profile_mfh_ww = profile_mfh_heating_ww - profile_mfh_heating
    return (profile_sfh_heating, profile_mfh_heating, profile_sfh_ww,
        profile_mfh_ww)


def create_hp_entities(region, year, demand, elec_bus, temp,
    share_sfh_hp, share_ww, share_air_hp, share_heating_rod, share_heat_storage,
    eta_th, eta_in, eta_out, cap_loss, opex_fix, holidays):

    # get profiles
    (profile_sfh_heating, profile_mfh_heating, profile_sfh_ww,
    profile_mfh_ww) = hp_load_profiles(
        region, year, demand, share_sfh_hp, share_ww, holidays)

    # create buses and sinks for each heat pump as well as heating and ww
    if share_air_hp != 0:
        # air heat pump heating
            # bus
        bus_hp = Bus(uid=('bus', region.name, 'residential', 'heat_pump_dec',
                          'air', 'heating'),
            type='heat_pump_dec', price=0, regions=[region], excess=False)
            # sink
        demand = sink.Simple(
            uid=('demand', region.name, 'residential', 'heat_pump_dec', 'air',
                 'heating'),
            inputs=[bus_hp],
            region=region)
        demand.val = (profile_sfh_heating + profile_mfh_heating) * share_air_hp
            # transformer heat pump
        cop = cop_heating(temp, 'air')
        transformer.Simple(
            uid=('transformer', region.name, 'hp', 'air', 'heating'),
            inputs=[elec_bus],
            outputs=[bus_hp],
            out_max=[(max(demand.val) * (1 - share_heating_rod)) * 1.01],
            eta=[cop],
            regions=[region])
        # transformer heating rod
        transformer.Simple(
            uid=('transformer', region.name, 'hp', 'air', 'heating', 'rod'),
            inputs=[elec_bus],
            outputs=[bus_hp],
            in_max=[None],
            out_max=[max(demand.val) * share_heating_rod],
            eta=[eta_th['heat_rod_dec']],
            regions=[region])
        # heat storage
        transformer.Storage(
            uid=('storage', region.name, 'hp', 'air', 'heating'),
            inputs=[bus_hp],
            outputs=[bus_hp],
            cap_max=max(demand.val) * 2 * share_heat_storage,
            out_max=[max(demand.val) * share_heat_storage],
            in_max=[max(demand.val) * share_heat_storage],
            eta_in=eta_in['heat_storage_dec'],
            eta_out=eta_out['heat_storage_dec'],
            cap_loss=cap_loss['heat_storage_dec'],
            opex_fix=opex_fix['heat_storage_dec'])

        # air heat pump warm water
            # bus
        if share_ww != 0:
            bus_hp = Bus(uid=('bus', region.name, 'residential', 'heat_pump_dec',
                          'air', 'ww'),
            type='heat_pump_dec', price=0, regions=[region], excess=False)
                # sink
            demand = sink.Simple(
                uid=('demand', region.name, 'residential', 'heat_pump_dec', 'air',
                     'ww'),
                inputs=[bus_hp],
                region=region)
            demand.val = (profile_sfh_ww + profile_mfh_ww) * share_air_hp
                # transformer heat pump
            cop = cop_ww(temp, profile_sfh_ww, profile_mfh_ww)
            transformer.Simple(
                uid=('transformer', region.name, 'hp', 'air', 'ww'),
                inputs=[elec_bus],
                outputs=[bus_hp],
                out_max=[(max(demand.val) * (1 - share_heating_rod)) * 1.01],
                eta=[cop],
                regions=[region])
                # transformer heating rod
            transformer.Simple(
                uid=('transformer', region.name, 'hp', 'air', 'ww', 'rod'),
                inputs=[elec_bus],
                outputs=[bus_hp],
                out_max=[max(demand.val) * share_heating_rod],
                eta=[eta_th['heat_rod_dec']],
                regions=[region])
                # heat storage
            transformer.Storage(
                uid=('storage', region.name, 'hp', 'air', 'ww'),
                inputs=[bus_hp],
                outputs=[bus_hp],
                cap_max=max(demand.val) * 2 * share_heat_storage,
                out_max=[max(demand.val) * share_heat_storage],
                in_max=[max(demand.val) * share_heat_storage],
                eta_in=eta_in['heat_storage_dec'],
                eta_out=eta_out['heat_storage_dec'],
                cap_loss=cap_loss['heat_storage_dec'],
                opex_fix=opex_fix['heat_storage_dec'])

    if share_air_hp != 1:
        # brine heat pump heating
            # bus
        bus_hp = Bus(uid=('bus', region.name, 'residential', 'heat_pump_dec',
                          'brine', 'heating'),
            type='heat_pump_dec', price=0, regions=[region], excess=False)
            # sink
        demand = sink.Simple(
            uid=('demand', region.name, 'residential', 'heat_pump_dec', 'brine',
                 'heating'),
            inputs=[bus_hp],
            region=region)
        demand.val = ((profile_sfh_heating + profile_mfh_heating) *
            (1 - share_air_hp))
            # transformer
        cop = cop_heating(temp, 'brine')
        transformer.Simple(
            uid=('transformer', region.name, 'hp', 'brine', 'heating'),
            inputs=[elec_bus],
            outputs=[bus_hp],
            out_max=[max(demand.val)],
            eta=[cop],
            regions=[region])
            # heat storage
        transformer.Storage(
            uid=('storage', region.name, 'hp', 'brine', 'heating'),
            inputs=[bus_hp],
            outputs=[bus_hp],
            cap_max=max(demand.val) * 2 * share_heat_storage,
            out_max=[max(demand.val) * share_heat_storage],
            in_max=[max(demand.val) * share_heat_storage],
            eta_in=eta_in['heat_storage_dec'],
            eta_out=eta_out['heat_storage_dec'],
            cap_loss=cap_loss['heat_storage_dec'],
            opex_fix=opex_fix['heat_storage_dec'])

        # brine heat pump warm water
            # bus
        if share_ww != 0:
            bus_hp = Bus(uid=('bus', region.name, 'residential', 'heat_pump_dec',
                              'brine', 'ww'),
                type='heat_pump_dec', price=0, regions=[region], excess=False)
                # demand
            demand = sink.Simple(
                uid=('demand', region.name, 'residential', 'heat_pump_dec', 'brine',
                     'ww'),
                inputs=[bus_hp],
                region=region)
            demand.val = (profile_sfh_ww + profile_mfh_ww) * (1 - share_air_hp)
                # transformer
            cop = cop_ww(temp, profile_sfh_ww, profile_mfh_ww, type_hp='brine')
            transformer.Simple(
                uid=('transformer', region.name, 'hp', 'brine', 'ww'),
                inputs=[elec_bus],
                outputs=[bus_hp],
                out_max=[max(demand.val)],
                eta=[cop],
                regions=[region])
                # heat storage
            transformer.Storage(
                uid=('storage', region.name, 'hp', 'brine', 'ww'),
                inputs=[bus_hp],
                outputs=[bus_hp],
                cap_max=max(demand.val) * 2 * share_heat_storage,
                out_max=[max(demand.val) * share_heat_storage],
                in_max=[max(demand.val) * share_heat_storage],
                eta_in=eta_in['heat_storage_dec'],
                eta_out=eta_out['heat_storage_dec'],
                cap_loss=cap_loss['heat_storage_dec'],
                opex_fix=opex_fix['heat_storage_dec'])
    return