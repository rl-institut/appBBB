# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 13:46:12 2016

@author: elisa
"""


import pandas as pd
import numpy as np
import os

#from oemof.db import coastdat
from oemof.core.network.entities import Bus
from oemof.core.network.entities.components import sinks as sink
from oemof.core.network.entities.components import transformers as transformer

import helper_SmEnOs as hls
import helper_heat_pump as hhp
import helper_BBB as hlsb


def create_decentral_entities(Regions, regionsBBB, demands_df, conn, year,
                              time_index, eta_th, eta_in, eta_out, cap_loss,
                              opex_fix, opex_var, eta_th_chp, eta_el_chp,
                              holidays):

    heating_system_commodity = {
        'hard_coal_dec': 'hard_coal',
        'hard_coal_dec_ind': 'hard_coal', #  bedarf und Bus
        'lignite_dec': 'lignite',
        'lignite_dec_ind': 'lignite',
        'oil_dec': 'oil',
        'natural_gas_dec': 'natural_gas',
        'other_dec': 'waste',
        'biomass_dec': 'biomass',
        'biomass_dec_ind': 'biomass'}
    heat_params = hls.get_bdew_heatprofile_parameters()
    (share_sfh_hp, share_ww, share_air_hp,
        share_heating_rod, share_heat_storage) = hls.get_hp_parameters()
    am, pm, profile_factors = hls.ind_profile_parameters()
    for region in Regions.regions:
        if region.name == 'BE':
            global_bus = 'BE'
        else:
            global_bus = 'BB'
        # get regional heat demand for different ressources
        regID = region.name
        demand_sectors = demands_df.query('region==@regID')
        ## get temperature of region as np array [°C]
        #multiWeather = coastdat.get_weather(conn, region.geom, year)
        #temp = np.zeros([len(multiWeather[0].data.index), ])
        #for weather in multiWeather:
            #temp += weather.data['temp_air'].as_matrix()
        #temp = pd.Series(temp / len(multiWeather) - 273.15)
        #region.temp = temp
        filename = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                        'temp'))
        temp = pd.read_pickle(filename)
        region.temp = temp
        # create empty dataframe for district heating demand
        dh_demand = pd.Series(0, index=time_index)
        # residential sector
        sec = 'HH'
        region.share_efh = heat_params.ix[global_bus]['share_EFH']
        region.building_class = heat_params.ix[global_bus]['building_class']
        region.wind_class = heat_params.ix[global_bus]['wind_class']
        for ressource in list(demand_sectors.query("sector==@sec")['type']):
            data_dem = float(demand_sectors.query(
                        "sector==@sec and type==@ressource")['demand'])
            if data_dem != 0:
                if ressource == 'electricity':
                    pass
                if ressource == 'heat_pump_dec':
                    elec_bus = [obj for obj in Regions.entities
                        if obj.uid == "('bus', '%s', 'elec')" % (regID)][0]
                    hhp.create_hp_entities(region, year, data_dem,
                        elec_bus, temp, share_sfh_hp, share_ww,
                        share_air_hp, share_heating_rod, share_heat_storage,
                        eta_th, eta_in, eta_out, cap_loss, opex_fix,
                        holidays)
                elif ressource in list(heating_system_commodity.keys()):
                    # create bus(bedarfsbus)
                    Bus(uid=('bus', region.name, sec, ressource), type=ressource,
                        price=0, regions=[region], excess=False)
                    # create sink
                    demand = sink.Simple(uid=('demand', region.name, sec,
                        ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == ('bus', region.name, sec, ressource)],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    profile_efh = hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem * region.share_efh),
                        shlp_type='EFH', ww_incl=True)
                    profile_mfh = hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem *
                            (1 - region.share_efh)),
                        shlp_type='MFH', ww_incl=True)
                    demand.val = (profile_efh + profile_mfh)
                    # create transformer
                    transformer.Simple(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+global_bus+"', '"+\
                                heating_system_commodity[ressource]+"')"],
                        outputs=[obj for obj in Regions.entities
                            if obj.uid == ('bus', region.name, sec, ressource)],
                        out_max=[max(demand.val)],
                        eta=[eta_th[ressource]],
                        opex_var=opex_var[ressource],
                        regions=[region])
                elif ressource == 'dh':
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    profile_efh = hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(
                            data_dem * region.share_efh),
                        shlp_type='EFH', ww_incl=True)
                    profile_mfh = hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem *
                            (1 - region.share_efh)),
                        shlp_type='MFH', ww_incl=True)
                    dh_demand += profile_efh + profile_mfh
                elif ressource == 'bhkw_bio':
                    # create bus(bedarfsbus)
                    Bus(uid="('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')", type=ressource,
                        price=0, regions=[region], excess=False)
                    # create sink
                    demand = sink.Simple(uid=('demand', region.name, sec,
                        ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    profile_efh = hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem * region.share_efh),
                        shlp_type='EFH', ww_incl=True)
                    profile_mfh = hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem *
                            (1 - region.share_efh)),
                        shlp_type='MFH', ww_incl=True)
                    demand.val = (profile_efh + profile_mfh)
                    # create transformer#
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities if obj.uid == 
                                "('bus', '"+global_bus+"', 'biomass')"],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '"+region.name+"', 'elec')"][0], 
                                [obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource], 
                                eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource], eta_th_chp[ressource]],
                        opex_var=opex_var[ressource],
                        regions=[region])
                elif ressource == 'bhkw_gas':
                    # create bus(bedarfsbus)
                    Bus(uid="('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')", type=ressource,
                        price=0, regions=[region], excess=False)
                    # create sink
                    demand = sink.Simple(uid=('demand', region.name, sec,
                        ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    profile_efh = hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem * region.share_efh),
                        shlp_type='EFH', ww_incl=True)
                    profile_mfh = hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem *
                            (1 - region.share_efh)),
                        shlp_type='MFH', ww_incl=True)
                    demand.val = (profile_efh + profile_mfh)
                    # create transformer
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+global_bus+"', 'natural_gas')"],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '"+region.name+"', 'elec')"][0],
                                [obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource], eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource], eta_th_chp[ressource]],
                        opex_var=opex_var[ressource],
                        regions=[region])
                else:
                    print('Folgender Bedarf wird nicht berücksichtigt:')
                    print(region.name)
                    print(sec)
                    print(ressource)
        # commercial sector
        sec = 'GHD'
        region.wind_class = heat_params.ix[global_bus]['wind_class']
        region.building_class = None
        for ressource in list(demand_sectors.query("sector==@sec")['type']):
            data_dem = float(demand_sectors.query(
                        "sector==@sec and type==@ressource")['demand'])
            if data_dem != 0:
                if ressource == 'electricity':
                    pass                
                if ressource in list(heating_system_commodity.keys()):
                    # create bus
                    Bus(uid=('bus', region.name, sec, ressource), type=ressource,
                        price=0, regions=[region], excess=False)
                    # create sink
                    demand = sink.Simple(uid=('demand', region.name, sec,
                        ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == ('bus', region.name, sec, ressource)],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=data_dem,
                        shlp_type='GHD', ww_incl=True)
                    # create transformer
                    transformer.Simple(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+global_bus+"', '"+heating_system_commodity[ressource]+"')"],
                        outputs=[obj for obj in Regions.entities
                            if obj.uid == ('bus', region.name, sec, ressource)],
                        out_max=[max(demand.val)],
                        eta=[eta_th[ressource]],
                        opex_var=opex_var[ressource],
                        regions=[region])
                elif ressource == 'dh':
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    dh_demand += hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=data_dem,
                        shlp_type='GHD', ww_incl=True)
                elif ressource == 'bhkw_bio':
                    # create bus
                    Bus(uid="('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')", type=ressource,
                        price=0, regions=[region], excess=False)
                    # create sink
                    demand = sink.Simple(uid=('demand', region.name, sec,
                        ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=data_dem,
                        shlp_type='GHD', ww_incl=True)
                    # create transformer
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+global_bus+"', 'biomass')"],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '"+region.name+"', 'elec')"][0], 
                                [obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource], eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource], eta_th_chp[ressource]],
                        opex_var=opex_var[ressource],
                        regions=[region])
                elif ressource == 'bhkw_gas':
                    # create bus
                    Bus(uid="('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')", type=ressource,
                        price=0, regions=[region], excess=False)
                    # create sink
                    demand = sink.Simple(uid=('demand', region.name, sec,
                        ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = hls.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=data_dem,
                        shlp_type='GHD', ww_incl=True)
                    # create transformer
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+global_bus+"', 'natural_gas')"],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '"+region.name+"', 'elec')"][0], 
                                [obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource], eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource], eta_th_chp[ressource]],
                        opex_var=opex_var[ressource],
                        regions=[region])
                else:
                    print('Folgender Bedarf wird nicht berücksichtigt:')
                    print(region.name)
                    print(sec)
                    print(ressource)
        # industrial sector
        sec = 'IND'
        for ressource in list(demand_sectors.query("sector==@sec")['type']):
            data_dem = float(demand_sectors.query(
                        "sector==@sec and type==@ressource")['demand'])
            if data_dem != 0:
                if ressource == 'electricity':
                    pass
                if ressource in list(heating_system_commodity.keys()):
                    # create bus
                    Bus(uid=('bus', region.name, sec, ressource), type=ressource,
                        price=0, regions=[region], excess=False)
                    # create sink
                    demand = sink.Simple(uid=('demand', region.name, sec,
                        ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == ('bus', region.name, sec, ressource)],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = (hls.call_ind_profile(
                        time_index, data_dem, holidays=holidays,
                        am=am, pm=pm, profile_factors=profile_factors))
                    # create transformer
                    transformer.Simple(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+global_bus+"', '"+\
                                heating_system_commodity[ressource]+"')"],
                        outputs=[obj for obj in Regions.entities
                            if obj.uid == ('bus', region.name, sec, ressource)],
                        out_max=[max(demand.val)],
                        eta=[eta_th[ressource]],
                        opex_var=opex_var[ressource],
                        regions=[region])
                elif ressource == 'bhkw_bio':
                    # create bus
                    Bus(uid="('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')", type=ressource,
                        price=0, regions=[region], excess=False)
                    # create sink
                    demand = sink.Simple(uid=('demand', region.name, sec,
                        ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = (hls.call_ind_profile(
                        time_index, data_dem, holidays=holidays,
                        am=am, pm=pm, profile_factors=profile_factors))
                    # create transformer
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+global_bus+"', 'biomass')"],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '"+region.name+"', 'elec')"][0], 
                                [obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource], eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource], eta_th_chp[ressource]],
                        opex_var=opex_var[ressource],
                        regions=[region])
                elif ressource == 'bhkw_gas':
                    # create bus
                    Bus(uid="('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')", type=ressource,
                        price=0, regions=[region], excess=False)
                    # create sink
                    demand = sink.Simple(uid=('demand', region.name, sec,
                        ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = (hls.call_ind_profile(
                        time_index, data_dem, holidays=holidays,
                        am=am, pm=pm, profile_factors=profile_factors))
                    # create transformer
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+global_bus+"', 'natural_gas')"],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '"+region.name+"', 'elec')"][0], 
                                [obj for obj in Regions.entities
                                if obj.uid == "('bus', '"+region.name+"', '"+sec+"', '"+ressource+"')"][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource], eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource], eta_th_chp[ressource]],
                        opex_var=opex_var[ressource],
                        regions=[region])
                elif ressource == 'dh':
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    dh_demand += hls.call_ind_profile(
                        time_index, data_dem, holidays=holidays,
                        am=am, pm=pm, profile_factors=profile_factors)
    
        # create dh sink
            # create sink
        demand = sink.Simple(uid=('demand', region.name, 'dh'),
            inputs=[obj for obj in Regions.entities
                    if obj.uid == "('bus', '"+region.name+"', 'dh')"],
            region=region)

        demand.val = dh_demand
        transformer.Simple(
            uid=('transformer', region.name, 'dh_peak_heating'),
            inputs=[obj for obj in Regions.entities
                    if obj.uid == "('bus', '"+global_bus+"', 'natural_gas')"],
            outputs=[obj for obj in Regions.entities
                    if obj.uid == "('bus', '"+region.name+"', 'dh')"],
            out_max=[max(demand.val)],
            eta=[eta_th['natural_gas']],
            opex_var=opex_var['natural_gas'],
            regions=[region])
