import os
import pandas as pd
import numpy as np
import pickle

from oemof.db import coastdat
from oemof.core.network.entities import Bus
from oemof.core.network.entities.components import sinks as sink
from oemof.core.network.entities.components import transformers as transformer

import helper_heat_pump as hhp
import helper_BBB as hlsb


def create_decentral_entities(Regions, regionsBBB, demands_df, year,
                              time_index, eta_th, eta_in, eta_out, cap_loss,
                              opex_fix, opex_var, eta_th_chp, eta_el_chp,
                              holidays):
                                  
    # dictionary to assign commodity to heating system
    heating_system_commodity = {
        'hard_coal_dec': 'hard_coal',
        'hard_coal_dec_ind': 'hard_coal',
        'lignite_dec': 'lignite',
        'lignite_dec_ind': 'lignite',
        'oil_dec': 'oil',
        'natural_gas_dec': 'natural_gas',
        'other_dec': 'waste',
        'biomass_dec': 'biomass',
        'biomass_dec_ind': 'biomass'}
    # get parameters for generation of heat load profiles 
    heat_params = hlsb.get_bdew_heatprofile_parameters()
    # get heat pump parameters
    (share_sfh_hp, share_ww, share_air_hp,
        share_heating_rod, share_heat_storage) = hlsb.get_hp_parameters()
    # get parameters for generation of industrial heat load profiles 
    am, pm, profile_factors = hlsb.ind_profile_parameters()
    
    for region in Regions.regions:
        if region.name == 'BE':
            global_bus = 'BE'
        else:
            global_bus = 'BB'
            
        # get regional heat demand for different ressources
        regID = region.name
        demand_sectors = demands_df.query('region==@regID')
        
        # get temperature of region as np array [°C]
        # used coastdat data cannot be provided, therefore only the used
        # temperature timeseries is provided
#        multiWeather = coastdat.get_weather(conn, region.geom, year)
#        temp = np.zeros([len(multiWeather[0].data.index), ])
#        for weather in multiWeather:
#            temp += weather.data['temp_air'].as_matrix()
#        temp = pd.Series(temp / len(multiWeather) - 273.15)
#        filename = os.path.abspath(
#            os.path.join(os.path.dirname(__file__), 'data',
#                        'temp_' + regID + '.pickle'))
#        pickle.dump(temp, open(filename, "wb" ))        
        filename = os.path.abspath(
            os.path.join(os.path.dirname(__file__), 'data',
                         'temp_' + regID + '.pickle'))
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
                    Bus(uid=('bus', region.name, sec, ressource),
                        type=ressource,
                        price=0, 
                        regions=[region], 
                        excess=False)
                    # create sink
                    demand = sink.Simple(
                        uid=('demand', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities if
                            obj.uid == ('bus', region.name, sec, ressource)],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    profile_efh = hlsb.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem * region.share_efh),
                        shlp_type='EFH', ww_incl=True)
                    profile_mfh = hlsb.call_heat_demandlib(region, time_index,
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
                            if obj.uid == (
                                'bus', region.name, sec, ressource)],
                        out_max=[max(demand.val)],
                        eta=[eta_th[ressource]],
                        regions=[region])
                elif ressource == 'dh':
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    profile_efh = hlsb.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(
                            data_dem * region.share_efh),
                        shlp_type='EFH', ww_incl=True)
                    profile_mfh = hlsb.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem *
                            (1 - region.share_efh)),
                        shlp_type='MFH', ww_incl=True)
                    dh_demand += profile_efh + profile_mfh
                elif ressource == 'bhkw_bio':
                    # create bus(bedarfsbus)
                    Bus(uid=("('bus', '" + region.name + "', '" + sec +
                            "', '" + ressource+"')"),
                        type=ressource,
                        price=0, 
                        regions=[region], 
                        excess=False)
                    # create sink
                    demand = sink.Simple(
                        uid=('demand', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == ("('bus', '" + region.name + "', '" +
                            sec + "', '" + ressource + "')")],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    profile_efh = hlsb.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem * region.share_efh),
                        shlp_type='EFH', ww_incl=True)
                    profile_mfh = hlsb.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem *
                            (1 - region.share_efh)),
                        shlp_type='MFH', ww_incl=True)
                    demand.val = (profile_efh + profile_mfh)
                    # create transformer#
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities if obj.uid == 
                                "('bus', '" + global_bus + "', 'biomass')"],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '" + region.name + "', 'elec')"][0], 
                                [obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + region.name + 
                                "', '" + sec + "', '" + ressource + "')")][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource], 
                                eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource], eta_th_chp[ressource]],
                        regions=[region])
                elif ressource == 'bhkw_gas':
                    # create bus(bedarfsbus)
                    Bus(uid=("('bus', '" + region.name + "', '" + sec + 
                             "', '" + ressource+"')"),
                        type=ressource,
                        price=0, 
                        regions=[region], 
                        excess=False)
                    # create sink
                    demand = sink.Simple(
                        uid=('demand', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == ("('bus', '" + region.name + "', '" + 
                            sec + "', '" + ressource + "')")],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    profile_efh = hlsb.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem * region.share_efh),
                        shlp_type='EFH', ww_incl=True)
                    profile_mfh = hlsb.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=(data_dem *
                            (1 - region.share_efh)),
                        shlp_type='MFH', ww_incl=True)
                    demand.val = (profile_efh + profile_mfh)
                    # create transformer
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + global_bus + 
                                "', 'natural_gas')")],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '" + region.name + "', 'elec')"][0],
                                [obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + region.name + 
                                "', '" + sec + "', '" + ressource + "')")][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource], 
                                eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource], eta_th_chp[ressource]],
                        regions=[region])
                        
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
                    Bus(
                        uid=('bus', region.name, sec, ressource),
                        type=ressource,
                        price=0, 
                        regions=[region], 
                        excess=False)
                    # create sink
                    demand = sink.Simple(
                        uid=('demand', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == (
                                'bus', region.name, sec, ressource)],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = hlsb.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=data_dem,
                        shlp_type='GHD', ww_incl=True)
                    # create transformer
                    transformer.Simple(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + global_bus +
                                "', '" + heating_system_commodity[ressource] +
                                "')")],
                        outputs=[obj for obj in Regions.entities
                            if obj.uid == (
                                'bus', region.name, sec, ressource)],
                        out_max=[max(demand.val)],
                        eta=[eta_th[ressource]],
                        regions=[region])
                elif ressource == 'dh':
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    dh_demand += hlsb.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=data_dem,
                        shlp_type='GHD', ww_incl=True)
                elif ressource == 'bhkw_bio':
                    # create bus
                    Bus(
                        uid=("('bus', '" + region.name + "', '" + sec +
                        "', '" + ressource + "')"),
                        type=ressource,
                        price=0,
                        regions=[region],
                        excess=False)
                    # create sink
                    demand = sink.Simple(
                        uid=('demand', region.name, sec,
                        ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == ("('bus', '" + region.name + "', '" +
                            sec + "', '" + ressource + "')")],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = hlsb.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=data_dem,
                        shlp_type='GHD', ww_incl=True)
                    # create transformer
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + global_bus + 
                                    "', 'biomass')")],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '" + region.name + "', 'elec')"][0], 
                                [obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + region.name + 
                                "', '" + sec + "', '" + ressource + "')")][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource],
                                eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource],
                        eta_th_chp[ressource]],
                        regions=[region])
                elif ressource == 'bhkw_gas':
                    # create bus
                    Bus(
                        uid="('bus', '" + region.name + "', '" + sec + 
                            "', '" + ressource + "')",
                        type=ressource,
                        price=0, 
                        regions=[region], 
                        excess=False)
                    # create sink
                    demand = sink.Simple(
                        uid=('demand', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == ("('bus', '" + region.name + "', '" +
                                sec + "', '" + ressource + "')")],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = hlsb.call_heat_demandlib(region, time_index,
                        holidays=holidays,
                        annual_heat_demand=data_dem,
                        shlp_type='GHD', ww_incl=True)
                    # create transformer
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + global_bus +
                                    "', 'natural_gas')")],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '" + region.name + "', 'elec')"][0], 
                                [obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + region.name +
                                "', '" + sec + "', '" + ressource + "')")][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource], 
                                eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource], 
                        eta_th_chp[ressource]],
                        regions=[region])
                        
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
                    Bus(uid=('bus', region.name, sec, ressource), 
                        type=ressource,
                        price=0, 
                        regions=[region], 
                        excess=False)
                    # create sink
                    demand = sink.Simple(
                        uid=('demand', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities if
                            obj.uid == ('bus', region.name, sec, ressource)],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = (hlsb.call_ind_profile(
                        time_index, data_dem, holidays=holidays,
                        am=am, pm=pm, profile_factors=profile_factors))
                    # create transformer
                    transformer.Simple(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + global_bus +
                                "', '" + heating_system_commodity[ressource] +
                                "')")],
                        outputs=[obj for obj in Regions.entities if
                            obj.uid == ('bus', region.name, sec, ressource)],
                        out_max=[max(demand.val)],
                        eta=[eta_th[ressource]],
                        regions=[region])
                elif ressource == 'bhkw_bio':
                    # create bus
                    Bus(uid=("('bus', '" + region.name + "', '" + sec + 
                            "', '" + ressource + "')"),
                        type=ressource,
                        price=0, 
                        regions=[region], 
                        excess=False)
                    # create sink
                    demand = sink.Simple(
                        uid=('demand', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == ("('bus', '" + region.name + "', '" +
                                sec + "', '" + ressource + "')")],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = (hlsb.call_ind_profile(
                        time_index, data_dem, holidays=holidays,
                        am=am, pm=pm, profile_factors=profile_factors))
                    # create transformer
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + global_bus +
                                    "', 'biomass')")],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '" + region.name + "', 'elec')"][0], 
                                [obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + region.name +
                                    "', '" + sec + "', '" + ressource + 
                                    "')")][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource], 
                                eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource], 
                        eta_th_chp[ressource]],
                        regions=[region])
                elif ressource == 'bhkw_gas':
                    # create bus
                    Bus(uid=("('bus', '" + region.name + "', '" + sec + 
                             "', '" + ressource + "')"),
                        type=ressource,
                        price=0, 
                        regions=[region], 
                        excess=False)
                    # create sink
                    demand = sink.Simple(
                        uid=('demand', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                            if obj.uid == ("('bus', '" + region.name + "', '" +
                                sec + "', '" + ressource + "')")],
                        region=region)
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    demand.val = (hlsb.call_ind_profile(
                        time_index, data_dem, holidays=holidays,
                        am=am, pm=pm, profile_factors=profile_factors))
                    # create transformer
                    transformer.CHP(
                        uid=('transformer', region.name, sec, ressource),
                        inputs=[obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + global_bus +
                                    "', 'natural_gas')")],
                        outputs=[[obj for obj in Regions.entities if obj.uid ==
                                 "('bus', '"+region.name+"', 'elec')"][0], 
                                [obj for obj in Regions.entities
                                if obj.uid == ("('bus', '" + region.name +
                                "', '" + sec + "', '" + ressource + "')")][0]],
                        out_max=hlsb.get_out_max_chp(
                                max(demand.val), eta_th_chp[ressource], 
                                eta_el_chp[ressource]),
                        eta=[eta_el_chp[ressource], eta_th_chp[ressource]],
                        regions=[region])
                elif ressource == 'dh':
                    # create heat load profile and write to sink object
                    # heat load in [MWh/a]
                    dh_demand += hlsb.call_ind_profile(
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
            regions=[region])
