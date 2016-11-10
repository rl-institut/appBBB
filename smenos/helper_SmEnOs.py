# -*- coding: utf-8 -*-
"""
@author: Elisa
"""
import logging
import os
import pandas as pd
import numpy as np
from datetime import time as settime

from oemof.core.network.entities.components import transformers as transformer
from oemof.core.network.entities.components import sources as source
import demandlib.bdew as bdew
import demandlib.particular_profiles as profiles
from oemof.tools import helpers
from oemof import db


def get_parameters():
    'returns emission and cost parameters'
    filename = os.path.abspath(os.path.join(os.path.dirname(__file__),
                         'simulation_parameter_App_SmEnOs_csv.csv'))
    read_parameter = pd.read_csv(filename, delimiter=';', index_col=0)
    parameters = {}

    for col in read_parameter.columns:
        parameters[col] = {}
        for row in read_parameter.index:
            try:
                parameters[col][row] = float(read_parameter.loc[row][col])
            except:
                parameters[col][row] = read_parameter.loc[row][col]

    # emission factors [t/MWh]
    co2_emissions = parameters['co2_var']
    co2_fix = parameters['co2_fix']

    # efficiencies [-]
    eta_elec = parameters['eta_elec']
    eta_th = parameters['eta_th']
    eta_el_chp = parameters['eta_el_chp']
    eta_th_chp = parameters['eta_th_chp']
    eta_chp_flex_el = parameters['eta_chp_flex_el']
    sigma_chp = parameters['sigma_chp']
    beta_chp = parameters['beta_chp']

    opex_var = parameters['opex_var']
    opex_fix = parameters['opex_fix']
    capex = parameters['capex']
    lifetime = parameters['lifetime']
    wacc = parameters['wacc']

    c_rate_in = parameters['c_rate_in']
    c_rate_out = parameters['c_rate_out']
    eta_in = parameters['eta_in']
    eta_out = parameters['eta_out']
    cap_loss = parameters['cap_loss']

    return(co2_emissions, co2_fix, eta_elec, eta_th, eta_th_chp, eta_el_chp,
           eta_chp_flex_el, sigma_chp, beta_chp, opex_var, opex_fix, capex,
           c_rate_in, c_rate_out, eta_in, eta_out, cap_loss, lifetime, wacc)


def get_data_from_csv(filename):
    '''returns "data". data is a dict which includes the columns of the csv
    as dicts: columns of csv: data[column]; cells of csv: data[column][row]
    '''
    read_data = pd.read_csv(filename, delimiter=';', index_col=0)
    data = {}

    for col in read_data.columns:
        data[col] = {}
        for row in read_data.index:
            try:
                data[col][row] = float(read_data.loc[row][col])
            except:
                data[col][row] = read_data.loc[row][col]
    return data


def get_res_parameters():
    site = {'module_name': 'Yingli_YL210__2008__E__',
            'azimuth': 0,
            'tilt': 0,
            'albedo': 0.2,
            'hoy': 8760,
            'h_hub': 135,
            'd_rotor': 127,
            'wka_model': 'ENERCON E 126 7500',
            'h_hub_dc': {
                1: 135,
                2: 135,
                3: 98,
                4: 78,
                0: 135},
            'd_rotor_dc': {
                1: 127,
                2: 82,
                3: 82,
                4: 82,
                0: 82},
            'wka_model_dc': {
                1: 'ENERCON E 82 2300',
                2: 'ENERCON E 82 2300',
                3: 'ENERCON E 82 3000',
                4: 'ENERCON E 82 3000',
                0: 'ENERCON E 82 2300'},
            }
    return site


def get_offshore_parameters(conn):
    site = {'h_hub': 80,
            'd_rotor': 120,
            'wka_model': 'SIEMENS SWT 3.6 120',
            }

    return site


def get_bdew_heatprofile_parameters():
    #TODO Werte recherchieren
    bdew_heatprofile_parameters = pd.DataFrame(
        [{'share_EFH': 0.5, 'wind_class': 1, 'building_class': 1},
            {'share_EFH': 0.5, 'wind_class': 1, 'building_class': 1},
            {'share_EFH': 0.5, 'wind_class': 1, 'building_class': 1},
            {'share_EFH': 0.5, 'wind_class': 1, 'building_class': 1},
            {'share_EFH': 0.5, 'wind_class': 1, 'building_class': 1},
            {'share_EFH': 0.5, 'wind_class': 1, 'building_class': 1}],
        index=['BE', 'BB', 'MV', 'SN', 'ST', 'TH'])
    return bdew_heatprofile_parameters


def get_hp_parameters():
    # share_hp_new_building und share_fbh_old_building können evtl auch
    # angepasst werden für ausbauszenarien (JAZ im Auge behalten)

    # share of single family houses of all residential buildings that have a
    # heat pump (share_mfh_hp = 1 - share_sfh_hp)
    share_sfh_hp = 0.5
    share_ww = 0  # share of warm water of total heating demand
    # share of air hp of all heat pumps (share_brine_hp = 1 - share_air_hp)
    share_air_hp = 0.6  # Anm.: Sole-WP hauptsächlich in Neubauten, sodass
                        # Anteil von Luft-WP bei Sanierungsszenarien steigt
    share_heating_rod = 0.42  # share of heating rod in monoenergetic hp
                              # system
    share_heat_storage = 0.85  # share of hp systems with heat storage
    return (share_sfh_hp, share_ww, share_air_hp, share_heating_rod,
        share_heat_storage)


def get_demand(conn, regions):
    sql = """
        SELECT nuts.nuts_id,
            ds.energy_sector AS energy,
            ds.consumption_sector AS sector,
            d.demand
        FROM oemof.demand AS d
        JOIN oemof.demand_sector AS ds ON d.sector=ds.id
        JOIN oemof.geo_nuts_rg_2013 AS nuts ON nuts.gid=d.region
        WHERE nuts.nuts_id IN
        """ + str(regions)
    df = pd.DataFrame(
        conn.execute(sql).fetchall(),
        columns=['nuts_id', 'energy', 'sector', 'demand'])
    return df


def get_biomass_between_5and10MW(conn, geometry):
    sql = """
    SELECT sum(p_nenn_kwp)
    FROM oemof_test.energy_map as ee
    WHERE anlagentyp='Biomasse'
    AND p_nenn_kwp >= 5000
    AND p_nenn_kwp < 10000
    AND st_contains(ST_GeomFromText('{wkt}',4326), ee.geom)
        """.format(wkt=geometry.wkt)
    cap = pd.DataFrame(conn.execute(sql).fetchall(), columns=['capacity'])
    return cap


def get_biomass_under_5MW(conn, geometry):
    sql = """
    SELECT sum(p_nenn_kwp)
    FROM oemof_test.energy_map as ee
    WHERE anlagentyp='Biomasse'
    AND p_nenn_kwp < 5000
    AND st_contains(ST_GeomFromText('{wkt}',4326), ee.geom)
        """.format(wkt=geometry.wkt)
    cap = pd.DataFrame(conn.execute(sql).fetchall(), columns=['capacity'])
    return cap


def get_opsd_pps(conn, geometry):
    de_en = {
        'Braunkohle': 'lignite',
        'lignite': 'lignite',
        'Steinkohle': 'hard_coal',
        'coal': 'hard_coal',
        'Erdgas': 'natural_gas',
        'Öl': 'oil',
        'oil': 'oil',
        'Solarstrom': 'solar_power',
        'Windkraft': 'wind_power',
        'Biomasse': 'biomass',
        'biomass': 'biomass',
        'Wasserkraft': 'hydro_power',
        'run_of_river': 'hydro_power',
        'Gas': 'methan',
        'Mineralölprodukte': 'mineral_oil',
        'Abfall': 'waste',
        'waste': 'waste',
        'Sonstige Energieträger\n(nicht erneuerbar) ': 'waste',
        'other_non_renewable': 'waste',
        'Pumpspeicher': 'pumped_storage',
        'pumped_storage': 'pumped_storage',
        'Erdwärme': 'geo_heat',
        'gas': 'natural_gas'}
    translator = lambda x: de_en[x]

    sql = """
        SELECT fuel, technology, status, chp, capacity, capacity_uba,
        chp_capacity_uba, efficiency_estimate
        FROM oemof_test.kraftwerke_de_opsd_2 as pp
        WHERE st_contains(
        ST_GeomFromText('{wkt}',4326), ST_Transform(pp.geom, 4326))
        """.format(wkt=geometry.wkt)
    df = pd.DataFrame(
        conn.execute(sql).fetchall(), columns=['type', 'technology', 'status',
                                                'chp', 'cap_el', 'cap_el_uba',
                                                'cap_th_uba', 'efficiency'])
    df['type'] = df['type'].apply(translator)

    #  rename cc-gasturbines in df
    for row in df.index:
        if df.loc[row]['technology'] == 'CC' and df.loc[
                                        row]['type'] == 'natural_gas':
                df.loc[row, 'type'] = 'natural_gas_cc'
    return df


def get_kwdb_pps(conn, geometry):
    de_en = {
        'lignite': 'lignite',
        'lignite': 'lignite',
        'hard_coal': 'hard_coal',
        'coal': 'hard_coal',
        'Erdgas': 'natural_gas',
        'Öl': 'oil',
        'mineral_oil': 'oil',
        'Solarstrom': 'solar_power',
        'solar':'solar_power',
        'uran':'uran',
        'Windkraft': 'wind_power',
        'refuse': 'waste',
        'Biomasse': 'biomass',
        'biomass': 'biomass',
        'Wasserkraft': 'hydro_power',
        'other_fuels': 'waste',
        'run_of_river': 'hydro_power',
        'gas': 'natural_gas',
        'mine_gas':'natural_gas',
        'landfill_gas':'natural_gas',
        'Mineralölprodukte': 'mineral_oil',
        'multiple_fuels': 'waste',
        'waste': 'waste',
        'Sonstige Energieträger\n(nicht erneuerbar) ': 'waste',
        'other_non_renewable': 'waste',
        'Pumpspeicher': 'pumped_storage',
        'pumped_storage': 'pumped_storage',
        'Erdwärme': 'geo_heat',
        'gas': 'natural_gas',
        'onshore': 'onshore'}
    translator = lambda x: de_en[x]

    sql = """
        SELECT fuel, kwk, status, rated_power
        FROM oemof.register_conventional_power_plants as pp
        WHERE st_contains(
        ST_GeomFromText('{wkt}',4326), ST_Transform(pp.geom, 4326))
        """.format(wkt=geometry.wkt)
    df = pd.DataFrame(
        conn.execute(sql).fetchall(), columns=['type', 'chp', 'status',
                                               'cap_el'])
    df['type'] = df['type'].apply(translator)

    return df


def get_hydro_energy(conn, regions):
    sql = """
        SELECT state_short, capacity_2013, energy_2013
        FROM oemof_test.hydro_energy_aee as ror
        WHERE state_short IN
        """ + str(regions)
    df = pd.DataFrame(
        conn.execute(sql).fetchall(), columns=[
            'state_short', 'capacity_mw', 'energy_mwh'])
    return df


def get_pumped_storage_pps(conn, regions):
    sql = """
        SELECT state_short, power_mw, capacity_mwh
        FROM oemof_test.pumped_storage_germany as pspp
        WHERE state_short IN
        """ + str(regions)
    df = pd.DataFrame(
        conn.execute(sql).fetchall(),
        columns=['state_short', 'power_mw', 'capacity_mwh'])
    return df


def get_offshore_pps(conn, schema, table, start_year):
    sql = """
        SELECT farm_capacity, st_astext(geom)
        FROM {schema}.{table}
        WHERE start_year <= {start_year}
        """.format(**{'schema': schema, 'table': table,
                      'start_year': start_year})
    pps = pd.DataFrame(conn.execute(sql).fetchall())
    return pps


def entity_exists(esystem, uid):
    return len([obj for obj in esystem.entities if obj.uid == uid]) > 0


def create_opsd_summed_objects(esystem, region, pp, **kwargs):

    'Creates entities for each type of generation'

    typeofgen = kwargs.get('typeofgen')
    ror_cap = kwargs.get('ror_cap')
    pumped_storage = kwargs.get('pumped_storage')
    chp_faktor_flex = kwargs.get('chp_faktor_flex', 0.84)
    cap_initial = kwargs.get('cap_initial', 0)

    (co2_emissions, co2_fix, eta_elec, eta_th, eta_th_chp, eta_el_chp,
     eta_chp_flex_el, sigma_chp, beta_chp, opex_var, opex_fix, capex,
     c_rate_in, c_rate_out, eta_in, eta_out,
     cap_loss, lifetime, wacc) = get_parameters()

    # replace NaN with 0
    mask = pd.isnull(pp)
    pp = pp.where(~mask, other=0)

    capacity = {}
    capacity_chp_el = {}

#_______________________________________get biomass under 10 MW from energymap
    conn = db.connection()
    p_bio = get_biomass_between_5and10MW(conn, region.geom)
    p_bio_bhkw = get_biomass_under_5MW(conn, region.geom)
    p_bio5to10 = float(p_bio['capacity']) / 1000  # from kW to MW
    p_el_bhkw = float(p_bio_bhkw['capacity']) / 1000  # from kW to MW
    p_th_bhkw = p_el_bhkw * eta_th_chp['bhkw'] / eta_el_chp['bhkw']

#___________________________________________________________

#    efficiency = {}

    for typ in typeofgen:
        # capacity for simple transformer (chp = no)
        capacity[typ] = sum(pp[pp['type'].isin([typ])][pp['status'].isin([
            'operating'])][pp['chp'].isin(['no'])]['cap_el'])
        # el capacity for chp transformer (cap_el_uba +
           # cap_el(where cap_th_uba=0 and chp=yes))
        capacity_chp_el[typ] = sum(pp[pp['type'].isin([typ])][pp[
            'status'].isin(['operating'])][pp['chp'].isin(['yes'])][
            'cap_el_uba']) + sum(pp[pp['type'].isin([typ])][pp['status'].isin([
            'operating'])][pp['chp'].isin(['yes'])][pp['cap_th_uba'].isin(
            [0])]['cap_el'])

        # efficiency only for simple transformer from OPSD-List

#        efficiency[typ] = np.mean(pp[pp['type'].isin([typ])][pp[
#        'status'].isin(['operating'])][pp['chp'].isin(['no'])]['efficiency'])

        # choose the right bus for type of generation (biomass: regional bus)
        # and add biomass capacity from energymap between 5 and 10 MW
        if typ == 'biomass':
            resourcebus = [obj for obj in esystem.entities if obj.uid ==
                "('bus', '"+region.name+"', '"+typ+"')"]

            capacity_chp_el[typ] = capacity_chp_el[typ] + p_bio5to10

            #create biomass bhkw transformer (under 5 MW)
            transformer.CHP(
                uid=('transformer', region.name, typ, 'bhkw'),
                # takes from resource bus
                inputs=resourcebus,
                # puts in electricity and heat bus
                outputs=[[obj for obj in region.entities if obj.uid ==
                        "('bus', '"+region.name+"', 'elec')"][0],
                        [obj for obj in region.entities if obj.uid ==
                        "('bus', '"+region.name+"', 'dh')"][0]],
                #TODO: Wärme auf den richtigen bus legen!!!!
                in_max=[None],
                out_max=[p_el_bhkw, p_th_bhkw],
                eta=[eta_el_chp['bhkw'], eta_th_chp['bhkw']],
                opex_var=opex_var[typ],
                regions=[region])
        else:  # not biomass
            resourcebus = [obj for obj in esystem.entities if obj.uid ==
                "('bus', 'global', '" +typ+ "')"]

        if capacity_chp_el[typ] > 0:

            transformer.CHP(
                uid=('transformer', region.name, typ, 'chp'),
                inputs=resourcebus,
                outputs=[[obj for obj in region.entities if obj.uid ==
                        "('bus', '"+region.name+"', 'elec')"][0],
                        [obj for obj in region.entities if obj.uid ==
                        "('bus', '"+region.name+"', 'dh')"][0]],
                in_max=[None],
                out_max=get_out_max_chp(capacity_chp_el[typ], chp_faktor_flex,
                                        eta_th_chp[typ], eta_el_chp[typ]),
                eta=[eta_el_chp[typ], eta_th_chp[typ]],
                opex_var=opex_var[typ],
                regions=[region])

            transformer.SimpleExtractionCHP(
                uid=('transformer', region.name, typ, 'SEchp'),
                inputs=resourcebus,
                outputs=[[obj for obj in region.entities if obj.uid ==
                        "('bus', '"+region.name+"', 'elec')"][0],
                        [obj for obj in region.entities if obj.uid ==
                        "('bus', '"+region.name+"', 'dh')"][0]],
                in_max=[None],
                out_max=get_out_max_chp_flex(capacity_chp_el[typ],
                                             chp_faktor_flex, sigma_chp[typ]),
                out_min=[0.0, 0.0],
                eta_el_cond=eta_chp_flex_el[typ],
                sigma=sigma_chp[typ],	 # power to heat ratio in backpr. mode
                beta=beta_chp[typ],		# power loss index
                opex_var=opex_var[typ],
                regions=[region])

        if capacity[typ] > 0:
            transformer.Simple(
                uid=('transformer', region.name, typ),
                inputs=resourcebus,
                outputs=[[obj for obj in region.entities if obj.uid ==
                        "('bus', '"+region.name+"', 'elec')"][0]],
                in_max=[None],
                out_max=[float(capacity[typ])],
                eta=[eta_elec[typ]],
                opex_var=opex_var[typ],
                regions=[region])

    # pumped storage
    typ = 'pumped_storage'
    power = sum(pumped_storage[pumped_storage[
        'state_short'].isin([region.name])]['power_mw'])
    if power > 0:
        transformer.Storage(
            uid=('Storage', region.name, typ),
            # nimmt von strombus
            inputs=[[obj for obj in esystem.entities if obj.uid ==
                    "('bus', '"+region.name+"', 'elec')"][0]],
            # speist auf strombus ein
            outputs=[[obj for obj in region.entities if obj.uid ==
                    "('bus', '"+region.name+"', 'elec')"][0]],
            cap_max=float(sum(pumped_storage[pumped_storage[
                'state_short'].isin([region.name])]['capacity_mwh'])),
            cap_min=0,
            in_max=[power],
            out_max=[power],
            eta_in=eta_in[typ],
            eta_out=eta_out[typ],
            c_rate_in=c_rate_in[typ],
            c_rate_out=c_rate_out[typ],
            #  opex_fix=opex_fix[typ],
            capex=capex[typ],
            cap_initial=cap_initial,
            regions=[region])

    # run of river
    typ = 'run_of_river'
    energy = sum(ror_cap[ror_cap['state_short'].isin(
        [region.name])]['energy_mwh'])

    capacity[typ] = sum(ror_cap[ror_cap[
        'state_short'].isin([region.name])]['capacity_mw'])
    if energy > 0:
        source.FixedSource(
            uid=('FixedSrc', region.name, 'hydro'),
            # speist auf strombus ein
            outputs=[obj for obj in region.entities if obj.uid ==
                    "('bus', '"+region.name+"', 'elec')"],
            val=scale_profile_to_sum_of_energy(
                filename=kwargs.get('filename_hydro'),
                energy=energy,
                capacity=capacity[typ]),
            out_max=[1],  # inst. Leistung!
            regions=[region])


def get_out_max_chp(capacity_chp_el, chp_faktor_flex,
                    eta_th_chp, eta_el_chp):
    out_max_el = float(capacity_chp_el) * (1-chp_faktor_flex)
    out_max_th = out_max_el * eta_th_chp / eta_el_chp
    out = [out_max_el, out_max_th]
    return out


def get_out_max_chp_flex(capacity_chp_el, chp_faktor_flex, sigma_chp):
    out_max_el = float(capacity_chp_el) * (chp_faktor_flex)
    out_max_th = out_max_el / sigma_chp
    out = [out_max_el, out_max_th]
    return out


def scale_profile_to_capacity(filename, capacity):
    profile = pd.read_csv(filename, sep=",")
    generation_profile = (profile.values / np.amax(profile.values) *
                          float(capacity))
    return generation_profile


def scale_profile_to_sum_of_energy(filename, energy, capacity):
    profile = pd.read_csv(filename)
    generation_profile = profile['profile'] * float(energy) / (
        sum(profile['profile']))
    return generation_profile


def el_load_profiles(demand, ann_el_demand_per_sector, year, **kwargs):

    # read standard load profiles
    e_slp = bdew.ElecSlp(year, holidays=kwargs.get('holidays', None))
    
    # multiply given annual demand with timeseries
    elec_demand = e_slp.get_profile(ann_el_demand_per_sector)
    
    # Add the slp for the industrial group
    ilp = profiles.IndustrialLoadProfile(
        e_slp.date_time_index, holidays=kwargs.get('holidays', None))
    
    # Beginning and end of workday, weekdays and weekend days, and scaling 
    # factors by default
    elec_demand['i0'] = ilp.simple_profile(ann_el_demand_per_sector['i0'],
            am=kwargs.get('am'), pm=kwargs.get('pm'),
            profile_factors=kwargs.get('profile_factors'))
    
    # Resample 15-minute values to hourly values.
    elec_demand = elec_demand.resample('H').mean()

    demand.val = elec_demand.sum(axis=1)
    return demand
    

def scale_profile(demand, year, filename, annual_demand):
    '''
    scale a given profile to a given annual demand, which is the sum
    of the single profile values
    '''
    profile = pd.read_csv(filename, sep=",")

    demand.val = (profile / profile.sum() * annual_demand)
    return demand


def call_heat_demandlib(region, time_index, **kwargs):
    '''
    Calls the demandlib and creates an object which includes the demand
    timeseries.

    Required Parameters
    -------------------
    demand : Sink object
    region : Region object
    '''
    load_profile = bdew.HeatBuilding(
        time_index,
        holidays=kwargs.get('holidays', None),
        temperature=region.temp,
        shlp_type=kwargs.get('shlp_type', None),
        building_class=(region.building_class
                        if region.building_class is not None else 0),
        wind_class=region.wind_class,
        annual_heat_demand=kwargs.get('annual_heat_demand', None),
        name=kwargs.get('shlp_type', None),
        ww_incl=kwargs.get('ww_incl', True)).get_bdew_profile()
    return load_profile


def ind_profile_parameters():
    am = settime(7, 0, 0)
    pm = settime(20, 00, 0)
    profile_factors = {'week': {'day': 0.8, 'night': 0.6},
                       'weekend': {'day': 0.9, 'night': 0.7}}
    return am, pm, profile_factors


def call_ind_profile(time_index, annual_demand, **kwargs):
    '''
    Creates an industrial load profile as step load profile.
    '''
    ilp = profiles.IndustrialLoadProfile(time_index,
                                         holidays=kwargs.get('holidays', None))
    load_profile = ilp.simple_profile(annual_demand,
        am=kwargs.get('am', None), pm=kwargs.get('pm', None),
        profile_factors=kwargs.get('profile_factors', None))
    load_profile = load_profile.resample('H').mean()
    return load_profile
