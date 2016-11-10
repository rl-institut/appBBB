# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 08:19:23 2016

@author: Elisa.Gaudchau
"""

import logging
import pandas as pd
import numpy as np
from datetime import time as settime
import pyomo.environ as po

from oemof.core.network.entities.components import transformers as transformer
from oemof.core.network.entities.components import sources as source
from oemof.tools import helpers
from oemof import db
import helper_SmEnOs as hls

def get_parameters(conn_oedb):
    'returns emission and cost parameters'

    sql = """
        SELECT technology, co2_var, co2_fix, eta_elec,
             eta_th, eta_el_chp,
             eta_th_chp, eta_chp_flex_el, sigma_chp, beta_chp,
             opex_var, opex_fix, capex, c_rate_in, c_rate_out,
             eta_in, eta_out, cap_loss, lifetime, wacc
        FROM model_draft.abbb_simulation_parameter AS d
        """
    read_parameter = pd.DataFrame(
        conn_oedb.execute(sql).fetchall(),
        columns=['technology', 'co2_emissions', 'co2_fix', 'eta_elec',
                 'eta_th', 'eta_el_chp',
                 'eta_th_chp', 'eta_chp_flex_el', 'sigma_chp', 'beta_chp',
                 'opex_var', 'opex_fix', 'capex', 'c_rate_in', 'c_rate_out',
                 'eta_in', 'eta_out', 'cap_loss', 'lifetime', 'wacc'])

    parameters = {}
    for col in read_parameter.columns:
        parameters[col] = {}
        for row in read_parameter.index:
            key = read_parameter.technology[row]
            try:
                parameters[col][key] = float(read_parameter.loc[row][col])
            except:
                parameters[col][key] = read_parameter.loc[row][col]

    # emission factors [t/MWh]
    co2_emissions = parameters['co2_emissions']
    co2_fix = parameters['co2_fix']

    # efficiencies [-]
    eta_elec = parameters['eta_elec']
    eta_th = parameters['eta_th']
    eta_el_chp = parameters['eta_el_chp']
    eta_th_chp = parameters['eta_th_chp']

    eta_chp_flex_el = parameters['eta_chp_flex_el']
    eta_chp_flex_el['jaenschwalde'] = 0.35
    eta_chp_flex_el['schwarzepumpe'] = 0.39

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


def get_transmission(conn_oedb, scenario_name):
    'transmission capacities between BBB regions'

    sql = """
        SELECT from_region, to_region, capacity
        FROM model_draft.abbb_transmission_capacity AS d
        WHERE scenario = '""" + str(scenario_name) + """'"""
    read_parameter = pd.DataFrame(
        conn_oedb.execute(sql).fetchall(),
        columns=['from', 'to', 'cap'])

    transmission = get_dict_from_df(read_parameter)

    return(transmission)


def get_demand(conn_oedb, scenario_name):
    'heat and electrical demands in BBB regions'

    sql = """
        SELECT region, sector, type, demand
        FROM model_draft.abbb_demand AS d
        WHERE scenario = '""" + str(scenario_name) + """'"""
    read_parameter = pd.DataFrame(
        conn_oedb.execute(sql).fetchall(),
        columns=['region', 'sector', 'type', 'demand'])

    return(read_parameter)


def get_transformer(conn_oedb, scenario_name):
    'transformers in BBB regions'

    sql = """
        SELECT region, ressource, transformer, power
        FROM model_draft.abbb_transformer AS d
        WHERE scenario = '""" + str(scenario_name) + """'"""
    read_parameter = pd.DataFrame(
        conn_oedb.execute(sql).fetchall(),
        columns=['region', 'ressource', 'transformer', 'power'])

    return(read_parameter)


def get_dict_from_df(data_frame):
    '''returns "data". data is a dict which includes the columns of the table
    as dicts: columns of table: data[column]; cells of table: data[column][row]
    '''
    data = {}

    for col in data_frame.columns:
        data[col] = {}
        for row in data_frame.index:
            try:
                data[col][row] = float(data_frame.loc[row][col])
            except:
                data[col][row] = data_frame.loc[row][col]
    return data

def get_st_timeline(conn, year):
    sql = """
        SELECT "EFH_Altbau_san_HWW_S", "EFH_Altbau_san_HWW_SO",
        "EFH_Altbau_san_HWW_SW",
        "EFH_Altbau_san_WW_S", "EFH_Altbau_san_WW_SO",
        "EFH_Altbau_san_WW_SW",
        "MFH_Altbau_san_HWW_S", "MFH_Altbau_san_HWW_SO",
        "MFH_Altbau_san_HWW_SW",
        "MFH_Altbau_san_WW_S", "MFH_Altbau_san_WW_SO",
        "MFH_Altbau_san_WW_SW"
        FROM wittenberg.trnsys_st_zeitreihen AS d
        """
    read_parameter = pd.DataFrame(
            conn.execute(sql).fetchall(),
            columns=['EFH_Altbau_san_HWW_S', 'EFH_Altbau_san_HWW_SO',
                     'EFH_Altbau_san_HWW_SW',
                     'EFH_Altbau_san_WW_S', 'EFH_Altbau_san_WW_SO',
                     'EFH_Altbau_san_WW_SW',
                     'MFH_Altbau_san_HWW_S', 'MFH_Altbau_san_HWW_SO',
                     'MFH_Altbau_san_HWW_SW',
                     'MFH_Altbau_san_WW_S', 'MFH_Altbau_san_WW_SO',
                     'MFH_Altbau_san_WW_SW'])

    #timeline_st = pd.DataFrame(columns='st')
    timeline_st = pd.DataFrame(
        index=pd.date_range(pd.datetime(year, 1, 1, 0), periods=8760,
                            freq='H'), columns=['st'])

    for row in range(8760):
            timeline_st['st'][row] =\
                        0.15 * read_parameter['EFH_Altbau_san_HWW_S'][row] /\
                        read_parameter['EFH_Altbau_san_HWW_S'].sum() +\
                        0.05 * read_parameter['EFH_Altbau_san_HWW_SO'][row] /\
                        read_parameter['EFH_Altbau_san_HWW_SO'].sum() +\
                        0.05 * read_parameter['EFH_Altbau_san_HWW_SW'][row] /\
                        read_parameter['EFH_Altbau_san_HWW_SW'].sum() +\
                        0.15 * read_parameter['EFH_Altbau_san_WW_S'][row] /\
                        read_parameter['EFH_Altbau_san_WW_S'].sum() +\
                        0.05 * read_parameter['EFH_Altbau_san_WW_SO'][row] /\
                        read_parameter['EFH_Altbau_san_WW_SO'].sum() +\
                        0.05 * read_parameter['EFH_Altbau_san_WW_SW'][row] /\
                        read_parameter['EFH_Altbau_san_WW_SW'].sum() +\
                        0.15 * read_parameter['MFH_Altbau_san_HWW_S'][row] /\
                        read_parameter['MFH_Altbau_san_HWW_S'].sum() +\
                        0.05 * read_parameter['MFH_Altbau_san_HWW_SO'][row] /\
                        read_parameter['MFH_Altbau_san_HWW_SO'].sum() +\
                        0.05 * read_parameter['MFH_Altbau_san_HWW_SW'][row] /\
                        read_parameter['MFH_Altbau_san_HWW_SW'].sum() +\
                        0.15 * read_parameter['MFH_Altbau_san_WW_S'][row] /\
                        read_parameter['MFH_Altbau_san_WW_S'].sum() +\
                        0.05 * read_parameter['MFH_Altbau_san_WW_SO'][row] / \
                        read_parameter['MFH_Altbau_san_WW_SO'].sum() +\
                        0.05 * read_parameter['MFH_Altbau_san_WW_SW'][row] /\
                        read_parameter['MFH_Altbau_san_WW_SW'].sum()
    return(timeline_st)


def create_transformer(esystem, region, pp, conn, **kwargs):

    'Creates entities for each type of generation'

    typeofgen = kwargs.get('typeofgen')
    chp_faktor_flex = kwargs.get('chp_faktor_flex', 0.84)
    cap_initial = kwargs.get('cap_initial', 0)

    (co2_emissions, co2_fix, eta_elec, eta_th, eta_th_chp, eta_el_chp,
     eta_chp_flex_el, sigma_chp, beta_chp, opex_var, opex_fix, capex,
     c_rate_in, c_rate_out, eta_in, eta_out,
     cap_loss, lifetime, wacc) = get_parameters(conn)

    for typ in typeofgen:
        if region.name == 'BE':
            resourcebus = [obj for obj in esystem.entities if obj.uid ==
                           "('bus', 'BE', '"+typ+"')"]
        else:
            resourcebus = [obj for obj in esystem.entities if obj.uid ==
                           "('bus', 'BB', '"+typ+"')"]

        ########################## BHKW ,30% Fernwärme, 70 nur strom  ########
        if typ == 'bhkw_bio':
            try:
                capacity = float(pp.query(
                    'region==@region.name and ressource==@typ and transformer=="bhkw"')[
                    'power'])
            except:
                capacity = 0
            if capacity > 0:
                cap = capacity * 0.3
                transformer.CHP(
                    uid=('transformer', region.name, typ, 'dh'),
                    inputs=[obj for obj in esystem.entities if obj.uid ==
                           "('bus', 'BB', 'biomass')"],
                    outputs=[[obj for obj in region.entities if obj.uid ==
                             "('bus', '"+region.name+"', 'elec')"][0],
                            [obj for obj in region.entities if obj.uid ==
                             "('bus', '"+region.name+"', 'dh')"][0]],
                    in_max=[None],
                    out_max=get_out_max_chp(
                            cap, eta_th_chp[typ], eta_el_chp[typ]),
                    eta=[eta_el_chp[typ], eta_th_chp[typ]],
                    opex_var=opex_var[typ],
                    co2_var=co2_emissions[typ],
                    regions=[region])

                cap_2 = capacity * 0.7
                transformer.Simple(
                    uid=('transformer', region.name, typ),
                    inputs=[obj for obj in esystem.entities if obj.uid ==
                           "('bus', 'BB', 'biomass')"],
                    outputs=[[obj for obj in region.entities if obj.uid ==
                         "('bus', '"+region.name+"', 'elec')"][0]],
                    in_max=[None],
                    out_max=[cap_2],
                    eta=[eta_el_chp[typ]],
                    opex_var=opex_var[typ],
                    co2_var=co2_emissions[typ],
                    regions=[region])

        ########################## CHP #################################
        try:
            capacity = float(pp.query(
                'region==@region.name and ressource==@typ and transformer=="chp"')[
                'power'])
        except:
            capacity = 0
        if capacity > 0:
            transformer.CHP(
                uid=('transformer', region.name, typ, 'chp'),
                inputs=resourcebus,
                outputs=[[obj for obj in region.entities if obj.uid ==
                         "('bus', '"+region.name+"', 'elec')"][0],
                        [obj for obj in region.entities if obj.uid ==
                         "('bus', '"+region.name+"', 'dh')"][0]],
                in_max=[None],
                out_max=get_out_max_chp(
                        capacity, eta_th_chp[typ], eta_el_chp[typ]),
                eta=[eta_el_chp[typ], eta_th_chp[typ]],
                opex_var=opex_var[typ],
                co2_var=co2_emissions[typ],
                regions=[region])

        ########################## SE_chp #################################
        try:
            capacity = float(pp.query(
                'region==@region.name and ressource==@typ and transformer=="SE_chp"')[
                'power'])
        except:
            capacity = 0
        if capacity > 0:
            transformer.SimpleExtractionCHP(
                uid=('transformer', region.name, typ, 'SEchp'),
                inputs=resourcebus,
                outputs=[[obj for obj in region.entities if obj.uid ==
                         "('bus', '"+region.name+"', 'elec')"][0],
                        [obj for obj in region.entities if obj.uid ==
                         "('bus', '"+region.name+"', 'dh')"][0]],
                in_max=[None],
                out_max=get_out_max_chp_flex(capacity, sigma_chp[typ]),
                out_min=[0.0, 0.0],
                eta_el_cond=eta_chp_flex_el[typ],
                sigma=sigma_chp[typ],	 # power to heat ratio in backpr. mode
                beta=beta_chp[typ],		# power loss index
                opex_var=opex_var[typ],
                co2_var=co2_emissions[typ],
                regions=[region])

        ########################## T_el #################################
        try:
            capacity = float(pp.query(
                'region==@region.name and ressource==@typ and transformer=="T_el"')[
                'power'])
        except:
            capacity = 0
        if capacity > 0:
            transformer.Simple(
                uid=('transformer', region.name, typ),
                inputs=resourcebus,
                outputs=[[obj for obj in region.entities if obj.uid ==
                         "('bus', '"+region.name+"', 'elec')"][0]],
                in_max=[None],
                out_max=[capacity],
                eta=[eta_elec[typ]],
                opex_var=opex_var[typ],
                co2_var=co2_emissions[typ],
                regions=[region])

        ########################## T_heat #################################
        try:
            capacity = float(pp.query(
                'region==@region.name and ressource==@typ and transformer=="T_heat"')[
                'power'])
        except:
            capacity = 0
        if capacity > 0:
            if typ == 'powertoheat':
                transformer.Simple(
                    uid=('transformer', region.name, typ),
                    inputs=[[obj for obj in region.entities if obj.uid ==
                             "('bus', '"+region.name+"', 'elec')"][0]],
                    outputs=[[obj for obj in region.entities if obj.uid ==
                             "('bus', '"+region.name+"', 'dh')"][0]],
                    in_max=[None],
                    out_max=[capacity],
                    eta=[eta_th['heat_rod']],
                    opex_var=0,
                    co2_var=co2_emissions['heat_rod'],
                    regions=[region])
            else:
                transformer.Simple(
                    uid=('heat_transformer', region.name, typ),
                    inputs=resourcebus,
                    outputs=[[obj for obj in region.entities if obj.uid ==
                             "('bus', '"+region.name+"', 'dh')"][0]],
                    in_max=[None],
                    out_max=[capacity],
                    eta=[eta_th[typ]],
                    opex_var=opex_var[typ],
                    co2_var=co2_emissions[typ],
                    regions=[region])

        ########################## lignite LS #################################
    if region.name == 'LS':
        typ = 'lignite'
        capacity = float(pp.query(
                         'region=="LS" and ressource=="lignite_sp"')[
                         'power'])
        transformer.SimpleExtractionCHP(
            uid=('transformer', region.name, 'lignite_sp', 'SEchp'),
            inputs=[obj for obj in esystem.entities if obj.uid ==
                    "('bus', 'BB', 'lignite')"],
            outputs=[[obj for obj in region.entities if obj.uid ==
                     "('bus', 'LS', 'elec')"][0],
                    [obj for obj in region.entities if obj.uid ==
                     "('bus', 'LS', 'dh')"][0]],
            in_max=[None],
            out_max=get_out_max_chp_flex(capacity, sigma_chp['lignite']),
            out_min=[0.0, 0.0],
            eta_el_cond=eta_chp_flex_el['schwarzepumpe'],
            sigma=sigma_chp[typ],	 # power to heat ratio in backpr. mode
            beta=beta_chp[typ],		# power loss index
            opex_var=opex_var[typ],
            co2_var=co2_emissions[typ],
            regions=[region])

        capacity = float(pp.query(
                         'region=="LS" and ressource=="lignite_jw"')[
                         'power'])
        transformer.SimpleExtractionCHP(
            uid=('transformer', region.name, 'lignite_jw', 'SEchp'),
            inputs=[obj for obj in esystem.entities if obj.uid ==
                    "('bus', 'BB', 'lignite')"],
            outputs=[[obj for obj in region.entities if obj.uid ==
                     "('bus', 'LS', 'elec')"][0],
                    [obj for obj in region.entities if obj.uid ==
                     "('bus', 'LS', 'dh')"][0]],
            in_max=[None],
            out_max=get_out_max_chp_flex(capacity, sigma_chp['lignite']),
            out_min=[0.0, 0.0],
            eta_el_cond=eta_chp_flex_el['jaenschwalde'],  #TODO: softcode
            sigma=sigma_chp[typ],	 # power to heat ratio in backpr. mode
            beta=beta_chp[typ],		# power loss index
            opex_var=opex_var[typ],
            co2_var=co2_emissions[typ]*0.08,  # 92% Abscheiderate  #TODO: softcode
            regions=[region])


def get_out_max_chp(capacity, eta_th_chp, eta_el_chp):
    out_max_th = capacity * eta_th_chp / eta_el_chp
    out = [capacity, out_max_th]
    return out


def get_out_max_chp_flex(capacity, sigma_chp):
    out_max_th = capacity / sigma_chp
    out = [capacity, out_max_th]
    return out


def get_constraint_values(conn_oedb, scenario_name):
    'values for constraints in scenario'

    sql = """
        SELECT constr, val
        FROM model_draft.abbb_constraints AS d
        WHERE scenario = '""" + str(scenario_name) + """'"""
    read_parameter = pd.DataFrame(
        conn_oedb.execute(sql).fetchall(),
        columns=['constr', 'val'])

    return(read_parameter)

def add_constraint_export_minimum(om, Export_Regions, constraints):
    # returns all transport entities
    tmp_entities = [obj for obj in om.energysystem.entities
        if 'transport' in obj.uid]
    # returns all transport entities with transport to region in Export_Regions
    exports = [obj for obj in tmp_entities
        if any(region in obj.outputs[0].uid for region in Export_Regions)]
    # returns all transport entities that import from Berlin
    imports = [obj for obj in tmp_entities
        if 'BE' in obj.inputs[0].uid]
    # write list to hand over to constraint
    transports_ex = []
    for export in exports:
        transports_ex += [(export.uid, export.outputs[0].uid)]
    # write list to hand over to constraint
    transports_im = []
    for imp in imports:
        transports_im += [(imp.uid, imp.outputs[0].uid)]
    # add new constraint
    om.export_minimum_constraint = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in transports_ex for t in om.timesteps) -
        sum(om.w[i, o, t] for i, o in transports_im for t in om.timesteps)
        >= float(constraints.query('constr=="export_min"')['val'])))
    return


def add_constraint_import_berlin(om, constraints):
    # returns all transport entities
    tmp_entities = [obj for obj in om.energysystem.entities
        if 'transport' in obj.uid]
    # returns all transport entities that import from BB to Berlin
    imports = [obj for obj in tmp_entities
        if 'BE' in obj.outputs[0].uid]
    # write list to hand over to constraint
    transports_im = []
    for imp in imports:
        transports_im += [(imp.uid, imp.outputs[0].uid)]
    # add new constraint
    om.import_constraint = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in transports_im for t in om.timesteps)
        <= float(constraints.query('constr=="import_max_be"')['val'])))
    return


def add_constraint_co2_emissions(om, co2_emissions, constraints):

    # fossil ressources
    global_ressources = ['natural_gas', 'natural_gas_cc', 'lignite',
        'oil', 'waste', 'hard_coal']

    # create list of global ressource buses BB
    list_global_ressource_buses = []
    for ressource in global_ressources:
        list_global_ressource_buses += ["('bus', 'BB', '" + ressource + "')"]
    # create list with entities of global ressource buses
    global_ressource_buses_bb = [obj for obj in om.energysystem.entities
        if any(bus in obj.uid for bus in list_global_ressource_buses)]

    # create list of global ressource buses BE
    list_global_ressource_buses = []
    for ressource in global_ressources:
        list_global_ressource_buses += ["('bus', 'BE', '" + ressource + "')"]
    # create list with entities of global ressource buses
    global_ressource_buses_be = [obj for obj in om.energysystem.entities
        if any(bus in obj.uid for bus in list_global_ressource_buses)]

    # biogas
    biogas_transformer = [obj for obj in om.energysystem.entities
        if 'bhkw_bio' in obj.uid and 'transformer' in obj.uid]
    bb_regions = ['PO', 'UB', 'HF', 'OS', 'LS']
    biogas_transformer_bb = [obj for obj in biogas_transformer
        if any(region in obj.uid for region in bb_regions)]
    biogas_transformer_be = [obj for obj in biogas_transformer
        if 'BE' in obj.uid]

    # write list to hand over to BB constraint
    co2_source_bb = []
    for bus in global_ressource_buses_bb:
        for output in bus.outputs:
            co2_source_bb += [(bus.uid, output.uid, bus)]
    for transformer in biogas_transformer_bb:
        co2_source_bb += [(
            transformer.inputs[0].uid, transformer.uid, transformer.inputs[0])]

    # write list to hand over to BE constraint
    co2_source_be = []
    for bus in global_ressource_buses_be:
        for output in bus.outputs:
            co2_source_be += [(bus.uid, output.uid, bus)]
    for transformer in biogas_transformer_be:
        co2_source_be += [(
            transformer.inputs[0].uid, transformer.uid, transformer.inputs[0])]

    # add new constraint BB
    om.co2_emissions_bb = po.Constraint(expr=(
        sum(om.w[i, o, t] * co2_emissions[b.type]
        for i, o, b in co2_source_bb for t in om.timesteps)
        <= float(constraints.query('constr=="co2_max_bb"')['val'])))

    # add new constraint BE
    om.co2_emissions_be = po.Constraint(expr=(
        sum(om.w[i, o, t] * co2_emissions[b.type]
        for i, o, b in co2_source_be for t in om.timesteps)
        <= float(constraints.query('constr=="co2_max_be"')['val'])))
    return


def add_constraint_entities_BE(om):
    # Primärenergieverbrauch
    transformer_uids = {
        # Heizwerke
        "('heat_transformer', 'BE', 'oil')": 346000,
        "('heat_transformer', 'BE', 'natural_gas')": 918000,
        "('heat_transformer', 'BE', 'biomass')": 113000,
        # Kraftwerke
        "('transformer', 'BE', 'oil')": 71000,
        "('transformer', 'BE', 'natural_gas')": 2939000,
        "('transformer', 'BE', 'biomass')": 932000,
        "('transformer', 'BE', 'powertoheat')": 1005000,
        # Heizkraftwerke KWK
        "('transformer', 'BE', 'oil', 'SEchp')": 621000,
        "('transformer', 'BE', 'natural_gas_cc', 'SEchp')": 23713000,
        "('transformer', 'BE', 'biomass', 'SEchp')": 425000
        }
    BE_entities = []
    for key in list(transformer_uids.keys()):
        print(key)
        BE_entities += [obj for obj in om.energysystem.entities
            if obj.uid == key]
    # write list to hand over to constraint
    BE_transformer = []
    for entity in BE_entities:
        BE_transformer += [(entity.inputs[0].uid, entity.uid)]
    # add new constraints
    transformer = BE_transformer[0]
    om.sum_transformer_1 = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in [transformer] for t in om.timesteps)
        <= transformer_uids[transformer[1]]))
    transformer = BE_transformer[1]
    om.sum_transformer_2 = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in [transformer] for t in om.timesteps)
        <= transformer_uids[transformer[1]]))
    transformer = BE_transformer[2]
    om.sum_transformer_3 = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in [transformer] for t in om.timesteps)
        <= transformer_uids[transformer[1]]))
    transformer = BE_transformer[3]
    om.sum_transformer_4 = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in [transformer] for t in om.timesteps)
        <= transformer_uids[transformer[1]]))
    transformer = BE_transformer[4]
    om.sum_transformer_5 = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in [transformer] for t in om.timesteps)
        <= transformer_uids[transformer[1]]))
    transformer = BE_transformer[5]
    om.sum_transformer_6 = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in [transformer] for t in om.timesteps)
        <= transformer_uids[transformer[1]]))
    transformer = BE_transformer[6]
    om.sum_transformer_7 = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in [transformer] for t in om.timesteps)
        <= transformer_uids[transformer[1]]))
    transformer = BE_transformer[7]
    om.sum_transformer_8 = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in [transformer] for t in om.timesteps)
        <= transformer_uids[transformer[1]]))
    transformer = BE_transformer[8]
    om.sum_transformer_9 = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in [transformer] for t in om.timesteps)
        <= transformer_uids[transformer[1]]))
    transformer = BE_transformer[9]
    om.sum_transformer_a = po.Constraint(expr=(
        sum(om.w[i, o, t] for i, o in [transformer] for t in om.timesteps)
        <= transformer_uids[transformer[1]]))
    return