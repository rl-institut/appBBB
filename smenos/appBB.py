#!/usr/bin/python3
# -*- coding: utf-8

import pandas as pd
import warnings
import logging
import os

from oemof import db
from oemof.tools import logger
from oemof.core import energy_system as es
from oemof.db import tools
#from oemof.db import feedin_pg
from oemof.solph import predefined_objectives as predefined_objectives
from oemof.core.network.entities import Bus
from oemof.core.network.entities.components import sinks as sink
from oemof.core.network.entities.components import transports as transport
from oemof.core.network.entities.components import sources as source
from oemof.solph.optimization_model import OptimizationModel

import helper_SmEnOs as hls
import helper_BBB as hlsb
import helper_dec_BBB as hlsd
import numpy as np
from workalendar.europe import Germany

# choose scenario
scenario = 'gruene2030_bk0'

# Basic inputs
warnings.simplefilter(action="ignore", category=RuntimeWarning)
logger.define_logging()
year = 2010
time_index = pd.date_range('1/1/{0}'.format(year), periods=8760, freq='H')
time_index_demandlib = pd.date_range(
    '1/1/{0}'.format(year), periods=8760, freq='H')
conn = db.connection()
conn_oedb = db.connection(section='open_edb')

########### get data ###########################################

cap_initial = 0.0
chp_faktor_flex = 0.84  # share of flexible generation of CHP
max_biomass = 16111111  # 58 PJ per year
energy_emob_BB = 2183000  # 7,86 PJ per year (10% of whole traffic)
energy_emob_BE = 1298333
share_emob = {}
share_emob['PO'] = 0.16
share_emob['UB'] = 0.12
share_emob['HF'] = 0.29
share_emob['OS'] = 0.17
share_emob['LS'] = 0.25

# parameters
(co2_emissions, co2_fix, eta_elec, eta_th, eta_th_chp, eta_el_chp,
 eta_chp_flex_el, sigma_chp, beta_chp, opex_var, opex_fix, capex,
 c_rate_in, c_rate_out, eta_in, eta_out,
 cap_loss, lifetime, wacc) = hlsb.get_parameters(conn_oedb)

transmission = hlsb.get_transmission(conn_oedb, scenario)
demands_df = hlsb.get_demand(conn_oedb, scenario)
transformer = hlsb.get_transformer(conn_oedb, scenario)
# st = hlsb.get_st_timeline(conn, year)  # timeline for solar heat
cal = Germany()
holidays = dict(cal.holidays(2010))
    
############## Create a simulation object ########################
simulation = es.Simulation(
    timesteps=list(range(len(time_index))), verbose=True, solver='cbc',
    stream_solver_output=True,
    objective_options={'function': predefined_objectives.minimize_cost})

############## Create an energy system ###########################
Regions = es.EnergySystem(time_idx=time_index, simulation=simulation)

regionsBBB = pd.DataFrame(
    [{'abbr': 'PO', 'nutsID': ['DE40F', 'DE40D', 'DE40A']},
        {'abbr': 'UB', 'nutsID': ['DE40I', 'DE405']},
        {'abbr': 'HF', 'nutsID': [
            'DE408', 'DE40E', 'DE40H', 'DE401', 'DE404']},
        {'abbr': 'OS', 'nutsID': ['DE409', 'DE40C', 'DE403']},
        {'abbr': 'LS', 'nutsID': [
            'DE406', 'DE407', 'DE40B', 'DE40G', 'DE402']},
        {'abbr': 'BE', 'nutsID': 'DE3'}],
    index=['Prignitz-Oberhavel', 'Uckermark-Barnim', u'Havelland-Fläming',
           'Oderland-Spree', 'Lausitz-Spreewald', 'Berlin'])

for index, row in regionsBBB.iterrows():
    Regions.regions.append(es.Region(
        geom=tools.get_polygon_from_nuts(conn, row['nutsID']),
        name=row['abbr']))

region_bb = []
for region in Regions.regions:
    if region.name == 'BE':
        region_ber = region
    else:
        region_bb.append(region)  # list

filename = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                        'data_mob_bev.csv'))

emob_DE = pd.read_csv(filename, delimiter=',')
emob_BB = emob_DE['sink_fixed'] * energy_emob_BB / emob_DE['sink_fixed'].sum()
emob_BE = emob_DE['sink_fixed'] * energy_emob_BE / emob_DE['sink_fixed'].sum()

# Add electricity sink and bus for each region
for region in Regions.regions:
    # create electricity bus
    Bus(uid="('bus', '"+region.name+"', 'elec')", type='elec', price=0,
        regions=[region], excess=True, shortage=True, shortage_costs=1000000.0)

    # create districtheat bus
    Bus(uid="('bus', '"+region.name+"', 'dh')", type='dh', price=0,
        regions=[region], excess=True)

    # create electricity sink
    demand = sink.Simple(uid=('demand', region.name, 'elec'),
                         inputs=[obj for obj in region.entities if obj.uid ==
                                 "('bus', '"+region.name+"', 'elec')"],
                         regions=[region])
    el_demands = {}
    el_demands['h0'] = float(demands_df.query(
        'region==@region.name and sector=="HH" and type=="electricity"')[
        'demand'])
    el_demands['g0'] = float(demands_df.query(
        'region==@region.name and sector=="GHD" and type=="electricity"')[
        'demand'])
    el_demands['i0'] = float(demands_df.query(
        'region==@region.name and sector=="IND" and type=="electricity"')[
        'demand'])
    am, pm, profile_factors = hls.ind_profile_parameters()
    hls.el_load_profiles(demand, el_demands, year, holidays=holidays,
                         am=am, pm=pm, profile_factors=profile_factors)
    if region.name != 'BE':
        demand = sink.Simple(uid=('demand', region.name, 'elec', 'mob'),
                         inputs=[obj for obj in region.entities if obj.uid ==
                                 "('bus', '"+region.name+"', 'elec')"],
                         regions=[region])
        emob = emob_BB * share_emob[region.name]
        demand.val = emob
    else:
        demand = sink.Simple(uid=('demand', region.name, 'elec', 'mob'),
                         inputs=[obj for obj in region.entities if obj.uid ==
                                 "('bus', '"+region.name+"', 'elec')"],
                         regions=[region])
        demand.val = emob_BE

# Add global buses for BB and BE
typeofgen_global = ['natural_gas', 'natural_gas_cc', 'lignite',
                    'oil', 'waste', 'hard_coal']
# Add biomass bus for Berlin and Brandenburg
for typ in typeofgen_global:
    Bus(uid="('bus', 'BB', '"+typ+"')",
        type=typ, price=0,
        excess=False, balanced=False,
        co2_var=co2_emissions[typ],
        regions=Regions.regions)
    Bus(uid="('bus', 'BE', '"+typ+"')",
        type=typ, price=0,
        excess=False, balanced=False,
        co2_var=co2_emissions[typ],
        regions=Regions.regions)

Bus(uid="('bus', 'BB', 'biomass')",
    type='biomass',
    price=0,
    balanced=False,
    sum_out_limit=max_biomass,
    co2_var=co2_emissions['biomass'],
    regions=region_bb,
    excess=False)

Bus(uid="('bus', 'BE', 'biomass')",
    type='biomass',
    price=0,
    balanced=False,
    co2_var=co2_emissions['biomass'],
    regions=[region_ber],
    excess=False)
print("('bus', 'BE', 'biomass')")

################# create transformers ######################
                ########### decentral #####################
hlsd.create_decentral_entities(Regions, regionsBBB, demands_df, conn, year,
                               time_index_demandlib, eta_th, eta_in, eta_out,
                               cap_loss,
                               opex_fix, opex_var, eta_th_chp, eta_el_chp,
                               holidays)

# renewable parameters
site = hls.get_res_parameters()

# add biomass and powertoheat in typeofgen
# because its needed to generate powerplants from db
typeofgen_global.append('biomass')
typeofgen_global.append('powertoheat')
typeofgen_global.append('bhkw_bio')

for region in Regions.regions:
    logging.info('Processing region: {0}'.format(region.name))

                ########### central #####################
    hlsb.create_transformer(
        Regions, region, transformer, conn=conn_oedb,
        cap_initial=cap_initial,
        chp_faktor_flex=chp_faktor_flex,  # share of flexible generation of CHP
        typeofgen=typeofgen_global)

    #TODO Problem mit Erdwärme??!!
                ########### renewables #####################
    filename = os.path.abspath(os.path.join(
                os.path.dirname(__file__),
                  'res_timeseries_' + region.name + '.csv'))
    feedin_df = pd.read_csv(filename, delimiter=',', index_col=0)
#    feedin_df, cap = feedin_pg.Feedin().aggregate_cap_val(
#        conn, region=region, year=year, bustype='elec', **site)
    ee_capacities = {}
    ee_capacities['pv_pwr'] = float(transformer.query(
        'region==@region.name and ressource=="pv"')['power'])
    ee_capacities['wind_pwr'] = float(transformer.query(
        'region==@region.name and ressource=="wind"')['power'])

    for stype in feedin_df.keys():
        source.FixedSource(
            uid=('FixedSrc', region.name, stype),
            outputs=[obj for obj in region.entities if obj.uid ==
                     "('bus', '"+region.name+"', 'elec')"],
            val=feedin_df[stype],
            out_max=[ee_capacities[stype]])

# Remove orphan buses
buses = [obj for obj in Regions.entities if isinstance(obj, Bus)]
for bus in buses:
    if len(bus.inputs) > 0 or len(bus.outputs) > 0:
        logging.debug('Bus {0} has connections.'.format(bus.type))
    else:
        logging.debug('Bus {0} has no connections and will be deleted.'.format(
            bus.type))
        Regions.entities.remove(bus)

Import_Regions = ('MV', 'ST', 'SN', 'KJ')
Export_Regions = ('MV', 'ST', 'SN', 'KJ', 'BE')

for region_name in Import_Regions:
    Regions.regions.append(es.Region(name=region_name))

for region in Regions.regions:
    if region.name in Import_Regions:
        Bus(uid="('bus', '"+region.name+"', 'elec')", type='elec',
            regions=[region],
            excess=True)
        Bus(uid="('bus', '"+region.name+"', 'import')", type='elec',
            shortage=True,
            shortage_costs=opex_var['import_el'],
            regions=[region])

## print all entities of every region
#for entity in Regions.entities:
#    print(entity.uid)
#    if entity.uid[0] == 'transformer' or entity.uid[0] == 'FixedSrc':
#        print('out_max')
#        print(entity.out_max)
#        print('type(out_max)')
#        print(type(entity.out_max))

# Connect the electrical buses of federal states

for con in transmission['from']:  # Zeilen in transmission-Tabelle
    capacity = transmission['cap'][con]
    if capacity != 0:
        reg1 = transmission['from'][con]  # zeile x,Spalte 'from'
        reg2 = transmission['to'][con]  # zeile x,Spalte 'from'
        eta_trans = eta_elec['transmission']
        for entity in Regions.entities:
            if entity.uid == "('bus', '" + reg1 + "', 'elec')":
                ebus_1 = entity
            if entity.uid == "('bus', '" + reg2 + "', 'elec')":
                ebus_2 = entity
            if entity.uid == "('bus', '" + reg2 + "', 'import')":
                ebus_import = entity
        # if reg2 is not BB or BE only create transport into the region
        if reg2 in Import_Regions:
            # if transport to MV or ST, limit export in times of wind
            # feedin > 0.9
            if reg2 in ['MV', 'ST']:
                filename = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                            'res_timeseries_' + reg1 + '.csv'))
                wind_feedin = pd.read_csv(filename)
                wind_feedin['transport_condition'] = np.where(
                    wind_feedin['wind_pwr'] >= 0.9, 0, 1)
                transport.Simple(  # export-connection MV/ST
                    uid='transport_' + ebus_1.uid + ebus_2.uid,
                    inputs=[ebus_1], outputs=[ebus_2],
                    out_max=[capacity], in_max=[capacity], eta=[eta_trans],
                    ub_out=[wind_feedin['transport_condition'] * capacity])
                transport.Simple(  # import-connection MV/ST
                    uid='transport_' + ebus_import.uid + ebus_1.uid,
                    outputs=[ebus_1], inputs=[ebus_import],
                    out_max=[capacity], in_max=[capacity], eta=[eta_trans])
            else:
                transport.Simple(  # export-connection SN/KJ/
                    uid='transport_' + ebus_1.uid + ebus_2.uid,
                    inputs=[ebus_1], outputs=[ebus_2],
                    out_max=[capacity], in_max=[capacity], eta=[eta_trans])
                transport.Simple(  # import-connection SN/KJ
                    uid='transport_' + ebus_import.uid + ebus_1.uid,
                    outputs=[ebus_1], inputs=[ebus_import],
                    out_max=[capacity], in_max=[capacity], eta=[eta_trans])
        else:
            Regions.connect(ebus_1, ebus_2,
                            in_max=capacity,
                            out_max=capacity,
                            eta=eta_trans,
                            transport_class=transport.Simple)

# change uid tuples to strings
for entity in Regions.entities:
    entity.uid = str(entity.uid)
    print(entity.uid)

# Optimize the energy system
om = OptimizationModel(energysystem=Regions)

constraints = hlsb.get_constraint_values(conn_oedb, scenario)

hlsb.add_constraint_export_minimum(om, Export_Regions, constraints)
hlsb.add_constraint_import_berlin(om, constraints)
hlsb.add_constraint_co2_emissions(om, co2_emissions, constraints)
hlsb.add_constraint_entities_BE(om)
om.write_lp_file()
om.solve()
Regions.results = om.results()
logging.info(Regions.dump())
