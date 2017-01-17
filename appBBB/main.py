import pandas as pd
import numpy as np
import warnings
import logging
import os
from workalendar.europe import Germany

from oemof import db
from oemof.tools import logger
from oemof.core import energy_system as es
from oemof.solph import predefined_objectives as predefined_objectives
from oemof.core.network.entities import Bus
from oemof.core.network.entities.components import sinks as sink
from oemof.core.network.entities.components import transports as transport
from oemof.core.network.entities.components import sources as source
from oemof.solph.optimization_model import OptimizationModel

import helper_BBB as hlsb
import helper_dec_BBB as hlsd


################################# CHOOSE SCENARIO ############################                                                         
scenario = 'ES2030' #TODO Welche anderen Szenarien können gewählt werden?
##############################################################################

################################# BASIC SETTINGS #############################
# set logging
warnings.simplefilter(action="ignore", category=RuntimeWarning)
logger.define_logging()

# establish database connections
conn_oedb = db.connection(section='open_edb')

# set solver
solver = 'cbc'
##############################################################################

################################# GET/SET DATA ###############################
# create time indexes
year = 2010
time_index = pd.date_range('1/1/{0}'.format(year), periods=2, freq='H')
time_index_demandlib = pd.date_range(
    '1/1/{0}'.format(year), periods=8760, freq='H')
# get German holidays
cal = Germany()
holidays = dict(cal.holidays(year))

# set regions to be considered along with their nuts ID and abbreviation
regionsBBB = pd.DataFrame(
    [{'abbr': 'PO', 'nutsID': ['DE40F', 'DE40D', 'DE40A']},
     {'abbr': 'UB', 'nutsID': ['DE40I', 'DE405']},
     {'abbr': 'HF', 'nutsID': [
            'DE408', 'DE40E', 'DE40H', 'DE401', 'DE404']},
     {'abbr': 'OS', 'nutsID': ['DE409', 'DE40C', 'DE403']},
     {'abbr': 'LS', 'nutsID': [
            'DE406', 'DE407', 'DE40B', 'DE40G', 'DE402']},
     {'abbr': 'BE', 'nutsID': 'DE300'}],
    index=['Prignitz-Oberhavel', 'Uckermark-Barnim', u'Havelland-Fläming',
           'Oderland-Spree', 'Lausitz-Spreewald', 'Berlin'])
           
# set maximum biomass availability for Brandenburg
maximum_biomass_availability = 16111111  # 58 PJ/a

# get electricity demand for electromobility for Brandenburg and Berlin 
# from database
emob_energy = hlsb.get_emob_values(conn_oedb, scenario)
energy_emob_BB = float(emob_energy.query('region=="BB"')['energy'])
energy_emob_BE = float(emob_energy.query('region=="BE"')['energy'])
# set shares of electromobility
share_emob = {}
share_emob['PO'] = 0.16
share_emob['UB'] = 0.12
share_emob['HF'] = 0.29
share_emob['OS'] = 0.17
share_emob['LS'] = 0.25
# get electromobility timeseries from data directory
emob_data = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         'data', 'data_mob_bev.csv'))
emob_DE = pd.read_csv(emob_data, delimiter=',')
                                       
# get powerplant parameters from database
(co2_emissions, co2_fix, eta_elec, eta_th, eta_th_chp, eta_el_chp,
 eta_chp_flex_el, sigma_chp, beta_chp, opex_var, opex_fix, capex,
 c_rate_in, c_rate_out, eta_in, eta_out,
 cap_loss, lifetime, wacc) = hlsb.get_parameters(conn_oedb)

# get transmission capacities from database
transmission = hlsb.get_transmission(conn_oedb, scenario)

# get demand timeseries from database
demands_df = hlsb.get_demand(conn_oedb, scenario)

# get powerplant capacities of Berlin and Brandenburg
transformer_capacities = hlsb.get_transformer(conn_oedb, scenario)

# get parameters for industrial load profiles
am, pm, profile_factors = hlsb.ind_profile_parameters()
##############################################################################

######################### INITIALISE OEMOF ENERGY SYSTEM #####################
logging.info('Initialising oemof energy system')                                   
# create oemof simulation object
simulation = es.Simulation(
    timesteps=list(range(len(time_index))), verbose=True, solver=solver,
    stream_solver_output=True,
    objective_options={'function': predefined_objectives.minimize_cost})

# create oemof energy system object
Regions = es.EnergySystem(time_idx=time_index, simulation=simulation)
          
# append regions to EnergySystem object
for index, row in regionsBBB.iterrows():
    Regions.regions.append(es.Region(
        geom=hlsb.get_polygon_from_nuts(conn_oedb, row['nutsID']),
        name=row['abbr']))
# create lists with Region objects of all regions in Berlin and Brandenburg
region_ber = []
region_bb = []
for region in Regions.regions:
    if region.name == 'BE':
        region_ber.append(region)
    else:
        region_bb.append(region)
##############################################################################

###################### ADD ELECTRICITY BUSES AND SINKS #######################
logging.info('Adding electricity buses and sinks')
for region in Regions.regions:
    # create electricity bus
    Bus(uid="('bus', '" + region.name + "', 'elec')",
        type='elec',
        regions=[region],
        excess=True,
        shortage=True,
        shortage_costs=1000000.0) # randomly high shortage costs to avoid
                                  # shortage
    # create electricity sink
    demand = sink.Simple(
        uid=('demand', region.name, 'elec'),
        inputs=[obj for obj in region.entities if obj.uid ==
            "('bus', '"+region.name+"', 'elec')"],
        regions=[region]) 
    # get electricity demands for the residential, commercial and
    # industrial sectors              
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
    # generate load profiles for all three sectors and write aggregated
    # timeseries to sink object    
    hlsb.el_load_profiles(demand, el_demands, year, holidays=holidays,
                          am=am, pm=pm, profile_factors=profile_factors)
##############################################################################

##################### ADD ELECTROMOBILITY BUS AND SINK #######################
logging.info('Adding electromobility bus and sink')
for region in Regions.regions: 
    demand = sink.Simple(uid=('demand', region.name, 'elec', 'mob'),
                         inputs=[obj for obj in region.entities if obj.uid ==
                            "('bus', '" + region.name + "', 'elec')"],
                         regions=[region])          
    if region.name != 'BE':      
        demand.val = (
            (emob_DE['sink_fixed'] * energy_emob_BB / 
                emob_DE['sink_fixed'].sum()) *
            share_emob[region.name])
    else:
        demand.val = (
            emob_DE['sink_fixed'] * energy_emob_BE / 
                emob_DE['sink_fixed'].sum())
##############################################################################

######################## ADD GLOBAL RESSOURCE BUSES ##########################
logging.info('Adding global ressource buses')
# fossil fuels and waste
typeofgen_global = ['natural_gas', 'natural_gas_cc', 'lignite',
                    'oil', 'waste', 'hard_coal']
for typ in typeofgen_global:
    # global ressource bus
    ressource_bus = Bus(
        uid="('bus', 'source', '" + typ + "')",
        type=typ,
        shortage=True,
        shortage_costs=opex_var[typ],
        excess=False,
        regions=Regions.regions)
    # ressource bus Brandenburg
    bus_brandenburg = Bus(
        uid="('bus', 'BB', '" + typ + "')",
        type=typ,
        co2_var=co2_emissions[typ],
        excess=False,
        regions=Regions.regions)
    # ressource bus Berlin
    bus_berlin = Bus(
        uid="('bus', 'BE', '" + typ + "')",
        type=typ,
        co2_var=co2_emissions[typ],
        excess=False,
        regions=Regions.regions)
    # transport from global ressource bus to ressource bus of Brandenburg
    transport.Simple(  
        uid='transport_ressource_' + typ,
        outputs=[bus_brandenburg],
        inputs=[ressource_bus],
        out_max=[1000000],
        in_max=[1000000],
        eta=[1.0])
    # transport from global ressource bus to ressource bus of Berlin
    transport.Simple(
        uid='transport_ressource_BE' + typ,
        outputs=[bus_berlin], 
        inputs=[ressource_bus],
        out_max=[1000000],
        in_max=[1000000],
        eta=[1.0])

# global biomass bus
Bus_bio = Bus(
    uid="('bus', 'source', 'biomass')",
    type='biomass',
    shortage=True,
    shortage_costs=opex_var['biomass'],
    co2_var=co2_emissions['biomass'],
    regions=region_bb,
    excess=False)
# biomass bus Brandenburg
BusBB_Bio = Bus(
    uid="('bus', 'BB', 'biomass')",
    type='biomass',
    sum_out_limit=maximum_biomass_availability,
    co2_var=co2_emissions['biomass'],
    regions=region_bb,
    excess=False)
# biomass bus Berlin
BusBE_Bio = Bus(
    uid="('bus', 'BE', 'biomass')",
    type='biomass',
    co2_var=co2_emissions['biomass'],
    regions=region_ber,
    excess=False)   
# transport from global biomass bus to biomass bus Brandenburg
transport.Simple(  
    uid='transport_ressource_biomass',
    outputs=[BusBB_Bio],
    inputs=[Bus_bio],
    out_max=[1000000],
    in_max=[1000000],
    eta=[1.0])
# transport from global biomass bus to biomass bus Berlin
transport.Simple(  
    uid='transport_ressource_biomass_BE',
    outputs=[BusBE_Bio],
    inputs=[Bus_bio],
    out_max=[1000000],
    in_max=[1000000],
    eta=[1.0])
##############################################################################
    
######################## CREATE DISTRICT HEATING BUSES #######################
logging.info('Adding district heating buses')
for region in Regions.regions:
    Bus(uid="('bus', '" + region.name + "', 'dh')",
        type='dh',
        regions=[region],
        excess=True)
##############################################################################
        
###################### CREATE DECENTRAL HEATING ENTITIES #####################
logging.info('Creating decentral heating entities')
hlsd.create_decentral_entities(Regions, regionsBBB, demands_df, year,
                               time_index_demandlib, eta_th, eta_in, eta_out,
                               cap_loss, opex_fix, opex_var, eta_th_chp,
                               eta_el_chp, holidays)
##############################################################################
                              
######################### CREATE POWER PLANTS ################################
# add biomass, powertoheat and bhkw_bio to list typeofgen_global to generate
# powerplants from database in hlsb.create_transformer function
typeofgen_global.append('biomass')
typeofgen_global.append('powertoheat')
typeofgen_global.append('bhkw_bio')

for region in Regions.regions:
    # PV + Wind
      # load feed-in timeseries data
    filename = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), 'data',
            'res_timeseries_' + region.name + '.csv'))
    feedin_df = pd.read_csv(filename, delimiter=',', index_col=0)
      # get capacities of PV and wind from transformer capacity dataframe
    ee_capacities = {}
    ee_capacities['pv_pwr'] = float(transformer_capacities.query(
        'region==@region.name and ressource=="pv"')['power'])
    ee_capacities['wind_pwr'] = float(transformer_capacities.query(
        'region==@region.name and ressource=="wind"')['power'])
      # create fixed source objects for PV and wind
    for stype in feedin_df.keys():
        source.FixedSource(
            uid=('FixedSrc', region.name, stype),
            outputs=[obj for obj in region.entities if obj.uid ==
                     "('bus', '" + region.name + "', 'elec')"],
            val=feedin_df[stype],
            out_max=[ee_capacities[stype]])
    
    # other power plants
    hlsb.create_transformer(
        Regions, region, transformer_capacities, conn=conn_oedb,
        typeofgen=typeofgen_global)
##############################################################################

########################### REMOVE ORPHAN BUSES ##############################
# get all buses
buses = [obj for obj in Regions.entities if isinstance(obj, Bus)]
for bus in buses:
    if len(bus.inputs) > 0 or len(bus.outputs) > 0:
        logging.debug('Bus {0} has connections.'.format(bus.type))
    else:
        logging.debug('Bus {0} has no connections and will be deleted.'.format(
            bus.type))
        Regions.entities.remove(bus)
##############################################################################

################ CREATE ELECTRICITY BUSES FOR IMPORT REGIONS #################
# import-regions are regions where electricity can be imported from to 
# Brandenburg (Mecklenburg-Vorpommern, Sachsen-Anhalt, Sachsen, Polen)
Import_Regions = ('MV', 'ST', 'SN', 'KJ')

# append to Regions instance
for region_name in Import_Regions:
    Regions.regions.append(es.Region(name=region_name))

# create separate electricity buses for import and export
for region in Regions.regions:
    if region.name in Import_Regions:
        Bus(uid="('bus', '" + region.name + "', 'export')",
            type='elec',
            regions=[region],
            excess=True)
        Bus(uid="('bus', '" + region.name + "', 'import')",
            type='elec',
            shortage=True,
            shortage_costs=opex_var['import_el'],
            regions=[region])
##############################################################################
            
####################### CONNECT ELECTRICITY BUSES ############################
# iterate through transmission lines
# entry is the number of the transmission line
for entry in transmission['from']:
    capacity = transmission['cap'][entry]
    if capacity != 0:
        # get from and to region
        from_reg = transmission['from'][entry]
        to_reg = transmission['to'][entry]
        eta_trans = eta_elec['transmission']
        
        # get electricity bus entities
        for entity in Regions.entities:
            if entity.uid == "('bus', '" + from_reg + "', 'elec')":
                ebus_from = entity
            if (entity.uid == "('bus', '" + to_reg + "', 'export')" or
                entity.uid == "('bus', '" + to_reg + "', 'elec')"):
                ebus_export = entity
            if entity.uid == "('bus', '" + to_reg + "', 'import')":
                ebus_import = entity
                
        # if to_reg is an import region create different transmissions for
        # import and export electricity buses
        if to_reg in Import_Regions:
            
            # if transport to MV or ST, limit export in times of wind
            # feedin > 0.9
            if to_reg in ['MV', 'ST']:               
                # get wind feed-in timeseries
                filename = os.path.abspath(
                    os.path.join(
                        os.path.dirname(__file__), 'data',
                        'res_timeseries_' + from_reg + '.csv'))
                wind_feedin = pd.read_csv(filename)
                # in case of wind power feed-in > 0.9 set import to BB to zero
                wind_feedin['transport_condition'] = np.where(
                    wind_feedin['wind_pwr'] >= 0.9, 0, 1)
                    
                # export connection
                transport.Simple(  
                    uid='transport_' + ebus_from.uid + ebus_export.uid,
                    inputs=[ebus_from],
                    outputs=[ebus_export],
                    out_max=[capacity],
                    in_max=[capacity],
                    eta=[eta_trans],
                    ub_out=[wind_feedin['transport_condition'] * capacity])
                # import connection
                transport.Simple(
                    uid='transport_' + ebus_import.uid + ebus_from.uid,
                    inputs=[ebus_import],
                    outputs=[ebus_from],
                    out_max=[capacity],
                    in_max=[capacity],
                    eta=[eta_trans])
            
            # if transport to SN or KJ
            else:
                # export connection
                transport.Simple(  
                    uid='transport_' + ebus_from.uid + ebus_export.uid,
                    inputs=[ebus_from],
                    outputs=[ebus_export],
                    out_max=[capacity],
                    in_max=[capacity],
                    eta=[eta_trans])
                # import connection
                transport.Simple( 
                    uid='transport_' + ebus_import.uid + ebus_from.uid,
                    inputs=[ebus_import],
                    outputs=[ebus_from], 
                    out_max=[capacity],
                    in_max=[capacity],
                    eta=[eta_trans])
                    
        # if to_reg is not an import region
        else:
            Regions.connect(ebus_from, ebus_export,
                            in_max=capacity,
                            out_max=capacity,
                            eta=eta_trans,
                            transport_class=transport.Simple)
##############################################################################

############################### FINAL STEPS ##################################
# change uid tuples to strings
for entity in Regions.entities:
    entity.uid = str(entity.uid)

# create optimisation model from Regions instance
om = OptimizationModel(energysystem=Regions)

# add constraints
constraints = hlsb.get_constraint_values(conn_oedb, scenario)
Export_Regions = ('MV', 'ST', 'SN', 'KJ', 'BE')
hlsb.add_constraint_export_minimum(om, Export_Regions, constraints)
hlsb.add_constraint_import_berlin(om, constraints)
hlsb.add_constraint_co2_emissions(om, co2_emissions, constraints)
hlsb.add_constraint_entities_BE(om)

# write LP file
om.write_lp_file()
# solve problem
om.solve()
# save results to Regions instance
Regions.results = om.results()
# dump Regions instance
logging.info(Regions.dump())
##############################################################################