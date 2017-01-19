from os.path import expanduser
import csv
import pandas as pd
import matplotlib.pyplot as plt
from oemof.core import energy_system as es
from oemof.solph.predefined_objectives import minimize_cost
from oemof.outputlib import to_pandas as tpd
from oemof import db
import helper_BBB as hlsb


def create_es(solver, timesteps, year):
    """ 
    Creates a default energy system to load results into.
    """
    simulation = es.Simulation(solver=solver, 
                               timesteps=timesteps,
                               debug=False, 
                               objective_options={"function": minimize_cost})

    # Adding a time index to the energy system
    time_index = pd.date_range('1/1/' + year,
                               periods=len(timesteps),
                               freq='H')
    energysystem = es.EnergySystem(time_idx=time_index,
                                   simulation=simulation)
    return energysystem


def color_dict(reg):
    """
    Sets colors for entities in plot of electricity sector.
    """
    cdict = {
             # transformer
             "('FixedSrc', '" + reg + "', 'wind_pwr')": 'lightblue',
             "('FixedSrc', '" + reg + "', 'pv_pwr')": 'yellow',
 
             "('transformer', '" + reg + "', 'oil')": 'black',
             "('transformer', '" + reg + "', 'oil', 'chp')": 'black',
             "('transformer', '" + reg + "', 'oil', 'SEchp')": 'black',

             "('transformer', '" + reg + "', 'natural_gas')": 'lightgrey',
             "('transformer', '" + reg + "', 'natural_gas', 'chp')":
                 'lightgrey',
             "('transformer', '" + reg + "', 'natural_gas', 'SEchp')":
                 'lightgrey',

             "('transformer', '" + reg + "', 'natural_gas_cc')": 'darkgrey',
             "('transformer', '" + reg + "', 'natural_gas_cc', 'chp')":
                 'darkgrey',
             "('transformer', '" + reg + "', 'natural_gas_cc', 'SEchp')":
                 'darkgrey',

             "('transformer', '" + reg + "', 'HH', 'bhkw_gas')": 'grey',
             "('transformer', '" + reg + "', 'GHD', 'bhkw_gas')": 'grey',

             "('transformer', '" + reg + "', 'biomass')": 'lightgreen',
             "('transformer', '" + reg + "', 'biomass', 'chp')": 'lightgreen',
             "('transformer', '" + reg + "', 'biomass', 'SEchp')":
                 'lightgreen',

             "('transformer', '" + reg + "', 'HH', 'bhkw_bio')": 'green',
             "('transformer', '" + reg + "', 'GHD', 'bhkw_bio')": 'green',

             "('transformer', '" + reg + "', 'powertoheat')": 'lightsalmon',

             "('transformer', '" + reg + "', 'lignite_jw', 'SEchp')": 'brown',
             "('transformer', '" + reg + "', 'lignite_sp', 'SEchp')": 'orange',

             # demand
             "('demand', '" + reg + "', 'elec')": 'red',
             "('demand', '" + reg + "', 'elec', 'mob')": 'red',

             # shortage / excess
             "('bus', '" + reg + "', 'elec')_excess": 'purple',
             "('bus', '" + reg + "', 'elec')_shortage": 'blueviolet',

             # heat pump
             "('transformer', '" + reg + "', 'hp', 'brine', 'ww')": 'blue',
             "('transformer', '" + reg + "', 'hp', 'brine', 'heating')":
                 'blue',
             "('transformer', '" + reg + "', 'hp', 'air', 'ww')": 'blue',
             "('transformer', '" + reg + "', 'hp', 'air', 'heating')": 'blue',
             "('transformer', '" + reg + "', 'hp', 'air', 'ww', 'rod')":
                 'blue',
             "('transformer', '" + reg + "', 'hp', 'air', 'heating', 'rod')":
                 'blue',

             # transport
             "transport_('bus', 'UB', 'elec')('bus', 'OS', 'elec')": 'salmon',
             "transport_('bus', 'OS', 'elec')('bus', 'UB', 'elec')": 'salmon',
             "transport_('bus', 'OS', 'elec')('bus', 'LS', 'elec')":
                 'chocolate',
             "transport_('bus', 'LS', 'elec')('bus', 'OS', 'elec')":
                 'chocolate',
             "transport_('bus', 'OS', 'elec')('bus', 'BE', 'elec')": 'peru',
             "transport_('bus', 'BE', 'elec')('bus', 'OS', 'elec')": 'peru',
             "transport_('bus', 'LS', 'elec')('bus', 'HF', 'elec')":
                 'burlywood',
             "transport_('bus', 'HF', 'elec')('bus', 'LS', 'elec')":
                 'burlywood',
             "transport_('bus', 'HF', 'elec')('bus', 'PO', 'elec')":
                 'goldenrod',
             "transport_('bus', 'PO', 'elec')('bus', 'HF', 'elec')":
                 'goldenrod',
             "transport_('bus', 'HF', 'elec')('bus', 'BE', 'elec')": 'khaki',
             "transport_('bus', 'BE', 'elec')('bus', 'HF', 'elec')": 'khaki',
             "transport_('bus', 'PO', 'elec')('bus', 'OS', 'elec')":
                 'indianred',
             "transport_('bus', 'OS', 'elec')('bus', 'PO', 'elec')":
                 'indianred',
             "transport_('bus', 'UB', 'elec')('bus', 'KJ', 'elec')": 'lime',
             "transport_('bus', 'UB', 'elec')('bus', 'MV', 'elec')": 'cyan',
             "transport_('bus', 'PO', 'elec')('bus', 'MV', 'elec')": 'teal',
             "transport_('bus', 'PO', 'elec')('bus', 'ST', 'elec')":
                 'seagreen',
             "transport_('bus', 'HF', 'elec')('bus', 'ST', 'elec')":
                 'yellowgreen',
             "transport_('bus', 'LS', 'elec')('bus', 'SN', 'elec')":
                 'turquoise',
             "transport_('bus', 'BE', 'elec')('bus', 'HF', 'elec')": 'olive',
             "transport_('bus', 'BE', 'elec')('bus', 'OS', 'elec')":
                 'lightseagreen',
             "transport_('bus', 'KJ', 'import')('bus', 'UB', 'elec')": 'lime',
             "transport_('bus', 'MV', 'import')('bus', 'UB', 'elec')": 'cyan',
             "transport_('bus', 'MV', 'import')('bus', 'PO', 'elec')": 'teal',
             "transport_('bus', 'ST', 'import')('bus', 'PO', 'elec')":
                 'seagreen',
             "transport_('bus', 'ST', 'import')('bus', 'HF', 'elec')":
                 'yellowgreen',
             "transport_('bus', 'SN', 'import')('bus', 'LS', 'elec')":
                 'turquoise',
             "transport_('bus', 'HF', 'elec')('bus', 'BE', 'elec')": 'olive',
             "transport_('bus', 'OS', 'elec')('bus', 'BE', 'elec')":
                 'lightseagreen'}
    return cdict

def color_dict_dh(reg):
    """
    Sets colors for entities in plot of district heating.
    """
    cdict = {
             # transformer
             "('transformer', '" + reg + "', 'oil', 'chp')": 'black',
             "('transformer', '" + reg + "', 'oil', 'SEchp')": 'black',
             "('heat_transformer', '" + reg + "', 'oil')": 'black',

             "('transformer', '" + reg + "', 'natural_gas', 'chp')":
                 'lightgrey',
             "('transformer', '" + reg + "', 'natural_gas', 'SEchp')":
                 'lightgrey',
             "('heat_transformer', '" + reg + "', 'natural_gas')":
                 'lightgrey',
             "('transformer', '" + reg + "', 'dh_peak_heating')": 'khaki',

             "('transformer', '" + reg + "', 'natural_gas_cc', 'chp')":
                 'darkgrey',
             "('transformer', '" + reg + "', 'natural_gas_cc', 'SEchp')":
                 'darkgrey',

             "('transformer', '" + reg + "', 'biomass', 'chp')": 'lightgreen',
             "('transformer', '" + reg + "', 'biomass', 'SEchp')":
                 'lightgreen',
             "('heat_transformer', '" + reg + "', 'biomass')": 'lightgreen',

             "('transformer', '" + reg + "', 'lignite_jw', 'SEchp')": 'brown',
             "('transformer', '" + reg + "', 'lignite_sp', 'SEchp')": 'orange',

             "('transformer', '" + reg + "', 'powertoheat')": 'lightsalmon',

             # demand
             "('demand', '" + reg + "', 'dh')": 'red',

             # shortag / excess 
             "('bus', '" + reg + "', 'dh')_excess": 'purple',
             "('bus', '" + reg + "', 'dh')_shortage": 'blue'}
    return cdict


def stack_plot(energysystem, reg, bus, date_from, date_to):
    """
    Creates a stack plot of the specified bus.
    """
    # initialize plot
    myplot = tpd.DataFramePlot(energy_system=energysystem)
    
    # get dictionary with color of each entity in plot
    if bus == 'elec':
        cdict = color_dict(reg)
    elif bus == 'dh':
        cdict = color_dict_dh(reg)

    # slice dataframe to prepare for plot function
    myplot.slice_unstacked(
        bus_uid="('bus', '" + reg + "', '" + bus + "')",
        type="input",
        date_from=date_from,
        date_to=date_to)
    myplot.color_from_dict(cdict)

    # set plot parameters
    fig = plt.figure(figsize=(40, 14))
    plt.rc('legend', **{'fontsize': 18})
    plt.rcParams.update({'font.size': 18})
    plt.style.use('grayscale')

    # plot bus
    handles, labels = myplot.io_plot(
        bus_uid="('bus', '" + reg + "', '" + bus + "')", 
        cdict=cdict, 
        line_kwa={'linewidth': 4},
        ax=fig.add_subplot(1, 1, 1),
        date_from=date_from,
        date_to=date_to,
        )
    myplot.ax.set_ylabel('Power in MW')
    myplot.ax.set_xlabel('Date')
    myplot.ax.set_title(bus+" bus")
    myplot.set_datetime_ticks(tick_distance=24, date_format='%d-%m-%Y')
    myplot.outside_legend(handles=handles, labels=labels)

    plt.show()
    return (fig)


def sum_max_output_of_component(energysystem, from_uid, to_uid):
    """
    Returns the sum and the maximum of the flow from entity with 'from_uid'
    to entity with 'to_uid'.
    """
    results_bus = energysystem.results[[obj for obj in energysystem.entities
        if obj.uid == (from_uid)][0]]
    results_bus_component = results_bus[[obj for obj in energysystem.entities
        if obj.uid == (to_uid)][0]]
    return sum(results_bus_component), max(results_bus_component)


def timeseries_of_component(energysystem, from_uid, to_uid):
    """
    Returns the flow from entity with 'from_uid' to entity with 'to_uid'.
    """
    results_bus = energysystem.results[[obj for obj in energysystem.entities
        if obj.uid == (from_uid)][0]]
    results_bus_component = results_bus[[obj for obj in energysystem.entities
        if obj.uid == (to_uid)][0]]
    return results_bus_component


def print_validation_outputs(energysystem, reg, results_dc):
    """
    Returns sums and maximums of flows as well as full load hours of
    transformers.
    """
    # connect to database
    conn_oedb = db.connection(section='open_edb')

    # get paremeters of transformers from database
    (co2_emissions, co2_fix, eta_elec, eta_th, eta_th_chp, eta_el_chp,
     eta_chp_flex_el, sigma_chp, beta_chp, opex_var, opex_fix, capex,
     c_rate_in, c_rate_out, eta_in, eta_out,
     cap_loss, lifetime, wacc) = hlsb.get_parameters(conn_oedb)

    # list of possible power plants in region
    pp = [
        "('FixedSrc', '" + reg + "', 'wind_pwr')",
        "('FixedSrc', '" + reg + "', 'pv_pwr')",

        "('transformer', '" + reg + "', 'oil')",
        "('transformer', '" + reg + "', 'oil', 'chp')",
        "('transformer', '" + reg + "', 'oil', 'SEchp')",

        "('transformer', '" + reg + "', 'natural_gas')",
        "('transformer', '" + reg + "', 'natural_gas', 'chp')",
        "('transformer', '" + reg + "', 'natural_gas', 'SEchp')",

        "('transformer', '" + reg + "', 'natural_gas_cc')",
        "('transformer', '" + reg + "', 'natural_gas_cc', 'chp')",
        "('transformer', '" + reg + "', 'natural_gas_cc', 'SEchp')",

        "('transformer', '" + reg + "', 'biomass')",
        "('transformer', '" + reg + "', 'biomass', 'chp')",
        "('transformer', '" + reg + "', 'biomass', 'SEchp')",

        "('transformer', '" + reg + "', 'HH', 'bhkw_gas')",
        "('transformer', '" + reg + "', 'GHD', 'bhkw_gas')",

        "('transformer', '" + reg + "', 'HH', 'bhkw_bio')",
        "('transformer', '" + reg + "', 'GHD', 'bhkw_bio')",

        "('transformer', '" + reg + "', 'bhkw_bio')",
        "('transformer', '" + reg + "', 'bhkw_bio', 'dh')",

        "('transformer', '" + reg + "', 'dh_peak_heating')",

        "('transformer', '" + reg + "', 'lignite_jw', 'SEchp')",
        "('transformer', '" + reg + "', 'lignite_sp', 'SEchp')",

        "('transformer', '" + reg + "', 'powertoheat')"]

    # list of efficiencies of the above transformers
    eta_el = [
            1,
            1,
            eta_elec['oil'],
            eta_el_chp['oil'],
            eta_chp_flex_el['oil'],

            eta_elec['natural_gas'],
            eta_el_chp['natural_gas'],
            eta_chp_flex_el['natural_gas'],

            eta_elec['natural_gas_cc'],
            eta_el_chp['natural_gas_cc'],
            eta_chp_flex_el['natural_gas_cc'],

            eta_elec['biomass'],
            eta_el_chp['biomass'],
            eta_chp_flex_el['biomass'],

            eta_el_chp['bhkw_gas'],
            eta_el_chp['bhkw_gas'],

            eta_el_chp['bhkw_bio'],
            eta_el_chp['bhkw_bio'],

            eta_el_chp['bhkw_bio'],
            eta_el_chp['bhkw_bio'],

            0,  # dh_peakheating
            eta_chp_flex_el['jaenschwalde'],       
            eta_chp_flex_el['schwarzepumpe'],
            0  # powertoheat
            ]

    # list of CO2 emissions of the above transformers   
    co2 = [
            0,
            0,
            co2_emissions['oil'],
            co2_emissions['oil'],
            co2_emissions['oil'],

            co2_emissions['natural_gas'],
            co2_emissions['natural_gas'],
            co2_emissions['natural_gas'],

            co2_emissions['natural_gas_cc'],
            co2_emissions['natural_gas_cc'],
            co2_emissions['natural_gas_cc'],

            co2_emissions['biomass'],
            co2_emissions['biomass'],
            co2_emissions['biomass'],

            co2_emissions['bhkw_gas'],
            co2_emissions['bhkw_gas'],

            co2_emissions['bhkw_bio'],
            co2_emissions['bhkw_bio'],

            co2_emissions['bhkw_bio'],
            co2_emissions['bhkw_bio'],

            0,  # dh_peakheating
            co2_emissions['lignite'],       
            co2_emissions['lignite'],
            0  # powertoheat
            ]
    # get sum and maximum of each flow from transformer to bus as well as 
    # full load hours of each transformer
    ebus = "('bus', '" + reg + "', 'elec')"
    dhbus = "('bus', '" + reg + "', 'dh')" 
    summe_plant_dict = {}
    el_energy = list()
    dh_energy = list()
    for p in pp:
        print(p)
        # if flow from transformer to electricity bus
        try:
            summe_plant_dict[p], maximum = sum_max_output_of_component(
                energysystem, p, ebus)
            print(('sum:' + str(summe_plant_dict[p])))
            print(('max:' + str(maximum)))
            results_dc['sum ' + reg + str(p)] = summe_plant_dict[p]
            results_dc['max ' + reg + str(p)] = maximum
            el_energy.append(summe_plant_dict[p])
        except:
            print('nicht vorhanden')
            results_dc['sum ' + reg + str(p)] = 0
            results_dc['max ' + reg + str(p)] = 0
            el_energy.append(0)
        try:
            print(('vlh:' + str(summe_plant_dict[p] / maximum)))
            results_dc['vlh ' + reg + str(p)] = summe_plant_dict[p] / maximum
        except:
            results_dc['vlh ' + reg + str(p)] = 0
        print('\n')
        # if flow from transformer to district heating bus
        try: 
            summe_plant_dict['dh' + p], maximum = sum_max_output_of_component(
                energysystem, p, dhbus)
            print(('sum:' + str(summe_plant_dict['dh' + p])))
            print(('max:' + str(maximum)))
            results_dc['sum '+ reg + str(p) + '_dh'] = \
                summe_plant_dict['dh' + p]
            results_dc['max '+ reg + str(p) + '_dh'] = maximum
            dh_energy.append(summe_plant_dict['dh' + p])
        except:
            print('nicht vorhanden')
            dh_energy.append(0)
            results_dc['sum '+ reg + str(p)+'_dh'] = 0
            results_dc['max '+ reg + str(p)+'_dh'] = 0
        try:
            print(('vls:' + str(summe_plant_dict[p] / maximum)))
            results_dc['vlh ' + reg + str(p)+'_dh'] = (summe_plant_dict[p] /
                maximum)
        except:
            results_dc['vlh ' + reg + str(p)+'_dh'] = 0
        print('\n')

    # get sum and maximum of electricity shortage
    shortage_bus = "('bus', '" + reg + "', 'elec')_shortage"
    summe_plant, maximum = sum_max_output_of_component(
        energysystem, shortage_bus, ebus)
    print(('el_shortage_sum:' + str(summe_plant)))
    print(('el_shortage_max:' + str(maximum)))
    results_dc['el_shortage ' + reg] = str(summe_plant)
    results_dc['el_shortage_max ' + reg] = maximum
    print('\n')

    # get sum and maximum of excess in district heating
    excess_dh = "('bus', '" + reg + "', 'dh')_excess"
    summe_plant, maximum = sum_max_output_of_component(
        energysystem, dhbus, excess_dh)
    print(('dh_excess_sum:' + str(summe_plant)))
    print(('dh_excess_max:' + str(maximum)))
    results_dc['dh_excess_sum ' + reg] = summe_plant
    results_dc['dh_excess_max ' + reg] = maximum
    
    # get sum and maximum of electricity excess
    excess = "('bus', '" + reg + "', 'elec')_excess"    
    summe_plant, maximum = sum_max_output_of_component(
        energysystem, ebus, excess)
    print(('el_excess_sum:' + str(summe_plant)))
    print(('el_excess_max:' + str(maximum)))
    results_dc['el_excess_sum ' + reg] = summe_plant
    results_dc['el_excess_max ' + reg] = maximum
    
    # get sum of flows from wind turbines and pv systems to electricity bus
    sum_fee = (summe_plant_dict["('FixedSrc', '" + reg + "', 'wind_pwr')"] +
               summe_plant_dict["('FixedSrc', '" + reg + "', 'pv_pwr')"])
    print(('share excess wind + pv:' + str((summe_plant / sum_fee) * 100)))
    
    # create dataframe with power output of each transformer, electrical
    # efficiency and CO2 per MWh
    frame = pd.DataFrame(index=pp)
    frame['dh_energy'] = dh_energy    
    frame['energy_sum'] = el_energy
    frame['eta_el'] = eta_el
    frame['co2'] = co2
    
    return (results_dc, frame)


def print_exports(energysystem, results_dc, year, path):
    """
    Get exports from Brandenburg to neighbor regions and imports from neighbor 
    regions to Brandenburg.
    """
    export_from = ["('bus', 'UB', 'elec')",
                   "('bus', 'UB', 'elec')",
                   "('bus', 'PO', 'elec')",
                   "('bus', 'PO', 'elec')",
                   "('bus', 'HF', 'elec')",
                   "('bus', 'LS', 'elec')",
                   "('bus', 'HF', 'elec')",
                   "('bus', 'OS', 'elec')"]
    import_to = export_from
    export_to = ["transport_('bus', 'UB', 'elec')('bus', 'KJ', 'elec')",
                 "transport_('bus', 'UB', 'elec')('bus', 'MV', 'elec')",
                 "transport_('bus', 'PO', 'elec')('bus', 'MV', 'elec')",
                 "transport_('bus', 'PO', 'elec')('bus', 'ST', 'elec')",
                 "transport_('bus', 'HF', 'elec')('bus', 'ST', 'elec')",
                 "transport_('bus', 'LS', 'elec')('bus', 'SN', 'elec')",
                 "transport_('bus', 'BE', 'elec')('bus', 'HF', 'elec')",
                 "transport_('bus', 'BE', 'elec')('bus', 'OS', 'elec')"]
    import_from = ["transport_('bus', 'KJ', 'import')('bus', 'UB', 'elec')",
                 "transport_('bus', 'MV', 'import')('bus', 'UB', 'elec')",
                 "transport_('bus', 'MV', 'import')('bus', 'PO', 'elec')",
                 "transport_('bus', 'ST', 'import')('bus', 'PO', 'elec')",
                 "transport_('bus', 'ST', 'import')('bus', 'HF', 'elec')",
                 "transport_('bus', 'SN', 'import')('bus', 'LS', 'elec')",
                 "transport_('bus', 'HF', 'elec')('bus', 'BE', 'elec')",
                 "transport_('bus', 'OS', 'elec')('bus', 'BE', 'elec')"]

    time_index = pd.date_range('1/1/{0}'.format(year), periods=8760, freq='H')
    time_no_export = pd.DataFrame(index=time_index)
    exports = pd.DataFrame(index=time_index)
    imports = pd.DataFrame(index=time_index)
    export_total = 0
    for i in range(len(export_from)):
        print(export_to[i])
        # sum of export
        summe_ex, maximum = sum_max_output_of_component(
            energysystem, export_from[i], export_to[i])
        export_total += summe_ex
        print('export:')
        print(summe_ex)
        results_dc['export ' + export_to[i] + ' summe'] = summe_ex
        # maximum of export
        print('max:')
        print(maximum)
        results_dc['export ' + export_to[i] + ' maximum'] = maximum
        # timeseries
        exports[export_to[i]] = timeseries_of_component(
            energysystem, export_from[i], export_to[i])
        imports[export_to[i]] = timeseries_of_component(
            energysystem, import_from[i], import_to[i])
        time_no_export[export_to[i]] = (exports[export_to[i]] -
            imports[export_to[i]])
    # total export
    print('export_gesamt:')
    print(export_total)
    results_dc['export gesamt: '] = export_total
    # save import and export timeseries to csv
    exports.to_csv(path + 'exports.csv')
    imports.to_csv(path + 'imports.csv')
    time_no_export.to_csv(path + 'no_export.csv')
    
    return (results_dc, time_no_export)


def print_im_exports(energysystem, results_dc, year, path):
    """
    Adds flows between regions in Brandenburg and between Brandenburg and
    Berlin to results_dc.
    """
    export_from = ["('bus', 'UB', 'elec')",
                   "('bus', 'PO', 'elec')",
                   "('bus', 'HF', 'elec')",
                   "('bus', 'LS', 'elec')",
                   "('bus', 'OS', 'elec')",
                   "('bus', 'BE', 'elec')"]
    export_to = [
             "transport_('bus', 'UB', 'elec')('bus', 'OS', 'elec')",
             "transport_('bus', 'OS', 'elec')('bus', 'UB', 'elec')",
             "transport_('bus', 'OS', 'elec')('bus', 'LS', 'elec')",
             "transport_('bus', 'LS', 'elec')('bus', 'OS', 'elec')",
             "transport_('bus', 'OS', 'elec')('bus', 'BE', 'elec')",
             "transport_('bus', 'BE', 'elec')('bus', 'OS', 'elec')",
             "transport_('bus', 'LS', 'elec')('bus', 'HF', 'elec')",
             "transport_('bus', 'HF', 'elec')('bus', 'LS', 'elec')",
             "transport_('bus', 'HF', 'elec')('bus', 'PO', 'elec')",
             "transport_('bus', 'PO', 'elec')('bus', 'HF', 'elec')",
             "transport_('bus', 'HF', 'elec')('bus', 'BE', 'elec')",
             "transport_('bus', 'BE', 'elec')('bus', 'HF', 'elec')",
             "transport_('bus', 'PO', 'elec')('bus', 'OS', 'elec')",
             "transport_('bus', 'OS', 'elec')('bus', 'PO', 'elec')"]

    time_index = pd.date_range('1/1/{0}'.format(year), periods=8760, freq='H')
    BBB_Kuppelstellen = pd.DataFrame(index=time_index)
        
    export_all = 0
    for i in export_from:
        print(i)
        for k in export_to:
            print(k)
            try:
                summe_ex, maximum = sum_max_output_of_component(
                    energysystem, i, k)
                export_all += summe_ex
                print('from '+ i + ' to '+ k)
                print(summe_ex)
                results_dc['export from ' + i + ' to ' + k] = summe_ex
                results_dc['export from ' + i + ' to ' + k + ' maximum'] = \
                    maximum
                BBB_Kuppelstellen['export from ' + i + ' to ' + k] = \
                    timeseries_of_component(energysystem, i, k)
            except:
                pass
    # total of flows
    print('export_in_BBB_gesamt:')
    print(export_all)
    results_dc['export in BBB gesamt: '] = export_all
    # timeseries to csv
    BBB_Kuppelstellen.to_csv(path + 'kuppelstellen.csv')

    return results_dc


def get_share_ee(energysystem, reg, results_dc):
    """
    Get shares of wind and pv on demand fulfillment.
    """
    # get feedin timeseries from wind and pv to electricity bus
    ebus = "('bus', '" + reg + "', 'elec')"
    pv_time = timeseries_of_component(
        energysystem, "('FixedSrc', '" + reg + "', 'pv_pwr')", ebus)
    wind_time = timeseries_of_component(
        energysystem, "('FixedSrc', '" + reg + "', 'wind_pwr')", ebus)
    # get electricity demand timeseries
    demand_time = timeseries_of_component(
        energysystem, ebus, "('demand', '" + reg + "', 'elec')")
    
    # calculate shares
    res = pd.DataFrame(index=range(len(demand_time)),
                       columns=['ee', 'pv', 'wind'])
    for i in range(len(demand_time)):
        fee = demand_time[i] - pv_time[i] - wind_time[i]
        if fee < 0:
            res['ee'][i] = demand_time[i]
            res['pv'][i] = demand_time[i] * pv_time[i] / (
                            pv_time[i] + wind_time[i])
            res['wind'][i] = demand_time[i] * wind_time[i] / (
                            pv_time[i] + wind_time[i])
        else:
            res['ee'][i] = pv_time[i] + wind_time[i]
            res['pv'][i] = pv_time[i]
            res['wind'][i] = wind_time[i]   
    ee_share = sum(res['ee']) / sum(demand_time)
    pv_share = sum(res['pv']) / sum(demand_time)
    wind_share = sum(res['wind']) / sum(demand_time)
    
    # print shares and add to results_dc
    print('ee share:')
    print(ee_share)
    results_dc['ee share ' + reg] = ee_share
    print('pv share:')
    print(pv_share)
    results_dc['pv share ' + reg] = pv_share
    print('wind share:')
    print(wind_share)    
    results_dc['wind share ' + reg] = wind_share

    return results_dc

def co2(energysystem):
    """
    Calculate total CO2 emissions.
    """
    # retrieve specific CO2 emissions from database
    conn_oedb = db.connection(section='open_edb')
    (co2_emissions, co2_fix, eta_elec, eta_th, eta_th_chp, eta_el_chp,
     eta_chp_flex_el, sigma_chp, beta_chp, opex_var, opex_fix, capex,
     c_rate_in, c_rate_out, eta_in, eta_out,
     cap_loss, lifetime, wacc) = hlsb.get_parameters(conn_oedb)

    # fossil ressources
    global_ressources = ['natural_gas', 'natural_gas_cc', 'lignite',
        'oil', 'waste', 'hard_coal']
    # create list of global ressource buses BB
    list_global_ressource_buses = []
    for ressource in global_ressources:
        list_global_ressource_buses += ["('bus', 'BB', '" + ressource + "')"]
    # create list with entities of global ressource buses
    global_ressource_buses_bb = [obj for obj in energysystem.entities
        if any(bus in obj.uid for bus in list_global_ressource_buses)]
    # get yearly energy
    co2 = 0
    for bus in global_ressource_buses_bb:
        for output in bus.outputs:
            summe, maximum = sum_max_output_of_component(
                energysystem, bus.uid, output.uid)
            co2 += summe * co2_emissions[bus.type]

    # biogas
    biogas_transformer = [obj for obj in energysystem.entities
        if 'bhkw_bio' in obj.uid and 'transformer' in obj.uid]
    bb_regions = ['PO', 'UB', 'HF', 'OS', 'LS']
    biogas_transformer_bb = [obj for obj in biogas_transformer
        if any(region in obj.uid for region in bb_regions)]

    # write list to hand over to BB constraint
    for transformer in biogas_transformer_bb:
        summe, maximum = sum_max_output_of_component(
            energysystem, transformer.inputs[0].uid, transformer.uid)
        co2 += summe * co2_emissions[transformer.inputs[0].type]
    print('Total CO2 emissions in BB:')
    print(co2)
    return co2


def get_supply_demand_timeseries(energysystem, year, path):
    """
    Writes timeseries of all inputs and outputs of the electricity bus of
    each region as well as their sums to dataframe and saves to csv.
    """
    time_index = pd.date_range('1/1/{0}'.format(year), periods=8760, freq='H')
    # create dataframe for timeseries sums of outputs and inputs of electricity
    # bus
    supply_demand_sum = pd.DataFrame(index=time_index)
    for region in energysystem.regions:
        reg = region.name
        # create dataframe for timeseries of outputs and inputs of electricity
        # bus
        elec_out = pd.DataFrame(index=time_index) 
        elec_in = pd.DataFrame(index=time_index) 
        # get electricity bus entity and its results
        elec_bus = [obj for obj in energysystem.entities
            if obj.uid == ("('bus', '" + reg + "', 'elec')")][0]
        elec_bus_results = energysystem.results[[obj for obj in 
            energysystem.entities
            if obj.uid == ("('bus', '" + reg + "', 'elec')")][0]]
        # get outputs of electricity bus
        for obj in energysystem.entities:
            if 'demand' in obj.uid or 'hp' in obj.uid:        
                try:
                    elec_out[obj.uid] = elec_bus_results[[obj][0]]
                except:
                    pass
        # get inputs of electricity bus
        for obj in energysystem.entities:
            if ('transformer' in obj.uid or 'transport' in obj.uid or
                'FixedSrc' in obj.uid):
                obj_in = energysystem.results[[obj][0]]
                try:
                    elec_in[obj.uid] = obj_in[[elec_bus][0]]
                except:
                    pass   

        # save to csv
        elec_in.to_csv(path + reg + '_all_times_in.csv')
        elec_out.to_csv(path + reg + '_all_times_out.csv')
        # get residual as well as sum of all inputs and all outputs
        supply_demand_sum[reg] = elec_in.sum(axis=1) - elec_out.sum(axis=1)   
        supply_demand_sum[reg + 'in'] = elec_in.sum(axis=1)
        supply_demand_sum[reg + 'out'] = elec_out.sum(axis=1)   
    # save to csv 
    supply_demand_sum.to_csv(path + 'supply_minus_demand.csv')
        
    return supply_demand_sum 
    

if __name__ == "__main__":

    # load results
    path_to_dump = expanduser("~") + '/.oemof/dumps/'
    year = 2010
      # create dummy energy system
    energysystem = create_es('cbc', [t for t in range(8760)], str(year))
      # load dumped energy system
    energysystem.restore(path_to_dump)
    
    # weeks for stack plot
    date_from = {}
    date_to = {}
    date_from['spring'] = "2010-03-17 00:00:00"
    date_to['spring'] = "2010-03-24 00:00:00"
    date_from['summer'] = "2010-06-17 00:00:00"
    date_to['summer'] = "2010-06-24 00:00:00"
    date_from['autumn'] = "2010-09-17 00:00:00"
    date_to['autumn'] = "2010-09-24 00:00:00"
    date_from['winter'] = "2010-12-17 00:00:00"
    date_to['winter'] = "2010-12-24 00:00:00"
    
    # empty results_dc dictionary to write results into
    results_dc = {}
    
    # get all inputs and outputs of electricity bus of each region
    get_supply_demand_timeseries(energysystem, year, path_to_dump)
    # get exports from Brandenburg to neighbor regions and imports from  
    # neighbor regions to Brandenburg
    print_exports(energysystem, results_dc, year, path_to_dump) ##
    # add flows between regions in Brandenburg and between Brandenburg and
    # Berlin to results_dc
    print_im_exports(energysystem, results_dc, year, path_to_dump)
    # calculates total CO2 emissions
    results_dc['co2_all_BB'] = co2(energysystem)
    
    
    transformer_results_df = pd.DataFrame()
    for reg in ('HF', 'LS', 'UB', 'PO', 'BE', 'OS'):
        # create stack plots for electricity bus and district heating bus for
        # winter week
        week = 'winter' 
        for bus in ('elec', 'dh'):      
            fig = stack_plot(
                energysystem, reg, bus, date_from[week], date_to[week])
            fig.savefig(path_to_dump + reg + '_' + bus + '_' + week + '.png')
        # add sums and maximums of flows as well as full load hours of
        # transformers
        # return value frame is a dataframe with power output of each 
        # transformer, electrical efficiency and CO2 per MWh
        results_dc, frame = print_validation_outputs(
            energysystem, reg, results_dc)
        transformer_results_df = transformer_results_df.append(frame)     
        # get shares of wind and pv of demand fulfillment
        get_share_ee(energysystem, reg, results_dc)
    
    # write to csv
    transformer_results_df.to_csv(path_to_dump + 'co2_el_energy.csv')
    
    keys = list(results_dc.keys())
    values = list(results_dc.values())
    f = open(path_to_dump + '_results.csv', 'w', newline='')
    w = csv.writer(f, delimiter=';')
    w.writerow(keys)
    w.writerow(values)
    f.close