# -*- coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
from oemof.core import energy_system as es
from oemof.solph.predefined_objectives import minimize_cost
from oemof.outputlib import to_pandas as tpd
import csv
from oemof import db
import helper_BBB as hlsb


def create_es(solver, timesteps, year):
    simulation = es.Simulation(solver=solver, timesteps=timesteps,
                               stream_solver_output=True,
                               debug=False, duals=True,
                               objective_options={"function": minimize_cost})

    # Adding a time index to the energy system
    time_index = pd.date_range('1/1/' + year,
                               periods=len(timesteps),
                               freq='H')
    energysystem = es.EnergySystem(time_idx=time_index,
                                   simulation=simulation)
    return energysystem


def color_dict(reg):
    cdict = {
             # renewables
             "('FixedSrc', '"+reg+"', 'wind_pwr')": 'lightblue',
             "('FixedSrc', '"+reg+"', 'pv_pwr')": 'yellow',

             "('transformer', '"+reg+"', 'oil')": 'black',
             "('transformer', '"+reg+"', 'oil', 'chp')": 'black',
             "('transformer', '"+reg+"', 'oil', 'SEchp')": 'black',

             "('transformer', '"+reg+"', 'natural_gas')": 'lightgrey',
             "('transformer', '"+reg+"', 'natural_gas', 'chp')": 'lightgrey',
             "('transformer', '"+reg+"', 'natural_gas', 'SEchp')": 'lightgrey',

             "('transformer', '"+reg+"', 'natural_gas_cc')": 'darkgrey',
             "('transformer', '"+reg+"', 'natural_gas_cc', 'chp')": 'darkgrey',
             "('transformer', '"+reg+"', 'natural_gas_cc', 'SEchp')": 'darkgrey',

             "('transformer', '"+reg+"', 'HH', 'bhkw_gas')": 'grey',
             "('transformer', '"+reg+"', 'GHD', 'bhkw_gas')": 'grey',

             "('transformer', '"+reg+"', 'biomass')": 'lightgreen',
             "('transformer', '"+reg+"', 'biomass', 'chp')": 'lightgreen',
             "('transformer', '"+reg+"', 'biomass', 'SEchp')": 'lightgreen',

             "('transformer', '"+reg+"', 'HH', 'bhkw_bio')": 'green',
             "('transformer', '"+reg+"', 'GHD', 'bhkw_bio')": 'green',

             "('transformer', '"+reg+"', 'powertoheat')": 'lightsalmon',

             "('transformer', '"+reg+"', 'lignite_jw', 'SEchp')": 'brown',
             "('transformer', '"+reg+"', 'lignite_sp', 'SEchp')": 'orange',
             "('demand', '"+reg+"', 'elec')": 'red',
             "('demand', '"+reg+"', 'elec', 'mob')": 'red',
             "('bus', '"+reg+"', 'elec')_excess": 'purple',
             "('bus', '"+reg+"', 'elec')_shortage": 'blueviolet',
             "('transformer', '"+reg+"', 'hp', 'brine', 'ww')": 'blue',
             "('transformer', '"+reg+"', 'hp', 'brine', 'heating')": 'blue',
             "('transformer', '"+reg+"', 'hp', 'air', 'ww')": 'blue',
             "('transformer', '"+reg+"', 'hp', 'air', 'heating')": 'blue',
             "('transformer', '"+reg+"', 'hp', 'air', 'ww', 'rod')": 'blue',
             "('transformer', '"+reg+"', 'hp', 'air', 'heating', 'rod')": 'blue',

             "transport_('bus', 'UB', 'elec')('bus', 'OS', 'elec')": 'salmon',
             "transport_('bus', 'OS', 'elec')('bus', 'UB', 'elec')": 'salmon',
             "transport_('bus', 'OS', 'elec')('bus', 'LS', 'elec')": 'chocolate',
             "transport_('bus', 'LS', 'elec')('bus', 'OS', 'elec')": 'chocolate',
             "transport_('bus', 'OS', 'elec')('bus', 'BE', 'elec')": 'peru',
             "transport_('bus', 'BE', 'elec')('bus', 'OS', 'elec')": 'peru',
             "transport_('bus', 'LS', 'elec')('bus', 'HF', 'elec')": 'burlywood',
             "transport_('bus', 'HF', 'elec')('bus', 'LS', 'elec')": 'burlywood',
             "transport_('bus', 'HF', 'elec')('bus', 'PO', 'elec')": 'goldenrod',
             "transport_('bus', 'PO', 'elec')('bus', 'HF', 'elec')": 'goldenrod',
             "transport_('bus', 'HF', 'elec')('bus', 'BE', 'elec')": 'khaki',
             "transport_('bus', 'BE', 'elec')('bus', 'HF', 'elec')": 'khaki',
             "transport_('bus', 'PO', 'elec')('bus', 'OS', 'elec')": 'indianred',
             "transport_('bus', 'OS', 'elec')('bus', 'PO', 'elec')": 'indianred',

             "transport_('bus', 'UB', 'elec')('bus', 'KJ', 'elec')": 'lime',
             "transport_('bus', 'UB', 'elec')('bus', 'MV', 'elec')": 'cyan',
             "transport_('bus', 'PO', 'elec')('bus', 'MV', 'elec')": 'teal',
             "transport_('bus', 'PO', 'elec')('bus', 'ST', 'elec')": 'seagreen',
             "transport_('bus', 'HF', 'elec')('bus', 'ST', 'elec')": 'yellowgreen',
             "transport_('bus', 'LS', 'elec')('bus', 'SN', 'elec')": 'turquoise',
             "transport_('bus', 'BE', 'elec')('bus', 'HF', 'elec')": 'olive',
             "transport_('bus', 'BE', 'elec')('bus', 'OS', 'elec')": 'lightseagreen',

             "transport_('bus', 'KJ', 'import')('bus', 'UB', 'elec')": 'lime',
             "transport_('bus', 'MV', 'import')('bus', 'UB', 'elec')": 'cyan',
             "transport_('bus', 'MV', 'import')('bus', 'PO', 'elec')": 'teal',
             "transport_('bus', 'ST', 'import')('bus', 'PO', 'elec')": 'seagreen',
             "transport_('bus', 'ST', 'import')('bus', 'HF', 'elec')": 'yellowgreen',
             "transport_('bus', 'SN', 'import')('bus', 'LS', 'elec')": 'turquoise',
             "transport_('bus', 'HF', 'elec')('bus', 'BE', 'elec')": 'olive',
             "transport_('bus', 'OS', 'elec')('bus', 'BE', 'elec')": 'lightseagreen'}
    return cdict

def color_dict_dh(reg):
    cdict = {
             "('transformer', '"+reg+"', 'oil', 'chp')": 'black',
             "('transformer', '"+reg+"', 'oil', 'SEchp')": 'black',
             "('heat_transformer', '"+reg+"', 'oil')": 'black',

             "('transformer', '"+reg+"', 'natural_gas', 'chp')": 'lightgrey',
             "('transformer', '"+reg+"', 'natural_gas', 'SEchp')": 'lightgrey',
             "('heat_transformer', '"+reg+"', 'natural_gas')": 'lightgrey',
             "('transformer', '"+reg+"', 'dh_peak_heating')": 'khaki',

             "('transformer', '"+reg+"', 'natural_gas_cc', 'chp')": 'darkgrey',
             "('transformer', '"+reg+"', 'natural_gas_cc', 'SEchp')": 'darkgrey',

             "('transformer', '"+reg+"', 'biomass', 'chp')": 'lightgreen',
             "('transformer', '"+reg+"', 'biomass', 'SEchp')": 'lightgreen',
             "('heat_transformer', '"+reg+"', 'biomass')": 'lightgreen',

             "('transformer', '"+reg+"', 'lignite_jw', 'SEchp')": 'brown',
             "('transformer', '"+reg+"', 'lignite_sp', 'SEchp')": 'orange',

             "('transformer', '"+reg+"', 'powertoheat')": 'lightsalmon',

             "('demand', '"+reg+"', 'dh')": 'red',

             "('bus', '"+reg+"', 'dh')_excess": 'purple',
             "('bus', '"+reg+"', 'dh')_shortage": 'blue'}
    return cdict


def stack_plot(energysystem, reg, bus, date_from, date_to):
    # Plotting a combined stacked plot
    myplot = tpd.DataFramePlot(energy_system=energysystem)

    if bus == 'elec':
        cdict = color_dict(reg)
    elif bus == 'dh':
        cdict = color_dict_dh(reg)

    ## Plotting the input flows of the electricity bus
    myplot.slice_unstacked(bus_uid="('bus', '"+reg+"', '"+bus+"')", type="input",
                           date_from=date_from,
                           date_to=date_to)
    myplot.color_from_dict(cdict)

    fig = plt.figure(figsize=(40, 14))
    plt.rc('legend', **{'fontsize': 18})
    plt.rcParams.update({'font.size': 18})
    plt.style.use('grayscale')

    handles, labels = myplot.io_plot(
        bus_uid="('bus', '"+reg+"', '"+bus+"')", cdict=cdict, line_kwa={'linewidth': 4},
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
    results_bus = energysystem.results[[obj for obj in energysystem.entities
        if obj.uid == (from_uid)][0]]
    results_bus_component = results_bus[[obj for obj in energysystem.entities
        if obj.uid == (to_uid)][0]]
    return sum(results_bus_component), max(results_bus_component)

def timeseries_of_component(energysystem, from_uid, to_uid):
    results_bus = energysystem.results[[obj for obj in energysystem.entities
        if obj.uid == (from_uid)][0]]
    results_bus_component = results_bus[[obj for obj in energysystem.entities
        if obj.uid == (to_uid)][0]]
    return results_bus_component

def timeseries_of_component_dh(energysystem, from_uid, to_uid):
    results_bus = energysystem.results[[obj for obj in energysystem.entities
        if obj.uid == (from_uid)][0]]
    results_bus_component = results_bus[[obj for obj in energysystem.entities
        if obj.uid == (to_uid)][0]]
    return results_bus_component


def print_validation_outputs(energysystem, reg, results_dc):
    conn_oedb = db.connection(section='open_edb')

    (co2_emissions, co2_fix, eta_elec, eta_th, eta_th_chp, eta_el_chp,
     eta_chp_flex_el, sigma_chp, beta_chp, opex_var, opex_fix, capex,
     c_rate_in, c_rate_out, eta_in, eta_out,
     cap_loss, lifetime, wacc) = hlsb.get_parameters(conn_oedb)

    # capacities of pp
    pp = [
        "('FixedSrc', '"+reg+"', 'wind_pwr')",
        "('FixedSrc', '"+reg+"', 'pv_pwr')",

        "('transformer', '"+reg+"', 'oil')",
        "('transformer', '"+reg+"', 'oil', 'chp')",
        "('transformer', '"+reg+"', 'oil', 'SEchp')",

        "('transformer', '"+reg+"', 'natural_gas')",
        "('transformer', '"+reg+"', 'natural_gas', 'chp')",
        "('transformer', '"+reg+"', 'natural_gas', 'SEchp')",

        "('transformer', '"+reg+"', 'natural_gas_cc')",
        "('transformer', '"+reg+"', 'natural_gas_cc', 'chp')",
        "('transformer', '"+reg+"', 'natural_gas_cc', 'SEchp')",

        "('transformer', '"+reg+"', 'biomass')",
        "('transformer', '"+reg+"', 'biomass', 'chp')",
        "('transformer', '"+reg+"', 'biomass', 'SEchp')",

        "('transformer', '"+reg+"', 'HH', 'bhkw_gas')",
        "('transformer', '"+reg+"', 'GHD', 'bhkw_gas')",

        "('transformer', '"+reg+"', 'HH', 'bhkw_bio')",
        "('transformer', '"+reg+"', 'GHD', 'bhkw_bio')",

        "('transformer', '"+reg+"', 'bhkw_bio')",
        "('transformer', '"+reg+"', 'bhkw_bio', 'dh')",

        "('transformer', '"+reg+"', 'dh_peak_heating')",

        "('transformer', '"+reg+"', 'lignite_jw', 'SEchp')",
        "('transformer', '"+reg+"', 'lignite_sp', 'SEchp')",

        "('transformer', '"+reg+"', 'powertoheat')"]

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

    ebus = "('bus', '"+reg+"', 'elec')"
    dhbus = "('bus', '"+reg+"', 'dh')"
    short = "('bus', '"+reg+"', 'elec')_shortage"
    short_dh = "('bus', '"+reg+"', 'dh')_shortage"
    excess = "('bus', '"+reg+"', 'elec')_excess"
    excess_dh = "('bus', '"+reg+"', 'dh')_excess"
    summe_plant_dict = {}
    frame = pd.DataFrame(index=pp)
    el_energy = list()
    dh_energy = list()
    
    for p in pp:
        print(p)
        try:  # el_transformer
            summe_plant_dict[p], maximum = sum_max_output_of_component(
                energysystem, p, ebus)
            print(('sum:' + str(summe_plant_dict[p])))
            print(('max:' + str(maximum)))
            results_dc['sum '+ reg + str(p)] = summe_plant_dict[p]
            results_dc['max '+ reg + str(p)] = maximum
            el_energy.append(summe_plant_dict[p])
        except:
            print('nicht vorhanden')
            el_energy.append(0)
            results_dc['sum '+ reg + str(p)] = 0
            results_dc['max '+ reg + str(p)] = 0
        try:
            print(('vls:' + str(summe_plant_dict[p] / maximum)))
            results_dc['vlh ' + reg + str(p)] = summe_plant_dict[p] / maximum
        except:
            results_dc['vlh ' + reg + str(p)] = 0
        print('\n')
        try:  # heat transformer
            summe_plant_dict['dh'+p], maximum = sum_max_output_of_component(
                energysystem, p, dhbus)
            print(('sum:' + str(summe_plant_dict['dh'+p])))
            print(('max:' + str(maximum)))
            results_dc['sum '+ reg + str(p)+'_dh'] = summe_plant_dict['dh'+p]
            results_dc['max '+ reg + str(p)+'_dh'] = maximum
            dh_energy.append(summe_plant_dict['dh'+p])
        except:
            print('nicht vorhanden')
            dh_energy.append(0)
            results_dc['sum '+ reg + str(p)+'_dh'] = 0
            results_dc['max '+ reg + str(p)+'_dh'] = 0
        try:
            print(('vls:' + str(summe_plant_dict[p] / maximum)))
            results_dc['vlh ' + reg + str(p)+'_dh'] = summe_plant_dict[p] / maximum
        except:
            results_dc['vlh ' + reg + str(p)+'_dh'] = 0
        print('\n')

    # shortage
    summe_plant, maximum = sum_max_output_of_component(
        energysystem, short, ebus)
    print(('el_shortage_sum:' + str(summe_plant)))
    results_dc['el_shortage ' + reg] = str(summe_plant)
    print(('el_shortage_max:' + str(maximum)))
    results_dc['el_shortage_max ' + reg] = maximum
    print('\n')
    
    frame['dh_energy'] = dh_energy    
    frame['energy_sum'] = el_energy
    frame['eta_el'] = eta_el
#    frame['PE'] = float(el_energy) / float(eta_el)
    frame['co2'] = co2
#    frame['emissions'] = float(el_energy) * float(co2) / float(eta_el)

#    summe_plant, maximum = sum_max_output_of_component(
#        energysystem, short_dh, dh)
#    print(('dh_shortage_sum:' + str(summe_plant)))
#    results_dc[reg + ' dh_shortage'] = summe_plant
#    print(('dh_shortage_max:' + str(maximum)))
#    results_dc[reg + ' dh_shortage_max'] = maximum
#    print('\n')

    # excess
    summe_plant, maximum = sum_max_output_of_component(
        energysystem, ebus, excess)
    print(('el_excess_sum:' + str(summe_plant)))
    results_dc['el_excess_sum ' + reg] = summe_plant
    print(('el_excess_max:' + str(maximum)))
    results_dc['el_excess_max ' + reg] = maximum
# excess_dh
    summe_plant, maximum = sum_max_output_of_component(
        energysystem, dhbus, excess_dh)
    print(('dh_excess_sum:' + str(summe_plant)))
    results_dc['dh_excess_sum ' + reg] = summe_plant
    print(('dh_excess_max:' + str(maximum)))
    results_dc['dh_excess_max ' + reg] = maximum

    sum_fee = (summe_plant_dict["('FixedSrc', '"+reg+"', 'wind_pwr')"] +
               summe_plant_dict["('FixedSrc', '"+reg+"', 'pv_pwr')"])
    print(('share excess:' + str((summe_plant / sum_fee) * 100)))
    return (results_dc, frame)


def print_exports(energysystem, results_dc):

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

    time_index = pd.date_range('1/1/{0}'.format(2010), periods=8760, freq='H')
    time_no_export = pd.DataFrame(index=time_index)
    exports = pd.DataFrame(index=time_index)
    imports = pd.DataFrame(index=time_index)
    export_all = 0
    for i in range(len(export_from)):
        print(export_to[i])
#        time = timeseries_of_component(
#                energysystem, export_from[i], export_to[i])
#        print(time)
        summe_ex, maximum = sum_max_output_of_component(
            energysystem, export_from[i], export_to[i])
        export_all += summe_ex
        print('export:')
        print(summe_ex)
        results_dc['export ' + export_to[i] + ' summe'] = summe_ex

        print('max:')
        print(maximum)
        results_dc['export ' + export_to[i] + ' maximum'] = maximum
        exports[export_to[i]] = timeseries_of_component(
            energysystem, export_from[i], export_to[i])
        imports[export_to[i]] = timeseries_of_component(
            energysystem, import_from[i], import_to[i])
        time_no_export[export_to[i]] = exports[export_to[i]] - imports[export_to[i]]
    print('export_gesamt:')
    print(export_all)
    exports.to_csv(path+'exports.csv')
    imports.to_csv(path+'imports.csv')
    time_no_export.to_csv(path+'no_export.csv')
    results_dc['export gesamt: '] = export_all
    

    return(results_dc, time_no_export)


def print_im_exports(energysystem, results_dc):

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

    time_index = pd.date_range('1/1/{0}'.format(2010), periods=8760, freq='H')
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
                print('from '+i+' to '+k)
                print(summe_ex)
                results_dc['export from ' + i + ' to ' + k] = summe_ex
                results_dc['export from ' + i + ' to ' + k + ' maximum'] = maximum
                BBB_Kuppelstellen['export from ' + i + ' to ' + k] = \
                    timeseries_of_component(energysystem, i, k)
            except:
                pass

#        time = timeseries_of_component(
#                energysystem, export_from[i], export_to[i])
#        print(time)
    print('export_in_BBB_gesamt:')
    print(export_all)
    results_dc['export in BBB gesamt: '] = export_all
    BBB_Kuppelstellen.to_csv(path+'kuppelstellen.csv')

    return(results_dc)


## capacities
def get_share_ee(energysystem, reg, results_dc):

    dh_transformer = (
            "('transformer', '"+reg+"', 'oil', 'chp')",
            "('transformer', '"+reg+"', 'oil', 'SEchp')",
            "('transformer', '"+reg+"', 'biomass', 'chp')",
            "('transformer', '"+reg+"', 'biomass', 'SEchp')",
            "('transformer', '"+reg+"', 'natural_gas', 'chp')",
            "('transformer', '"+reg+"', 'natural_gas', 'SEchp')",
            "('transformer', '"+reg+"', 'natural_gas_cc', 'chp')",
            "('transformer', '"+reg+"', 'natural_gas_cc', 'SEchp')",
            "('transformer', '"+reg+"', 'lignite_sp', 'SEchp')",
            "('transformer', '"+reg+"', 'lignite_jw', 'SEchp')")

    ebus = "('bus', '"+reg+"', 'elec')"
    pv_time = timeseries_of_component(
                energysystem, "('FixedSrc', '"+reg+"', 'pv_pwr')", ebus)
    wind_time = timeseries_of_component(
                energysystem, "('FixedSrc', '"+reg+"', 'wind_pwr')", ebus)
    demand_time = timeseries_of_component(
                energysystem, ebus, "('demand', '"+reg+"', 'elec')")
                
    res = pd.DataFrame(index=range(len(demand_time)), columns=['ee', 'pv', 'wind'])
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
    print('ee share:')
    print(ee_share)
    results_dc['ee share ' + reg] = ee_share
    print('pv share:')
    print(pv_share)
    results_dc['pv share ' + reg] = pv_share
    print('wind share:')
    print(wind_share)    
    results_dc['wind share ' + reg] = wind_share

    return(results_dc)

def co2(energysystem):
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


def get_supply_demand_timeseries(energysystem):
    time_index = pd.date_range('1/1/{0}'.format(2010), periods=8760, freq='H')
    supply_demand_time = pd.DataFrame(index=time_index)
    for regio in energysystem.regions:
        reg = regio.name
        all_times_out = pd.DataFrame(index=time_index)
        all_times_in = pd.DataFrame(index=time_index)
    
        elec_bus = energysystem.results[[obj for obj in energysystem.entities
            if obj.uid == ("('bus', '"+reg+"', 'elec')")][0]]
        e_bus = [obj for obj in energysystem.entities
            if obj.uid == ("('bus', '"+reg+"', 'elec')")][0]
        for obj in energysystem.entities:
            if 'demand' in obj.uid or 'hp' in obj.uid:        
                try:
                    all_times_out[obj.uid] = elec_bus[[obj][0]]
                except:
                    pass
        for obj in energysystem.entities:
            if 'transformer' in obj.uid or 'transport' in obj.uid or 'FixedSrc' in obj.uid:
                obj_in = energysystem.results[[obj][0]]
                try:
                    all_times_in[obj.uid] = obj_in[[e_bus][0]]
                except:
                    pass
        time_in_sum = all_times_in.sum(axis=1)
        time_out_sum = all_times_out.sum(axis=1)    

        all_times_in.to_csv(path+reg+'_all_times_in.csv')
        all_times_out.to_csv(path+reg+'_all_times_out.csv')
    
        supply_demand_time[reg] = time_in_sum - time_out_sum
        supply_demand_time[reg+'in'] = time_in_sum
        supply_demand_time[reg+'out'] = time_out_sum

    return(supply_demand_time)
    

################# get results ############################

path = '/home/hendrik/UserShares/Elisa.Gaudchau/Oemof/dumps/Szenario_gruene2030_ohne_Braunkohle/'
# load dumped energy system
year = 2050
energysystem = create_es(
    'cbc', [t for t in range(8760)], str(year))
energysystem.restore(path)

buses = ('elec', 'dh')
regions_BBB = ('HF', 'LS', 'UB', 'PO', 'BE', 'OS')

#for week in ('spring', 'summer', 'autumn', 'winter'):

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

#results_dc = {}
#results_dc['co2_all_BB'] = co2(energysystem)
#
#
#supply_demand_time = get_supply_demand_timeseries(energysystem)
#supply_demand_time.to_csv(path+'supply_minus_demand.csv')
#
#print_exports(energysystem, results_dc)
#
#print_im_exports(energysystem, results_dc)
#frame_base = pd.DataFrame()
for reg in regions_BBB:
    week = 'winter' 
    for bus in buses:      
        fig = stack_plot(energysystem, reg, bus, date_from[week], date_to[week])
        fig.savefig(path+reg+'_'+bus+'_'+week+'.png')

#    results_dc, frame = print_validation_outputs(energysystem, reg, results_dc)
#    frame_base = frame_base.append(frame)       
#    get_share_ee(energysystem, reg, results_dc)
#
#frame_base.to_csv(path+'co2_el_energy.csv')
#
#x = list(results_dc.keys())
#y = list(results_dc.values())
#f = open(path + '_results.csv', 'w', newline='')
#w = csv.writer(f, delimiter=';')
#w.writerow(x)
#w.writerow(y)
#f.close
#
#f = open(path + '_results.csv', 'w', newline='')
#w = csv.writer(f, delimiter=';')
#w.writerow(x)
#w.writerow(y)
#f.close
