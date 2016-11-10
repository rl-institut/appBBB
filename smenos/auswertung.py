# -*- coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
from oemof.core import energy_system as es
from oemof.solph.predefined_objectives import minimize_cost
from oemof.outputlib import to_pandas as tpd


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


def color_dict():
    cdict = {
             # renewables
             'wind_onshore': 'lightblue',
             'wind_offshore': 'darkblue',
             'runofriver': 'blue',
             'pv': 'yellow',
             'geothermal': 'red',
             # storages
             'battery': 'cyan',
             'pumped_hydro': '#ffde32',
             # thermal pp
             'hard_coal_st': 'black',
             'hard_coal_chp': 'black',
             'waste_chp': '#42c77a',
             'lignite_st': 'brown',
             'lignite_chp': 'brown',
             'gas_gt': 'lightgrey',
             'gas_cc': 'grey',
             'gas_chp': 'darkgrey',
             'biomass_st': 'darkgreen',
             'uranium_st': 'pink',
             # else
             'el_demand': '#ce4aff',
             'electricity_shortage': 'purple',
             'electricity_excess': 'black',
             # h2
             'ely': 'darkgoldenrod',
             'methanisation': 'blue',
             'h2_storage': 'orange',
             'ch4_storage': 'orange'}
    return cdict


def stack_plot_gas(energysystem):
    # Plotting a combined stacked plot
    myplot = tpd.DataFramePlot(energy_system=energysystem)

    cdict = color_dict()

    ## Plotting the input flows of the electricity bus
    myplot.slice_unstacked(bus_uid="gas", type="input",
                           date_from="2050-01-01 00:00:00",
                           date_to="2050-01-14 00:00:00")
    myplot.color_from_dict(cdict)

    fig = plt.figure(figsize=(24, 14))
    plt.rc('legend', **{'fontsize': 19})
    plt.rcParams.update({'font.size': 19})
    plt.style.use('grayscale')

    handles, labels = myplot.io_plot(
        bus_uid="gas", cdict=cdict, line_kwa={'linewidth': 4},
        ax=fig.add_subplot(1, 1, 1),
        date_from="2050-01-01 00:00:00",
        date_to="2050-01-8 00:00:00",
        )
    myplot.ax.set_ylabel('Power in MW')
    myplot.ax.set_xlabel('Date')
    myplot.ax.set_title("Electricity bus")
    myplot.set_datetime_ticks(tick_distance=24, date_format='%d-%m-%Y')
    myplot.outside_legend(handles=handles, labels=labels)

    plt.show()
    return


def stack_plot(energysystem):
    # Plotting a combined stacked plot
    myplot = tpd.DataFramePlot(energy_system=energysystem)

    cdict = color_dict()

    ## Plotting the input flows of the electricity bus
    myplot.slice_unstacked(bus_uid="electricity", type="input",
                           date_from="2050-01-01 00:00:00",
                           date_to="2050-01-14 00:00:00")
    myplot.color_from_dict(cdict)

    fig = plt.figure(figsize=(24, 14))
    plt.rc('legend', **{'fontsize': 19})
    plt.rcParams.update({'font.size': 19})
    plt.style.use('grayscale')

    handles, labels = myplot.io_plot(
        bus_uid="electricity", cdict=cdict, line_kwa={'linewidth': 4},
        ax=fig.add_subplot(1, 1, 1),
        date_from="2050-01-01 00:00:00",
        date_to="2050-01-8 00:00:00",
        )
    myplot.ax.set_ylabel('Power in MW')
    myplot.ax.set_xlabel('Date')
    myplot.ax.set_title("Electricity bus")
    myplot.set_datetime_ticks(tick_distance=24, date_format='%d-%m-%Y')
    myplot.outside_legend(handles=handles, labels=labels)

    plt.show()
    return


def sum_max_output_of_component(energysystem, from_uid, to_uid):
    results_bus = energysystem.results[[obj for obj in energysystem.entities
        if obj.uid == (from_uid)][0]]
    results_bus_component = results_bus[[obj for obj in energysystem.entities
        if obj.uid == (to_uid)][0]]
    return sum(results_bus_component), max(results_bus_component)


def res_share(energysystem):
    # conventional
    from_uid = ['commodity_hard_coal', 'commodity_gas', 'commodity_waste',
        'commodity_uranium', 'commodity_lignite']
    to_uid = ['hard_coal', 'gas', 'waste', 'uranium', 'lignite']

    summe_conv = 0
    for i in range(len(from_uid)):
        summe_plant, maximum = sum_max_output_of_component(
            energysystem, from_uid[i], to_uid[i])
        summe_conv += summe_plant

    # renewables
    from_uid = ['biomass_st', 'wind_onshore', 'wind_offshore',
        'pv', 'geothermal', 'runofriver']
    summe_res = 0
    for i in range(len(from_uid)):
        summe_plant, maximum = sum_max_output_of_component(
            energysystem, from_uid[i], 'electricity')
        summe_res += summe_plant

    # shortage
    summe_shortage, maximum = sum_max_output_of_component(
            energysystem, 'electricity_shortage', 'electricity')

    return ((summe_res - summe_shortage) /
        (summe_conv + summe_res - summe_shortage))


def print_validation_outputs(energysystem):

    # capacities of pp
    pp = ['wind_onshore', 'wind_offshore', 'runofriver', 'pv',
            'geothermal', 'hard_coal_st', 'gas_gt', 'gas_chp', 'gas_cc',
            'waste_chp', 'uranium_st', 'lignite_st', 'lignite_chp',
            'hard_coal_chp', 'biomass_st']
    summe_plant_dict = {}
    for p in pp:
        print(p)
        summe_plant_dict[p], maximum = sum_max_output_of_component(
            energysystem, p, 'electricity')
        print(('sum:' + str(summe_plant_dict[p])))
        print(('max:' + str(maximum)))
        try:
            print(('vls:' + str(summe_plant_dict[p] / maximum)))
        except:
            pass
        print('\n')
    # ely
    summe_plant, maximum = sum_max_output_of_component(
            energysystem, 'electricity', 'ely')
    print('ely')
    print(('sum:' + str(summe_plant)))
    print(('max:' + str(maximum)))
    try:
        print(('vls:' + str(summe_plant / maximum)))
    except:
        pass
    print('\n')
    # methanisation
    summe_plant, maximum = sum_max_output_of_component(
            energysystem, 'methanisation', 'gas')
    print('methanisation')
    print(('sum:' + str(summe_plant)))
    print(('max:' + str(maximum)))
    try:
        print(('vls:' + str(summe_plant / maximum)))
    except:
        pass
    print('\n')

    # storages
    # battery
    summe_plant, maximum = sum_max_output_of_component(
            energysystem, 'battery', 'electricity')
    print('battery')
    print(('sum:' + str(summe_plant)))
    print(('max:' + str(maximum)))
    print('\n')
    # pumped hydro
    summe_plant, maximum = sum_max_output_of_component(
            energysystem, 'pumped_hydro', 'electricity')
    print('pumped_hydro')
    print(('sum:' + str(summe_plant)))
    print(('max:' + str(maximum)))
    print('\n')
    # h2_storage
    summe_plant, maximum = sum_max_output_of_component(
            energysystem, 'h2_storage', 'h2')
    print('h2_storage')
    print(('sum:' + str(summe_plant)))
    print(('max:' + str(maximum)))
    print('\n')

    # shortage
    summe_plant, maximum = sum_max_output_of_component(
            energysystem, 'electricity_shortage', 'electricity')
    print(('el_shortage_sum:' + str(summe_plant)))
    print(('el_shortage_max:' + str(maximum)))
    print('\n')

    # excess
    summe_plant, maximum = sum_max_output_of_component(
            energysystem, 'electricity', 'electricity_excess')
    print(('el_excess_sum:' + str(summe_plant)))
    print(('el_excess_max:' + str(maximum)))
    sum_fee = (summe_plant_dict['wind_onshore'] +
        summe_plant_dict['wind_offshore'] +
        summe_plant_dict['pv'])
    print(('share excess:' + str((summe_plant / sum_fee) * 100)))
    return


# load dumped energy system
year = 2050
energysystem = create_es(
    'cbc', [t for t in range(8760)], str(year))
energysystem.restore()

results_bus = energysystem.results[[obj for obj in energysystem.entities
        if obj.uid == ("('bus', 'BE', 'elec')_shortage")][0]]
print(results_bus)
## anteil ee
#print(res_share(energysystem))

## capacities
#print_validation_outputs(energysystem)

# stack el
# stack h2
# stack gas
#stack_plot(energysystem)

## print all buses
#for entity in energysystem.entities:
    #try:
        #print(entity.uid)
        #print(entity.crf)
    #except:
        #pass