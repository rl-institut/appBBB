#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Created on Thu Jun 16 10:35:19 2016

changed from feedin_pg in order to create fixedSource Offshore
"""


import pandas as pd
import os.path as path

from oemof.db import coastdat
from feedinlib import powerplants as pp
import helper_SmEnOs as hls
from shapely.wkt import loads as wkt_loads


class Feedin:
    ''

    def __init__(self):
        ''
        pass

    def aggregate_cap_val(self, conn, **kwargs):
        '''
        Returns the normalised feedin profile and installed capacity for
        a given region.

        Parameters
        ----------
        region : Region instance
        region.geom : shapely.geometry object
            Geo-spatial data with information for location/region-shape. The
            geometry can be a polygon/multi-polygon or a point.

        Returns
        -------
        feedin_df : pandas dataframe
            Dataframe containing the normalised feedin profile of the given
            windpower plants. Index of the dataframe is the hour of the year; 
        cap : pandas series
            Series containing the summed installed capacity of given wind
            turbines.
        '''
        
        [wind_df, cap] = self.get_timeseries(
            conn,
            **kwargs)


        if kwargs.get('store', False):
            self.store_full_df(wind_df, **kwargs)

        # Summerize the results to one column
        cap = cap.sum()
        df = (wind_df.sum(axis=1) / cap)
        feedin_df = df

        return feedin_df, cap

    def get_timeseries(self, conn, **kwargs):
        # get capacities and geometries from db 
        os_pps = hls.get_offshore_pps(conn, kwargs['schema'], 
                                      kwargs['table'], kwargs['start_year'])      

        wind_df = 0
        wind_cap = {}

        series_name_0 = 'Park'
        counter = 0
        for opp in os_pps.iterrows():    
            counter = counter + 1
            series_name = series_name_0 + str(counter)
            geometry = wkt_loads(opp[1][1])
            
            # get weather for each wind park
            w_cell = coastdat.get_weather(conn, geometry, 
                                          kwargs['year'])

            kwargs['wind_conv_type'] = kwargs['wka_model']
            kwargs['d_rotor'] = kwargs['d_rotor']
            kwargs['h_hub'] = kwargs['h_hub']
            
#                # Determine the feedin time series for the weather cell
#                # Wind energy
            wind_peak_power = opp[1][0]
            wind_power_plant = pp.WindPowerPlant(**kwargs)
            wind_series = wind_power_plant.feedin(
                    weather=w_cell, installed_capacity=wind_peak_power)
            wind_series.name = series_name
            wind_cap[series_name] = wind_peak_power
#
#                # Combine the results to a DataFrame
            try:
                wind_df = pd.concat([wind_df, wind_series], axis=1)
            except:
                wind_df = wind_series.to_frame()

        # Write capacity into a dataframe
        capw = pd.Series(pd.DataFrame.from_dict(wind_cap, orient='index')[0])
        capw.name = 'wind_pwr'
        print(wind_peak_power)

        return wind_df, capw

    def store_full_df(self, pv_df, wind_df, **kwargs):
        ''
        dpath = kwargs.get(
            'dpath', path.join(path.expanduser("~"), '.oemof'))
        filename = kwargs.get('filename', 'feedin_offshore')
        fullpath = path.join(dpath, filename)

        if kwargs['store'] == 'hf5':
            pv_df.to_hdf(fullpath + '.hf5', 'pv_pwr')
            wind_df.to_hdf(fullpath + '.hf5', 'wind_pwr')

        if kwargs['store'] == 'csv':
            pv_df.to_csv(fullpath + '_pv.csv')
            wind_df.to_csv(fullpath + '_wind.csv')

