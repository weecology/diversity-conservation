from __future__ import division
from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import colorsys
from numpy import log10
import pandas as pd
import numpy
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
import macroecotools
import seaborn as sns
import random
import glob
import os

included_species = pd.read_csv('data/mapping_data/included_species_ids.csv')
#clean up excluded families
included_species = included_species[(included_species['AOU'] > 2880)]
included_species = included_species[(included_species['AOU'] < 3650) | (included_species['AOU'] > 3810)]
included_species = included_species[(included_species['AOU'] < 3900) | (included_species['AOU'] > 3910)]
included_species = included_species[(included_species['AOU'] < 4160) | (included_species['AOU'] > 4210)]
included_species = included_species[(included_species['AOU'] != 7010)]
included_species.rename(columns= {'AOU':'species'}, inplace = True)
AOU_list = pd.DataFrame(included_species['species'])
sisid_list = pd.DataFrame(included_species['sisid'])

data = pd.read_csv('data/mapping_data/bbs_species_2016.csv', delimiter=',', sep='\s*,\s*')
data.rename(columns = {'site_id':'site'}, inplace = True)
data.rename(columns = {'species_id':'species'}, inplace = True)
data = pd.merge(data, AOU_list, how='inner', on=['species']) #exclude species whose ranges are not mostly in north america


tenyear_subset = data[(data['year'] <= 2015) & (data['year'] >= 2005)]
fiveyear_subset = data[(data['year'] <= 2015) & (data['year'] >= 2010)]

tenyear_count = []
for site, site_data in tenyear_subset.groupby('site'):
    count = len(site_data['year'].unique())
    tenyear_count.append([site, count])
tenyear = pd.DataFrame(tenyear_count, columns = ['site', 'counts'])
    
fiveyear_count = []
for site, site_data in fiveyear_subset.groupby('site'):
    count = len(site_data['year'].unique())
    fiveyear_count.append([site, count])
fiveyear = pd.DataFrame(fiveyear_count, columns = ['site', 'counts'])

plt.hist(tenyear['counts'])
plt.savefig('figures/sensitivity/tenyear_counts.png')

plt.hist(fiveyear['counts'])
plt.savefig('figures/sensitivity/fiveyear_counts.png')

def filter_survey_data (data, start_year, end_year, min_years_surveyed, year_column, site_column):
    data = data[(data[year_column] <= end_year) & (data[year_column] >= start_year)]
    sites = []
    for site, site_data in data.groupby(site_column):
        count = len(site_data[year_column].unique())
        if count >= min_years_surveyed:
            sites.append(site)
    return data[data[site_column].isin(sites)]

test = filter_survey_data(data, 2005, 2015, 5, 'year', 'site')
len(tenyear[tenyear.counts >= 5]['site'].unique())


def plot_sites_by_characteristic(dataframe, lat_col, long_col, title=None, char_column=None, bins=None, dataframe2=None, lat_col2=None, long_col2=None):
    map = Basemap(projection='merc',llcrnrlat=23.5,urcrnrlat=57, llcrnrlon=-140,urcrnrlon=-50,lat_ts=20,resolution='l')
    map.drawcoastlines(linewidth = 1.25)
    plt.title(title)
    
    if not char_column:    
        lats = dataframe[lat_col]
        longs = dataframe[long_col]
        x,y = map(longs.values,lats.values)
        map.plot(x, y, ls='', marker='o', markersize=4)

    if char_column:
        blues = sns.color_palette("Blues", n_colors=bins)
        dataframe['quantile'] = pd.qcut(dataframe[char_column], bins)
        grouped = dataframe.groupby('quantile')
        
        i= -1
        for groupname, groupdata, in grouped:
            i = i + 1
            colors = blues[i]
            lats = groupdata["lat"]
            longs = groupdata["long"]
            x,y = map(longs.values,lats.values)
            map.plot(x, y, ls='', marker='o', color=colors, markersize=4)
    plt.hold(True)
    if lat_col2:    
        lats = dataframe2[lat_col2]
        longs = dataframe2[long_col2]
        x,y = map(longs.values,lats.values)
        map.plot(x, y, ls='', marker='o', markersize=4, color='brown')    

def get_hotspots(data, richness_column, cell=False):
    sort = data.sort([richness_column], ascending=False)
    if cell is False:
        hotspots = sort.head(int(round(len(sort)*0.05)))
    else:
        num_hotspots = int(round(0.05 * (len(sort)-sort[richness_column].isnull().sum())))
        hotspots = sort.head(num_hotspots)
        hotspots['hotspot'] = [1]*num_hotspots   
    return hotspots

richness_by_site = macroecotools.richness_in_group(data, 
                                                   ['site', 'lat', 'long'], ['species'])

#plot according to richness at site
hotspot_sites = get_hotspots(richness_by_site, 'richness')
plot_sites_by_characteristic(richness_by_site, lat_col='lat', long_col='long', title='survey richness', 
                             char_column='richness', bins=10, dataframe2=hotspot_sites, lat_col2='lat', long_col2='long')

def plot_range_of_survey_years(dataframe, range_years_surveyed, start_year, end_year, year_column, site_column, species_column, lat_col, long_col):
    j = 1
    plot_num = len(range_years_surveyed)        
    fig = plt.figure()
    for i in range_years_surveyed:
        filtered = filter_survey_data(data= dataframe, start_year = start_year, end_year = end_year, min_years_surveyed= i, year_column = year_column, site_column = site_column)
        filtered_richness = macroecotools.richness_in_group(filtered, [site_column, lat_col, long_col], [species_column])
        hotspots = get_hotspots(filtered_richness, 'richness')

        ax = fig.add_subplot(plot_num, 1, j)
        plot_sites_by_characteristic(filtered_richness, lat_col = lat_col, long_col = long_col, char_column='richness', bins=10, dataframe2=hotspots, lat_col2='lat', long_col2='long')
        j = j+1
    plt.show()
    
survey_years = [5, 9, 11]
plot_range_of_survey_years(data, survey_years, 2005, 2015, 'year', 'site', 'species', 'lat', 'long')