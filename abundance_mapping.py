from __future__ import division
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import colorsys
from numpy import log10
import pandas as pd
import numpy
from mpl_toolkits.basemap import Basemap
import macroecotools
import seaborn as sns
import random
import glob
import os

data = pd.read_csv('bbs_abundances_by_site.csv', delimiter=',')

#plot sites
def plot_sites_by_characteristic(dataframe, lat_col, long_col, char_column=None, bins=None):
    plt.figure()
    map = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=71, llcrnrlon=-170,urcrnrlon=-50,lat_ts=20,resolution='l')
    map.drawcoastlines(linewidth = 1.25)
    
    if not char_column:    
        lats = dataframe[lat_col]
        longs = dataframe[long_col]
        x,y = map(longs.values,lats.values)
        map.plot(x, y, ls='', marker='o')

    if char_column:
        blues = sns.color_palette("Blues", n_colors=bins)
        dataframe['quantile'] = pd.qcut(dataframe[char_column], bins)
        grouped = dataframe.groupby('quantile')
        
        i=-1
        for groupname, groupdata, in grouped:
            i = i + 1
            colors = blues[i]
            lats = groupdata["lat"]
            longs = groupdata["long"]
            x,y = map(longs.values,lats.values)
            map.plot(x, y, ls='', marker='o', color=colors)


plot_sites_by_characteristic(data, 'lat', 'long')

#plot according to richness at site
richness_by_site = macroecotools.richness_in_group(data, ['site', 'lat', 'long'], ['species'])

plot_sites_by_characteristic(richness_by_site, lat_col='lat', long_col='long', char_column='richness', bins=10)

#plot rare species
data_species = data.groupby('species')

total_sites = len(np.unique(data['site']))

rarity_prop = []

for species, species_data in data_species:
    occurence_sites = len(species_data['site'])
    proportion = occurence_sites/total_sites
    rarity_prop.append([species, proportion])
sp_rarity = pd.DataFrame(rarity_prop, columns=['species', 'proportion'])
data_w_proportion = pd.merge(sp_rarity, data, on='species')

data_w_proportion_species = data_w_proportion.groupby('species')

uniq_prop = []
for species, species_data in data_w_proportion_species:
    mean=np.mean(species_data['proportion'])
    uniq_prop.append(mean)
med = np.median(uniq_prop)
    
data_rare = data_w_proportion[data_w_proportion['proportion'] < med]

plot_sites_by_characteristic(data_rare, lat_col='lat', long_col='long')

plt.show()

#grid sampling

def get_sites_by_grid(dataframe, site_col, lat_col, long_col, band_width, sites_in_cell):
    dataframe = dataframe[['site', 'lat', 'long']].drop_duplicates()
    min_lat = min(dataframe[lat_col])
    max_lat = max(dataframe[lat_col])
    min_long = min(dataframe[long_col])
    max_long = max(dataframe[long_col])
    band_degrees = (band_width/40000)*360
    lat_start = min_lat - 0.001
    lat_end = lat_start
    data_selection = pd.DataFrame()
    while lat_end < max_lat:
        long_start = min_long - 0.001
        long_end = long_start
        lat_end = lat_end + band_degrees
        while long_end < max_long:
            long_end = long_end + band_degrees
            data_sub = dataframe[(dataframe[lat_col] > lat_start) & (dataframe[lat_col] < lat_end) & (dataframe[long_col] > long_start) & (dataframe[long_col] < long_end)]
            long_start = long_end
            if len(data_sub['site']) >= sites_in_cell:
                selection = data_sub.ix[random.sample(data_sub.index, sites_in_cell)]
                data_selection = data_selection.append(selection)
        lat_start = lat_end
    return data_selection

if os.path.isfile('selected_sites.csv') == True:
    selected_sites = pd.read_csv('selected_sites.csv', delimiter=',')
    print ('yes')
else:
    selected_sites = get_sites_by_grid(data_w_proportion, 'site', 'lat', 'long', 100, 3)
    selected_sites.to_csv('selected_sites.csv')