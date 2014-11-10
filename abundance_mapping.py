from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import colorsys
from numpy import log10
import pandas as pd
import numpy
from mpl_toolkits.basemap import Basemap
import macroecotools
import seaborn as sns

data = pd.read_csv('bbs_abundances_by_site.csv', delimiter=',')

#plot sites
map = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=71, llcrnrlon=-170,urcrnrlon=-50,lat_ts=20,resolution='l')
map.drawcoastlines(linewidth = 1.25)
lats = data["lat"]
longs = data["long"]
x,y = map(longs.values,lats.values)
map.plot(x, y, ls='', marker='o')
plt.show()

#plot according to richness at site
richness_by_site = macroecotools.richness_in_group(data, ['site', 'lat', 'long'], ['species'])

def plot_sites_by_characteristic(dataframe, lat_col, long_col, char_column=None, bins=None):
    map = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=71, llcrnrlon=-170,urcrnrlon=-50,lat_ts=20,resolution='l')
    map.drawcoastlines(linewidth = 1.25)
    
    if not char_column:    
        lats = dataframe[lat_col]
        longs = dataframe[long_col]
        x,y = map(longs.values,lats.values)

    if char_column:
        blues = sns.color_palette("Blues", n_colors=bins)
        dataframe['quantile'] = pd.qcut(dataframe['char_column'], bins)
        grouped = dataframe.groupby('quantile')
        
        i=-1
        for groupname, groupdata, in grouped:
            i = i + 1
            colors = blues[i]
            lats = groupdata["lat"]
            longs = groupdata["long"]
            x,y = map(longs.values,lats.values)
            print colors
            map.plot(x, y, ls='', marker='o', color=colors)
    plt.show()

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

def sum_or_not(list_of_numbers, tosum=True):
    if tosum:
        return sum(list_of_numbers)