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

plt.figure()
plt.hist(uniq_prop, bins=20)

plt.show()