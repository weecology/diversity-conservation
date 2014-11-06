from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import colorsys
from numpy import log10
import pandas as pd
import numpy
from mpl_toolkits.basemap import Basemap
import macroecotools

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

map = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=71, llcrnrlon=-170,urcrnrlon=-50,lat_ts=20,resolution='l')
map.drawcoastlines(linewidth = 1.25)
lats = richness_by_site["lat"]
longs = richness_by_site["long"]
x,y = map(longs.values,lats.values)
blues = np.linspace(0, 225, num=10)
richness_by_site['quantile'] = pd.qcut(richness_by_site['richness'], 10)
grouped = richness_by_site.groupby('quantile')

i=-1
for groupname, groupdata, in grouped:
    i = i + 1
    colors = (0, 0, blues[i])
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

