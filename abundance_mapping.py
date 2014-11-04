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

#plot rare species
data_rare = data.groupby(['species', 'site']).count()
data_species = data.groupby('species')

total_sites = len(np.unique(data['site']))

rarity_prop = []

for species, species_data in data_species:
    occurence_sites = len(species_data['site'])
    proportion = occurence_sites/total_sites
    rarity_prop.append([species, proportion])
sp_rarity = pd.DataFrame(rarity_prop, columns=['species', 'proportion'])
data_w_proportion = pd.merge(sp_rarity, data, on='species')
