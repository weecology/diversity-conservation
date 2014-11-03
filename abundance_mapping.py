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

#determine number of sites
data_site = data.groupby('site')
sites = []
for site, site_data in data_site:
    sites.append(site)
print len(sites)

#plot sites
map = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=71, llcrnrlon=-170,urcrnrlon=-50,lat_ts=20,resolution='l')
map.drawcoastlines(linewidth = 1.25)
lats = data["lat"]
longs = data["long"]
x,y = map(longs.values,lats.values)
map.plot(x, y, ls='', marker='o')
plt.show()

#plot according to richness at site
richness_by_site = macroecotools.richness_in_group(data, ['site'], ['species'])