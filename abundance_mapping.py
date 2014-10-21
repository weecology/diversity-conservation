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

data_site = data.groupby('site')

#determine number of sites
sites = []
for site, site_data in data_site:
    sites.append(site)
print len(sites)

richness_by_site = macroecotools.richness_in_group(data, ['site'], ['species'])

map = Basemap(projection='merc',llcrnrlat=-57,urcrnrlat=71, llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='l')
map.drawcoastlines(linewidth = 1.25)
lats = data["lat"]
longs = data["long"]
x,y = map(longs,lats)
map.plot(x,y, ls='', marker=markers[i], markerfacecolor=colors[i], markeredgewidth=0.25, markersize=markersizes)