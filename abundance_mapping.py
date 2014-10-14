from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import colorsys
from numpy import log10
import pandas as pd
import numpy
from mpl_toolkits.basemap import Basemap
import os
import macroecotools

data = pd.read_csv('bbs_abundances_by_site.csv', delimiter=',')

data_site = data.groupby('site')

#determine number of sites
sites = []
for site, site_data in data_site:
    sites.append(site)
print len(sites)

macroecotools.richness_in_group(data, ['site'], ['species'])
