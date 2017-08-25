from __future__ import division
from __future__ import print_function
import numpy as np
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

included_species = pd.read_csv('mapping_data/included_species_ids.csv')
#clean up excluded families
included_species = included_species[(included_species['AOU'] > 2880)]
included_species = included_species[(included_species['AOU'] < 3650) | (included_species['AOU'] > 3810)]
included_species = included_species[(included_species['AOU'] < 3900) | (included_species['AOU'] > 3910)]
included_species = included_species[(included_species['AOU'] < 4160) | (included_species['AOU'] > 4210)]
included_species = included_species[(included_species['AOU'] != 7010)]
included_species.rename(columns= {'AOU':'species'}, inplace = True)
AOU_list = pd.DataFrame(included_species['species'])
sisid_list = pd.DataFrame(included_species['sisid'])

data = pd.read_csv('mapping_data/bbs_species_2016.csv', delimiter=',', sep='\s*,\s*')
data.rename(columns = {'site_id':'site'}, inplace = True)
data.rename(columns = {'species_id':'species'}, inplace = True)
data = pd.merge(data, AOU_list, how='inner', on=['species']) #exclude species whose ranges are not mostly in north america


tenyear_subset = data[(data['year'] <= 2015) & (data['year'] >= 2005)]
fiveyear_subset = data[(data['year'] <= 2015) & (data['year'] >= 2010)]

tenyear_count = []
for site, site_data in tenyear_subset.groupby('site'):
    count = len(site_data['year'].unique())
    tenyear_count.append([site, count])
tenyear = pd.DataFrame(tenyear_count, columns = ['site', 'count'])
    
fiveyear_count = []
for site, site_data in fiveyear_subset.groupby('site'):
    count = len(site_data['year'].unique())
    fiveyear_count.append([site, count])
fiveyear = pd.DataFrame(fiveyear_count, columns = ['site', 'count'])

plt.hist(tenyear['count'])
plt.hist(fiveyear['count'])