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

#plot sites
def plot_sites_by_characteristic(dataframe, lat_col, long_col, title=None, char_column=None, bins=None):
    plt.figure()
    map = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=71, llcrnrlon=-170,urcrnrlon=-50,lat_ts=20,resolution='l')
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

#plot rare species
def get_rarity_proportion(dataframe, species_column, site_column):
    data_species = dataframe.groupby(species_column)
    total_sites = len(np.unique(dataframe[site_column]))
    rarity_prop = []
    for species, species_data in data_species:
        occurence_sites = len(species_data[site_column])
        proportion = occurence_sites/total_sites
        rarity_prop.append([species, proportion])
    sp_rarity = pd.DataFrame(rarity_prop, columns=[species_column, 'proportion'])
    data_w_proportion = pd.merge(sp_rarity, dataframe, on=species_column)
    return data_w_proportion


def get_median_rarity_proportion(dataframe, species_column, proportion_column):
    dataframe_species = dataframe.groupby(species_column)
    uniq_prop = []
    for species, species_data in dataframe_species:
        mean=np.mean(species_data[proportion_column])
        uniq_prop.append(mean)
    med = np.median(uniq_prop)
    return med

#find centroid of cell
def get_centroid(points):  
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    cent_lat = sum(x) / len(points)
    cent_long = sum(y) / len(points)
    return cent_lat, cent_long


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
    cellid = 0
    data_selection = pd.DataFrame()
    centroid_coordinates = pd.DataFrame()
    while lat_end < max_lat:
        long_start = min_long - 0.001
        long_end = long_start
        lat_end = lat_end + band_degrees
        while long_end < max_long:
            long_end = long_end + band_degrees
            data_sub = dataframe[(dataframe[lat_col] > lat_start) & (dataframe[lat_col] < lat_end) & (dataframe[long_col] > long_start) & (dataframe[long_col] < long_end)]
            corners = [(lat_start, long_start), (lat_start, long_end), (lat_end, long_start), (lat_end, long_end)]
            cent_lat, cent_long = get_centroid(corners)
            cellid = cellid + 1
            centroid_coordinates = centroid_coordinates.append([(cent_lat, cent_long, cellid)])
            long_start = long_end
            if len(data_sub['site']) >= sites_in_cell:
                selection = data_sub.ix[random.sample(data_sub.index, sites_in_cell)]
                selection['cellid'] = cellid                
                data_selection = data_selection.append(selection)
        lat_start = lat_end
    centroid_coordinates.columns = ['cent_lat', 'cent_long', 'cellid']
    cell_info = pd.merge(centroid_coordinates, data_selection, how = 'left', on = ['cellid'])
    return cell_info

data = pd.read_csv('bbs_abundances_by_site.csv', delimiter=',')
cell_info = get_sites_by_grid(data, 'site', 'lat', 'long', 100, 3)

#SURVEY DATA
data = pd.read_csv('bbs_abundances_by_site.csv', delimiter=',')

#plot sites
plot_sites_by_characteristic(data, 'lat', 'long', title='sites')

#plot according to richness at site
richness_by_site = macroecotools.richness_in_group(data, ['site', 'lat', 'long'], ['species'])
plot_sites_by_characteristic(richness_by_site, lat_col='lat', long_col='long', title='survey richness', char_column='richness', bins=10)

#plot sites with rare species, not adjusted for spatial bias
data_w_proportion = get_rarity_proportion(data, 'species', 'site')
median_rarity = get_median_rarity_proportion(data_w_proportion, 'species', 'proportion')
data_rare = data_w_proportion[data_w_proportion['proportion'] < median_rarity]
#plot_sites_by_characteristic(data_rare, lat_col='lat', long_col='long')

#plot selected sites to eliminate spatial bias
if os.path.isfile('selected_sites.csv') == True:
    selected_sites = pd.read_csv('selected_sites.csv', delimiter=',')
    print ('yes')
else:
    selected_sites = get_sites_by_grid(data, 'site', 'lat', 'long', 100, 3)
    selected_sites.to_csv('selected_sites.csv')
    
plot_sites_by_characteristic(selected_sites, lat_col='lat', long_col='long', title='grid-selected sites')

#plot sites with rare species
data_from_selected_sites = pd.merge(selected_sites, data, how='left', on=['site', 'lat', 'long'])
selected_w_proportion = get_rarity_proportion(data_from_selected_sites, 'species', 'site')
selected_median = get_median_rarity_proportion(selected_w_proportion, 'species', 'proportion')
selected_rare = selected_w_proportion[selected_w_proportion['proportion'] < selected_median]

plot_sites_by_characteristic(selected_rare, 'lat', 'long', title='sites with rare species (survey)')

#plot sites according to richness of rare species
selected_rare = selected_rare.drop('proportion', 1)
rarity_richness_by_site = macroecotools.richness_in_group(selected_rare, ['site', 'lat', 'long'], ['species'])
plot_sites_by_characteristic(rarity_richness_by_site, lat_col='lat', long_col='long', char_column='richness', bins=2, title='richness of rare species (survey)')

#75th percentile richness
data_rare_high = rarity_richness_by_site[rarity_richness_by_site['richness'] > 5]
#plot_sites_by_characteristic(data_rare_high, lat_col='lat', long_col='long', char_column='richness', bins=4)


#RANGE DATA
range_map = pd.read_csv('rangemap_species.csv')
range_map = range_map.sort('site')
range_abun = macroecotools.richness_in_group(range_map, ['site', 'lat', 'long'], ['sisid'])

range_selected = pd.merge(selected_sites, range_map, how='left', on=['site', 'lat', 'long'])
range_prop = get_rarity_proportion(range_selected, 'sisid', 'site')
range_median = get_median_rarity_proportion(range_prop, 'sisid', 'proportion')
range_rare = range_prop[range_prop['proportion'] < range_median]

#plot range map abundance at site points
plot_sites_by_characteristic(range_abun, lat_col='lat', long_col='long', char_column='richness', bins=10, title="range map richness")

#plot sites with rare species
plot_sites_by_characteristic(range_rare, 'lat', 'long', title='sites with rare species (range)')

#plot sites according to richness of rare species
range_rare = range_rare.drop('proportion', 1)
range_rarity_richness = macroecotools.richness_in_group(range_rare, ['site', 'lat', 'long'], ['sisid'])
plot_sites_by_characteristic(range_rarity_richness, lat_col='lat', long_col='long', char_column='richness', bins=2, title='rare species richenss (range)')

#75th percentile richness
#range_rare_high = range_rarity_richness[range_rarity_richness['richness'] > 15]
#plot_sites_by_characteristic(range_rare_high, lat_col='lat', long_col='long', char_column='richness', bins=4)


#Range map rarity definition
range_area = pd.read_csv('species_area.csv')
range_area_uniq = range_area.groupby('sisid', as_index=False).sum()
range_area_full = pd.merge(range_area_uniq, range_map, on=['sisid'])
rare_range = range_area_full[range_area_full['shape_area'] < np.median(np.unique(range_area_full['shape_area']))]
rare_range_full = pd.merge(rare_range, range_abun, on=['site', 'lat', 'long'])

plot_sites_by_characteristic(rare_range_full, 'lat', 'long', char_column='richness', bins=10, title='sites with small range species')