from __future__ import division
from __future__ import print_function
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
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
import pdb

#plot sites
def plot_sites_by_characteristic(dataframe, lat_col, long_col, title=None, char_column=None, bins=None, dataframe2=None, lat_col2=None, long_col2=None):
    plt.figure()
    map = Basemap(projection='merc',llcrnrlat=23.5,urcrnrlat=57, llcrnrlon=-140,urcrnrlon=-50,lat_ts=20,resolution='l')
    map.drawcoastlines(linewidth = 1.25)
    plt.title(title)

    if not char_column:
        lats = dataframe[lat_col]
        longs = dataframe[long_col]
        x,y = map(longs.values,lats.values)
        map.plot(x, y, ls='', marker='o', markersize=4)

    if char_column:
        blues = sns.color_palette("Blues", n_colors=bins)
        dataframe['quantile'] = pd.qcut(dataframe[char_column], bins, duplicates='drop')
        grouped = dataframe.groupby('quantile')

        i= -1
        for groupname, groupdata, in grouped:
            i = i + 1
            colors = blues[i]
            lats = groupdata["lat"]
            longs = groupdata["long"]
            x,y = map(longs.values,lats.values)
            map.plot(x, y, ls='', marker='o', color=colors, markersize=2)
    plt.hold(True)
    if lat_col2:
        lats = dataframe2[lat_col2]
        longs = dataframe2[long_col2]
        x,y = map(longs.values,lats.values)
        map.plot(x, y, ls='', marker='o', markersize=2, color='brown')
    plt.show()

#plot rare species
def get_rarity_proportion(dataframe, species_column, site_column):
    data_species = dataframe.groupby(species_column)
    total_sites = len(np.unique(dataframe[site_column]))
    rarity_prop = []
    for species, species_data in data_species:
        occurence_sites = len(np.unique(species_data[site_column]))
        proportion = occurence_sites/total_sites
        rarity_prop.append([species, proportion])
    sp_rarity = pd.DataFrame(rarity_prop, columns=[species_column, 'proportion'])
    data_w_proportion = pd.merge(sp_rarity, dataframe, how='right', on=species_column)
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
    dataframe = dataframe[[site_col, lat_col, long_col]].drop_duplicates()
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
            if len(data_sub[site_col]) >= sites_in_cell:
                selection = data_sub.ix[random.sample(data_sub.index, sites_in_cell)]
                selection['cellid'] = cellid
                data_selection = data_selection.append(selection)
        lat_start = lat_end
    centroid_coordinates.columns = ['cent_lat', 'cent_long', 'cellid']
    cell_info = pd.merge(centroid_coordinates, data_selection, how = 'left', on = ['cellid'])
    return cell_info

def get_hotspots(data, richness_column, cell=False):
    sort = data.sort_values([richness_column], ascending=False)
    if cell is False:
        hotspots = sort.head(int(round(len(sort)*0.05)))
    else:
        num_hotspots = int(round(0.05 * (len(sort)-sort[richness_column].isnull().sum())))
        hotspots = sort.head(num_hotspots)
        hotspots['hotspot'] = [1]*num_hotspots
    return hotspots



#species for whome the majority of their range is in North America
included_species = pd.read_csv('data/mapping_data/included_species_ids.csv')


#SURVEY DATA
#formatting
data_non_exclude = pd.read_csv('data/mapping_data/bbs_species_2016.csv', delimiter=',', sep='\s*,\s*')
data_non_exclude.rename(columns = {'site_id':'site'}, inplace = True)
data_non_exclude.rename(columns = {'species_id':'species'}, inplace = True)
data_non_exclude = data_non_exclude[(data_non_exclude['year'] <= 2015) & (data_non_exclude['year'] >= 2005)]

def plot_hotspots(data, title, bins = 10, species_col = 'species', type = "plot"):
  richness_by_site = macroecotools.richness_in_group(data, ['site', 'lat', 'long'], [species_col])
  hotspot_sites = get_hotspots(richness_by_site, 'richness')

  if type == "plot":
    plot_sites_by_characteristic(richness_by_site, lat_col='lat', long_col='long',
                                title=title, char_column='richness', bins=bins, dataframe2=hotspot_sites, lat_col2='lat', long_col2='long')
  if type == "hotspots":
      return(hotspot_sites)

#plot richness map with non-North America majority species included
plot_hotspots(data_non_exclude, 'survey richness all species')
plt.savefig('figures/survey_site_non_exclude.png')

#exclude non-North America majority
data = pd.merge(data_non_exclude, included_species[['AOU']], how='inner', left_on=['species'], right_on=['AOU']) #exclude species whose ranges are not mostly in north america
plot_hotspots(data, 'survey richness species excluded')
plt.savefig('figures/survey_site.png')






##estimates
#data_est = pd.read_csv('site_biodiversity_estimates.csv', delimiter=',')
#data_est['site'] = richness_by_site['site']
#data_est = pd.merge(richness_by_site[['site', 'lat', 'long']], data_est, how='left', on='site')
#hotspot_sites_est = get_hotspots(data_est, 'Jack1ab') 
#plot_sites_by_characteristic(data_est, lat_col='lat', long_col='long', title='survey richness estimates', char_column='Jack1ab', bins=10, dataframe2=hotspot_sites_est, lat_col2='lat', long_col2='long')

#plot sites with rare species, not adjusted for spatial bias
data_w_proportion = get_rarity_proportion(data, 'species', 'site')
median_rarity = get_median_rarity_proportion(data_w_proportion, 'species', 'proportion')
data_rare = data_w_proportion[data_w_proportion['proportion'] < median_rarity]
#plot_sites_by_characteristic(data_rare, lat_col='lat', long_col='long')

if os.path.isfile('selected_sites.csv') == True:
    selected_sites = pd.read_csv('data/mapping_data/selected_sites.csv', delimiter=',')
    print ('yes')
else:
    selected_sites = get_sites_by_grid(data, 'site', 'lat', 'long', 100, 3)
    selected_sites.to_csv('selected_sites.csv')
    
#find rare species
data_from_selected_sites = pd.merge(selected_sites, data, how='left', on=['site', 'lat', 'long'])
selected_w_proportion = get_rarity_proportion(data_from_selected_sites, 'species', 'site')
selected_median = get_median_rarity_proportion(selected_w_proportion, 'species', 'proportion')
selected_rare = selected_w_proportion[selected_w_proportion['proportion'] < selected_median]
selected_rare = selected_rare.drop('proportion', 1)

#plot sites according to richness of rare species
plot_hotspots(selected_rare, bins = 5, title = "survey rarity")
plt.savefig('figures/survey_site_rare_new.png')

#RANGE DATA
#richness
range_map_inc = pd.read_csv('data/mapping_data/rangemap_species_2016-order.csv', usecols=['site', 'sisid'])
range_map_inc = pd.merge(range_map_inc, data[['site', 'lat', 'long']].drop_duplicates(), on = 'site', how = 'left')
range_map_inc = range_map_inc.sort_values('site')
range_map = pd.merge(range_map_inc, included_species[["sisid"]], how='inner', on=['sisid'])
#range_abun = macroecotools.richness_in_group(range_map, ['site', 'lat', 'long'], ['sisid'])

plot_hotspots(range_map, title = "range map richness", species_col = 'sisid')
plt.savefig('figures/range_site.png')

#rarity
range_selected = pd.merge(selected_sites, range_map, how='left', on=['site', 'lat', 'long'])
range_prop = get_rarity_proportion(range_selected, 'sisid', 'site')
range_median = get_median_rarity_proportion(range_prop, 'sisid', 'proportion')
range_rare = range_prop[range_prop['proportion'] < range_median]
range_rare = range_rare.drop('proportion', 1)

#plot sites according to richness of rare species
plot_hotspots(range_rare, title = 'range rarity', species_col = 'sisid')
plt.savefig('figures/range_site_rare.png')

##Range map rarity definition
#range_area = pd.read_csv('species_area.csv')
#range_area_uniq = range_area.groupby('sisid', as_index=False).sum()
#range_area_full = pd.merge(range_area_uniq, range_map, on=['sisid'])
#rare_range = range_area_full[range_area_full['shape_area'] < np.median(np.unique(range_area_full['shape_area']))]
#rare_range_full = pd.merge(rare_range, range_abun, on=['site', 'lat', 'long'])

#plot_sites_by_characteristic(rare_range_full, 'lat', 'long', char_column='richness', bins=10, title='sites with small range species')

#CELL MAPPING
def get_unique_cell_richness (data, cell_id_column, cell_lat_column, cell_long_column, speciesid_column):
    uniq_cell_abun = pd.DataFrame()
    for cell, cell_data in data.groupby(cell_id_column):
        count = len(np.unique(cell_data[speciesid_column]))
        uniq_cell_abun = uniq_cell_abun.append([(cell, count)])
    uniq_cell_abun.columns = [cell_id_column, 'total_richness']
    uniq_cell_abun = uniq_cell_abun[uniq_cell_abun['total_richness'] > 1]
    uniq_cell_abun_loc = pd.merge(selected_sites[[cell_lat_column, cell_long_column, cell_id_column]].drop_duplicates(), uniq_cell_abun, how='left', on=[cell_id_column])
    return uniq_cell_abun_loc

def plot_cell_feature (data, cell_id_column, cell_lat_column, cell_long_column, richness_column, title=None, second_feature_data=None, second_feature_column=None):
    lats = np.asarray(np.unique(data[cell_lat_column]))
    lons = np.asarray(np.unique(data[cell_long_column]))
    lons, lats = np.meshgrid(lons,lats)
    
    richness = np.array(data[richness_column])
    richness.shape = (len(np.unique(lats)), len(np.unique(lons)))
    richness_mask = ma.masked_where(np.isnan(richness),richness)
    
    if second_feature_column:
        second_feature = np.array(second_feature_data[second_feature_column])
        second_feature.shape = (len(np.unique(lats)), len(np.unique(lons)))
        second_feature_data_mask = ma.masked_where(np.isnan(second_feature),second_feature)
        
    fig = plt.figure()
    m = Basemap(projection='merc',llcrnrlat=23.5,urcrnrlat=57, llcrnrlon=-140,urcrnrlon=-50,lat_ts=20,resolution='l')
    m.drawcoastlines(linewidth = 1.25)
    if np.nanmin(richness) < 20:
        vmin=0
    else:
        vmin=round(np.nanmin(richness)-20, -1)
    im1 = m.pcolormesh(lons,lats,richness_mask,shading='flat',cmap=plt.cm.Blues,latlon=True, vmin=vmin)
    if second_feature_column:
        im2 = m.pcolormesh(lons,lats,second_feature_data_mask,shading='flat',cmap=plt.cm.RdYlBu,latlon=True)
    cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
    plt.title(title)

#observed richness
cell_site_species = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid', 'site']], data, how='left', on=['site'])
uniq_cell_abun = get_unique_cell_richness(cell_site_species, 'cellid', 'cent_lat', 'cent_long', 'species')   
obs_hotspot_cell = get_hotspots(uniq_cell_abun, 'total_richness', cell=True)
all_hotspot_cell = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid']].drop_duplicates(), obs_hotspot_cell, how='left', on=['cellid', 'cent_lat', 'cent_long'])
plot_cell_feature(uniq_cell_abun, 'cellid', 'cent_lat', 'cent_long', 'total_richness', title='Observed Survey Richness with Hotspots', second_feature_data=all_hotspot_cell, second_feature_column='hotspot')
plt.savefig('figures/survey_cell.png')

#observed richness including all observed species
cell_site_species_inc = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid', 'site']], data, how='left', on=['site'])
uniq_cell_abun_inc = get_unique_cell_richness(cell_site_species_inc, 'cellid', 'cent_lat', 'cent_long', 'species')
obs_hotspot_cell_inc = get_hotspots(uniq_cell_abun_inc, 'total_richness', cell=True)
all_hotspot_cell_inc = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid']].drop_duplicates(), obs_hotspot_cell_inc, how='left', on=['cellid', 'cent_lat', 'cent_long'])
plot_cell_feature(uniq_cell_abun_inc, 'cellid', 'cent_lat', 'cent_long', 'total_richness', title='Observed Survey Richness with Hotspots non exclude', second_feature_data=all_hotspot_cell_inc, second_feature_column='hotspot')
plt.savefig('figures/survey_cell_non_exclude.png')

##estimated richness
#cell_bio_est = pd.read_csv("cell_estimates.csv", delimiter=",")
#uniq_cell = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid', 'site']], richness_by_site[['site', 'richness']], how='right', on=['site'])
#cells = np.unique(uniq_cell['cellid'].dropna())
#cells = pd.DataFrame(cells)
#cells.columns=['cellid']
#cell_bio_est = cell_bio_est.join(cells, lsuffix='_left', rsuffix='_right')
#cell_abun_est = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid']].drop_duplicates(), cell_bio_est, how='left', on=['cellid'])
#est_hotspot_cells = get_hotspots(cell_abun_est, 'Jack1ab', cell=True)
#all_est_hotspot_cells = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid']].drop_duplicates(), est_hotspot_cells, how='left', on=['cellid', 'cent_lat', 'cent_long'])
#plot_cell_feature(cell_abun_est, 'cellid', 'cent_lat', 'cent_long', 'Jack1ab', title='Estimated Survey Richness with Hotspots', second_feature_data=all_est_hotspot_cells, second_feature_column='hotspot')

#range map
cell_range_species = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid', 'site']], range_map, how='left', on=['site'])
uniq_range_cell = get_unique_cell_richness(cell_range_species, 'cellid', 'cent_lat', 'cent_long', 'sisid')
range_hotspot_cells = get_hotspots(uniq_range_cell, 'total_richness', cell=True)
all_range_hotspot_cells = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid']].drop_duplicates(), range_hotspot_cells, how='left', on=['cellid', 'cent_lat', 'cent_long'])
plot_cell_feature(uniq_range_cell, 'cellid', 'cent_lat', 'cent_long', 'total_richness', title='Range Richness Hotspots', second_feature_data=all_range_hotspot_cells, second_feature_column='hotspot')
plt.savefig('figures/range_cell.png')

#range map richness for all species
cell_range_species_inc = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid', 'site']], range_map_inc, how='left', on=['site'])
uniq_range_cell_inc = get_unique_cell_richness(cell_range_species_inc, 'cellid', 'cent_lat', 'cent_long', 'sisid')
range_hotspot_cells_inc = get_hotspots(uniq_range_cell_inc, 'total_richness', cell=True)
all_range_hotspot_cells_inc = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid']].drop_duplicates(), range_hotspot_cells_inc, how='left', on=['cellid', 'cent_lat', 'cent_long'])
plot_cell_feature(uniq_range_cell, 'cellid', 'cent_lat', 'cent_long', 'total_richness', title='Range Richness Hotspots non exclude', second_feature_data=all_range_hotspot_cells_inc, second_feature_column='hotspot')
plt.savefig('figures/range_cell_non_exclude.png')

#rare species
rare_survey_cell = get_unique_cell_richness(selected_rare[['cent_lat', 'cent_long', 'cellid', 'site', 'lat', 'long', 'abundance', 'species']], 'cellid', 'cent_lat', 'cent_long', 'species')
rare_survey_hotspot_cells = get_hotspots(rare_survey_cell, 'total_richness', cell=True)
all_rare_survey_hotspot_cells = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid']].drop_duplicates(), rare_survey_hotspot_cells, how='left', on=['cellid', 'cent_lat', 'cent_long'])
plot_cell_feature(rare_survey_cell, 'cellid', 'cent_lat', 'cent_long', 'total_richness', 'Rare Survey Richness with Hotspots', second_feature_data=all_rare_survey_hotspot_cells, second_feature_column='hotspot')
plt.savefig('figures/survey_cell_rare.png')

rare_range_cell = get_unique_cell_richness(range_rare, 'cellid', 'cent_lat', 'cent_long', 'sisid')
rare_range_hotspot_cells = get_hotspots(rare_range_cell, 'total_richness', cell=True)
all_rare_range_hotspot_cells = pd.merge(selected_sites[['cent_lat', 'cent_long', 'cellid']].drop_duplicates(), rare_range_hotspot_cells, how='left', on=['cellid', 'cent_lat', 'cent_long'])
plot_cell_feature(rare_range_cell,'cellid', 'cent_lat', 'cent_long', 'total_richness', 'Rare Range Map Richness', second_feature_data=all_rare_range_hotspot_cells, second_feature_column='hotspot')
plt.savefig('figures/range_cell_rare.png')

#one to one plotting
uniq_cell_abun.columns = ['cent_lat', 'cent_long', 'cellid', 'survey_richness']
uniq_range_cell.columns = ['cent_lat', 'cent_long', 'cellid', 'range_richness']
rich_comp = pd.merge(uniq_cell_abun, uniq_range_cell, how='left', on=['cent_lat', 'cent_long', 'cellid'])
rich_comp['line'] = rich_comp['survey_richness']
ax = rich_comp.plot(kind='scatter', x='survey_richness', y='range_richness')
plt.plot(rich_comp['survey_richness'], rich_comp['line'], 'k-')

#RARITY WEIGHTED RICHNESS
#survey data
def rwr_priority_sites(data, species_col):
    rwr_data = pd.DataFrame()
    for species, species_data in data.groupby(species_col):
        rwr_data = rwr_data.append([(species, 1/(len(species_data)))])
    rwr_data.columns = [species_col,'rwr']
    rwr_data = pd.merge(rwr_data, data, how='right', on=species_col)
    
    rwr_sites = pd.DataFrame()
    for site, site_data in rwr_data.groupby('site'):
        rwr_sum = site_data['rwr'].sum()
        rwr_sites = rwr_sites.append([(site, rwr_sum)])
    rwr_sites.columns = ['site', 'rwr']
    rwr_sites = rwr_sites.sort_values('rwr', ascending = False)
    
    num_rwrs = int(round(0.05 * len(rwr_sites)))
    rwr_sites = rwr_sites[:num_rwrs]
    rwr_sites_data = pd.merge(rwr_sites[['site']], data[['site', 'lat', 'long']], how='left', on='site')
    rwr_sites_data = rwr_sites_data.drop_duplicates()
    return rwr_sites_data

rwr_sites_survey = rwr_priority_sites(data_from_selected_sites, 'species')
rwr_sites_range = rwr_priority_sites(range_selected, 'sisid')

plot_sites_by_characteristic(rwr_sites_survey, 'lat', 'long', title='RWR Sites Survey Data')
plot_sites_by_characteristic(rwr_sites_range, 'lat', 'long', title='RWR Sites Range Data')



#bar plot of data type comparisons

#get hotspots for site data
hotspot_sites = plot_hotspots(data, 'richness', type = "hotspots")
rare_hotspots_sites = plot_hotspots(selected_rare, 'richness', type = "hotspots")
range_rich_hotspot = plot_hotspots(range_map, 'richness', species_col = 'sisid', type = "hotspots")
range_rare_hotspot = plot_hotspots(range_rare, 'richness', species_col = 'sisid', type = "hotspots")

#hotspots for non-exclude figures
survey_inc_hotspot = plot_hotspots(data_non_exclude, 'richness', type = "hotspots")
range_inc_hotspot = plot_hotspots(range_map_inc, 'richness', species_col = 'sisid', type = "hotspots")

site_rich_comp = (len(pd.merge(hotspot_sites, range_rich_hotspot, how='inner', on=['site', 'lat', 'long']))*2)/(len(hotspot_sites)+len(range_rich_hotspot))
site_rare_comp = (len(pd.merge(rare_hotspots_sites, range_rare_hotspot, how='inner', on=['site', 'lat', 'long']))*2)/(len(range_rare_hotspot)+len(rare_hotspots_sites))
cell_rich_comp = (len(pd.merge(obs_hotspot_cell, range_hotspot_cells, how='inner', on=['cellid', 'cent_lat', 'cent_long']))*2)/(len(obs_hotspot_cell)+len(range_hotspot_cells))
cell_rare_comp = (len(pd.merge(rare_survey_hotspot_cells, rare_range_hotspot_cells, how='inner', on=['cellid', 'cent_lat', 'cent_long']))*2)/(len(rare_survey_hotspot_cells)+len(rare_range_hotspot_cells))
rwr_comp = (len(pd.merge(rwr_sites_survey, rwr_sites_range, how='inner', on=['site', 'lat', 'long']))*2)/(len(rwr_sites_survey)+len(rwr_sites_range))
survey_ex_comp = (len(pd.merge(hotspot_sites, survey_inc_hotspot, how='inner', on=['site', 'lat', 'long']))*2)/(len(hotspot_sites)+len(survey_inc_hotspot))
range_ex_comp = (len(pd.merge(range_rich_hotspot, range_inc_hotspot, how='inner', on=['site', 'lat', 'long']))*2)/(len(range_rich_hotspot)+len(range_inc_hotspot))

perc = [site_rich_comp, cell_rich_comp, site_rare_comp, cell_rare_comp]
N = len(perc)
ind = np.arange(N)
width = 0.35

fig, ax = plt.subplots()
rects1 = ax.bar(ind, perc, width, color='brown')


ax.bar(ind, perc, width, color='maroon')
ax.set_ylim([0,1])
ax.set_ylabel('Hotspot Similarity Percentage')
tick_labels = ['site richness', 'cell richness', 'site rarity', 'cell rarity']
ax.set_xticks(ind+width)
ax.set_xticklabels(tick_labels)


plt.bar(ind, perc, color='maroon')
plt.xticks(ind, tick_labels)
plt.savefig('figures/comparison_barplot.png')

def bbs_exclude_families(species_list, aou_column, return_cols = None):
    # pass a dataframe with species ID (AOU) column and return that dataframe (or a subset of it)
    # with observations of bad BBS species excluded

    included_species = species_list[(species_list[aou_column] > 2880)]
    included_species = included_species[(included_species[aou_column] < 3650) | (included_species[aou_column] > 3810)]
    included_species = included_species[(included_species[aou_column] < 3900) | (included_species[aou_column] > 3910)]
    included_species = included_species[(included_species[aou_column] < 4160) | (included_species[aou_column] > 4210)]
    included_species = included_species[(included_species[aou_column] != 7010)]
    if aou_column == "AOU":
        included_species.rename(columns= {aou_column:'species'}, inplace = True)
    if return_cols:
        print("some columns")
        return included_species[return_cols]
    else:
        print("all columns")
        return included_species

analysis_species = bbs_exclude_families(included_species, "AOU", ['species', 'sisid'])
