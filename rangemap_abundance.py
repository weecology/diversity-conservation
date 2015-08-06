import numpy as np
from osgeo import ogr
import os
import pandas as pd

os.chdir('.../Documents/reserve_selection/data/BOTW')
driver = ogr.GetDriverByName('OpenFileGDB')
gdb = driver.Open("BOTW.gdb", 0)
    
# NOTES:
# There is one layer, which includes polygons for all species
# We can get a layer using
all_sp_layer = gdb.GetLayer("All_Spp")
# There are ~17000 features, which are one or more polygons for each species
# We can get a feature by using
first_feature = all_sp_layer.GetFeature(1) #indexing appears to start at 1
# Each feature (which is a polygon) has 16 "items" or "attributes"
# We can get a list of features by
first_feature.items()
# We can get a particular value by
first_feature.GetField('SCINAME')

# Filter features to only include those where Seasonal == 1 or 2 (i.e., where the polygon is either resident or breeding season)
all_sp_layer.SetAttributeFilter("Seasonal = 1 or Seasonal = 2")

# Get species names for every polygon
#species = []
#for i in range(all_sp_layer.GetFeatureCount()):
    #feature = all_sp_layer.GetFeature(i + 1)
    #sciname = feature.GetField('SCINAME')
    #species.append(sciname)

#import sites
inDataSource = ogr.Open("bbs_sites_coordinates_wrapper.vrt")
site_lyr = inDataSource.GetLayer('bbs_sites')
for feat in site_lyr:
    geom = feat.GetGeometryRef()
    print geom.ExportToWkt() 
 
species_count = pd.DataFrame()
for i in range(site_lyr.GetFeatureCount()):
    print "i %d" %i
    num_species = []
    indv_site = site_lyr.GetFeature(i+1)
    indv_site_geom = indv_site.GetGeometryRef()
    for j in range(all_sp_layer.GetFeatureCount()):
        print "j %d" %j
        indv_species = all_sp_layer.GetFeature(j+1)
        indv_species_geom = indv_species.GetGeometryRef()
        inter = indv_site_geom.Intersection(indv_species_geom)
        if inter != "GEOMETRYCOLLECTION EMPTY":
            num_species.append(1)
    species_count.append([indv_site.GetField('site'), len(num_species)])
