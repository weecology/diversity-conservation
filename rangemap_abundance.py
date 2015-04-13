import numpy as np
from osgeo import ogr
import os
import pandas as pd

os.chdir('/Users/karinorman/Documents/reserve_selection/data/BOTW')
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

#example intersection    
wkt1 = "POLYGON ((1208064.271243039 624154.6783778917, 1208064.271243039 601260.9785661874, 1231345.9998651114 601260.9785661874, 1231345.9998651114 624154.6783778917, 1208064.271243039 624154.6783778917))"
wkt2 = "POLYGON ((1199915.6662253144 633079.3410163528, 1199915.6662253144 614453.958118695, 1219317.1067437078 614453.958118695, 1219317.1067437078 633079.3410163528, 1199915.6662253144 633079.3410163528)))"

poly1 = ogr.CreateGeometryFromWkt(wkt1)
poly2 = ogr.CreateGeometryFromWkt(wkt2)

intersection = poly1.Intersection(poly2)

print intersection.ExportToWkt()    

#intersection of two species
first_species = all_sp_layer.GetFeature(1)
second_species = all_sp_layer.GetFeature(4)

first_species_geom = first_species.GetGeometryRef()
second_species_geom = second_species.GetGeometryRef()

inter = first_species_geom.Intersection(second_species_geom)
print inter

for i in range(all_sp_layer.GetFeatureCount()):
    second_species = all_sp_layer.GetFeature(i+1)
    second_species_geom = second_species.GetGeometryRef()
    inter = first_species_geom.Intersection(second_species_geom)
    print inter
    
# Filter features to only include those where Seasonal == 1 or 2 (i.e., where the polygon is either resident or breeding season)
all_sp_layer.SetAttributeFilter("Seasonal = 1 or Seasonal = 2")

# Get species names for every polygon
species = []
for i in range(all_sp_layer.GetFeatureCount()):
    feature = all_sp_layer.GetFeature(i + 1)
    sciname = feature.GetField('SCINAME')
    species.append(sciname)

# unique list of site coordinates
# loop over polygons for each point 
# list of species names for each site

inDataSource = ogr.Open("bbs_sites_coordinates_wrapper.vrt")
site_lyr = inDataSource.GetLayer('bbs_sites')
for feat in site_lyr:
    geom = feat.GetGeometryRef()
    print geom.ExportToWkt() 


for i in range(site_lyr.GetFeatureCount()):
    indv_site = site_lyr.GetFeature(i+1)
    indv_site_geom = indv_site.GetGeometryRef()
    inter = first_species_geom.Intersection(indv_site_geom)
    print inter