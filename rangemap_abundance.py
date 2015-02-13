import numpy as np
from osgeo import ogr
import os

os.chdir(/Users/karinorman/Documents/reserve_selection/data/BOTW)
driver = ogr.GetDriverByName('OpenFileGDB')
gdb = driver.Open("BOTW.gdb", 0)

# list to store layers'names
featsClassList = []

# parsing layers by index
for featsClass_idx in range(gdb.GetLayerCount()):
    featsClass = gdb.GetLayerByIndex(featsClass_idx)
    featsClassList.append(featsClass.GetName())

# sorting
featsClassList.sort()

# printing
for featsClass in featsClassList:
    print featsClass