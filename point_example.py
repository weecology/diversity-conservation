import ogr
inDataSource = ogr.Open("bbs_sites_coordinates_wrapper.vrt")
lyr = inDataSource.GetLayer("bbs_sites")
for feat in lyr:
    geom = feat.GetGeometryRef()
    print geom.ExportToWkt()

