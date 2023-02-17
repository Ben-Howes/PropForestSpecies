## Calculate total forested area at each degree latitude on the globe

import ee
import numpy
import geehydro
import folium
import pandas as pd
import os 
from multiprocessing import Pool

## Initialise earth engine
ee.Initialize()

image = ee.Image("UMD/hansen/global_forest_change_2021_v1_9")
esa = ee.ImageCollection("ESA/WorldCover/v100").first()

region = ee.Geometry.BBox(-180, -90, 180, 90)
image_clipped = image.clip(region)

forest = image.select('treecover2000').gte(70)

task_config = {
    'scale': 10000
}

task = ee.batch.Export.image(forest, 'forest_map', task_config)

task.start()

land = image.select('datamask').eq(1)

task_config = {
    'scale': 10000
}

task = ee.batch.Export.image(land, 'land_map', task_config)

task.start()

max_latitudes = range(-89, 91)
min_latitudes = range(-90, 90)

# ## We are also interested in breaking the world into 3 general longitudinal bands
# ## to see if any patterns hold across general regions
# ## These are as follows:
#     ## Americas: -180 to -30
#     ## Europe, Africa and ME: -30 to 60
#     ## Asia and Australia: 60 to 180

# max_longitudes = range(-179, 181)
# min_longitudes = range(-180, 180)

# final = pd.DataFrame()

# for i in range(len(max_latitudes)):
#     for j in range(len(max_longitudes)):
        
#         print(i)

#         lat_area = ee.Geometry.BBox(min_longitudes[j], min_latitudes[i], max_longitudes[j], max_latitudes[i])

#         clipped = image.clip(lat_area)
#         esa_clipped = esa.clip(lat_area)

#         area = ee.Image.pixelArea().updateMask(esa.select('Map').neq(80))

#         total_area = area.reduceRegion(**{
#                         'reducer':ee.Reducer.sum(),
#                         'geometry':esa_clipped.geometry(),
#                         'scale':1000,
#                         })

#         for k in [30, 50, 70]:
            
#             df = pd.DataFrame({'latitude_max': max_latitudes[i], 'latitude_min':min_latitudes[i], 
#                                'longitude_max': max_longitudes[j], 'longitude_min' : min_longitudes[j], 'forest_area':0, 
#                     'land_area':0, 'forest_def':k}, index = [i])
            
#             forest = ee.Image.pixelArea().updateMask(image.select('treecover2000').gte(k))

#             forest_area = forest.reduceRegion(**{
#                             'reducer':ee.Reducer.sum(),
#                             'geometry':clipped.geometry(),
#                             'scale':1000,
#                             })

#             df.forest_area[df.index[df['latitude_max'] == max_latitudes[i]]] = forest_area.getInfo()['area']
#             df.land_area[df.index[df['latitude_max'] == max_latitudes[i]]] = total_area.getInfo()['area']
            
#             final = pd.concat([final,df])
        
# final.to_csv("Documents/PhD/matrix_Response/Data/Ranges/forest_by_latitude.csv")

# ## Download land area map

# land = esa.select('Map').neq(80)
# region = ee.Geometry.BBox(-180, -90, 180, 90)
# land_clipped = land.clip(region)

# task_config = {
#     'scale': 5000,  
#     # 'region': land_clipped.geometry()
#     }

# task = ee.batch.Export.image(land_clipped, 'land_map', task_config)

# task.start()

# def forest_area(i):
    
#     final = pd.DataFrame()
#     print(i)
           
#     for j in range(len(max_longitudes)):

#         lat_area = ee.Geometry.BBox(min_longitudes[j], min_latitudes[i], max_longitudes[j], max_latitudes[i])

#         clipped = image.clip(lat_area)
#         esa_clipped = esa.clip(lat_area)

#         esa_area = ee.Image.pixelArea().updateMask(esa.select('Map').neq(80))

#         esa_total_area = esa_area.reduceRegion(**{
#                         'reducer':ee.Reducer.sum(),
#                         'geometry':esa_clipped.geometry(),
#                         'scale':1000,
#                         })
        
#         hansen_area = ee.Image.pixelArea().updateMask(image.select('datamask').eq(1))

#         hansen_total_area = hansen_area.reduceRegion(**{
#                 'reducer':ee.Reducer.sum(),
#                 'geometry':clipped.geometry(),
#                 'scale':1000,
#                 })

#         for k in [30, 50, 70]:
            
#             df = pd.DataFrame({'latitude_max': max_latitudes[i], 'latitude_min':min_latitudes[i], 
#                                'longitude_max': max_longitudes[j], 'longitude_min' : min_longitudes[j], 'forest_area':0, 
#                     'esa_land_area':0, 'hansen_land_area':0, 'forest_def':k}, index = [i])
            
#             forest = ee.Image.pixelArea().updateMask(image.select('treecover2000').gte(k))

#             forest_area = forest.reduceRegion(**{
#                             'reducer':ee.Reducer.sum(),
#                             'geometry':clipped.geometry(),
#                             'scale':1000,
#                             })

#             df.forest_area[df.index[df['latitude_max'] == max_latitudes[i]]] = forest_area.getInfo()['area']
#             df.esa_land_area[df.index[df['latitude_max'] == max_latitudes[i]]] = esa_total_area.getInfo()['area']
#             df.hansen_land_area[df.index[df['latitude_max'] == max_latitudes[i]]] = hansen_total_area.getInfo()['area']
            
#             final = pd.concat([final,df])
            
#     return(final)
        
# forest_area(120)
# map(test_func,[66,53,0,94])

# with Pool(processes=12) as P:
#     out = P.map(forest_area, range(len(max_latitudes)))
    
# final = pd.concat(out)

# final.to_csv("forest_by_latlong1.csv")