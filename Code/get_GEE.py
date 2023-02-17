## Sample % cover of land types around all sites in each study at 200m and 2000m scales
## Data taken from the ESA 10m land cover dataset

import ee
import numpy
import geehydro
import folium
import pandas as pd
import os 

pd.options.mode.chained_assignment = None  # default='warn'

## Initialise earth engine
ee.Initialize()

## Read in csv with coordiantes for every site in every study
sites = pd.read_csv("Documents/PhD/matrix_Response/Data/all_sites.csv", encoding = "ISO-8859-1")[['pid', 'plot', 'latitude', 'longitude']]
sites = sites.rename(columns = {"plot" : "site"})

final = pd.DataFrame()

for i in range(len(sites)):
    
    print(i)
    
    ## Define site coordinates
    site = ee.Geometry.Point([sites.longitude[i], sites.latitude[i]])

    ## Load in ESA map within bounds of site and filter to first image: 2020
    image = ee.Image(ee.ImageCollection("ESA/WorldCover/v100")
                    .filterBounds(site)
                    .first())

    ## Set buffer size for measuring classes within
    ## Then buffer and clip image to just get that area
    buffer = 2000
    buff = site.buffer(buffer)
    clipped = image.clip(buff)

    ## Make data frame of classes and their codes to be filled in later with area
    codes = [10,20,30,40,50,60,70,80,90,95,100]
    landuses = ('trees', 'shrub', 'grass', 'crops', 'urban', 'barren', 'snow',
            'water', 'wetland', 'mangrove','moss')
    df = pd.DataFrame({'pid': sites.pid[i], 'plot':sites.site[i], 'landuse':landuses, 'code': codes, 'area':0, 'buffer':buffer})

    for code in df.code:
        
        ## Select image and mask with just the pixels that are of teh class (code) we're interested in
        area = ee.Image.pixelArea().updateMask(image.select('Map').eq(code))

        ## Sum all pixels of that class using the area, this means we don't need to worry about crs
        out = area.reduceRegion(**{
            'reducer':ee.Reducer.sum(),
            'geometry':clipped.geometry(),
            'scale':10,
            })

        ## Put area output into our df
        df.area[df.index[df['code'] == code]] = out.getInfo()['area']

    ## calculate the percentage of the area taken up by each class
    df['percentage'] = df['area']/sum(df.area)

    final = pd.concat([final,df])