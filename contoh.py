import os
from os import path
import glob
from pykml import parser
from shapely.geometry import Point, Polygon
import ast
import pandas as pd
from datetime import datetime
import ctypes
import sys
import geopandas as gpd
import numpy as np

def Random_Points_in_Bounds(polygon, number):   
    minx, miny, maxx, maxy = polygon.bounds
    x = np.round(np.random.uniform(minx, maxx, number),10)
    y = np.round(np.random.uniform(miny, maxy, number),10) 
    return x, y

# read coordinates points file
df = pd.read_excel("polygons_and_coordinates/coordinates.xlsx")
df.columns = df.columns.str.lower()

# sort directory and quantity
lstdir = (glob.glob(os.path.join('polygons_and_coordinates', '*.kml')))
lstdirsorted = sorted(lstdir)
lendir = len(lstdirsorted)
print(str(lendir) + " .kml files detected")

for filepath in lstdirsorted:
    
    # open .kml
    with open(filepath) as f:
        
        root = parser.parse(f).getroot()
        
        print("Running... " + root.Document.name.text)
        
        #get coordinates tag
        strpolygon = root.Document.Placemark.MultiGeometry.LineString.coordinates.text
        
        # cleaning
        strpolygon = strpolygon.strip()
        strpolygon = "[("+ strpolygon + ")]"
        strpolygon = strpolygon.replace(" ","),(")
        strpolygon = strpolygon.replace(" ","),(")
        strpolygon = strpolygon.replace(",0)",")")
        
        #input
        lstpolygon = ast.literal_eval(strpolygon)
        polygon = Polygon(lstpolygon)
        
        gdf_poly = gpd.GeoDataFrame(index=["myPoly"], geometry=[polygon])
        
        x,y = Random_Points_in_Bounds(polygon, 100000) #number of desired points
        
        df['latitude'] = pd.Series(y.tolist())
        df['longitude'] = pd.Series(x.tolist())

        print(root.Document.name.text + " successfully processed")
        
#save
df.to_excel (r"polygons_and_coordinates/coordinates.xlsx", index = False, header=True)
print("Finished!")

MessageBox = ctypes.windll.user32.MessageBoxW
MessageBox(None, 'Finished', 'Finished!', 0)