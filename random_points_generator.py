import os
import glob
import pandas as pd
import numpy as np
from shapely.geometry import Polygon, Point
from pykml import parser
import geopandas as gpd
import ctypes

def random_points_in_polygon(polygon, number):
    minx, miny, maxx, maxy = polygon.bounds
    x, y = [], []
    while len(x) < number:
        px = np.random.uniform(minx, maxx)
        py = np.random.uniform(miny, maxy)
        point = Point(px, py)
        if polygon.contains(point):
            x.append(round(px, 8))
            y.append(round(py, 8))
    return np.array(x), np.array(y)

input_folder = "polygons_and_coordinates"
output_file = os.path.join(input_folder, "coordinates.xlsx")

if not os.path.exists(input_folder):
    os.makedirs(input_folder)

kml_files = sorted(glob.glob(os.path.join(input_folder, "*.kml")))
if not kml_files:
    print("❌ Tidak ada file .kml di folder:", input_folder)
    exit()

print(f"✅ {len(kml_files)} file KML terdeteksi")

all_points = []

for filepath in kml_files:
    with open(filepath, "r", encoding="utf-8") as f:
        root = parser.parse(f).getroot()

        name = root.Document.name.text if hasattr(root.Document, "name") else os.path.basename(filepath)
        print(f"🔹 Memproses {name} ...")

        coords_text = root.Document.Placemark.MultiGeometry.LineString.coordinates.text.strip()

        coords_list = []
        for c in coords_text.split():
            lon, lat, *_ = map(float, c.split(","))
            coords_list.append((lon, lat))

        polygon = Polygon(coords_list)
        if not polygon.is_valid:
            print(f"⚠️ Polygon tidak valid di {name}, dilewati.")
            continue

        x, y = random_points_in_polygon(polygon, number=1000)

        df_temp = pd.DataFrame({
            "Polygon": name,
            "Longitude": x,
            "Latitude": y
        })
        all_points.append(df_temp)

if all_points:
    df_all = pd.concat(all_points, ignore_index=True)
    df_all.to_excel(output_file, index=False)
    print(f"✅ Titik acak berhasil disimpan di: {output_file}")
else:
    print("⚠️ Tidak ada titik yang dihasilkan.")

ctypes.windll.user32.MessageBoxW(0, "Proses selesai!", "Random Point Generator", 0)
