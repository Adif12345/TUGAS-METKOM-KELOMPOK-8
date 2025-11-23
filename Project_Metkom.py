"""
Koordinat Point Generator — PRO / GIS MINI
Versi untuk VSCode (single-file app)

Fitur utama:
- GUI modern (ttk/ttkbootstrap bila tersedia)
- Tab: Generate | Transform | Map Preview | Export | Settings
- Banyak mode generate: square, rectangle, circle, ellipse, hex, line, polyline, ring, random uniform, random gaussian
- Input: start coord, offsets, total points or points per side
- Rotate without scaling (around origin)
- Convert LatLon <-> UTM (pyproj)
- Export: TXT, CSV, XLSX, GeoJSON, Shapefile (via geopandas)
- Map preview via geopandas/matplotlib
"""

import os
import math
import json
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter.scrolledtext import ScrolledText
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# GIS libs
try:
    import geopandas as gpd
    from shapely.geometry import Point, LineString, Polygon
    from pyproj import Transformer, CRS
    HAS_GIS = True
except Exception:
    HAS_GIS = False

# Try ttkbootstrap for nicer theme; fallback to ttk
try:
    import ttkbootstrap as tb
    from ttkbootstrap.constants import *
    THEME = "flatly"
    USE_BOOTSTRAP = True
except Exception:
    USE_BOOTSTRAP = False

# -------------------------
# Helper utilities
# -------------------------
def ensure_gis_available():
    if not HAS_GIS:
        raise RuntimeError("GIS libraries not found. Install geopandas, shapely, pyproj, fiona.")

def to_geodf(points_list, crs_epsg=4326):
    """points_list: iterable of (x,y). Returns GeoDataFrame of Points."""
    ensure_gis_available()
    geoms = [Point(x, y) for x, y in points_list]
    gdf = gpd.GeoDataFrame({"x":[p[0] for p in points_list], "y":[p[1] for p in points_list]}, geometry=geoms, crs=f"EPSG:{crs_epsg}")
    return gdf

def export_gdf_to_shapefile(gdf, path):
    ensure_gis_available()
    gdf.to_file(path, driver="ESRI Shapefile")

def export_gdf_to_geojson(gdf, path):
    ensure_gis_available()
    gdf.to_file(path, driver="GeoJSON")

def safe_float(v, default=0.0):
    try:
        return float(v)
    except Exception:
        return default

# -------------------------
# Geometry generation routines
# -------------------------
def generate_square_grid(n_side, spacing, x0, y0, total_mode=False):
    """Generate square grid (nxn) anchored at x0,y0 (lower-left)."""
    if total_mode:
        # n_side is total points: compute nearest square
        Ntot = int(n_side)
        n = int(np.ceil(np.sqrt(Ntot)))
    else:
        n = int(n_side)
    xs = np.arange(0, n * spacing, spacing)
    ys = np.arange(0, n * spacing, spacing)
    pts = [(x0 + x, y0 + y) for x in xs for y in ys]
    return pts

def generate_rect_grid(nx, ny, spacing_x, spacing_y, x0, y0):
    xs = np.arange(0, nx * spacing_x, spacing_x)
    ys = np.arange(0, ny * spacing_y, spacing_y)
    return [(x0 + x, y0 + y) for x in xs for y in ys]

def generate_circle(center_x, center_y, radius, n_points):
    thetas = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    return [(center_x + radius*np.cos(t), center_y + radius*np.sin(t)) for t in thetas]

def generate_ellipse(center_x, center_y, a, b, n_points, rotation_deg=0):
    thetas = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    rot = math.radians(rotation_deg)
    pts = []
    for t in thetas:
        x = a * math.cos(t)
        y = b * math.sin(t)
        xr = x*math.cos(rot) - y*math.sin(rot)
        yr = x*math.sin(rot) + y*math.cos(rot)
        pts.append((center_x + xr, center_y + yr))
    return pts

def generate_hex_grid(n_side, spacing, x0, y0):
    # axial hex grid building, roughly covering n_side x n_side
    pts = []
    for row in range(n_side):
        for col in range(n_side):
            x = col * spacing * 0.75
            y = row * spacing * (math.sqrt(3)/2)
            # offset every other column
            if col % 2 == 1:
                y += (math.sqrt(3)/4) * spacing
            pts.append((x0 + x, y0 + y))
    return pts

def generate_line(x0, y0, x1, y1, n_points):
    xs = np.linspace(x0, x1, n_points)
    ys = np.linspace(y0, y1, n_points)
    return list(zip(xs, ys))

def generate_polyline(coords_list, density_per_segment=10):
    # coords_list: [(x1,y1),(x2,y2),...]
    pts = []
    for i in range(len(coords_list)-1):
        (x0,y0),(x1,y1) = coords_list[i], coords_list[i+1]
        segment = generate_line(x0,y0,x1,y1,density_per_segment)
        pts.extend(segment[:-1])  # avoid double-count
    pts.append(coords_list[-1])
    return pts

def generate_ring(center_x, center_y, inner_r, outer_r, n_points):
    # uniformly sample annulus ring by angle + radius
    thetas = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    pts = []
    for t in thetas:
        r = (inner_r + outer_r) / 2.0
        pts.append((center_x + r*math.cos(t), center_y + r*math.sin(t)))
    return pts

def generate_random_uniform(xmin, xmax, ymin, ymax, n_points):
    xs = np.random.uniform(xmin, xmax, n_points)
    ys = np.random.uniform(ymin, ymax, n_points)
    return list(zip(xs, ys))

def generate_random_gaussian(cx, cy, sigma_x, sigma_y, n_points):
    xs = np.random.normal(cx, sigma_x, n_points)
    ys = np.random.normal(cy, sigma_y, n_points)
    return list(zip(xs, ys))

# -------------------------
# Coordinate transforms (LatLon <-> UTM)
# -------------------------
def lonlat_to_utm(lons, lats, zone, hemisphere='south'):
    ensure_gis_available()
    # select EPSG code
    epsg = 32700 + int(zone) if hemisphere == 'south' else 32600 + int(zone)
    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    ux, uy = transformer.transform(lons, lats)
    return ux, uy

def utm_to_lonlat(xs, ys, zone, hemisphere='south'):
    ensure_gis_available()
    epsg = 32700 + int(zone) if hemisphere == 'south' else 32600 + int(zone)
    transformer = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
    lon, lat = transformer.transform(xs, ys)
    return lon, lat

# -------------------------
# Main Application GUI
# -------------------------
class PROGISApp:
    def __init__(self, master):
        if USE_BOOTSTRAP:
            self.root = tb.Window(themename=THEME, title="Koordinat Point Generator - PRO GIS")
            self.master = self.root
        else:
            self.root = master
            self.master = master
            master.title("Koordinat Point Generator - PRO GIS")
            master.geometry("1200x820")

        self.points = []            # current points in map coordinates (units: meter)
        self.current_crs_epsg = 32748  # default example (UTM zone 48S) — EPSG code if needed
        self.transformed = []       # transformed points after ops

        self._build_ui()

    def _build_ui(self):
        # Top frame for tabs
        self.tab_control = ttk.Notebook(self.master)
        # Tabs
        self.tab_generate = ttk.Frame(self.tab_control)
        self.tab_transform = ttk.Frame(self.tab_control)
        self.tab_map = ttk.Frame(self.tab_control)
        self.tab_export = ttk.Frame(self.tab_control)
        self.tab_settings = ttk.Frame(self.tab_control)

        self.tab_control.add(self.tab_generate, text="Generate")
        self.tab_control.add(self.tab_transform, text="Transform")
        self.tab_control.add(self.tab_map, text="Map Preview")
        self.tab_control.add(self.tab_export, text="Export")
        self.tab_control.add(self.tab_settings, text="Settings")
        self.tab_control.pack(expand=1, fill="both")

        # ----------------- GENERATE TAB -----------------
        g = self.tab_generate
        lf = ttk.LabelFrame(g, text="Options", padding=8)
        lf.pack(side="left", fill="y", padx=8, pady=8)

        ttk.Label(lf, text="Mode:").grid(row=0, column=0, sticky="w")
        self.mode_var = tk.StringVar(value="square")
        ttk.Combobox(lf, textvariable=self.mode_var, values=[
            "square","rectangle","circle","ellipse","hex","line","polyline","ring","random_uniform","random_gaussian"
        ], width=18).grid(row=0, column=1, padx=4, pady=2)

        # Start coordinate
        ttk.Label(lf, text="Start X:").grid(row=1, column=0, sticky="w")
        self.start_x = tk.DoubleVar(value=500000.0)
        ttk.Entry(lf, textvariable=self.start_x, width=12).grid(row=1, column=1)
        ttk.Label(lf, text="Start Y:").grid(row=2, column=0, sticky="w")
        self.start_y = tk.DoubleVar(value=908000.0)
        ttk.Entry(lf, textvariable=self.start_y, width=12).grid(row=2, column=1)

        # Grid size and points
        ttk.Label(lf, text="Spacing (m):").grid(row=3, column=0, sticky="w")
        self.spacing = tk.DoubleVar(value=10.0)
        ttk.Entry(lf, textvariable=self.spacing, width=12).grid(row=3, column=1)

        ttk.Label(lf, text="Points per side / Total:").grid(row=4, column=0, sticky="w")
        self.count_side = tk.IntVar(value=10)
        ttk.Entry(lf, textvariable=self.count_side, width=12).grid(row=4, column=1)

        # Rectangle specific
        ttk.Label(lf, text="Rect NX:").grid(row=5, column=0, sticky="w")
        self.rect_nx = tk.IntVar(value=10)
        ttk.Entry(lf, textvariable=self.rect_nx, width=8).grid(row=5, column=1, sticky="w")
        ttk.Label(lf, text="Rect NY:").grid(row=6, column=0, sticky="w")
        self.rect_ny = tk.IntVar(value=6)
        ttk.Entry(lf, textvariable=self.rect_ny, width=8).grid(row=6, column=1, sticky="w")

        # Circle/Ellipse specific
        ttk.Label(lf, text="Circle radius (m):").grid(row=7, column=0, sticky="w")
        self.circle_r = tk.DoubleVar(value=50.0)
        ttk.Entry(lf, textvariable=self.circle_r, width=12).grid(row=7, column=1)
        ttk.Label(lf, text="Ellipse a,b (m):").grid(row=8, column=0, sticky="w")
        self.ellipse_a = tk.DoubleVar(value=60.0)
        self.ellipse_b = tk.DoubleVar(value=30.0)
        ttk.Entry(lf, textvariable=self.ellipse_a, width=6).grid(row=8, column=1, sticky="w")
        ttk.Entry(lf, textvariable=self.ellipse_b, width=6).grid(row=8, column=1, sticky="e")

        # Line & polyline
        ttk.Label(lf, text="Line x1,y1,x2,y2:").grid(row=9, column=0, sticky="w")
        self.line_x1 = tk.DoubleVar(value=500000.0); self.line_y1 = tk.DoubleVar(value=908000.0)
        self.line_x2 = tk.DoubleVar(value=500200.0); self.line_y2 = tk.DoubleVar(value=908200.0)
        ttk.Entry(lf, textvariable=self.line_x1, width=6).grid(row=9, column=1, sticky="w")
        ttk.Entry(lf, textvariable=self.line_y1, width=6).grid(row=9, column=1)
        ttk.Entry(lf, textvariable=self.line_x2, width=6).grid(row=10, column=1, sticky="w")
        ttk.Entry(lf, textvariable=self.line_y2, width=6).grid(row=10, column=1)

        # Random gaussian params
        ttk.Label(lf, text="Gaussian sigma X,Y:").grid(row=11, column=0, sticky="w")
        self.gauss_sx = tk.DoubleVar(value=30.0)
        self.gauss_sy = tk.DoubleVar(value=30.0)
        ttk.Entry(lf, textvariable=self.gauss_sx, width=6).grid(row=11, column=1, sticky="w")
        ttk.Entry(lf, textvariable=self.gauss_sy, width=6).grid(row=11, column=1)

        # Hex params
        ttk.Label(lf, text="Hex side n:").grid(row=12, column=0, sticky="w")
        self.hex_n = tk.IntVar(value=8)
        ttk.Entry(lf, textvariable=self.hex_n, width=6).grid(row=12, column=1, sticky="w")

        # Controls
        ttk.Button(lf, text="Generate", command=self.on_generate).grid(row=13, column=0, pady=6)
        ttk.Button(lf, text="Clear", command=self.on_clear).grid(row=13, column=1, pady=6)

        # Right side: plotting area & table
        right_frame = ttk.Frame(self.tab_generate)
        right_frame.pack(side="left", fill="both", expand=True, padx=8, pady=8)

        self.fig, self.ax = plt.subplots(figsize=(7,6))
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        self.canvas_plot = FigureCanvasTkAgg(self.fig, master=right_frame)
        self.canvas_plot.get_tk_widget().pack(fill="both", expand=True)

        self.table = ttk.Treeview(right_frame, columns=("x","y"), show="headings", height=8)
        self.table.heading("x", text="X")
        self.table.heading("y", text="Y")
        self.table.pack(fill="x")

        # ----------------- TRANSFORM TAB -----------------
        t = self.tab_transform
        tf = ttk.LabelFrame(t, text="Transform & Projection", padding=8)
        tf.pack(side="left", fill="y", padx=8, pady=8)

        # Rotate center
        ttk.Label(tf, text="Rotate center X:").grid(row=0, column=0)
        self.rot_cx = tk.DoubleVar(value=500000.0)
        ttk.Entry(tf, textvariable=self.rot_cx, width=12).grid(row=0, column=1)
        ttk.Label(tf, text="Rotate center Y:").grid(row=1, column=0)
        self.rot_cy = tk.DoubleVar(value=908000.0)
        ttk.Entry(tf, textvariable=self.rot_cy, width=12).grid(row=1, column=1)
        ttk.Label(tf, text="Angle (deg):").grid(row=2, column=0)
        self.rot_angle = tk.DoubleVar(value=0.0)
        ttk.Entry(tf, textvariable=self.rot_angle, width=12).grid(row=2, column=1)
        ttk.Button(tf, text="Rotate", command=self.on_rotate).grid(row=3, column=0, columnspan=2, pady=6)

        ttk.Separator(tf, orient="horizontal").grid(row=4, column=0, columnspan=2, sticky="we", pady=6)
        ttk.Label(tf, text="UTM Zone:").grid(row=5, column=0)
        self.utm_zone = tk.IntVar(value=48)
        ttk.Entry(tf, textvariable=self.utm_zone, width=8).grid(row=5, column=1)
        ttk.Label(tf, text="Hemisphere:").grid(row=6, column=0)
        self.utm_hemi = tk.StringVar(value="south")
        ttk.Combobox(tf, textvariable=self.utm_hemi, values=["north","south"], width=8).grid(row=6, column=1)
        ttk.Button(tf, text="To LonLat (UTM->LL)", command=self.on_utm_to_ll).grid(row=7, column=0, columnspan=2, pady=6)
        ttk.Button(tf, text="To UTM (LL->UTM)", command=self.on_ll_to_utm).grid(row=8, column=0, columnspan=2, pady=6)

        # ----------------- MAP TAB -----------------
        m = self.tab_map
        mf = ttk.LabelFrame(m, text="Map Preview & Overlays", padding=8)
        mf.pack(fill="both", expand=True, padx=8, pady=8)
        # small plot reuse
        self.fig_map, self.ax_map = plt.subplots(figsize=(8,6))
        self.canvas_map = FigureCanvasTkAgg(self.fig_map, master=mf)
        self.canvas_map.get_tk_widget().pack(fill="both", expand=True)
        # overlay controls
        ovf = ttk.Frame(mf)
        ovf.pack(fill="x")
        ttk.Button(ovf, text="Plot Points", command=self.plot_map_points).pack(side="left", padx=4)
        ttk.Button(ovf, text="Load Shapefile / GeoJSON", command=self.load_vector).pack(side="left", padx=4)

        # ----------------- EXPORT TAB -----------------
        e = self.tab_export
        ef = ttk.LabelFrame(e, text="Export Options", padding=8)
        ef.pack(fill="both", expand=True, padx=8, pady=8)

        ttk.Button(ef, text="Export TXT", command=self.export_txt).pack(anchor="w")
        ttk.Button(ef, text="Export CSV", command=self.export_csv).pack(anchor="w")
        ttk.Button(ef, text="Export XLSX", command=self.export_xlsx).pack(anchor="w")
        ttk.Button(ef, text="Export GeoJSON", command=self.export_geojson).pack(anchor="w")
        ttk.Button(ef, text="Export Shapefile (.shp)", command=self.export_shapefile).pack(anchor="w")

        # ----------------- SETTINGS TAB -----------------
        s = self.tab_settings
        sf = ttk.LabelFrame(s, text="Settings & Info", padding=8)
        sf.pack(fill="both", expand=True, padx=8, pady=8)
        txt = ScrolledText(sf, height=20)
        txt.insert("1.0", "Koordinat Point Generator - PRO GIS MINI\\n\\nRequirements:\\n- Python 3.8+\\n- geopandas, shapely, pyproj, fiona, pandas, numpy, matplotlib\\n\\nWorkflow tips:\\n- Generate grid or other shapes in GENERATE tab\\- Transform / rotate in TRANSFORM tab\\- Preview on MAP tab\\- Export from EXPORT tab\\")
        txt.configure(state="disabled")
        txt.pack(fill="both", expand=True)

    # -------------------------
    # Event handlers
    # -------------------------
    def on_generate(self):
        mode = self.mode_var.get()
        x0 = self.start_x.get()
        y0 = self.start_y.get()
        spacing = self.spacing.get()
        cnt = self.count_side.get()

        pts = []
        try:
            if mode == "square":
                pts = generate_square_grid(cnt, spacing, x0, y0, total_mode=False)
            elif mode == "rectangle":
                nx = self.rect_nx.get(); ny = self.rect_ny.get()
                pts = generate_rect_grid(nx, ny, spacing, spacing, x0, y0)
            elif mode == "circle":
                pts = generate_circle(x0, y0, self.circle_r.get(), int(cnt))
            elif mode == "ellipse":
                pts = generate_ellipse(x0, y0, self.ellipse_a.get(), self.ellipse_b.get(), int(cnt))
            elif mode == "hex":
                pts = generate_hex_grid(self.hex_n.get(), spacing, x0, y0)
            elif mode == "line":
                pts = generate_line(self.line_x1.get(), self.line_y1.get(), self.line_x2.get(), self.line_y2.get(), int(cnt))
            elif mode == "polyline":
                # for demo, create simple polyline from several points set relative to x0,y0
                coords = [(x0, y0), (x0+spacing*5, y0+spacing*2), (x0+spacing*8, y0-spacing*3)]
                pts = generate_polyline(coords, density_per_segment=int(cnt))
            elif mode == "ring":
                pts = generate_ring(x0, y0, self.circle_r.get()*0.6, self.circle_r.get()*1.2, int(cnt))
            elif mode == "random_uniform":
                pts = generate_random_uniform(x0, x0+spacing*cnt, y0, y0+spacing*cnt, int(cnt))
            elif mode == "random_gaussian":
                pts = generate_random_gaussian(x0 + spacing*cnt/2.0, y0 + spacing*cnt/2.0, self.gauss_sx.get(), self.gauss_sy.get(), int(cnt))
            else:
                pts = []
        except Exception as e:
            messagebox.showerror("Error Generate", str(e))
            return

        self.points = pts
        self.transformed = []
        self._update_table(self.points)
        self._plot_points(self.points, title=f"Generated: {mode}")

    def on_clear(self):
        self.points = []
        self.transformed = []
        self._update_table([])
        self._plot_points([], title="Cleared")

    def _update_table(self, pts):
        # clear
        for item in self.table.get_children():
            self.table.delete(item)
        for x,y in pts:
            self.table.insert("", "end", values=(f"{x:.3f}", f"{y:.3f}"))

    def _plot_points(self, pts, title="Plot"):
        self.ax.clear()
        if pts:
            xs, ys = zip(*pts)
            self.ax.scatter(xs, ys, s=10)
        self.ax.set_title(title)
        self.ax.grid(True)
        self.canvas_plot.draw()

    def on_rotate(self):
        if not self.points:
            messagebox.showwarning("No Points", "Generate points first.")
            return
        cx = self.rot_cx.get(); cy = self.rot_cy.get()
        angle = math.radians(self.rot_angle.get())
        c, s = math.cos(angle), math.sin(angle)
        rotated = []
        for x,y in self.points:
            dx = x - cx; dy = y - cy
            xr = dx*c - dy*s + cx
            yr = dx*s + dy*c + cy
            rotated.append((xr, yr))
        self.transformed = rotated
        self._update_table(self.transformed)
        self._plot_points(self.transformed, title="Rotated Points")

    def on_ll_to_utm(self):
        if not HAS_GIS:
            messagebox.showerror("Missing libs", "Geopandas/pyproj required.")
            return
        # expects points in lon, lat currently (user must place lon/lat in start_x/start_y)
        lons = [p[0] for p in self.points]
        lats = [p[1] for p in self.points]
        zone = self.utm_zone.get(); hemi = self.utm_hemi.get()
        try:
            ux, uy = lonlat_to_utm(lons, lats, zone, hemi)
            self.transformed = list(zip(ux, uy))
            self._update_table(self.transformed)
            self._plot_points(self.transformed, title="LL -> UTM")
        except Exception as e:
            messagebox.showerror("Transform Error", str(e))

    def on_utm_to_ll(self):
        if not HAS_GIS:
            messagebox.showerror("Missing libs", "Geopandas/pyproj required.")
            return
        xs = [p[0] for p in self.points]
        ys = [p[1] for p in self.points]
        zone = self.utm_zone.get(); hemi = self.utm_hemi.get()
        try:
            lon, lat = utm_to_lonlat(xs, ys, zone, hemi)
            self.transformed = list(zip(lon, lat))
            self._update_table(self.transformed)
            self._plot_points(self.transformed, title="UTM -> LL")
        except Exception as e:
            messagebox.showerror("Transform Error", str(e))

    def plot_map_points(self):
        if not HAS_GIS:
            messagebox.showerror("Missing libs", "Geopandas & shapely required for map preview.")
            return
        self.ax_map.clear()
        # plot base if any
        try:
            if hasattr(self, "last_vector_gdf") and self.last_vector_gdf is not None:
                self.last_vector_gdf.plot(ax=self.ax_map, facecolor="none", edgecolor="gray")
        except Exception:
            pass
        if self.points:
            gdf = to_geodf(self.points, crs_epsg=32748)  # temporary EPSG assume UTM48S for plotting ease
            gdf.plot(ax=self.ax_map, markersize=8, color="red")
        if self.transformed:
            gdf2 = to_geodf(self.transformed, crs_epsg=32748)
            gdf2.plot(ax=self.ax_map, markersize=8, color="blue")
        self.ax_map.set_title("Map Preview")
        self.canvas_map.draw()

    def load_vector(self):
        if not HAS_GIS:
            messagebox.showerror("Missing libs", "Geopandas required.")
            return
        path = filedialog.askopenfilename(title="Open vector file", filetypes=[("Shapefile","*.shp"),("GeoJSON","*.geojson;*.json"),("All","*.*")])
        if not path:
            return
        try:
            gdf = gpd.read_file(path)
            self.last_vector_gdf = gdf
            messagebox.showinfo("Loaded", f"Loaded vector with {len(gdf)} features.")
            self.plot_map_points()
        except Exception as e:
            messagebox.showerror("Load error", str(e))

    # -------------------------
    # Export functions
    # -------------------------
    def export_txt(self):
        pts = self.transformed if self.transformed else self.points
        if not pts:
            messagebox.showwarning("No Data", "No point to export.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text file","*.txt")])
        if not path: return
        with open(path, "w") as f:
            f.write("X\tY\n")
            for x,y in pts:
                f.write(f"{x:.6f}\t{y:.6f}\n")
        messagebox.showinfo("Saved", f"Saved {len(pts)} points to {path}")

    def export_csv(self):
        pts = self.transformed if self.transformed else self.points
        if not pts:
            messagebox.showwarning("No Data","No point to export.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV","*.csv")])
        if not path: return
        df = pd.DataFrame(pts, columns=["X","Y"])
        df.to_csv(path, index=False)
        messagebox.showinfo("Saved", f"Saved {len(pts)} points to {path}")

    def export_xlsx(self):
        pts = self.transformed if self.transformed else self.points
        if not pts:
            messagebox.showwarning("No Data","No point to export.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel","*.xlsx")])
        if not path: return
        df = pd.DataFrame(pts, columns=["X","Y"])
        df.to_excel(path, index=False)
        messagebox.showinfo("Saved", f"Saved {len(pts)} points to {path}")

    def export_geojson(self):
        if not HAS_GIS:
            messagebox.showerror("Missing libs","Geopandas required.")
            return
        pts = self.transformed if self.transformed else self.points
        if not pts:
            messagebox.showwarning("No Data","No point to export.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".geojson", filetypes=[("GeoJSON","*.geojson")])
        if not path: return
        gdf = to_geodf(pts, crs_epsg=32748)  # choose EPSG accordingly or ask user
        gdf.to_file(path, driver="GeoJSON")
        messagebox.showinfo("Saved", f"Saved {len(pts)} points to {path}")

    def export_shapefile(self):
        if not HAS_GIS:
            messagebox.showerror("Missing libs","Geopandas required.")
            return
        pts = self.transformed if self.transformed else self.points
        if not pts:
            messagebox.showwarning("No Data","No point to export.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".shp", filetypes=[("Shapefile","*.shp")])
        if not path: return
        gdf = to_geodf(pts, crs_epsg=32748)
        gdf.to_file(path, driver="ESRI Shapefile")
        messagebox.showinfo("Saved", f"Saved {len(pts)} points to {path}")

# -------------------------
# Run
# -------------------------
def main():
    if USE_BOOTSTRAP:
        app = PROGISApp(None)
        app.root.mainloop()
    else:
        root = tk.Tk()
        app = PROGISApp(root)
        root.mainloop()

if __name__ == "__main__":
    main()
