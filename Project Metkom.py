import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import numpy as np
import random
from pyproj import Proj, transform

class KoordinatPointGeneratorAdvanced:
    def __init__(self, root):
        self.root = root
        self.root.title("Koordinat Point Generator - Advanced (Metode Komputasi UGM)")
        self.root.geometry("1100x750")
        self.points = []
        self.transformed_points = []

        frame_param = ttk.LabelFrame(root, text="Parameter Generator", padding=10)
        frame_param.pack(fill="x", padx=10, pady=10)

        ttk.Label(frame_param, text="Jumlah Titik (N):").grid(row=0, column=0)
        self.num_points = tk.IntVar(value=25)
        ttk.Entry(frame_param, textvariable=self.num_points, width=8).grid(row=0, column=1)

        ttk.Label(frame_param, text="X Max:").grid(row=0, column=2)
        self.x_max = tk.DoubleVar(value=50)
        ttk.Entry(frame_param, textvariable=self.x_max, width=8).grid(row=0, column=3)

        ttk.Label(frame_param, text="Y Max:").grid(row=0, column=4)
        self.y_max = tk.DoubleVar(value=50)
        ttk.Entry(frame_param, textvariable=self.y_max, width=8).grid(row=0, column=5)

        ttk.Label(frame_param, text="Mode:").grid(row=0, column=6)
        self.mode = tk.StringVar(value="Random")
        ttk.Combobox(frame_param, textvariable=self.mode, values=["Random", "Grid"], width=10).grid(row=0, column=7)

        ttk.Button(frame_param, text="Generate", command=self.generate_points).grid(row=0, column=8, padx=10)
        ttk.Button(frame_param, text="Clear", command=self.clear_all).grid(row=0, column=9)
        ttk.Button(frame_param, text="Save as TXT", command=self.save_txt).grid(row=0, column=10)
        ttk.Button(frame_param, text="Export CSV", command=self.save_csv).grid(row=0, column=11)

        frame_transform = ttk.LabelFrame(root, text="Transformasi Koordinat", padding=10)
        frame_transform.pack(fill="x", padx=10, pady=5)

        ttk.Label(frame_transform, text="Rotasi (derajat):").grid(row=0, column=0)
        self.angle = tk.DoubleVar(value=30.0)
        ttk.Entry(frame_transform, textvariable=self.angle, width=8).grid(row=0, column=1)

        ttk.Button(frame_transform, text="Terapkan Rotasi", command=self.apply_rotation).grid(row=0, column=2, padx=10)

        ttk.Label(frame_transform, text="Konversi ke UTM Zona:").grid(row=0, column=3)
        self.utm_zone = tk.IntVar(value=48)  # default Indonesia tengah
        ttk.Entry(frame_transform, textvariable=self.utm_zone, width=5).grid(row=0, column=4)
        ttk.Button(frame_transform, text="Konversi UTM", command=self.apply_utm).grid(row=0, column=5, padx=10)

        frame_table = ttk.LabelFrame(root, text="Tabel Koordinat (X, Y)", padding=10)
        frame_table.pack(fill="both", expand=True, padx=10, pady=10)

        self.tree = ttk.Treeview(frame_table, columns=("x", "y"), show="headings")
        self.tree.heading("x", text="X")
        self.tree.heading("y", text="Y")
        self.tree.column("x", width=150, anchor="center")
        self.tree.column("y", width=150, anchor="center")
        self.tree.pack(fill="both", expand=True)

        frame_plot = ttk.LabelFrame(root, text="Scatter Plot", padding=10)
        frame_plot.pack(fill="both", expand=True, padx=10, pady=10)

        self.figure, self.ax = plt.subplots(figsize=(6, 5))
        self.canvas = FigureCanvasTkAgg(self.figure, master=frame_plot)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def generate_points(self):
        N = self.num_points.get()
        xmax = self.x_max.get()
        ymax = self.y_max.get()
        mode = self.mode.get()

        if mode == "Random":
            self.points = [(random.uniform(0, xmax), random.uniform(0, ymax)) for _ in range(N)]
        else:
            nx = int(np.sqrt(N))
            ny = nx
            x = np.linspace(0, xmax, nx)
            y = np.linspace(0, ymax, ny)
            self.points = [(xi, yi) for xi in x for yi in y]

        self.transformed_points = []
        self.update_table(self.points)
        self.update_plot()

    def update_table(self, points):
        for i in self.tree.get_children():
            self.tree.delete(i)
        for x, y in points:
            self.tree.insert("", "end", values=(f"{x:.3f}", f"{y:.3f}"))

    def update_plot(self):
        self.ax.clear()
        if self.points:
            xs, ys = zip(*self.points)
            self.ax.scatter(xs, ys, color='dodgerblue', label='Asli')
        if self.transformed_points:
            xs2, ys2 = zip(*self.transformed_points)
            self.ax.scatter(xs2, ys2, color='red', label='Transformasi')
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_title("Scatter Plot Koordinat")
        self.ax.legend()
        self.ax.grid(True)
        self.canvas.draw()

    def clear_all(self):
        self.points = []
        self.transformed_points = []
        for i in self.tree.get_children():
            self.tree.delete(i)
        self.ax.clear()
        self.canvas.draw()

    def apply_rotation(self):
        if not self.points:
            messagebox.showwarning("Peringatan", "Generate titik dulu sebelum rotasi!")
            return
        angle_rad = np.radians(self.angle.get())
        cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
        self.transformed_points = [(x * cos_a - y * sin_a, x * sin_a + y * cos_a) for x, y in self.points]
        self.update_plot()
        messagebox.showinfo("Transformasi", f"Titik berhasil dirotasi {self.angle.get()}°")

    def apply_utm(self):
        if not self.points:
            messagebox.showwarning("Peringatan", "Generate titik dulu sebelum konversi!")
            return
        zone = self.utm_zone.get()
        proj_geo = Proj(proj='latlong', datum='WGS84')
        proj_utm = Proj(proj='utm', zone=zone, datum='WGS84')

        longitudes = [x for x, _ in self.points]
        latitudes = [y for _, y in self.points]

        x_utm, y_utm = transform(proj_geo, proj_utm, longitudes, latitudes)
        self.transformed_points = list(zip(x_utm, y_utm))
        self.update_plot()
        messagebox.showinfo("Transformasi", f"Koordinat berhasil dikonversi ke UTM zona {zone}")

    def save_txt(self):
        if not self.points:
            messagebox.showwarning("Peringatan", "Tidak ada data untuk disimpan.")
            return
        file = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text file", "*.txt")])
        if file:
            with open(file, "w") as f:
                f.write("X\tY\n")
                for x, y in self.transformed_points or self.points:
                    f.write(f"{x:.3f}\t{y:.3f}\n")
            messagebox.showinfo("Sukses", "Data disimpan ke file TXT.")

    def save_csv(self):
        if not self.points:
            messagebox.showwarning("Peringatan", "Tidak ada data untuk disimpan.")
            return
        file = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV file", "*.csv")])
        if file:
            df = pd.DataFrame(self.transformed_points or self.points, columns=["X", "Y"])
            df.to_csv(file, index=False)
            messagebox.showinfo("Sukses", "Data disimpan ke file CSV.")

if __name__ == "__main__":
    root = tk.Tk()
    app = KoordinatPointGeneratorAdvanced(root)
    root.mainloop()
