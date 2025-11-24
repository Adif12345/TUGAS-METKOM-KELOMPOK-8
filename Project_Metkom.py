"""
Koordinat Point Generator — PRO / GIS MINI (Square & Azimuth Only)
Versi untuk VSCode (single-file app)

Fitur utama:
- GUI dengan mode generate: square dan line by azimuth
- Input azimuth hanya muncul jika mode azimuth dipilih
- Visualisasi scatter plot
- Export: TXT, CSV, XLSX
"""

import math
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tkinter.scrolledtext import ScrolledText

try:
    import ttkbootstrap as tb
    from ttkbootstrap.constants import *
    THEME = "flatly"
    USE_BOOTSTRAP = True
except Exception:
    USE_BOOTSTRAP = False

# -------------------------
# Geometry generation
# -------------------------
def generate_square_grid(n_side, spacing, x0, y0):
    n = int(n_side)
    xs = np.arange(0, n * spacing, spacing)
    ys = np.arange(0, n * spacing, spacing)
    pts = [(x0 + x, y0 + y) for x in xs for y in ys]
    return pts

def generate_line_azimuth(x0, y0, spacing, n_points, azimuth_deg):
    az = math.radians(azimuth_deg)
    pts = [(x0 + i * spacing * math.sin(az), y0 + i * spacing * math.cos(az)) for i in range(n_points)]
    return pts

# -------------------------
# Main GUI Application
# -------------------------
class PROGISApp:
    def __init__(self, master):
        if USE_BOOTSTRAP:
            self.root = tb.Window(themename=THEME, title="Koordinat Point Generator - PRO GIS (Square & Azimuth)")
            self.master = self.root
        else:
            self.root = master
            self.master = master
            master.title("Koordinat Point Generator - PRO GIS (Square & Azimuth)")
            master.geometry("1000x700")

        self.points = []
        self._build_ui()

    def _build_ui(self):
        self.tab_control = ttk.Notebook(self.master)
        self.tab_generate = ttk.Frame(self.tab_control)
        self.tab_export = ttk.Frame(self.tab_control)
        self.tab_about = ttk.Frame(self.tab_control)

        self.tab_control.add(self.tab_generate, text="Generate")
        self.tab_control.add(self.tab_export, text="Export")
        self.tab_control.add(self.tab_about, text="About")
        self.tab_control.pack(expand=1, fill="both")

        # GENERATE TAB
        lf = ttk.LabelFrame(self.tab_generate, text="Options", padding=8)
        lf.pack(side="left", fill="y", padx=8, pady=8)

        ttk.Label(lf, text="Mode:").grid(row=0, column=0, sticky="w")
        self.mode_var = tk.StringVar(value="square")
        mode_box = ttk.Combobox(lf, textvariable=self.mode_var, values=["square", "line_azimuth"], width=18)
        mode_box.grid(row=0, column=1, padx=4, pady=2)
        mode_box.bind("<<ComboboxSelected>>", self.toggle_azimuth_option)

        ttk.Label(lf, text="Start X:").grid(row=1, column=0, sticky="w")
        self.start_x = tk.DoubleVar(value=500000.0)
        ttk.Entry(lf, textvariable=self.start_x, width=12).grid(row=1, column=1)

        ttk.Label(lf, text="Start Y:").grid(row=2, column=0, sticky="w")
        self.start_y = tk.DoubleVar(value=908000.0)
        ttk.Entry(lf, textvariable=self.start_y, width=12).grid(row=2, column=1)

        ttk.Label(lf, text="Spacing (m):").grid(row=3, column=0, sticky="w")
        self.spacing = tk.DoubleVar(value=10.0)
        ttk.Entry(lf, textvariable=self.spacing, width=12).grid(row=3, column=1)

        ttk.Label(lf, text="Total Points:").grid(row=4, column=0, sticky="w")
        self.count_side = tk.IntVar(value=10)
        ttk.Entry(lf, textvariable=self.count_side, width=12).grid(row=4, column=1)

        # azimuth field (hidden by default)
        self.azimuth_label = ttk.Label(lf, text="Azimuth (°):")
        self.azimuth_entry_var = tk.DoubleVar(value=45.0)
        self.azimuth_entry = ttk.Entry(lf, textvariable=self.azimuth_entry_var, width=12)

        ttk.Button(lf, text="Generate", command=self.on_generate).grid(row=6, column=0, pady=6)
        ttk.Button(lf, text="Clear", command=self.on_clear).grid(row=6, column=1, pady=6)

        right_frame = ttk.Frame(self.tab_generate)
        right_frame.pack(side="left", fill="both", expand=True, padx=8, pady=8)

        self.fig, self.ax = plt.subplots(figsize=(7, 6))
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        self.canvas_plot = FigureCanvasTkAgg(self.fig, master=right_frame)
        self.canvas_plot.get_tk_widget().pack(fill="both", expand=True)

        self.table = ttk.Treeview(right_frame, columns=("x", "y"), show="headings", height=8)
        self.table.heading("x", text="X")
        self.table.heading("y", text="Y")
        self.table.pack(fill="x")

        # EXPORT TAB
        ef = ttk.LabelFrame(self.tab_export, text="Export Options", padding=8)
        ef.pack(fill="both", expand=True, padx=8, pady=8)

        ttk.Button(ef, text="Export TXT", command=self.export_txt).pack(anchor="w")
        ttk.Button(ef, text="Export CSV", command=self.export_csv).pack(anchor="w")
        ttk.Button(ef, text="Export XLSX", command=self.export_xlsx).pack(anchor="w")

        # ABOUT TAB
        sf = ttk.LabelFrame(self.tab_about, text="About", padding=8)
        sf.pack(fill="both", expand=True, padx=8, pady=8)
        txt = ScrolledText(sf, height=15)
        txt.insert("1.0", "Koordinat Point Generator - PRO GIS MINI (Square & Azimuth)\n\nAplikasi sederhana untuk menghasilkan titik koordinat otomatis.\n\nFitur:\n- Generate titik persegi dan garis berdasarkan azimuth\n- Input azimuth hanya muncul jika mode azimuth dipilih\n- Visualisasi scatter plot\n- Export hasil ke TXT, CSV, XLSX\n")
        txt.configure(state="disabled")
        txt.pack(fill="both", expand=True)

    def toggle_azimuth_option(self, event=None):
        mode = self.mode_var.get()
        if mode == "line_azimuth":
            self.azimuth_label.grid(row=5, column=0, sticky="w")
            self.azimuth_entry.grid(row=5, column=1)
        else:
            self.azimuth_label.grid_remove()
            self.azimuth_entry.grid_remove()

    def on_generate(self):
        mode = self.mode_var.get()
        x0 = self.start_x.get()
        y0 = self.start_y.get()
        spacing = self.spacing.get()
        cnt = self.count_side.get()

        if mode == "square":
            pts = generate_square_grid(cnt, spacing, x0, y0)
        elif mode == "line_azimuth":
            azimuth = self.azimuth_entry_var.get()
            pts = generate_line_azimuth(x0, y0, spacing, cnt, azimuth)
        else:
            pts = []

        self.points = pts
        self._update_table(pts)
        self._plot_points(pts, title=f"Generated: {mode}")

    def on_clear(self):
        self.points = []
        self._update_table([])
        self._plot_points([], title="Cleared")

    def _update_table(self, pts):
        for item in self.table.get_children():
            self.table.delete(item)
        for x, y in pts:
            self.table.insert("", "end", values=(f"{x:.3f}", f"{y:.3f}"))

    def _plot_points(self, pts, title="Plot"):
        self.ax.clear()
        if pts:
            xs, ys = zip(*pts)
            self.ax.scatter(xs, ys, s=10)
        self.ax.set_title(title)
        self.ax.grid(True)
        self.canvas_plot.draw()

    def export_txt(self):
        if not self.points:
            messagebox.showwarning("No Data", "No points to export.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text file", "*.txt")])
        if not path:
            return
        with open(path, "w") as f:
            f.write("X\tY\n")
            for x, y in self.points:
                f.write(f"{x:.6f}\t{y:.6f}\n")
        messagebox.showinfo("Saved", f"Saved {len(self.points)} points to {path}")

    def export_csv(self):
        if not self.points:
            messagebox.showwarning("No Data", "No points to export.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not path:
            return
        df = pd.DataFrame(self.points, columns=["X", "Y"])
        df.to_csv(path, index=False)
        messagebox.showinfo("Saved", f"Saved {len(self.points)} points to {path}")

    def export_xlsx(self):
        if not self.points:
            messagebox.showwarning("No Data", "No points to export.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel", "*.xlsx")])
        if not path:
            return
        df = pd.DataFrame(self.points, columns=["X", "Y"])
        df.to_excel(path, index=False)
        messagebox.showinfo("Saved", f"Saved {len(self.points)} points to {path}")

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
