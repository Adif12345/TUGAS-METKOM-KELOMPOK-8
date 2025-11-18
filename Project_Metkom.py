# Berikut versi kode Anda dengan penjelasan per baris di samping (komentar)

import tkinter as tk                  # Import modul GUI Tkinter
from tkinter import ttk, filedialog, messagebox  # Import widget tambahan Tkinter
import matplotlib.pyplot as plt        # Import Matplotlib untuk plotting
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg  # Embed plot dalam Tkinter
import pandas as pd                    # Import Pandas untuk ekspor CSV
import numpy as np                     # Import NumPy untuk perhitungan numerik
import random                          # Import random untuk generate titik acak
from pyproj import Proj, transform     # Import pyproj untuk transformasi UTM

class KoordinatPointGeneratorAdvanced:            # Deklarasi kelas utama aplikasi
    def __init__(self, root):                    # Konstruktor kelas
        self.root = root                        # Simpan referensi root Tkinter
        self.root.title("Koordinat Point Generator - Advanced (Metode Komputasi UGM)")  # Judul window
        self.root.geometry("1100x750")          # Ukuran window
        self.points = []                         # List titik asli
        self.transformed_points = []             # List titik hasil transformasi

        frame_param = ttk.LabelFrame(root, text="Parameter Generator", padding=10)  # Frame parameter
        frame_param.pack(fill="x", padx=10, pady=10)  # Layout frame

        ttk.Label(frame_param, text="Jumlah Titik (N):").grid(row=0, column=0)  # Label N
        self.num_points = tk.IntVar(value=25)             # Variabel jumlah titik
        ttk.Entry(frame_param, textvariable=self.num_points, width=8).grid(row=0, column=1)  # Input N

        ttk.Label(frame_param, text="X Max:").grid(row=0, column=2)   # Label X max
        self.x_max = tk.DoubleVar(value=50)         # Variabel X max
        ttk.Entry(frame_param, textvariable=self.x_max, width=8).grid(row=0, column=3)

        ttk.Label(frame_param, text="Y Max:").grid(row=0, column=4)   # Label Y max
        self.y_max = tk.DoubleVar(value=50)         # Variabel Y max
        ttk.Entry(frame_param, textvariable=self.y_max, width=8).grid(row=0, column=5)

        ttk.Label(frame_param, text="Mode:").grid(row=0, column=6)   # Label mode generate
        self.mode = tk.StringVar(value="Random")   # Mode awal Random
        ttk.Combobox(frame_param, textvariable=self.mode, values=["Random", "Grid"], width=10).grid(row=0, column=7)  # ComboBox

        ttk.Button(frame_param, text="Generate", command=self.generate_points).grid(row=0, column=8, padx=10)  # Tombol generate
        ttk.Button(frame_param, text="Clear", command=self.clear_all).grid(row=0, column=9)                # Tombol clear
        ttk.Button(frame_param, text="Save as TXT", command=self.save_txt).grid(row=0, column=10)         # Tombol save txt
        ttk.Button(frame_param, text="Export CSV", command=self.save_csv).grid(row=0, column=11)          # Tombol save csv

        frame_transform = ttk.LabelFrame(root, text="Transformasi Koordinat", padding=10)   # Frame transformasi
        frame_transform.pack(fill="x", padx=10, pady=5)

        ttk.Label(frame_transform, text="Rotasi (derajat):").grid(row=0, column=0)    # Label rotasi
        self.angle = tk.DoubleVar(value=30.0)                # Variabel sudut rotasi
        ttk.Entry(frame_transform, textvariable=self.angle, width=8).grid(row=0, column=1)

        ttk.Button(frame_transform, text="Terapkan Rotasi", command=self.apply_rotation).grid(row=0, column=2, padx=10)  # Tombol rotasi

        ttk.Label(frame_transform, text="Konv