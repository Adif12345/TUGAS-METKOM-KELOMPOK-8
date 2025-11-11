
1. Tujuan Program

Program ini dibuat untuk **menghasilkan dan menampilkan titik-titik koordinat (X dan Y)** secara otomatis.
Program juga bisa **menyimpan hasil koordinat ke file Excel/CSV** dan **menampilkan grafik sebaran titik (scatter plot)** secara interaktif melalui GUI.

---

2. Library yang Digunakan

Dalam kode ini kita menggunakan beberapa library penting untuk membuat program:

* **`tkinter`** → untuk membuat GUI (antarmuka pengguna).
* **`matplotlib`** → untuk menampilkan grafik atau plot koordinat.
* **`pandas`** → untuk menyimpan data ke dalam file CSV/Excel.
* **`numpy`** → untuk perhitungan matematis, terutama ketika membuat titik-titik grid.
* **`random`** → untuk menghasilkan nilai koordinat secara acak.

> Jadi, kombinasi ini membuat program kita punya antarmuka, visualisasi, dan kemampuan penyimpanan data.

---

3. Struktur Program

Program ini berbasis **class** bernama `KoordinatPointGenerator`.

Di dalam class ini ada tiga bagian utama:

1. Bagian Input dan Parameter

   * User bisa mengatur jumlah titik, batas maksimum sumbu X dan Y, serta mode (“Random” atau “Grid”).
   * Tombol “Generate” digunakan untuk membuat titik.
   * Tombol “Save as CSV” digunakan untuk menyimpan hasilnya.

2. Bagian Tabel Data

   * Menampilkan hasil koordinat (X, Y) dalam bentuk tabel.

3. Bagian Grafik (Scatter Plot)

   * Menampilkan visualisasi dari titik-titik yang dihasilkan.

---

4. Penjelasan Tiap Fungsi Utama
🔹 Fungsi `generate_points()`

* Mengambil input dari user (jumlah titik, Xmax, Ymax, mode).
* Jika mode **Random**, maka titik dibuat secara acak menggunakan `random.uniform()`.
* Jika mode **Grid**, titik dibuat secara merata menggunakan `numpy.linspace()`.

> Fungsi ini jadi “otak” utama program, yang menentukan data apa yang akan ditampilkan.

---
🔹 Fungsi `update_table()`

* Membersihkan tabel lama.
* Menambahkan data titik hasil generate ke dalam tabel GUI.

> Tujuannya agar user bisa melihat nilai koordinat setiap titik.

---

🔹 Fungsi `update_plot()`

* Menghapus grafik lama, lalu menggambar titik-titik baru dengan `matplotlib.scatter`.
* Memberikan label pada sumbu X dan Y serta grid agar terlihat rapi.

> Ini membantu pengguna memahami sebaran titik secara visual.

---
🔹 Fungsi `save_csv()`

* Menyimpan hasil koordinat ke dalam file CSV menggunakan `pandas.DataFrame`.
* Menampilkan notifikasi sukses menggunakan `messagebox.showinfo()`.

> Jadi, data yang sudah digenerate tidak hilang — bisa dipakai untuk analisis lebih lanjut di Excel.

---

5. Cara Kerja Program (Alur Eksekusi)

1. Program dijalankan → GUI muncul.
2. User mengisi parameter (N, Xmax, Ymax, mode).
3. Klik tombol **Generate** → titik dibuat dan tampil di tabel serta grafik.
4. Jika ingin disimpan, klik **Save as CSV** → hasil tersimpan otomatis.

---

6. Pengembangan Lanjutan (Opsional)

Program ini bisa dikembangkan menjadi:

* Konversi koordinat ke **UTM / sistem geospasial**.
* Penambahan fitur **rotasi titik**.
* Menampilkan **peta dasar** agar hasilnya lebih geofisika-realistik.

---
7. Kesimpulan

Program ini menunjukkan **Python dalam metode komputasi**, khususnya untuk:

* Pembuatan data koordinat secara otomatis.
* Visualisasi data dalam bentuk grafik.
* Integrasi antarmuka GUI dan penyimpanan data.