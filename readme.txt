Readme of 2D and 3D LSMPS Python v1.0 by Devvskiii
Important notes : 
-If anything happened, ask Dew AE17. 
line : dalfarisy
WA : 081932038690

Cara Pemakaian : 
1. Cek kelengkapan file yang ada di code. Pastikan lengkap dan tidak corrupt.
List file 2D : 
-Seperangkat file Geom (Geometri)
-ANSYSPlot
-ResultViewer
-plotter
-LSMPSElasticSolver
-LSMPSGeneral
-neighbourfind
List File 3D : 
-Seperangkat file Geom (Geometri)
-LSMPSElasticSolver3D
-LSMPSgeneral
-neighbourfind3d
-Viewer_Geom (File matlab)
-Viewer_Result (File matlab)
-Viewer_ResultANSYS (File matlab)
2. Pastikan anda sudah menginstall python library yg dibutuhkan. Installnya caranya
ada di internet, cari aja.
List library yg dibutuhkan : 
-math
-scipy
-numpy
-array
-math
-datetime
3. Running file geometri apapun dengan python compiler kesukaan anda. File geometri 
ditandai dengan awalan "Geom" atau "Geom3D". Sebelumnya perlu dicek untuk diubah boundary condition
atau tidak. settingan E, v, dan jumlah partikel ada di bagian atas file geom. File geom dapat diubah
datanya sebagai berikut : 
a : panjang geometri
Nx : jumlah partikel yg disebar pada sumbu x
Ny : jumlah partikel yg disebar pada sumbu y
E : Modulus Young
v : Poisson ratio
rho : density
go = konstanta gravitasi
xmid : menentukan posisi titik tengah lingkaran
ymid : menentukan posisi titik tengah lingkaran
radius : radius lingkaran
untuk posisi boundary sudah jelas mana top, right, bottom, left.
cara buat boundary : 
-fixed support boundary disp di semua arah = 0, force segala arah = "nan"
-wilayah yg diberi pembebanan gaya, disp = "nan", force segala arah (berikan nilai sesuai gaya)
4. Setelah file geometri berhasil dirun, seharusnya akan ada file "geom.txt". Apabila tidak ada,
berarti perlu sedikit modifikasi dalam code. ubah bagian bawah sehingga hanya menghasilkan file 
geom.txt.
5. Setelah ada file "geom.txt", ada beberapa settings yang dapat diubah pada code ini, ada di
LSMPSElasticSolver3D atau LSMPSElasticSolver. Settingsnya adalah : 
-ukuran scatter
-neighboursearch method
-scaled cut off ratio
-target neighbour number
-boundary plot method
-skala displacement
-skema kontur,
-data type
-konversi unit dari kN dan mm ke N dan m
-pilihan untuk membuat result.txt atau tidak
6. Setelah pakai settings, silakan run LSMPSElasticSolver atau LSMPSElasticSolver3D dan tunggu
7. Untuk post processing, dapat memakai software berikut : 
2D : Buka resultviewer, nantinya akan dihasilkan beberapa settings seperti apa 
hasil analytic yg akan dibandingkan
3D : buka Viewer_Result, nantinya akan dihasilkan beberapa settings 
seperti apa hasil analytic yg akan dibandingkan, file 3D run pakai matlab karena postproc 
3D Python sangat lambat dan hasilnya sama saja.

