Interpolation of ERA5 fields on a gaussian grid (to generate initial fields for the SWE model):

1) Unzip the 5 files in "data_in" directory : gunzip field_DDMMYYYY_00_025res.dat.gz
2) Edit the file "read_era5_field_for_gg.f90" and choose the resolution of target Gaussian grid
3) Create a binary : gfortran -o main.exe read_era5_field_for_gg.f90 gauleg.f90
4) Execute the binary file main.exe to create the necessary files for run the shallow water model 
