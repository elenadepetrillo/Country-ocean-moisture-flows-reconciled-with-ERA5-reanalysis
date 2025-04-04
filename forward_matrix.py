#FORWARD BILATERAL MATRIX OF ATMOSPHERIC MOISTURE FLOWS BETWEEN COUNTRIES/OCEANS from the UTrack dataset (https://doi.org/10.1594/PANGAEA.912710, https://github.com/ObbeTuinenburg/UTrack_global_database)

# Created by Elena De Petrillo, Simon Felix Fahrlaender, Luca Monaco
# elena.depetrillo@polito.it, simon.fahrlaender@pik-potsdam.de, luca.monaco@polito.it

# Released 28 March 2024
''''
The Utrack dataset is a 4D netcdf gridded world of monthly moisture connections (given as probabilities ranging from 0 to 1) between source and target cells. 
The first two coordinates of the dataset (source_lat, source_lon) describe the location of a source cell in the world, while the last two (target_lat, target_lon) describe the location of a target in the world.

This script extracts from the dataset the probability that an evapotranspiration particle from a source cell precipitates into a target t, i.e. a 'forward footprint' (fw_fp).
The fw_fp is given as a probability grid over the whole world, for each source cell s and for each month m. 
For each specific month, the forward footprint fw_fp is multiplied by the days of the month and the evapotranspiration (ET) from ERA5 of source cell s, to obtain the monthly moisture flux. Further multiplying 
by the area of the source cell, the moisture volumes are obtained. The moisture volumes are Precipitation (TP) in m^3 at target cell t and Evapotranspiration (ET) in m^3 at source cell s.
To obtain the annual moisture flows, the volumes for each month are cumulated. Subsequently, the cell-scale volumes are aggregated to country/ocean scale.

The reference year is the average year between 2008-2017, so ERA5 ET data are averaged over this time window when extracted.

The script runs using cell-scale parallel programming.

The output is a square matrix in which each row and column index is associated with a country/ocean code (available in the input list).
The row shows the annual evaporation between one country/ocean and another, while the column shows the annual precipitation between one country/ocean and another.

The sum over the row is the total evaporation of a country/ocean, while the sum over the column is its total precipitation in the average year.

'''

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4 as nc
import multiprocessing

#Get the days in the m month
def get_days_in_month(month):
    days_in_month = {
        '01': 31,
        '02': 28,
        '03': 31,
        '04': 30,
        '05': 31,
        '06': 30,
        '07': 31,
        '08': 31,
        '09': 30,
        '10': 31,
        '11': 30,
        '12': 31
    }
    return days_in_month.get(month, None)

#UTrack spatial resolution: choose between 1.0 or 0.5
res = "0.5"

#Number of processors to run the multiprocessing pool at the cell-scale
num_processors = 

root_path = "your root path_"
era5_path = root_path
input_nc_path = root_path + "/input_nc"
input_csv_path = root_path + "/input_csv"
script_path = root_path + "/Script_"
output_path = root_path + "/your output folder"
file_path = root_path + "/your file path"

if res == "1.0": 
    step = 1
elif res == "0.5":
    step = 0.5

lats = np.arange(90, -90, -step)
lons = np.arange(0, 360, step)

#rearrange longitudes to (-180 -- 180 ) for plotting
lons_new = np.arange(-180, 180, step) 

lons_mat, lats_mat = np.meshgrid(lons_new, lats)
months = [str(i).zfill(2) for i in range(1, 13)]

#load countries and oceans codes and names list
country_ocean_file = pd.read_csv(input_csv_path + "/country_list_cntr.csv", sep=";", keep_default_na=False) 
country_ocean_iso = country_ocean_file ['FID'].astype(str).to_numpy()
country_ocean_list = country_ocean_file ['NUM_ID'].astype(int).to_numpy()

# to run the code for selected countries
    # iso_to_index = {iso: index for index, iso in enumerate(country_iso)} 
    # country_oceans_indices = []
    # country_oceans_indices.append(iso_to_index['ISO_country_1'])
    # country_oceans_indices.append(iso_to_index['ISO_country_1'])


#load the countries and oceans raster mask
mask = xr.open_dataset(input_nc_path + '/CNTR_OC_mask_EUROSTAT_GOaS.nc')
#elige the Band
country_ocean_mask = mask['Band1'].values.astype(int) 

#function to retrive UTrack fraction of evaporation to precipitation in the forward case at the cell scale
#the function will go in the parallel run
def country_cell_loop(c_cell):
    lat = lats[c_cell[0]]
    lon = lons[c_cell[1]]
    fw_fp = mf.sel(sourcelat=lat, sourcelon=lon).values
    mask = np.where(fw_fp == 0)
    fw_fp = fw_fp * (-0.1)
    fw_fp = np.exp(fw_fp)
    fw_fp[mask] = 0
    fw_fp = fw_fp / np.sum(fw_fp)
    et = EvapoTransp.sel(lat=lat, lon=lon).values * days * (-1) 
    area=cell_area.sel(lat=lat, lon=lon).values #extract the area for the specific cell
    fw_fp = et * fw_fp * area
    return fw_fp, area

matrix_fw = np.zeros((len(country_ocean_iso), len(country_ocean_iso))) #define the bilateral forward moisture volumes matrix to be given in output by this file

#begin the cycle over months
for m in months:
    #for each month extract the utrack file, Evapotranspiration and Total Precipitation
    utrack_clim = xr.open_dataset(input_nc_path + "/utrack_climatology_" + res + "_" + m + ".nc")
    mf = utrack_clim["moisture_flow"]
    utrack_clim.close()
    ERA5 = xr.open_dataset(input_nc_path + "/ERA5_" + m + "_" + res + ".nc")
    EvapoTransp = ERA5["e"].mean(dim="time")
    TotalPrec = ERA5["tp"].mean(dim="time")
    ERA5.close()

    #Extract the world cell areas
    Map_area = xr.open_dataset(input_nc_path + "/area_" + res + ".nc")
    cell_area = Map_area["cell_area"]
    area_map = cell_area.values
    Map_area.close()

    days = get_days_in_month(m) #get the days in the m mounth

    # in case the analysis is perfomerd only for some countries/oceans use:
    #for i, c_idx in enumerate(country_indices): 
    
    #in case the analysis runs over all the countries/oceans
    for i, c in enumerate(country_ocean_list): #begin cycle over countries/oceans
        iso = country_ocean_iso[i]  # Get the ISO code of the current country/ocean
        if iso in country_ocean_iso:
            country_mask_indices = np.where(country_ocean_mask == int(c)) #get the indices of the cells who are contained by the country/ocean c
            country_cell = list(zip(country_mask_indices[0], country_mask_indices[1])) #reshape them as lists of points (lat, lon)
            pool = multiprocessing.Pool(num_processors) #define the pool to be used in parallel programming
            country_output = pool.map(country_cell_loop, country_cell) #parallel programming, feeding country_cell to country_cell_loop
            pool.close()
            pool.join()
            cumulated_flux = np.zeros((lats.shape[0], lons.shape[0])) #define a variable to cumulate the world cell_scale volumes from country_cell from the pool

            for o in country_output: #cumulate
                cumulated_flux += o[0]

	    #Now cumulated_flux must be subdivided into contributions from -potentally- every country/ocean in the world
            cumulated_flux_by_country = np.zeros(len(country_ocean_iso))
            for c1 in range(len(country_ocean_iso)):#begin another cycle over countries/oceans
                country_cell1 = np.where(country_ocean_mask == country_ocean_list[c1])
                flux_cellbycell = cumulated_flux[country_cell1] #extract the contribution of every cell from the country/ocean having index c1
                flux_cellbycell = np.where(np.isnan(flux_cellbycell), 0, flux_cellbycell) #correct any NaN
                somma = np.nansum(flux_cellbycell) #sum the contribution of every cell from the country/ocean having index c1
                cumulated_flux_by_country[c1] = somma #assign it

            #fill the row of the matrix as the evaporation of the i- country
            matrix_fw[i, :] += cumulated_flux_by_country #assign the row sum to the output matrix

np.savetxt(file_path, matrix_fw, delimiter=';') #save the output matrix

print("Saved")
