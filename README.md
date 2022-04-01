# EarthAlbedo.jl

This package holds code used for the calculation of Earth's albedo for a given satellite, sun location, and set of reflectivity data. Note that this has been ported over from MATLAB's 'Earth Albedo Toolbox', and uses reflectivity data downloaded from https://bhanderi.dk/downloads/ .

## Downloading reflectivity data

The reflectivity data is based off of a set measured reflectivity data known as the TOMS dataset. The surface of the Earth is divided into cells and the TOMS data corresponding to each cell are averaged and used to estimate the reflectivity of a given area of the Earth's surface. 

A version of this data has been processed and saved in the 'data/' folder. The source data can be downloaded by selecting the MATLAB converted TOMS data corresponding to the desired year from https://bhanderi.dk/downloads/  (e.g., refl = matread("tomsdata2005/2005/ga050101-051231.mat"), and then creating a struct.

Note that the 'Earth Albedo Toolbox' does not need to be downloaded for this to work.

## Use 

The primary function is 'albedo(sat, sun, refl_data)' which takes in a vector containing the satellite and sun positions in the inertial frame (centered on Earth), as well as a set of reflectivity data. The provided data from the above link is from the TOMS dataset and is the average measured reflectivity for each section of Earth over a specified time period.