# EarthAlbedo

Code used for the calculation of Earth's albedo for a given satellite, sun location, and set of reflectivity data. Note that this has been ported over from MATLAB's 'Earth Albedo Toolbox', and uses reflectivity data downloaded from https://bhanderi.dk/downloads/ .

## Use 

The primary function is 'albedo(sat, sun, refl)' which takes in a vector containing the satellite and sun positions in the inertial frame (centered on Earth), as well as a set of reflectivity data. The provided data from the above link is from the TOMS dataset and is the average measured reflectivity for each section of Earth over a specified time period.

Note that two additional functions have been provided ('get_albedo_cell_centers' and 'get_diode_albedo'); the first generates the center of each cell in ECEF frame, and the second determines how much of an effect the total albedo has on a given photodiode on the satellite.
