# ReadV2

Read and plot COSMOS formate v2 seismic ground motion records such as posted on CESMD website (Center of Strong Motion Data https://www.strongmotioncenter.org/). The code can also be used to view building instrument records on the CESMD site and those also available at the HCAI website (https://hcai.ca.gov/construction-finance/facility-detail/ - navigate to a hospital that has instrumented buildings and look under the Instrumented Buildings Tab). Python code reads a .v2 file that contains one or three channels (Free-Field instruments have 3 channels, Instrumented Buildings have records for individual channels.)

User interface uses tkinter.

This code functions as a viewer for the ground motion data contained in the v2 file.  
Code shows location of seismic instrument that originated the record on a map. Plots acceleration, integrated velocity and displacement time history for each component. Plot response spectra in SA vs Time Period or ADRS format - computes energy content of each component. Plots 3D orbit plots. Plots rotated resultant in the maximum acceleration, velocity or displacement directions. Create RotD50, RotD00, RotD100 response spectra. Compare to Geomean Spectra. Plot resultant spectra in a Tripartite format. Compare to ASCE 7-22 design spectra using any coordinates in the US - default is the location of the instrument (that is, using the latitude and longitude of the instrument). 

Save time vs. acceleration in a text format.

Example input v2 file from the Ferndale earthquake is included.  See png files to see sample outputs.

Needs many Python packages: numpy matplotlib itertools tkintner scipy


Changes 5/3/2023
  *Added option to plot Arias Intensity
  *Changed ASCE7-22 url per changes by USGS
  
Changes 5/19/2023
 *Added option to rotate to a specified angle

 Changes 10/17/2024
 *Added D5-75 and D5-95 calculations reported on the Arias Intensity plots - easily changed to any interval.
