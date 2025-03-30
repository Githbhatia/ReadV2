# ReadV2

Read and plot CISMIP formated v2 OR COSMOS formated V2c seismic ground motion records such as posted on CESMD website (Center of Strong Motion Data https://www.strongmotioncenter.org/). The code can also be used to view building instrument records on the CESMD site and those also available at the HCAI website (https://hcai.ca.gov/construction-finance/facility-detail/ - navigate to a hospital that has instrumented buildings and look under the Instrumented Buildings Tab). Python code reads a .v2 or V2c files that contains one or three channels (Free-Field instruments have 3 channels, Instrumented Buildings have records for individual channels.)

User interface uses tkinter.

This code functions as a viewer for the ground motion data contained in the v2 file (1 file containing 3 channels with acceleration, velocity or displacment) or V2c files (9 containing each containing acceleration, velocity or displacement for 3 channels). Zip files containing the v2 or V2c files can be directly read including when they are double zipped as is the case when downloaded from the CESMD website (Center of Strong Motion Data https://www.strongmotioncenter.org/).

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
 *Added D5-75 and D5-95 calculations reported on the Arias Intensity plots - easily changed to any interval in code.

 Changes 12/7/2024
 *Revised the 3D-orbit plot to have equal axis in all three directions.  Added option to create response spectra for velocity and displacement in addition to acceleration as requested by a user.

Changes 12/8/2024
*Revised 3D-orbit plot to have colors based on z displacements (up-down).  Changed the main gui window to be scrollable for users with lower resolution monitors.

Changes 12/22/2024
*Animated 3D-orbit plot  with play, stop and reverse buttons with accompanying side plots that can be acceleration, velocity or displacement.  Circular gridlines added to orbit plots for rotated directions.  v2 files can be zip archives as downloaded from CESMD or CSMIP sites. Other usability improvements.

Changes 12/23/2024
*As some downloaded CESMD records are zipped once and others zipped twice, added code to accept either.

Changes 12/23/2024
*Added markers for animations.  Fixed all the window titles with record time and date.  Other usability improvements. Standardize directions for 3D orbit plot as in other orbit plots.

Changes 12/29/2024
*Added option to read COSMOS V2c files.  Other improvements to speed up animation and spectra calculations (now uses only selected portion of the record).

Changes 3/29/2025
*Added revised file name schema used by USGS for v2C files for some internation stations.
