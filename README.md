hydro_analysis
==============

Analyze the hydrodynamic evolution file

* SurfaceFinder class use Cornelius algorithm to find iso-thermal hyper-surface
in the hydrodynamic medium. Cornelius is developed and written by Pasi Huovinen
and Hannah Peterson. It is an open source code through OSCAR project. The project
webpage can be found at https://karman.physics.purdue.edu/OSCAR/index.php/CORNELIUS.
The copyright of this subroutine can be found there.

* FluidcellStatistic class include various kinds of flow analysis for the hydrodynamic
medium, including compute the evolution of the flow velocity, spatial eccentricity,
momentum anisotropy, Knudsen and inverse Reynold's number for hydrodynamic fluid cells.

* Hydroinfo_MUSIC class is added to support read in hydrodynamic evolution from 
MUSIC outputs.

To do:
* add cmake compile files
