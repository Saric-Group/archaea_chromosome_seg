#############################################
#############################################
Requirements
#############################################
#############################################

- A Bourne shell compatible “Unix” shell program (frequently this is bash)
- Open MPI version 4.1.2 (with libraries)
- GCC version 11.4.0

- Python version 3.10.12 
- Numpy version 1.25.2
- Scipy version 1.11.1
- scikit-learn version is 1.4.1

- LAMMPS version 02/08/2023 (stable release)

#############################################
Compiling LAMMPS (for Unix-based systems)
#############################################

(See also: https://docs.lammps.org/Build_make.html)

The LAMMPS simulations require to slightly modify two files in the LAMMPS source code.
These are:

- pair_cosine_squared.cpp
- pair_cosine_squared.h

The modified files can be found in the attached lammps directory.

The following instructions allow to install LAMMPS using make on Ubuntu or other similar Unix-based systems.

1) Download LAMMPS stable version 02/08/2023 tarball from the following page: https://download.lammps.org/tars/index.html

2) Uncompress tarball: 

	tar -xzvf lammps*.tar.gz

3) Copy modified files pair_cosine_squared.cpp and pair_cosine_squared.h into lammps-2Aug2023/src/EXTRA-PAIR, replacing existing files.

4) Move into the directory lammps-2Aug2023/src

3) Install additional packages (for examples: make yes-rigid)

ASPHERE  
EXTRA-COMPUTE 
EXTRA-MOLECULE 
EXTRA-PAIR  
MOLECULE 
RIGID

5) Install LAMMPS using make: 
	
	make mpi

#############################################
#############################################
LAMMPS simulations
#############################################
#############################################

In the lammps directory, you will find a compiled LAMMPS executable lmp_mpi (compiled on Ubuntu 22.04.5).
We however recommend to recompile the executable, following the steps listed in "Compiling LAMMPS"

#############################################
Initial state creation
#############################################

1) Move into the directory initial_state_creation

2) Run python script create_2rings_noangles_connected.py

3) You will obtain a LAMMPS data file 2rings_connected_N500_box140.06.lammpsdata

4) Run LAMMPS script rings_squeeze_noangles_soft_harmonic_connected.in

5) You will obtain a LAMMPS data file data.compress

6) Run python script data_compress_add_mesh_ellipsoid_atom_longbox_upscale_patchypoly.py

7) You will obtain a LAMMPS data file data_mesh_ylz_patchy_poly_diameter2.00_natoms_sphere6750_lx123.24_ly82.16_lz82.16.lammpsdata - This is the initial configuration for the mobile phase simulation.

#############################################
Mobile phase
#############################################

1) Move into the directory mobile_phase

2) Copy the LAMMPS data file data_mesh_ylz_patchy_poly_diameter2.00_natoms_sphere6750_lx123.24_ly82.16_lz82.16.lammpsdata (initial state) into the directory.

3) Run LAMMPS script ylz_patchypoly_hc_mobile.in:

	Example: ../../lammps/lmp_mpi -in ylz_patchypoly_hc_mobile.in

4) You will obtain a LAMMPS data file data.ylz_patchy_hc_mobile

#############################################
Compaction phase
#############################################

1) Move into one of the directories in the compaction_phase directory (for example, fast_compaction)

2) Copy the LAMMPS data file data.ylz_patchy_hc_mobile, generated at the end of the mobile phase simulation, into the directory and rename it data (an example file is already present in the directory).

3) Run the LAMMPS input file, for example ylz_patchypoly_hc_compactFast_wetUnif.in:

	Example: ../../lammps/lmp_mpi -in ylz_patchypoly_hc_compactFast_wetUnif.in

4) You will obtain a .lammpsdata file (final state after compaction).

#############################################
#############################################
Data Analysis
#############################################
#############################################

The analysis folder contains three python scripts to perform analysis of the LAMMPS simulations.

1) membrane_volume_radius_patchypoly.py

This scripts analyzes the membrane configurations. It approximates the membrane surface using the convex hull method, and returns several quantities:

-natoms.dat 					Number of membrane atoms
-hull_radius_av.dat  			Mean membrane radius
-hull_volume_vs_time.dat  		Membrane volume vs time
-hull_area_vs_time.dat  		Membrane area vs time
-hull_radius_vs_time.dat  		Membrane radius vs time 
-neigh_distance_vs_time.dat 	Mean distance between neighbouring membrane particles vs time
-area_per_bead_vs_time.dat 		Mean area per membrane bead vs time (estimate)

2) lda_segregation_patchy.py

This script performs LDA analysis on the coordinates of the two polymers and returns error fraction f.
The segregation efficiency can be obtaind as s=1-2*f, as detailed in the Methods (see Supporting Information).
It returns:

-lda_frac_errors_vs_time.dat 	Error fraction f vs time
-lda_normal_vs_time.dat 		Normal to the LDA plane (segregation axis) vs time (3 components, x,y,z)
-last_lda_normal.dat  			Last LDA normal of the simulation

3) cm_dist_patchy.py

This script computes the distance between centers of mass (CoMs) of the two polymers. It returns:

-cm_dist_vs_time.dat 			Distance between the CoMs vs time
-cm_dist_resc_rmem_vs_time.dat 	Normalized (by mean membrane radius) distance between the CoMs vs time
-cm_vec_vs_time.dat 			Vector connecting the CoMs vs time (3 components, x,y,z)
