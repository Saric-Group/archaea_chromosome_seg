"""
Creates a LAMMPS data file with patchy beads arranged on a cubic lattice, for testing purposes with the following scripts:

	- ylz_patchypoly_hc_compactFast_wetUnif.in
	- ylz_patchypoly_hc_compactSlow_wetUnif.in

The present script does not generate membrane beads, and the beads that would normally be connected into 2 polymers are here disconnected.
The purpose of this script is to test the behavior resulting from combining dynamic groups with the command neigh_modify exclude.

It is assumed that the LAMMPS data have the following atom_type, atom_ids and mol_ids:

atom_types:

Polymer beads: atom_type=1
Membrane beads: atom_type=2
Attractive patches: atom_type=3

atom_ids:

Beads of polymer n.1: 	1 <= atom_id <= polymer_length (default 500)
Beads of polymer n.2: 	polymer_length+1 (default 501) <= atom_id <= 2*polymer_length (default 1000)
Patches of polymer n.1: 	2*polymer_length+1 (default 1001) <= atom_id <= 3*polymer_length (default 1500)
Patches of polymer n.2: 	3*polymer_length+1 (default 1501) <= atom_id <= 4*polymer_length (default 2000)
Membrane beads: 			atom_id > 2000

mol_ids:

Each polymer bead has the same mol_id as the patch it is attached to, as (polymer bead + patch) is treated as a single rigid body in the simulation.

Polymer 1: 				1 <= mol_id <= polymer_length (default 500)
Polymer 2: 				polymer_length+1 (default 501) <= mol_id <= 2*polymer_length (default 1000)
"""

def cubic_lattice_points(n_points, density, origin=(0.0, 0.0, 0.0)):
	"""
	Arrange points on a simple cubic lattice, using a number of points
	approximated to the nearest perfect cube.

	Parameters
	----------
	n_points : int
		Target number of points.
	density : float
		Number density (points per unit volume).
	origin : tuple of 3 floats, optional
		Coordinates of the lattice origin.

	Returns
	-------
	coords : ndarray of shape (m**3, 3)
		Cartesian coordinates of the lattice points, where m**3 is the
		perfect cube closest to n_points.
	n_used : int
		Actual number of points used (= m**3).
	spacing : float
		Lattice spacing corresponding to the requested density.
	"""

	if n_points <= 0:
		raise ValueError("n_points must be positive")
	if density <= 0:
		raise ValueError("density must be positive")

	# Find the integer m such that m^3 is closest to n_points
	m = int(rint(n_points ** (1/3)))
	n_used = m**3

	# Lattice spacing from density
	spacing = density ** (-1/3)

	#Box size
	box=spacing*m

	# Generate cubic lattice
	x, y, z = meshgrid(
		arange(m),
		arange(m),
		arange(m),
		indexing="ij"
	)

	coords = column_stack((x.ravel(), y.ravel(), z.ravel()))
	coords = spacing * coords + array(origin)

	return box,n_used,coords 

import os, re, sys, time
from numpy import *
from numpy.random import *
#from scipy.spatial import KDTree
#import trimesh

#########################################################################
#DO NOT CHANGE THESE VALUES (unless you have a very good reason...)
ellipsoidflag_polymer=0
ellipsoidflag_patch=0
ellipsoidflag_sphere=1
#########################################################################

npoints=int(input("Enter number of points:"))
density=float(input("Enter density (points per unit volume):"))

poly_diameter=2.0
dist_neighbors=0.5 		#average distance between neighboring points on the sphere (value should be a bit smaller than the desired one)
density_ellipsoid=6/pi 	#density must be set to 6/pi for ellipsoids of mass 1 and diameter 1
mass_polymer=1.0
mass_patch=0.01

shift=1e-8 	#A small shift of coordinates is needed, otherwise LAMMPS will complain that the rigid bodies particles are the same particle

box,natoms,coors=cubic_lattice_points(npoints,density)
natoms_sphere=0
natoms_tot=2*natoms+natoms_sphere
nbonds=0


with open('patchy_beads_diameter%.2f_natoms%d_lx%.2f_ly%.2f_lz%.2f.lammpsdata'%(poly_diameter,natoms_tot,box,box,box),"w") as fout:
	fout.write("LAMMPS data file\n\n%d atoms\n3 atom types\n%d ellipsoids\n%d bonds\n1 bond types\n\n0.0 %f xlo xhi\n0.0 %f ylo yhi\n0.0 %f zlo zhi\n\nAtoms\n\n"%(natoms_tot,natoms_sphere,nbonds,box,box,box))

	#NB: assuming atom_style hybrid ellispoid bond molecular in LAMMPS input file
	#The format of the Atoms section must be:
	#atom-ID atom-type x y z ellipsoidflag density molecule-ID
	#ellipsoidflag=1 for ylz particles, 0 for point particles
	#IMPORTANT: density must be set to 6/pi for particles of mass 1 and diameter 1
	#IMPORTANT: for particles with ellipsoidflag=0, the density field is interpreted as the particle's mass!

	#mol_id_all=[]
	for k,r in enumerate(coors):
		atom_id=k+1
		x=r[0]
		y=r[1]
		z=r[2]
		xflag=0
		yflag=0
		zflag=0

		mol_id=atom_id #The polymer bead and its patch have the same mol id

		#Write polymer bead coordinates
		atom_type=1
		fout.write("%d %d %.10f %.10f %.10f %d %.6f %d %d %d %d\n"%(atom_id,atom_type,x,y,z,ellipsoidflag_polymer,mass_polymer,mol_id,xflag,yflag,zflag))
		atom_type=3
		fout.write("%d %d %.10f %.10f %.10f %d %.6f %d %d %d %d\n"%(atom_id+natoms,atom_type,x+rand()*shift,y+rand()*shift,z+rand()*shift,ellipsoidflag_patch,mass_patch,mol_id,xflag,yflag,zflag))

