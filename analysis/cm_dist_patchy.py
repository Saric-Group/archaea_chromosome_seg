#!/usr/bin/env python3
"""
This script computes the distance between centers of mass (CoMs) of the two polymers. It returns:
-cm_dist_vs_time.dat 			Distance between the CoMs vs time
-cm_dist_resc_rmem_vs_time.dat 	Normalized (by mean membrane radius) distance between the CoMs vs time
-cm_vec_vs_time.dat 			Vector connecting the CoMs vs time (3 components, x,y,z)

LAMMPS data file format
-----------------------
It is assumed that the LAMMPS data have the following atom_type, atom_ids and mol_ids:

atom_types
----------
Polymer beads: atom_type=1
Membrane beads: atom_type=2
Attractive patches: atom_type=3

atom_ids
----------
Beads of polymer n.1: 	1 <= atom_id <= polymer_length (default 500)
Beads of polymer n.2: 	polymer_length+1 (default 501) <= atom_id <= 2*polymer_length (default 1000)
Patches of polymer n.1: 	2*polymer_length+1 (default 1001) <= atom_id <= 3*polymer_length (default 1500)
Patches of polymer n.2: 	3*polymer_length+1 (default 1501) <= atom_id <= 4*polymer_length (default 2000)
Membrane beads: 			atom_id > 2000

mol_ids
----------
Each polymer bead has the same mol_id as the patch it is attached to, as (polymer bead + patch) is treated as a single rigid body in the simulation.

Polymer 1: 				1 <= mol_id <= polymer_length (default 500)
Polymer 2: 				polymer_length+1 (default 501) <= mol_id <= 2*polymer_length (default 1000)
"""

def compute_cm(coors):
	xcm=average(coors[:,0])
	ycm=average(coors[:,1])
	zcm=average(coors[:,2])
	return array([xcm,ycm,zcm])

import os, sys, time
from numpy import *

folder="dumplin"
polymer_length=500
print('Using polymer length = %d'%polymer_length)

if((folder!='dump') and (folder!='dumplin')):
	print('Folder does not exist (%s)'%folder)
	sys.exit()

path='.'
r_mem=float(loadtxt("%s/membrane_radius_equil.dat"%path))
files=[f for f in os.listdir("%s/%s"%(path,folder)) if ("lammpstrj" in f)]

cmdist_vs_time=array([0,0])
cmdist_vs_time_rescaled=array([0,0])
cmvec_vs_time=array(zeros(4))
for fin in files:
	with open("%s/%s/%s"%(path,folder,fin),'r') as myfile:
		head=[next(myfile) for i in range(9)]
	time=int(head[1].split()[0])				
	xmin=float(head[5].split()[0])
	xmax=float(head[5].split()[1])
	ymin=float(head[6].split()[0])
	ymax=float(head[6].split()[1])
	zmin=float(head[7].split()[0])
	zmax=float(head[7].split()[1])
	my_format=head[8]
	box_x=xmax-xmin
	box_y=ymax-ymin
	box_z=zmax-zmin

	data_tmp=loadtxt("%s/%s/%s"%(path,folder,fin),skiprows=9)	
	data=data_tmp[data_tmp[:,0]==1]  #only take polymers coors

	if('type id mol xu yu zu' in my_format):
		mol_ids=data[:,2].astype('int32')
		poly_ids=where(mol_ids <= polymer_length, 1, 2)
		coordinates=data[:,3:6]
	elif('type mol xu yu zu' in my_format):
		mol_ids=data[:,1].astype('int32')
		poly_ids=where(mol_ids <= polymer_length, 1, 2)
		coordinates=data[:,2:5]
	else:
		print('ERROR: unrecognized format!')

	poly1_coor=coordinates[poly_ids==1]
	poly2_coor=coordinates[poly_ids==2]
	poly1_cm=compute_cm(poly1_coor)
	poly2_cm=compute_cm(poly2_coor)

	cmvec=poly1_cm-poly2_cm

	#dist=sqrt((poly1_cm[0]-poly2_cm[0])**2+(poly1_cm[1]-poly2_cm[1])**2+(poly1_cm[2]-poly2_cm[2])**2)
	dist=linalg.norm(cmvec)
	dist_rescaled=dist/r_mem #rescale with membrane diameter

	cmvec_vs_time=vstack((cmvec_vs_time,concatenate(([time],cmvec/dist)))) #normalized vector connecting the centers of mass
	cmdist_vs_time=vstack((cmdist_vs_time,array([time,dist])))
	cmdist_vs_time_rescaled=vstack((cmdist_vs_time_rescaled,array([time,dist_rescaled])))

cmdist_vs_time_sorted=cmdist_vs_time[cmdist_vs_time[:,0].argsort()][1:]
savetxt("%s/cm_dist_vs_time.dat"%path,cmdist_vs_time_sorted,fmt='%d %.6f')
cmdist_vs_time_rescaled_sorted=cmdist_vs_time_rescaled[cmdist_vs_time_rescaled[:,0].argsort()][1:]
savetxt("%s/cm_dist_resc_rmem_vs_time.dat"%path,cmdist_vs_time_rescaled_sorted,fmt='%d %.6f')
cmvec_vs_time_sorted=cmvec_vs_time[cmvec_vs_time[:,0].argsort()][1:]
savetxt("%s/cm_vec_vs_time.dat"%path,cmvec_vs_time_sorted,fmt='%d %.6f %.6f %.6f')
