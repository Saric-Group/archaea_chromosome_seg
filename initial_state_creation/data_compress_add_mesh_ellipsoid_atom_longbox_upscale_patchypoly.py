#Adds a spherical mesh membrane to a configurations of polymers 
#initially enclosed in a sphere of radius r_wall.
#
#This version additionally upscales the size of the polymer beads,
#so that the polymer bead diameter is larger than the membrane bead one.
#
#In this version, the polymer beads are rigid bodies with the same center of mass.
#The rigid body is made of a type 1 particle and a type 3 particle, whereas the membrane beads are type 2.
#Particles belonging to the same rigid bodies must have the same molecular id
#
#Using rigid bodies is required to split the repulsive and attractive contributions of the interaction 
#between two different particles

import os, re, sys, time
from numpy import *
from numpy.random import *
from scipy.spatial import KDTree
import trimesh

#=================================
#Quaternion stuff
#=================================

def unit_vector(vector):
   """ Returns the unit vector of the vector. """
   if max(abs(vector))>1e-12:
      return vector / linalg.norm(vector)
   else:
      return vector

def angle_between(v1, v2):
   """ Returns the angle between vectors 'v1' and 'v2' """
   v1_u = unit_vector(v1)
   v2_u = unit_vector(v2)
   return arccos(clip(dot(v1_u, v2_u), -1.0, 1.0))

def distance_between(v1, v2):
   """ Returns the distance between position vectors 'v1' and 'v2' """
   return sqrt(sum(square(array(v1)-array(v2))))

def get_vrot_ang(orient):
   """ Get rotation axis and angle for quaternion
   Quaternions are by default pointing along the x directions. This function
   first finds the axis along which it has to rotate the quaternion and then
   the angle by which it has to rotate it in the right hand rule sense
   """
   xhat = array([1,0,0])
   if   isclose( dot(unit_vector(orient),xhat), 1,rtol=0, atol=1e-6 ):
      vrot = xhat
      ang = 0
   elif isclose( dot(unit_vector(orient),xhat),-1,rtol=0, atol=1e-6 ):
      vrot = array([0,1,0])
      ang = pi
   else:
      vrot = unit_vector(cross(xhat,orient))
      ang = angle_between(xhat,orient)
   return vrot,ang

def average_nearest_neighbor_distance(points):
    tree = KDTree(points)
    distances, _ = tree.query(points, k=2)  # k=2 because the nearest point includes the point itself
    return mean(distances[:, 1])  # Skip distance to itself which is 0

#=================================
#Code for quaternions
#=================================
# for i in range(0, data.particles.count):
#   z, x, y = data.particles.positions[i]  # Some routine that gives you the position of particles
#   orient = unit_vector(array([x,y,z]))
#   vrot, theta = get_vrot_ang(orient)  # vrot in cartesians, theta in radians
#   a, b, c = vrot
#   quatw, quati, quatj, quatk = cos(theta / 2), a * sin(theta / 2), b * sin(theta / 2), c * sin(theta / 2) # these is the quaternion you need to put in the LAMMPS input file

def generate_quaternions(vertices,box_x,box_y,box_z):
	#Takes coordinates of the vertices and uses them to generate quaternions components
	quaternions=empty((0,4))
	for r in vertices:
		r_shift=array([r[0]-0.5*box_x,r[1]-0.5*box_y,r[2]-0.5*box_z])
		orient = unit_vector(r_shift)
		vrot, theta = get_vrot_ang(orient)  # vrot in cartesians, theta in radians
		a, b, c = vrot
		quatw, quati, quatj, quatk = cos(theta / 2), a * sin(theta / 2), b * sin(theta / 2), c * sin(theta / 2) # these is the quaternion you need to put in the LAMMPS input file
		quaternions=vstack((quaternions,[quatw, quati, quatj, quatk]))
	return quaternions

def uniform_sphere_points(N, radius):
	#Generates uniformly distributed points on the surface of a sphere of given radius.
	#Note that if the area per point ~pi*r^2, the radius for N points is R~0.5*r*sqrt(N).
	points = zeros((N, 3))
	indices = arange(0, N, dtype=float) + 0.5
	phi = arccos(1 - 2*indices/N)
	theta = pi * (1 + 5**0.5) * indices
	x = radius * cos(theta) * sin(phi)
	y = radius * sin(theta) * sin(phi)
	z = radius * cos(phi)
	points[:, 0] = x
	points[:, 1] = y
	points[:, 2] = z
	current_distance = average_nearest_neighbor_distance(points)
	return points,current_distance

natoms_sphere=6750
poly_diameter=2.0

ellipsoidflag_polymer=0
ellipsoidflag_sphere=1
dist_neighbors=0.5 #average distance between neighboring points on the sphere (value should be a bit smaller than the desired one)
density=6/pi 	#density must be set to 6/pi for particles of mass 1 and diameter 1

shift=1e-8 	#a small shift of coordinates is needed or LAMMPS will complain that the rigid bodies are the same particle

#GENERATE MESH
radius=0.5*dist_neighbors*sqrt(natoms_sphere)
box_new_x=6*radius
box_new_y=4*radius
box_new_z=4*radius
vertices, achieved_distance = uniform_sphere_points(natoms_sphere,radius)
print("N = %d, R = %.2f, dist = %.3f"%(natoms_sphere,radius,achieved_distance))

vertices[:,0]=vertices[:,0]+0.5*box_new_x
vertices[:,1]=vertices[:,1]+0.5*box_new_y
vertices[:,2]=vertices[:,2]+0.5*box_new_z

#GENERATE QUATERNIONS
quaternions_vertices=generate_quaternions(vertices,box_new_x,box_new_y,box_new_z)

fname='data.compress'

print(fname)
fname_noext=os.path.splitext(fname)[0]
with open('%s'%fname,"r") as f:
	lines=[l.rstrip().split() for l in f.readlines()]
natoms=int(lines[2][0])
nbonds_init=int(lines[4][0])
box_x=float(lines[7][1])
box_y=float(lines[8][1])
box_z=float(lines[9][1])
natoms_tot=2*natoms+natoms_sphere

atom_data=lines[26:26+natoms]
bonds_init=array(lines[32+2*natoms:32+2*natoms+nbonds_init],dtype='int32')
bonds=bonds_init[bonds_init[:,1]==1] #only keep type 1 bonds
coors=array(atom_data,dtype='float32')[:,3:6]
nbonds=len(bonds)

#Shift coordinates due to box rescaling and upscale them
coors[:,0]=(coors[:,0]-0.5*box_x)*poly_diameter+0.5*box_new_x
coors[:,1]=(coors[:,1]-0.5*box_y)*poly_diameter+0.5*box_new_y
coors[:,2]=(coors[:,2]-0.5*box_z)*poly_diameter+0.5*box_new_z

with open('%s_mesh_ylz_patchy_poly_diameter%.2f_natoms_sphere%d_lx%.2f_ly%.2f_lz%.2f.lammpsdata'%(fname_noext,poly_diameter,natoms_sphere,box_new_x,box_new_y,box_new_z),"w") as fout:
	fout.write("LAMMPS data file\n\n%d atoms\n3 atom types\n%d ellipsoids\n%d bonds\n1 bond types\n\n0.0 %f xlo xhi\n0.0 %f ylo yhi\n0.0 %f zlo zhi\n\nAtoms\n\n"%(natoms_tot,natoms_sphere,nbonds,box_new_x,box_new_y,box_new_z))

	#NB: assuming atom_style hybrid ellispoid bond molecular in LAMMPS input file
	#The format of the Atoms section must be:
	#atom-ID atom-type x y z ellipsoidflag density molecule-ID
	#ellipsoidflag=1 for ylz particles, 0 for point particles
	#density must be set to 6/pi for particles of mass 1 and diameter 1

	#mol_id_all=[]
	for k,l in enumerate(atom_data):
		atom_id=int(l[0])
		#poly_id=int(l[1])
		atom_type=int(l[2])
		x=coors[k][0]
		y=coors[k][1]
		z=coors[k][2]
		xflag=0
		yflag=0
		zflag=0

		fout.write("%d %d %.10f %.10f %.10f %d %.6f %d %d %d %d\n"%(atom_id,atom_type,x,y,z,ellipsoidflag_polymer,density,atom_id,xflag,yflag,zflag))
		fout.write("%d 3 %.10f %.10f %.10f %d %.6f %d %d %d %d\n"%(atom_id+natoms,x+rand()*shift,y+rand()*shift,z+rand()*shift,ellipsoidflag_polymer,density,atom_id,xflag,yflag,zflag))

	mol_id_poly_max=2*natoms

	for k,r in enumerate(vertices):
		atom_id=2*natoms+k+1
		#atom_type=2
		#mol_id=mol_id_poly_max+1
		x=r[0]
		y=r[1]
		z=r[2]
		xflag=0
		yflag=0
		zflag=0

		fout.write("%d 2 %.8f %.8f %.8f %d %.6f %d %d %d %d\n"%(atom_id,x,y,z,ellipsoidflag_sphere,density,mol_id_poly_max+1,xflag,yflag,zflag))

	fout.write("\nEllipsoids\n\n")
	#line syntax: atom-ID shapex shapey shapez quatw quati quatj quatk
	#shapex,shapey,shapez = 3 diameters of ellipsoid (distance units)
	#quatw,quati,quatj,quatk = quaternion components for orientation of atom

	for k, q in enumerate(quaternions_vertices):
		atom_id=2*natoms+k+1
		fout.write("%d 0.999999 1.000000 1.000000 %.6f %.6f %.6f %.6f\n"%(atom_id,q[0],q[1],q[2],q[3]))

	fout.write("\nBonds\n\n")

	for r in bonds:
		fout.write("%d 1 %d %d\n"%(r[0],r[2],r[3]))

