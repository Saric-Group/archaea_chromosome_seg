#Takes a configuration of patchy polymers + YLZ membrane, and rotates it so that the normal to the plane 
#obtained from the LDA analysis becomes the new x axis.

def compute_cm(coors):
	xcm=average(coors[:,0])
	ycm=average(coors[:,1])
	zcm=average(coors[:,2])
	return array([xcm,ycm,zcm])

def find_perpendicular_vectors(v):
    """Given a unit vector v, find two perpendicular unit vectors."""
    v = v / np.linalg.norm(v)  # Ensure v is normalized
    # Choose an arbitrary vector not parallel to v
    if np.allclose(v, [0, 0, 1]):
        temp = np.array([0, 1, 0])  # Use y-axis if v is close to z-axis
    else:
        temp = np.array([0, 0, 1])  # Default choice is the z-axis
    # Compute first perpendicular vector
    v1 = np.cross(v, temp)
    v1 /= np.linalg.norm(v1)  # Normalize
    # Compute second perpendicular vector
    v2 = np.cross(v, v1)
    v2 /= np.linalg.norm(v2)  # Normalize
    return v1, v2

############################
# Coordinates transformation
############################

def change_basis_3d(points, u, v, w):
    """
    Transforms a set of 3D points to a new coordinate system defined by (u, v, w).
    
    Parameters:
    points (numpy array): Nx3 array of points in the original basis
    u, v, w (numpy arrays): Basis vectors defining the new coordinate system
    
    Returns:
    transformed_points: Nx3 array of points in the new coordinate system
    """
    # Construct the transformation matrix (each column is a basis vector)
    R = np.column_stack((u, v, w))  # Equivalent to np.array([u, v, w]).T
    
    # Transform points to the new basis
    transformed_points = np.dot(points, np.linalg.inv(R).T)  # Compute P' = R^T * P
    return transformed_points

############################
# Quaternions stuff
############################

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

#=================================
#Ivan's code for quaternions
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

import os, re, sys, time
from numpy import *
from numpy.random import *
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.svm import LinearSVC

input_file=input("Enter input file name (paths to the folders containing the dump folders):")
fname='data.ylz_patchy_hc_mobile'

polymer_length=500
print('Using polymer length = %d'%polymer_length)

shift=1e-8 	#a small shift of coordinates is needed or LAMMPS will complain that the rigid bodies are the same particle

with open(input_file,'r') as f:
	paths=[l.rstrip() for l in f.readlines()]

for path in paths:
	print (path)

	lda_norm_vec=loadtxt("%s/last_lda_normal.dat"%path)
	perpendicular_vectors=find_perpendicular_vectors(lda_norm_vec)
	orthonormal_basis=[lda_norm_vec,perpendicular_vectors[0],perpendicular_vectors[1]]

	with open('%s/%s'%(path,fname),"r") as fin:
		lines=[l.rstrip().split() for l in fin.readlines()]
	natoms=int(lines[2][0])
	nbonds=int(lines[4][0])
	box_x=float(lines[8][1])
	box_y=float(lines[9][1])
	box_z=float(lines[10][1])

	#The data file format should be:
	#Atoms section:
	#atom-ID atom-type x y z ellipsoidflag density molecule-ID xflag yflag zflag
	#Ellipsoids section:
	#atom-ID shapex shapey shapez quatw quati quatj quatk

	atoms_tmp=array(lines[24:24+natoms],dtype='float32')
	bonds=array(lines[30+2*natoms:30+2*natoms+nbonds],dtype='int32')
	ellipsoids_tmp=array(lines[33+2*natoms+nbonds:],dtype='float32')
	nmem=len(ellipsoids_tmp)

	atoms=atoms_tmp[argsort(atoms_tmp[:,0])] #sort by atom id
	ellipsoids=ellipsoids_tmp[argsort(ellipsoids_tmp[:,0])] #sort by atom id

	atom_ids=atoms[:,0].astype(int)
	atom_types=atoms[:,1].astype(int)
	coors=atoms[:,2:5]-compute_cm(atoms[:,2:5]) 
	ellipsoidflags=atoms[:,5].astype(int)
	densities=atoms[:,6]
	mol_ids=atoms[:,7].astype(int)
	xyz_flags=atoms[:,8:]

	rotated_coors=change_basis_3d(coors,orthonormal_basis[0],orthonormal_basis[1],orthonormal_basis[2])
	rotated_coors=rotated_coors+0.5*array([box_x,box_y,box_z]) #Recenter in middle of the box
	rotated_coors=rotated_coors%array([box_x,box_y,box_z]) #Wrap inside box (for beads evaporated from membrane)

	################################################################################################
	#Use LDA to verify that the normal to the LDA plane is indeed the x axis

	#Using LDA: https://scikit-learn.org/0.16/modules/generated/sklearn.lda.LDA.html
	poly_coors=rotated_coors[atom_types==1]
	mol_ids_poly=mol_ids[atom_types==1]
	poly_ids=where(mol_ids_poly <= polymer_length, 1, 2)
	lda=LinearDiscriminantAnalysis(n_components=1)
	lda=lda.fit(poly_coors, poly_ids) #Fit LDA model according to the given training data and parameters.
	predicted=lda.predict(poly_coors) #Predict class labels for samples in X.

	#Write normal to separation plane
	lda_coef = lda.coef_[0]
	lda_intercept = lda.intercept_
	lda_norm_vec_new_nonorm=array([lda_coef[0],lda_coef[1],lda_coef[2]])
	lda_norm_vec_new=lda_norm_vec_new_nonorm/linalg.norm(lda_norm_vec_new_nonorm)
	#print(lda_norm_vec)
	print(lda_norm_vec_new)
	################################################################################################

	ellips_atom_ids=ellipsoids[:,0]
	ellips_shapes=ellipsoids[:,1:4]
	ellips_quaternions=ellipsoids[:,4:]
	rotated_membrane_coors=rotated_coors[atom_types==2]

	#Generate new quaternions for rotated vescicle
	new_quaternions=generate_quaternions(rotated_membrane_coors,box_x,box_y,box_z)

	#Structure of atoms section:
	#IDs from 1 to 2*polymer_length (1000) are polymer coordinates
	#IDs from 2*polymer_length+1 (1001) to 4*polymer_length (2000) are patch coordinates
	#IDs from 4*polymer_length+1 onwards are membrane coordinates

	with open('%s/%s_rotate_xaxis_lda_norm'%(path,fname),"w") as fout:
		fout.write("LAMMPS data file\n\n%d atoms\n3 atom types\n%d ellipsoids\n%d bonds\n1 bond types\n\n0.0 %f xlo xhi\n0.0 %f ylo yhi\n0.0 %f zlo zhi\n\nAtoms\n\n"%(natoms,len(ellipsoids),nbonds,box_x,box_y,box_z))

		#NB: assuming atom_style hybrid ellispoid bond molecular in LAMMPS input file
		#The format of the Atoms section must be:
		#atom-ID atom-type x y z ellipsoidflag density molecule-ID
		#ellipsoidflag=1 for ylz particles, 0 for point particles
		#density must be set to 6/pi for particles of mass 1 and diameter 1

		for n in range(len(atoms)):
			atom_id=atom_ids[n]
			atom_type=atom_types[n]
			x_rot=rotated_coors[n][0]
			y_rot=rotated_coors[n][1]
			z_rot=rotated_coors[n][2]
			ellipsoidflag=ellipsoidflags[n]
			density=densities[n]
			mol_id=mol_ids[n]
			xflag=0
			yflag=0
			zflag=0

			if(atom_type==3):
				fout.write("%d %d %.10f %.10f %.10f %d %.6f %d %d %d %d\n"%(atom_id,atom_type,x_rot+rand()*shift,y_rot+rand()*shift,z_rot+rand()*shift,ellipsoidflag,density,mol_id,xflag,yflag,zflag))
			else:
				fout.write("%d %d %.10f %.10f %.10f %d %.6f %d %d %d %d\n"%(atom_id,atom_type,x_rot,y_rot,z_rot,ellipsoidflag,density,mol_id,xflag,yflag,zflag))

		fout.write("\nBonds\n\n")

		for r in bonds:
			fout.write("%d 1 %d %d\n"%(r[0],r[2],r[3]))

		fout.write("\nEllipsoids\n\n")
		#line syntax: atom-ID shapex shapey shapez quatw quati quatj quatk
		#shapex,shapey,shapez = 3 diameters of ellipsoid (distance units)
		#quatw,quati,quatj,quatk = quaternion components for orientation of atom

		for k in range(len(ellipsoids)):
			atom_id=ellips_atom_ids[k]
			shapex=ellips_shapes[k][0]
			shapey=ellips_shapes[k][1]
			shapez=ellips_shapes[k][2]
			qx=new_quaternions[k][0]
			qy=new_quaternions[k][1]
			qz=new_quaternions[k][2]
			qw=new_quaternions[k][3]
			fout.write("%d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n"%(atom_id,shapex,shapey,shapez,qx,qy,qz,qw))



