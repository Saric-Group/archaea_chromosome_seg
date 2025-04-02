#Creates 2 rings which are initially connected to each other (like a circular stair...)

from numpy import *
from numpy.random import *

def rotate2d(vector,angle): #Rotate vector by angle
	r=array([cos(angle),-sin(angle),sin(angle),cos(angle)])
	rot=reshape(r,(2,2)) #2d rotation matrix
	return rot.dot(vector)

bond_length=0.8

ring_size=500
theta=2*pi/ring_size
radius=bond_length/(2*sin(theta/2))

box_x=1.1*(radius*2)
box_y=1.1*(radius*2)
box_z=1.1*(radius*2)

nrings=2
natoms=nrings*ring_size
nbonds=(1+nrings)*ring_size

vec=array([1,0])

with open('2rings_connected_N%d_box%.2f.lammpsdata'%(ring_size,box_x),'w') as fout:
	fout.write("LAMMPS data file\n\n%d atoms\n1 atom types\n%d bonds\n2 bond types\n\n0.0 %f xlo xhi\n0.0 %f ylo yhi\n0.0 %f zlo zhi\n\nMasses\n\n1 1\n\nAtoms\n\n"%(natoms,nbonds,box_x,box_y,box_z))
	
	atom_id=0
	mol_id=0
	bonds_intra=[]
	bonds_inter=[]

	ix0=box_x*0.5
	iy0=box_y*0.5

	iz0=0.5*(box_z-1)

	mol_id+=1
	r0=array([ix0,iy0,iz0])

	for n in range(ring_size):
		atom_id+=1
		rotvec=rotate2d(vec,n*theta)
		r=r0+radius*array([rotvec[0],rotvec[1],0])
		fout.write("%d %d 1 %.15f %.15f %.15f\n"%(atom_id,mol_id,r[0],r[1],r[2]))
		if(n>0):
			bonds_intra.append((atom_id,atom_id-1))
	bonds_intra.append((atom_id,atom_id-ring_size+1))

	iz0=0.5*(box_z+bond_length)
	
	mol_id+=1
	r0=array([ix0,iy0,iz0])

	for n in range(ring_size):
		atom_id+=1
		rotvec=rotate2d(vec,n*theta)
		r=r0+radius*array([rotvec[0],rotvec[1],0])
		fout.write("%d %d 1 %.15f %.15f %.15f\n"%(atom_id,mol_id,r[0],r[1],r[2]))
		bonds_inter.append((atom_id,atom_id-ring_size))
		if(n>0):
			bonds_intra.append((atom_id,atom_id-1))
	bonds_intra.append((atom_id,atom_id-ring_size+1))

	if((len(bonds_intra)+len(bonds_inter))!=nbonds):
		print("WARNING: Number of bonds different from expected number of bonds (have %d, expect %d)"%(len(bonds),nbonds))

	#Print bonds
	fout.write("\nBonds\n\n")
	for nb, bond in enumerate(bonds_intra):
		fout.write("%d 1 %d %d\n"%(nb+1,bond[0],bond[1]))
	for nb, bond in enumerate(bonds_inter):
		fout.write("%d 2 %d %d\n"%(nb+1,bond[0],bond[1]))