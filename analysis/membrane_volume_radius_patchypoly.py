#This scripts analyzes the membrane configurations. It approximates the membrane surface using the convex hull method, and returns several quantities:
#-natoms.dat 					Number of membrane atoms
#-hull_radius_av.dat  			Mean membrane radius
#-hull_volume_vs_time.dat  		Membrane volume vs time
#-hull_area_vs_time.dat  		Membrane area vs time
#-hull_radius_vs_time.dat  		Membrane radius vs time 
#-neigh_distance_vs_time.dat 	Mean distance between neighbouring membrane particles vs time
#-area_per_bead_vs_time.dat 		Mean area per membrane bead vs time (estimate)


def average_nearest_neighbor_distance(points):
    tree = KDTree(points)
    distances, _ = tree.query(points, k=2)  # k=2 because the nearest point includes the point itself
    return mean(distances[:, 1])  # Skip distance to itself which is 0

import os, sys, time
from numpy import * 
#import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from scipy.spatial import ConvexHull

folder='dumplin'
tmin_av_radius=3e7
ptype=2

print("\nUsing tmin_av_radius=%e\n"%tmin_av_radius)

path='.'
files=[f for f in os.listdir(path+'/'+folder) if ("lammpstrj" in f)]

hull_volume_vs_time=empty((0,2))
hull_area_vs_time=empty((0,2))
hull_radius_vs_time=empty((0,2))
neigh_distance_vs_time=empty((0,2))
area_per_bead_vs_time=empty((0,2))
hull_radius_av=0
n_radius_values=0
for n,fname in enumerate(files):
	with open('%s/%s/%s'%(path,folder,fname),'r') as myfile:
		head=[next(myfile) for i in range(8)]
	time=int(head[1].split()[0])
	xmin=float(head[5].split()[0])
	xmax=float(head[5].split()[1])
	ymin=float(head[6].split()[0])
	ymax=float(head[6].split()[1])
	zmin=float(head[7].split()[0])
	zmax=float(head[7].split()[1])
	box_x=xmax-xmin
	box_y=ymax-ymin
	box_z=zmax-zmin

	data=loadtxt("%s/%s/%s"%(path,folder,fname),skiprows=9)
	coors=data[data[:,0]==ptype][:,3:6]
	#Sometimes, particles fly away from the membrane; we try to neglect them.
	coors=coors[coors[:,0]<box_x]
	coors=coors[coors[:,1]<box_y]
	coors=coors[coors[:,2]<box_z]
	natoms=len(coors)
	neigh_distance=average_nearest_neighbor_distance(coors)
	area_per_bead=pi*(0.5*neigh_distance)**2 
	#NB: this tentative estimate for the area per bead does not correspond to the one we get by plotting
	#R=(N*a/4*pi)^(1/2), with R = membrane radius, N = n. of particles, a = area per bead.
	hull=ConvexHull(coors)
	hull_volume=hull.volume
	hull_area=hull.area
	hull_radius=(3*hull_volume/(4*pi))**(1/3.)
	if(time>tmin_av_radius):
		hull_radius_av+=hull_radius
		n_radius_values+=1
	hull_volume_vs_time=vstack((hull_volume_vs_time,array([time,hull_volume])))
	hull_area_vs_time=vstack((hull_area_vs_time,array([time,hull_area])))
	hull_radius_vs_time=vstack((hull_radius_vs_time,array([time,hull_radius])))
	neigh_distance_vs_time=vstack((neigh_distance_vs_time,array([time,neigh_distance])))
	area_per_bead_vs_time=vstack((area_per_bead_vs_time,array([time,area_per_bead])))

hull_radius_av/=n_radius_values

hull_volume_vs_time=hull_volume_vs_time[argsort(hull_volume_vs_time[:,0])]
hull_area_vs_time=hull_area_vs_time[argsort(hull_area_vs_time[:,0])]
hull_radius_vs_time=hull_radius_vs_time[argsort(hull_radius_vs_time[:,0])]
neigh_distance_vs_time=neigh_distance_vs_time[argsort(neigh_distance_vs_time[:,0])]
area_per_bead_vs_time=area_per_bead_vs_time[argsort(area_per_bead_vs_time[:,0])]

savetxt('%s/natoms.dat'%path,array([natoms]))
savetxt('%s/hull_radius_av.dat'%path,array([hull_radius_av]))
savetxt('%s/hull_volume_vs_time.dat'%path,hull_volume_vs_time)
savetxt('%s/hull_area_vs_time.dat'%path,hull_area_vs_time)
savetxt('%s/hull_radius_vs_time.dat'%path,hull_radius_vs_time)
savetxt('%s/neigh_distance_vs_time.dat'%path,neigh_distance_vs_time)
savetxt('%s/area_per_bead_vs_time.dat'%path,area_per_bead_vs_time)


