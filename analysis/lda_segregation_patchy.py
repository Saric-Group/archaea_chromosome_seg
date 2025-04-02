#This script performs LDA analysis on the coordinates of the two polymers and returns error fraction f.
#The segregation efficiency can be obtaind as s=1-2*f, as detailed in the Methods (see Supporting Information).
#It returns:
#-lda_frac_errors_vs_time.dat 	Error fraction f vs time
#-lda_normal_vs_time.dat 		Normal to the LDA plane (segregation axis) vs time (3 components, x,y,z)
#-last_lda_normal.dat  			Last LDA normal of the simulation

import os, sys, time
from numpy import *
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.svm import LinearSVC

folder='dumplin'
polymer_length=500
print('Using polymer length = %d'%polymer_length)

if((folder!='dump') and (folder!='dumplin')):
	print('Folder does not exist (%s)'%folder)
	sys.exit()

path='.'
files=[f for f in os.listdir("%s/%s"%(path,folder)) if ("lammpstrj" in f)]

ferrors_vs_time=array([0,0])
lda_normal_vs_time=array([0,0,0,0])

for fin in files:
	with open("%s/%s/%s"%(path,folder,fin),'r') as myfile:
		head=[next(myfile) for i in range(9)]
	my_format=head[8]
	time=int(head[1].split()[0]) 

	data_tmp=loadtxt("%s/%s/%s"%(path,folder,fin),skiprows=9) 
	data=data_tmp[data_tmp[:,0]==1] #only take polymers coors 

	#mol_ids = loadtxt(target_file, usecols=1, dtype=int16)
	#coordinates = loadtxt(target_file, usecols=(2, 3, 4))
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

	natoms=len(mol_ids)

	#Using LDA: https://scikit-learn.org/0.16/modules/generated/sklearn.lda.LDA.html
	lda=LinearDiscriminantAnalysis(n_components=1)
	lda=lda.fit(coordinates, poly_ids) #Fit LDA model according to the given training data and parameters.
	predicted=lda.predict(coordinates) #Predict class labels for samples in X.
	ferrors=count_nonzero(predicted != poly_ids)/natoms
	ferrors_vs_time=vstack((ferrors_vs_time,array([time,ferrors])))
	#print('time=%d f_errors=%f'%(time,ferrors))

	#Write normal to separation plane 
	lda_normal_nonorm = lda.coef_.flatten()
	lda_normal = lda_normal_nonorm/linalg.norm(lda_normal_nonorm) 
	lda_normal_vs_time=vstack((lda_normal_vs_time,array([time,lda_normal[0],lda_normal[1],lda_normal[2]])))


ferrors_vs_time_sorted=ferrors_vs_time[ferrors_vs_time[:,0].argsort()][1:]
savetxt("%s/lda_frac_errors_vs_time.dat"%path,ferrors_vs_time_sorted,fmt='%d %.3e')

lda_normal_vs_time_sorted=lda_normal_vs_time[lda_normal_vs_time[:,0].argsort()][1:]
savetxt("%s/lda_normal_vs_time.dat"%path,lda_normal_vs_time_sorted,fmt='%d %.3f %.3f %.3f')

last_lda_normal=lda_normal_vs_time_sorted[-1][1:]
savetxt("%s/last_lda_normal.dat"%path,last_lda_normal)
#print(f"Number of errors = {count_nonzero(predict != poly_ids)} out of {len(poly_ids)}")



