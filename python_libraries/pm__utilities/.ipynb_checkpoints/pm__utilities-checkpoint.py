import numpy as np
import math
import copy
import scipy.optimize as opti
#import transformations as tr

def getTetrahedronVolume(x,y,z):
	"""
	Returns volume of a tetrahedron defined by vectors x,y,z
	"""
	return abs(1/6.*np.dot(x,np.cross(y,z)))

def getParallelepipedVolume(x,y,z):
	"""
	Returns volume of a parallelepiped
	"""
	return abs(np.dot(x,np.cross(y,z)))


def angle(v1, v2):
	"""
	Returns angle between vectors v1 and v2 (in radians)
	"""
	cosang = np.dot(v1, v2)
	sinang = np.linalg.norm(np.cross(v1, v2))
	return 180./math.pi*np.arctan2(sinang, cosang)
	
def toFractional(pos_carth,a,b,c):
	"""
	Converts carthesian coordinates in the form [ax, ay, az] or set of carthesian coordinates
	in the form [[ax, ay, az], [bx, by, bz],..., [cx, cy, cz]] into fractional coordinates with 
	respect to lattice vectors a, b, c """

	latvec_vertical=np.array([a,b,c]).transpose()	#vertically ordered lattice vectors
	#print "toFractional: ","************"
	#print "toFractional: ",a,b,c,latvec_vertical
	#print "toFractional: ",pos_carth
	#print "toFractional: ","************"
	res=[]
	for ii in range(len(pos_carth)):
		res.append(np.linalg.solve(latvec_vertical,np.array(pos_carth[ii]).transpose()).tolist())
	return copy.deepcopy(np.array(res))

#	return copy.deepcopy(np.linalg.solve(latvec_vertical,pos_carth.transpose()).tolist())

def toCarthesian(pos_frac,a,b,c):
	return copy.deepcopy(np.dot(pos_frac,[a,b,c]).tolist())

def chop(array,threshold=1e-10):
	"""
	Set all values in the array, which are smaller than the threshold, to zero.
	"""
	shape=np.shape(array)
#	print "chop function:", shape
	aux_array=(array.reshape(1,-1))[0]
	for ii in range(len(aux_array)): 
		if abs(aux_array[ii])<threshold:
			aux_array[ii]=0.

#	print "chop function (before reshape):",aux_array
	#aux_array.reshape(shape)
	aux_array=(aux_array.reshape(shape))
#	print "chop function (after reshape):",aux_array
	return aux_array

def projectToDirection(vector, direction=np.array([1,0,0])):
	"""Returns projection of n-dimensional vector to given n-dimensional direction.
	Returned array is numpy.array
	Both vector and direction can be given as simple n-dimensional list of numbers
	[a,b,c, ... , ] or as N-dimensional set of triplets [[a,b,c],[d,e,f], ... ],
	and can even differ in shape, but must have same number of compoents"""
	
	#First reshape the direction to array 1xN
	direction_resh=np.reshape(direction,(1,-1))[0]
	#Reshape also the vector to 1xN
	vector_resh=np.reshape(vector,(1,-1))[0]
	if len(direction_resh) != len(vector_resh):
		print("Projection: projected vector and direction differ in number of components, returning 0")
		return 0.
	#get norm of direction
	norm_direction=np.linalg.norm(direction_resh)
	
	return (1./norm_direction**2)*np.dot(vector_resh,direction_resh)

#def projectToSubspace(vectors,subspace,perpendicular='no')
#	"""Project vector or set of vectors to subspace """
#	space_dim=1	#dimensionality of the space in which we do the projection
#
#	#Record original shape
#	vector_shape=np.shape(vector)
#
#	#Derive dimensionality of the vector space
#	if len(vector_shape)>3: print "Projection: maximum allowed dimension is 3, terminating..."; exit(1)
#	if len(vector_shape) == 3:
#		#interpret of set of vectors, each of which is composed of "subvectors" (e.g. displacements associated with individual atoms)
#		print "Projection: converting multiple vectors"
#		space_dim=vector.shape[1]*vector.shape[2]
#	else if len(vector_shape) == 2:
#		#interpret as one vector, each of which is composed of "subvectors" (e.g. displacements associated with individual atoms 
#		space_dim=vector_shape[0]*vector.shape[1]
#	else if len(vector_shape) == 1:
#		#interpret as one vecto
#		space_dim=vector_shape[0]
#	print "Projection: vector_shape = ",vector_shape, "dimension of space = ", space_dim
#
#	#Reshape vectorst to standard form
#	vector_resh=np.reshape(vector,(-1,space_dim))
#	if len(vector_resh[-1]) != space_dim: print "Projection: projected vector has wrong number of components, terminating..."; exit(1)
#	print "Projection: vector_resh = ",vector_resh
#
#	direction_resh=np.reshape(direction,(-1,space_dim))
#	if len(direction_resh[-1]) != space_dim: print "Projection: projection direction has wrong number of components, terminating..."; exit(1)
#	print "Projection: direction_resh = ",direction_resh
#
#	#Test orthogonality of the subspace base vectors. If not mutually perpendicular, issue warning
#	
#	#Perform projection
#
#	#If perpendicular == 'yes', subtract projction from original vector (projecting to space perpendicular to those given by parameter subspace
#
#	#Recover original form of vectors
#

def fitCoordinateSystem(fit_data,euler_angles_init=(0.001,0.001,0.001)):
	"""
	Finds the coordinate system, which approximates best 
	vectors in the fit_data (n*3) array. These vectors are supposed to be perturbed direction x,y,z of the coordinate system.
	Vectors in fit_data are normalized before fitting.
	Note: all vectors are taken with the same weight (sigma)\n
	:parameter fit_data: matrix of vectors to which the coordinate system will be fitted.\n
	:rtype: [coordinate system vectors, transformation matrix, vector of Eulre angles]
	"""

	def aux_errorFunctionForCoordinateSystemFitting(arg,alpha,beta,gamma):
		""" Function accepts (array) of vectors and evaluates array of distances from vectors of a coordinate system characterized by the Euler angles alpha, beta, gamma """
		
		carth_x=np.array([1.,0.,0.])
		carth_y=np.array([0.,1.,0.])
		carth_z=np.array([0.,0.,1.])

		argT=arg.transpose()
		argT_normalized=copy.deepcopy(argT)
		
		""" Normalize vectors to which is the coordinate system fitted """
		for ii in range(len(argT)):
			argT_normalized[ii,:]/=np.linalg.norm(argT[ii,:])

		R=tr.euler_matrix(alpha,beta,gamma,'rxyz')[:3,:3]

		""" Transform carthesian coordinates according to alpha,beta,gamma """
		rot_carth_x=np.dot(R,carth_x)
		rot_carth_y=np.dot(R,carth_y)
		rot_carth_z=np.dot(R,carth_z)

		error=[]
		
		""" Evaluate distance between rot_carth_x, rot_carth_y, rot_carth_z and +- vector from arg and save the smallest distance """
		for ii in range(len(argT)):
			totalmin=1000
			actmin=np.linalg.norm(   argT_normalized[ii]-rot_carth_x)
			if actmin<totalmin:totalmin=actmin
			actmin=np.linalg.norm(   argT_normalized[ii]-rot_carth_y);
			if actmin<totalmin:totalmin=actmin
			actmin=np.linalg.norm(   argT_normalized[ii]-rot_carth_z);
			if actmin<totalmin:totalmin=actmin

			error.append(totalmin)

		return np.array(error)
	
	init    = np.array(euler_angles_init)		#initial_values for the fit
	sigma	= np.ones_like(range(len(fit_data)))	#Weights are the same for all vectors
	target	= np.zeros_like(sigma)			#Target difference between fitted coordinate system and the provided vectors is zero

	res=opti.curve_fit(aux_errorFunctionForCoordinateSystemFitting,fit_data.transpose(),target,init,sigma)

	print( "Euler angles of the coordinate system are: (%.2f,%.2f,%.2f) deg."%tuple(180/math.pi*res[0]) )
	print( "Coordinate system is defined bythe basis:" )
	
	R=tr.euler_matrix(res[0][0],res[0][1],res[0][2],'rxyz')[:3,:3]
	print( "x'="+str(np.dot(R,[1,0,0]) ))
	print( "y'="+str(np.dot(R,[0,1,0]) ))
	print( "z'="+str(np.dot(R,[0,0,1]) ))
	
	ret_coords=np.array([np.dot(R,[1,0,0]),np.dot(R,[0,1,0]),np.dot(R,[0,0,1])])
	ret_transf_matrix=R
	ret_angles=res[0]

	print( "c"+str(ret_coords ))
	print( "t"+str(ret_transf_matrix ))
	print( "a"+str(ret_angles ))

	return [ret_coords,ret_transf_matrix,ret_angles]

def fracToCarthesianFromAngles(a_len, b_len, c_len, alpha, beta, gamma, frac1, frac2, frac3):
	carthx=a_len*frac1+b_len*math.cos(math.pi*gamma/180.)*frac2+c_len*math.cos(math.pi*beta/180.)*frac3
	carthy=b_len*math.sin(math.pi*gamma/180.)*frac2+c_len*(math.cos(math.pi*alpha/180.)-math.cos(math.pi*beta/180.)*math.cos(math.pi*gamma/180.))/math.sin(math.pi*gamma/180.)*frac3
	carthz=c_len*math.sqrt(1.-(math.cos(math.pi*alpha/180.))*(math.cos(math.pi*alpha/180.))-(math.cos(math.pi*beta/180.))*(math.cos(math.pi*beta/180.))-(math.cos(math.pi*gamma/180.))*(math.cos(math.pi*gamma/180.))+2*math.cos(math.pi*alpha/180.)*math.cos(math.pi*beta/180.)*math.cos(math.pi*gamma/180.))/math.sin(math.pi*gamma/180.)*frac3
	
	return[carthx,carthy,carthz]
