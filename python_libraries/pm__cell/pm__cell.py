#!/usr/bin/env python

import	sys
import	numpy as np
import	re
import	os
import	shutil
import	glob
import	copy
import	math
import	pm__constants as pmco
import	pm__utilities as pmut
import	pm__data_from_file as pmdf
#import	pm__latex_aid as pmla
import	logging as log
import time 

"""
History
5.11.2019: Procedure setLatticeVectors(self,abc) was altered. When the lattice vectors are newly set, the fractional coordinates of atoms are recalculated with new vectors
20.5.2019: changes in changeElementNames(self,element_names): given atom names can be interpretted as species names or list of names of all atoms in the system
"""

class Atom:
	"""
	Class for storing information about individual atom
	"""
	def __init__(self):
		self.position_frac=[0.,0.,0.,]
		self.name=""
		self.charge=0.
		self.coreshell="core"
		self.wyckoff=""
		self.spin=0.
	

	def __str__(self):
		aux_string= "Information about atom:\n   Name = %s\n   Charge = %.3f\n   Coreshell = %s\n   Wyckoff = %s\n   Spin = %.3f\n   Coordinates (fractional) = (%.3f,%.3f,%.3f)"%(
					self.name,
					self.charge,
					self.coreshell,
					self.wyckoff,
					self.spin,
					self.position_frac[0],
					self.position_frac[1],
					self.position_frac[2],
					)
		return aux_string

	#def __eq__(self, other):
	#	if isinstance(other, Atom):
	#		return self.a == other.a and self.b == other.b
	#	return False

	
class Lattice:
	"""
	Class for storing information about crystallographic lattice
	"""
	def __init__(self):
		self.a=[1,0,0]		#cell vector a
		self.b=[0,1,0]		#cell vector b
		self.c=[0,0,1]		#cell vector c
		
		self.alpha=0.		#angle between b and c
		self.beta=0.		#angle between a and c
		self.gamma=0.		#angle between a and b
		self.mult=1.
	
	def __str__(self):
		self.getAnglesFromVectors()
		str_string="Lattice vectors:\n"
		str_string=str_string + \
				"     a: % .6f % .6f % .6f\n"   % (self.a[0],self.a[1],self.a[2]) + \
				"     b: % .6f % .6f % .6f\n"   % (self.b[0],self.b[1],self.b[2]) + \
				"     c: % .6f % .6f % .6f\n"   % (self.c[0],self.c[1],self.c[2]) + \
				"     alpha=%10.6f beta=%10.6f gamma=%10.6f\n" % (self.alpha,self.beta,self.gamma) + \
				"     |a|  =%10.6f |b| =%10.6f |c|  =%10.6f, V=%10.6f"   % (np.linalg.norm(self.a),np.linalg.norm(self.b),np.linalg.norm(self.c),np.dot(self.c,np.cross(self.a,self.b)))
		return str_string
	
	def getAnglesFromVectors(self):
		self.alpha=np.arccos(np.dot(self.b,self.c)/(np.linalg.norm(self.b)*np.linalg.norm(self.c)))*180./math.pi
		self.beta =np.arccos(np.dot(self.a,self.c)/(np.linalg.norm(self.a)*np.linalg.norm(self.c)))*180./math.pi
		self.gamma=np.arccos(np.dot(self.a,self.b)/(np.linalg.norm(self.a)*np.linalg.norm(self.b)))*180./math.pi
	
	def getLatticeArray(self): return np.array([self.a,self.b,self.c])
	def getLatticeArrayCol(self): return np.array([self.a,self.b,self.c]).T
	def getLattRecArrayCol(self): 
		aux_V		=abs(np.dot(self.a,np.cross(self.b,self.c)))
		aux_factor	=2*math.pi/aux_V
		aux_b1=		np.cross(self.b,self.c)
		aux_b2=		np.cross(self.c,self.a)
		aux_b3=		np.cross(self.a,self.b)
		return aux_factor*np.transpose(np.array([aux_b1,aux_b2,aux_b3]))

class Cell:
	"""
	Class for storing information about an atomic unit cell\n
	>>> import pm__cell as pmc
	>>> cell=pmc.Cell("POSCAR")
	"""
	

#	shell_models = {
#    "PTO_shimada": {
#        "Pb": {
#            	'core': {'mass': 207.2, 'charge': 5.49595},
#            	'shell': {'mass': 0.0001, 'charge': -3.63322}
#        	}
#    			},
#    		"Ti": {
#        			'core': {'mass': 47.867, 'charge': 19.36901},
#        			'shell': {'mass': 0.0001, 'charge': -16.27952}
#   				 },
#    		"O": {
#        			'core': {'mass': 15.999, 'charge': 2.54843},
#        			'shell': {'mass': 0.0001, 'charge': -4.19917}
#    			}
#					}

	def __init__(self,filename="default_file",convention=None,prescribe_N=0,comment="no comment"):
		"""
		__init__ file for the Cell class\n
		:param filename: Structural input file
		:param convention: Type of the structural file if it cannot be recognized from filename
		"""
		self.atom=[]
		self.lattice=None
		self.comment=comment
		self.position_flag=""
		self.species_count=[]
		self.species_name=[]
		#self.shell_model_name=None
		#self.shell_model_values={}
		self.shell_models ={}
		self.N=0

		def zeroList(self,perscribe_N):
			self.position_flag='Direct'
			self.lattice=Lattice()
			self.N=prescribe_N
			aux_atom=Atom()
			for ii in range(self.N):
				aux_atom.position_frac=[0.,0.,0.,]
				aux_atom.name="NN"
				self.atom.append(copy.deepcopy(aux_atom))
	

		def readFromPOSCAR(self,filename):
			"""
			Reads atomic structure from POSCAR (CONTCAR) VASP-file\n
			:param filename: Name of the VASP structural file
			:type filename: str
			"""
			log.info("Reading file in VASP format")
			fin = open(filename)
			
			#get lattice
			aux_lat=Lattice()
			aux_is_names_line=0
			
			self.comment=	fin.readline() # skip first line with comment
			
			line_mult=	fin.readline() # line with multiplication factor
			aux_lat.mult=	float(line_mult)
			
#			aux_line_spl=fin.readline().split()
#			if(len(aux_line_spl) != 3): print("Three inputs expected for lattice vector, terminating..."); exit(1);
#			aux_lat.a=		[float(aux_line_spl[0]),float(aux_line_spl[1]),float(aux_line_spl[2])] # translation vector a
#			aux_line_spl=fin.readline().split()
#			if(len(aux_line_spl) != 3): print("Three inputs expected for lattice vector, terminating..."); exit(1);
#			aux_lat.b=		[float(aux_line_spl[0]),float(aux_line_spl[1]),float(aux_line_spl[2])] # translation vector a
#			aux_line_spl=fin.readline().split()
#			if(len(aux_line_spl) != 3): print("Three inputs expected for lattice vector, terminating..."); exit(1);
#			aux_lat.c=		[float(aux_line_spl[0]),float(aux_line_spl[1]),float(aux_line_spl[2])] # translation vector a
			aux_lat.a=list(map(float,fin.readline().split()))
			aux_lat.b=list(map(float,fin.readline().split()))
			aux_lat.c=list(map(float,fin.readline().split()))

			if(np.linalg.det([aux_lat.a,aux_lat.b,aux_lat.c])==0):
				print("Lattice vectors do not span 3D space, terminating...")
				exit(1)
			

			if aux_lat.mult !=1.:
				aux_lat.a=[aux_lat.mult*xx for xx in aux_lat.a ]
				aux_lat.b=[aux_lat.mult*xx for xx in aux_lat.b ]
				aux_lat.c=[aux_lat.mult*xx for xx in aux_lat.c ]
				aux_lat.mult=1.
			
			aux_lat.getAnglesFromVectors()
			self.lattice=aux_lat

			#get atoms
			aux_atom=Atom()
		
			aux_fin=fin.readline().split()
			try:
				int(aux_fin[0])
			except Exception:
				#additional line with element names
				aux_is_names_line=1
				self.species_name=aux_fin

				aux_fin=fin.readline().split()	#and read the next line


			self.species_count=list(map(int,aux_fin))
			self.N=sum(self.species_count)
			self.position_flag=fin.readline().strip()

			for ii in range(self.N):
				#aux_line=fin.readline().split()
				aux_line=re.sub("[, !#%:\t]+",' ',fin.readline()).split() #splits read line and 
											#removes comments which may 
											#precede name of elements (e.g.
											#in POSCARs from Bilbao
											#crystallographic server)
				pos=list(map(float,aux_line[0:3]))
				aux_atom.position_frac=pos

				if len(aux_line) > 3:
					aux_atom.name=aux_line[3]
				self.atom.append(copy.deepcopy(aux_atom))


			fin.close()
			if aux_is_names_line:
				self.changeElementNames(self.species_name)

		def readFromGulp(self,filename):
			"""
			Reads atomic structure from GULP file\n
			:param filename: Name of the GULP structural file
			"""
			self.mult=1.
			self.position_flag='Direct'
			log.info("Reading file in GULP format")
			fin = open(filename)
			
			#get lattice
			aux_lat=Lattice()
				
			aux_line=fin.readline()
			while aux_line and aux_line.split() != [] and aux_line.split()[0] !="cell":
				aux_line=fin.readline()
			
			[a_len,b_len,c_len,aux_lat.alpha,aux_lat.beta,aux_lat.gamma]=list(map(float,fin.readline().split()))
			#print("  a length: %lf\n  b length: %lf\n  c length: %lf\n  alpha   : %lf\n  beta    : %lf\n  gamma   : %lf" % tuple([a_len,b_len,c_len,aux_lat.alpha,aux_lat.beta,aux_lat.gamma]))
			
			aux_lat.a=pmut.fracToCarthesianFromAngles(a_len,b_len,c_len,aux_lat.alpha,aux_lat.beta,aux_lat.gamma,1,0,0)
			aux_lat.b=pmut.fracToCarthesianFromAngles(a_len,b_len,c_len,aux_lat.alpha,aux_lat.beta,aux_lat.gamma,0,1,0)
			aux_lat.c=pmut.fracToCarthesianFromAngles(a_len,b_len,c_len,aux_lat.alpha,aux_lat.beta,aux_lat.gamma,0,0,1)

			self.lattice=aux_lat

			#get atoms
			aux_atom=Atom()
			aux_line=fin.readline().split()
			neco=aux_line[1]#zjistit vyznam...

			aux_line=fin.readline()
			aux_split=aux_line.split()
			try:
				#PM20151107: In the following condition "aux_split != []" prevents from reading empty line
				while aux_line and aux_split != [] and pmco.isAtom(aux_split[0]):
					if aux_split[1] == "core" or aux_split[1] == "shel":
						pos=list(map(float,aux_split[2:5]))		#second position on the line is usually core/shell
						if aux_split[1] == "shel": aux_atom.coreshell="shell"
					else:
						pos=list(map(float,aux_split[1:4]))
					aux_atom.position_frac=pos
					aux_atom.name=aux_split[0]
					self.N+=1
					self.atom.append(copy.deepcopy(aux_atom))


					aux_line=fin.readline()
					aux_split=aux_line.split()
			except StopIteration:
				{}	#catching end of file

			
			fin.close()

		def readFromDLPoly(self,filename):

			"""
			Reads atomic structure from DLPOLY file\n
			:param filename: Name of the DLPOLY structural file
			"""
			self.mult=1.
			self.position_flag='Direct'
			log.info("Reading file in DLPOLY format")
			fin = open(filename)
			
			#get lattice
			aux_lat=Lattice()
			
			aux_lat.a=		list(map(float,fin.readline().split())) # translation vector a
			aux_lat.b=		list(map(float,fin.readline().split())) # translation vector b
			aux_lat.c=		list(map(float,fin.readline().split())) # translation vector c
		
			aux_lat.getAnglesFromVectors()
			self.lattice=aux_lat

			#get atoms
			aux_atom=Atom()

			aux_line=fin.readline()
			aux_split=aux_line.split()
			try:
				#PM20151107: In the following condition "aux_split != []" prevents from reading empty line
				while aux_line and aux_split != []:
					#read position x y z
					pos=list(map(float,aux_split[1:4]))
					#convert atomic positions to fractional:
					latvec_vertical=np.array([self.lattice.a,self.lattice.b,self.lattice.c]).transpose()	#vertically ordered lattice vectors
					aux_atom.position_frac=copy.deepcopy(np.linalg.solve(latvec_vertical,pos))
					
					aux_atom.name=aux_split[0]
					self.N+=1
					self.atom.append(copy.deepcopy(aux_atom))
			


					aux_line=fin.readline()
					aux_split=aux_line.split()
			except StopIteration:
				{}	#catching end of file
			
			fin.close()
		
		def readFromXYZ(self,filename):

			"""
			Reads atomic structure from XYZ file\n
			:param filename: Name of the structural file
			"""
			self.mult=1.
			self.position_flag='Direct'
			log.info("Reading file in XYZ format")
			fin = open(filename)
			
			aux_lat=Lattice()
				
			aux_line=fin.readline()
			while aux_line and aux_line.split() != [] and aux_line.split()[0] !="lattice":
				aux_line=fin.readline()
			
			#get lattice
			aux_lat=Lattice()
			aux_is_names_line=0
			aux_lat.mult=	1.
			
			aux_lat.a=		list(map(float,fin.readline().split())) # translation vector a
			aux_lat.b=		list(map(float,fin.readline().split())) # translation vector b
			aux_lat.c=		list(map(float,fin.readline().split())) # translation vector c

			if aux_lat.mult !=1.:
				aux_lat.a=[aux_lat.mult*xx for xx in aux_lat.a ]
				aux_lat.b=[aux_lat.mult*xx for xx in aux_lat.b ]
				aux_lat.c=[aux_lat.mult*xx for xx in aux_lat.c ]
				aux_lat.mult=1.
			
			aux_lat.getAnglesFromVectors()
			self.lattice=aux_lat

			#get atoms
			aux_atom=Atom()
			aux_line=fin.readline()
			print(aux_line)
			while aux_line and aux_line.split() != [] and aux_line.split()[0] !="basis":
				aux_line=fin.readline()

			aux_line=fin.readline()
			aux_split=aux_line.split()
			try:
				#PM20151107: In the following condition "aux_split != []" prevents from reading empty line
				while aux_line and aux_split != []:
					#read position x y z
					pos=list(map(float,aux_split[1:4]))
					#convert atomic positions to fractional:
					latvec_vertical=np.array([self.lattice.a,self.lattice.b,self.lattice.c]).transpose()	#vertically ordered lattice vectors
					aux_atom.position_frac=copy.deepcopy(np.linalg.solve(latvec_vertical,pos))
					
					aux_atom.name=aux_split[0]
					self.N+=1
					self.atom.append(copy.deepcopy(aux_atom))
			


					aux_line=fin.readline()
					aux_split=aux_line.split()
			except StopIteration:
				{}	#catching end of file

			
			fin.close()
		
		def readFromLAMMPSStructure(self,filename):#,model_name):
					"""
					Reads atomic structure from LAMMPSStructure file\n
					:param filename: Name of the LAMMPS structural file
					"""
					self.mult=1.
					self.position_flag='Direct'
					log.info("Reading file in LAMMPS format")
					fin = open(filename,"r")		

					with fin as file:
						aux_line = file.readlines()								
					#get lattice
					aux_lat=Lattice()
					aux_atom=Atom()				
					data = {}
					cart_pos_lst = []
					model_data = {}
					model_data['PTO_shimada'] = {}
					for i, lines in enumerate(aux_line):
						words = lines.split()
						#print(words)
						#print(model_data['PTO_shimada'],'at start')
						if "ATOMS" in words and 'NUMBER' in words:
							#print(int(aux_line[i+1]))							
							atm_line = aux_line[i+1]
							temp_line = atm_line.split()
							#print(temp_line[0])
							data['atoms'] =  int(temp_line[0])	
							# hard coded the values of model['atom types'] = 6 since we have AB03  shell model so cadidate rememmber to change it later
							data['atom types'] = 6
							data['bond types'] = 3	
							for jj in range(1,data['atom types']+1):
								atom_type = jj								
								model_data['PTO_shimada'][atom_type] = {}							
								model_data['PTO_shimada'][atom_type]['core'] = {"mass": None, "charge": None }
								model_data['PTO_shimada'][atom_type]['shell'] = {"mass": None,  "charge": None}													
						if "BOUNDS" in words and 'BOX' in words:
							bound_line =  aux_line[i+1]
							temp_line = bound_line.split()
							xlo_bound = float(temp_line[0])
							xhi_bound = float(temp_line[1])
							xy = float(temp_line[2])
							bound_line =  aux_line[i+2]
							temp_line = bound_line.split()
							ylo_bound = float(temp_line[0])
							yhi_bound = float(temp_line[1])
							xz = float(temp_line[2])
							bound_line =  aux_line[i+3]
							temp_line = bound_line.split()
							zlo = float(temp_line[0])
							zhi = float(temp_line[1])
							yz = float(temp_line[2])
							aux_lat.a[0] = (xhi_bound - max(0.0,xy,xz,xy+xz))- (xlo_bound-min(0.0,xy,xz, xy+xz))
							aux_lat.b[0] = xy
							aux_lat.b[1] = (yhi_bound-max(0.0,yz)) - (ylo_bound-min(0.0,xy))
							aux_lat.c[0] = xz
							aux_lat.c[1] = yz
							aux_lat.c[2] = zhi-zlo
							#print(temp_line)
						#print(model_data['PTO_shimada'])	
						if 'atoms' in words and  'ITEM' not in words:
							data['atoms'] = int(words[words.index('atoms') - 1])
						if 'bonds' in words and  'ITEM' not in  words:
							data['bonds'] = int(words[words.index('bonds') - 1])
						if 'atom' in words and 'types' in  words:
							data['atom types'] = int(words[words.index('atom') - 1])
						if 'bond' in words and 'types' in  words:
							data['bond types'] = int(words[words.index('bond') - 1])
						if "xlo" in words and  'ITEM:' not in  words:							
							lx = float(words[words.index('xlo') - 1]) - float(words[words.index('xlo') - 2])
							aux_lat.a[0] =lx
						if "ylo" in words:
							ly = float(words[words.index('ylo') - 1]) - float(words[words.index('ylo') - 2])
							aux_lat.b[1] = ly
						if "zlo" in words:
							lz = float(words[words.index('zlo') - 1]) - float(words[words.index('zlo') - 2])						
							aux_lat.c[2] = lz						
						# for Triclinic (non-orthogonal) simulation boxes
						if ("xy"  in words) and ("BOX" not in words): 
							yz = float(words[words.index('xy') - 1])
							xz = float(words[words.index('xy') - 2])
							xy = float(words[words.index('xy') - 3])
							a = lx
							b = np.sqrt(ly**2 +xy**2)
							c = np.sqrt(lz**2 + xz**2 + yz**2)
							alpha = np.arccos((xy*yz +lx*yz)/(b*c))
							beta = np.arccos(xz/c)
							gamma = np.arccos(xy/b)
							aux_lat.a[0] = a # a_x
							aux_lat.b[0] = b * np.cos(gamma)# b_x
							aux_lat.b[1] = b * np.sin(gamma)# b_y
							aux_lat.c[0] = c * np.cos(beta)# c_x
							aux_lat.c[1] = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma) #c_y
							aux_lat.c[2] = np.sqrt(c**2 - aux_lat.c[0]**2 - aux_lat.c[1]**2) #c_z
							#print(alpha,beta,gamma)
						aux_lat.getAnglesFromVectors()
						self.lattice=aux_lat
						#print(model_data)	
						if "Masses" in words:
							for j in range(2,data['atom types']+2):
								mass_line = aux_line[i+j]
								masses = mass_line.split()
								atom_type = int(masses[0])
								atom_mass = float(masses[1])								
								# Ensure that aux_model_data['PTO_shimada'][atom_type] exists
								if 'PTO_shimada' not in model_data:
									model_data['PTO_shimada'] = {}
								if atom_type not in model_data['PTO_shimada']:
									model_data['PTO_shimada'][atom_type] = {}							
								model_data['PTO_shimada'][atom_type]['core'] = {"mass": None, "charge": None }
								model_data['PTO_shimada'][atom_type]['shell'] = {"mass": None,  "charge": None}
								if(atom_type<=int(data['atom types']/2)):
									model_data['PTO_shimada'][atom_type]["core"]["mass"] = atom_mass
								else:
									model_data['PTO_shimada'][atom_type]["shell"]["mass"] = atom_mass

						#print(data['atoms'],"before_atoms")
						model ={}
						model['total_atoms'] = {}
						if ("Atoms" in words) or ("x" in words and "y" in words and "z" in words) :
							ii =2 
							if ("x" in words and "y" in words and "z" in words):
								ii=1							
							 
							#for j in range(2,model['atoms']+2):
							for j in range(ii,data['atoms']+ii):
								atoms_line = aux_line[i+j]
								atoms_details = atoms_line.split()
								if ii == 2:
									atoms_tag = int(atoms_details[0])
									bond_tag = int(atoms_details[1])
									atom_type = int(atoms_details[2])
									atom_charge = float(atoms_details[3]) 
									pos = [float(atoms_details[4]),float(atoms_details[5]),float(atoms_details[6])]
								else:
									atoms_tag = int(atoms_details[0])
									bond_tag = 0
									atom_charge = 0.0
									atom_type = int(atoms_details[1])
									pos = [float(atoms_details[3]),float(atoms_details[4]),float(atoms_details[5])]														
								cart_pos_lst.append(pos)
								if atoms_tag not in model['total_atoms']:
									model['total_atoms'][atoms_tag] = {}
									model['total_atoms'][atoms_tag]["bond_tag"] = bond_tag
									model['total_atoms'][atoms_tag]["atom_type"] = atom_type
									model['total_atoms'][atoms_tag]["atom_charge"] = atom_charge
									model['total_atoms'][atoms_tag]["pos"] = pos										
								#print(j,atom_type,atoms_tag)
								if(atom_type>int(data['atom types']/2)):
									aux_atom.coreshell = "shell"
									aux_atom.position_frac = [0,0,0]
									aux_atom.name = str(atom_type)
									aux_atom.charge = atom_charge
									#model_data['PTO_shimada'][atom_type] = atom_type
									model_data['PTO_shimada'][atom_type]["shell"]["charge"] = atom_charge
									self.N+=1
									self.atom.append(copy.deepcopy(aux_atom))
								else:
									aux_atom.coreshell = "core"
									aux_atom.position_frac =  [0,0,0]
									aux_atom.name = str(atom_type)
									aux_atom.charge = atom_charge
									self.N+=1
									self.atom.append(copy.deepcopy(aux_atom))
									#print(model_data)
									model_data['PTO_shimada'][atom_type]["core"]["charge"] = atom_charge									
					cart_pos_lst = np.array(cart_pos_lst)			
					self.setCartesian(cart_pos_lst)
					self.shell_models = copy.deepcopy(model_data)
					

		if 	convention == 'POSCAR' or \
			re.findall('VASP',filename.upper()) or \
			re.findall('CONTCAR',filename.upper()) or \
			re.findall('POSCAR',filename.upper()) or \
			re.findall('CHGCAR',filename.upper()) :
			#re.findall('vasp',filename)    or re.findall('VASP',filename) or \
			#re.findall('contcar',filename) or re.findall('CONTCAR',filename) or \
			#re.findall('poscar',filename)  or re.findall('POSCAR',filename) or \
			#re.findall('chgcar',filename)  or re.findall('CHGCAR',filename) :
			readFromPOSCAR(self,filename)
		elif convention == 'GULP' or os.path.splitext(filename)[1] == '.GULP' or os.path.splitext(filename)[1] == '.gulp':
			readFromGulp(self,filename)
		elif convention == 'DLPOLY' or os.path.splitext(filename)[1] == '.DLPOLY' or os.path.splitext(filename)[1] == '.dlpoly':
			readFromDLPoly(self,filename)
		elif convention == 'XYZ' or os.path.splitext(filename)[1] == '.XYZ' or os.path.splitext(filename)[1] == '.xyz':
			readFromXYZ(self,filename)
		elif convention == 'LAMMPSStructure' or os.path.splitext(filename)[1] == '.LAMMPSStructure' or os.path.splitext(filename)[1] == '.lammpssructure':
			readFromLAMMPSStructure(self,filename)
		elif convention == 'zerolist':
			zeroList(self,prescribe_N)
		else:
			
			log.warning("Unknown convention \""+str(convention)+"\" or file type \""+filename+"\" for file with atomic structure. Constructing empty cell")
			zeroList(self,prescribe_N)
			
			#log.warning("Unknown file type, file "+filename+" cannot be used to get atomic structure")

		self.countPresentSpecies()
		#if self.species_name: self.orderSpecies()	#if the atom names are already assigned and species evaluated

	def __str__(self):
		""" Prints information about cell in a standard format """
		str_string=str(self.lattice)+"\n"
		str_string+="Number of atoms: %d\n"%self.N
		str_string+="Atoms:\n"
		str_string=str_string+"   %-4s %-3s %-32s %-7s %-4s %-4s"%("Id","At","Position (fractional)","Charge","c/s","Wyc")+"\n"
		for ii in range(self.N):
  			str_string=str_string+"   % 3d: %3s (% .6f % .6f % .6f) %+7.2f %3s %3s\n" % (	ii,
  												self.atom[ii].name,
  												self.atom[ii].position_frac[0],
  												self.atom[ii].position_frac[1],
  												self.atom[ii].position_frac[2],
  												self.atom[ii].charge,
  												self.atom[ii].coreshell,
  												str(self.atom[ii].wyckoff),
  											    )
		if self.N==0:
			str_string=str_string+"      No atoms prescribed"
		str_string=str_string+"\n"

		str_string=str_string+"Species count: "+str(self.species_count)+"\n"
		str_string=str_string+"Species names: "+str(self.species_name)+"\n"
		str_string.strip()
		return str_string

	def orderSpecies(self):
		""" This goes through the list of atoms and groups the species according to species_name list """

		new_list=[]
		for ii in self.species_name:
			for jj in self.atom:
				if jj.name == ii:
					new_list.append(copy.deepcopy(jj))

		del self.atom
		self.atom=new_list

	def replicate(self,R):
		"""R is list of [n_x,n_y,n_z] where n idicates the how many times this cell will repeat in respective direction"""
		if len(R)!=3:
			log.error("wrong number of elements in list R")
			return
		if all(not isinstance(ii, int) and ii < 0 for ii in R):
			log.error("please provide positive intergers ")
			return
		x = self.lattice.a
		y = self.lattice.b
		z = self.lattice.c
		xr = [xx * R[0] for xx in x]
		yr = [yy * R[1] for yy in y]
		zr = [zz * R[2] for zz in z]
		cart_cord = self.getArrayCarthesian()
		repllica  = Cell(prescribe_N=0)
		repllica.lattice.a = xr
		repllica.lattice.b = yr
		repllica.lattice.c = zr
		repllica.lattice.getAnglesFromVectors()
		nx = R[0]
		ny = R[1]
		nz = R[2]
		cart_cord_repllica =[]
		for ii in range(0,nz):
			for jj in range(0,ny):
				for kk in range(0,nx):					
					for aa in range(len(cart_cord)):
						ll = cart_cord[aa]
						pp = ll.tolist()                                    
						pp[0] = pp[0] + (kk * x[0])
						pp[1] = pp[1] + (jj * y[1])
						pp[2] = pp[2] + (ii * z[2])
						cart_cord_repllica.append(pp)
						repllica.appendAtom([0,0,0],self.atom[aa].name)
		cart_cord_repllica = np.array(cart_cord_repllica)
		repllica.setCartesian(cart_cord_repllica)
		repllica.shell_models = copy.deepcopy(self.shell_models)
		return repllica
	

	def printLatex(self,filename="structure.tex"):
		""" Print information about cell in Latex format """
		aaa=10
		f=open(filename,"w")

		f.write("\\documentclass[aps,prb,10pt,twocolumn,amsmath,amssymb,showpacs,letterpaper,showkeys]{revtex4-1}\n")
		#f.write("\\documentclass{article}\n")
		f.write("\\usepackage{graphicx}\n")
		f.write("\\usepackage{epsfig}\n")
		f.write("\\begin{document}\n")

		f.write("\\begin{table}\n")
		f.write("\\begin{tabular}{ll}\n")
		f.write("\\hline \hline\n")

		f.write(r"$|a|$    & %s\\"%pmla.writeNumber(np.linalg.norm(self.lattice.a),format="f",prec=4));
		f.write(r"$|b|$    & %s\\"%pmla.writeNumber(np.linalg.norm(self.lattice.b),format="f",prec=4));
		f.write(r"$|c|$    & %s\\"%pmla.writeNumber(np.linalg.norm(self.lattice.c),format="f",prec=4));
		f.write(r"$\alpha$ & %s\\"%pmla.writeNumber(self.lattice.alpha,format="f",prec=4));
		f.write(r"$\beta$  & %s\\"%pmla.writeNumber(self.lattice.beta, format="f",prec=4));
		f.write(r"$\gamma$ & %s\\"%pmla.writeNumber(self.lattice.gamma,format="f",prec=4));
		f.write("\n")
		f.write("\hline\n")
		for ii in range(self.N):
			f.write(r"%s    & (%s,%s,%s)\\"%(	self.atom[ii].name,
								pmla.writeNumber(self.atom[ii].position_frac[0],format="f",prec=4),
								pmla.writeNumber(self.atom[ii].position_frac[1],format="f",prec=4),
								pmla.writeNumber(self.atom[ii].position_frac[2],format="f",prec=4)))
			f.write("\n")
												

			
		f.write("\\hline \hline\n")
		f.write("\\end{tabular}\n")
		f.write(r"\caption{CAPTION}")
		f.write("\n")
		f.write("\\label{tab_properties}\n")
		f.write("\\end{table}\n")

		f.write("\\bibliography{/home/vader/jabref/marton}\n")
		f.write("\\end{document}\n")
		f.close()

	def setElementSpins(self,element_spins):
		log.info("Changing spins to: %s"%str(element_spins))
		print("Changing spins to: %s"%str(element_spins))
		if not self.N == len(element_spins):
			log.warning("Number of provided spins in \"setElementSpins()\" is not equal to number of elements")

		for ii in range(self.N):
			if ii < len(element_spins):
				self.atom[ii].spin=float(element_spins[ii])
	
	def setElementNames(self,element_names):
		"""
		Another name for :func:`~pm__cell.Cell.changeElementNames`
		"""
		self.changeElementNames(element_names)

	def changeElementNames(self,element_names):
		"""
		Change or set names of element in the cell.\n
		:param element_names: contains list of element names
		:type element_names: str
		"""
		#print(element_names)
		#print(self.species_count)
		#exit(1)
		if len(element_names) == len(self.species_count):
			log.info("Changing names of species to: " + str(element_names))
			self.species_name=[element_names[ii] for ii in range(len(self.species_count))]
			aux_cnt1=0; aux_cnt2=0
			for ii in range(len(self.species_count)):
				aux_cnt1=copy.deepcopy(aux_cnt2)
				aux_cnt2=copy.deepcopy(aux_cnt2+self.species_count[ii])
				for jj in range(aux_cnt1,aux_cnt2):
					self.atom[jj].name=self.species_name[ii]

		elif len(element_names) == self.N:
			#print("List is utilized to name individual atoms")
			for ii in range(self.N):
				self.atom[ii].name=element_names[ii]

		else:			
			print("Number of given strings (%s) is not equal to number of species (%s) or total number of atoms in the system (%s), terminating"%(len(element_names),len(self.species_count),self.N))
			sys.exit(1)

		self.countPresentSpecies()

	
	def countPresentSpecies(self):
		"""
		Fills structures species_count and species_name with number of species and their names, respectively.
		"""
		aux_species_dictionary={}
		aux_species_list=[]

		
		if 	(not self.species_count and not self.species_name) or \
			(    self.species_count and     self.species_name):
			#in the first case first counting ever
			#in the second case rocounting
			self.species_count=[]
			self.species_name=[]
			
			for ii in range(self.N):
				if self.atom[ii].name in aux_species_dictionary.keys():	
					aux_species_dictionary[self.atom[ii].name]+=1
				else:
					aux_species_dictionary[self.atom[ii].name] =1
					aux_species_list.append(self.atom[ii].name)

			for ii in range(len(aux_species_list)):
				sp=aux_species_list[ii]
				self.species_count.append(aux_species_dictionary[sp])	
				self.species_name.append(sp)

		if self.species_count and not self.species_name:
			#number of species are prescribed, but not names (e.g. from older POSCARs)
			aux_store_species_count=copy.deepcopy(self.species_count)
			for ii in range(self.N):
				aux_species_dictionary[self.atom[ii].name]=aux_species_dictionary.get(self.atom[ii].name,0)+1
				if not self.atom[ii].name in self.species_name:
					self.species_name.append(self.atom[ii].name)
				
			self.species_count=[0 for ii in self.species_name]
			for ii in range(len(self.species_name)):
				self.species_count[ii]=aux_species_dictionary[self.species_name[ii]]

			if not self.species_count == aux_store_species_count:
				log.warning("Species cannot be determined from structure file")
				self.species_count=copy.deepcopy(aux_store_species_count)
				self.species_name=[]
		
		if not self.species_count and self.species_name:
			self.species_count=[0 for ii in self.species_name]
			log.warning("Situation when species_name are known, but counts is not yet implemented")

	def removeAtom(self,atomid):
		if atomid<self.N:
			#print "Removing atom number %d"%atomid
			#print self.atom[atomid]
			del(self.atom[atomid])
			self.N-=1
			
			self.countPresentSpecies()
		else:
			log.warning("Atom not removed, index out of range")

	def appendAtom(self,position_frac=[0.,0.,0.,],name="undef",charge=0.0,coreshell="core",wyckoff="",spin=0.):
		
		
		aux_atom=Atom()
		self.atom.append(copy.deepcopy(aux_atom))
	
		self.atom[-1].name=name
		self.atom[-1].position_frac[0]=position_frac[0]
		self.atom[-1].position_frac[1]=position_frac[1]
		self.atom[-1].position_frac[2]=position_frac[2]
		self.atom[-1].charge=charge
		self.atom[-1].coreshell=coreshell
		self.atom[-1].wyckoff=wyckoff
		self.atom[-1].spin=spin

		self.N+=1
		self.countPresentSpecies()

	def rearrangeSpecialSystem(self):
		"""
		For some known systems which we frequently work with this procedure allow automatic reordering of lattice vectors
		to a common setting.\n
		Known systems are:\n
		BFO_12: transformation to almost 90deg. setting, and order lattice vectors from longest to shortest\n
		"""	
		if self.getSystemType() == 'unknown':
			log.warning("Unknown system, no rearrangement is applied")
		else:
			import operator
			#sort lattice vectors from longest to shortest
			aux=[[np.linalg.norm(self.lattice.a),0],[np.linalg.norm(self.lattice.b),1],[np.linalg.norm(self.lattice.c),2]]
			aux.sort(key=operator.itemgetter(0),reverse=True)
			aux_order=[aux[ii][1] for ii in range(len(aux))]

			if(self.getSystemType() == 'BFO_12'):
				#Project BFO under pressure 
				log.info("System BFO 12 found, choosing special arrangement of unticell (a=longest, b=sqrt(2)a, c=2a)")
				
				#swap the last two indices in array
				aux_swap=copy.deepcopy(aux_order[1])
				aux_order[1]=aux_order[2]
				aux_order[2]=aux_swap
				self.permuteDirections(aux_order)
				self.reorderElements(['Bi','Fe','O'])

				#check angles whether they are close to 90 degrees, which is the preferred settin
				if abs(self.lattice.alpha-90) > 5 or abs(self.lattice.beta-90) > 5 or abs(self.lattice.gamma-90) > 5:
					if (self.lattice.gamma-108)<5:
						log.info("108 degree problem, using a'=a+b, b'=b, c'=c transformation")
						log.info("angles: alpha=%f, beta=%f, gamma=%f"%(self.lattice.alpha,self.lattice.beta,self.lattice.gamma))
						
						self.tranform(np.array([[1.,1.,0.],[0.,1.,0.],[0.,0.,1.]]))
						
						self.lattice.getAnglesFromVectors()
						log.info("angles: alpha=%f, beta=%f, gamma=%f"%(self.lattice.alpha,self.lattice.beta,self.lattice.gamma))


	def tranform(self,M):
		""" M is transformation matrix """
		if np.shape(M) != (3,3):
			log.error("wrong shape of transformation matrix, transformation is not performed")
			return
		if np.linalg.det(M) == 0:
			log.error("Transformation matrix is not regular, transformation is not performed")
			return
		if np.linalg.det(M) != 1:
			log.error("Transformation mcaatrix is not unitary (changes volume), transformation is not performed")
			return

		log.info("Tranformation matrix (det=%lf):\n"%np.linalg.det(M) + str(M))
		lat_orig=np.array([self.lattice.a,self.lattice.b,self.lattice.c])
		lat_prime=np.dot(M,lat_orig)
		for ii in range(3):
			for jj in range (3):
				if abs(lat_prime[ii,jj])<1E-10: lat_prime[ii,jj]=0.

		self.lattice.a=copy.deepcopy(lat_prime[0,:])
		self.lattice.b=copy.deepcopy(lat_prime[1,:])
		self.lattice.c=copy.deepcopy(lat_prime[2,:])

		iM=np.linalg.inv(M)
		log.info("Inversed tranformation matrix:\n" + str(iM))
		TiM=iM.transpose()
		log.info("Inversed and transposed tranformation matrix:\n" + str(TiM))

		#Tranfrom all the fractional coordinates
		for ii in range(self.N):
			aux=self.atom[ii].position_frac
			self.atom[ii].position_frac=np.dot(TiM,aux)

	def bringAtomsInsideCell(self):
		""" Add integer multiple of lattice constant to the current position of the atom in order that all its fractional
		coordinates are in range <0,1) """
		for ii in range(self.N):
			for jj in range(3):
				self.atom[ii].position_frac[jj]=self.atom[ii].position_frac[jj] % 1
				if abs(self.atom[ii].position_frac[jj]-1.) < 1E-10: self.atom[ii].position_frac[jj] = 0.	#added 30.5.2019 in order to avoid eE-17 -> 1.0000000000


	def writeToXYZ(self, filename):
		nfile=os.path.splitext(filename)[0] + ".xyz"
		log.info("Writing file to " + nfile + " in xyz format")
		fout = open(nfile,"w")

		fout.write("%d\n\n"%self.N)

		#fout.write("#   % 15.10f % 15.10f % 15.10f\n" % tuple(self.lattice.a))
		#fout.write("#   % 15.10f % 15.10f % 15.10f\n" % tuple(self.lattice.b))
		#fout.write("#   % 15.10f % 15.10f % 15.10f\n" % tuple(self.lattice.c))
		#fout.write("\n")
		
		#tisk v kartezskych souradnicich
		for ii in range(self.N):
			#get carthesian coordintes
			x1=[self.atom[ii].position_frac[0]*self.lattice.a[jj] for jj in range(3)]
			x2=[self.atom[ii].position_frac[1]*self.lattice.b[jj] for jj in range(3)]
			x3=[self.atom[ii].position_frac[2]*self.lattice.c[jj] for jj in range(3)]
			pos=list(map(sum,zip(x1,x2,x3)))
			fout.write("%2s  " % self.atom[ii].name + "% 15.10f % 15.10f % 15.10f\n" % tuple(pos))

		fout.close()

	def writeToCIF(self, filename):
		nfile=os.path.splitext(filename)[0] + ".cif"
		log.info("Writing file to " + nfile + " in cif format")
		fout = open(nfile,"w")

		fout.write("data_I\n\n")
		fout.write("_cell_length_a %.6f\n" % np.linalg.norm(self.lattice.a))
		fout.write("_cell_length_b %.6f\n" % np.linalg.norm(self.lattice.b))
		fout.write("_cell_length_c %.6f\n" % np.linalg.norm(self.lattice.c))
		fout.write("_cell_angle_alpha %.6f\n" % self.lattice.alpha)
		fout.write("_cell_angle_beta  %.6f\n" % self.lattice.beta )
		fout.write("_cell_angle_gamma %.6f\n" % self.lattice.gamma)
		fout.write("\nloop_\n")
		fout.write("\t_atom_site_type_symbol\n")
		fout.write("\t_atom_site_fract_x\n")
		fout.write("\t_atom_site_fract_y\n")
		fout.write("\t_atom_site_fract_z\n")
		for ii in range(self.N):
			fout.write("\t%2s  " % self.atom[ii].name + "% 15.10f % 15.10f % 15.10f\n" % tuple(self.atom[ii].position_frac))
		
		#tisk v kartezskych souradnicich, toto je spravne, ale nepouzite
		#for ii in range(self.N):
		#	#get carthesian coordintes
		#	x1=[self.atom[ii].position_frac[0]*self.lattice.a[jj] for jj in range(3)]
		#	x2=[self.atom[ii].position_frac[1]*self.lattice.b[jj] for jj in range(3)]
		#	x3=[self.atom[ii].position_frac[2]*self.lattice.c[jj] for jj in range(3)]
		#	pos=list(map(sum,zip(x1,x2,x3)))
		#	fout.write("%2s  " % self.atom[ii].name + "% 15.10f % 15.10f % 15.10f\n" % tuple(pos))

		fout.close()


	def writeToLAMMPSStructure(self, filename):
			""" writes the current cell to LAMMPSStructure file """
			nfile=os.path.splitext(filename)[0] + ".LAMMPSStructure"
			log.info("Writing file to " + nfile + " in LAMMPSStructure format")
			fout = open(nfile,"w")
			fout.write(f"\n")
			fout.write(f"\n")
			fout.write(f"{int(self.N)}\tatoms\n")
			fout.write(f"{int(self.N/2)}\tbonds\n\n")
			fout.write(f"{int(len(self.species_count))}\tatom types\n")
			fout.write(f"{int(len(self.species_count)/2)}\tbond types\n\n")
			a = np.linalg.norm(self.lattice.a)
			b = np.linalg.norm(self.lattice.b)
			c = np.linalg.norm(self.lattice.c)
			lx = a
			xy = b * np.cos((self.lattice.gamma)*(np.pi/180))
			xz = c * np.cos(self.lattice.beta*(np.pi/180))
			ly = np.sqrt(b**2 - xy**2)
			yz = ((b*c*np.cos(self.lattice.alpha*(np.pi/180))) - (xy*xz))/ly
			lz = np.sqrt(c**2-xz**2-yz**2)
			fout.write(f"{0.0}\t {lx:.6f}\t xlo xhi\n")
			fout.write(f"{0.0}\t {ly:.6f}\t ylo yhi\n")
			fout.write(f"{0.0}\t {lz:.6f}\t zlo zhi\n")
			fout.write(f"{xy:.6f}\t {xz:.6f}\t {yz:.6f}\t xy xz yz\n\n")

			fout.write("Masses\n\n")
			for ii in range(len(self.species_count)):
				if(self.shell_models['PTO_shimada'][ii+1]['core']['mass']!= None):
					fout.write(f"{ii+1}\t {self.shell_models['PTO_shimada'][ii+1]['core']['mass']} \n")
				else:
					fout.write(f"{ii+1}\t {self.shell_models['PTO_shimada'][ii+1]['shell']['mass']} \n")	
			#fout.write(f"\n")
			fout.write(f"\nAtoms\n\n")
			jj=0
			for ii in range(self.N):
				#pos = tuple(self.atom[ii].position_frac)				
				#get carthesian coordintes
				x1=[self.atom[ii].position_frac[0]*self.lattice.a[jj] for jj in range(3)]
				x2=[self.atom[ii].position_frac[1]*self.lattice.b[jj] for jj in range(3)]
				x3=[self.atom[ii].position_frac[2]*self.lattice.c[jj] for jj in range(3)]
				pos=list(map(sum,zip(x1,x2,x3)))
				if (ii %2) == 0:
					jj=jj +1 								
				for kk in range(len(self.species_count)):
					atom_tag =	int(self.species_name[kk])
					if self.atom[ii].name == self.species_name[kk]:
						c_charge =  self.shell_models['PTO_shimada'][atom_tag]['core']['charge'] 
						s_charge =  self.shell_models['PTO_shimada'][atom_tag]['shell']['charge']
						if(c_charge!= None):
							fout.write("% 3d % 3d  " % (ii+1,jj) +  self.atom[ii].name + " % 15.10f" % c_charge + "% 15.10f % 15.10f % 15.10f\n" % tuple(pos))
						else:
							fout.write("% 3d % 3d  " % (ii+1,jj) +  self.atom[ii].name + " % 15.10f" % s_charge + "% 15.10f % 15.10f % 15.10f\n" % tuple(pos))
			fout.write(f"\nBonds\n\n")
			kk=1
			for ii in range(int(self.N)):
						if int(self.atom[ii].name) <= int(len(self.species_count)/2): 
							fout.write( "%3d \t %3d \t %3d \t %3d \n" % (kk, int(self.atom[ii].name),ii+1,ii+2))
							kk = kk + 1

			fout.write(f"\nCS-Info\n\n")
			jj=0
			for ii in range(int(self.N)):
				if (ii %2) == 0:
					jj=jj +1 
				fout.write("% 3d \t % 3d \n" % (ii+1,jj))	
			fout.close()




	def writeToFINDSYM(self,filename):
		if not self.species_name or not self.species_count:
			log.warning("You need to specify elements in order to write to findsym format")
			return

		nfile=os.path.splitext(filename)[0] + ".findsym"
		log.info("Writing file to " + nfile + " in findsym format")
		fout = open(nfile,"w")

		fout.write("#Designed to provide inputs for online form FINDSYM on: http://stokes.byu.edu/iso/findsym.php\n\n" % np.linalg.norm(self.lattice.a))

		fout.write("a     %.10f\n" % np.linalg.norm(self.lattice.a))
		fout.write("b     %.10f\n" % np.linalg.norm(self.lattice.b))
		fout.write("c     %.10f\n" % np.linalg.norm(self.lattice.c))
		fout.write("alpha %.10f\n" % self.lattice.alpha)
		fout.write("beta  %.10f\n" % self.lattice.beta )
		fout.write("gamma %.10f\n" % self.lattice.gamma)
		
		for ii in range(len(self.species_name)):
			for jj in range(self.species_count[ii]):
				fout.write(self.species_name[ii]+" ")
		
		fout.write("\n")
		
		for ii in range(self.N):
			fout.write("  % 15.10f % 15.10f % 15.10f\n" % tuple(self.atom[ii].position_frac))

		fout.close()

	def rotateCell(self,matrix):
		assert	(np.linalg.det(np.linalg.matrix_power(matrix,4))-1) < 1E-8 or\
			(np.linalg.det(np.linalg.matrix_power(matrix,3))-1) < 1E-8,\
			"Matrix must define operation which converts a cube, centered in the origin of Cartesian coordinates with faces in direction of axes, to itself" 


		latt_ideal	=np.array([[0,1,1,],[1,0,1,],[1,1,0,]]).T	#cell shape in terms of pseudocubic vectors of the crystal

		self_latt	=self.lattice.getLatticeArrayCol()		#lattice vectors of the initial cell
		self_cart	=self.getArrayCarthesian()			#cartesian coordinates of atoms in initial cell

		cart	=np.dot(matrix,self_cart.T).T				#rotated cartesian positions
		latt_aux=np.dot(matrix,self_latt)				#rotated lattice vectors
		
		coefs	=np.linalg.solve(np.dot(matrix,latt_ideal),latt_ideal)	#calculation of coefficients, which must be used in order to bring te 
										#rotated lattice vector to the orignal setting
		latt	=np.array(np.dot(latt_aux,coefs))			#using the coefficients to define new lattice vectors consistent with
										#those before rotation
		aux_cell=copy.deepcopy(self)
		self.setCartesian(cart)						#setting of rotated cartesian positions of atoms
		self.setLatticeVectors(latt.T)					#setting lattice vectors (transposition: row vectors are expected)
		
		self.bringAtomsInsideCell()					
		del(aux_cell)

	#def writeToPOSCAR(self,filename="POSCAR_no_name",nowrite=False,decimal=10):
	def writeToPOSCAR(self,filename="POSCAR_no_name",nowrite=False,decimal=10):
		"""
		Writes structure POSCAR (CONTCAR) VASP file
		:param filename: Name of the output file
		:type filename: str
		:return string
		"""
		if not (re.findall("POSCAR",filename) or re.findall("CONTCAR",filename)):
			nfile=os.path.splitext(filename)[0] + ".POSCAR"
		else:
			nfile=filename
		log.info("Writing to " + nfile + " in VASP format")
		#TODO: remove negative zero
		#TODO: allow prescribing of decimal places for output
		aux_string=""
		aux_string=("#  %s" % self.comment) # removed \n tag getting error in VASP 6.4.2 when deepcopy existing cell 
		aux_string+=("% 15.10f\n" % self.lattice.mult)		
		aux_string+=("% 15.10f % 15.10f % 15.10f\n" % tuple(self.lattice.a))
		aux_string+=("% 15.10f % 15.10f % 15.10f\n" % tuple(self.lattice.b))
		aux_string+=("% 15.10f % 15.10f % 15.10f\n" % tuple(self.lattice.c))

		for ii in range(len(self.species_name)):
			aux_string+=(" %3s" % self.species_name[ii])
		aux_string+=("\n")
		for ii in range(len(self.species_count)):
			aux_string+=(" %3d" % self.species_count[ii])
		aux_string+=("\n")

		aux_string+=("%s\n" % self.position_flag)

		for ii in range(len(self.species_count)):
			for jj in range(self.N):
				if(self.species_name[ii]==self.atom[jj].name):
					aux_string+=("% 15.10f % 15.10f % 15.10f" % tuple(self.atom[jj].position_frac) + " %2s" % self.atom[jj].name + "\n")
					#aux_string+=("% 15.*f % 15.*f % 15.*f" % ((decimal,self.atom[jj].position_frac[0]),(decimal,self.atom[jj].position_frac[1]),(decimal,self.atom[jj].position_frac[2]),) + " %2s" % self.atom[jj].name + "\n")
		

		if nowrite == False:
			fout = open(nfile,"w")
			fout.write(aux_string)
			fout.close()

		return aux_string


	def writeToISOTROPY_findsym(self,filename="FINDSYM_no_name",tolerance=0.001):
		"""
		Writes structure for the findsym script from the ISOTROPY suit
		:param filename: Name of the outpu file
		:type filename: str
		"""
		nfile=filename
		log.info("Writing to " + nfile + " in findsym format")

		fout = open(nfile,"w")
		
		fout.write("#%s\n" % self.comment)
		fout.write("%d   tolerance for positions of atoms\n"%tolerance)
		fout.write("2    form of lattice parameters: to be entered as lengths and angles\n")
		fout.write("% 15.10f % 15.10f % 15.10f   % 15.10f % 15.10f % 15.10f  a,b,c,alpha,beta,gamma (in angstroms and degrees)\n" % (np.linalg.norm(self.lattice.a),np.linalg.norm(self.lattice.b),np.linalg.norm(self.lattice.c),self.lattice.alpha,self.lattice.beta,self.lattice.gamma))
		
		fout.write("2      form of primitive lattice vectors\n")
		fout.write("P      primitive\n")

		#fout.write("1      form of primitive lattice vectors\n")
		#fout.write("% 15.10f % 15.10f % 15.10f  primitive lattice vectors\n" % 	tuple(self.lattice.a))
		#fout.write("% 15.10f % 15.10f % 15.10f\n" % 				tuple(self.lattice.b))
		#fout.write("% 15.10f % 15.10f % 15.10f\n" % 				tuple(self.lattice.c))

		fout.write("%d            number of atoms in the primitive unit cell\n"%self.N)
		for ii in range(len(self.species_count)):
			for jj in range(self.species_count[ii]):
				fout.write("%d "%(ii+1))
		fout.write(" type of each atom\n")

		for ii in range(self.N):
			fout.write("% 15.10f % 15.10f % 15.10f" % tuple(self.atom[ii].position_frac) + " %2s" % self.atom[ii].name + "\n")
		
		fout.close()


	def writeToGULP(self,filename,decimal=6):
		""" Writes structure to GULP format """
		nfile=os.path.splitext(filename)[0] + ".GULP"
		log.info("Writing file to " + nfile + " in GULP format")
		fout = open(nfile,"w")
		
		fout.write("cell\n")
		
		fout.write("% 15.10f % 15.10f % 15.10f % 10.6f % 10.6f % 10.6f\n" % (	np.linalg.norm(self.lattice.a),
											np.linalg.norm(self.lattice.b),
											np.linalg.norm(self.lattice.c),
											self.lattice.alpha,
											self.lattice.beta,
											self.lattice.gamma))
		
		fout.write("fractional 1\n")
		
		for ii in range(self.N):
			fout.write(" % 2s" % self.atom[ii].name + "% 15.10f % 15.10f % 15.10f" % tuple(self.atom[ii].position_frac) + "\n")
		
		fout.close()

	def displaceCarthesian(self,displacement_list):
		"""
		Displaces atoms in the cell by amount specified in displacement_list.\n
		:param displacement_list: list of three-component lists representing displacement for each atom in Carthesian coordinates. 
		"""
		aux_natom=len(self.atom)
		aux_ndisp=len(displacement_list)
		if aux_natom != aux_ndisp:
			log.warning("Number of displacements is not equal to number of atoms incell")

	#	aux_displ_frac=np.array([self.atom[ii].position_frac for ii in range(self.N)])
	#	aux_carth=pmut.toCarthesian(aux_displ_frac,self.lattice.a,self.lattice.b,self.lattice.c)
	#	aux_carth_displaced=displacement_list
		displacement_list_fractional=pmut.toFractional(displacement_list,self.lattice.a,self.lattice.b,self.lattice.c)
		for ii in range(self.N):
			self.atom[ii].position_frac=np.add(self.atom[ii].position_frac,displacement_list_fractional[ii]).tolist()

	def getArrayFractional(self):
		"""
		Returns a (N,3)-dimensional numpy-array with fractional displacements of all atoms in the cell.
		"""
		aux_arr=[]						#auxiliary empty list
		for ii in range(self.N):				#for all atoms in the cell
			aux_arr.append(self.atom[ii].position_frac)	#append fractional coordinates to the list
		return np.array(aux_arr)				#and return the list converted to numpy array

	def getArrayCarthesian(self):
		ar_frac=self.getArrayFractional()
		return np.dot(ar_frac,[self.lattice.a,self.lattice.b,self.lattice.c])

	def displaceFractional(self,displacement_list):
		aux_natom=len(self.atom)
		aux_ndisp=len(displacement_list)
		if aux_natom != aux_ndisp:
			if aux_ndisp != 3:
				log.warning("Number of displacements is not equal to number of atoms incell")
			else:
				aux_displacement_list=np.ones((aux_natom,3))
				for ii in range(aux_natom):
					aux_displacement_list[ii][0]=displacement_list[0]
					aux_displacement_list[ii][1]=displacement_list[1]
					aux_displacement_list[ii][2]=displacement_list[2]
				displacement_list=aux_displacement_list#.tolist()

		for ii in range(aux_natom):
			self.atom[ii].position_frac=np.add(self.atom[ii].position_frac,displacement_list[ii]).tolist()

	def displaceRandom(self,maxdispl=0.01):
		"""
		Displaces all atoms from their actual position by a random vector. Displacement is done in Carthesian coordinates
		:param displ: maximum displacement in Angstroems
		"""
		aux_displ=maxdispl*2*(np.random.rand(self.N,3)-1.)
		aux_displ_frac=pmut.toFractional(aux_displ,self.lattice.a,self.lattice.b,self.lattice.c)

		for ii in range(self.N):
			self.atom[ii].position_frac=np.add(self.atom[ii].position_frac,aux_displ[ii]).tolist()
		
	def permuteDirections(self,permutation):
		#permutation specifies which direction will move to first position, second and third: [1,0,2] permutes directins 0 and 1
		log.info("Permuting traslation vectors of the lattice: %d->%d, %d->%d, %d->%d" % (permutation[0],0,permutation[1],1,permutation[2],2))
		aux=[copy.deepcopy(self.lattice.a),copy.deepcopy(self.lattice.b),copy.deepcopy(self.lattice.c)]
		vec=[self.lattice.a,self.lattice.b,self.lattice.c]
		for ii in range(3):
			for jj in range(3):
				vec[ii][jj]=aux[permutation[ii]][permutation[jj]]
		self.lattice.getAnglesFromVectors()	#update angles
		
		#update fractional coordinates
		for ii in range(self.N):
			aux=copy.deepcopy(self.atom[ii].position_frac)
			for jj in range(3):
				self.atom[ii].position_frac[jj]=aux[permutation[jj]]
		
	def reorderElements(self,order):
		""" Change order of elements in the record
		:param order: E.g. order=['Bi','Fe','O'] 
		"""
		log.info("Changing order or elements to: " + str(order))

		aux_permut_atoms=[]
		aux_permut_present=[]
		for spc in order:
			for ii in range(self.N):
				if self.atom[ii].name == spc:
					aux_permut_atoms.append(ii)
			for ii in range(len(self.species_name)):
				if self.species_name[ii] == spc:
					aux_permut_present.append(ii)

		if self.N != len(aux_permut_atoms):
			log.error("Error in permuting atoms, no permutation is performed")
			return
			
		aux1=copy.deepcopy(self.atom)
		for ii in range(self.N):
			self.atom[ii]=aux1[aux_permut_atoms[ii]]
	
		aux1=copy.deepcopy(self.species_count)
		aux2=copy.deepcopy(self.species_name)
		for ii in range(len(self.species_name)):
			self.species_count[ii]   =aux1[aux_permut_present[ii]]
			self.species_name[ii]=aux2[aux_permut_present[ii]]
	
	def getSystemType(self):
		if	len(self.species_count) == 3 and \
			'Bi' in self.species_name and \
			'Fe' in self.species_name and \
			'O'  in self.species_name and \
			self.N == 12*5:
			return 'BFO_12'
		else:
			return 'unknown'

	def getCellVolume(self):
		"""
		Returns volume of the cell (parallelepiped shape)
		"""
		return abs(np.dot(self.lattice.a,np.cross(self.lattice.b,self.lattice.c)))

	def getSymmetry(self,tolerance=0.0001):
		self.writeToISOTROPY_findsym("findsym_aux.in",tolerance=tolerance)
		os.system("~/source/ISOTROPY_software_suit/findsym < findsym_aux.in > findsym_aux.out")
		line=pmdf.find_line_with_pattern("findsym_aux.out","Space Group")
		os.system("rm findsym* -f")
		line_split=line.split()

		return [line_split[2],line_split[3],line_split[4]]
	
	def shiftToZero(self,whichatom):
		shift=copy.deepcopy([-self.atom[whichatom].position_frac[ii] for ii in range(3)])
			
		if shift[0]>0.5:	shift[0]=-1.+shift[0]
		if shift[1]>0.5:	shift[1]=-1.+shift[1]
		if shift[2]>0.5:	shift[2]=-1.+shift[2]
		
		if shift[0]<-0.5:	shift[0]=1.+shift[0]
		if shift[1]<-0.5:	shift[1]=1.+shift[1]
		if shift[2]<-0.5:	shift[2]=1.+shift[2]
		
		for ii in range(self.N):
			for jj in range(3):
				self.atom[ii].position_frac[jj]+=shift[jj]
				if abs(self.atom[ii].position_frac[jj]) < 1E-10: self.atom[ii].position_frac[jj]=+0.#remove negative zeros
		
		for ii in range(3): self.atom[whichatom].position_frac[ii]=0.

	def setCartesian(self,coord_list,atom_ids=None):
		""" sets Carthesian positions of atoms. atom_ids can be used to specify for which atoms the changes apply """
		if atom_ids == None: atom_ids=range(len(coord_list))	#if we didn't choose particular indices, it is for all atoms in the cell

		coord_list_fractional=pmut.toFractional(coord_list,self.lattice.a,self.lattice.b,self.lattice.c)
		for ix,ii in zip(atom_ids,range(len(atom_ids))):
			self.atom[ix].position_frac=coord_list_fractional[ii].tolist()

	def setFractional(self,coord_list,atom_ids=None):
		aux_natom=len(self.atom)
		aux_ndisp=len(coord_list)
		if aux_natom != aux_ndisp:
			log.warning("Number of coords is not equal to number of atoms incell")

		if not atom_ids == None:
			for ii in range(atom_ids):
				self.atom[atom_ids[ii]].position_frac=coord_list[ii].tolist()

		else:
			for ii in range(min(aux_natom,aux_ndisp)):
				self.atom[ii].position_frac=coord_list[ii].tolist()
	
	def setLatticeVectors(self,abc):
		#convert position of atoms to carthesian and store
		if(self.N !=0):
			coord_list_carthesian=np.dot(self.getArrayFractional(),[self.lattice.a,self.lattice.b,self.lattice.c])
		else:
			coord_list_carthesian=np.array([])
		
		#set new lattice vectors
		self.lattice.a=abc[0]
		self.lattice.b=abc[1]
		self.lattice.c=abc[2]
		self.lattice.getAnglesFromVectors()
		
		#calculate fractional coordinates from new lattice vectors
		if(self.N !=0):
			coord_list_fractional=pmut.toFractional(coord_list_carthesian,self.lattice.a,self.lattice.b,self.lattice.c)
		else:
			coord_list_fractional=np.array([])

		for ii in range(len(coord_list_fractional)):
			self.atom[ii].position_frac=coord_list_fractional[ii].tolist()

	def match(self,refCell):
		"""
		Displace atoms in the cell by a integer combination of lattice vectors in order that it
		matches position of some atom of the same kind (name) from refCell as good as possible.\n
		Both structures must have identical number of particles\n
		Routine does not change element order.\n
		Shape of both cells is not required identical, but it must be similar (only works with the fractional coordinates)\n
		Notice that different atoms may be moved by different amount\n
		:param 	refCell	: Object containing the reference cell\n
		:return		: Cell with shifted atoms
		"""
	
		"================================================================================================================"					
		if not refCell.N == self.N:		#number of atoms must be equal
			raise AssertionError("match() requires same number of atoms in both atomic structures")
		"================================================================================================================"					
							#lattice vectors are required similar (but not necessarily identical); this
							#condition is not necessary, but might help to avoid mistakes

		if not (np.linalg.norm(np.array(refCell.lattice.a) - np.array(self.lattice.a)) < 0.5 and
			np.linalg.norm(np.array(refCell.lattice.b) - np.array(self.lattice.b)) < 0.5 and
			np.linalg.norm(np.array(refCell.lattice.c) - np.array(self.lattice.c)) < 0.5):
			raise AssertionError("match() requires that both structures have similar lattice vectors")
							
		"================================================================================================================"					
							#The atoms of same kinds and identical number from each kind are required
		atom_self={}
		atom_ref={}
		for ii in range(refCell.N):
			atom_self[self   .atom[ii].name]=atom_self.get(self   .atom[ii].name,0)+1
			atom_ref [refCell.atom[ii].name]=atom_ref .get(refCell.atom[ii].name,0)+1
		matching_atoms=1
		for ii in atom_ref.keys():
			if not atom_self[ii] == atom_ref[ii]: same_atom=0
		if not matching_atoms:
			log.error("match: names of the elements in either cell are not matching (check element names), terminating...")
			exit(1)
		"================================================================================================================"					

		ref_displ_frac=np.array([refCell.atom[ii].position_frac for ii in range(refCell.N)])
		act_displ_frac=np.array([self.atom[ii].position_frac for ii in range(self.N)])
		ref_carth=pmut.toCarthesian(ref_displ_frac,refCell.lattice.a,refCell.lattice.b,refCell.lattice.c)
		act_carth=pmut.toCarthesian(act_displ_frac,self.lattice.a,self.lattice.b,self.lattice.c)

		reordering_atoms=np.zeros(refCell.N)
		search_range=[-1,0,1]
		#search_range=[-2,-1,0,1,2]
		shifts=refCell.N*[[0,0,0]]
		if(refCell.N>100):
			print("Matching two atomic structures, this may take some time")
			N10percent=int(refCell.N/10.)
		for ii in range(refCell.N):											#for all atoms in reference cell
			if refCell.N>100 and ii%N10percent==0: print(str((ii*10.)/N10percent)+"% ... ")
			dist=1000.												#set very large distance (in Angstroems)
			for jj in range(self.N):										#for every atom in the actual structure
				if self.atom[jj].name == refCell.atom[ii].name and self.atom[jj].coreshell == refCell.atom[ii].coreshell:	#but only for the atom of the same kind
					for aa in search_range:									#apply displacement in three directions in the search_range
						for bb in search_range:
							for cc in search_range:
								act_shift=aa*np.array(self.lattice.a)+bb*np.array(self.lattice.b)+cc*np.array(self.lattice.c)
								act_carth_shifted=act_carth[jj]+act_shift
								dist_act_ref=np.linalg.norm(act_carth_shifted-ref_carth[ii])
								if dist_act_ref<dist:
									dist=dist_act_ref
									#print "%3d"%ii, ": smallest distance so far: ", dist, ", for atom ", self.atom[jj].name, ", id: ", jj
									reordering_atoms[ii]=jj
									shifts[ii]=[aa,bb,cc]

		copyself=copy.deepcopy(self)
		for ii in range(self.N):
			self.atom[ii].position_frac=	[copyself.atom[int(reordering_atoms[ii])].position_frac[jj]+shifts[ii][jj] for jj in [0,1,2]]
			self.atom[ii].name=		copyself.atom[int(reordering_atoms[ii])].name
			self.atom[ii].charge=           copyself.atom[int(reordering_atoms[ii])].charge
			self.atom[ii].coreshell=        copyself.atom[int(reordering_atoms[ii])].coreshell
			self.atom[ii].wyckoff=          copyself.atom[int(reordering_atoms[ii])].wyckoff
			self.atom[ii].spin=             copyself.atom[int(reordering_atoms[ii])].spin
			
		print("Matching atoms: reordering of atoms: "+str( reordering_atoms)+"\n")
	
	def structureHash(self,distance=5):
		aux_cart=self.getArrayCarthesian();

		distance_dict={}
		Nextension=1
		
		res_atom=0.
		#res_latt=0.
		for ii in range(-Nextension,Nextension+1):
			for jj in range(-Nextension,Nextension+1):
				for kk in range(-Nextension,Nextension+1):
					for aa in range(len(aux_cart)):
						aux_pos=aux_cart[aa]+ii*np.array(self.lattice.a)+jj*np.array(self.lattice.b)+kk*np.array(self.lattice.c)
						aux_dist=np.linalg.norm(aux_cart[0]-aux_pos)
						if aux_dist <= distance:
							res_atom+=math.trunc(1000000*aux_dist)/1000000.

		print("Final sum for the structure with distance = ", distance, "from the 0-th atom is: ")
		#print "  Lattice: ", res_latt
		print("  Atoms  : %.6f"%res_atom)
							

		
