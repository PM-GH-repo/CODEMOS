#!/usr/bin/env python
import time
import glob
import os
import pm__task_management as pmtm
import hashlib
import datetime

#History of Changes
#10.09.2018: changed delay(1) to delay(0.5)

class Emperor:
	"""
	Class for managing project: plenty of simulations 
	"""

	def __init__(	self,
			name="UNKNOWN_PROJECT",
			dirs={"dir_OUTPUT":"OUTPUT/","dir_RUNDIR":"RUNDIR/","dir_PASS":"PASS/",},
			files={"file_log":"LOGFILE.dat",},
			nPROC=1):
		self.name=name
		self.dir_OUTPUT=	dirs["dir_OUTPUT"]
		self.dir_RUNDIR=	dirs["dir_RUNDIR"]
		self.dir_PASS=		dirs["dir_PASS"]
		self.file_log=		files["file_log"]
		self.simulations=	[]
		self.dir_ROOT=		os.getcwd()
		self.nPROC=		nPROC
		self.strategy=		"from_first"
		self.delay_interval=	10		#seconds
		self.are_running=	[]

		self.last_simulation=None
		self.simulation_counter=0

		self.simulations_to_be_launched=self.simulations

		for ii in [self.dir_OUTPUT,self.dir_RUNDIR,self.dir_PASS]:
			if not os.path.exists(ii):
				print "Creating directory ",ii
				os.system("mkdir "+ii)
			
		pmtm.timestamp(self.file_log)

	def __str__(self):
		aux_string=""
		return aux_string

	def addSimulation	(	
				self,
				name=["unknown_simulation"],
				files={	"files_copy":[],
					"files_string":{},
					"files_store":[],
					"files_pass":[],
					"files_adopt":{},
				},
				depends_on_sims=[],
				command="",
				hashbase=""
				):
		""" 	name: 		array of simulation names, which will be _ separated
			files_copy: 	array of paths (with respect to the parent directory) to files,
					which are supposed to be copied to the calculation directory
			files_string:	dictionary (filename:string) with strings for input files
			files_store:	array of filenames for files, which should be stored
			files_pass:	files which will be passed to following simulations (copied to directory PASS)
			files_adopt:	files which will be adopted/taken from previous simulations (copied from directory PASS)
			command:	path to the program which conducts calculations"""

		aux_name=self.name+"_"
		if type(name) is list:
			for ii in name:
				aux_name+=("_"+ii)
		else:		
			aux_name+=name

		aux_sim=Simulation(name=aux_name,files=files,
				depends_on_sims=depends_on_sims,
				command=command,
				project=self,
				hashbase=hashbase,
				identificator=self.simulation_counter,
				emperor=self)

		self.simulations.append(aux_sim)
		
		self.last_simulation=aux_sim			#update value for the last simulation
		self.simulation_counter+=1			#increment simulation counter

		return aux_sim


	def printAll(self):
		for ii in self.simulations:
			self.printSeparator()
			print "ID = ",ii.identificator
			print ii

	def printSeparator(self):
		print "===================================================================================================="
		
				
	def runSingleSimulation(self):
		print ("runSingleSimulation is not yet implemented, terminating")
		exit(1)
	
	def chooseSimulationToLaunch(self):
		aux_array=[]

		#print "simulations_to_be_launched",self.simulations_to_be_launched
		for ii in self.simulations_to_be_launched+self.are_running:
			aux_found_all=1
			aux_string="../"+self.dir_OUTPUT+"*"+ii.name+"*"
			for jj in ii.files_store:
				aux_outfiles=glob.glob(aux_string+jj)
				if not aux_outfiles:
					aux_found_all=0
			#print "aux_found_all",aux_found_all
			if	aux_found_all:
				ii.is_finished="yes"
				ii.is_running="no"
				self.are_running=[x for x in self.are_running if x != ii]
			#print "is_running",ii.is_running
			#print "is_finished",ii.is_finished
			#print "is_identical_to",ii.is_identical_to

			if ii.is_running=="no" and ii.is_finished=="no" and not ii.is_identical_to:
				aux_array.append(ii)
		#print "aux_array",aux_array
		if self.strategy == "from_first":
			self.simulations_to_be_launched=aux_array
		elif self.strategy == "from_last":
			self.simulations_to_be_launched=aux_array[::-1]#invert the array 
		else:
			print("Unknown strategy for launching simulations, terminating ...")
			exit(1)
	
		
		#print "DEBUG: simulations to be launched are: ",self.simulations_to_be_launched

		sim=None
		for ii in self.simulations_to_be_launched:
			#does depend on unfinished simulation?
			aux_depends_on_unfinished="no"
			for jj in ii.depends_on_sims:
				if jj.is_finished == "no":							#This simulation depends on other simulation, which needs tobe evaluated first
					aux_depends_on_unfinished="yes"
			if aux_depends_on_unfinished == "yes": 
				#print ii.name," depend on unfinished ",jj.name
				pass
			if 	not aux_depends_on_unfinished == "yes" and sim == None:				#Does not depend on not-calculated task, no simulation was yet selected for calculataion
				#print "Simulation ",ii.name," does not depend on unfinished and no simulation was yet chosen, I pick it!"
				sim=ii									#Pick actual simulation for calculation
			if	not aux_depends_on_unfinished == "yes" and ii.is_prerequisite_for:		#Does not depend on not-calculated task, is prerequisite for other simulations
				#print "Simulation ",ii.name," does not depend on unfinished and moreover it is prerequisite for other simulations, I pick it!"
				sim=ii									#Select actual simulation for calculation (can overwrite already picked simulation, which is not prerequisite of other simulations
				break									#And continue, because this cannot be owerwritten
		if sim:		print("Simulation chosen to be launched is: "+sim.name)
		else:		print("No simulation is selected for launch")
		
		return sim

	def runAllSimulations(self,nPROC=0,strategy="from_first",functions=[]):
		



		local_nPROC=(nPROC if nPROC>0 else self.nPROC)						#override choice made at creation of the Emperor class	
		self.strategy=strategy
		
		os.chdir(self.dir_RUNDIR)								#change directory to RUNDIR
		
		chosen_simulation=self.chooseSimulationToLaunch()					
		while(chosen_simulation):
			time.sleep(0.5)									#waiting to ensure that the "running" locks are properly written
			pmtm.wait_until_maximum_running_processes(maximum=local_nPROC-1,delay=self.delay_interval)		#waiting for enough processors
		
			chosen_simulation.createCalculationDir()
			#print("SIMULATION TO BE LAUNCHED:\n"+str(chosen_simulation))
			
			os.chdir(chosen_simulation.rundir)
			aux_command=r"( touch ../running_"+chosen_simulation.name+" ; nohup "+chosen_simulation.command+" 2>out ; "		#2>/dev/null
			for jj in chosen_simulation.files_store:
				aux_command+="cp "+jj+" ../../"+self.dir_OUTPUT+chosen_simulation.name+"__"+jj+";  "
			for jj in chosen_simulation.files_pass:
				aux_command+="cp "+jj+" ../../"+self.dir_PASS+chosen_simulation.name+"__"+jj+";  "
			aux_command+=" rm ../running_"+chosen_simulation.name+"; rm ./* ; rmdir ../"+chosen_simulation.name+"; ) &"
			#aux_command+=" rm ../running_"+chosen_simulation.name+"; ) &"
			
			#print("LAUNCH COMMAND:\n"+aux_command
			print("Launching simulation: "+chosen_simulation.name)

			chosen_simulation.is_running="yes"
			self.are_running.append(chosen_simulation)
			
			#exit(1)
			os.system(aux_command)

			print("")

			os.chdir("..")
			chosen_simulation=self.chooseSimulationToLaunch()
			while self.simulations_to_be_launched and not chosen_simulation:		#no simulation is selected although there are simulations to be calculated:
				#HERE I SHOUD CHECK THAT SOME SIMULATION IS RUNNING AND EXIT IF NOT!!!	#the only plausible reason for this is that simulations, on which all the waiting simulatinons depend, are still running.
				time.sleep(self.delay_interval)						#waiting to ensure that the "running" locks are properly written
				print("Waiting to finish simulation, which is prerequisite for other simulations ...")
				chosen_simulation=self.chooseSimulationToLaunch()
			for ii in functions:
				print "=================================",os.getcwd()
				#print "function = ",str(functions),str(ii)
				ii

		pmtm.wait_until_maximum_running_processes(maximum=0,delay=self.delay_interval)

		os.chdir("..")
		
		print("All simulations are finished")

	def deleteUnusedSimulations(self):
		a="none"

			
class Simulation:	
	def __init__(self,	name="unknown_simulation",
				files={},
				depends_on_sims=[],
				command="",
				project="",
				hashbase="",
				identificator="undefined",
				emperor=None):
		"""
		Class Simulations summarizes all information about a given simulation.\n
		:param name: name of the simulation. The string is then used to store the data\n
		:type name: str
		"""
		def unify_array(arg):
#			print "JSEM TU",arg
			if type(arg) is list:
#				print "jsem list"
				if not (arg == [] or arg == [None,] or arg == ["",]):
#					print "-------not",arg
					return arg

				else:
#					print "-------yes",arg
					return []
			elif type(arg) is dict:
#				print "jsem dict"
				to_be_deleted=[]
				for ii in arg:
					if arg[ii] == None or arg[ii] == "" or arg[ii] == False: to_be_deleted.append(ii)
				for ii in to_be_deleted:
					del[arg[ii]]
				return arg
			else:
				print "jsem neco jineho, wrong type passed, terminating..."
				exit(1)
			

		#from argument
		self.name=			name
		self.project=			project
		self.files_copy=		unify_array(files["files_copy"])
		self.files_string=		unify_array(files["files_string"])
		self.files_store=		unify_array(files["files_store"])
		self.files_pass=		unify_array(files["files_pass"])
		self.files_adopt=		unify_array(files["files_adopt"])
#		self.files_adopt=		files["files_adopt"]
		self.depends_on_sims=		unify_array(depends_on_sims)
#		if not (depends_on_sims == [] or depends_on_sims == [None] or depends_on_sims == ["",]): self.depends_on_sims=		depends_on_sims
#		else: self.depends_on_sims = []
		self.command=			command
		self.identificator=		identificator
		self.emperor=			emperor
		
		#derived
		self.rundir=			self.name
		self.is_finished=		"no"
		self.is_running=		"no"
		self.is_prerequisite_for=	[]

		aux_string=hashbase
		if not aux_string:
			""" based on input strings """
			for ii in self.files_string.keys():
				aux_string+=self.files_string[ii]
			""" based on copied files """
			for ii in self.files_copy:
				if not os.path.exists(ii[0]):
					print "Warning, file indexed for copy does not yet exist. This can make sence if it is outcome of a simulation, on which present simulation depend"
					aux_string+=str(datetime.datetime.now()) #Because the file is not know, add a string, which is unique (based on date)
				else:
					with open(ii[0], 'r') as infile:
						aux_cnt=0
						for line in infile:
							if aux_cnt<100 and line:
								aux_string+=line
							else:
								break

			self.hashstring=hashlib.sha224(aux_string).hexdigest()

		self.is_identical_to=	None
		for ii in project.simulations:
			if ii.hashstring==self.hashstring:
				self.is_identical_to=ii
				break

		#print self.depends_on_sims
		for ii in self.depends_on_sims:
			#print self.depends_on_sims
			#print "i am depending on :",ii.name
			ii.is_prerequisite_for.append(self)

	def getResultName(self):
		"""
		Returns instance of the Simulation class, which contains relevant results. Sometimes, the simulation is
		identified as identical to other simulation. In such case, the simulation is not launched, but results 
		from the identical simulations are utilized.\n
		"""
		sim=self
		while sim.is_identical_to:
			sim=sim.is_identical_to

		return sim

	def createCalculationDir(self):
		#create own directory and create/copy input files
		if not os.path.exists(self.name):os.mkdir(self.name)
		os.chdir(self.rundir)
		for ii in self.files_string.keys():		#Write all files which are provided as strings
			f=open(ii,"w")
			f.write(self.files_string[ii])
			f.close()
		
		for ii in self.files_copy:			#Copy all files (path is based on parent directory)
			if type(ii) is list or type(ii) is tuple:				#two values are provided meaning that the file must be renamed
				if not len(ii) == 2: 
					print ("length of array for file rename must be 2, terminating...")
					exit(1)
				file_origin=ii[0]
				file_rename=ii[1]
			else:
				file_origin=ii
				file_rename=""
			if os.path.exists("../../"+file_origin):
				os.system("cp ../../"+file_origin+" ./"+file_rename)
			else: 
				print("File ../../"+file_origin+" was not found, cannot copy ...")

		for ii in self.files_adopt.keys():
			aux_string="../../"+self.emperor.dir_PASS+self.files_adopt[ii].name+"__"+ii
			if os.path.exists(aux_string):
				print "Copying file "+aux_string
				os.system("cp "+aux_string+" ./"+ii)
			else:
				print "File "+aux_string+" not found, PROBLEM !!!!!!!!!!!!!!!!!!!"
				exit(1)


		os.chdir("..")

	def deleteCalculationDir(self):
		os.system("rm -f "+self.name)

	def __str__(self):
		aux_string=""+\
		"                     SIMULATION BASIC INFORMATION\n"+\
		"---------------------------------------------------------------------------------------\n"+\
		"name               : "+self.name+"\n"+\
		"identificator      : "+str(self.identificator)+"\n"+\
		"project            : "+self.project.name+"\n"+\
		"files_copy         : "+str(self.files_copy)+"\n"+\
		"files_strings      : "+str(self.files_string.keys())+"\n"+\
		"files_store        : "+str(self.files_store)+"\n"+\
		"files_pass         : "+str(self.files_pass)+"\n"+\
		"files_adopt        : "
		for ii in self.files_adopt.keys():
			aux_string+="\'"+ii+"\'"+" from simulation "+str(self.files_adopt[ii].name)+"\n"
		if self.depends_on_sims:	aux_string+="depends_on_sims    : None\n"
		else:				aux_string+="depends_on_sims    : "+str([self.depends_on_sims[ii].name for ii in range(len(self.depends_on_sims))])+"\n"
		aux_string+=\
		"is_prerequisite_for: "+str([self.is_prerequisite_for[ii].name for ii in range(len(self.is_prerequisite_for))])+"\n"+\
		"command            : "+self.command+"\n"+\
		"is_finished        : "+self.is_finished+"\n"+\
		"is_running         : "+self.is_running+"\n"+\
		"is_identical_to    : "+(self.is_identical_to.name if self.is_identical_to else "no")+"\n"+\
		"rundir             : "+self.rundir+"\n"+\
		"hash               : "+self.hashstring+"\n"+\
		"emperor            : "+str(self.emperor)+"\n"

		return aux_string
