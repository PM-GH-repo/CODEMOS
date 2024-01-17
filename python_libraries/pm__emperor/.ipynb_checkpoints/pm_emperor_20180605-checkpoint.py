#!/usr/bin/env python
import time
import glob
import os
import pm__task_management as pmtm
import hashlib

class Emperor:
	"""
	Class for managing project: plenty of simulations 
	"""

	def __init__(self,name="UNKNOWN_PROJECT",dir_OUTPUT="OUTPUT/",dir_RUNDIR="RUNDIR/",file_log="LOGFILE.dat",nPROC=1):
		self.name=name
		self.dir_OUTPUT=dir_OUTPUT
		self.dir_RUNDIR=dir_RUNDIR
		self.file_log=file_log
		self.simulations=[]
		self.dir_ROOT=os.getcwd()
		self.nPROC=nPROC
		self.strategy="from_first"
		self.delay_interval=10

		self.simulations_to_be_launched=[]

		for ii in [self.dir_OUTPUT,self.dir_RUNDIR]:
			if not os.path.exists(ii):os.system("mkdir "+ii)
		pmtm.timestamp(self.file_log)

	def __str__(self):
		aux_string=""
		return aux_string

	def addSimulation(	self,
				name="unknown_simulation",
				files_copy=[],
				files_string={},
				files_store=[],
				depend_on_sims=[],
				command="",
				hashbase=""):
		""" 	name: 		simulation name
			files_copy: 	array of paths (with respect to the parent directory) to files,
					which are supposed to be copied to the calculation directory
			files_string:	dictionary (filename:string) with strings for input files
			files_store:	array of filenames for files, which should be stored
			command:	path to the program which conducts calculations"""

		aux_sim=Simulation(name=self.name+"__"+name,
				files_copy=files_copy,
				files_string=files_string,
				files_store=files_store,
				depend_on_sims=depend_on_sims,
				command=command,
				project=self,
				hashbase=hashbase,
				identificator=len(self.simulations))

		self.simulations.append(aux_sim)

		return aux_sim


	def printAll(self):
		for ii in self.simulations:
			self.printSeparator()
			print "ID = ",ii.identificator
			print ii

	def printSeparator(self):
		print "===================================================================================================="
		
				
	def runSingleSimulation(self):
		a="none"
	
	def chooseSimulationToLaunch(self):
		#initiate array simulations_to_be_launched
		aux_array=[]
		for ii in self.simulations:
			if not ii.is_running and not ii.is_finished and not ii.is_identical_to:
				#was it calculated previously?
				aux_found_all=1
				aux_string="../"+self.dir_OUTPUT+"*"+ii.name+"*"
				for jj in ii.files_store:
					aux_outfiles=glob.glob(aux_string+jj)
					if not aux_outfiles:
						aux_found_all=0
						#print aux_string+jj+" not found in "+self.dir_OUTPUT
				if	aux_found_all:
					ii.is_finished="yes"
					ii.is_running=None


				if 	not ii.is_identical_to and not ii.is_finished:
					aux_array.append(ii)

		if self.strategy == "from_first":
			self.simulations_to_be_launched=aux_array
		elif self.strategy == "from_last":
			self.simulations_to_be_launched=aux_array[::-1]#invert the array 
		else:
			print("Unknown strategy for launching simulations, terminating ...")
			exit(1)
	
		
		print self.simulations_to_be_launched
		#for ii in self.simulations_to_be_launched:
		#	print ii



		sim=None
		for ii in self.simulations_to_be_launched:
			#does depend on unfinished simulation?
			aux_depend=None
			for jj in ii.depend_on_sims:
				if jj in self.simulations_to_be_launched:				#This simulation depends on other simulation, which needs tobe evaluated first
					aux_depend="yes"
					break;
			if 	not aux_depend and not ii.is_identical_to and sim == None:		#Not depend on not-calculated task, is not same as any from the group, no simulation was yet selected for calculataion
				sim=ii									#Pick actual simulation for calculation
			if	not aux_depend and not ii.is_identical_to and ii.prerequisit_of_sims:	#Not depend on not-calculated task, is not same as any from the group, is prerequisite for other simulations
				sim=ii									#Select actual simulation for calculation (can overwrite already picked simulation, which is not prerequisite of other simulations
				continue								#And continue, because this cannot be owerwritten
		if sim:		print("Chosen simulation to be launched is: "+sim.name)
		else:		print("No simulation is selected for simulation")
		
		return sim

	def runAllSimulations(self,nPROC=0,strategy="from_first"):
		local_nPROC=(nPROC if nPROC>0 else self.nPROC)						#override choice made at creation of the Emperor class	
		self.strategy=strategy
		
		os.chdir(self.dir_RUNDIR)								#change directory to RUNDIR
		
		chosen_simulation=self.chooseSimulationToLaunch()					
		while(chosen_simulation):
			time.sleep(1)									#waiting to ensure that the "running" locks are properly written
			pmtm.wait_until_maximum_running_processes(maximum=1,delay=self.delay_interval)		#waiting for enough processors
		
			chosen_simulation.createCalculationDir()
			print("SIMULATION TO BE LAUNCHED:\n"+str(chosen_simulation))
			os.chdir(chosen_simulation.directory)
			aux_command=r"( touch ../running_"+chosen_simulation.name+" ; nohup "+chosen_simulation.command+" 2>out ; "		#2>/dev/null
			for jj in chosen_simulation.files_store:
				aux_command+="cp "+jj+" ../../"+self.dir_OUTPUT+chosen_simulation.name+"__"+jj+";  "
			aux_command+=" rm ../running_"+chosen_simulation.name+"; rm ./* ; rmdir ../"+chosen_simulation.name+"; ) &"
			print("LAUNCH COMMAND:")
			self.printSeparator()
			print(aux_command)

			chosen_simulation.is_running="yes"
			os.system(aux_command)

			print("")

			os.chdir("..")
			chosen_simulation=self.chooseSimulationToLaunch()
			while self.simulations_to_be_launched and not chosen_simulation:		#no simulation is selected although there are simulations to be calculated:
				#HERE I SHOUD CHECK THAT SOME SIMULATION IS RUNNING AND EXIT IF NOT!!!	#the only plausible reason for this is that simulations, on which all the waiting simulatinons depend, are still running.
				time.sleep(self.delay_interval)						#waiting to ensure that the "running" locks are properly written
				print("Waiting to finis simulation, on which all other simulations depend...")
				chosen_simulation=self.chooseSimulationToLaunch()

		pmtm.wait_until_maximum_running_processes(maximum=0,delay=self.delay_interval)

		os.chdir("..")

	def deleteUnusedSimulations(self):
		a="none"

			
class Simulation:	
	def __init__(self,	name="unknown_simulation",
				files_copy=[],
				files_string={},
				files_store=[],
				depend_on_sims=[],
				command="",
				project="",
				hashbase="",
				identificator=0):
		#from argument
		self.name=			name
		self.project=			project
		self.files_copy=		files_copy
		self.files_string=		files_string
		self.files_store=		files_store
		self.depend_on_sims=		depend_on_sims
		self.command=			command
		self.identificator=		identificator
		
		#derived
		self.directory=			self.name
		self.is_finished=		None
		self.is_running=		None
		self.prerequisit_of_sims=	[]

		aux_string=hashbase
		if not aux_string:
			for ii in self.files_string.keys():
				aux_string+=self.files_string[ii]
			self.hashstring=hashlib.sha224(aux_string).hexdigest()

		self.is_identical_to=	None
		for ii in project.simulations:
			if ii.hashstring==self.hashstring:
				self.is_identical_to=ii
				self.is_finished="yes"
				break

		for ii in self.depend_on_sims:
			print "i am depending on :",ii.name
			ii.prerequisit_of_sims.append(self)

	def createCalculationDir(self):
		#create own directory and create/copy input files
		if not os.path.exists(self.name):os.mkdir(self.name)
		os.chdir(self.directory)
		for ii in self.files_string.keys():		#Write all files which are provided as strings
			f=open(ii,"w")
			f.write(self.files_string[ii])
			f.close()
		
		for ii in self.files_copy:			#Copy all files (path is based on parent directory)
			os.system("cp ../../"+ii+" .")
		os.chdir("..")

	def deleteCalculationDir(self):
		os.system("rm -f "+self.name)

	def __str__(self):
		aux_string="====================================================================================================\n"+\
		"Simulation "+self.name+": basic information\n"+\
		"name               : "+self.name+"\n"+\
		"identificator      : "+str(self.identificator)+"\n"+\
		"project            : "+self.project.name+"\n"+\
		"input_copy         : "+str(self.files_copy)+"\n"+\
		"input_strings      : "+str(self.files_string.keys())+"\n"+\
		"store              : "+str(self.files_store)+"\n"+\
		"depend_on_sims     : "+str([self.depend_on_sims[ii].name for ii in range(len(self.depend_on_sims))])+"\n"+\
		"prerequisit_of_sims: "+str([self.prerequisit_of_sims[ii].name for ii in range(len(self.prerequisit_of_sims))])+"\n"+\
		"command            : "+self.command+"\n"+\
		"is_finished        : "+("yes" if self.is_finished else "no")+"\n"+\
		"is_running         : "+("yes" if self.is_running else "no")+"\n"+\
		"is_identical_to    : "+(self.is_identical_to.name if self.is_identical_to else "no")+"\n"+\
		"directory          : "+self.directory+"\n"+\
		"hash               : "+self.hashstring+"\n"

		return aux_string
