#!/usr/bin/env python

import glob
import time
import datetime
import sys
import os


def getRunningCPU():
	return len(glob.glob('./running_*'))


def wait_until_NCPU_less_than(NCPU,delay=10):
	running_cores=len(glob.glob('./running_*'))
	while running_cores >= NCPU:
		print(" Running %d processes, waiting until it is %d, sleeping ..."%(running_cores,NCPU))
		time.sleep(delay)
		running_cores=len(glob.glob('./running*'))

def wait_until_maximum_running_processes(maximum=0,delay=10):
	time.sleep(1)							#give the chance to the previously launched process to write its "running__" file. Sometimes this take some time
	running_processes=len(glob.glob('./running_*'))
	while running_processes > maximum:
		print(" %d processes running, waiting ..."%(running_processes,))
		time.sleep(delay)
		running_processes=len(glob.glob('./running*'))

def stamp():
	return str(datetime.datetime.now())

def timestamp(file,flag="a"):
	fout=open(file,flag)
	timestamp=stamp()
	fout.write("# %s\n" % timestamp)
	fout.close()

def recordcommand(file):
	fout=open(file,"a+")
	fout.write("%s\n"%' '.join(sys.argv))
	fout.close()

def printCommand(file="RESULT_command.dat"):
	f=open(file,"a+")
	f.write("%s\n"%str(sys.argv))
	f.close()
