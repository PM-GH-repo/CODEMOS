#!/usr/bin/env python

import sys
import numpy as np
import re
import os
import shutil
import glob
import copy


def parse_properties_file(filename):
	class property:
		def __init__(self,line,metainfo):
			self.metainfo=metainfo
			#print "Constructor line: %s"%line
			#line_split=re.split('\W+',line)
			line_split=line.split()
			self.name=line_split[0]
			print(line_split)
			if line_split[1][0]=='~':
				self.value=float(line_split[1][1:])
				self.metainfo["approximate"]="yes"
				#print self.metainfo["approximate"]
			else:
				print(line_split)
				self.value=float(line_split[1])
				self.metainfo["approximate"]="no"

		def __str__(self):
			aux_approx=""
			aux_aux1=""
			aux_aux2=""
			if self.metainfo["approximate"]=="yes":	aux_approx="~"
			if self.metainfo["aux1"]!="--":	aux_aux1=self.metainfo["aux1"]
			if self.metainfo["aux2"]!="--":	aux_aux2=self.metainfo["aux2"]

			aux_string=	" %10s"%self.name+\
					" %1s%14.6f"%(aux_approx,self.value)+\
					" %8s"%self.metainfo["material"]+\
					" %8s"%self.metainfo["symmetry"]+\
					" %25s"%self.metainfo["reference"]+\
					" %s"%aux_aux1+\
					" %s"%aux_aux2
			return aux_string
		def writeLatex(self,ndecimal=0):
			return r"%.*f"%(ndecimal,self.value)+"\cite{%s}"%self.metainfo["reference"]

	properties=[]
	dict_head={}
	with open(filename) as f:
		for line in f:
			#line_split=re.split('\W+',line)
			line_split=line.split()
			if line_split:
				#print "line_split=",line_split
				if line_split[0][0]=="#":		#skip comments
					continue	
				if line_split[0] == "HEAD":		#read header
					dict_head={}
					print("READING HEADER")
					dict_head["material"]=		line_split[1]
					dict_head["method"]=		line_split[2]
					dict_head["temperature"]=	line_split[3]
					dict_head["symmetry"]=		line_split[4]
					dict_head["domain"]=		line_split[5]
					dict_head["crystal"]=		line_split[6]
					dict_head["source"]=		line_split[7]
					dict_head["relevance"]=		line_split[8]
					dict_head["aux1"]=		line_split[9]
					dict_head["aux2"]=		line_split[10]
					dict_head["reference"]=		line_split[11]
					continue
				properties.append(property(line=line,metainfo=dict_head))
			
	for ii in range(len(properties)):
		print(str(properties[ii]))

	f.close()

	return properties

def find_line_with_pattern(filename,pattern,whichline='last'):
	loc_pattern=re.compile(pattern)
	
	if whichline == 'first':
		with open(filename) as f:
			for line in f:
				if loc_pattern.search(line):
					f.close()
					return line
		f.close()
		return
	if whichline == 'last':
		#print "searching pattern: ", pattern
		lastline=""
		with open(filename) as f:
			for line in f:
			#	print "line=%s" % line
				if loc_pattern.search(line):
					lastline=copy.deepcopy(line)
#					print "=================== pattern found: ", lastline.strip()
#				print line.strip(),"LASTLINE: ",lastline.strip()
		f.close()
		return lastline
	if whichline == 'all':
		lines=[]
		with open(filename) as f:
			for line in f:
				if loc_pattern.search(line):
					lines.append(copy.deepcopy(line))
#					print "=================== pattern found: ", lastline.strip()
#				print line.strip(),"LASTLINE: ",lastline.strip()
		f.close()
		return lines
	if isinstance(int(whichline),int) == True:
		nn=int(whichline)
		counter=0
		with open(filename) as f:
			for line in f:
			#	print "line=%s" % line
				if loc_pattern.search(line):
					counter+=1
					if counter == nn:
						f.close()
						return line
		f.close()
		return
	
	return 

def find_line_after_pattern(filename,pattern,whichline='last',pluslines=0,nlines=1):
	loc_pattern=re.compile(pattern)
	
	if whichline == 'last':
		#print "searching pattern: ", pattern
		lastline_single=""
		lastline_multiple=[]
		fin=open(filename,"r")
		line=fin.readline()
		while line != "":
			if loc_pattern.search(line):
				for ii in range(pluslines):
					line=fin.readline()
				if nlines==1:
					lastline_single=copy.deepcopy(line)
				if nlines>1:
					lastline_multiple=[]
					lastline_multiple.append(copy.deepcopy(line))
					for jj in range(nlines-1):
						line=fin.readline()
						lastline_multiple.append(copy.deepcopy(line))


			line=fin.readline()
		fin.close()
		
	if nlines==1:return lastline_single
	if nlines>1: return lastline_multiple
	

def replace_in_file(file_template,file_final,pattern,subst):
#	print "input file = %s" % file_template
#	print "output file = %s" % file_final 
#	print "pattern = %s" % pattern
#	print "substitute with = %s" % subst

	filenames_are_equal=0
	if file_final == file_final:
		#print("Same files, make backup copy")
		filenames_are_equal=1
		shutil.copy(file_final,"tmp_final.txt")
		file_final='tmp_final.txt'

	fin = open(file_template)
	fout = open(file_final, "wt")
	for line in fin:
#		print line.strip(), "subst = ",subst
		line =  (re.sub(pattern, subst, line))
		fout.write(line)
#		print line.strip(), "subst = ",subst
	fin.close()
	fout.close()

	if filenames_are_equal:
		shutil.move('tmp_final.txt',file_template)

def wait_until_NCPU_less_or_equal(NCPU):
	running_cores=len(glob.glob('./running*'))
	while running_cores > NCPU:
		print("Max number of processes reached, sleeping ...")
		time.sleep(10)
		running_cores=len(glob.glob('./running*'))


"""read data from file *.arr ----------------------------------------------------------------------"""
def read_arr(input_file):
	#set default values
	X=1; Y=1; Z=1; DIM=3; PDIM=1;

	fin=open(input_file);
	#read header
	line=re.split('#| |\t|\n',fin.next())
	line=filter(None,line)
	if re.search('SIM',line[0]):
		#simulation timestamp found on first line, read dimensions of the box on next line
		line=re.split('#| |\t|\n',fin.next())
		line=filter(None,line)

	#print line

	#TODO: generalize to arbitrary dimension
	pointer=0
	DIM=int(line[pointer])
	pointer+=1
	if(DIM >= 1):
		X=int(line[pointer])
		pointer+=1
	if(DIM >= 2):
		Y=int(line[pointer])
		pointer+=1
	if(DIM >= 3):
		Z=int(line[pointer])
		pointer+=1
	PDIM=int(line[pointer])
	fin.close()

	print("Reading file \"%s\", DIM=%d, X=%d, Y=%d, Z=%d, PDIM=%d ..." % (input_file,DIM,X,Y,Z,PDIM))


	P_from_file=np.loadtxt(input_file).reshape(Z,X,Y,PDIM)		#reshape to Z,X,Y (format of *.arr)
	P_resh=np.zeros((X,Y,Z,PDIM))
	for ii in xrange(X):
		for jj in xrange(Y):
			for kk in xrange(Z):
				for pp in xrange(PDIM):		#reshape to "expected order" X,Y,Z 
					P_resh[ii][jj][kk][pp]=P_from_file[kk][ii][jj][pp]
	Px=np.zeros((X,Y,Z))					#getting Px,Py,Pz to separate arrays
	Py=np.zeros((X,Y,Z))
	Pz=np.zeros((X,Y,Z))
	for ii in xrange(X):
		for jj in xrange(Y):
			for kk in xrange(Z):
				if PDIM >=1:
					Px[ii][jj][kk]=P_resh[ii][jj][kk][0]
				if PDIM >=2:
					Py[ii][jj][kk]=P_resh[ii][jj][kk][1]
				if PDIM >=3:
					Pz[ii][jj][kk]=P_resh[ii][jj][kk][2]

	return DIM,X,Y,Z,PDIM,Px,Py,Pz

"""save data to *.arr file ------------------------------------------------------------------------"""
def write_arr(output_file="save_array.arr",DIM=1,X=1,Y=1,Z=1,PDIM=1,Px=[1.0],Py=None,Pz=None):
	print("Writing file \"%s\"" % output_file)
	if X==None:
		X=1
	if Y==None:
		Y=1
	if Z==None:
		Z=1

	fout=open(output_file,'w');
	if DIM == 1:
		fout.write("#1\t%d\t%d\n"%(X,PDIM))
	if DIM == 2:
		fout.write("#2\t%d\t%d\t%d\n"%(X,Y,PDIM))
	if DIM == 3:
		fout.write("#3\t%d\t%d\t%d\t%d\n"%(X,Y,Z,PDIM))

	if DIM == 1:
		if PDIM==1:
			for xx in xrange(X):
				fout.write(" % .5f\n" % (Px[xx]))
		if PDIM==2:
			for xx in xrange(X):
				fout.write(" % .5f % .5f\n" % (Px[xx],Py[xx]))
		if PDIM==3:
			for xx in xrange(X):
				fout.write(" % .5f % .5f % .5f\n" % (Px[xx],Py[xx],Pz[xx]))
	if DIM == 2:
		if PDIM==1:
			for xx in xrange(X):
				for yy in xrange(Y):
					fout.write(" % .5f" % (Px[xx,yy]))
				fout.write("\n")
		if PDIM==2:
			for xx in xrange(X):
				for yy in xrange(Y):
					fout.write(" % .5f % .5f" % (Px[xx,yy],Py[xx,yy]))
				fout.write("\n")
		if PDIM==3:
			for xx in xrange(X):
				for yy in xrange(Y):
					fout.write(" % .5f % .5f % .5f" % (Px[xx,yy],Py[xx,yy],Pz[xx,yy]))
				fout.write("\n")
	if DIM == 3:
		if PDIM==1:
			for zz in xrange(Z):
				for xx in xrange(X):
					for yy in xrange(Y):
						fout.write(" % .5f" % (Px[xx,yy,zz]))
					fout.write("\n")
				fout.write("\n")
		if PDIM==2:
			for zz in xrange(Z):
				for xx in xrange(X):
					for yy in xrange(Y):
						fout.write(" % .5f % .5f" % (Px[xx,yy,zz],Py[xx,yy,zz]))
					fout.write("\n")
				fout.write("\n")
		if PDIM==3:
			for zz in xrange(Z):
				for xx in xrange(X):
					for yy in xrange(Y):
						fout.write(" % .5f % .5f % .5f" % (Px[xx,yy,zz],Py[xx,yy,zz],Pz[xx,yy,zz]))
					fout.write("\n")
				fout.write("\n")
	fout.close()
