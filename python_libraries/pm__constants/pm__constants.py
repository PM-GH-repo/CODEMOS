#!/usr/bin/env python

import sys
import numpy as np
import re
import os
import shutil
import glob
import copy

elementary_charge=1.602176565e-19

perovskite_typical_a = 4;	#typical perovskite lattice constants 4 Angstroems

def assignChargeTensor(name,coshe="core",chargeset="NOMINAL"):
	if chargeset == "PZT_SHIMADA" or chargeset == "PTO_SHIMADA" or chargeset == "NOMINAL":
		return assignCharge(name,coshe,chargeset)*np.diag((1,1,1))
	""" Born effective charges need to be assigned differently, as two identical atoms might have same
	charge depending on symmetry of its position within crystal """

def assignCharge(name,coshe="core",chargeset="NOMINAL"):

	if chargeset == "PZT_SHIMADA":
		if(name=="Pb" and coshe=="core" ): return  11.850027
		if(name=="Pb" and coshe=="shell"): return  -9.809102
		if(name=="Zr" and coshe=="core" ): return -10.199181
		if(name=="Zr" and coshe=="shell"): return  14.566555
		if(name=="Ti" and coshe=="core" ): return   1.019923
		if(name=="Ti" and coshe=="shell"): return   3.347449
		if(name=="O"  and coshe=="core" ): return   1.288949
		if(name=="O"  and coshe=="shell"): return  -3.425048

	if chargeset == "PTO_SHIMADA":
		if(name=="Pb" and coshe=="core" ): return   5.495849
		if(name=="Pb" and coshe=="shell"): return  -3.633216
		if(name=="Ti" and coshe=="core" ): return  19.369090
		if(name=="Ti" and coshe=="shell"): return -16.279515
		if(name=="O"  and coshe=="core" ): return   2.548431
		if(name=="O"  and coshe=="shell"): return  -4.199143

	if chargeset == "NOMINAL":
		""" Perovskite A-atoms """
		if(name=="Pb" and coshe=="core" ): return   2.00
		if(name=="Pb" and coshe=="shell"): return   0.00
		if(name=="Bi" and coshe=="core" ): return   3.00
		if(name=="Bi" and coshe=="shell"): return   0.00
		if(name=="Ba" and coshe=="core" ): return   2.00
		if(name=="Ba" and coshe=="shell"): return   0.00
		if(name=="Li" and coshe=="core" ): return   1.00
		if(name=="Li" and coshe=="shell"): return   0.00
		if(name=="K"  and coshe=="core" ): return   1.00
		if(name=="K"  and coshe=="shell"): return   0.00
		if(name=="Sr" and coshe=="core" ): return   2.00
		if(name=="Sr" and coshe=="shell"): return   0.00
		
		""" Perovskite B-atoms """
		if(name=="Zr" and coshe=="core" ): return   4.00
		if(name=="Zr" and coshe=="shell"): return   0.00
		if(name=="Ti" and coshe=="core" ): return   4.00
		if(name=="Ti" and coshe=="shell"): return   0.00
		if(name=="Fe" and coshe=="core" ): return   3.00
		if(name=="Fe" and coshe=="shell"): return   0.00
		if(name=="Nb" and coshe=="core" ): return   5.00
		if(name=="Nb" and coshe=="shell"): return   0.00
		if(name=="Ta" and coshe=="core" ): return   5.00
		if(name=="Ta" and coshe=="shell"): return   0.00
		if(name=="Cd" and coshe=="core" ): return   2.00
		if(name=="Cd" and coshe=="shell"): return   0.00
		
		""" Oxygen """
		if(name=="O"  and coshe=="core" ): return  -2.00
		if(name=="O"  and coshe=="shell"): return   0.00
		
		""" Others """
		if(name=="Sc" and coshe=="core" ): return   3.00
		if(name=="Sc" and coshe=="shell"): return   0.00
	
		""" GVS atoms """		
		if(name=="Ga" and coshe=="core" ): return  +4.00
		if(name=="Ga" and coshe=="shell"): return   0.00
		if(name=="V"  and coshe=="core" ): return  +3.00
		if(name=="V"  and coshe=="shell"): return   0.00
		if(name=="Mo" and coshe=="core" ): return  +3.00	#Allowed according to wiki
		if(name=="Mo" and coshe=="shell"): return   0.00
		if(name=="S"  and coshe=="core" ): return  -2.00
		if(name=="S"  and coshe=="shell"): return   0.00
		if(name=="Se" and coshe=="core" ): return  -2.00
		if(name=="Se" and coshe=="shell"): return   0.00


def getABO(name):
	if name == "Pb": return "A"
	if name == "Ba": return "A"
	if name == "Bi": return "A"
	if name == "Li": return "A"
	if name == "K" : return "A"
	if name == "Sr": return "A"

	if name == "Fe": return "B"
	if name == "Zr": return "B"
	if name == "Ti": return "B"
	if name == "Sc": return "B"
	if name == "Nb": return "B"
	if name == "Ta": return "B"
	if name == "Cd": return "B"

	if name == "O" : return "O"
	return None

def isAtom(name):
	""" Returns 1 if the parameter represent name of an element, 0 otherwise """
	if name == "Pb": return 1
	if name == "Ba": return 1
	if name == "Bi": return 1
	if name == "Sr": return 1
	if name == "Fe": return 1
	if name == "Zr": return 1
	if name == "Ti": return 1
	if name == "Sc": return 1
	if name == "Nb": return 1
	if name == "K" : return 1
	if name == "Li": return 1
	if name == "Ta": return 1
	if name == "O" : return 1
	if name == "Cd": return 1
	return 0
