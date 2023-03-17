#!/usr/bin/env python
#-*- coding: utf-8 -*-

'''
Script para convertir una archivo sdf a pdb
'''
import os
import sys
from openbabel import openbabel


def SDFtoPDB(molecule):
	'''Convertir compuestos de formato con openbabel'''
	obConversion = openbabel.OBConversion()
	obConversion.SetInAndOutFormats("sdf", "pdb")
	mol = openbabel.OBMol()
	notatend = obConversion.ReadFile(mol, molecule) 
	while notatend:
		name = molecule.replace('.sdf', '')
		obConversion.WriteFile(mol,name+'.pdb')
		mol = openbabel.OBMol()
		notatend = obConversion.Read(mol)


SDFtoPDB(sys.argv[1])
	
	
	

	 
	
		
