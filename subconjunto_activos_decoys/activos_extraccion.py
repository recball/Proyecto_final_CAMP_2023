#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys

def extractActives(lista_picked):
	'''Extraer ligandos activos de un archivo .picked'''
	ligandos_activos = []
	with open(lista_picked, 'r') as listaPicked:
		for archivo in listaPicked:
			with open(archivo.strip(), 'r') as decoys_activos:
				ligandos_activos.append(decoys_activos.readlines()[0])
	
	with open('lista_activos.txt', 'a') as activos:
		for activo in ligandos_activos:
			activos.write(activo)
			



extractActives(sys.argv[1])
