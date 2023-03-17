#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri March 10 13:04:37 2023

@author: recball
"""

import pandas as pd 
import sys

def SMINA_SCORE(lista_log, output):
    '''Extrae las energias de los archivos .log obtenidos desde VINA/SMINA
    imput: archivo.txt con los archivos.log
    output'''
    nombres = list()
    scores =list()
    with open (lista_log, 'r') as f:
        for compuesto in f:
            nombres.append(compuesto.replace('.log\n', ''))
            print(compuesto)
            with open(compuesto.strip(), 'r') as logfile:
                for p, line in enumerate(logfile):
                    if p == 25:
                        score = line.split()[1]
                        scores.append(float(score))
        informacion = { 'nombre': nombres, 'SMINA vinardo': scores}
        data = pd.DataFrame(informacion)
        print(data.sort_values(by=['SMINA vinardo'], ascending=True))
        data.sort_values(by=['SMINA vinardo'], ascending=True).to_csv(output+'.csv', index=True)

if len(sys.argv) == 3:
	SMINA_SCORE(sys.argv[1], sys.argv[2])
