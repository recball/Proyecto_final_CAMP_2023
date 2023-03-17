#!/usr/bin/env python
# coding: utf-8

# In[301]:


#validacion de la metodología de docking: construccion de una curva ROC
'''Descarga de archivos predeterminados de compuestos activos y decoys en formato smile, los archivos de formato .picked
 y busqueda de el ligando activo utilizado para la construccion de los decoys'''


# In[1]:


import csv
import pandas as pd
import os 
from openbabel import openbabel
import numpy as np
import random 
import subprocess


# In[303]:


os.system('wget https://dude.docking.org//targets/cdk2/dudgen_clustered/search')


# In[298]:


picked = []
with open('search.txt', 'r') as csvfile:
    datareader = csv.reader(csvfile)
    for row in datareader:
        if row:
            if row[0].endswith('.picked</a></li>'):
                picked.append(row)


# In[297]:


picked_limpios = [(picked[p][0].replace('    <li><img src="https://dude.docking.org//icons/text.gif"> <a href="','')).replace('</a></li>','') for p,picked_file in enumerate(picked)]
picked_limpios_completo = [picked_limpios[p].split('">')[0] for p,picked in enumerate(picked_limpios)]
picked_limpios_completo


# In[295]:


for picked in picked_limpios_completo:
    os.system(f'wget {picked}')


# In[305]:


# El numero de moleculas a generar y luego largar el docking es demasiado grande, por lo que para esta prueba piloto
# seleccionaré un subconjunto de ligandos: 50 ligandos activos y 2500 decoys, seleccionados de manera aleatoria

picked_seleccionados_al_azar= [random.choice(picked_limpios_completo) for numero_ligandos in range(50)]
picked_seleccionados_al_azar


# In[306]:


#Abrir los archivos para extraer ligandos, y obtener una lista de los compuestos activos para la quinasa
ligandos_activos_smile = []
ligandos_activos_ID = []
for picked in picked_seleccionados_al_azar:
    picked = picked.replace('https://dude.docking.org//targets/cdk2/dudgen_clustered/search/','')
    with open(picked, 'r') as contenido:
        for p,row in enumerate(contenido):
            if p==0:
                ligandos_activos_smile.append(row.split('\t')[1])
                ligandos_activos_ID.append(row.split('\t')[2])


# In[307]:


#Extraccion de los smiles de los decoys desde cada archivo .picked
decoys = [picked_file.replace('https://dude.docking.org//targets/cdk2/dudgen_clustered/search/','') for picked_file in picked_seleccionados_al_azar]
smiles_decoys = []
for decoy in decoys:
    with open(decoy, 'r') as decoys_smiles:
        for p,row in enumerate(decoys_smiles):
            if p!=0:
                smiles_decoys.append(row)


# In[313]:


#generar las conformaciones desde los smiles
def conformeros3D(lista):
    '''Generar estructuras 3D a partir de smiles en una lista utilizando openbabel'''
    decoys_smile = []
    for picked in decoys:
        with open(picked, 'r') as picked_file, open('smiles_fail.txt', 'a') as smile_fail,open('lista_completa.txt', 'a') as nombres:
            for p, smile in enumerate(picked_file): # ej: COc1ccc(-c2cc(=O)c3c(O)cc(O)cc3o2)cc1  ZINC ID
                    ID_ZINC = smile.strip().split('\t')[1] # Modificar de acuerdo al separador utilizado
                    ID_sdf_zinc = '-O'+ID_ZINC+'.sdf'
                    print(smile.strip().split('\t')[0])
                    fails=0
                    smile = '-:'+smile.strip()
                    try:
                        procesamiento_obabel = subprocess.run(['obabel', '-:',smile, '-osdf', ID_sdf_zinc, '--gen3D'], stdout = subprocess.PIPE,stderr=subprocess.STDOUT, text=True, timeout= 60)
                        for line in procesamiento_obabel.stdout:
                            fails +=1
                        if fails > 21:
                            print(f'La molécula {ID_ZINC} tuvo un error de procesamiento con Openbabel, verificar en el archivo, n = {p}')   
                            smile_fail.write(smile.strip().split('\t')[0]+','+'error quimico'+'\n')
                        else:
                            print(f'La molecula {ID_ZINC} fue correctamente procesada con Openbabel, verificar en el archivo, n = {p}')
                    except subprocess.TimeoutExpired:
                        print(f'La molecula {ID_ZINC} no corrió adecuadamente y el proceso fue detenido abruptamente, n = {p}')
                        nombres.write(ID_ZINC+','+'error de ejecución'+'\n')


# In[314]:


conformeros3D(smiles_decoys)


# In[ ]:


# 47 Moleculas fallaron en generar conformaciones 


# In[ ]:


#Para cada molecula minimizar y graficar las energias de moleculas minimizadas vs no minimizadas para ver si existe una correlacion
#entre la energia de la molecula luego de la minimización y antes: de esta forma veo si es mejor
#minimizar antes de hacer docking en este sistema o es indistinto.


# In[ ]:


def OBMolMinimize(molecule, formato_input, formato_out,steps = 2500,):
        '''Minimiza una molecula con Steepest. Descendent y el campo de fuerza mmff94
        default: steps = 2500, e = 1e-6f'''
        counter = 0
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats(formato_input, formato_out)
        mol = openbabel.OBMol()
        forcefield = openbabel.OBForceField.FindForceField("MMFF94")
        notatend = obConversion.ReadFile(mol, molecule)
        while notatend:
                counter+= 1
                name = molecule.strip().replace('.sdf', '')+'_mmff94.sdf'
                forcefield.Setup(mol)
                print(f'Minimizando {name}')
                forcefield.SteepestDescent(steps)
                forcefield.GetCoordinates(mol)
                obConversion.WriteFile(mol,name)
                mol = openbabel.OBMol()
                notatend = obConversion.Read(mol)


# In[323]:


os.system('ls *.sdf > decoys_activos.txt')


# In[324]:


#Minimizacion de ligandos con openbabel utilizando el campo de fuerza MMFF94
with open('decoys_activos.txt', 'r') as decoys_a_minimizar:
    for molecula in decoys_a_minimizar:
        OBMolMinimize(molecula.strip(), 'sdf', 'sdf')        


# In[ ]:


#Calculo de energias para las estructuras sin minimizar vs minimizadas


# In[ ]:


# Correlacion entre energias para ver si vale la pena mimimizar


# In[ ]:


# Docking


# In[ ]:


# Construcción de una curva ROC

