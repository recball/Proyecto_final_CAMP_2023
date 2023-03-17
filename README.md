
Objetivo del trabajo

Efectuar la búsqueda asistida por computadora de de inhibidores de la proteina CDK2, utilizando herramientas como el dcoking y la dinamica molecular, con un postprocesamiento utilizando herramientas de python.
La descripcion de lo hecho hasta el momento corresponde a la parte inicial del proyecto que consiste en la validacion del programa de docking a utilizar
basado en la performance para disriminar entre compuestos activos y decoys, lo que se cuantifica a traves de la contruccion de una curva ROC y el 
calculo del area bajo la misma. 


PAQUETES UTILIZADOS
-OPENBABEL : https://openbabel.org/docs/dev/Installation/install.html
-PANDAS: pip install pandas 
-MATPLOTLIB: pip install matplotlib
-OS - parte de la libreria estandar de python
-SYS - parte de la libreria estandar de python
-NUMPY - pip install numpy
-RANDOM - parte de la libreria estandar de python
-SUBPROCESS - parte de la libreria estandar de python
-AUTODOCKTOOLS- descargar desde https://autodock.scripps.edu/download-autodock4/
-smina.statics- docking() - https://sourceforge.net/projects/smina/



CARPETAS

ligandos_de_partida: ESTRUCTURAS DE LOS COMPUESTOS ACTIVOS Y DECOYS GENERADAS CON OPENBABEL
ligandos_minimizados: ESTRUCTURAS DE LOS DECOYS EN FORMATO PDB, PDBQT MINIMIZADOS CON EL CAMPO DE FUERZA MMFF94
ligandos_no_minimizados: ESTRUCTURAS DE LOS DECOYS Y LIGANDOS SIN MINIMIZAR
receptor: ESTRUCTURA DE LA PROTEINA CDK2
GRAFICOS: RELACION ENTRE NUMERO DE COMPUESTOS ACTIVOS Y VALOR DE AUC Y ENRICHMENT FACTOR
scripts_utilizados: SCRIPTS UTILIZADOS EN EL ANALISIS 
construccionCurvaRoc: CURVAS ROC Y DATASETS EVALUADOS
subconjunto_activos_decoys:ARCHIVOS .PICKED SELECCIONADOS ALEATORIAMENTE PARA CONSTRUIR LA CURVA ROC A PARTIR DE 50 LIGANDOS Y 2500 DECOYS


CORRIDA
 
1) Obtencion del dataset de moleculas para procesar (activos + decoys) para el armado de la curva ROC

A) 'procesamiento_conformeros.ipynb' localizado dentro de la carpeta proyecto_final_camp
Las funciones de este script son las de descargar los decoys, generar las conformaciones 3D y minimizar las estructuras generadas

celda #1: DESCARGA DE ARCHIVOS .PICKED DESDE LA BASE DE DATOS DE DECOYS DUDE: https://dude.docking.org/ PARA EL TARGET CDK2 (en el directorio raiz).
La salida de esta celda es un archivo de tipo html con la información del total de decoys disponibles en la base de datos.

celda #2: Extraccion de los archivos .picked donde se almacenan los smiles de decoys y compuestos activos desde el archivo html

celda #3:Obtencion de las direcciones de los decoys en DUDE

celda #4: El numero de moleculas a generar y luego largar el docking es demasiado grande, por lo que para esta prueba piloto seleccionaré un subconjunto de ligandos: 50 ligandos activos y 2500 decoys, seleccionados de manera aleatoria.

celda #5: Abrir los archivos .picked para extraer ligandos, y obtener una lista de los compuestos activos para la quinasa que se encuentran en la primera linea.

celda #6: Extraccion de los smiles de los decoys desde cada archivo .picked

celda #7: Generar las conformaciones desde los smiles utilizando openbabel como generador de estructuras 3D

celda #8: Función para minimizar ligandos con el campo de fuerza MMFF94.

celda #10: Minimizacion de ligandos con openbabel utilizando el campo de fuerza MMFF94.

celda #11: Mover los archivos .sdf minimizados y no minimizados a sus respectivas carpetas.


2) Selección del receptor a utilizar: Estructura de CDK2 disponible en la base de datos protein data bank (PDB)

https://www.rcsb.org/structure/4kd1
Resolution: 1.70 Å  R-Value Free: 0.232 
ligando: 	3-[({3-ethyl-5-[(2S)-2-(2-hydroxyethyl)piperidin-1-yl]pyrazolo[1,5-a]pyrimidin-7-yl}amino)methyl]-1-hydroxypyridinium
Cadena A seleccionada para estudios posteriores

3) Seleccion-estudio del sitio activo, tomando todos los residuos a 5A del ligando cristalizado 
protonación con el servidor molprobity, pH 7.4:_http://molprobity.biochem.duke.edu/

4) Dentro de la carpeta de ligandos minimizado ejecutar el script SDFtoPDB.py convertir los SDF a PDB ya que es el formato reconocido por prepare_ligand4.py, este recibe como argumento una molecula en formato sdf y devuelve su equivalente en PDB

bash: for ligando in $(ls *.sdf); do python3 SDFtoPDB.py $ligando; done

5) Obtencion de los archivos de ligandos y decoys en formato .pdbqt utilizando el script en python de ambertools (prepare_ligand4.py)

bash: for ligando in $(ls *.pdb); do prepare_ligand4.py -l $ligando; done

6) Con prepare_receptor.py -r receptor

7) Preparación de los archivos del dataset con smina (preparación de archivos de entrada con el script en bash de otra persona, que prefiere que no lo suba a git)

8) Correr el script .sh localizado dentro de docking: para cada conjunto: activos y decoys. La salida de

9) postprocesamiento de resultados de docking:
COn el fin de avaluar como varia la ROC con el tamaño del dataset (ligandos + activos) con el numero de ligandos activos y decoys,
en dude por cada activo se generan 50 decoys, por lo que se evaluaron conjuntos de 5,10,20,30 y 50 ligandos activos y los valores de AUC obtenidos 
para cada caso. 
Dentro de las carpetas docking_decoys, docking_activos  existen carpetas que reflejan el numero de decoys para cada caso, y de activos.
Las salidas del docking son las conformaciones de cada ligando en formato pdbqt y un archivo .log que contiene los valores de score para las diferentes conformaciones obtenidas

10) para cada carpeta dentro de las carpetas Con el script smina_score.py extraer los scores obtenido para cada molecula por docking a partir de un archivo que contenga los nombres de los archivos .log y el nombre del archivo final con estos scores. 


docking_decoys
250
500
1000
1500
4kd1FH (2500 decoys)

docking_activos
5
10
20
30
50
4kd1FH (50 activos)

bash: smina_score.py lista_log nombre_del_archivo_final

11) Mover los archivos con los scores a la carpeta construccionCurvaRoc y definir subcarpetas con el siguiente formato
numeroDEactivos_activos_numeroDEdecoys_decoys.

12) Dentro de cada carpeta se encuentra un archivo ROC que contiene scores para activos y decoys distinguir entre ellos marcando los activos con el nombre 'activo' antes del final del nombre,
correr el script rocker.py (externo para calcular el area de la curva): python3 rocker.py -an 'activos' -EF -r -c 2 ROC.csv dentro de cada una para obtener ROC, AUC_value EF y el grafico de la ROC para cada uno (salidas)

an activos: marcaje que distingue activos
EF: calculo del factor de enriquecimiento
r: score reverso (scores mas negativos como mejores)
c 2: utilizar columna 2 con scores
ROC: archivo con datos

Guardar el grafico y desde consola los valores de AUC, y EF en la carpeta graficos 

13) Utilizar el script relacion_parametros_ROC.ipynb para graficar la AUC (en cada caso) vs el numero de compuestos activos (AUC_vs_activos.txt)

NUMERO,AUC,ENRINCHMENT_FACTOR
5,0.851199999999999,17.0 
10,0.7444000000000031, 10.200000000000001
20,0.7697905671907174,5.1000000000000005
30,0.7844222222222307,10.200000000000001
50,0.7697905671907174,8.677150786308973
(en la carpeta se encuentra un archivo AUC_vs_activos que contiene los valores de AUC para cada conjunto de activos: 5,10,15,20,50)

14) En la visualizacion del grafico queria observar si existia una relacion lineal entre aumentar el numero de ligandos activos y el incremento del area bajo
la curva ROC

CONCLUSIONES Y OBSERVACIONES
- lA finalidad de este proyecto (que enmarca dentro del diseño de inhibidores de CDK2) es evaluar la eficiencia de smina para clasificar entre compuestos activos
e inactivos hacia esta quinasa, esto se realiza sobre conjuntos de diferente tamaño para ver si el area de la curva ROC (Area Under Curve) sigue un comportamiento constante en valor, aumenta o disminuye al aumentar o disminuir el numero de activos y decoys, en una proporcion de 1 activo cada 50 decoys.

- A grandes rasgos, siguiendo lo explicado arriba para el codigo y salidas el primer paso consiste en obtener los decoys y activos desde DUDE en fornato smiles
reconstruirlos a 3D utilizando openbabel, minimizar las estructuras y prepararlas en formato pdbqt para corridas de docking. Por otro lado se hace la seleccion del receptor a utilizar de CDK2 desde protein data Bank en base a la resolucion de la estructura y los parametros de calidad (menos de 2 anstrongs).
El docking fue utilizado para obtener las conformaciones para cada molecula dentro del sitio activo de esta quinasa, y los scores de cada una en un archivo .log, y estos archivos fueron posteriormente procesados haciendo una seleccion aleatoria de compuestos activos y decoys con las relaciones (5 activos:250 decoys, 10 activos: 500 decoys, 20 activos: 1000 decoys, 30 activos: 1500 decoys, 50 activos: 2500 decoys), ya para cada caso se graficaron las curvas ROC correspondientes para obtener el area de cada curva. El analisis exploratorio de los datos mostró que a menor numero de activos la curva adquiere un mayor valor, lo cual tiene sentido debido a que existe menos diversidad de activos y si clasifica bien por ejemplo 5 (en la de 5 activos), el area es mayor, mientras que al aumentar el numero de activos esta area se reduce debido a que el clasificador esta puesto a prueba con un dataset mayor y comete mas errores en la clasificación, por lo que lo mejor para evaluar la metodologia de docking es utilizar un dataset con una alta diversidad de compuestos y alto en numero. Entre los errores de mi estudio se encuentran justamente la selección aleatoria de cada decoy y activo, debido a que no se sabe si se eligen compuestos muy parecidos entre si, o hay una diversidad estructural alta entre ellos. En los proximos analisis que deberian revisarse hace falta revisar los criterios de seleccion para activos y decoyS.
Algo que no mencioné es que un AUC = 1 habla de una clasificacion perfecta, un AUC = 0.5 o menor habla de una seleccion al azar entre activos y decoys y una mala performance del programa de docking.

 

