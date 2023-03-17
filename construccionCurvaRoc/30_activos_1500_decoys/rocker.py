#!/usr/bin/env python
#-*- coding: utf-8 -*-

'''
Script for plotting ROC-curves (Receiver Operator Characteristic) and 
calculation of AUC (Area Under Curve) and two kinds of enrichment factors
'''

from __future__ import division
import sys
import os
import re
from operator import itemgetter
import matplotlib as mpl
from matplotlib import pyplot
import argparse
import numpy as np
import matplotlib.font_manager as font_manager

__version__ = '0.1.4'

class CenteredFormatter(mpl.ticker.ScalarFormatter):
    """Acts exactly like the default Scalar Formatter, but yields an empty
    label for ticks at "center"."""
    center = 0
    def __call__(self, value, pos=None):
        if value == self.center:
            return ''
        else:
            return mpl.ticker.ScalarFormatter.__call__(self, value, pos)

def plotter(**kwargs):
    '''Runtime plotter function.'''
    
    if kwargs.get('font', None):
        mpl.rc('font', family=str(kwargs.get('font')))
    
    if kwargs.get('axeswidth', None):
        for thing in ('axes.linewidth', 'xtick.major.pad', 'xtick.major.size', 'xtick.major.width', 'xtick.minor.pad', 'xtick.minor.size', 'xtick.minor.width', 'ytick.major.pad', 'ytick.major.size', 'ytick.major.width', 'ytick.minor.pad', 'ytick.minor.size', 'ytick.minor.width'):
            mpl.rcParams[thing] = mpl.rcParams[thing]*kwargs.get('axeswidth')
    
    sortreverse = not kwargs.get('score_reverse',False)
    
    ofname = kwargs.get('pic',None)
    if ofname:
        filename, ftype = os.path.splitext(ofname)
        if not imagetest(ftype):
            print( "\"{}\"-format is not supported!".format(ftype) )
            sys.exit(0)
        
    actives = []
    label = []
    scores = []
    
    #Handling of special kwargs.
    specidic = matchbox( kwargs.get('kwargs','') )
    
    for i in range(0,len(kwargs['input'])):
        if kwargs.get('activelist', False):
            afname = kwargs['activelist'][i] #sys.argv[i]
            sfname = kwargs['input'][i] #sys.argv[i+1]
            actives.append( readActives(afname) )
            print("Loaded {} actives from {}".format(len(actives[-1]), afname) )
            
        elif kwargs.get('activename', False):
            sfname = kwargs['input'][i]
            actives.append( readActivesByRegex(sfname, kwargs.get('activename') ) )
            print("Loaded {} actives from {} starting with {}".format(len(actives[-1]), sfname, kwargs.get('activename')))
        
        columnN = kwargs['column'][i] if len(kwargs['column']) > 1 else kwargs['column'][0]
        
        lab, sco = readScores(sfname,columnN)
        label.append(lab)
        scores.append(sorted(sco, key = lambda item: item[1],reverse=sortreverse))
        print("Loaded {} {} scores from {}".format(len(scores[-1]), label[-1], sfname))
        
    
    #ROC / AUC
    if True:    
        if kwargs.get('legenditems',False):
            label = kwargs.get('legenditems')
        #label = kwargs.get('legenditems', label)
        
        colors = kwargs.get('colors')
        linest = kwargs.get('styles')
        fg = kwargs.get('figsize')
        dperi = kwargs.get('dotsperinch')
        legend = kwargs.get('legend')
        
        
        linew = kwargs.get('linewidth', 2)
        if not kwargs.get('noROC',False):
            fig = mpl.pyplot.figure(figsize=fg, dpi=dperi, **specidic.get('figure',{}))
            ax = fig.add_subplot(1,1,1)
            
        setdict = {
        'randomline'    :   not kwargs.get('norandom', False), 
        'randomleg'     :   kwargs.get('legenditems')[-1] if kwargs.get('legenditems') else 'random',
        'loglimit'      :   kwargs.get('logplot', 0.0), 
        'XYlab'         :   kwargs.get('labels'), 
        'labelsize'     :   kwargs.get('labelsize', False), 
        'ticksize'      :   kwargs.get('ticksize', False),
        'annotations'   :   kwargs.get('annotations', tuple()),
        'annotationsize':   kwargs.get('annotationsize', False),
        'linewidth'     :   linew,
        #'fontname'      :   kwargs.get('font', None),
        }
        
        if not kwargs.get('noROC',False):
            initializePlot(**setdict)
            
            logval = kwargs.get('logplot', False)
            logval = logval if logval!=0 else False
        
        for i, item in enumerate(actives):
            tpr, fpr = calcRates(actives[i], scores[i])
            
            if not kwargs.get('noROC',False):
                print("Plotting ROC Curve ...")
            
                if logval:
                    fpr=list( map(lambda x: x if x>0 else logval, fpr) )
                
                mpl.pyplot.plot(fpr, tpr, color=colors[i%len(colors)], linestyle=linest[i%len(linest)], linewidth=linew, label= str(label[i%len(label)]), **specidic.get('plot', {}))
                
        if not kwargs.get('noROC',False):
            randomPlot(**setdict)
        
            if logval:
                ax.set_xscale('log')
                ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x,pos:str(x)))
            elif kwargs.get('noorigin',False):
                formatter = CenteredFormatter()
                formatter.center = 0
                mpl.pyplot.gca().xaxis.set_major_formatter(formatter)
                mpl.pyplot.gca().yaxis.set_major_formatter(formatter)
            legkwar = specidic.get('legend',{})
            savePlot(kwargs.get('logplot', False), legend, kwargs.get('legendsize',False), ofname, dperi, specidic.get('legend',{}) )
        
    if kwargs.get('enrichmentfactor', False):
        #per = kwargs.get('enrichmentfactor') if isinstance(kwargs.get('enrichmentfactor'), (float, int)) else 1.0
        per = kwargs.get('enrichmentfactor')
        for i in range(len(actives)):
            calcEF(actives[i],scores[i], per)
        
    if kwargs.get('enrichmentfactordec', False):
        #perd = kwargs.get('enrichmentfactordec') if isinstance(kwargs.get('enrichmentfactordec'), (float, int)) else 1.0
        perd = kwargs.get('enrichmentfactordec')
        for i in range(len(actives)):
            calcEFdec(actives[i],scores[i], perd)
    
    if kwargs.get('BEDROC', False):
        #perd = kwargs.get('enrichmentfactordec') if isinstance(kwargs.get('enrichmentfactordec'), (float, int)) else 1.0
        alpha = kwargs.get('BEDROC')
        for i in range(len(actives)):
            calcBEDROC(actives[i],scores[i], alpha)
            
    

def calcEF(actives, scores, percent):
    '''
    Calculate Enrichment factor.
    EF = (act_pro/mol_pro)/(act_all/mol_all)
    '''
    length = int(np.round(len(scores)*percent/100.0))
    sub = scores[:length]
    suba = [1 if a[0] in actives else 0 for a in sub]
    EF = (float(sum(suba))/len(sub))/(float(len(actives))/len(scores))
    
    print('EF_{} = {}'.format(percent,EF))
    
def calcEFdec(actives, scores, percent):
    '''
    Calculate Enrichment factor.
    n = mols when n in_acts found
    EF = act_n / act_all
    '''
    length = int(np.round((len(scores)-len(actives))*percent/100.0))
    newsco = 0
    decs = 0
    for item in scores:
        if decs >= length:
            break
        elif item[0] in actives:
            newsco += 1
        else:
            decs += 1
    
    EF = 100.0*newsco/len(actives)
    print('EF_{}%decs = {}'.format(percent,EF))
    
def calcBEDROC(actives, scores, alpha):
    '''
    Calculate Boltzmann-Enhanced Discrimination of ROC
    '''
    
    N = len(scores)
    Ra = len(actives)/N
    #Ri = 1 - Ra
    
    RIEmin = (1 - np.exp(alpha*Ra))/(Ra*(1-np.exp(alpha)))
    RIEmax = (1 - np.exp(-alpha*Ra))/(Ra*(1-np.exp(-alpha)))
    
    acts = [(i+1)/N for i, (name, score) in enumerate(scores) if name in actives]
    
    RIE = (1/Ra)*(np.sum([np.exp(-alpha*y) for y in acts]))*(np.exp(alpha/len(scores))-1)/(1-np.exp(-alpha))
    BEDROC = (RIE - RIEmin)/(RIEmax - RIEmin) #(float(sum(suba))/len(sub))/(float(len(actives))/len(scores))
    
    
    print('BEDROC_{} = {}'.format(alpha,BEDROC))
    
def printfonts():
    for font in listfonts():
        print(font)
    

def readActives(fname):
    '''
    Read active names from file
    '''
    actives = []
    with open(fname, 'r') as f:
        for line in f:
            id = line.strip()
            actives.append(id)
    return actives
    
def readActivesByRegex(fname, regexes):
    '''
    Read active names from scores file by regex
    '''
    myregs = [re.compile(regex) for regex in regexes]
    actives = []
    with open(fname, 'r') as f:
        for line in f:
            iden = re.split('[\s,;]', line)[0].strip() #line.strip()
            for myreg in myregs:
                if myreg.match(iden):
                    actives.append(iden)
                    break
    return actives

def readScores(path,column):
    '''
    Read scores file and return names and scores
    '''
    with open(path, 'r') as f:
        data = [[cell.strip() for cell in re.split('[\t ]', line)] for line in f]
    label = data[0][column-1].strip()
    
    scores = []
    for row in data[1:]:
        if len(row)==0:
            continue
        scores.append((row[0],float(row[column-1])))
    return label, scores

def calcRates(actives, scores, AUC=True, EF=False, dec=False):
    '''
    Calculate True positive and false positive rates and print enrichment factors.
    '''
    tpf = [0]  # true positives found
    fpf = [0]  # false positives found
	
    foundactives = 0
    founddecoys = 0
    for idx, (id, score) in enumerate(scores):
        if id in actives:
            foundactives += 1
        else:
            founddecoys += 1
        tpf.append(foundactives)
        fpf.append(founddecoys)
    nractives = foundactives
    nrdecoys = founddecoys
    tpr = list( map(lambda x: float(x)/foundactives, tpf) )
    fpr = list( map(lambda x: float(x)/founddecoys, fpf) )
    
    area = 0
    for i in range(len(tpr)-1):
        if fpr[i] == fpr[i+1]:
            continue
        else:
            area += tpr[i]
    
    narea = area/float(nrdecoys)
    
    q1=narea/(2-narea)
    q2=2*narea**2/(1+narea)
    
    ernarea = (((narea-narea**2)+(nractives-1)*(q1-narea**2)+(nrdecoys-1)*(q2-narea**2))/(nrdecoys*nractives))**0.5
    
    print('AUC={}+-{}'.format(narea,ernarea))
    return tpr, fpr

def kwargParse(string):
    params = string.strip().split(',')
    kwargs = {}
    for param in params:
        a, b = param.split('=')
        kwargs[a.strip()] = eval( b.strip() )
    return kwargs

def initializePlot(**kwargs):
    '''
    Initialize a ROC plot
    '''
    
    XYlab= kwargs.get('XYlab', ('FPR','TPR')) #('FPR','TPR') if not kwargs.get('XYlab',False) else kwargs.get('XYlab')
    '''
    fontsize = kwargs.get('labelsize', 14) #14 if not kwargs.get('fontsize',False) else kwargs.get('fontsize')[1]
    ticksize = kwargs.get('ticksize', 14)
    '''
    
    annotations = kwargs.get('annotations', tuple() )
    annotationsize = kwargs.get('annotationsize', 14)
    
    labels = {'fontsize': kwargs.get('labelsize', 14)}
    ticks = {'labelsize': kwargs.get('ticksize', 14)}
    
    anno = {}
    if 'fontname' in kwargs:
        labels['fontname'] = kwargs['fontname']
        anno['fontname'] = kwargs['fontname']

    #pyplot.xlabel( str(XYlab[0]), fontsize=fontsize)
    #pyplot.ylabel( str(XYlab[1]), fontsize=fontsize)
    
    mpl.pyplot.xlabel( str(XYlab[0]), **labels )
    mpl.pyplot.ylabel( str(XYlab[1]), **labels )
    
    #pyplot.tick_params(axis='x', labelsize=ticksize)
    #pyplot.tick_params(axis='y', labelsize=ticksize)    
    mpl.pyplot.tick_params(axis='x', **ticks)
    mpl.pyplot.tick_params(axis='y', **ticks)
    
    for item in annotations:
        splite = [ cell.strip() for cell in item.split(',') ]
        xc,yc,text = (a.strip('\\') for a in splite[:3])
        asize = float(splite[3]) if len(splite)>3 else annotationsize
        mpl.pyplot.text(float(xc),float(yc), str(text), fontsize=asize, **anno)
        
def randomPlot(**kwargs):
    '''
    Initialize a ROC plot
    '''
    randomline = kwargs.get('randomline',False)  #True if not kwargs.get('randomline',False) else kwargs.get('randomline')
    randomleg = kwargs.get('randomleg','random')
    limit = kwargs.get('loglimit', 0) #0 if not kwargs.get('loglimit',False) else kwargs.get('loglimit')

    if randomline:
        #x = [0.0, 1.0]
        x = list(np.linspace(limit,1,101))
        mpl.pyplot.plot(x, x, linestyle='dashed', color='gray', linewidth=kwargs.get('linewidth', 2), label=randomleg)

def savePlot(log=False, legend=False, legsize=12, fname=None, dpi=300, legkwargs={}, savekwargs={}):
    '''
    Save created plot
    '''
    if log:
        mpl.pyplot.xlim(log, 1.0)
    else:
        mpl.pyplot.xlim(0.0, 1.0)
    mpl.pyplot.ylim(0.0, 1.0)
    
    if legend or legend is 0:
        #if isinstance(legend, int):
        try:
            tryleg = int(legend)
            legend = {  0  : 'best'        , #(only implemented for axis legends)
                        1  : 'upper right' ,
                        2  : 'upper left'  ,
                        3  : 'lower left'  ,
                        4  : 'lower right' ,
                        5  : 'right'       ,
                        6  : 'center left' ,
                        7  : 'center right',
                        8  : 'lower center',
                        9  : 'upper center',
                        10 : 'center'      }.get(tryleg, 'best')
        except ValueError:
            pass
        #frameon = legkwargs.get('frameon', False)
        #labelspacing = legkwargs.get('labelspacing', 0.2) 
        mpl.pyplot.legend(fontsize=legsize, loc=legend, **legkwargs) #frameon=frameon, labelspacing=labelspacing, handletextpad=0.1, borderaxespad=0, handlelength=1.5 )
    #if kwargs.get('tight_layout', True):
    if savekwargs.get('tight_layout', True):
        mpl.pyplot.tight_layout()
    if fname:
        mpl.pyplot.savefig(fname, dpi=dpi, **savekwargs)
    else:
        mpl.pyplot.show()

#def plotROC(pyplot, actives, scores, color, label, style='solid'):
def plotROC(actives, scores, color, label, style='solid'):
    '''Plot a single ROC curve'''
    tpr, fpr = calcRates(actives, scores)
    mpl.pyplot.plot(fpr, tpr, color=color, linestyle=style, linewidth=2, label=label)

def imagetest(ftype):
    '''Test if requested image type is supported.'''
    fig = mpl.pyplot.figure()
    return (ftype.strip('. ') in fig.canvas.get_supported_filetypes())

def bobo(x):
    if x in ('True', 'T'):
        return True
    elif x in ('False','F'):
        return False
    else:
        raise ValueError('String not True, T, False or F')
    
def numer(stri):
    for met in (int, float, bobo, str):
        try:
            return met(stri)
        except ValueError:
            continue
        
def matchbox(kwa):
     box = dict()
     if isinstance(kwa, str):
         for mat in re.finditer('([^=:, ]+)\s*[:=]\s*\{(.+?)\}',kwa):
             box[mat.groups()[0]] = dict()
             #print mat.groups()
             for mato in re.finditer('\s*([^=:, ]+)\s*[:=]\s*([^,]+)\s*', mat.groups()[1]):
                 box[mat.groups()[0]][mato.groups()[0]] = numer(mato.groups()[1])
     return box

def listfonts():
    listed = font_manager.get_fontconfig_fonts()
    names = [font_manager.FontProperties(fname=fname).get_name() for fname in listed]
    return sorted(list(set(names)))

class InputError(Exception):
    pass

def main():
    handles(sys.argv[1:])
    
def handles(arguments):
    arger = argparse.ArgumentParser(description= 'Draw ROC curve from molecule score list and defined actives. \n\nVersion {} \n\nIf you publish your work using rocker, please refer to it: \n\nLätti, S., Niinivehmas, S., & Pentikäinen, O. (2016). Rocker : Open source, easy-to-use tool for AUC and enrichment calculations and ROC visualization. Journal of Cheminformatics, 8 (45). doi:10.1186/s13321-016-0158-y'.format(__version__), formatter_class=argparse.RawTextHelpFormatter )
    
    arger.add_argument("input", metavar = 'scores.csv', nargs ='*', type = str, help="Specify a set of score files.")
    
    active = arger.add_mutually_exclusive_group(required = False)
    active.add_argument("-an", "--activename", metavar='name_regex', nargs='+', type=str, default=None, help="Specify regex for name of active ligands. If this argument is given, do not give actives.")
    active.add_argument("-al", "--activelist", metavar = 'actives.csv', nargs ='+', type = str, default = None, help="Specify a set of active files.")
    
    arger.add_argument("-c", "--column", metavar='N', nargs='+', type=int, default=(2,), help="Specify them column that is includes the desired score. Give 1 or 1 for each set of actives/hits. default=2.")
    arger.add_argument("-r", "--score_reverse", action='store_true', default=False, help="Reverses scores so that small (negative) scores are better.")
    
    arger.add_argument("-nro", "--noROC", action='store_true', help="Don't plot ROC (receiver operator charasteristic) curve.")
    arger.add_argument("-p",  "--pic", type=str, metavar='path.png', help="Specify path to output figure.")
    arger.add_argument("-dpi","--dotsperinch", type = int, default=300,    help="Dots per inch. Resolution of output pic.")
    
    arger.add_argument("-s",  "--figsize", type = float, metavar=('Xsize','Ysize'), nargs=2, default=(5,5),  help="Image size in inches. Give 2 numbers.")
    
    arger.add_argument("-la", "--labels", type = str, metavar=('Xlabel','Ylabel'), nargs=2, default=('FPR','TPR'),  help="Labels for X and Y axels.")
    arger.add_argument("-las","--labelsize", type = float, metavar='fontsize', default= 12,  help="Fontsize of x and y axis labels.")
    
    arger.add_argument("-lp", "--logplot", nargs='?', type=float, metavar='X_min', default=0.0, const=0.001, help="Draw logarithmic plot instead of linear.")
    arger.add_argument("-cl", "--colors", type = str, metavar='color', nargs='+', default=('blue','red','black','green','purple','cyan'),  help="Colors of the plots.")
    arger.add_argument("-st", "--styles", type = str, metavar='style', nargs='+', default=('solid',),  help="Linestyle. solid, dashed.")
    arger.add_argument("-lw", "--linewidth", type = float, metavar='width', default = 2,  help="Linewidth.")
    arger.add_argument("-aw", "--axeswidth", type = float, metavar='width',  help="Width of the axes.")
    
    arger.add_argument("-l",  "--legend", metavar='location', nargs='?', const='best', default=False, help="Location of legend. Number from 0 to 10, or one of 'best','upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center' or 'center'. If no location is provided, best is assumed.")
    arger.add_argument("-li", "--legenditems", type = str, metavar='legend', nargs='+', default=False,  help="Texts appearing in legend. Remember to specify random if it's plotted, random should be the last entry.") #not implemented
    arger.add_argument("-les","--legendsize", type = float, metavar='fontsize', default= 12,  help="Fontsize of legend texts.")
    arger.add_argument("-lf", "--listfonts", action = "store_true",  help="List availeable fonts.")
    arger.add_argument("-f",  "--font", type = str, metavar='fontname', default = None,  help="Name of font used.")
    
    arger.add_argument("-nra","--norandom", default=False, action='store_true', help="Don't plot random line.") #not implemented
    
    arger.add_argument("-ts", "--ticksize", type = float, metavar='fontsize', default= 12,  help="Fontsize of tick numbers.")
    
    arger.add_argument("-a",  "--annotations", type = str, metavar="'X,Y,text'", nargs='+', default=tuple(),  help="Make text annotations on your plot. Give a number of strings with dot delimitering X-coordinate,Y-coordinate,text.")
    arger.add_argument("-as", "--annotationsize", type = float, metavar='fontsize', default= 12,  help="Fontsize of annotations.")
    
    arger.add_argument("-EF", "--enrichmentfactor", type=float, metavar='percent', nargs='?', const=1, default=False, help="Calculate enrichment factor. Optionally give percentage. Dafaults to 1.")
    arger.add_argument("-EFd","--enrichmentfactordec", type=float, metavar='percent', nargs='?', const=1, default=False, help="Calculate EFdec, that is fraction of actives found when a given percentage of decoys are found. Optionally give percentage. Dafaults to 1.")
    arger.add_argument("-BR", "--BEDROC", type=float, metavar='alpha', nargs='?', const=20.0, default=False, help="Calculate BEDROC. Alpha defaults to 20.0.")
    
    arger.add_argument("-no", "--noorigin", action="store_true", help="Do not print origo.")
    
    arger.add_argument("-kw",  "--kwargs", type=str, metavar='dict', default={}, help="Give python style kwargs to methods. Only numeric, boolean and string values accepted. Format is 'figure:{dpi:600, rightpadding:1.2},legend:{frameon:False}'. Accepts arguments for figure, plot, legend and savefig.")
    
    args = arger.parse_args(arguments)
    
    if vars(args).get('listfonts', False):
        printfonts()
        sys.exit()
    
    if len(vars(args).get('input')) == 0:
        print('You must give at least one input file if you don\'t call for font names.')
        arger.print_usage()
        sys.exit()
    
        #raise InputError('You must give at least one input file.')
    
    if not ( vars(args).get('activename') or vars(args).get('activelist') ):
        print('You must define actives with --activename or --activelist.')
        arger.print_usage()
        sys.exit()
        #raise InputError('You must define actives with --activename or --activelist.')
    
    
    sys.exit(plotter(**vars(args)))
    
    
if __name__ == "__main__":
    main()
    
