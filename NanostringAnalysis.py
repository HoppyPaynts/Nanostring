#!usr/bin/python

##############################################################################
# Benjamin Rambo-Martin
# October, 2014
#
#   This script produces Copy Number counts from Nanostring samples. Run from 
# directory containing .RCC nanostring raw data files. Also requires the
# SampleManifest from Nanostring. Create a text file of sample manifest
# containing only ID rows. Reference is set as '20140922_cnv11_em11_12.RCC'. 
# Two output files are produced, one with copy number counts for each sample 
# for each locus, and another with all data and calculations.
##############################################################################



import numpy as np
from scipy import stats
import glob
import re

# Build sample key dictionary
SampleKey = {}
keyfile = open('SampleManifest.txt', 'r')
for line in keyfile:
    if 'cnv' in line:
        array = line.split()[9][:-3]
        em = 'em' + str(line.split()[9].split('_')[1][3:]) + '_'
        num = str(line.split()[11])
        array += em + str(0) + num + '.RCC'
        if '20140922_cnv11_em11_12.RCC' in array:                       # SET REFERENCE HERE
            array = 'Reference_' + array
        SampleKey[array] = line.split()[1]
        arrayroot = array[:-6]
    else:
        if int(line.split()[9]) < 10:
            array = arrayroot + str(0) + str(line.split()[9]) + '.RCC'
            if '20140922_cnv11_em11_12.RCC' in array:                   # SET REFERENCE HERE
                array = 'Reference_' + array
            SampleKey[array] = line.split()[1]
        else:
            array = arrayroot + str(line.split()[9]) + '.RCC'
            if '20140922_cnv11_em11_12.RCC' in array:                   # SET REFERENCE HERE
                array = 'Reference_' + array
            SampleKey[array] = line.split()[1]
keyfile.close()

# Build Main Dictionary with filenames as Keys and empty dictionary as Values
files = glob.glob('*RCC')
MainD = {}
for i in files:
    MainD[i] = {}

# Reference file check
if '20140922_cnv11_em11_12.RCC' not in files:
    print '\n\nReference file \'20140922_cnv11_em11_12.RCC\' required!\n\nOPERATION TERMINATED\n\n'
    import sys
    sys.exit(0)

# Existing data Dictionary builder function
def DictMaker(filehandle, line, MainD=MainD):
    if 'Positive' in line:
        MainD[filehandle]['Positive_conc'].append(re.split('\(|\)',line)[1])
        MainD[filehandle]['Positive_count'].append(line.split(',')[3].rstrip())
    if 'RESTRICTIONSITE+B' in line:
        MainD[filehandle]['b'] = float(line.split(',')[3].rstrip())
    if 'RESTRICTIONSITE-C' in line:
        MainD[filehandle]['c'] = float(line.split(',')[3].rstrip())
    if 'RESTRICTIONSITE+A' in line:
        MainD[filehandle]['a'] = float(line.split(',')[3].rstrip())
    if 'RESTRICTIONSITE-D' in line:
        MainD[filehandle]['d'] = float(line.split(',')[3].rstrip())
    if 'Invariant' in line:
        MainD[filehandle]['Invariant'].append(line.split(',')[3].rstrip())
        MainD[filehandle]['Invariant_probe_names'].append(line.split(',')[1])        
    if 'Endogenous' in line:
        MainD[filehandle]['Endogenous'].append(line.split(',')[3].rstrip())
        MainD[filehandle]['Endogenous_probe_names'].append(line.split(',')[1])

# Fill individual sample dictionaries with available data
for k in MainD.iterkeys():
    MainD[k]['Positive_conc'] = []
    MainD[k]['Positive_count'] = []
    MainD[k]['Invariant'] = []
    MainD[k]['Invariant_probe_names'] = []
    MainD[k]['Invariant_norm'] = []
    MainD[k]['Invariant_CN'] = []
    MainD[k]['Endogenous'] = []
    MainD[k]['Endogenous_probe_names'] = []
    MainD[k]['Endogenous_norm'] = []
    MainD[k]['Endogenous_CN'] = []
    f = open(k, 'r')
    for line in f:
        DictMaker(k, line)
    f.close()
    posconc = np.log10(map(float, MainD[k]['Positive_conc']))
    poscount = np.log10(map(float, MainD[k]['Positive_count']))    
    MainD[k]['r2'] = np.square(stats.linregress(posconc, poscount)[2])
    MainD[k]['Restriction_control'] = [int(MainD[k]['c']/MainD[k]['b']), int(MainD[k]['d']/MainD[k]['a'])]
    MainD[k]['Invariant_mean'] = np.mean(map(float,MainD[k]['Invariant']))

# Calculate invariant probe mean's mean
invariant_means = []
for k in MainD.iterkeys():
    invariant_means.append(MainD[k]['Invariant_mean'])
invarmean = np.mean(invariant_means)

# Calculate normalization factor
for k in MainD.iterkeys():
    MainD[k]['Normalization_factor'] = invarmean / MainD[k]['Invariant_mean']

# Calculate normalized probe counts
for k in MainD.iterkeys():
    MainD[k]['Endogenous'] = map(float, MainD[k]['Endogenous'])
    MainD[k]['Invariant'] = map(float, MainD[k]['Invariant'])
    for j in MainD[k]['Endogenous']:
        MainD[k]['Endogenous_norm'].append(j*MainD[k]['Normalization_factor'])
    for l in MainD[k]['Invariant']:
        MainD[k]['Invariant_norm'].append(l*MainD[k]['Normalization_factor'])

# Set reference sample
MainD['Reference_20140922_cnv11_em11_12.RCC'] = MainD.pop('20140922_cnv11_em11_12.RCC') # SET REFERENCE HERE

# Calculate Copy Number Values
for k in MainD.iterkeys():
    for j in range(0,len(MainD[k]['Endogenous_norm'])):
        MainD[k]['Endogenous_CN'].append(MainD[k]['Endogenous_norm'][j] / MainD['Reference_20140922_cnv11_em11_12.RCC']['Endogenous_norm'][j])
    for l in range(0,len(MainD[k]['Invariant_norm'])):
        MainD[k]['Invariant_CN'].append(MainD[k]['Invariant_norm'][l] / MainD['Reference_20140922_cnv11_em11_12.RCC']['Invariant_norm'][l])

# Print output
invariant_probe_names = MainD['Reference_20140922_cnv11_em11_12.RCC']['Invariant_probe_names']      # SET REFERENCE HERE
endogenous_probe_names = MainD['Reference_20140922_cnv11_em11_12.RCC']['Endogenous_probe_names']    # SET REFERENCE HERE
CNout = open('NanostringCopyNumberCounts.txt', 'a+')
print>>CNout, 'Probe',
for k in MainD.iterkeys(): # Header
    print>>CNout, SampleKey[k],
print>>CNout
for i in range(0,len(invariant_probe_names)): # Invariant CNs
    print>>CNout, invariant_probe_names[i].split('|')[3],
    for k in MainD.iterkeys():
        print>>CNout, '%.2f' % MainD[k]['Invariant_CN'][i],
    print>>CNout
for i in range(0, len(endogenous_probe_names)): # Endogenous CNs
    print>>CNout, endogenous_probe_names[i].split('|')[3],
    for k in MainD.iterkeys():
        print>>CNout, '%.2f' % MainD[k]['Endogenous_CN'][i],
    print>>CNout
CNout.close()

# Print Integer CN output
def CN(i):
	i = i*3
	if i < 0.4:
		i = 0
	elif 0.6 < i < 1.4:
		i = 1
	elif 1.6 < i < 2.4:
		i = 2
	elif 2.6 < i < 3.4:
		i = 3
	elif 3.6 < i < 4.4:
		i = 4
	elif 4.6 < i < 5.4:
		i = 5
	elif 5.6 < i < 6.4:
		i = 6
	return i

CNint = open('NanostringCopNumberCounts.integer.txt', 'a+')
print>>CNint, 'Probe',
for k in MainD.iterkeys(): # Header
    print>>CNint, SampleKey[k],
print>>CNint
for i in range(0,len(invariant_probe_names)): # Invariant CNs
    print>>CNint, invariant_probe_names[i].split('|')[3],
    for k in MainD.iterkeys():
    	c = CN(MainD[k]['Invariant_CN'][i])
        print>>CNint, '%.2f' % c,
    print>>CNint
for i in range(0, len(endogenous_probe_names)): # Endogenous CNs
    print>>CNint, endogenous_probe_names[i].split('|')[3],
    for k in MainD.iterkeys():
    	c = CN(MainD[k]['Endogenous_CN'][i])
        print>>CNint, '%.2f' % c,
    print>>CNint
CNint.close()

#Print full output
CompleteOut = open('Nanostring_all_calculations.txt', 'a+')
print>>CompleteOut, 'Feature',
for k in MainD.iterkeys():
    print>>CompleteOut, k,
print>>CompleteOut
print>>CompleteOut, 'SampleID',
for k in MainD.iterkeys():
    print>>CompleteOut, SampleKey[k],
print>>CompleteOut
for i in range(0,6):
    print>>CompleteOut, 'Positive_concentration_' + str(i), #positive concentration
    for k in MainD.iterkeys():
        print>>CompleteOut, MainD[k]['Positive_conc'][i],
    print>>CompleteOut
for i in range(0,6):
    print>>CompleteOut, 'Positive_count_' + str(i), #positive count
    for k in MainD.iterkeys():
        print>>CompleteOut, MainD[k]['Positive_count'][i],
    print>>CompleteOut
print>>CompleteOut, 'Rsquared', #r2
for k in MainD.iterkeys():
    print>>CompleteOut, '%.2f' % MainD[k]['r2'],
print>>CompleteOut
for i in range(0,2): #restriction control
    print>>CompleteOut, 'RestrictionQC_' + str(i),
    for k in MainD.iterkeys():
        print>>CompleteOut, MainD[k]['Restriction_control'][i],
    print>>CompleteOut
for i in range(0, len(invariant_probe_names)): #invariant raw count
    print>>CompleteOut, 'RawCount_' + invariant_probe_names[i],
    for k in MainD.iterkeys():
        print>>CompleteOut, MainD[k]['Invariant'][i],
    print>>CompleteOut
print>>CompleteOut, 'Invariant_mean', #invariant mean
for k in MainD.iterkeys():
    print>>CompleteOut, MainD[k]['Invariant_mean'],
print>>CompleteOut
print>>CompleteOut, 'Mean_of_Invariant_Means', #mean of invariant means
for i in range(0,len(MainD)):
    print>>CompleteOut, invarmean,
print>>CompleteOut
print>>CompleteOut, 'Normalization_factor', #normalization factor
for k in MainD.iterkeys():
    print>>CompleteOut, MainD[k]['Normalization_factor'],
print>>CompleteOut
for i in range(0, len(invariant_probe_names)): #invariant normalized count
    print>>CompleteOut, 'NormalizedCount_' + invariant_probe_names[i],
    for k in MainD.iterkeys():
        print>>CompleteOut, MainD[k]['Invariant_norm'][i],
    print>>CompleteOut
for i in range(0, len(invariant_probe_names)): #invariant CN
    print>>CompleteOut, 'CN_' + invariant_probe_names[i],
    for k in MainD.iterkeys():
        print>>CompleteOut, MainD[k]['Invariant_CN'][i],
    print>>CompleteOut
for i in range(0, len(endogenous_probe_names)): #endogenous raw count
    print>>CompleteOut, 'RawCount_' + endogenous_probe_names[i],
    for k in MainD.iterkeys():
        print>>CompleteOut, MainD[k]['Endogenous'][i],
    print>>CompleteOut
for i in range(0, len(endogenous_probe_names)): #endogenous normalized count
    print>>CompleteOut, 'NormalizedCount_' + endogenous_probe_names[i],
    for k in MainD.iterkeys():
        print>>CompleteOut, MainD[k]['Endogenous_norm'][i],
    print>>CompleteOut
for i in range(0, len(endogenous_probe_names)): #endogenous CN
    print>>CompleteOut, 'CN_' + endogenous_probe_names[i],
    for k in MainD.iterkeys():
        print>>CompleteOut, MainD[k]['Endogenous_CN'][i],
    print>>CompleteOut

CompleteOut.close()
