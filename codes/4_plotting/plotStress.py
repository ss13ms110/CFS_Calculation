# script to plot stress maps for an event
# user input required: srcmodID of event

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from shapely.geometry import Polygon, LinearRing, Point
import subprocess as sp
import pandas as pd
import os


# FUNCTIONS =============================================================

def sigFilter(cfs, a, b):
    dummy1 = a*cfs - b
    dummy2 = np.exp(-dummy1)
    sig = 1/(1 + dummy2)
    return sig

def dist3D(lat1, lon1, z1,  lat2, lon2, z2):
    """
    Surface distance (in km) between points given as [lat,lon]
    """
    R0 = 6367.3        
    D = R0 * np.arccos(
        np.sin(np.radians(lat1)) * np.sin(np.radians(lat2)) +
        np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.cos(np.radians(lon1-lon2)))
    r = np.sqrt(np.square(D) + np.square(z1-z2))
    return r
# =========================================================================


# MAIN
# paths =======================================
stressPath = "./workingData/stressData/"
ASPath = "./workingData/interData/ASFiles/"
cataPath = "./workingData/catalog.dat"
figPath = "./figs/stress"
# =============================================

# prams =======================================
sMin = 0.2
sMax = 0.8
a = 10
b = 1
dd = 2.5
sType = ['MAS0', 'MAS', 'OOP', 'VM', 'MS', 'VMS']
# =============================================

# ask user for stress metric =====================================
smResp = raw_input("MAS0, MAS, OOP, VM, MS, VMS? [1/2/3/4/5/6]: ")
if not smResp in ['1', '2', '3', '4', '5', '6']:
    print "Wrong input... [1,2,3,4,5,6]"
    quit()
print

srcmodID = raw_input("Enter SRCMOD ID of slip model: ")
evPkl = stressPath + '/' + sType[int(smResp)-1] + '/' + srcmodID + '.pkl'
if not os.path.exists(evPkl):
    print "Event %s does not exist. Quiting..." %(srcmodID)
    quit()
print

depth = float(raw_input("Enter depth for stress map: "))
if not depth in np.arange(2.5,50,5):
    print "Wrong depth. Quiting..."
    quit()
print
# ================================================================

# get lat lon of main shock
cataRow = sp.check_output("cat %s | grep %s.fsp" %(cataPath, srcmodID), shell=True)
Mlat = float(cataRow.split()[3])
Mlon = float(cataRow.split()[4])

ASFile = ASPath + '/' + srcmodID + '.txt'

# check for AS and poly files
if not os.path.exists(ASFile):
    print "Aftersock file does not found. Quiting..."
    quit()
print

# load stress pickle file
strDf = pd.read_pickle(evPkl)
sCol = 'D'+str(depth)
strDf = strDf.filter(['lat', 'lon', sCol])

stressData = sigFilter(np.array(strDf[sCol]), a, b)

strDf[sCol] = stressData

# filter between min and max
sVals = np.array(strDf[sCol])
sVals[np.where(sVals < sMin)] = sMin + 0.001
sVals[np.where(sVals >=sMax)] = sMax - 0.001
strDf[sCol] = sVals


# get aftershock grids
# load AS data
ASdata = np.loadtxt(ASFile)
ASrows = ASdata[(ASdata[:,3] >= depth-dd) & (ASdata[:,3] <= depth+dd)]
indx = np.zeros(len(strDf))

ASlat = []
ASlon = []
for ASrow in ASrows:
    dist = dist3D(ASrow[0], ASrow[1], depth, strDf['lat'], strDf['lon'], depth)
    
    imin = np.argmin(np.array(dist))
    ASlat.append(strDf.iloc[imin]['lat'])
    ASlon.append(strDf.iloc[imin]['lon'])


# plotting CFS_sigmoid
plt.figure(figsize=(8, 8))

# activate latex text rendering
rc('text', usetex=True)
rc('axes', linewidth=2)
rc('font', weight='bold')

rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

# CS2 = plt.tricontourf(strDf['lon'], strDf['lat'], strDf[sCol], np.linspace(0.2, 0.8, 100), vmin = 0.2, vmax = 0.8, cmap="Reds", linewidth=0.1)
CS2 = plt.tricontourf(strDf['lon'], strDf['lat'], strDf[sCol], np.linspace(0.2, 0.8, 100), vmin = 0.2, vmax = 0.8, cmap="Reds")

m2 = plt.cm.ScalarMappable(cmap="Reds")
m2.set_array(strDf[sCol])
m2.set_clim(0.0,0.8)
cb2=plt.colorbar(m2, ticks=[0.0, 0.4, 0.8], fraction=0.047, pad=0.04)
cb2.set_label(r'$sig(aS_n \thinspace - \thinspace b)$', fontsize=20)
plt.scatter(ASlon, ASlat, marker='s', c='black', s=20)
plt.scatter(Mlon, Mlat, marker='*', c='yellow', s=400)
plt.xlim(min(strDf['lon'])-0.2, max(strDf['lon'])+0.2)
plt.ylim(min(strDf['lat'])-0.2, max(strDf['lat'])+0.2)
plt.tick_params(axis="x", labelsize=18)
plt.tick_params(axis="y", labelsize=18)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("Longitude")
plt.ylabel("Latitute")
# plt.show()
plt.savefig("%s/%s_%s.pdf" %(figPath,  sType[int(smResp)-1], srcmodID), dpi=500)