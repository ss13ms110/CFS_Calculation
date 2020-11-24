# By Shubham Sharma
# Date - 19-11-2019
# Script to calculate ROC curves for different stress metrics

import numpy as np
import os
import funcFile
from funcFile import bcol
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt

resp=raw_input(bcol.BOLD + "Take aftersocks from saved data or from pickle? [d/p]\nIf running first time choose [P]: " + bcol.ENDC)
if not resp in ['D', 'd', 'P', 'p']:
    print bcol.FAIL + "Wrong input... [D/d or P/p]" + bcol.ENDC
    quit()
print
resp2=raw_input(bcol.BOLD + "Do you want to save polygons? [y/n]\nChoose [Y] if running first time: " + bcol.ENDC)
if not resp2 in ['Y', 'y', 'N', 'n']:
    print bcol.FAIL + "Wrong input... [Y/y or N/n]" + bcol.ENDC
    quit()
print
resp3=raw_input(bcol.BOLD + "MAS0, MAS, OOP, VM, MS, VMS? [1/2/3/4/5/6]: " + bcol.ENDC)
if not resp3 in ['1', '2', '3', '4', '5', '6']:
    print bcol.FAIL + "Wrong input... [1,2,3,4,5,6]" + bcol.ENDC
    quit()
print

# paths to directories
pscmp_dir_path = "./workingData/psgrn+pscmp_out/pscmp_out"
slip_dir_path = "./workingData/srcmod"
AS_dir_path = "./workingData/interData/ASFiles"
ROC_dir_path = "./workingData/ROCData/ROCout"
poly_dir_path = "./workingData/interData/polys"
eventFile = "./workingData/events.txt"
pickleFile = "./../rawData/isc_rev.pkl"
stress_dir_path = "./workingData/stressData"

# declare parameters
dist = 100
days = 365
depth_range = np.arange(2.5,50,5)
Mcut = -1
f = 0.4
skempton = 0
pscmpOGF = 'GF'
mas_oop = 8
sType = ['MAS0', 'MAS', 'OOP', 'VM', 'MS', 'VMS']
# Master fault or OOP CFS
if resp3 == '1':
    pscmpOGF = 'okada'
if resp3 == '3':
    mas_oop = 9

# start by listing SRCMOD files
slip_list = os.listdir(slip_dir_path)

for slip_file in slip_list:
    print bcol.OKBLUE + 'Working on ' + bcol.ENDC + slip_file
    slip_file_path = '%s/%s' %(slip_dir_path, slip_file)

    lat_hypo, lon_hypo, z_hypo, yy, mm, dd, M, Nfault, strike, dip, rake, lat, lon, z = funcFile.read_fsp_file(slip_file_path)

    # calculate lat lon of fault corners
    la_corners, lo_corners = funcFile.getCorners(Nfault, lat, lon, z)

    # get buffer region around the fault
    long2km = funcFile.dist(lat_hypo, lon_hypo-0.5, lat_hypo, lon_hypo+0.5)
    distDeg = dist/long2km

    [xBuffer, yBuffer], polyBuffer = funcFile.createBuffer(la_corners, lo_corners, distDeg)

    # paths for pscmp_out files and AS data
    psc_dir_name = slip_file.split('.')[0]
    psc_dir = '%s/%s/%s' %(pscmp_dir_path, pscmpOGF, psc_dir_name)
    as_file_name = "%s/%s.txt" %(AS_dir_path, psc_dir_name)

    # if opted, save polygon
    if resp2 == "y" or resp2 == "Y":
        polyFname = '%s/%s.poly' %(poly_dir_path, psc_dir_name)
        polyid = open(polyFname,'w')
        for i in range(len(xBuffer)):
            polyid.write('%8.3f  %8.3f\n' %(xBuffer[i], yBuffer[i]))


    # load stress tensors, CFS and OOP stress from pscmp_out files ===========================
    latAll = []
    lonAll = []
    depAll = []
    CFSall = []
    sTensor = []
    for depth in depth_range:
        psc_file_path = '%s/coseis_pscmp_%s_%s_km.inp' %(psc_dir, psc_dir_name, depth)
        
        data = np.loadtxt(psc_file_path, skiprows=1, usecols=(0, 1, 5, 6, 7, 8, 9, 10, 17, 21))
        if depth == 2.5:
            areaIndex = funcFile.getRegion(data[:,0:2], polyBuffer)
        
        latAll = np.append(latAll, data[areaIndex,0])
        lonAll = np.append(lonAll, data[areaIndex,1])
        depAll = np.append(depAll, np.full(len(data[areaIndex,0]),depth))
        if resp3 in ['1','2','3']:
            CFSall = np.append(CFSall, data[areaIndex,mas_oop])
        else:
            sTensor = np.append(sTensor, data[areaIndex,2:8])

    if resp3 in ['4','5','6']:
        sTensor = np.reshape(sTensor, (len(latAll),6))
    #==========================================================================================
    
    # continue to next event if there is faulty data
    if np.isnan(np.sum(sTensor)):
        continue

    # get start and end time for aftershock catalog
    start_dateTime, end_dateTime, EVresp = funcFile.getDateTime(eventFile, slip_file, days)

    # check if event is present in events file
    if EVresp:
        # get aftershocks in buffer region

        AScataBuffer, ASresp = funcFile.getAScataBuff(resp, as_file_name, pickleFile, start_dateTime, end_dateTime, polyBuffer)

        # check if aftershock data is present then calculate ROC values==================================================
        if ASresp:

            # load ROC calc class
            preROC = funcFile.loadROCdata(latAll, lonAll, depAll, AScataBuffer, Mcut)

            # depending on response calculate ROC values
            if resp3 in ['1', '2', '3']:
                ROCdata, stress = preROC.ROCcfs(CFSall)
                # convert stress to Mpa
                stress = stress/1e6
            
            elif resp3 == '4':
                ROCdata, stress = preROC.ROCvm(np.array(sTensor), [strike[0], dip[0], rake[0], f], skempton)
            
            elif resp3 == '5':
                ROCdata, stress = preROC.ROCms(np.array(sTensor))

            elif resp3 == '6':
                ROCdata, stress = preROC.ROCvms(np.array(sTensor))

            else:
                print "errrr..."
                quit()
        #==================================================================================================================

            # save ROC and stress values==========================================

            # ROC
            ROCfname = "%s/%s/%s.roc" %(ROC_dir_path, sType[int(resp3)-1], slip_file.split('.')[0])
            fROC = open(ROCfname, 'w')

            for row in ROCdata:
                fROC.write('%9.6f  %9.6f  %5.2f\n' %(row[0], row[1], row[2]))

            # Stress
            # stressDir = '%s/%s/%s' %(stress_dir_path, sType[int(resp3)-1], slip_file.split('.')[0])
            stressFile = '%s/%s/%s.pkl' %(stress_dir_path, sType[int(resp3)-1], slip_file.split('.')[0])
            strDf = pd.DataFrame()
            # if not os.path.exists(stressDir):
            #     os.makedirs(stressDir)
            
            n = 0
            for d in np.unique(depAll):

                latb = []
                lonb = []
                stressb = []
                for i in range(sum(areaIndex)):
                    if n == 0:
                        latb.append(latAll[n+i])
                        lonb.append(lonAll[n+i])
                    stressb.append(stress[n+i])
                if n == 0:
                    strDf['lat'] = latb
                    strDf['lon'] = lonb
                strDf['D'+str(d)] = stressb

                n += i
            strDf.to_pickle(stressFile)

print
print bcol.OKGREEN + "ROC analysis completed for " + bcol.ENDC + sType[int(resp3)-1]
print