# script to plot ROC curves
# user inputs required to plot
# specific ROC curves

import numpy as np
import matplotlib.pyplot as plt
import os

# FUNCTION
def pltFunc(fname, ttl):
    plt.tick_params(axis="x", labelsize=18)
    plt.tick_params(axis="y", labelsize=18)
    plt.xlabel("False positive rate", fontsize=20)
    plt.ylabel("True positive rate", fontsize=20)
    plt.legend(loc='best', prop={'size': 14})
    plt.title('ROC for %s' %(ttl))
    plt.savefig("%s/%s.pdf" %(figPath, fname), dpi=500)




# PATHS ==================================
ROCPath = './workingData/ROCData/ROCout'
figPath = './figs/ROC'
# ========================================

# PRAMS ==================================
sType = ['MAS0', 'MAS', 'OOP', 'VM', 'MS', 'VMS']

# ========================================

# user inputs ============================
uResp = raw_input('Single event or ROC average [1/2]: ')
if not uResp in ['1', '2']:
    print "Wrong input %s. Quiting..." %(uResp)
    quit()
print

if uResp == '1':
    srcmodID = raw_input("Enter SRCMOD ID of slip model: ")
    for sT in sType:
        ROCFile = ROCPath + '/' + sT + '/' + srcmodID + '.roc'
        if not os.path.exists(ROCFile):
            print "ROC %s does not exist. Quiting..." %(srcmodID)
            quit()
    print
else:
    plotMode = raw_input('Metric (e.g. MAS0, VM, ...) or leave empty for all plots: ')
    if not plotMode in sType+[""]:
        print "Wrong input %s, Quiting..." %(plotMode)
        quit()
    print
# ========================================

#MAINS
if uResp == '1':
    # initiate figure
    plt.figure(figsize=(8,8))
    
    # loop into stress metric
    for sT in sType:
        ROCFile = ROCPath + '/' + sT + '/' + srcmodID + '.roc'

        FP, TP = np.loadtxt(ROCFile, usecols=(0,1), unpack=True)

        AUC = abs(np.trapz(TP, FP))
        # plot step
        plt.step(FP, TP, label='%s  | AUC = %5.3f' %(sT, AUC), lw=2.0)
    
    pltFunc(srcmodID, srcmodID)

else:
    metList = sType
    allTag = True
    if plotMode in sType:
        metList = [plotMode]
        allTag = False
    
    # for combined plot of all ROC
    allTP = []
    allFP = []
    allAUC = []
    for met in metList:
        # initiate figure
        plt.figure(figsize=(8,8))
        ROCfPath = ROCPath + '/' + met

        ROCfList = os.listdir(ROCfPath)

        TFP = []
        TTP = []
        bins = np.linspace(0, 1, 500)
        cnt = 0
        for ROCf in ROCfList:
            ROCfName = ROCfPath + '/' + ROCf

            FP, TP = np.loadtxt(ROCfName, usecols=(0,1), unpack=True)

            TFP = np.append(TFP, FP)
            TTP = np.append(TTP, TP)
            plt.step(FP,TP, color="grey")
        
        TFPdigi = np.digitize(TFP, bins)
        TPbin_means = [TTP[TFPdigi == i].mean() for i in range(1, len(bins))]
        new_bins = [(a+b)/2.0 for a, b in zip(bins[:len(bins)-1], bins[1:])]
        AUC = np.trapz(TPbin_means, new_bins)

        if allTag:
            allTP.append(TPbin_means)
            allFP.append(new_bins)
            allAUC.append(AUC)

        plt.plot(new_bins, TPbin_means, color="blue", linewidth=3.0, label="Average ROC\nAUC = %5.3f" %(abs(AUC)))
        plt.plot(np.array([0,1]), np.array([0,1]), '--', color="black", linewidth=2.0)

        pltFunc(met, met)

    if allTag:
        plt.figure(figsize=(8,8))
        
        for i in range(len(allTP)):

            plt.plot(allFP[i], allTP[i], label='%s  | AUC = %5.3f' %(sType[i], allAUC[i]))

        pltFunc('ROCmean', 'all metric')

