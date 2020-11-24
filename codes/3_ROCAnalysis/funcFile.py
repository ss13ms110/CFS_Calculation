# functions used in calcROC.py

import numpy as np
from shapely.geometry import Polygon, LinearRing, Point
import matplotlib.pyplot as plt
import pandas as pd
from numpy import linalg as LA
import random as rn
import os
import datetime as dt

# For coloured outputs in terminal
class bcol:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def dist(lat1, lon1, lat2, lon2):
        """
        Surface distance (in km) between points given as [lat,lon]
        """
        R0 = 6367.3        
        D = R0 * np.arccos(
            np.sin(np.radians(lat1)) * np.sin(np.radians(lat2)) +
            np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.cos(np.radians(lon1-lon2)))
        return D

def dist3D(lat1, lon1, z1,  lat2, lon2, z2):
    """
    distance (in km) between points given as [lat,lon, depth]
    """
    R0 = 6367.3        
    D = R0 * np.arccos(
        np.sin(np.radians(lat1)) * np.sin(np.radians(lat2)) +
        np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.cos(np.radians(lon1-lon2)))
    r = np.sqrt(np.square(D) + np.square(z1-z2))
    return r

def subfinder(mylist, col, pattern):
    return [x for x in mylist if x[col] == pattern]

# function to read srcmod file
def read_fsp_file(inname):

    Nfault = []
    lat = []
    lon = []
    z = []
    slip = []
    rake = []
    strike = []
    dip = []
    L = []
    W = []
    fin = open(inname, "r")
    for line in fin:
        column = line.split()
        Ncol = len(column)
        if column:
            Ncol = len(column)
            if column[0] == "%":
                for n in range(Ncol):
                    if column[n].find('/', 5) > 0:
                        date = column[n].split('/')
                        mm = int(date[0])
                        dd = int(date[1])
                        yy = int(date[2])
                    if column[n] == 'Loc' and Ncol >= n+4 and column[n+3] == '=':
                        lat_hypo = float(column[n+4])
                        lon_hypo = float(column[n+7])
                        z_hypo = float(column[n+10])
                    if column[n] == 'Mw' and Ncol >= n+2 and column[n+1] == '=':
                        M = float(column[n+2])
                    if column[n] == 'STRK' and Ncol >= n+2 and column[n+1] == '=':
                        strike_all = float(column[n+2])
                    if column[n] == 'DIP' and Ncol >= n+2 and column[n+1] == '=':
                        dip_all = float(column[n+2])
                    if column[n] == 'RAKE' and Ncol >= n+2 and column[n+1] == '=':
                        rake_all = float(column[n+2])
                    if column[n] == 'STRIKE' and Ncol >= n+2 and column[n+1] == '=' :
                        strike_all = float(column[n+2])
                    if column[n] == 'DIP' and Ncol >= n+2 and column[n+1] == '=' :
                        dip_all = float(column[n+2])
                    if column[n] == 'Dx' and Ncol >= n+2 and column[n+1] == '=' :
                        Dx = float(column[n+2])
                    if column[n] == 'Dz' and Ncol >= n+2 and column[n+1] == '=' :
                        Dz = float(column[n+2])
                    if column[n] == 'Nsbfs' and column[n+1] == '=':
                        Nfault.append(int(column[n+2]))
            else:
                strike.append(strike_all)
                dip.append(dip_all)
                lat.append(float(column[0]))
                lon.append(float(column[1]))
                z.append(float(column[4]))
                if Ncol >= 7:
                    rake.append(float(column[6]))
                else:
                    rake.append(rake_all)
    fin.close()
    return lat_hypo, lon_hypo, z_hypo, yy, mm, dd, M, Nfault, strike, dip, rake, lat, lon, z

# function to get corners of a fault
def getCorners(Nfault, lat, lon, z):
    na = 0
    nb = 0
    la_corners = [[]]*len(Nfault)
    lo_corners = [[]]*len(Nfault)
    i = 0
    for nf in Nfault:
        nb += nf
        n_depths = len(np.unique(z[na:nb]))
        n_strk = int(nf/float(n_depths))

        la_corners_tmp = [lat[na+0], lat[na+n_strk-1], lat[na+nf-1], lat[na+nf-n_strk]]
        lo_corners_tmp = [lon[na+0], lon[na+n_strk-1], lon[na+nf-1], lon[na+nf-n_strk]]
        
         
        la_corners[i] = la_corners[i] + la_corners_tmp
        lo_corners[i] = lo_corners[i] + lo_corners_tmp
        i += 1
        na += nf
    return la_corners, lo_corners

# get buffer around the fault
def createBuffer(las, los, distDeg):
    poly = Polygon()
    for Nplane in range(len(las)):
        coords = [(los[Nplane][0],las[Nplane][0]), (los[Nplane][1],las[Nplane][1]), (los[Nplane][2],las[Nplane][2]), (los[Nplane][3],las[Nplane][3])]

        lr = LinearRing(coords)
        polyTmp = Polygon(lr).buffer(distDeg)
        poly = poly.union(polyTmp)

    return np.array(poly.exterior.xy), poly

# get lat lon inside the buffer
def getRegion(latlon, polyBuffer):

    areaIndex = []
    for i in range(len(latlon)):
        pt = Point(latlon[i,1], latlon[i,0])

        areaIndex.append(polyBuffer.contains(pt))
    return areaIndex

# get start and end time for aftrshocks
def getDateTime(eventFile, slip_file, days):
    EVresp = False
    start_dateTime, end_dateTime = None, None
     # load event file
    eventfileID = open(eventFile, 'r')
    eventrows = eventfileID.readlines()[1:]
    eventfileID.close()

    eventData = []
    for line in eventrows:
        eventData.append(line.split()[0:13])
    
    # find mainshock stats
    substring = subfinder(eventData, 12, slip_file)

    if len(substring) != 0:
        evYR = int(substring[0][0])
        evMN = int(substring[0][2])
        evDY = int(substring[0][3])
        evHR = int(substring[0][4])
        evMI = int(substring[0][5])
        evSE = float(substring[0][6])
        
        # change here for aftershock begin and end time [days - seconds]
        start_dateTime = dt.datetime(evYR, evMN, evDY, evHR, evMI, int(evSE)) + dt.timedelta(seconds=1)
        end_dateTime = start_dateTime + dt.timedelta(days=days) - dt.timedelta(seconds=1)

        EVresp = True
    
    return start_dateTime, end_dateTime, EVresp

# get aftershocks in given duration and buffer zone
def getISCcatalog(pickleFile, startTime, endTime, polyBuffer):
    print "Loading ISC picke file..."
    df = pd.read_pickle(pickleFile)
    print "Pickle file loaded!!"

    #convert strings to numbers and fill in blanks with NaNs
    df[['magnitude', 'depth', 'latitude', 'longitude']]=df[['magnitude', 'depth', 'latitude', 'longitude']].apply(pd.to_numeric, errors='coerce')

    #quality filtering
    df = df.replace(r'^\s*$', np.nan, regex=True)
    df = df[(df.magnitude.notnull()) & (df.depth.notnull()) & (df.latitude.notnull()) & (df.longitude.notnull()) & (df.magnitude_type.notnull()) & (df.magnitude_author.notnull())]

    #make sure UTM works
    df = df[(df.latitude < 84) & (df.latitude > -80)]
    df = df[(df['depth']<>0.0)]

    # time filtering
    df = df[(df.datetime>startTime) & (df.datetime<endTime)]
    df.replace('\s+', '',regex=True,inplace=True)

    #magnitude type filtering 
    df = df[df['magnitude_type'].isin(['Mb','mb','MB', 'mB', 'ML','Ml','ml','mL','Mw','MW','mW', 'mw', 'Ms','MS', 'mS', 'ms'])]

    aslist = []
    for i, row in df.iterrows():
        Ala = row['latitude']
        Alo = row['longitude']
        Amg = row['magnitude']
        Ade = row['depth']

        pt = Point(Alo, Ala)

        if polyBuffer.contains(pt):
            aslist.append([Ala, Alo, Amg, Ade])
    
    asarray = np.array(aslist)


    return asarray

# depending on response read or load AS catalog
def getAScataBuff(resp, as_file_name, pickleFile, start_dateTime, end_dateTime, polyBuffer):
    AScataBufferNew = []
    ASresp = False
    if resp in ["p", "P"]:
        ASid = open(as_file_name, 'w')
        AScataBuffer = getISCcatalog(pickleFile, start_dateTime, end_dateTime, polyBuffer)
        
        if len(AScataBuffer) != 0:
            ASresp = True
            for i in range(len(AScataBuffer)):
                ASid.write('%8.3f  %8.3f  %4.1f  %5.1f\n' %(AScataBuffer[i,0], AScataBuffer[i,1], AScataBuffer[i,2], AScataBuffer[i,3]))
        else:
            print ("No aftershock data received for this event...")

    else:
        if os.path.isfile(as_file_name) and os.stat(as_file_name).st_size != 0:
            AScataBuffer = np.loadtxt(as_file_name)
            if len(AScataBuffer.shape) == 1:
                ASresp = False
            else:
                ASresp = True
        else:
            AScataBuffer = []
            print ("No aftershock data received for this event...")

    return AScataBuffer, ASresp


# define class to preload lat, lon, and aftershock data
class loadROCdata:
    '''
    AScata = [lat, lon, mag, depth]
    '''

    def __init__(self, lat, lon, depth, AScata, Mcut):
        
        self.Mcut = Mcut
        self.Mmax = np.full(len(lat), -10)
        # calculated the position of AS
        for i in range(len(AScata[:,0])):
            r = dist3D(AScata[i,0], AScata[i,1], AScata[i,3], lat , lon, depth)

            indx = np.argmin(r)

            if AScata[i,2] > self.Mmax[indx]:
                self.Mmax[indx] = AScata[i,2]

    def _ROCvalues(self, CFS):
        ROCdata = []

        for CFSth in np.sort(CFS):
            TP = len(CFS[((CFS>=CFSth) & (self.Mmax >=self.Mcut))])
            FP = len(CFS[((CFS>=CFSth) & (self.Mmax < self.Mcut))])
            FN = len(CFS[((CFS <CFSth) & (self.Mmax >=self.Mcut))])
            TN = len(CFS[((CFS <CFSth) & (self.Mmax < self.Mcut))])    

            if FN > 0 or TP > 0:
                TPR = TP/float(TP + FN)
                FPR = FP/float(TN + FP)
                ROCdata.append([FPR, TPR, CFSth])
        
        return ROCdata

    def _cmbfix(self, strike, dip, rake, sxx, syy, szz, sxy, syz, szx, f, skempton):
        '''
        calculate Coulomb Stress for fixed receiver (based on cmbfix.f) 
        
        Input:	
        stress tensor, friction coefficient
        receiver orientation parameter (strike, dip and rake)
        return:												    
        Coulomb stress (cmb) 
        '''

        DEG2RAD = np.pi/180.0
        
        s11 = sxx
        s12 = sxy
        s13 = szx
        s21 = sxy
        s22 = syy
        s23 = syz
        s31 = szx
        s32 = syz
        s33 = szz

        strike_radian = strike * DEG2RAD
        dip_radian    = dip * DEG2RAD
        rake_radian   = rake * DEG2RAD

        ns1 =  np.sin(dip_radian) * np.cos(strike_radian + 0.5 * np.pi)
        ns2 =  np.sin(dip_radian) * np.sin(strike_radian + 0.5 * np.pi)
        ns3 = -np.cos(dip_radian)

        rst1 = np.cos(strike_radian)
        rst2 = np.sin(strike_radian)
        rst3 = 0.0

        rdi1 = np.cos(dip_radian) * np.cos(strike_radian + 0.5 * np.pi)
        rdi2 = np.cos(dip_radian) * np.sin(strike_radian + 0.5 * np.pi)
        rdi3 = np.sin(dip_radian)
        
        ts1 = rst1 * np.cos(rake_radian) - rdi1 * np.sin(rake_radian)
        ts2 = rst2 * np.cos(rake_radian) - rdi2 * np.sin(rake_radian)
        ts3 = rst3 * np.cos(rake_radian) - rdi3 * np.sin(rake_radian)

        sigg = 0.0
        tau  = 0.0
        sigg += ns1 * s11 * ns1
        tau  += ts1 * s11 * ns1    
        sigg += ns2 * s21 * ns1
        tau  += ts2 * s21 * ns1
        sigg += ns3 * s31 * ns1
        tau  += ts3 * s31 * ns1
        
        sigg += ns1 * s12 * ns2
        tau  += ts1 * s12 * ns2    
        sigg += ns2 * s22 * ns2
        tau  += ts2 * s22 * ns2
        sigg += ns3 * s32 * ns2
        tau  += ts3 * s32 * ns2
        
        sigg += ns1 * s13 * ns3
        tau  += ts1 * s13 * ns3    
        sigg += ns2 * s23 * ns3
        tau  += ts2 * s23 * ns3
        sigg += ns3 * s33 * ns3
        tau  += ts3 * s33 * ns3

        p=-skempton*(sxx+syy+szz)/3.0
        cmb = tau + f * (sigg + p) 
        return cmb

    def _sVM(self, sTensor, sdrf0, sk):
        '''
        sdrf0 = strike, dip, rake, friction list
        sTensor[0,1,2,3,4,5] = Sxx, Syy, Szz, Sxy, Syz, Szx
        S11, S12, S13, S21, S22, S23, S31, S32, S33 = Sxx, Sxy, Szx, Sxy, Syy, Syz, Szx, Syz, Szz
        '''
        sTensor = sTensor/1e6
        Sxx, Syy, Szz, Sxy, Syz, Szx = sTensor.T
        
        # # sd in strike, dip, rake and f
        # std = (25, 25, 25, 0.1)
        # Nsample = 1500
        std = (30, 30, 30, 0.1)
        Nsample = 1500

        # get gaussian distribution of strike, dip, rake and friction
        sdrfDist = []
        for i, sdrf in enumerate(sdrf0):
            sdrfDist.append(np.random.normal(sdrf, std[i], Nsample))
        sdrfDist = np.array(sdrfDist)

        svm = []
        for i in range(len(Sxx)):
            CFS = self._cmbfix(sdrfDist[0], sdrfDist[1], sdrfDist[2], Sxx[i], Syy[i], Szz[i], Sxy[i], Syz[i], Szx[i], sdrfDist[3], sk)
            svm.append(np.sum(CFS[(CFS>0.0)])/(float(Nsample)))
            # CFSn = CFS[(CFS>0.0)]
            # svm.append(np.sum(CFSn)/(len(CFSn)))
        
        return np.array(svm)

    def _sMS(self, sTensor):
        sTensor = sTensor/1e6
        Sxx, Syy, Szz, Sxy, Syz, Szx = sTensor.T

        sms = []
        for i in range(len(Sxx)):
            s11, s12, s13, s21, s22, s23, s31, s32, s33 = Sxx[i], Sxy[i], Szx[i], Sxy[i], Syy[i], Syz[i], Szx[i], Syz[i], Szz[i]

            w, v = LA.eig(np.array([[s11, s12, s13], [s21, s22, s23], [s31, s32, s33]]))

            x1 = w[0]
            x3 = w[2]
            sms.append(np.absolute(x1-x3)/2.)

        return np.array(sms)

    def _sVMS(self, sTensor):
        sTensor = sTensor/1e6
        Sxx, Syy, Szz, Sxy, Syz, Szx = sTensor.T

        svms = []
        for i in range(len(Sxx)):
            s11, s12, s13, s21, s22, s23, s31, s32, s33 = Sxx[i], Sxy[i], Szx[i], Sxy[i], Syy[i], Syz[i], Szx[i], Syz[i], Szz[i]

            w, v = LA.eig(np.array([[s11, s12, s13], [s21, s22, s23], [s31, s32, s33]]))

            IN1 = w[0] + w[1] + w[2]
            IN2 = w[0]*w[1] + w[1]*w[2] + w[0]*w[2]

            svms.append(np.sqrt(IN1*IN1 - 3*IN2))

        return np.array(svms)

    # ROC calculating main attributes here==============
    def ROCcfs(self, CFSval):
        ROCdata = self._ROCvalues(CFSval)

        return ROCdata, CFSval

    def ROCvm(self, sTensor, sdrf0, sk):
        stressVal = self._sVM(sTensor, sdrf0, sk)

        ROCdata = self._ROCvalues(stressVal)

        return ROCdata, stressVal

    def ROCms(self, sTensor):
        stressVal = self._sMS(sTensor)

        ROCdata = self._ROCvalues(stressVal)

        return ROCdata, stressVal

    def ROCvms(self, sTensor):
        stressVal = self._sVMS(sTensor)

        ROCdata = self._ROCvalues(stressVal)

        return ROCdata, stressVal
    #====================================================