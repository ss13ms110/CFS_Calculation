import numpy as np
import sys
import os
import shapely.geometry
import utm
import pyproj
import matplotlib.pyplot as plt

def dist(lat1, lon1, lat2, lon2):
        """
        Surface distance (in km) between points given as [lat,lon]
        """
        R0 = 6367.3        
        D = R0 * np.arccos(
            np.sin(np.radians(lat1)) * np.sin(np.radians(lat2)) +
            np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.cos(np.radians(lon1-lon2)))
        return D


def read_fsp_file(inname):
    # fin = open(inname, "r")
    # Nfault = 0
    # for line in fin:
    #     column = line.split()
    #     if column:
    #         if column[0] == "%":
    #             continue                    
    #         else:
    #             Nfault += 1
    # fin.close()

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
                L.append(Dx)
                W.append(Dz)
                lat.append(float(column[0]))
                lon.append(float(column[1]))
                z.append(float(column[4]))
                slip.append(float(column[5]))
                if Ncol >= 7:
                    rake.append(float(column[6]))
                else:
                    rake.append(rake_all)
    fin.close()
    return lat_hypo, lon_hypo, z_hypo, yy, mm, dd, M, Nfault, strike, dip, rake, rake_all, L, W, lat, lon, z, slip 

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
    return np.array(la_corners), np.array(lo_corners)


# function to get lat-lon range 
def getRange(la_corners, lo_corners, dxy, laHy, loHy):
    la1, la2 = np.min(la_corners), np.max(la_corners)
    lo1, lo2 = np.min(lo_corners), np.max(lo_corners)
    
    # extend lat lon to 200 kms from corners
    la1, la2 = la1 - 200/lat2km, la2 + 200/lat2km
    lo1, lo2 = lo1 - 200/lon2km, lo2 + 200/lon2km


    # get zone from lat-lon Hypo
    _,_,zone,_ = utm.from_latlon(laHy, loHy)
    p = pyproj.Proj(proj='utm', zone=zone, ellps='WGS84')
    x1, y1 = p(lo1, la1)
    x2, y2 = p(lo2, la2)

    stepsize = dxy*1000 # 5 km grid step size

    
    
    # Iterate over 2D area

    lats = []
    lons = []
    x = x1
    i = 0
    while x < x2:
        y = y1
        j = 0
        while y < y2:
            p1, p2 = p(x, y, inverse=True)
            lons.append(p1)
            lats.append(p2)
            y += stepsize
            j += 1
        x += stepsize
        i += 1

    
    return i, j, min(lats), max(lats), min(lons), max(lons)


# inputs to generate Okada vs GF tensors
OGF = raw_input("Generate inputs for Okada solution or GF solution [O/G]: ")
print 
# main:
slipmodeldir = './workingData/srcmod'
pscminputdir = './workingData/PSCMP-INPUTFILES'

# parameters
eta = 0.4
dxy=5
B = 0.00
S1 = 0.1E+07   
S2 = 0.0E+07    
S3 = -0.2E+07
lbda = 3.0E+10
mu = 3.0E+10
lat2km = 110.57
depth = np.arange(2.5, 50, 5)

# Okada or GF
if OGF in ['O', 'o']:
    OGFDir = 'okada'
elif OGF in ['G', 'g']:
    OGFDir = 'GF'
else:
    print "Wrong input [O/G or o/g]: %s" %(OGF)
    quit()



# open file to write out_file names

print "Converting SRCMOD to inputs for PSCMP code..."
print 

sliplist = os.listdir(slipmodeldir)
for i in range(len(sliplist)):
    slipfile = '%s/%s' % (slipmodeldir, sliplist[i])
    lat_hypo, lon_hypo, z_hypo, yy, mm, dd, M, Nfault, strike, dip, rake, rake0, L, W, lat, lon, z, slip = read_fsp_file(slipfile)
    
    
    # get corners -------------------------------------------------------------------------
    la_corners, lo_corners = getCorners(Nfault, lat, lon, z)
    
    lon2km = dist(lat_hypo, lon_hypo-0.5, lat_hypo, lon_hypo+0.5)

    ndx, ndy, laty1, laty2, lonx1, lonx2 = getRange(la_corners, lo_corners, dxy, lat_hypo, lon_hypo)

    # --------------------------------------------------------------------------------------
    
    z=np.array(z)
    nA = 0
    for nB in Nfault:
        nB=int(nB)
        nB = nA+nB
        if z[nA] < 0:
            z[nA:nB] = z[nA:nB] - z[nA]
        nA = nB
    name = sliplist[i][:-4]
    Nfault = np.array(Nfault)

    
    for dp in depth:
        outname = '%s/%s/pscmp_%s_%s_km.inp' % (pscminputdir, OGFDir, name, dp)
        fout=open(outname,"w")
        fout.write('# OBSERVATION ARRAY\n')
        fout.write('# =================\n')
        fout.write(" 2\n")
        fout.write(' %d    %.1f  %.1f\n' % (ndy, laty1, laty2))
        fout.write(' %d    %.1f  %.1f\n' % (ndx, lonx1, lonx2))
        fout.write('# OUTPUT\n')
        fout.write('# ======\n')
        fout.write(' 0   -0.0068  0.3848  -0.9205    ![option turned off] cosines of LOS direction from ground to Evisat satellite orbit\n')
        fout.write('#\n')
        fout.write('# =>    select (1/0) output for Coulomb failure stress changes (only for snapshots, see below): icfs,\n')
        fout.write('#       friction, Skempton ratio, strike, dip, and rake angles [deg] describing the uniform regional master\n')
        fout.write('#       fault mechanism, the uniform regional principal stresses: sigma1, sigma2 and sigma3 [Pa] in\n')
        fout.write('#       arbitrary order (the orietation of the pre-stress field will be derived by assuming that the master\n')
        fout.write('#       fault is optimally oriented according to Coulomb failure criterion)\n')
        fout.write('#\n')
        fout.write('#       if this option is selected (icfs = 1), the snapshots will include additional data:\n')
        fout.write('#\n')
        fout.write('#       CFS_Max = increase of maximum CFS (can be oriented differently before and after the earthquake)\n')
        fout.write('#       CFS_Mas = increase of CFS for the given master fault mechanism\n')
        fout.write('#       CFS_Mas_Opt = increase of CFS on the master fault at the postseismic optimal rake direction\n')
        fout.write('#       Sigma_Mas = increase of normal stress on the master fault\n')
        fout.write('#       Rake_Mas_Opt = postseismic optimal rake on the master fault\n')
        fout.write('#       CFS_Opt = increase of CFS for the postseismic optimal mechanism (in 3D space)\n')
        fout.write('#       Sigma_Opt_1/2 = increase of normal stress on the first/second postseismic optimal faults\n')
        fout.write('#       Strike_Opt_1/2, Dip_Opt_1/2, Rake_Opt_1/2 = the first/second postseismic optimal focal mechanisms\n')
        fout.write('#\n')
        fout.write('#    Note: the first 3D optimally orieted fault is the closest to the master fault.\n')
        fout.write('#\n')
        fout.write('#=======================================================================================================\n')
        # print ' %d     %5.3f  %5.3f  %7.3f   %6.3f   %7.3f    %3.1E   %3.1E   %3.1E\n' %(1,eta,B,strike[0],dip[0],rake[0],S1,S2,S3)
        fout.write(' %d     %5.3f  %5.3f  %7.3f   %6.3f   %7.3f    %3.1E   %3.1E   %3.1E\n' %(1,eta,B,strike[0],dip[0],rake0,S1,S2,S3))
        fout.write("'./workingData/psgrn+pscmp_out/pscmp_out/%s/'\n" %(OGFDir))
        fout.write(' %d                %d                %d\n' %(0,0,0))
        fout.write(' %s    %s    %s\n' %("'U_north.dat'","'U_east.dat'","'U_down.dat'"))
        fout.write(' %d            %d            %d            %d            %d            %d\n' %(0,0,0,0,0,0))
        fout.write(' %s   %s   %s   %s   %s   %s\n' %("'S_nn.dat'","'S_ee.dat'","'S_dd.dat'","'S_ne.dat'","'S_ed.dat'","'S_dn.dat'"))
        fout.write(' %d               %d               %d               %d               %d\n' %(0,0,0,0,0))
        fout.write(' %s    %s    %s    %s    %s\n' %("'Tilt_n.dat'","'Tilt_e.dat'","'Rotation.dat'","'geoid.dat'","'Gravity.dat'"))
        fout.write(' %d\n' % (1))
        fout.write('      %4.2f  %s              %s\n' %(0.00,"'coseis.dat'","|coseismic"))
        fout.write('#    %5.2f  %s             %s\n' %(10.00,"'days_10.dat'","|coseismic + postseismic of the first 10 days"))
        fout.write('#    %5.2f  %s            %s\n' %(30.00,"'months_1.dat'","|coseismic + postseismic of the first 1 month"))
        fout.write('#    %5.2f  %s            %s\n' %(90.00,"'months_3.dat'","|coseismic + postseismic of the first 3 months"))
        fout.write('#   %6.2f  %s             %s\n' %(365.00,"'years_1.dat'","|coseismic + postseismic of the first 1 year"))
        fout.write('#  %7.2f  %s            %s\n' %(3650.00,"'years_10.dat'","|coseismic + postseismic of the first 10 years"))
        fout.write('#=======================================================================================================\n')
        fout.write("# GREEN'S FUNCTION DATABASE\n")
        fout.write('# =========================\n')
        fout.write('# 1. selection (0/1) of earth structure model (0 = homogeneous elastic halfspace, i.e., only for\n')
        fout.write('#    co-seismic elastic Okada solutions, 1 = multi-layered viscoelastic halfspace, for which numerical\n')
        fout.write('#    co- and post-seismic Greens functions calculated with psgrn2019 are needed): iesmodel\n')
        fout.write('#\n')
        fout.write('#    IF (iesmodel = 0 for analytical Okada solutions) THEN\n')
        fout.write('#  \n')
        fout.write('# 2. observation depth [km] and the two Lame coefficients lambda and mue [Pa]\n')
        fout.write('#\n')
        fout.write('#    ELSE IF (iesmodel = 1 for numerical psgrn2019 Greens functions) THEN\n')
        fout.write('#\n')
        fout.write('# 2. directory where the Greens functions are sts2005KASHMI01KONC.fspored: grndir\n')
        fout.write('# 3. file names (without extensions!) for the 13 Greens functions:        \n')
        fout.write('#=======================================================================================================\n')
        
        if OGF in ['O', 'o']:
            fout.write(' %d\n' %(0))
            fout.write(' %4.1f  %3.1E  %3.1E\n' %(dp,lbda,mu))
            fout.write("#'./workingData/psgrn+pscmp_out/psgrn_out/%s_%s/'\n" %(name, dp))
            fout.write('# %s  %s  %s\n' %("'uz'","'ur'","'ut'"))
            fout.write('# %s  %s  %s  %s  %s  %s\n' % ("'szz'","'srr'","'stt'","'szr'","'srt'","'stz'"))
            fout.write('# %s  %s  %s  %s  %s\n' %("'tr'","'tt'","'rot'","'gd'","'gr'"))

        elif OGF in ['G', 'g']:
            fout.write(' %d\n' %(1))
            fout.write('# %4.1f  %3.1E  %3.1E\n' %(dp,lbda,mu))
            fout.write("'./workingData/psgrn+pscmp_out/psgrn_out/%s_%s/'\n" %(name, dp))
            fout.write(' %s  %s  %s\n' %("'uz'","'ur'","'ut'"))
            fout.write(' %s  %s  %s  %s  %s  %s\n' % ("'szz'","'srr'","'stt'","'szr'","'srt'","'stz'"))
            fout.write(' %s  %s  %s  %s  %s\n' %("'tr'","'tt'","'rot'","'gd'","'gr'"))

        fout.write('#=======================================================================================================\n')
        fout.write("# RECTANGULAR SUBFAULTS\n")
        fout.write('# =====================\n')
        fout.write("#===============================================================================\n")
        fout.write("# n_faults\n")
        fout.write("#-------------------------------------------------------------------------------\n")
        fout.write("%d\n" % (sum(Nfault)))
        fout.write("#-------------------------------------------------------------------------------\n")
        fout.write("# n   O_lat   O_lon    O_depth length  width strike dip   np_st np_di start_time\n")
        fout.write("# [-] [deg]   [deg]    [km]    [km]     [km] [deg]  [deg] [-]   [-]   [day]\n")
        fout.write("#     pos_s   pos_d    slp_stk slp_dip open\n")
        fout.write("#     [km]    [km]     [m]     [m]     [m]\n")
        fout.write("#-------------------------------------------------------------------------------\n")
        for i in range(len(z)):
            fout.write(" %1d    %8.4f %8.4f %7.4f %6.2f %6.2f %6.2f %6.2f %3d %3d %5.2f\n" % (i+1,lat[i], lon[i], z[i], L[i], W[i], strike[i], dip[i], 1, 1, 0.0))
            slip_di = (-1)*slip[i]*np.sin(np.deg2rad(rake[i]))
            slip_st = slip[i]*np.cos(np.deg2rad(rake[i]))
            fout.write(" %6.2f  %6.2f  %6.3f  %6.3f  %6.3f\n" % (0.0, 0.5*W[i], slip_st, slip_di, 0.0))
        fout.write('#================================end of input===================================')
        fout.close()

print "Done!"
print "----------Output saved in 'workingData/PSCMP-INPUTFILES'----------"
print 