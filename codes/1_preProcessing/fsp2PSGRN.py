import numpy as np
import os
import mpl_toolkits
mpl_toolkits.__path__.append('/usr/lib/python2.7/dist-packages/mpl_toolkits/')
from mpl_toolkits.basemap import Basemap

# FUNCTIONS
def read_fsp_file(inname):

    Nfault = []
    lat = []
    lon = []
    z = []
    dip = []
    fin = open(inname, "r")
    for line in fin:
        column = line.split()
        Ncol = len(column)
        if column:
            Ncol = len(column)
            if column[0] == "%":
                for n in range(Ncol):
                    if column[n] == 'DIP' and Ncol >= n+2 and column[n+1] == '=':
                        dip_all = float(column[n+2])
            else:
                dip.append(dip_all)
                z.append(float(column[4]))
    fin.close()
    return dip, z

# MAIN
slp_path = './workingData/srcmod'
psgrn_path = './workingData/PSGRN-INPUTFILES'
crust2_path = './../codes/packages/CRUST_2.0'
pwd = os.getcwd()
cataFile = './workingData/catalog.dat'

catalog = open(cataFile, 'r')
bm = Basemap()

#prams
nr = 61
r1 = 0.0
r2 = 600.0
zs1 = 0.0
sr = 4.0
dp1 = 2.5
dp2 = 50.0
ddp = 5.0
dr_min = ((r2-r1)/float(nr-1))*(2/(1+sr))

# vel prams
etak = 0.0E+00
etakL = 1.0E+17
etam = 0.0E+00
etamL = 1.0E+19
alpha = 1.0
alphaL = 0.5


print "Converting SRCMOD to inputs for PSGRN code..."
print 
for line in catalog:
    la = float(line.split()[3])
    lo = float(line.split()[4])
    fname = line.split()[-1]
     
    slpfile = '%s/%s' %(slp_path, fname)

    dip, z = read_fsp_file(slpfile)
    zs2 = int(max(z)+10)
    nzs = int((zs2-zs1)/(dr_min*np.sin(np.deg2rad(dip[0]))))

    # flag for continent or oceanic coordinates
    flg = int(bm.is_land(lo, la))

    # get velocity model
    os.chdir(crust2_path)
    ftmp = open('tmp.sh', 'w')
    ftmp.write('%s\n' %("dummy=`./getCN2point <<EOF"))
    ftmp.write('%f,%f\n' %(la, lo))
    ftmp.write('%s' %("EOF`"))
    ftmp.close()
    # run this script
    os.system("./tmp.sh")

    # open file "outcr"
    outcr = open('outcr','r')
    dummy = outcr.readline()
    dummy = outcr.readline()
    mantleRow = outcr.readline()
    outrows = outcr.readlines()[2:]
    outcr.close()

    # change dir
    os.chdir(pwd)

    i = 0
    vel_arr = np.empty([7,4])
    for row in outrows:
        vel_arr[i][0] = float(row.split()[0])
        vel_arr[i][1] = float(row.split()[1])
        vel_arr[i][2] = float(row.split()[2])
        vel_arr[i][3] = float(row.split()[2])
        i+=1
    

    vel_tmp = []
    j = 0
    for i in range(7):
        
        if i == 2:
            if vel_arr[2][0] != 0 or vel_arr[i+1][0] != 0:
                vel_tmp.append(float(vel_arr[i][0]+vel_arr[i+1][0]))
                vel_tmp.append(float((vel_arr[i][1]+vel_arr[i+1][1])/2.0))
                vel_tmp.append(float((vel_arr[i][2]+vel_arr[i+1][2])/2.0))
                vel_tmp.append(float((vel_arr[i][3]+vel_arr[i+1][3])/2.0))
                j+=1
        elif i == 4 or i == 5 or i == 6:
            vel_tmp.append(float(vel_arr[i][0]))
            vel_tmp.append(float(vel_arr[i][1]))
            vel_tmp.append(float(vel_arr[i][2]))
            vel_tmp.append(float(vel_arr[i][3]))
            j+=1

    vel_mdl = np.reshape(vel_tmp, (j,4))

    dp = dp1
    while dp <= dp2:

        # open file to write output
        outfile = '%s/psgrn_%s_%s_km.inp' %(psgrn_path, fname[:-4], dp)
        fout = open(outfile, 'w')

        fout.write('%s\n' %("#============================================================================="))
        fout.write('%s\n' %("# This is input file of FORTRAN program \"psgrn2019\" for computing responses"))
        fout.write('%s\n' %("# (Green's functions) of a multi-layered viscoelastic halfspace to point"))
        fout.write('%s\n' %("# dislocation sources buried at different depths. All results will be stored in"))
        fout.write('%s\n' %("# the given directory and provide the necessary data base for the program"))
        fout.write('%s\n' %("# \"pscmp2019\" for computing time-dependent deformation, geoid and gravity changes"))
        fout.write('%s\n' %("# induced by an earthquake with extended fault planes via linear superposition."))
        fout.write('%s\n' %("#"))
        fout.write('%s\n' %("# written by Rongjiang Wang"))
        fout.write('%s\n' %("# Helmholtz Centre Potsdam"))
        fout.write('%s\n' %("# GFZ German Research Centre for Geosciences"))
        fout.write('%s\n' %("# e-mail: wang@gfz-potsdam.de"))
        fout.write('%s\n' %("#"))
        fout.write('%s\n' %("# Last modified: Potsdam, Feb, 2019"))
        fout.write("%s\n" %("#"))

        fout.write('%s\n' %("# PARAMETERS FOR SOURCE-OBSERVATION CONFIGURATIONS"))
        fout.write('%s\n' %("# ================================================"))

        fout.write('     %5.1f     %1d\n' %(dp, flg))
        fout.write(' %4d %7.1f %7.1f %5.1f\n' %(nr, r1, r2, sr))
        fout.write(' %4d %7.1f %7.1f\n' %(nzs, zs1, zs2))

        fout.write('%s\n' %("# ================================================"))
        fout.write("%s\n" %("#"))

        fout.write('%s\n' %("# PARAMETERS FOR TIME SAMPLING"))
        fout.write('%s\n' %("# ============================"))
        fout.write(' %5d  %7.1f\n' %(1, 1))
        fout.write('%s\n' %("# ============================"))
        fout.write("%s\n" %("#"))

        fout.write('%s\n' %("# PARAMETERS FOR WAVENUMBER INTEGRATION"))
        fout.write('%s\n' %("# ====================================="))
        fout.write(' %5.2f\n' %(0.05))
        fout.write(' %5.2f\n' %(1.00))
        fout.write('%s\n' %("# ====================================="))
        fout.write("%s\n" %("#"))

        fout.write('%s\n' %("# PARAMETERS FOR OUTPUT FILES"))
        fout.write('%s\n' %("# ==========================="))
        fout.write('%s\n' %("'./workingData/psgrn+pscmp_out/tmp_out/'"))
        fout.write('%s\n' %("'uz'  'ur'  'ut'"))
        fout.write('%s\n' %("'szz' 'srr' 'stt' 'szr' 'srt' 'stz'"))
        fout.write('%s\n' %("'tr'  'tt'  'rot' 'gd'  'gr'"))
        fout.write('%s\n' %("# ==========================="))
        fout.write("%s\n" %("#"))

        # write header for model
        fout.write('%s\n' %("#------------------------------------------------------------------------------"))
        fout.write('%d             %s\n' %(len(vel_mdl)*2+1, "|int: no_model_lines;"))
        fout.write('%s\n' %("#------------------------------------------------------------------------------"))
        fout.write('%s\n' %("# no  depth[km]  vp[km/s]  vs[km/s]  rho[kg/m^3] etak[Pa*s] etam[Pa*s] alpha"))
        fout.write('%s\n' %("#------------------------------------------------------------------------------"))

        j = 1
        for i in range(len(vel_mdl)):
            if i == 0:
                fout.write('%d      %4.1f       %6.4f    %6.4f     %6.1f     %0.1E    %0.1E   %5.3f\n' %(j, vel_mdl[i][0]-vel_mdl[i][0], vel_mdl[i][1], vel_mdl[i][2], vel_mdl[i][3]*1000, etak, etam, alpha))
                j+=1
                fout.write('%d      %4.1f       %6.4f    %6.4f     %6.1f     %0.1E    %0.1E   %5.3f\n' %(j, vel_mdl[i][0], vel_mdl[i][1], vel_mdl[i][2], vel_mdl[i][3]*1000, etak, etam, alpha))
                vTmp = vel_mdl[i][0]
                j+=1        
            else:
                fout.write('%d      %4.1f       %6.4f    %6.4f     %6.1f     %0.1E    %0.1E   %5.3f\n' %(j, vTmp, vel_mdl[i][1], vel_mdl[i][2], vel_mdl[i][3]*1000, etak, etam, alpha))
                j+=1
                fout.write('%d      %4.1f       %6.4f    %6.4f     %6.1f     %0.1E    %0.1E   %5.3f\n' %(j, sum(vel_mdl[:i+1,0]), vel_mdl[i][1], vel_mdl[i][2], vel_mdl[i][3]*1000, etak, etam, alpha))
                vTmp = sum(vel_mdl[:i+1,0])
                j+=1
        fout.write('%d      %4.1f       %6.4f    %6.4f     %6.1f     %0.1E    %0.1E   %5.3f\n' %(j, vTmp, float(mantleRow.split()[7]), float(mantleRow.split()[8]), float(mantleRow.split()[9])*1000, etakL, etamL, alphaL))
        
        fout.write('%s' %("#=======================end of input==========================================="))

        fout.close()

        dp += ddp

print "Done!"
print "----------Output saved in 'workingData/PSGRN-INPUTFILES'----------"
print 

    

    
    
    
