#!/bin/bash

echo "Run PSGRN(1) or PSCMP(2)? [1/2]"
read resp

# PSGRN=================================================
if [ $resp -eq "1" ]
then
inDir='workingData/PSGRN-INPUTFILES'
outDir='workingData/psgrn+pscmp_out'

rm -f workingData/not_worked_psgrn.list
for f in ${inDir}/*
do
file=`echo $f | awk -F"/" '{print $3}'`
id=`echo $file | awk -F"_" '{print $2"_"$3}'`
depth=`echo $file | awk -F"_" '{print $3}'`

cp $f $file

./../codes/packages/PSGRN+PSCMP/bin/psgrn2019 <<EOF
$file
EOF

if [ -f ${outDir}/tmp_out/szz.ss ]
then

if [ ! -d ${outDir}/psgrn_out/${id} ]
then
mkdir ${outDir}/psgrn_out/${id}
fi
mv ${outDir}/tmp_out/* ${outDir}/psgrn_out/${id}

else
echo $file >> workingData/not_worked_psgrn.list
fi
rm $file
done
#=========================================================

# PSCMP===================================================
elif [ $resp -eq "2" ]
then

# okada or GF ----------------------------------------
echo "Run Okada or GF [O/G]: "
read OGF



if [ $OGF == 'O' ] || [ $OGF == 'o' ]
then
inDir='workingData/PSCMP-INPUTFILES/okada'
outDir='workingData/psgrn+pscmp_out/pscmp_out/okada'
elif [ $OGF == 'G' ] || [ $OGF == 'g' ]
then
inDir='workingData/PSCMP-INPUTFILES/GF'
outDir='workingData/psgrn+pscmp_out/pscmp_out/GF'
else
echo "Wrong Input... [O/G or o/g]"
exit
fi
#-----------------------------------------------------

rm -f workingData/not_worked_pscmp.list
for f in ${inDir}/*
do
file=`echo $f | awk -F"/" '{print $4}'`
depth=`echo $file | awk -F"_" '{print $3}'`
id=`echo $file | awk -F"_" '{print $2}'`

cp $f $file

./../codes/packages/PSGRN+PSCMP/bin/pscmp2019 <<EOF
$file
EOF

if [ -f ${outDir}/coseis.dat ]
then

if [ ! -d ${outDir}/${id} ]
then
mkdir ${outDir}/${id}
fi

mv ${outDir}/coseis.dat ${outDir}/${id}/coseis_${file}
else
echo $file >> workingData/not_worked_pscmp.list
fi
rm $file

done
#============================================================

else
echo "Wrong Input... [1/2]"
fi