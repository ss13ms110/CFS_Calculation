#!/bin/bash

echo "Run PSGRN(1) or PSCMP(2)? [1/2]"
read resp

if [ $resp -eq "1" ]
then

rm -f not_worked_psgrn.list
for f in inputs/PSGRN-INPUTFILES/*
do
file=`echo $f | awk -F"/" '{print $3}'`
id=`echo $f | awk -F"_" '{print $2"_"$3}'`
depth=`echo $f | awk -F"_" '{print $3}'`

cp $f $file

./bin/psgrn2019 <<EOF
$file
EOF

if [ -f outputs/tmp_out/szz.ss ]
then

if [ ! -d outputs/psgrn_out/${id} ]
then
mkdir outputs/psgrn_out/${id}
fi
mv outputs/tmp_out/* outputs/psgrn_out/${id}

else
echo $file >> not_worked_psgrn.list
fi
rm $file
done

elif [ $resp -eq "2" ]
then

rm -f not_worked_pscmp.list
for f in inputs/PSCMP-INPUTFILES/*
do
file=`echo $f | awk -F"/" '{print $3}'`
depth=`echo $f | awk -F"_" '{print $3}'`

cp $f $file

./bin/pscmp2019 <<EOF
$file
EOF

if [ -f outputs/pscmp_out/coseis.dat ]
then
mv outputs/pscmp_out/coseis.dat outputs/pscmp_out/coseis_${file}
else
echo $file >> not_worked_pscmp.list
fi
rm $file

done

else
echo "Wrong Input..."
fi