#!/bin/bash

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