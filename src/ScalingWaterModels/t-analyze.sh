#!/bin/bash
rm -f All.dat
python test-10.py
mkdir -p Combined-10
#cp Combined-05/scri*.R Combined-10/.
mv All.dat Combined-10/.
cd Combined-10/.
sed -i -e 's/Madrid/Madrid+Princeton/g' All.dat
sed -i -e 's/Evan/Madrid+Princeton/g' All.dat
mv All.dat All-1.dat
R CMD BATCH script-trick.R
sleep 5
open RescaledTSEOS.pdf
cd ..
