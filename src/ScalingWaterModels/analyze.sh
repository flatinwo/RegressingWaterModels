#!/bin/bash
rm -f All.dat
python test-09.py
mkdir -p Combined-09
#cp Combined-05/scri*.R Combined-06/.
mv All.dat Combined-09/.
cd Combined-09/.
sed -i -e 's/Madrid/Madrid+Princeton/g' All.dat
sed -i -e 's/Evan/Madrid+Princeton/g' All.dat
mv All.dat All-1.dat
R CMD BATCH script-trick.R
sleep 3
open RescaledTSEOS.pdf
cd ..
