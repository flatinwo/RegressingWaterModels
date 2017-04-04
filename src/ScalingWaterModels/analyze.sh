#!/bin/bash
rm -f All.dat
python test-$1.py
mkdir -p Combined-$1
cp Combined-05/scri*.R Combined-$1/.
mv All.dat Combined-$1/.
cd Combined-$1/.
sed -i -e 's/Madrid/Madrid+Princeton/g' All.dat
sed -i -e 's/Evan/Madrid+Princeton/g' All.dat
mv All.dat All-1.dat
R CMD BATCH script-trick.R
sleep 3
open RescaledTSEOS.pdf
cd ..
