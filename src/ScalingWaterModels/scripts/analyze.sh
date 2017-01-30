#!/bin/bash
cd ..
python test-04.py
mv All.dat Combined-03/.
cd Combined-03/.
sed -i -e 's/Madrid/Madrid+Princeton/g' All.dat
sed -i -e 's/Evan/Madrid+Princeton/g' All.dat
mv All.dat All-1.dat
R CMD BATCH script-trick.R
open RescaledTSEOS.pdf
cd ..
