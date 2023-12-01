#! /bin/bash

for f in small14-preprocess small14-valid-iso small14-valid-aniso small14_fig_iso small14_fig_aniso HR46-preprocess HR45-NumChk HR46-ref-process T144-preprocess T144-ref-process
do
    jupyter nbconvert --to notebook --execute --inplace $f.ipynb
done
