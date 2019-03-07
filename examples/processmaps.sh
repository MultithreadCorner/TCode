#!/bin/bash
'path to executables'/loadmaps_cuda  -E 'path to input map'/my_electron_mobility.txt -e 'path to input map'/my_electric_field.txt -w 'path to input map'/my_weighting_field.txt -H 'path to input map'/my_hole_mobility.txt --efield_sc 2 --wfield_sc 2 --hmob_sc 2 --emob_sc 2 --emob_scale 1e8  --hmob_scale 1e8 --efield_scale 1e-4 --wfield_scale 1e-4
