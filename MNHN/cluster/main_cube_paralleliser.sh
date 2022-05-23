#!/bin/bash

pathFiles="/trinity/home/pturk/Projet_MNHN/MNHN/cluster/script_cube_auto_generated"

for file in $pathFiles/cube*.sh
do

sbatch $file

done