#!/bin/bash

pathFiles="/trinity/home/pturk/scriptCube"

for file in $pathFiles/cube*.sh
do

sbatch $file

done