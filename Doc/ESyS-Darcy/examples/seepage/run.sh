#!/bin/bash

#a bash file to run the simulation and create vtk files for visualization in Paraview
#use the command "sh run.sh" in terminal for execution 


if [ -d "results" ]; then
  rm -r results
fi

mkdir results

mpirun -np 2 `which esysparticle` upward_seepage.py

dump2vtk -i results/snap -o results/particles_ -rot -t 0 6 100000

printf "\nSimulation done, results saved in /results\n"
