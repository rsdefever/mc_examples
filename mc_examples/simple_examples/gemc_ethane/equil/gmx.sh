#!/bin/bash

python gmx.py

gmx grompp -f em.mdp -c in.gro -p topol.top -o em
gmx mdrun -v -deffnm em

gmx grompp -f eq.mdp -c em.gro -p topol.top -o eq -maxwarn 1
gmx mdrun -v -deffnm eq


