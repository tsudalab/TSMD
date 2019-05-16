# TSMD
Python implementations of TS-MD and PaCS-MD [Harada, R. and Kitao, A., "Parallel cascade selection molecular dynamics (PaCS-MD) to generate conformational transition pathway," J. Chem. Phys., 139, 07B611 1. (2013)].
Our implementations use GROMACS command for all of the manipulation about MD simulation.
The required version is GROMACS 2016.5.
## USAGE
### equivration
Before running TS-MD or PaCS-MD, we have to do adding ion, energy minimization, nvt equilibration and npt equilibration.
We need the .mdp files to do each manipulation of such preparations. The .mdp files here are intended only for use with this chignolin(1uao) folding sample. You have to set the parameters for each of your tasks.  
```
initialize.py() 
```
### Requrements
Both TS-MD and PaCS-MD program require the files shown below.
- reactant structure (.gro)
- checkpoint file after the equilibration (.cpt)
- topology file (.top)
- target structure (.gro)
- mdp file for short md (.mdp)

### TS-MD
- Options
  - -r : Reactant structure file (.gro).
  - -t : Target structure file (.gro).
  - -top : Topology file (.top).
  - -s : Step size. default is 1000.
  - -c : A parameter of ucb socre. Default is 0.05.
  - -cn : Resume from the previous state. Default is 0.
  - -ntmpi : The number of mpi.
  - -ntomp : The number of open MP
  - -del   : If you want to delete the intermediate trajectory files, set this to 1. default is 0.
  - -thresh : A thereshold (Å) that decide whether one structure is similar to others.
  - -alp : A penalty parameter.
  - -ctype : Choose the way how to calclate UCB. (normal, adaptive, adaptive2)

If you have the reactant structure file as reactant.gro, target structure file as target.gro, and topology file as topol.top
your command is 
```
python pats_md.py -r reactant -t target -top topol
```
If TS-MD doesn't reach the enough RMSD, the intermediate state is preserve as vars.pikcle automatically.
You can resume the process by setting the -cn to 1.

### PaCS-MD
```
python pacs_md.py -r reactant -t target -top topl
```
- Options
  - -r : Reactant structure file (.gro).
  - -t : Target structure file (.gro).
  - -top : Topology file (.top).
  - -s : Step size. Default is 1000.
  - -k : The number of parallelization. Default is 5
  - -cn : Resume from the previous state. Default is 0.
  - -ntmpi : The number of mpi.
  - -ntomp : The number of open MP.
  - -del   : If you want to delete the intermediate trajectory files, set this to 1. default is 0.
