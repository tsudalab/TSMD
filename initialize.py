import os

def initialize_reactant():
    os.system('gmx_mpi pdb2gmx -f 1uao_unfold.pdb -o 1uao_processed.gro -ignh -water spce -ff amber99sb')

    # Set the simulation config
    os.system('gmx_mpi editconf -f 1uao_processed.gro -o 1uao_newbox.gro -c -d 1.0 -bt cubic')

    # Add solvate
    os.system('gmx_mpi solvate -cp 1uao_newbox.gro -cs spc216.gro -o 1uao_solv.gro -p topol.top')

    # Add ions
    os.system('gmx_mpi grompp -f ions.mdp -c 1uao_solv.gro -p topol.top -o ions.tpr -maxwarn 5')
    os.system("echo 13 | gmx_mpi genion -s ions.tpr -o 1uao_solv_ions.gro -p topol.top -pname NA -nname CL -np 2")

    # Energy minimization
    os.system('gmx_mpi grompp -f minim.mdp -c 1uao_solv_ions.gro -p topol.top -o em.tpr')
    os.system('gmx_mpi mdrun -deffnm em')

    # Equilibration
    os.system('gmx_mpi grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -r em.gro')
    os.system('gmx_mpi mdrun -deffnm nvt')

    os.system('gmx_mpi grompp -f npt.mdp -c nvt.gro -p topol.top -o 0.tpr -r nvt.gro')
    os.system('gmx_mpi mdrun -deffnm 0')

def initialize_target():
    os.system('gmx_mpi pdb2gmx -f 1uao_model1.pdb -o target_processed.gro -ignh -water spce -ff amber99sb -p target_topol.top')

    # Set the simulation config
    os.system('gmx_mpi editconf -f target_processed.gro -o target_newbox.gro -c -d 1.0 -bt cubic')

    # Add solvate
    os.system('gmx_mpi solvate -cp target_newbox.gro -cs spc216.gro -o target_solv.gro -p target_topol.top')

    # Add ions
    os.system('gmx_mpi grompp -f ions.mdp -c target_solv.gro -p target_topol.top -o target_ions.tpr -maxwarn 5')
    os.system("echo 13 | gmx_mpi genion -s target_ions.tpr -o target_solv_ions.gro -p target_topol.top -pname NA -nname CL -np 2")

    # Energy minimization
    os.system('gmx_mpi grompp -f minim.mdp -c target_solv_ions.gro -p target_topol.top -o target_em.tpr -maxwarn 5')
    os.system('gmx_mpi mdrun -deffnm target_em')

    # Equilibration
    os.system('gmx_mpi grompp -f nvt.mdp -c target_em.gro -p target_topol.top -o target_nvt.tpr -r target_em.gro -maxwarn 5')
    os.system('gmx_mpi mdrun -deffnm target_nvt')

    os.system('gmx_mpi grompp -f npt.mdp -c target_nvt.gro -p target_topol.top -o target_npt.tpr -r target_nvt.gro -maxwarn 5')
    os.system('gmx_mpi mdrun -deffnm target_npt')

initialize_reactant()
initialize_target()
