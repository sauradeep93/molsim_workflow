# molsim_workflow
A computational workflow which takes as input a CIF (Crystallographic Information File) of a MOF (Metal Organic Framework) and then performs the following operations:

1. Preparation of lammps input file (for energy minimization) using lammps_interface (https://github.com/peteboyd/lammps_interface)
2. Energy minimzation of the MOF using lammps (https://github.com/lammps/lammps)
3. Partial charge assignment of the MOF using EQeq method (https://github.com/danieleongari/EQeq)
4. GCMC (Grand Canonical Monte Carlo) simulations using RASPA (https://github.com/iRASPA/RASPA2)

Please make sure that all the respective programs are installed in your system. 
