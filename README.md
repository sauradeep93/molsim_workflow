# material_design_&_screening_workflow
A computational workflow to screen materials for carbon capture, hydrogen storage, methane storage and other energy related applications. The workflow takes as input a CIF (Crystallographic Information File) of a MOF (Metal Organic Framework) and then performs the following operations:

(0.) generation of in-silico structures from building block files using ToBaCCo (https://github.com/tobacco-mofs/tobacco_3.0)

1. Preparation of lammps input file (for energy minimization) using lammps_interface (https://github.com/peteboyd/lammps_interface)
2. Energy minimzation of the MOF using lammps (https://github.com/lammps/lammps)
3. Partial charge assignment of the MOF using EQeq method (https://github.com/danieleongari/EQeq)
4. Pore geometry property calculation by Zeo++ (https://github.com/danieleongari/aiida-zeopp)
5. GCMC (Grand Canonical Monte Carlo) simulations using RASPA (https://github.com/iRASPA/RASPA2)

Please make sure that all the respective programs are installed in your system. 

Publications that have used this workflow:
https://pubs.acs.org/doi/10.1021/acsami.1c16220

Workflow contributors: Sauradeep Majumdar, Elias Moubarak

Contact: sauradeep.majumdar@epfl.ch
