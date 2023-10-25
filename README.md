# material_design_&_screening_workflow
A computational workflow to screen materials for carbon capture, hydrogen storage, methane storage and other energy related applications. The workflow takes as input a CIF (Crystallographic Information File) of a MOF (Metal Organic Framework) and then performs the following operations (Step 0 can be used to build a cif file and then pass it on to Step 1):

0. (Optional) Generation of in-silico structures from building block files using ToBaCCo (https://github.com/tobacco-mofs/tobacco_3.0)

1. Preparation of lammps input file (for energy minimization) using lammps_interface (https://github.com/peteboyd/lammps_interface)
2. Energy minimzation of the MOF using lammps (https://github.com/lammps/lammps)
3. Partial charge assignment of the MOF using EQeq method (https://github.com/danieleongari/EQeq)
4. Pore geometry property calculation by Zeo++ (https://github.com/danieleongari/aiida-zeopp)
5. GCMC (Grand Canonical Monte Carlo) simulations using RASPA (https://github.com/iRASPA/RASPA2)

Please make sure that all the respective programs are installed in your system. 


# Usage
Each of the python files can be run independently or together.

'python opt_localpc.py' will take in a batch of cif files and compute steps 1 to 5. This version is ideal to run on local computers for number of structures in the order of 100's.

For evaluating 1000's of structures, it might be benefecial to run these scripts individually on a remote cluster, and for that we suggest looking at the examples provided in the examples_run.zip folder





If you find our work useful, please cite us.

Publications that have used this workflow:
1. https://pubs.acs.org/doi/10.1021/acsami.1c16220
2. https://doi.org/10.26434/chemrxiv-2023-71mjq-v2

Workflow contributors: Sauradeep Majumdar, Elias Moubarak 

Contact: sauradeep.majumdar@epfl.ch
