import os
import shutil
import subprocess
import numpy as np
import math





# Calculate the minimum unit cells needed for GCMC simulation
def cell_units(lens, angs, co):
    # lens and cutoff (co) are given in A and angs in degrees
    # unpack parameters
    a = lens[0]
    b = lens[1]
    c = lens[2]
    alpha = math.radians(angs[0])
    beta = math.radians(angs[1])
    gamma = math.radians(angs[2])
    # Start the calculations (https://en.wikipedia.org/wiki/Fractional_coordinates)
    K = c*(math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma)
    V = math.sqrt(1-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2+2*math.cos(alpha)*math.cos(beta)*math.cos(gamma))
    Minv = [[a, 0, 0], [b*math.cos(gamma), b*math.sin(gamma), 0] , [c*math.cos(beta), K, c*V/math.sin(gamma)]]
    axb = np.cross(Minv[0],Minv[1])
    bxc = np.cross(Minv[1],Minv[2])
    cxa = np.cross(Minv[2],Minv[0])
    x = np.dot(Minv[0],bxc)/np.linalg.norm(bxc)
    y = np.dot(Minv[1],cxa)/np.linalg.norm(cxa)
    z = np.dot(Minv[2],axb)/np.linalg.norm(axb)
    xmin = math.ceil(2*co/x)
    ymin = math.ceil(2*co/y)
    zmin = math.ceil(2*co/z)
    return xmin, ymin, zmin



# Extract the edge lengths and angles from the .cif file
def extract_geometry(structure):
    with open (structure + '.cif', 'r') as fi:
        data = fi.readlines()
        for line in data:
            if "_cell_length_a" in line:
                a = float(line.split()[1])
            elif "_cell_length_b" in line:
                b = float(line.split()[1])
            elif "_cell_length_c" in line:
                c = float(line.split()[1])
            elif "_cell_angle_alpha" in line:
                alpha = float(line.split()[1])
            elif "_cell_angle_beta" in line:
                beta = float(line.split()[1])
            elif "_cell_angle_gamma" in line:
                gamma = float(line.split()[1])
        lens = [a, b, c]
        angs = [alpha, beta, gamma]
        fi.close()
    unitcell = cell_units(lens, angs, 12)
    return unitcell


# Create a simulation.input file for GCMC calculations
def GCMC(cycles, structure, unitcell, T, P, molecule, block):
    # T and P are in (K) and (Pa) respectively
    with open("simulation.input","w") as fo:
        str_out = ""
        str_out += "SimulationType               MonteCarlo\n"
        str_out += "NumberOfCycles               %i\n"%(cycles)
        str_out += "NumberOfInitializationCycles %i\n"%(cycles/10)
        str_out += "PrintEvery                   %i\n"%(cycles/10)
        str_out += "PrintPropertiesEvery         %i\n"%(cycles/10)
        str_out += "RestartFile                  no\n\n"
        str_out += "Forcefield                   local\n\n"
        str_out += "Framework               0\n"
        str_out += "FrameworkName           %s\n"%(structure)
        str_out += "UseChargesFromCIFFile   yes\n"
        str_out += "UnitCells               %i %i %i\n"%(unitcell[0],unitcell[1],unitcell[2])
        str_out += "ExternalTemperature     %s\n"%str(T)
        str_out += "ExternalPressure        %s\n\n"%str(P)
        str_out += "Component 0 MoleculeName                 %s\n"%(molecule)
        str_out += "            MoleculeDefinition           local\n"
        if (block == "yes"):
            str_out += "            BlockPockets                 %s\n"%(block)
            str_out += "            BlockPocketsFilename         %s\n"%(structure)
        str_out += "            TranslationProbability       1.0\n"
        str_out += "            RotationProbability          1.0\n"
        str_out += "            ReinsertionProbability       1.0\n"
        str_out += "            SwapProbability              1.0\n"
        str_out += "            CreateNumberOfMolecules      0\n"
        fo.write(str_out)
        fo.close()

# Extract loading and heat of adsorption from a GCMC output file
def extract_GCMC(structure, unitcell, T, P):
    with open ("output_"+ structure + "_" + str(unitcell[0]) + "." +
               str(unitcell[1]) + "." + str(unitcell[2]) + "_" +
               str(T) + ".000000_" + str(P) + ".data", 'r') as fi:
        data = fi.readlines()
        for line in data:
            if "Enthalpy of adsorption:" in line:
                Q_line = data[data.index(line) + 10]
                Q = float(Q_line.split()[0])
            elif "Average loading absolute [mol/kg framework]" in line:
                L = float(line.split()[5])
    fi.close()
    # L (mmol/g) and Q (J/mmol)
    return L, Q

# Create job.bash file to run RASPA
def slurm(structure):
    with open("job.bash","w") as fo:
        str_out = ""
        str_out += "#!/bin/bash\n\n"
        str_out += "#SBATCH --no-requeue\n"
        str_out += "#SBATCH --wait\n"
        str_out += "#SBATCH --job-name   %s\n"%(structure)
        str_out += "#SBATCH --get-user-env\n"
        str_out += "#SBATCH --output     jobout.txt\n"
        str_out += "#SBATCH --error      joberr.txt\n"
        str_out += "#SBATCH --nodes      1\n"
        str_out += "#SBATCH --ntasks     1\n"
        str_out += "#SBATCH --partition  serial\n"
        str_out += "#SBATCH --time       05:00:00\n\n"
        str_out += "export RASPA_DIR=/home/kjablonk/RASPA/simulations/\n"
        str_out += "export DYLD_LIBRARY_PATH=/home/kjablonk/RASPA/simulations/lib\n"
        str_out += "export LD_LIBRARY_PATH=/home/kjablonk/RASPA/simulations/lib\n\n"
        str_out += "/home/kjablonk/RASPA/simulations/bin/simulate 'simulation.input'\n"
        fo.write(str_out)
        fo.close()

# Run RASPA calculations
def run_RASPA(structure, unitcell, T, molecule, block, cycles, P):
    GCMC(cycles, structure, unitcell, T, P, molecule, block)
    #slurm(structure)
    subprocess.call('simulate')

# Read RASPA's output
def read_RASPA(structure, unitcell, T, molecule, P):
    current_dir = os.getcwd()
    new_dir = current_dir + "/Output/System_0"
    os.chdir(new_dir)

    L, Q = extract_GCMC(structure, unitcell, T, P)
    os.chdir(current_dir)
    return L, Q

# Obtain blocking spheres using zeo++
def zeo(structure):
    subprocess.call(["network", "-ha", "-block", "1.525", "500", structure + ".cif"])
