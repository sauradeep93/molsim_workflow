from gcmc import *


def cifs(extension):
    directory = os.getcwd()
    structure_list = []

    for f in os.listdir(directory):
        if f.endswith('.' + extension):
            structure_list.append(f)

    return structure_list



structure_list= cifs("cif")


for structure in structure_list:
    #structure=structure_list[i].rsplit('.', 1)[0]
    print(structure)
    print(structure.rsplit('.', 1)[0])


    #Get the cell parameters
    unitcell=extract_geometry(structure.rsplit('.',1)[0])
    print(unitcell)
   
   
   
    #Run zeo ++ calculations - blocking spheres
    zeo(structure.rsplit('.', 1)[0])

    # Check if the block file is empty or not
    with open(structure.rsplit('.',1)[0] + ".block") as f9:
        if "0" in f9.readline():
            block = "no"
        else:
            block = "yes"
        f9.close()





    molecule = "CO2"       # You can change it to N2
    T_ads = 298            # Adsorption temperature (K)
    T_des = 363            # Desorption temperature (K)
    P_ads = 100000         # Adsorption pressure (Pa)
    P_des = 10000          # Desorption pressure (Pa)
    cycles = 10000          # Number of cycles

    run_RASPA(structure.rsplit('.',1)[0], unitcell, T_ads, molecule, block, cycles, P_ads)
    #run_RASPA(structure.rsplit('.',1)[0], unitcell, T_des, molecule, block, cycles, P_des)

    L_ads, Q_ads = read_RASPA(structure.rsplit('.',1)[0], unitcell, T_ads, molecule, P_ads)
    #L_des, Q_des = read_RASPA(structure.rsplit('.',1)[0], unitcell, T_des, molecule, P_des)

    with open("output_data","a+") as fo:
        str_out = ""
        str_out += "Structure:" + structure.rsplit('.',1)[0]+".cif\n"
        #str_out += "T_ads: %s\n"%str(T_ads)
        #str_out += "P_ads: %s\n"%str(P_ads)
        str_out += "L_ads: %s\n"%str(L_ads)
        str_out += "Q_ads: %s\n"%str(Q_ads)
        #str_out += "T_des: %s\n"%str(T_des)
        #str_out += "P_des: %s\n"%str(P_des)
        #str_out += "L_des: %s\n"%str(L_des)
        #str_out += "Q_des: %s\n"%str(Q_des)
        fo.write(str_out)
        fo.close()
