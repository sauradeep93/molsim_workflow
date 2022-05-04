#!/home/sauradeepmajumdar/anaconda3/bin/python

#import argparse
import math
#import json
#import numpy as np 
#import networkx as nx
#from scipy.spatial import distance
#import sys
#import math
import os
import subprocess
import shutil
#from collections import deque
#import matplotlib.pyplot as plt



def cifs(extension):
    curr_directory = os.getcwd()
   # new_directory = os.path.join(curr_directory, "cifs")
    list_cifs = []

    for f in os.listdir(curr_directory):
        if f.endswith('.' + extension):
            list_cifs.append(f)

    return list_cifs


list_cifs = cifs("cif")


# running lammps_interface to generate input and data files for structure optimization in lammps..Suitable editing/modifications of the generated file shave been done for convenience


path = os.getcwd()
print(path)


for i in list_cifs:
    print(i)
    subprocess.call(["lammps-interface", "--minimize", "--cutoff=6", "--fix-metal", i])
    print("in." + i.split('.')[0])
        
   

    
    with open("data." + i.split('.')[0], "r") as f0:
        data = f0.read().splitlines()
        c=c1=c2=0
        for line in data:
            c = c+1
            if "Masses" in line:
                c1 = c
            if "Bond Coeffs" in line:
                c2 =c

        c4=(c2-3)-(c1+1)+1 
        a=c1+1
       
        atom_type=" "
        while(a>c1 and a<(c2-2)):
            atom_type +=(data[a].split('#')[1])
            a=a+1

    with open("in." + i.split('.')[0], "r") as f1:
        data = f1.readlines()
        data1 =  data[0:22]
        data2 =  data[22:len(data)]
        f1.close() 
    
    os.remove("in." + i.split('.')[0])
    with open("in." + i.split('.')[0], "w") as f2:
        str_out = ""
        for j in data1:
            str_out+=j
        #str_out+="dump" + "            " + i.split('.')[0] + "_xyzmov"+ " " +"all xyz 100" + " " +  i.split('.')[0] + "_mov.xyz\n"
        #str_out+="dump_modify" + "     "  +   i.split('.')[0] + "_xyzmov element" + " " + atom_type + "\n"
        str_out+="dump            d1 all custom 100" + " " + "dump." + i.split('.')[0] + " "+ "element xs ys zs\n"
        str_out+="dump_modify     d1 element" + atom_type +"\n"
       
        for k in data2:
            str_out+=k
        str_out+="undump"+ "          "+ i.split('.')[0] + "_xyzmov"
        str_out+="\n"
        f2.write(str_out)
        f2.close()
       
    
    # running energy minimization in lammps using input and data files generated from lammps_interface


    with open("trial1.sh", "w+") as f3:
         str_out=""
         str_out+="#!/bin/bash\n"
         str_out+="lmp_serial < in." + i.split('.')[0]
         f3.write(str_out)
         f3.close()

    
    
    subprocess.call(["sh" , "trial1.sh"])

    #subprocess.call(["lmp_serial", "<", "in." + i.split('.')[0]])



    #EXTRACTING THE COORDINATES OF THE LAST FRAME FROM DUMP FILES OF LAMMPS AND CONVERTING IT TO CIF FILES


    with open("dump." + i.split('.')[0], "r") as f4, open(i.split('.')[0]+".dump","w+") as f5:
        data=f4.readlines()
        f4.close()
        n_atoms= int(data[3])
        #print(n_atoms)
        #for row in f4.readlines()[-(n_atoms+9):]:
        for row in data[-(n_atoms+9):]:
            f5.write(row)
        f4.close()
        f5.close()

    file_path_opt = os.path.join(path, "charge_cifs", i.split('.')[0]+"_opt.cif")
    with open(i.split('.')[0]+".dump","r") as f6, open(file_path_opt,"w+") as f7:
        data = f6.read().splitlines()
        f6.close()

        for line in data:
            if line=='':
                data.remove(line)


        yy=' '.join(data[5:8])
        #for i in yy.split():
        #  print(i,float(i))
        values=[float(i) for i in yy.split()]
        print(values)
        print(values[0])


        xlo_bound=values[0]
        xhi_bound=values[1]
        xy=values[2]
        ylo_bound=values[3]
        yhi_bound=values[4]
        xz=values[5]
        zlo_bound=values[6]
        zhi_bound=values[7]
        yz=values[8]
        xlo=xlo_bound-min(0.0,xy,xz,xy+xz)
        xhi=xhi_bound-max(0.0,xy,xz,xy+xz)
        ylo=ylo_bound-min(0.0,yz)
        yhi=yhi_bound-max(0.0,yz)
        zlo=zlo_bound
        zhi=zhi_bound
        lx=xhi-xlo
        ly=yhi-ylo
        lz=zhi-zlo
        a=lx
        b= math.sqrt(math.pow(ly,2)+ math.pow(xy,2))
        c= math.sqrt(math.pow(lz,2)+ math.pow(xz,2) + math.pow(yz,2))
        alph= math.acos(((xy*xz)+(ly*yz))/(b*c))
        bet= math.acos(xz/c)
        gamm = math.acos(xy/b)
        alpha =math.degrees(alph)
        beta =math.degrees(bet)
        gamma =math.degrees(gamm)
        print(a,b,c,alpha,beta,gamma)
   #end of lattice parameters calculation



        count = 0
        count_lim = 0
        loop_list = []
        atom_site_list = []

        for line in data:
            count_lim += 1
        print(count_lim)
        data = data [0:count_lim]## need to put correct line number
        flags = data[9:len(data)] #first 9 lines(0 -8) do  not contain the coordinates
        print(len(data))
        print(len(flags))

        str_out=""
        str_out+=("CIF file\n")
        str_out+=("loop_\n")
        str_out+=("_symmetry_equiv_pos_as_xyz\n")
        str_out+=("\'x, y, z\'\n")
        str_out+=("_symmetry_space_group_name_H-M    P1\n")
        str_out+=("_symmetry_Int_Tables_number       1\n")
        str_out+=("_cell_length_a"+ " " + str(a) +"\n")
        str_out+=("_cell_length_b" + " " + str(b) + "\n")
        str_out+=("_cell_length_c" + " " + str(c) + "\n")
        str_out+=("_cell_angle_alpha" + " " + str(alpha) +"\n")
        str_out+=("_cell_angle_beta" + " " + str(beta) + "\n")
        str_out+=("_cell_angle_gamma" + " " + str(gamma) +"\n")
        str_out+=("loop_\n")
        str_out+=("_atom_site_type_symbol\n")
        str_out+=("_atom_site_label\n")
        str_out+=("_atom_site_fract_x\n")
        str_out+=("_atom_site_fract_y\n")
        str_out+=("_atom_site_fract_z\n")



        for j in range(0, len(flags)):
            a = flags[j].split()
    
            if "Ti6+4" in a[0]:
                a.append("Ti")

            if "Co6+3" in a[0]:
                a.append("Co")

            if "Cu3+1" in a[0]:
                a.append("Cu")

            if "Ni4+2" in a[0]:
                a.append("Ni")

            if "O_2" in a[0]:
                a.append("O")

            if "O_R" in a[0]:
                a.append("O")
  
            if "O_3" in a[0]:
                a.append("O")
           
            if "H_" in a[0]:
                a.append("H")
    
            if "N_R" in a[0]:
               a.append("N")

            if "C_3" in a[0]:
               a.append("C")
 
            if "C_R" in a[0]:
               a.append("C")

            if "C_2" in a[0]:
               a.append("C")

            if "C_1" in a[0]:
               a.append("C")
    
            if "N_1" in a[0]:
               a.append("N")
            
            if "F_" in a[0]:
               a.append("F")
            
            if "Br" in a[0]:
               a.append("Br")
            
            if "S_3+6" in a[0]:
               a.append("S")
            
            if "S_R" in a[0]:
               a.append("S")
            
            
            
         
            str_out+=((a[4])+ "  " +(a[4])+ "  "+(a[1])+ "  "+(a[2])+"  "+ (a[3]))
            str_out+="\n"
        
        str_out+=("_end")
        str_out+="\n"
        f7.write(str_out)
        f7.close()
    

    
    #file_path_charge = os.path.join(path, "charge_cifs", i.split('.')[0]+"_opt.cif")
    
    #charge calculation using eqeq method

  #  subprocess.call(["./eqeq", file_path_opt])
