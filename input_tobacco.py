import os
import argparse
import math
import pandas as pd
import numpy as np
import subprocess







def files(extension):
    curr_directory = os.getcwd()
   # new_directory = os.path.join(curr_directory, "cifs")
    list_files = []

    for f in os.listdir(curr_directory):
        if f.endswith('.' + extension):
            list_files.append(f)

    return list_files


list_xyz = files("xyz")




path = os.getcwd()
print(path)


for i in list_xyz:
    print(i)
    #subprocess.call(["lammps-interface", "--minimize", "--cutoff=6", "--fix-metal", i])
    #print("in." + i.split('.')[0])


    molecule = pd.read_csv(i, skiprows=2, delim_whitespace=True,names=['atom', 'x', 'y', 'z'])
    s=[]
    c=1
    for j in range(0, len(molecule.x)):
        s.append(molecule.atom[j] + str(c))
        c=c+1
    print(len(molecule.x))
    print(min(molecule.x), max(molecule.x))


    xl = (max(molecule.x)-min(molecule.x))
    yl = (max(molecule.y)-min(molecule.y))
    zl = (max(molecule.z)-min(molecule.z))
    xc = min(molecule.x) + xl/2
    yc = min(molecule.y) + yl/2
    zc = min(molecule.z) + zl/2


    centre=[xc,yc,zc]

    print(xl,yl,zl)
    a =3 * max(xl,yl,zl)
    print(a)
    b=c=a
    print(a,b,c)
    print(c)
    #alpha=beta=gamma= 90 degrees


    molecule.x = molecule.x - xc
    molecule.y = molecule.y - yc
    molecule.z = molecule.z - zc
 
    u = (1/a)*molecule.x
    v =(1/b)*molecule.y
    w =(1/c)*molecule.z

    ofile = open(i.split('.xyz')[0]+'.cif', 'w+')

    print('CIF file', file = ofile)
    
    print('loop_', file = ofile)
    print('_symmetry_equiv_pos_as_xyz', file = ofile)
    print('\'x, y, z\'', file = ofile)
    
    print('_cell_length_a', a, file = ofile)
    print('_cell_length_b',b, file = ofile)
    print('_cell_length_c', c,file = ofile)
    #print(a,b,c)
    print('_cell_angle_alpha','90', file = ofile)
    print('_cell_angle_beta','90' ,file = ofile)
    print('_cell_angle_gamma','90' ,file = ofile)
    
    print('loop_', file = ofile)
    
    print('_atom_site_type_symbol', file = ofile)
    print('_atom_site_label', file = ofile)
    print('_atom_site_fract_x', file = ofile)
    print('_atom_site_fract_y', file = ofile)
    print('_atom_site_fract_z', file = ofile)
    
    
    
    
    s1=[]
    s2=[]
    s3=[]
    s4=[]
    s5=[]
    count=1
    for k in range(0, len(molecule.x)):
        s1.append(molecule.atom[k])
        s2.append(molecule.atom[k])
        s3.append(u[k])
        s4.append(v[k])
        s5.append(w[k])
        count=count+1
    
    #s5[0]
        print("%s %8s %12.6f %10.6f %10.6f" %(s1[k], s2[k], u[k]+0.5, v[k]+0.5, w[k]+0.5), file =ofile)
    #print(s1,s2,s3,s4,s5, file=ofile)

