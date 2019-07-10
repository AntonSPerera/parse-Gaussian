# Updated to the 28 of March, 2018

import os
import numpy as np

# This script reads the a set of output of gaussian force calculation and gives as output forces, coordinates,
# dipole, quadrupole and an energies files that will be read in the next module by a matlab code.

# The script needs an input file: "input_matlab.txt"

inputfile = "input_matlab.txt"
if not (os.path.exists('./' + inputfile)):
    r = 'Error!\n"{}" not found.'.format(inputfile)
    print(r)
    exit()

with open(inputfile,'r') as f:
    lines = f.readlines()
    # Path to the log files
    path = lines[2].split()[0]
    # Path for the output
    out_path = lines[3].split()[0]
    # Flag force
    flag_force = lines[4].split()[0]
    # Flag coordinate
    flag_coordinate = lines[5].split()[0]
    # Flag energy 
    flag_energy = lines[6].split()[0]
    # Flag Counterpoise
    flag_counter = lines[7].split()[0]
    # flag_dipole
    flag_dipole = lines[8].split()[0]
    # flag quadrupole
    flag_quadrupole = lines[9].split()[0]
    # Number of atoms
    n_atoms = int(lines[10].split()[0])
    # Number of water molecules
    n_wat = int(lines[11].split()[0])
    # Topology flag

# Check *.log files extracted from the Gaussian calculation path
if not (os.path.exists(path)):
    r = 'Error!\n"{}" not found.'.format(path)
    print(r)
    exit()

if n_atoms <= 0:
   print('Please write a correct number of atoms as input : {}'.format(str(n_atoms)))
   exit()

diz_flag = {'Y','y','n','N','yes','no','No','Yes'}
flag_force = flag_force.lower()
flag_coordinate = flag_coordinate.lower()
flag_energy = flag_energy.lower()
flag_counter = flag_counter.lower()
flag_dipole = flag_dipole.lower()
flag_quadrupole = flag_quadrupole.lower()

sys = {'polymer','water_solution','hetero_solution'}

if flag_force not in diz_flag or flag_coordinate not in diz_flag or flag_energy not in diz_flag or flag_counter not in diz_flag or flag_dipole not in diz_flag or flag_quadrupole not in diz_flag:
    print('Please insert a correct value for the files flag')
    exit()
# ------------------------------------------------------------------------------------
h2kj= 2625.5;
  # Conversion from Harthree to kiloJoule/mol
d = {'30': 'Zn', '8': 'O', '1': 'H', '6': 'C', '7': 'N', '16': 'S', '12': 'Mg', '11': 'Na', '28': 'Ni', '20': 'Ca', '':'', '17':'Cl','35':'Br' };
d_atms = ['Zn','C','H','O','Mg','N','Na','Ni','Ca','S','Fe','Cl','Br']

tip3p_charges = {'OW': -0.834, 'HW': 0.417}
tip3p_sigma = {'OW': 3.15061e-01, 'HW': 0.0}
tip3p_eps = {'OW': 6.36386e-01, 'HW': 0.0}

forces_atomic_numbers = []
force_coords = []
force = []
forces = False
energy = []
complexation = []
coordinate = []
coord =  []
coord_atomic_numbers = []
quadrupole = []
quadrupole_A = []
quadrupole_B = []
dipole = []
#--------------------------------------------------------------------------------------

os.chdir(out_path)

if os.path.exists('forces'+ str(n_wat) +'.txt'):
    os.remove('forces'+ str(n_wat) +'.txt')
if os.path.exists('coordinates'+ str(n_wat) +'.txt'):
    os.remove('coordinates'+ str(n_wat) +'.txt')
if os.path.exists('energy'+ str(n_wat) +'.txt'):
    os.remove('energy'+ str(n_wat) +'.txt')
if os.path.exists('complexation'+ str(n_wat) +'.txt'):
    os.remove('complexation'+ str(n_wat) +'.txt')
if os.path.exists('dipole'+str(n_wat)+'.txt'):
    os.remove('dipole'+str(n_wat)+'.txt')
if os.path.exists('quadrupole_TrA'+str(n_wat)+'.txt'):
    os.remove('quadrupole_TrA'+str(n_wat)+'.txt')
if os.path.exists('quadrupole_TrB'+str(n_wat)+'.txt'):
    os.remove('quadrupole_TrB'+str(n_wat)+'.txt') 
if os.path.exists('quadrupole' + str(n_wat) + '.txt'):
    os.remove('quadrupole' + str(n_wat) + '.txt')

os.chdir(path)
for filename in os.listdir(path):
        if filename[-4:]=='.log':
                print(filename)
                with open(filename,'r') as f1:
                        lines = f1.readlines();
                forces = False
                flag_coord = False
                flg_quad = False
                flg_dip = False
                flg_trace = False
                if 'Normal termination of Gaussian' not in lines[-1].strip():
                    print('the file {} didn t converged'.format(filename))
                    exit()
                for line in lines:
                    if line.strip()[:26]  == 'ONIOM: extrapolated energy':
                        energy.append(line.strip().split()[4])
                    if line.strip()[:9] == 'SCF Done:':
                        energy.append(line.strip().split()[4])
                    if 'complexation' in line:
                        complexation.append(line.strip().split()[-3])
                    if line.strip()[36:] == 'Forces (Hartrees/Bohr)':
                        forces = True
                        force.append(str(n_atoms) +  '\n')
                        force.append('Title forces forces forces forces (Hartrees/Bohr) \n')
                    if 'Search for a local minimum.' in line:
                        forces = False
                    if forces == True:
                         if len(line.strip().split()) == 5:
                             #print(line.strip().split())
                             if line.strip().split()[1] in d:
                                 prv = [];
                                 forces_atomic_numbers.append(line.strip().split()[1]);
                                 prv.append(line.strip().split()[1])
                                 prv.append(line.strip().split()[2]);
                                 prv.append(line.strip().split()[3]);
                                 prv.append(line.strip().split()[4]);
                                 force_coords.append(prv);
                                 #print(d[prv[0]])
                                 #print(prv[1])
                                 #print(forces_atomic_numbers[-1])
                                 force.append('{:>5}{:>15}{:>15}{:>15}\n'.format(d[prv[0]],prv[1],prv[2],prv[3]))
                    if 'Symbolic Z-matrix:' in line:
                        flag_coord = True;
                        coordinate.append(str(n_atoms)+ '\n')
                        coordinate.append('Title  cooordinates  cooordinates  cooordinates  cooordinates  cooordinates:\n')
                    if 'NAtoms=' in line:
                        flag_coord = False; 
                    if flag_coord == True:
                        if len(line.strip().split()) == 4:
                            tmp_at = line.strip().split()[0]
                            for index_let, letter in enumerate(tmp_at):
                                if letter == '(':
                                    lett_atm = tmp_at[0:index_let]
                                    
                            if lett_atm in d_atms:
                                coordinate.append('{:>5}{:>15}{:>15}{:>15}\n'.format(lett_atm, line.strip().split()[1],line.strip().split()[2],line.strip().split()[3]))
                    if 'Dipole moment (field-independent basis, Debye)' in line:
                        flg_dip = True
                    if flg_dip == True and len(line.split()) == 8 :
                        dipole.append([line.strip().split()[1], line.strip().split()[3], line.strip().split()[5]])  
                        flg_dip = False
                    if 'Quadrupole moment (field-independent basis, Debye-Ang)' in line and 'Traceless' not in line:
                        flg_quad = True
                    if flg_quad == True and len(line.split()) == 6 and 'XX' in line:
                        quadrupole.append([line.strip().split()[1],line.strip().split()[3],line.strip().split()[5]])
                        flg_quad = False
                    if 'Traceless Quadrupole moment (field-independent basis, Debye-Ang):' in line:
                        flg_trace = True   
                    if flg_trace == True and len(line.strip().split()) == 6:
                        if 'XX=' in line:
                            quadrupole_A.append([line.strip().split()[1], line.strip().split()[3], line.strip().split()[5]])
                        if 'XY=' in line:
                            quadrupole_B.append([line.strip().split()[1],line.strip().split()[3],line.strip().split()[5]])
                            flg_trace = False



os.chdir(out_path)
if flag_force[0] == 'y':
    with open('forces'+ str(n_wat) +'.txt','x') as f:
        for ff in force:
            f.write(ff)
if flag_coordinate[0] =='y':
    with open('coordinates'+ str(n_wat) +'.txt','x') as c:
        for cc in coordinate:
            c.write(cc)
if flag_energy[0] == 'y' and flag_counter[0] == 'n':
    with open('energy'+ str(n_wat) +'.txt','x') as e:
        #e.write('Values are in Hartrees\n')
        for ee in energy:
            e.write(ee + '\n')
if flag_energy[0] == 'y' and flag_counter[0] == 'y':
    n = len(energy)
    l_energy = []
    sys = []
    w_1 = []
    i_1 = []
    w_2 = []
    i_2 = []
    
    sys = energy[0:n:5]
    i_1 = energy[1:n:5]
    w_1 = energy[2:n:5]
    i_2 = energy[3:n:5]
    w_2 = energy[4:n:5]
    if len(sys) !=  len(i_1) or  len(i_1)!= len(i_2) or len(i_2) != len(w_1) or len(w_1) != len(w_2):
        print('check the files you have, someone might not have converged')
        exit()
    with open('energy'+ str(n_wat) +'.txt','x') as e:
        e.write('{:>15}{:>15}{:>15}{:>15}{:>15}{:>15}{:>20}\n'.format('sys','water','ion', 'water corr','ion corr','diff','diff corr [Hartrees]'))
        for jj in range(len(sys)):
            diff = float(sys[jj]) - float(w_1[jj]) - float(i_1[jj])
            diff_corr = float(sys[jj]) - float(w_2[jj]) - float(i_2[jj])
            if diff > 0 or diff_corr > 0:
                print('Difference between the global system and each fragment is positive: {} corr {} '.format(str(diff), str(diff_corr)))
            e.write('{:>15}{:>15}{:>15}{:>15}{:>15}{:>15}{:>15}\n'.format(sys[jj],w_1[jj],i_1[jj],w_2[jj],i_2[jj],str(round(diff,8)),str(round(diff_corr,8))))
       


if flag_counter[0] == 'y':
    n = len(complexation)
    with open('complexation'+ str(n_wat) +'.txt','x') as c:
        c.write('{:15}{:15}\n'.format('raw','corrected [kcal/mol]'))
        for ii in range(0,n,2):
            c.write('{:15}{:15}\n'.format(complexation[ii], complexation[ii+1]))

if flag_dipole[0] == 'y':
    with open('dipole' + str(n_wat) + '.txt','x') as ddipol:
        for dip in dipole:
            ddipol.write('{:>15}{:>15}{:>15}\n'.format(dip[0],dip[1],dip[2]))    

if flag_quadrupole[0] == 'y':
    with open('quadrupole_TrA' + str(n_wat) + '.txt','x') as Qa:
        for qqa in quadrupole_A:
            Qa.write('{:>15}{:>15}{:>15}\n'.format(qqa[0],qqa[1],qqa[2]))
    with open('quadrupole_TrB' +str(n_wat) + '.txt','x') as Qb:
        for qqb in quadrupole_B:
            Qb.write('{:>15}{:>15}{:>15}\n'.format(qqb[0],qqb[1],qqb[2]))
    with open('quadrupole' +str(n_wat) + '.txt','x') as Qq:
        for qqq in quadrupole:
            Qq.write('{:>15}{:>15}{:>15}\n'.format(qqq[0],qqq[1],qqq[2]))

























