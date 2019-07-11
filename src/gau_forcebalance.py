import os
import pandas as pd
import numpy as np

# Costants
d       = {'30': 'Zn', '8': 'O', '1': 'H', '6': 'C', '7': 'N', '16': 'S', '12': 'Mg',
           '11': 'Na', '28': 'Ni', '20': 'Ca', '':'', '17':'Cl','35':'Br' };

d_atms  = ['Zn','C','H','O','Mg','N','Na','Ni','Ca','S','Fe','Cl','Br']



# ------------------- FORCE BALANCE PARSER ----------------------------------------#
def Gau_Forcebalance(filepath, ener_flag, job, name_file = ['qdata.txt', 'all.gro']):
    """
        Function which parse a Gaussian output file and extract 
        Energy, Forces, Coordinate and Moments from the log file
        
        Inputs
            filename
            name_file
            ener_flag
            job
        
        Returns
            file-dataframe for each type
    """
    energy = -1
    force, gro, moment, c, fa                      = [], [], [], [], []
    dipole, quadrupole, quadrupole_A, quadrupole_B  = [], [], [], []
    
    flag_coord, flag_forces      = False, False
    flg_dip, flg_trace, flg_quad = False, False,False
    
    flags = {'ONIOM':'ONIOM: extrapolated energy', 
             'Counterpoise':'complexation',
             'SCF':'SCF Done',
             'BSSE':'Counterpoise corrected energy'}
    
    if ener_flag not in flags.keys():
        print("Choose a proper flag for the energies {}".format(ener_flag))
        return -1
    
    
    
    if os.path.isfile(filepath) == False:
        print("The file {} doesnt exists, please check it".format(filepath[filepath.rfind('/') + 1:]))
        return -1
    else:
        with open(filepath, 'r') as f:
            for line in f:
                # Energy
                if flags[ener_flag] in line:
                    if ener_flag == 'Counterpoise':
                        energy = float(line.strip().split()[-3])
                    elif ener_flag == 'BSSE':
                        energy = float(line.strip().split()[-1])
                    else:
                        energy = float(line.strip().split()[4])
                
                # Start Parsing Forces
                elif line.strip()[36:] == 'Forces (Hartrees/Bohr)':
                    flag_forces = True
                    
                # Stop Parsing Forces
                elif 'Search for a local minimum.' in line:
                    flag_forces = False
                
                # Start Parsing Coordinate
                elif 'Symbolic Z-matrix:' in line:
                    flag_coord = True;
                    
                # Stop Parsing Coordinate
                elif 'NAtoms=' in line:
                    flag_coord = False; 
                    
                if flag_forces == True:
                    if len(line.strip().split()) == 5 and line.strip().split()[1] in d:  # Check the atomic number:
                        force.append([ float(line.strip().split()[2]),
                                       float(line.strip().split()[3]),
                                       float(line.strip().split()[4])])
                            
                
                if flag_coord == True:
                        if len(line.strip().split()) == 4:
                            tmp_at = line.strip().split()[0]
                            for index_let, letter in enumerate(tmp_at):
                                if letter == '(':
                                    lett_atm = tmp_at[0:index_let]        
                            if lett_atm in d_atms:
                                c.append([float(line.strip().split()[1]),
                                          float(line.strip().split()[2]),
                                          float(line.strip().split()[3])])
                                gro.append([lett_atm,
                                              float(line.strip().split()[1]),
                                              float(line.strip().split()[2]),
                                              float(line.strip().split()[3])])

        try:
            os.chdir('./resultsfb')
                
        except:
            os.mkdir('./resultsfb')
            os.chdir('./resultsfb')
        
        ## Check qdata.txt
        if os.path.isfile(name_file[0]) == False:
                with open(name_file[0], 'x') as f:
                    f.write("")
         
        # Check all.gro
        if os.path.isfile(name_file[1]) == False:
                with open(name_file[1], 'x') as f:
                    f.write("")
                    
        if ener_flag == 'Counterpoise' :
            convert_energy = 627.509 # from kcal/mol to Hartree
        else :
            convert_energy = 1
        
        # Write qdata.txt
        with open(name_file[0], 'a') as qd:
            qd.write("JOB {}\n".format(job))
            
            c     = np.array(c)
            force = np.array(force)
            cq    = c.flatten()
            force = force.flatten()
            
            qd.write("COORDS  ")
            for _c in cq:
                 qd.write("{:.4f} ".format(_c*0.1))
        
            qd.write("\nENERGY {} \n".format(str(energy*convert_energy)))
            
            qd.write("FORCES ")
            for _f in force:
                 qd.write("{:.9f} ".format(_f))
            qd.write("\n \n")
         
        conta_atom = 1
        with open(name_file[1], 'a') as allgro:
            allgro.write("Coordinate from {}\n".format(filepath[filepath.rfind('/') + 1:]))
            allgro.write("   {}\n".format(int(len(cq)/3)))
           
            for _gro in gro:
                allgro.write("    1SOL {:>4}{:>4}    {:.4f}    {:.4f}     {:.4f}\n".format(
                                                 _gro[0].replace('O','OW').replace('H','HW'),
                                                 conta_atom,
                                                 _gro[1]*0.1,
                                                 _gro[2]*0.1,
                                                 _gro[3]*0.1))
                conta_atom +=1
            allgro.write("  0.000000   0.000000   0.000000\n")
        
        os.chdir(filepath[:filepath.rfind('/')])