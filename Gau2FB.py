#!/usr/bin/env python


import numpy as np
import os
import argparse

"""
Gau2FB.py is a script that takes in input a set of Gaussian output
".log" files and parse it building the files qdata.txt, all.gro, 
necessary to the ForceBalance procedure to perform a fitting on 
energy and forces evaluated at a quantum mechanical level.



"""

# List of atoms
d_atms  = ["None",'H','He',
            'Li','Be','B','C','N','O','F','Ne',
            'Na','Mg','Al','Si','P','S','Cl','Ar',
            'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
            'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
            'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
            'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
            'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt']

# Dictionary atom number : atom
d       = {'1':'H','2':'He',
           '3':'Li','4':'Be','5':'B','6':'C','7':'N','8':'O','9':'F','10':'Ne',
           '11':'Na','12':'Mg','13':'Al','14':'Si','15':'P','16':'S','17':'Cl','18':'Ar',
           
           '19':'K','20':'Ca','21':'Sc','22':'Ti','23':'V','24':'Cr','25':'Mn','26':'Fe','27':'Co','28':'Ni','29':'Cu',
           '30':'Zn','31':'Ga','32':'Ge','33':'As','34':'Se','35':'Br','36':'Kr',
           
           '37':'Rb','38':'Sr','39':'Y','40':'Zr','41':'Nb','42':'Mo','43':'Tc','44':'Ru','45':'Rh','46':'Pd','47':'Ag','48':'Cd',
           '49':'In','50':'Sn','51':'Sb','52':'Te','53':'I','54':'Xe',
           
           '55':'Cs','56':'Ba','57':'La','58':'Ce','59':'Pr','60':'Nd','61':'Pm','62':'Sm','63':'Eu','64':'Gd','65':'Tb','66':'Dy',
           '67':'Ho','68':'Er','69':'Tm','70':'Yb',
           
           '71':'Lu','72':'Hf','73':'Ta','74':'W','75':'Re','76':'Os','77':'Ir','78':'Pt','79':'Au','80':'Hg','81':'Tl','82':'Pb',
           '83':'Bi','84':'Po','85':'At','86':'Rn',
           
           '87':'Fr','88':'Ra','89':'Ac','90':'Th','91':'Pa','92':'U','93':'Np','94':'Pu','95':'Am','96':'Cm','97':'Bk','98':'Cf',
           '99':'Es','100':'Fm','101':'Md','102':'No','103':'Lr','104':'Rf','105':'Db','106':'Sg','107':'Bh','108':'Hs','109':'Mt'}
    
  
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
    force, gro, c, fa            = [], [], [], []
    flag_coord, flag_forces      = False, False
    
    
    
    
    
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
        
        ## Check if qdata.txt exists
        if os.path.isfile(name_file[0]) == False:
                with open(name_file[0], 'x') as f:
                    f.write("")
         
        # Check if all.gro exists
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
            c     = c.flatten()
            
            force = np.array(force)
            force = force.flatten()
            
            qd.write("COORDS  ")
            for _c in c:
                 qd.write("{:.4f} ".format(_c*0.1))             # From nm to A
        
            qd.write("\nENERGY {} \n".format(str(energy*convert_energy)))
            
            qd.write("FORCES ")
            for _f in force:
                 qd.write("{:.9f} ".format(_f))
            qd.write("\n \n")
         
        conta_atom = 1
        with open(name_file[1], 'a') as allgro:
            allgro.write("Coordinate from {}\n".format(filepath[filepath.rfind('/') + 1:]))
            allgro.write("   {}\n".format(int(len(c)/3)))
            for _gro in gro:
                allgro.write("    1SOL {:>4}{:>4}    {:.4f}    {:.4f}     {:.4f}\n".format(
                                                 _gro[0],
                                                 conta_atom,
                                                 _gro[1]*0.1,
                                                 _gro[2]*0.1,
                                                 _gro[3]*0.1))
                conta_atom +=1
            allgro.write("  0.000000   0.000000   0.000000\n")
        
        os.chdir(filepath[:filepath.rfind('/')])

        
def walk_Gau_Forcebalance(path, ener_flag):
    
    job = 0
    for filename in os.listdir(path):
        if filename[-4:]=='.log':
                Gau_Forcebalance(path + '/'+ filename, ener_flag, job=job)
                job +=1
                

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)





def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--path', type=dir_path)
    parser.add_argument('--ener_flag')
    
    args = parser.parse_args()
    
    path      = args.path
    ener_flag = args.ener_flag
    
    flags = {'ONIOM':'ONIOM: extrapolated energy', 
             'Counterpoise':'complexation',
             'SCF':'SCF Done',
             'BSSE':'Counterpoise corrected energy'}
    
    if ener_flag not in flags.keys():
        print("Choose a proper flag for the energies {}".format(ener_flag))
        return -1
    
    walk_Gau_Forcebalance(path, ener_flag)

if __name__ == "__main__":
    main()