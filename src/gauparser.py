import os
import pandas as pd
import numpy as np

from gaussian_class import *


# Costants
d       = {'30': 'Zn', '8': 'O', '1': 'H', '6': 'C', '7': 'N', '16': 'S', '12': 'Mg',
           '11': 'Na', '28': 'Ni', '20': 'Ca', '':'', '17':'Cl','35':'Br' };

d_atms  = ['Zn','C','H','O','Mg','N','Na','Ni','Ca','S','Fe','Cl','Br']



def extract_forces(filepath, outfile):
    """
        Function which takes as input the filepath of a Gaussian output file .log,
        extract the forces and write it to a file 'outfile'
        
        Inputs:
            filepath
            outfile
        
        Returns
            None
    """
    forces = False
    force  = []
    
    if os.path.isfile(filepath) == False:
        print("The file {} doesnt exists, please check it".format(filepath))
        return -1
    else:
        with open(filepath, 'r') as f:
            for line in f:
                if line.strip()[36:] == 'Forces (Hartrees/Bohr)':
                    forces = True
                    force.append('Forces (Hartrees/Bohr)  {} \n'.format(filepath[filepath.rfind('/') + 1:]))
                
                elif 'Search for a local minimum.' in line:
                    forces = False
                        
                if forces == True:
                    if len(line.strip().split()) == 5:
                        lforce, aforce = force_line(line)
                        force.append(lforce)
        
    
        with open(outfile, 'w') as f:
            for line in force:
                f.write(line)

                
def extract_energy(filepath, outfile, flag):
    """
        Function which takes as input the filepath of a Gaussian output file .log,
        extract the energy and write it to a file 'outfile'
        
        Inputs:
            filepath
            outfile
            flag      = [ONIOM, Counterpoise, SCF]
        
        Returns
            None
    """
    flags = {'ONIOM':'ONIOM: extrapolated energy',
             'Counterpoise':'complexation', 
             'SCF':'SCF Done', 
             'BSSE':'Counterpoise corrected energy'}

    if flag not in flags.keys():
        print("Error, the flag must be one of the following: {} ".format(flags.keys()))
        return -1
    
    energy  = -1
    
    if os.path.isfile(filepath) == False:
        print("The file {} doesnt exists, please check it".format(filepath))
        return -1
    else:
        with open(filepath, 'r') as f:
            for line in f:
                if flags[flag] in line:
                    if flag == 'Counterpoise':
                        energy = line.strip().split()[-3]
                    elif flag == 'BSSE':
                        energy = line.strip().split()[-1]
                    else:
                        energy = line.strip().split()[4]
    
        with open(outfile, 'w') as f:
            f.write(str(energy))
            

def extract_coord(filepath, outfile):
    """
        Function which takes as input the filepath of a Gaussian output file .log,
        extract the coordinates and write it to a file 'outfile'
        
        Inputs:
            filepath
            outfile
        
        Returns
            None
    """
  
    coord  = []
    flag_coord = False
    
    if os.path.isfile(filepath) == False:
        print("The file {} doesnt exists, please check it".format(filepath[filepath.rfind('/') + 1:]))
        return -1
    else:
        with open(filepath, 'r') as f:
            for line in f:
                if 'Symbolic Z-matrix:' in line:
                    flag_coord = True;
                    coord.append('Coordinate {}:\n'.format(filepath[filepath.rfind('/') + 1:]))
                elif 'NAtoms=' in line:
                    flag_coord = False; 
                
                if flag_coord == True:
                        if len(line.strip().split()) == 4:
                            
                            tmp_at = line.strip().split()[0]
                            for index_let, letter in enumerate(tmp_at):
                                if letter == '(':
                                    lett_atm = tmp_at[0:index_let]        
                            if lett_atm in d_atms:
                                coord.append('{:>5}{:>15}{:>15}{:>15}\n'.format(lett_atm, 
                                                                                line.strip().split()[1],
                                                                                line.strip().split()[2],
                                                                                line.strip().split()[3]))
                                
        with open(outfile, 'w') as f:
            for line in coord:
                f.write(line)
                
                
def extract_moment(filepath, outfile):
    """
        Function which takes as input the filepath of a Gaussian output file .log,
        extract the coordinates and write it to a file 'outfile'
        
        Inputs:
            filepath
            outfile
        
        Returns
            None
    """
  
    dipole, quadrupole, quadrupole_A, quadrupole_B  = [], [], [], []
    flg_dip, flg_trace, flg_quad = False, False,False
    
    if os.path.isfile(filepath) == False:
        print("The file {} doesnt exists, please check it".format(filepath[filepath.rfind('/') + 1:]))
        return -1
    else:
        with open(filepath, 'r') as f:
            for line in f:
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
                            
        with open(outfile, 'w') as f:
            f.write("Dipole {}\n".format(filepath[filepath.rfind('/') + 1:]))
            for dip in dipole:
                f.write('{:>15}{:>15}{:>15}\n'.format(dip[0],dip[1],dip[2]))    
                
            f.write("Quadrupole A {}\n".format(filepath[filepath.rfind('/') + 1:]))
            for qqa in quadrupole_A:
                f.write('{:>15}{:>15}{:>15}\n'.format(qqa[0],qqa[1],qqa[2]))
                
            f.write("Quadrupole B {}\n".format(filepath[filepath.rfind('/') + 1:]))
            for qqb in quadrupole_B:
                f.write('{:>15}{:>15}{:>15}\n'.format(qqb[0],qqb[1],qqb[2]))
            
            f.write("Quadrupole {}\n".format(filepath[filepath.rfind('/') + 1:]))
            for qqq in quadrupole:
                f.write('{:>15}{:>15}{:>15}\n'.format(qqq[0],qqq[1],qqq[2]))                       
                    
                        
def force_line(line):
    """
        extract the information from a force line of Gaussian logfile
    """
    
    if line.strip().split()[1] in d:  # Check the atomic number
        lforce = '{:>5}{:>15}{:>15}{:>15}\n'.format(
            d[line.strip().split()[1]],
            line.strip().split()[2],
            line.strip().split()[3],
            line.strip().split()[4]
        )
        
        aforce = [ d[line.strip().split()[1]],
                   float(line.strip().split()[2]),
                   float(line.strip().split()[3]),
                   float(line.strip().split()[4])]
        return lforce, aforce
    else: 
        return  " ", " "

    

    
def extract_all(filepath, name_files, ener_flag, flag_class = True):
    """
        Function which parse a Gaussian output file and extract 
        Energy, Forces, Coordinate and Moments from the log file
        
        Inputs
            filename
            name_files
            ener_flag
        
        Returns
            file-dataframe for each type
    """
    energy = -1
    force, coord, moment, c, fa                      = [], [], [], [], []
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
                        energy = line.strip().split()[-3]
                    elif ener_flag == 'BSSE':
                        energy = line.strip().split()[-1]
                    else:
                        energy = line.strip().split()[4]
                
                # Start Parsing Forces
                elif line.strip()[36:] == 'Forces (Hartrees/Bohr)':
                    flag_forces = True
                    force.append('Forces (Hartrees/Bohr)  {} \n'.format(filepath[filepath.rfind('/') + 1:]))
                    
                # Stop Parsing Forces
                elif 'Search for a local minimum.' in line:
                    flag_forces = False
                
                # Start Parsing Coordinate
                elif 'Symbolic Z-matrix:' in line:
                    flag_coord = True;
                    coord.append('Coordinate {}:\n'.format(filepath[filepath.rfind('/') + 1:]))
                    
                # Stop Parsing Coordinate
                elif 'NAtoms=' in line:
                    flag_coord = False; 
                    
                # Start Parsing Dipole Moment    
                elif 'Dipole moment (field-independent basis, Debye)' in line:
                    flg_dip = True
                
                # Start Parsing Quadrupole Moment
                elif 'Quadrupole moment (field-independent basis, Debye-Ang)' in line and 'Traceless' not in line:
                    flg_quad = True
                    
                # Start Parsing Quadrupole Moment
                elif 'Traceless Quadrupole moment (field-independent basis, Debye-Ang):' in line:
                    flg_trace = True   
                    
                if flg_trace == True and len(line.strip().split()) == 6:
                    
                    if 'XX=' in line:
                        quadrupole_A.append([line.strip().split()[1],
                                             line.strip().split()[3],
                                             line.strip().split()[5]])
                    if 'XY=' in line:
                        quadrupole_B.append([line.strip().split()[1],
                                             line.strip().split()[3],
                                             line.strip().split()[5]])
                        flg_trace = False
                
                
                if flg_dip == True and len(line.split()) == 8 :
                        dipole.append([line.strip().split()[1],
                                       line.strip().split()[3],
                                       line.strip().split()[5]])  
                        flg_dip = False
                        
                if flg_quad == True and len(line.split()) == 6 and 'XX' in line:
                    quadrupole.append([line.strip().split()[1],
                                       line.strip().split()[3],
                                       line.strip().split()[5]])
                    flg_quad = False
                    
                if flag_forces == True:
                    if len(line.strip().split()) == 5:
                        lforce, aforce = force_line(line)
                        force.append(lforce)
                        fa.append(aforce)
                
                if flag_coord == True:
                        if len(line.strip().split()) == 4:
                            tmp_at = line.strip().split()[0]
                            for index_let, letter in enumerate(tmp_at):
                                if letter == '(':
                                    lett_atm = tmp_at[0:index_let]        
                            if lett_atm in d_atms:
                                coord.append([lett_atm,
                                              float(line.strip().split()[1]),
                                              float(line.strip().split()[2]),
                                              float(line.strip().split()[3])])
                                c.append([lett_atm,
                                              float(line.strip().split()[1]),
                                              float(line.strip().split()[2]),
                                              float(line.strip().split()[3])])
        try:
            os.chdir('./results')
        except:
            os.mkdir('./results')
            os.chdir('./results')
        
        with open(name_files['moment'], 'w') as f:
            f.write("Dipole {}\n".format(filepath[filepath.rfind('/') + 1:]))
            for dip in dipole:
                f.write('{:>15}{:>15}{:>15}\n'.format(dip[0],dip[1],dip[2]))    
                
            f.write("Quadrupole A {}\n".format(filepath[filepath.rfind('/') + 1:]))
            for qqa in quadrupole_A:
                f.write('{:>15}{:>15}{:>15}\n'.format(qqa[0],qqa[1],qqa[2]))
                
            f.write("Quadrupole B {}\n".format(filepath[filepath.rfind('/') + 1:]))
            for qqb in quadrupole_B:
                f.write('{:>15}{:>15}{:>15}\n'.format(qqb[0],qqb[1],qqb[2]))
            
            f.write("Quadrupole {}\n".format(filepath[filepath.rfind('/') + 1:]))
            for qqq in quadrupole:
                f.write('{:>15}{:>15}{:>15}\n'.format(qqq[0],qqq[1],qqq[2]))
        
        with open(name_files['coordinate'], 'w') as f:
            for line in coord:
                f.write('{:>5}{:>15}{:>15}{:>15}\n'.format(line[0],line[1],line[2],line[3]))
        
        with open(name_files['energy'], 'w') as f:
            f.write(str(energy))
        
        with open(name_files['force'], 'w') as f:
            for line in force:
                f.write(line)
        
        if flag_class:
            # Moments
            df_moment = pd.DataFrame()
            df_moment['dipole-x'] = dipole[:][0]
            df_moment['dipole-y'] = dipole[:][1]
            df_moment['dipole-z'] = dipole[:][2]

            df_moment['quadru-x'] = quadrupole[:][0]
            df_moment['quadru-y'] = quadrupole[:][1]
            df_moment['quadru-z'] = quadrupole[:][2]

            df_moment['quadruA-x'] = quadrupole_A[:][0]
            df_moment['quadruA-y'] = quadrupole_A[:][1]
            df_moment['quadruA-z'] = quadrupole_A[:][2]

            df_moment['quadruB-x'] = quadrupole_B[:][0]
            df_moment['quadruB-y'] = quadrupole_B[:][1]
            df_moment['quadruB-z'] = quadrupole_B[:][2]

            # Coordinates
            c = np.array(c)
            df_coord = pd.DataFrame()
            df_coord['atom-type'] = c[:,0]
            df_coord['x']         = c[:,1]
            df_coord['y']         = c[:,2]
            df_coord['z']         = c[:,3]

            # Forces
            fa = np.array(fa[1:])
            df_force = pd.DataFrame()
            df_force['atom-type'] = fa[:,0]
            df_force['x']         = fa[:,1]
            df_force['y']         = fa[:,2]
            df_force['z']         = fa[:,3]

            # Energies
            df_energy = pd.DataFrame()
            df_energy['energy'] = np.ones(np.shape(fa[:,0]))*float(energy)


            return gau_file(df_moment, df_coord, df_force, df_energy)
    
