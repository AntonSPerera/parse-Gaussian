import os

d       = {'30': 'Zn', '8': 'O', '1': 'H', '6': 'C', '7': 'N', '16': 'S', '12': 'Mg', '11': 'Na', '28': 'Ni', '20': 'Ca', '':'', '17':'Cl','35':'Br' };
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
                        force.append(force_line(line))
        
    
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
    flags = {'ONIOM':'ONIOM: extrapolated energy', 'Counterpoise':'complexation', 'SCF':'SCF Done', 'BSSE':'Counterpoise corrected energy'}

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
                                coord.append('{:>5}{:>15}{:>15}{:>15}\n'.format(lett_atm, line.strip().split()[1],line.strip().split()[2],line.strip().split()[3]))
                                
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
        force = '{:>5}{:>15}{:>15}{:>15}\n'.format(
            d[line.strip().split()[1]],
            line.strip().split()[2],
            line.strip().split()[3],
            line.strip().split()[4]
        )
        return force
    else: 
        return  " "
    
