class gau_file():
    """
        Class gau_file --> Object derived by the parsing of a gaussian output file
        contains moment, energy, coordinate and forces
        
        Inputs:
            coordinate
            moments
            forces
            energy
        
        Returns
            self
    """
    
    def __init__(self, df_moment, df_coord, df_force, df_ener):
        self.moment = df_moment
        self.coord  = df_coord
        self.energy = df_ener
        self.force  = df_force

    