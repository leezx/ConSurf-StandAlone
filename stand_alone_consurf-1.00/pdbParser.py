
import re

class pdbParser:

    def __init__(self):

        self.SEQRES = {}
        self.ATOM = ""
        self.ATOM_withoutX = {}
        self.type = ""
        self.MODIFIED_COUNT = 0
        self.MODIFIED_LIST = ""




    def read(self, file, query_chain, DNA_AA, atom_position_filename):

        #conversion_table = {"ALA" : "A", "ARG" : "R", "ASN" : "N", "ASP" : "D", "CYS" : "C", "GLN" : "Q", "GLU" : "E", "GLY" : "G", "HIS" : "H", "ILE" : "I", "LEU" : "L", "LYS" : "K", "MET" : "M", "PHE" : "F", "PRO" : "P", "SER" : "S", "THR" : "T", "TRP" : "W", "TYR" : "Y", "VAL" : "V", "A" : "a", "T" : "t", "C" : "c", "G" : "g", "U" : "u", "I" : "i", "DA" : "a", "DT" : "t", "DC" : "c", "DG" : "g", "DU" : "u", "DI" : "i", "5CM" : "c", "5MU" : "t", "N" : "n"}
        conversion_table = {"ALA" : "A", "ARG" : "R", "ASN" : "N", "ASP" : "D", "CYS" : "C", "GLN" : "Q", "GLU" : "E", "GLY" : "G", "HIS" : "H", "ILE" : "I", "LEU" : "L", "LYS" : "K", "MET" : "M", "PHE" : "F", "PRO" : "P", "SER" : "S", "THR" : "T", "TRP" : "W", "TYR" : "Y", "VAL" : "V", "A" : "a", "T" : "t", "C" : "c", "G" : "g", "U" : "u", "I" : "i", "DA" : "a", "DT" : "t", "DC" : "c", "DG" : "g", "DU" : "u", "DI" : "i", "5CM" : "c", "N" : "n"}
        modified_residues = {"MSE" : "MET", "MLY" : "LYS", "HYP" : "PRO", "CME" : "CYS", "CGU" : "GLU", "SEP" : "SER", "KCX" : "LYS", "MLE" : "LEU", "TPO" : "THR", "CSO" : "CYS", "PTR" : "TYR", "DLE" : "LEU", "LLP" : "LYS", "DVA" : "VAL", "TYS" : "TYR", "AIB" : "ALA", "OCS" : "CYS", "NLE" : "LEU", "MVA" : "VAL", "SEC" : "CYS", "PYL" : "LYS"}
        localMODRES = {}
        FIRST = {} # first residue in chain
        fas_pos = 0
        QUERY = False
        TER = False # reached termination flag 

        if DNA_AA == "Nuc":

            UnknownChar = "N"

        else:

            UnknownChar = "X"
                        
        # open file to read MODRES
        try:

            PDBFILE = open(file, 'r')

        except:

            return 0

        try:

            MODRES_FILE = open(file + ".MODRES", 'w')

        except:

            return 0

        # read the MODRES
        line = PDBFILE.readline()
        while line != "" and not re.match(r'^ATOM', line):

            if re.match(r'^MODRES', line):

                MODRES = line[12:15].strip() # strip spaces to support NUC
                CHAIN = line[16:17]
                if CHAIN == " ":
                    
                    CHAIN = "NONE"
                    
                # we only look at the query chain
                if CHAIN != query_chain:
                    
                    line = PDBFILE.readline()
                    continue
                
                RES = line[24:27].strip() # strip spaces to support NUC

                if not MODRES in localMODRES:

                    localMODRES[MODRES] = RES
                    MODRES_FILE.write(MODRES + "\t" + RES + "\n")

                elif localMODRES[MODRES] != RES:

                    localMODRES[MODRES] = "" # two different values to the same residue

            line = PDBFILE.readline()


        MODRES_FILE.close()
        PDBFILE.close()
        
        # open position file for writing 
        try:

            CORR = open(atom_position_filename, 'w')

        except:

            return 0
        
        # reopen file to read all the file
        try:

            PDBFILE = open(file, 'r')

        except:

            return 0

        line = PDBFILE.readline()
        while line != "":

            line = line.strip()

            if re.search(r'^SEQRES', line): # SEQRES record

                chain = line[11:12] # get chain id
                if chain == " ":

                    chain = "NONE"
                                        
                # check if this chain was already processed
                if not chain in self.SEQRES:

                    self.SEQRES[chain] = ""

                # we skip the chain if it is not the query
                if query_chain != chain:
                    
                    line = PDBFILE.readline()
                    continue

                # convert to one letter format
                for acid in line[19:70].split():

                    if chain == query_chain and self.type == "":

                        if len(acid) < 3:

                            self.type = "Nuc"

                        elif len(acid) == 3:

                            if acid in modified_residues:

                                if len(modified_residues[acid]) < 3:

                                    self.type = "Nuc"

                                else:

                                    self.type = "AA"

                            else:

                                self.type = "AA"

                    # regular conversion
                    if acid in conversion_table:

                        # add to chain
                        self.SEQRES[chain] += conversion_table[acid]

                    # modified residue
                    else:

                        # count this modified residue
                        self.MODIFIED_COUNT += 1

                        # check if residue is identified
                        if acid in modified_residues and modified_residues[acid] in conversion_table:

                            self.SEQRES[chain] += conversion_table[modified_residues[acid]]

                            # add to modified residue list
                            if not acid + " > " in self.MODIFIED_LIST:

                                self.MODIFIED_LIST += acid + " > " + conversion_table[modified_residues[acid]] + "\n"

                        elif acid in localMODRES and localMODRES[acid] != "" and localMODRES[acid] in conversion_table:

                            self.SEQRES[chain] += conversion_table[localMODRES[acid]]

                            # add to modified residue list
                            if not acid + " > " in self.MODIFIED_LIST:

                                self.MODIFIED_LIST += acid + " > " + conversion_table[localMODRES[acid]] + "\n"

                        else:

                            # set residue name to X or N
                            self.SEQRES[chain] += UnknownChar

                            # add message to front of modified residue list
                            modified_changed_to_X_or_N_msg = "Modified residue(s) in this chain were converted to the one letter representation '" + UnknownChar + "'\n"

                            if not "Modified residue" in self.MODIFIED_LIST:

                                self.MODIFIED_LIST = modified_changed_to_X_or_N_msg + self.MODIFIED_LIST

            elif re.search(r'^ATOM', line):

                # extract atom data
                res = line[17:20].strip() # for DNA files there is only one or two letter code
                chain = line[21:22]
                pos = line[22:27]

                if chain == " ":

                    chain = "NONE"

                # find modified residue
                if res in modified_residues:
                    
                    mod_res = modified_residues[res]
                    
                elif res in localMODRES and localMODRES[res] != "" and localMODRES[res] in conversion_table:
                    
                    mod_res = localMODRES[res]

                else:
                    
                    mod_res = res
                  
                # convert residue to one letter
                if res in conversion_table:
                    
                    oneLetter = conversion_table[res]
                    
                else:
                    
                    oneLetter = UnknownChar
                  
                # check if we reached a new residue
                if not chain in FIRST:
                    
                    FIRST[chain] = True
                    last_pos = pos
                    new_res = True
                    self.ATOM_withoutX[chain] = oneLetter
                    
                elif pos != last_pos:
                    
                    last_pos = pos
                    new_res = True
                    self.ATOM_withoutX[chain] += oneLetter

                else:

                     new_res = False
                     
                if query_chain != chain:
                         
                    line = PDBFILE.readline()
                    continue 
                
                else:
                    
                    QUERY = True
                
                # writing atom position file
                if new_res:

                    fas_pos += 1
                    CORR.write("%s\t%d\t%s\n" %(mod_res, fas_pos, pos))
                        
                residue_number = int(line[22:26].strip())
                        
                # check type 
                if self.type == "":

                    if len(res) < 3:

                        self.type = "Nuc"

                    elif len(res) == 3:

                        if res in modified_residues:

                            if len(modified_residues[res]) < 3:

                                self.type = "Nuc"

                            else:

                                self.type = "AA"

                        else:

                            self.type = "AA"
                            
                if FIRST[chain]:
                    
                    last_residue_number = residue_number
                    FIRST[chain] = False

                elif last_residue_number < residue_number:
                    
                    while residue_number != last_residue_number + 1: # For Disorder regions
                    
                            self.ATOM += UnknownChar
                            last_residue_number += 1
                            
                if new_res:
                    
                    if mod_res in conversion_table:
                        
                        self.ATOM += conversion_table[mod_res]
                        
                    else:
                        
                        self.ATOM += UnknownChar
                                            
                last_residue_number = residue_number

            elif re.search(r'^HETATM', line):

                # extract hetatm data
                res = line[17:20].strip() # for DNA files there is only one or two letter code
                chain = line[21:22]
                pos = line[22:27]

                if chain == " ":

                    chain = "NONE"
                    
                if query_chain != chain or TER:
                         
                    line = PDBFILE.readline()
                    continue 
                
                else:
                    
                    QUERY = True
                    
                # find modified residue
                if res in modified_residues:
                    
                    mod_res = modified_residues[res]
                    
                elif res in localMODRES and localMODRES[res] != "" and localMODRES[res] in conversion_table:
                    
                    mod_res = localMODRES[res]

                else:
                    
                    mod_res = res
                  
                # convert residue to one letter
                if res in conversion_table:
                    
                    oneLetter = conversion_table[res]
                    
                else:
                    
                    oneLetter = UnknownChar
                  
                # check if we reached a new residue
                if not chain in FIRST:
                    
                    FIRST[chain] = True
                    last_pos = pos
                    new_res = True
                    self.ATOM_withoutX[chain] = oneLetter
                    
                elif pos != last_pos:
                    
                    last_pos = pos
                    new_res = True
                    self.ATOM_withoutX[chain] += oneLetter

                else:

                     new_res = False
                
                # writing atom position file
                if new_res:

                    fas_pos += 1
                    CORR.write("%s\t%d\t%s\n" %(mod_res, fas_pos, pos))
                        
                residue_number = int(line[22:26].strip())
                  
                """
                # check type 
                if self.type == "":

                    if len(res) < 3:

                        self.type = "Nuc"

                    elif len(res) == 3:

                        if res in modified_residues:

                            if len(modified_residues[res]) < 3:

                                self.type = "Nuc"

                            else:

                                self.type = "AA"

                        else:

                            self.type = "AA"
                """
                
                if FIRST[chain]:
                    
                    last_residue_number = residue_number
                    FIRST[chain] = False

                elif last_residue_number < residue_number:
                    
                    while residue_number != last_residue_number + 1: # For Disorder regions
                    
                            self.ATOM += UnknownChar
                            last_residue_number += 1
                            
                if new_res:
                    
                    if mod_res in conversion_table:
                        
                        self.ATOM += conversion_table[mod_res]
                        
                    else:
                        
                        self.ATOM += UnknownChar
                                            
                last_residue_number = residue_number   

            elif re.search(r'^TER', line):
                
                if QUERY:
                    
                    TER = True
                
            line = PDBFILE.readline()

        PDBFILE.close()
        CORR.close()
        return 1





    def get_type(self):

        return self.type

    def get_SEQRES(self):

        return self.SEQRES


    def get_ATOM(self):

        return self.ATOM

    def get_ATOM_withoutX(self):

        return self.ATOM_withoutX

    def get_MODIFIED_COUNT(self):
        

        return self.MODIFIED_COUNT
        


    def get_MODIFIED_LIST(self):

        return self.MODIFIED_LIST
























































