


import re
import os
import time
import GENERAL_CONSTANTS

bayesInterval = 3
ColorScale = {0 : 9, 1 : 8, 2 : 7, 3 : 6, 4 : 5, 5 : 4, 6 : 3, 7 : 2, 8 : 1}

def assign_colors_according_to_r4s_layers(Output, rate4site_filename, algorithm):

    # read the rate4site output into an array.
    # calculates layers according to the max and min grades from the output
    # Output : a pointer to array ; assigns each position in the MSA (array index) with a grade (arracy value)


    try:

        RATE4SITE = open(rate4site_filename, 'r')

    except:

        return("assign_colors_according_to_r4s_layers : can't open " + str(rate4site_filename))

    line = RATE4SITE.readline()
    while line != "":

        line.rstrip()

        if algorithm == "Bayes":

            # baysean
            match1 = re.match(r'^\s+(\d+)\s+(\w)\s+(\S+)\s+\[\s*(\S+),\s*(\S+)\]\s+\S+\s+(\d+)\/(\d+)', line)
            if match1:

                Output.append({'POS' : int(match1.group(1)), 'SEQ' : match1.group(2), 'GRADE' : float(match1.group(3)), 'INTERVAL_LOW' : float(match1.group(4)), 'INTERVAL_HIGH' : float(match1.group(5)), 'MSA_NUM' : int(match1.group(6)), 'MSA_DENUM' : match1.group(7)})

        else:

            # Maximum likelihood
            match2 = re.match(r'^\s*(\d+)\s+(\w)\s+(\S+)\s+(\d+)\/(\d+)', line)
            if match2:

                Output.append({'POS' : int(match2.group(1)), 'SEQ' : match2.group(2), 'GRADE' : float(match2.group(3)), 'INTERVAL_LOW' : float(match2.group(3)), 'INTERVAL_HIGH' : float(match2.group(3)), 'MSA_NUM' : int(match2.group(4)), 'MSA_DENUM' : match2.group(5)})

        line = RATE4SITE.readline()

    RATE4SITE.close()

    # the lower the number the more the conservation 
    
    max_cons = Output[0]['GRADE']
    min_cons = Output[0]['GRADE']
    for element in Output:

        if element['GRADE'] < max_cons:

            max_cons = element['GRADE']
            
        if element['GRADE'] > min_cons:

            min_cons = element['GRADE']

    # unity of conservation to be colored
    ConsColorUnityLeft = -1 * (max_cons / 4.5)
    ConsColorUnityRight = min_cons / 4.5
    """
    if max_cons >= 0:

        ConsColorUnity = max_cons
    """
    # calculates the grades for each color
    NoLayers = 10
    LeftLayers = 5
    RightLayers = 5
    ColorLayers = []
    i = 0
    while i < LeftLayers:

        ColorLayers.append(max_cons + i * ConsColorUnityLeft)
        i += 1
        
    i = 0
    while i < RightLayers:

        ColorLayers.append(ConsColorUnityRight / 2 + i * ConsColorUnityRight)
        i += 1

    # gives the color to the interval
    for element in Output:

        i = 0
        while not 'INTERVAL_LOW_COLOR' in element or not 'INTERVAL_HIGH_COLOR' in element or not 'COLOR' in element:

            if not 'INTERVAL_LOW_COLOR' in element:

                if i == NoLayers - 1:

                    element['INTERVAL_LOW_COLOR'] = 8

                elif element['INTERVAL_LOW'] >= ColorLayers[i] and element['INTERVAL_LOW'] < ColorLayers[i + 1]:

                    element['INTERVAL_LOW_COLOR'] = i

                elif element['INTERVAL_LOW'] < ColorLayers[0]:

                    element['INTERVAL_LOW_COLOR'] = 0

            if not 'INTERVAL_HIGH_COLOR' in element:

                if i == NoLayers - 1:

                    element['INTERVAL_HIGH_COLOR'] = 8

                elif element['INTERVAL_HIGH'] >= ColorLayers[i] and element['INTERVAL_HIGH'] < ColorLayers[i + 1]:

                    element['INTERVAL_HIGH_COLOR'] = i

                elif element['INTERVAL_HIGH'] < ColorLayers[0]:

                    element['INTERVAL_HIGH_COLOR'] = 0

            if not 'COLOR' in element:

                if i == NoLayers - 1:

                    element['COLOR'] = ColorScale[i - 1]

                elif element['GRADE'] >= ColorLayers[i] and element['GRADE'] < ColorLayers[i + 1]:

                    element['COLOR'] = ColorScale[i]

            i += 1

        if element['INTERVAL_HIGH_COLOR'] - element['INTERVAL_LOW_COLOR'] > bayesInterval or element['MSA_NUM'] <= 5:

            element['ISD'] = 1

        else:

            element['ISD'] = 0

    return("OK", ColorLayers)

def print_precentage_nuc(ref_nucleotide_freq, ref_position_totalNuc, out_file):

    nuc_arr = ["A", "C", "G", "T", "U", "OTHER"]

    try:

        OUT = open(out_file, 'w')

    except:

        return("print_precentage_nuc : Could not open the file " + out_file + " for writing.")

    OUT.write("\"The table details the nucleic acid variety in % for each position in the query sequence.\"\n\"Each column shows the % for that nucleic-acid, found in position ('pos') in the MSA.\"\n\"In case there are nucleotides  which are not a standard nucleic-acid in the MSA, they are represented under column 'OTHER'\"\n\n")
    OUT.write("pos")
    for nuc in nuc_arr:

        OUT.write(nuc)

    OUT.write("\n")
    # order the lines according to nucleotide position in the sequence
    for position in sorted(ref_nucleotide_freq.keys()):

        total = 0
        nuc_found = "no"
        other = ""
        other_val = 0
        max_precent = 0
        max_nuc = ""
        OUT.write(str(position))
        # the total number of nucleic acids found in the MSA for that position
        total += ref_position_totalNuc[position]
        # for each position , sort the nucleic acids variety
        nuc_in_position = sorted(ref_nucleotide_freq[position].keys())
        OUT.write("POS:%d\tNUC_IN_POS:" %position)
        for NUC_IN_POS in nuc_in_position:

            OUT.write("," + NUC_IN_POS)

        OUT.write("\n")
        j = 0
        nuc = nuc_in_position[j]
        OUT.write("NUC:" + nuc + "\n")
        val = 100 * (ref_nucleotide_freq[position][nuc] / total)
        # in order to print a table, we go by the sorted array
        for char in nuc_arr:

            # if there are non standart nuc, we calculate its total value seperately and add it under column OTHER
            i = 0
            while not nuc.upper() in "ACTGU" and j < len(nuc_in_position):

                other_val += val
                j += 1
                if j < len(nuc_in_position):

                    nuc = nuc_in_position[j]
                    val = 100 * (ref_nucleotide_freq[position][nuc] / total)

            if char == "OTHER" and other_val != 0:

                other = other_val
                OUT.write(",%.3f" %other)

            elif char == nuc.upper():

                if val > max_precent:

                    max_precent = val
                    max_nuc = nuc

                elif val == max_precent:

                    max_nuc += nuc

                OUT.write(",%.3f" %val)

                j += 1
                if j < len(nuc_in_position):

                    nuc = nuc_in_position[j]
                    val = 100 * (ref_nucleotide_freq[position][nuc] / total)

            else:

                OUT.write(",")

            OUT.write("POS\t%d\tI:%d\t%s\tVAL:%.3f\n" %(position, i, char, val))
            i += 1

        OUT.write("\n")

    OUT.close()
    return("OK")

def print_precentage(ref_residue_freq, ref_position_totalAA, out_file, ref_ConSurfGrades):

    # for each position, calculate the % for each residue which is found in the MSA in that position
    # the output is written in "Cvs" format, meaning there are ',' signs for each tab carachter

    aa_arr = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "OTHER"]

    try:

        OUT = open(out_file, 'w')

    except:

        return("print_precentage : Could not open the file " + out_file + " for writing.", 'PANIC')

    OUT.write("\"The table details the residue variety in % for each position in the query sequence.\"\n\"Each column shows the % for that amino-acid, found in position ('pos') in the MSA.\"\n\"In case there are residues which are not a standard amino-acid in the MSA, they are represented under column 'OTHER'\"\n\n")
    OUT.write("pos")
    for a in aa_arr:

        OUT.write("," + a)

    OUT.write(",MAX AA,ConSurf Grade\n")
    # order the lines according to AA position in the sequence
    for position in sorted(ref_residue_freq.keys()):

        total = 0
        aa_found = "no"
        other = ""
        other_val = 0
        max_precent = 0
        max_AA = ""
        OUT.write(str(position))
        # the total number of animo acids found in the MSA for that position
        total += ref_position_totalAA[position]
        # For each position , sort the amino acids variety
        aa_in_position = sorted(ref_residue_freq[position].keys())
        j = 0
        aa = aa_in_position[j]
        val = 100 * (ref_residue_freq[position][aa] / total)
        # in order to print a table, we go by the sorted array
        for char in aa_arr:

            # if there are non standart aa, we calculate its total value seperately and add it under column OTHER
            while not aa.upper() in "KRHDEYWSTFLIMCNQAVPG" and j < len(aa_in_position):

                other_val += val
                j += 1
                if j < len(aa_in_position):

                    aa = aa_in_position[j]
                    val = 100 * (ref_residue_freq[position][aa] / total)

            if char == "OTHER" and other_val != 0:

                other = other_val
                OUT.write(",%.3f" %other)

            elif char == aa.upper():

                if val > max_precent:

                    max_precent = val
                    max_AA = aa

                elif val == max_precent:

                    max_AA += aa

                OUT.write(",%.3f" %val)

                j += 1
                if j < len(aa_in_position):

                    aa = aa_in_position[j]
                    val = 100 * (ref_residue_freq[position][aa] / total)

            else:

                OUT.write(",")

        OUT.write(",%s %.3f" %(max_AA, max_precent))

        if position -1 < len(ref_ConSurfGrades):

            OUT.write(",%d" %ref_ConSurfGrades[position -1]['COLOR'])
            if ref_ConSurfGrades[position -1]['ISD'] == 1:

                OUT.write("*")

        OUT.write("\n")

    OUT.close()
    return(["OK"])

def create_atom_position_file(pdb_filename , atom_position_filename, chain_id, ref_to_return_hash, Type):

    tr_aa = {'LYS' : "K", 'TYR' : "Y", 'LEU' : "L", 'GLN' : "Q", 'ARG' : "R", 'TRP' : "W", 'ILE' : "I", 'ALA' : "A", 'HIS' : "H", 'SER' : "S", 'MET' : "M", 'VAL' : "V", 'ASP' : "D", 'THR' : "T", 'CYS' : "C", 'PRO' : "P", 'GLU' : "E", 'PHE' : "F", 'ASN' : "N", 'GLY' : "G"}
    modified_residues = {"MSE" : "MET", "MLY" : "LYS", "HYP" : "PRO", "CME" : "CYS", "CGU" : "GLU", "SEP" : "SER", "KCX" : "LYS", "MLE" : "LEU", "TPO" : "THR", "CSO" : "CYS", "PTR" : "TYR", "DLE" : "LEU", "LLP" : "LYS", "DVA" : "VAL", "TYS" : "TYR", "AIB" : "ALA", "OCS" : "CYS", "NLE" : "LEU", "MVA" : "VAL", "SEC" : "CYS", "PYL" : "LYS"}
    last_pos = ""
    fas = ""
    fas_pos = 0
    first = 1

    if not os.path.exists(pdb_filename) or os.path.getsize(pdb_filename) == 0:

        ref_to_return_hash['ERROR'] = "cp_rasmol_gradesPE_and_pipe.create_atom_position_file : the file " + pdb_filename + " does not exist ot is empty."
        return

    try:

        PDB = open(pdb_filename, 'r')

    except:

        ref_to_return_hash['ERROR'] = "cp_rasmol_gradesPE_and_pipe.create_atom_position_file : could not open the file " + pdb_filename + " for reading."
        return

    try:

        CORR = open(atom_position_filename, 'w')

    except:

        ref_to_return_hash['ERROR'] = "cp_rasmol_gradesPE_and_pipe.create_atom_position_file : could not open the file " + atom_position_filename + " for writing."
        return

    # going over the PDB file and extracting the chain sequence.
    line = PDB.readline()
    pdb_line = 0
    while line != "":

        pdb_line += 1
        if (re.match(r'^ATOM', line) or re.match(r'^HETATM', line)) and line[21:22] == chain_id:

            residue = line[17:20]
            if residue in modified_residues:
                
                residue = modified_residues[residue]
                
            pos = line[22:27]
            if Type == "AA":

                if residue in tr_aa:

                    if first == 1:

                        fas += tr_aa[residue]
                        last_pos = pos
                        fas_pos = 1
                        CORR.write("%s\t%d\t%s\n" %(residue, fas_pos, pos))

                    first = 0
                    if pos != last_pos:

                        fas += tr_aa[residue]
                        last_pos = pos
                        fas_pos += 1
                        CORR.write("%s\t%d\t%s\n" %(residue, fas_pos, pos))

                else:

                    if not 'WARNING' in ref_to_return_hash:

                        ref_to_return_hash['WARNING'] = "line %s : residue %s is not legal\n" %(pdb_line, residue)

                    else:

                        ref_to_return_hash['WARNING'] += "line %s : residue %s is not legal\n" %(pdb_line, residue)

            else:

                residue = residue.strip()
                if first == 1:

                    fas += residue
                    last_pos = pos
                    fas_pos += 1
                    CORR.write("%s\t%d\t%s\n" %(residue, fas_pos, pos))

                first = 0
                if pos != last_pos:

                    fas += residue
                    fas_pos += 1
                    CORR.write("%s\t%d\t%s\n" %(residue, fas_pos, pos))

                last_pos = pos

        line = PDB.readline()

    PDB.close()
    CORR.close()

    ref_to_return_hash['ATOM_AA_SEQ'] = fas
    ref_to_return_hash['AA_SEQ'] = fas

def match_seqres_pdb(pdbseq, query_seq, atom_position, chain, ref_fas2pdb, Type):

    # matches the position in the seqres/msa sequence to the position in the pdb
    """
    pdbseq = ""
    query_seq = ""

    try:

        ALN = open(seqres_atom_aln, 'r')

    except:

        return("match_seqres_pdb : Could not open the file " + seqres_atom_aln + " for reading.", "PANIC")

    line = ALN.readline()
    while line != "":

        match1 = re.match(r'^ATOM_\S+\s+(\S+)', line)
        if match1:

            pdbseq += match1.group(1)
            print("match_seqres_pdb\tATOM\t" + match1.group(1))

        else:

            match2 = re.match(r'^SEQRES_\S+\s+(\S+)', line)
            if match2:

                query_seq += match2.group(1)
                print("match_seqres_pdb\tSEQRES\t" + match2.group(1))

            else:

                match3 = re.match(r'^MSA_\S+\s+(\S+)', line)
                if match3:

                    query_seq += match3.group(1)
                    print("match_seqres_pdb\tMSA\t" + match3.group(1))

                else:

                    match4 = re.match(r'^QUERY\s+(\S+)', line)
                    if match4:

                        query_seq += match4.group(1)
                        print("match_seqres_pdb\tQUERY\t" + match4.group(1))

        line = ALN.readline()

    ALN.close()
    """
    UnKnownChar = ""
    if Type == "AA":

        UnKnownChar = "X"

    else:

        UnKnownChar = "N"

    # creating the hash that matches the position in the ATOM fasta to its position
    # in the pdb file and also the fasta ATOM position to the correct residue
    try:

        MATCH = open(atom_position, 'r')

    except:

        return ("match_seqres_pdb : Could not open the file " + atom_position + " for reading.", "PANIC")

    AtomPos_Regex = ""
    if Type == "AA":

        AtomPos_Regex = "(\w{3})\s+(\d+)\s+(\S+)"

    else:

        AtomPos_Regex = "(\w)\s+(\d+)\s+(\S+)"

    match_ATOM = {}
    res_ATOM = {}
    line = MATCH.readline()
    while line != "":

        match = re.match(AtomPos_Regex, line)
        if match:

            match_ATOM[int(match.group(2))] = match.group(1) + ":" + match.group(3) + ":" + chain

        line = MATCH.readline()

    MATCH.close()

    # creating a hash in which the key is the position in the aln (i.e the
    # position in the seqres) and the value is the correct position in the fas.
    pos_count = 0
    pos = 0
    aln_fas_seqres = {}
    for char in query_seq:

        if char != '-' and char != UnKnownChar:

            pos_count += 1
            aln_fas_seqres[pos] = pos_count

        pos += 1

    length_of_seqres = pos_count

    # creating a hash in which the key is the position in the aln (i.e the
    # position in the atoms) and the value is the correct position in the fas.
    pos_count = 0
    pos = 0
    aln_fas_atoms = {}
    for char in pdbseq:

        if char != '-' and char != UnKnownChar:

            pos_count += 1
            aln_fas_atoms[pos] = pos_count

        pos += 1

    length_of_atoms = pos_count

    pos = 0
    for char in query_seq:

        if char != '-' and char != UnKnownChar:

            fas_pos = int(aln_fas_seqres[pos])
            if pdbseq[pos] == '-' or pdbseq[pos] == UnKnownChar:

                ref_fas2pdb[fas_pos] = '-'

            else:

                ref_fas2pdb[fas_pos] = match_ATOM[aln_fas_atoms[pos]]

        pos += 1

    return("OK", length_of_seqres, length_of_atoms)

def print_rasmol(out_file, rasmol_isd, ref_colors_array, chain, proteopedia, scale = "legacy"):

    # print out new format of rasmol


    consurf_rasmol_colors = ["", "[16,200,209]", "[140,255,255]", "[215,255,255]", "[234,255,255]", "[255,255,255]", "[252,237,244]", "[250,201,222]", "[240,125,171]", "[160,37,96]", "[255,255,150]"]
    consurf_rasmol_colors_CBS = ["", "[27,120,55]", "[90,174,97]", "[166,219,160]", "[217,240,211]", "[247,247,247]", "[231,212,232]", "[194,165,207]", "[153,112,171]", "[118,42,131]", "[255,255,150]"]

    try:

        OUT = open(out_file, 'w')

    except:

        return ("print_rasmol : Could not open the file " + out_file + " for writing.", "PANIC")

    if proteopedia != "yes":

        OUT.write("select all\ncolor [200,200,200]\n\n")

    i = len(ref_colors_array) - 1
    while i > 0:

        if i == 10 and rasmol_isd == "no":

            i -= 1
            continue

        if len(ref_colors_array[i]) > 0:

            OUT.write(print_selected(ref_colors_array[i], "no"))
            OUT.write("\nselect selected and :%s\n" %chain)
            if scale == "legacy":

                OUT.write("color %s\nspacefill\n" %consurf_rasmol_colors[i])

            else:

                OUT.write("color %s\nspacefill\n" %consurf_rasmol_colors_CBS[i])

            OUT.write("define CON%d selected\n\n" %i)

        i -= 1

    OUT.close()
    return(["OK"])

def print_selected(arr_ref, print_for_pipe):

    total_text = ""
    string = ""
    if print_for_pipe == "yes":

        string = "! select "

    else:

        string = "select "

    total_length = len(string)

    if len(arr_ref) > 0:

        for aa in arr_ref:

            aa = aa.replace(":", "")
            total_length += len(aa)
            if total_length > 80:

                if re.search(r', $', string):

                    string = string[:-2]

                total_text += string + "\n"
                if print_for_pipe == "yes":

                    string = "! select selected or %s, " %aa

                else:

                    string = "select selected or %s, " %aa

                total_length = len(string)

            else:

                string += aa + ", "
                total_length += 2

    else:

        total_text += string + "none"


    if re.search(r', $', string):

        string = string[:-2]
        total_text += string

    return total_text

def create_gradesPE_ConSurf(Output, ref_match, ref_residue_freq, no_isd_residue_color, isd_residue_color, gradesPE_file, Type, ref_Solv_Acc_Pred, ALGORITHM, layers_array):


    # printing the the ConSurf gradesPE file

    seq3d_grades_isd = ""
    seq3d_grades = ""

    try:

        PE = open(gradesPE_file, 'w')

    except:

        return("create_gradesPE_ConSurf : can't open '" + gradesPE_file + "' for writing.", "PANIC")

    if Type == "AA":

        PE.write("\t Amino Acid Conservation Scores\n")

    else:

        PE.write("\t Nucleic Acid Conservation Scores\n")

    PE.write("\t=======================================\n\n")
    PE.write("The layers for assigning grades are as follows.\n")
    for i in range(1, len(layers_array)):
        
        if layers_array[i - 1] < 0:
        
            left_end = "%.3f" %layers_array[i - 1]
            
        else:
            
            left_end = " %.3f" %layers_array[i - 1]
    
            
        if layers_array[i] < 0:
        
            right_end = "%.3f" %layers_array[i]
            
        else:
            
            right_end = " %.3f" %layers_array[i]
            
        PE.write("from %s to %s the grade is %d\n" %(left_end, right_end, 10 - i))
        
    PE.write("\nIf the difference between the colors of the CONFIDENCE INTERVAL COLORS is more than 3 or the msa number (under the column titled MSA) is less than 6, there is insufficient data and an * appears in the COLOR column.\n")

    if Type == "AA":

        PE.write("\n- POS: The position of the AA in the SEQRES derived sequence.\n")

    else:

        PE.write("\n- POS: The position of the nucleotide in the SEQRES derived sequence.\n")

    PE.write("- SEQ: The SEQRES derived sequence in one letter code.\n")
    if Type == "AA":

        PE.write("- 3LATOM: The ATOM derived sequence in three letter code, including the AA's positions as they appear in the PDB file and the chain identifier.\n")

    else:

        PE.write("- ATOM: The ATOM derived sequence, including the nucleotide's positions as they appear in the PDB file and the chain identifier.\n")

    PE.write("- SCORE: The normalized conservation scores.\n")
    PE.write("- COLOR: The color scale representing the conservation scores (9 - conserved, 1 - variable).\n")
    PE.write("- CONFIDENCE INTERVAL: When using the bayesian method for calculating rates, a confidence interval is assigned to each of the inferred evolutionary conservation scores.\n")
    PE.write("- CONFIDENCE INTERVAL COLORS: When using the bayesian method for calculating rates. The color scale representing the lower and upper bounds of the confidence interval.\n")
    if Type == "AA":

        PE.write("- MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.\n")
        PE.write("- RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.\n\n")

    else:

        PE.write("- MSA DATA: The number of aligned sequences having an nucleotide (non-gapped) from the overall number of sequences at each position.\n")
        PE.write("- NUCLEOTIDE VARIETY: The nucleotide variety at each position of the multiple sequence alignment.\n\n")

    if ALGORITHM == "Bayes":
        
        CONFIDENCE_INTERVAL = "CONFIDENCE INTERVAL\tCONFIDENCE INTERVAL COLORS \t"
        
    else:
        
        CONFIDENCE_INTERVAL = ""
  
    if ref_Solv_Acc_Pred != "":
        
        PE.write(" POS\t SEQ\t    3LATOM\tSCORE\t\tCOLOR\t%sB/E    FUNCTION \tMSA DATA\t" %CONFIDENCE_INTERVAL)

    else:
    
        PE.write(" POS\t SEQ\t    ATOM\tSCORE\t\tCOLOR\t%s\tMSA DATA\t" %CONFIDENCE_INTERVAL)

    if Type == "AA":
        
        PE.write("RESIDUE VARIETY\n")
        
    else:
        
        PE.write("NUCLEOTIDE VARIETY\n")
        
    PE.write("    \t    \t        \t(normalized)\t        \t               \n")

    for elem in Output:

        pos = elem['POS']
        var = ""
        atom_3L = ref_match[pos]
        
        if ref_Solv_Acc_Pred != "" and atom_3L != "-" and atom_3L.split(":")[1] in ref_Solv_Acc_Pred:

            Solv_Acc_Pred = ref_Solv_Acc_Pred[atom_3L.split(":")[1]]

        else:

            Solv_Acc_Pred = ""
            
        PE.write("%4d" %pos)
        PE.write("\t%4s" %elem['SEQ'])
        PE.write("\t%10s" %atom_3L)
        PE.write("\t%6.3f" %elem['GRADE'])
        if elem['ISD'] == 1:

            PE.write("\t\t%3d*" %elem['COLOR'])

        else:

            PE.write("\t\t%3d" %elem['COLOR'])

        if ALGORITHM == "Bayes":
            
            PE.write("\t%6.3f, " %elem['INTERVAL_LOW'])
            PE.write("%6.3f" %elem['INTERVAL_HIGH'])
            PE.write("\t\t\t%5d," %ColorScale[elem['INTERVAL_LOW_COLOR']])
            PE.write("%1d\t\t" %ColorScale[elem['INTERVAL_HIGH_COLOR']])
        
        if Solv_Acc_Pred != "":

            PE.write("\t%3s" %Solv_Acc_Pred)

            # FUNCT/STRUCT COL
            if Solv_Acc_Pred == "e":

                if elem['COLOR'] == 9 or elem['COLOR'] == 8:

                    PE.write("\tf\t")

                else:

                    PE.write("\t\t ")

            elif elem['COLOR'] == 9:

                PE.write("\ts\t")

            else:

                PE.write("\t\t ")
                
        elif ref_Solv_Acc_Pred != "":
            
            PE.write("\t\t\t ")
            
        else:

            PE.write("\t ")

        PE.write("\t%8s" %(str(elem['MSA_NUM']) + "/" + elem['MSA_DENUM']))
        for aa in ref_residue_freq[pos]:

            var += aa + ","

        if re.search(r',$', var):

            var = var[:-1]

        PE.write("\t" + var + "\n")
        # the amino-acid in that position, must be part of the residue variety in this column
        if not re.search(elem['SEQ'], var, re.IGNORECASE):

            PE.closee()
            return("create_gradesPE_ConSurf : in position %s, the amino-acid %s does not match the residue variety: %s." %(pos, elem['SEQ'], var), "PANIC")

        # printing the residue to the rasmol script
        # assigning grades to seq3d strings
        if not '-' in atom_3L:

            atom_3L = re.search(r'(.+):', atom_3L).group(1)
            if Type == "Nuc":

                atom_3L = "D" + atom_3L

            color = elem['COLOR']
            no_isd_residue_color[color].append(atom_3L)
            if elem['ISD'] == 1:

                isd_residue_color[10].append(atom_3L)
                seq3d_grades_isd += "0"

            else:

                isd_residue_color[color].append(atom_3L)
                seq3d_grades_isd += str(color)

            seq3d_grades += str(color)

        else:

            seq3d_grades_isd += "."
            seq3d_grades += "."

    PE.write("\n\n*Below the confidence cut-off - The calculations for this site were performed on less than 6 non-gaped homologue sequences,\n")
    PE.write("or the confidence interval for the estimated score is equal to- or larger than- 4 color grades.\n")
    PE.close()

    return("OK", seq3d_grades_isd, seq3d_grades)

def create_score_file(Output, scoresFile, ref_match):

    # Scores file for NGL

    try:

        FH = open(scoresFile, 'w')

    except:

        return("create_gradesPE_ConSurf : can't open '" + scoresFile + "' for writing.", "PANIC")

    FH.write("ResInd,Score\n")

    for elem in Output:

        pos = elem['POS']
        residueStr = ref_match[pos]
        #match = re.search(r'^[a-z]{1,3}(\d+):', residueStr, re.IGNORECASE)
        if residueStr != "-":

            resIndex = residueStr.split(":")[1]
            if elem['ISD'] == 1:

                color = 0

            else:

                color = elem['COLOR']

            FH.write("%s,%d\n" %(resIndex, color))

    FH.close()



def create_part_of_pipe_new(pipe_file, unique_seqs, db, seq3d_grades_isd, seq3d_grades, length_of_seqres, length_of_atom, ref_isd_residue_color, ref_no_isd_residue_color, E_score, iterations, max_num_homol, MSAprogram, algorithm, matrix, scale = "legacy"):

    # creating part of the pipe file, which contains all the non-unique information.
    # each chain will use this file to construct the final pdb_pipe file, to be viewed with FGiJ

    if scale == "legacy":

        scale_block = "!color color_grade0 FFFF96 insufficient data yellow\n!color color_grade1 10C8D1 turquoise variable\n!color color_grade2 8CFFFF\n!color color_grade3 D7FFFF\n!color color_grade4 EAFFFF\n!color color_grade5 FFFFFF\n!color color_grade6 FCEDF4\n!color color_grade7 FAC9DE\n!color color_grade8 F07DAB\n!color color_grade9 A02560 burgundy conserved"

    else:

        scale_block = "!color color_grade0 FFFF96 insufficient data yellow\n!color color_grade1 1b7837 variable\n!color color_grade2 5aae61\n!color color_grade3 a6dba0\n!color color_grade4 d9f0d3\n!color color_grade5 f7f7f7\n!color color_grade6 e7d4e8\n!color color_grade7 c2a5cf\n!color color_grade8 9970ab\n!color color_grade9 762a83 conserved"

    # design the seq3d to be printed out to the pipe file
    seq3d_grades_isd = design_string_for_pipe(seq3d_grades_isd)
    seq3d_grades = design_string_for_pipe(seq3d_grades)

    # creating the frequencies array which corresponds the number of residues in each grade
    [consurf_grade_freqs_isd, consurf_grade_freqs] = freq_array(ref_isd_residue_color, ref_no_isd_residue_color)

    # Taking Care of Strings
    if max_num_homol == "all":

        max_num_homol = "\"all\""

    # write to the pipe file
    try:

        PIPE = open(pipe_file, 'w')

    except:

        return("cannot open the file " + pipe_file + " for writing.", 'PANIC')

    PIPE.write("""! consurf_psi_blast_e_value = %s;
! consurf_psi_blast_database = "%s";
! consurf_psi_blast_iterations = %s;
! consurf_max_seqs = %s;
! consurf_apd = 0.99;
! consurf_alignment = "%s";
! consurf_method = "%s";
! consurf_substitution_model =  "%s";
!
! consurf_seqres_length = %s;
! consurf_atom_seq_length = %s;
! consurf_unique_seqs = %s;
! consurf_grade_freqs_isd = %s;
! consurf_grade_freqs = %s;
!
! seq3d_grades_isd =
%s
!
! seq3d_grades =
%s
!
!
!! ====== CONTROL PANEL OPTIONS SECTION ======
!js.init
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! pipe_title_enlarged = false;
! pipe_background_color = "white";
!
!! Specify the custom consurf control panel
!!
! pipe_cp1 = "consurf/consurf.htm";
!
!! If you want the frontispiece to be reset every time you enter this
!! page, use false. If this is a one-page presentation (no contents)
!! and you want to be able to return from QuickViews without resetting
!! the view, use true.
!!
! frontispiece_conditional_on_return = true;
!
!! Open the command input slot/message box to 30%% of window height.
!!
! pipe_show_commands = true;
! pipe_show_commands_pct = 30;
!
!! Don't show the PiPE presentation controls in the lower left frame.
!!
! pipe_hide_controls = true;
!
!! Hide development viewing mode links at the bottom of the control panel.
!!
! pipe_tech_info = false;
!
!! pipe_start_spinning = true; // default is PE's Preference setting.
!! top.nonStopSpin = true; // default: spinning stops after 3 min.
!!
!! ====== COLORS SECTION ======
!!
!color color_carbon C8C8C8
!color color_sulfur FFC832
!
!! Ten ConSurf color grades follow:
!!
%s
!
!
!! ====== SCRIPTS SECTION ======
!!----------------------------------------
!!
!spt #name=select_and_chain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!!----------------------------------------
!!
!spt #name=view01
! @spt consurf_view_isd
!
!!----------------------------------------
!!
!spt #name=hide_all
! restrict none
! ssbonds off
! hbonds off
! dots off
! list * delete
!
!!----------------------------------------
!! common_spt uses CPK carbon gray (or phosphorus yellow) for backbones.
!!
!spt #name=common_spt
! @spt hide_all
! select all
! color [xC8C8C8] # rasmol/chime carbon gray
! select nucleic
! color [xFFA500] # phosphorus orange
! select hetero
! color cpk
! select not hetero
! backbone 0.4
! javascript top.water=0
!
! ssbonds 0.3
! set ssbonds backbone
! color ssbonds @color_sulfur
!
! select hetero and not water
! spacefill 0.45
! wireframe 0.15
! dots 50
!
! select protein
! center selected
!
!!----------------------------------------
!!
!spt #name=consurf_view_isd
! @spt common_spt
! @for $=0, 9
! @spt select_isd_grade$
! @spt select_and_chain
! color @color_grade$
! spacefill
! @endfor
! zoom 115
!
!!----------------------------------------
""" %(E_score, db, iterations, max_num_homol, MSAprogram, algorithm, matrix, length_of_seqres, length_of_atom, unique_seqs, consurf_grade_freqs_isd, consurf_grade_freqs, seq3d_grades_isd, seq3d_grades, scale_block))

    lineToPrint = ""
    i = 9
    while i > 0:

        PIPE.write("!!\n!spt #name=select_isd_grade%d\n!\n" %i)
        lineToPrint = print_selected(ref_isd_residue_color[i], "yes")
        if re.search(r'select', lineToPrint):

            PIPE.write(lineToPrint + "\n")

        PIPE.write("!\n!\n!!----------------------------------------\n")
        i -= 1

    PIPE.write("!!\n!spt #name=select_isd_grade0\n")
    lineToPrint = print_selected(ref_isd_residue_color[10], "yes")
    if re.search(r'select', lineToPrint):

        PIPE.write(lineToPrint + "\n")

    PIPE.write("!\n!\n!!----------------------------------------\n")

    i = 9
    while i > 0:

        PIPE.write("!!\n!spt #name=select_grade%d\n!\n" %i)
        lineToPrint = print_selected(ref_no_isd_residue_color[i], "yes")
        if re.search(r'select', lineToPrint):

            PIPE.write(lineToPrint + "\n")

        PIPE.write("!\n!\n!!----------------------------------------\n")
        i -= 1

    PIPE.write("!!\n!spt #name=select_grade0\n! select none\n!!\n")
    PIPE.write("!! ====== END OF CONSURF PiPE BLOCK ======\n")
    PIPE.close()
    return(["OK"])

def extract_data_from_pdb(input_pdb_file):

    header = ""
    title = ""
    compnd = ""

    try:

        PDB = open(input_pdb_file, 'r')

    except:

        return("extract_data_from_pdb : Could not open the file " + input_pdb_file + " for reading.", "PANIC")

    line = PDB.readline()
    while line != "":

        match1 = re.match(r'^HEADER', line)
        if match1:

            header = line.rstrip()

        else:

            match2 =re.match(r'^TITLE\s+\d*\s(.*)', line)
            if match2:

                title += match2.group(1) + " "

            else:

                match3 = re.match(r'^COMPND\s+\d*\s(.*)', line)
                if match3:

                    compnd += match3.group(1) + " "

                elif re.match(r'^SOURCE', line) or re.match(r'^KEYWDS', line) or re.match(r'^AUTHOR', line) or re.match(r'^SEQRES', line) or re.match(r'^ATOM', line):

                    break # no nead to go over all the pdb

        line = PDB.readline()

    PDB.close()
    return("OK", header, title, compnd)

def design_string_for_pipe(string_to_format):

    # take a string aaaaaaa and returns it in this format: ! "aaa" +\n! "aa";\n

    part = string_to_format
    newPart = ""

    while len(part) > 73:

        newPart += "! \"" + part[:73] + "\" +\n"
        part = part[73:]

    newPart += "! \"" + part + "\" ;"

    return newPart

def freq_array(isd_residue_color, no_isd_residue_color):

    # design the frequencies array

    consurf_grade_freqs_isd = "Array(" + str(len(isd_residue_color[10]))
    i = 1
    while i < 10:

        consurf_grade_freqs_isd += "," + str(len(isd_residue_color[i]))
        i += 1

    consurf_grade_freqs_isd += ")"

    consurf_grade_freqs = "Array(0"
    i = 1
    while i < 10:

        consurf_grade_freqs += "," + str(len(no_isd_residue_color[i]))
        i += 1

    consurf_grade_freqs += ")"

    return(consurf_grade_freqs_isd, consurf_grade_freqs)

def create_consurf_pipe_new(results_dir, IN_pdb_id_capital, chain, ref_header_title, final_pipe_file, identical_chains, partOfPipe, current_dir, run_number, msa_filename, query_name_in_msa = "", tree_filename = "", submission_time = "", completion_time = "", run_date = ""):


    # Create the pipe file for FGiJ

    if chain == 'NONE':

        chain = ""
        identical_chains = ""

    pdb_dir = results_dir + IN_pdb_id_capital + "/"

    # read info from the pdb file
    ans = ref_header_title

    header_line = ""
    if len(ans) > 1:

        header_line = ans[1]

    title_line = ""
    if len(ans) > 2:

        title_line = ans[2]

    elif len(ans) > 3:

        title_line = ans[3]

    if title_line == "":

        title_line = "! \"No title or compound description was found in the PDB file\";"

    else:

        title_line = design_string_with_spaces_for_pipe(title_line)

    # design the identical chains line
    identical_chains_line = "! consurf_identical_chains = \"%s\";" %identical_chains.replace(" ", "")

    current_dir += "_sourcedir"
    # in case there is a source dir - we determine the var query_name_in_msa
    if os.path.exists(current_dir):

        try:

            SOURCEDIR = open(current_dir, 'r')

        except:

            return("create_consurf_pipe : cannot open " + current_dir + " for reading.")

        match = re.match(r'(\d[\d\w]{3})\/(\w)', SOURCEDIR.readline())
        if match:

            query_name_in_msa = SOURCEDIR.group(1) + SOURCEDIR.group(2)
            SOURCEDIR.close()

    if query_name_in_msa == "":

        query_name_in_msa = IN_pdb_id_capital + chain.upper()

    # write to the pipe file
    try:

        PIPE_PART = open(partOfPipe, 'r')

    except:

        return("create_consurf_pipe : cannot open " + partOfPipe + " for reading.")

    try:

        PIPE = open(final_pipe_file, 'w')

    except:

        return("create_consurf_pipe : cannot open " + final_pipe_file + " for writing.")

    if header_line != "":

        PIPE.write(header_line + "\n")

    else:

        PIPE.write("HEADER                                 [THIS LINE ADDED FOR JMOL COMPATIBILITY]\n")

    PIPE.write("""!! ====== IDENTIFICATION SECTION ======
!js.init
! consurf_server = "consurf";
! consurf_version = "3.0";
! consurf_run_number = \"%s\";
! consurf_run_date = \"%s\";
! consurf_run_submission_time = \"%s\";
! consurf_run_completion_time = \"%s\";
!
! consurf_pdb_id = \"%s\";
! consurf_chain = \"%s\";
%s
! consurf_msa_filename = \"%s\";
! consurf_msa_query_seq_name = \"%s\";
! consurf_tree_filename = \"%s\";
!
""" %(run_number, run_date, submission_time, completion_time, IN_pdb_id_capital, chain, identical_chains_line, msa_filename, query_name_in_msa, tree_filename))

    titleFlag = 0
    line = PIPE_PART.readline()
    while line != "":

        if re.match(r'^~~~+', line):

            if titleFlag == 0:

                PIPE.write("! pipe_title = \"<i>ConSurf View:</i> %s chain %s.\"\n!! pipe_subtitle is from TITLE else COMPND\n!!\n! pipe_subtitle =\n%s\n" %(IN_pdb_id_capital, chain, title_line))
                titleFlag = 1

            elif chain != "":

                PIPE.write("! select selected and :%s\n" %chain)

            else:

                PIPE.write("! select selected and protein\n")

        else:

            PIPE.write(line)

        line = PIPE_PART.readline()

    PIPE_PART.close()
    PIPE.close()

    return("OK")

def design_string_with_spaces_for_pipe(part_input):

    words = part_input.split()
    newPart = "! \"" +words[0]
    part = ""
    for word in words[1:]:

        # if adding another word to the string will yeild a too long string - we cut it.
        if len(word) + 1 + len(newPart) > 76:

            part += newPart + " \" +\n"
            newPart = "! \"" + word

        else:

            newPart += " " + word

    part += newPart + "\" ;"

    return part

def add_pdb_data_to_pipe(pdb_file, pipe_file):

    # create the file to be shown using FGiJ. read the pdb file and concat header pipe to it.

    try:

        PIPE = open(pipe_file, 'a')

    except:

        return("add_pdb_data_to_pipe: cannot open " + pipe_file + " for writing.")

    try:

        PDB_FILE = open(pdb_file, 'r')

    except:

        return("add_pdb_data_to_pipe: cannot open the " + pdb_file + " for reading.")

    line = PDB_FILE.readline()
    while line != "":

        if not re.match(r'^HEADER', line):

            PIPE.write(line)

        line = PDB_FILE.readline()

    PIPE.close()
    PDB_FILE.close()

    return("OK")


def ReplaceTempFactConSurfScore(query_chain, input, PE, out1, out2):

    # Creates The ATOM section with ConSurf grades instead of the TempFactor column, creates PDB file with ConSurf grades

    grades = {}
    gradesPE_info = {}
    gradesPE_info_with0 = {}

    insufficient = read_ConSurf_gradesPE(PE, gradesPE_info, gradesPE_info_with0)

    if insufficient != "yes" and insufficient != "no":

        # error in read_ConSurf_gradesPE
        return(insufficient)

    try:

        PDB = open(input, 'r')

    except:

        return("could not open file '" + input + "' for reading.\n")

    try:

        OUTP = open(out1, 'w')

    except:

        return("could not open the file '" + out1 + "' for writing.\n")

    if insufficient == "yes":

        try:

            OUT_ISD = open(out2, 'w')

        except:

            return("could not open the file '" + out2 + "' for writing.\n")

    line = PDB.readline()
    while line != "":

        if re.match(r'^\s*$', line):

            # line is empty
            line = PDB.readline()
            continue

        if re.match(r'^ATOM', line):

            # delete new line character
            line = line[:-1]

            # if ATOM rec too short, right pad it with spaces till it is long enough to contain occupancy & tempFactor
            while len(line) < 72:

                line += " "

            tempFactor = line[60:72]
            tempFactor_isd = tempFactor
            chain = line[21:22]
            if chain == query_chain or (chain == " " and re.match(r'none', query_chain, re.IGNORECASE)):

                residue = (line[22:27]).strip()
                if residue in gradesPE_info:

                    # the TF is updated with the number from gradesPE
                    tempFactor = "     " + gradesPE_info[residue] + "      "

                if insufficient == "yes" and residue in gradesPE_info:

                    tempFactor_isd = "     " + gradesPE_info_with0[residue] + "      "

            else:

                tempFactor = "            "
                tempFactor_isd = "            "

            OUTP.write(line[:60] + tempFactor + "\n")
            if insufficient == "yes":

                OUT_ISD.write(line[:60] + tempFactor_isd + "\n")

        else:

            OUTP.write(line)
            if insufficient == "yes":

                OUT_ISD.write(line)

        line = PDB.readline()

    OUTP.close()
    if insufficient == "yes":

        OUT_ISD.close()

    return("OK", insufficient)

def read_ConSurf_gradesPE(gradesPE_file, gradesPE_hash_ref, gradesPE_0_hash_ref):

    # the routine matches each position in the gradesPE file its ConSurf grade. In case there was a grade mark with *, we put it in a seperate hash with the grade 0.
    # the routine returns "yes" if a * was found and "no" otherwise

    insufficient = "no"

    try:

        GRADES = open(gradesPE_file, 'r')

    except:

        return("cp_rasmol_gradesPE_and_pipe.read_ConSurf_gradesPE can't open '" + gradesPE_file + "' for reading.")

    line = GRADES.readline()
    while line != "":

        if re.match(r'^\s*\d+\s+\w', line):

            grades = line.split()
            if grades[2] != "-":
            
                 grades[2] = ((grades[2]).split(':'))[1]
                 
            #grades[2] = re.sub(r'[a-z]', "", grades[2], flags=re.IGNORECASE)
            if re.match(r'\d\*?', grades[4]):

                # if it is insufficient color - we change its grade to 0, which will be read as light yellow
                if re.match(r'(\d)\*', grades[4]):

                    gradesPE_hash_ref[grades[2]] = grades[4]
                    gradesPE_0_hash_ref[grades[2]] = "0"
                    insufficient = "yes"

                else:

                    gradesPE_hash_ref[grades[2]] = grades[4]
                    gradesPE_0_hash_ref[grades[2]] = grades[4]

        line = GRADES.readline()

    GRADES.close()

    return(insufficient)

def read_Rate4Site_gradesPE(gradesPE_file, gradesPE_hash_ref):

    # the routine matches each position in the gradesPE file its Rate4Site grade.

    try:

        GRADES = open(gradesPE_file, 'r')

    except:

        return("cp_rasmol_gradesPE_and_pipe.read_Rate4Site_gradesPE: Can't open '" + gradesPE_file + "' for reading.")

    line = GRADES.readline()
    while line != "":

        if re.match(r'^\s*\d+\s+\w', line):

            grades = line.split()
            #ResNum = re.sub(r'[a-z]', "", grades[2], flags=re.IGNORECASE)
            if grades[2] != "-":
            
                 grades[2] = (grades[2]).split(":")[1]
                 
            gradesPE_hash_ref[grades[2]] = grades[3]

        line = GRADES.readline()

    GRADES.close()

    return("OK")

def replace_tempFactor(PdbFile, chain, HashRef, Out):

    # replace the tempFactor column in the PDB file

    # read the PDB file to an array
    try:

        READ_PDB = open(PdbFile, 'r')

    except:

        return("cp_rasmol_gradesPE_and_pipe.replace_tempFactor: Can't open '" + PdbFile + "' for reading.")

    # write the PDB file and replace the tempFactor column
    # with the new one.
    try:

        WRITE_PDB = open(Out, 'w')

    except:

        return("cp_rasmol_gradesPE_and_pipe.replace_tempFactor: Can't open '" + Out + "' for writing.")

    line = READ_PDB.readline()
    while line != "":

        if re.match(r'^ATOM', line):

            PDBchain = line[21:22]
            ResNum = (line[22:27]).strip()

            if PDBchain == chain and ResNum in HashRef:

                while len(HashRef[ResNum]) < 6:

                    HashRef[ResNum] = " " + HashRef[ResNum]

                WRITE_PDB.write(line[:60] + HashRef[ResNum] + line[66:])

            else:

                WRITE_PDB.write(line[:60] + "      " + line[66:])

        else:

            WRITE_PDB.write(line)

        line = READ_PDB.readline()

    READ_PDB.close()
    WRITE_PDB.close()

    return("OK")

def create_rasmol_page(fileName, identical_chains, run_number, working_dir, chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, chimerax_script_for_figure, chimerax_script_for_figure_isd, rasmolFILE, rasmol_isdFILE, pdb_file_name, Used_PDB_Name, rasmolFILE_CBS, rasmol_isdFILE_CBS):

    # Create the html describing how to create rasmol Hige resolution images

    server_html_dir = GENERAL_CONSTANTS.CONSURF_HTML_DIR
    consurf2016_url = GENERAL_CONSTANTS.CONSURF_URL + "2016/"

    try:

        RASMOL_HTML = open(fileName, 'w')

    except:

        return('sys_error', "cp_rasmol_gradesPE_and_pipe.create_rasmol_page: Can't open the file " + fileName + " for writing.")

    RASMOL_HTML.write("""<?php
\$TITLE_LINE = "ConSurf Image using RasMol";

include (\"%stemplates/definitions.tpl\");
include (\"%stemplates/header.tpl\");
?>
    <center><h2>Create a high resolution RasMol figure</h2></center>""" %(server_html_dir, server_html_dir))

    RASMOL_HTML.write("<ol>")
    RASMOL_HTML.write("<li>Download the <a href = \"%s\" download=\"%s.pdb\">PDB file</a></li>" %(pdb_file_name, Used_PDB_Name))
    RASMOL_HTML.write("<li>Download the RasMol Coloring Script:</li>")
    RASMOL_HTML.write("(a) <A HREF='%s' download> Showing Insufficient Data</A>" %rasmol_isdFILE)
    RASMOL_HTML.write(" [or the <A HREF='%s'>color-blind friendly version</A>]<br>\n" %rasmol_isdFILE_CBS)
    RASMOL_HTML.write("(b) <A HREF='%s' download> Hiding Insufficient Data</A>" %rasmolFILE)
    RASMOL_HTML.write(" [or the <A HREF='%s'>color-blind friendly version</A>]<br>\n" %rasmolFILE_CBS)
    RASMOL_HTML.write("<li>Start RasMol</li>")
    RASMOL_HTML.write("<li>Load the PDB file (File → Open → select \'Protein Databank\' → select the PDB file)</li>")
    RASMOL_HTML.write("<li>Load the RasMol Coloring Script(File → Open →  select \'RasMol Script\' → select the script) </li>")
    RASMOL_HTML.write("</ol>")
    RASMOL_HTML.write("<br>")

    if identical_chains != chain and chain != "NONE":

        RASMOL_HTML.write("<form action=\"%sproject_on_identical_chains.php\" method=\"post\">" %consurf2016_url)
        RASMOL_HTML.write("<br>Act on identical chains: ")

        for identical_chain in identical_chains.split():

            if identical_chain == chain:

                RASMOL_HTML.write("&nbsp<input type=checkbox name=\"identicalChainToProjectOn[]\" value=%s id=\"identicalChainToProjectOn%s\" checked=\"checked\">%s  " %(identical_chain, identical_chain, identical_chain))

            else:

                RASMOL_HTML.write("&nbsp<input type=checkbox name=\"identicalChainToProjectOn[]\" value=%s id=\"identicalChainToProjectOn%s\">%s  " %(identical_chain, identical_chain, identical_chain))

        # hidden params
        RASMOL_HTML.write("<input type=\"hidden\" name=\"run_number\" value=%s id=\"run_number\">" %run_number)
        RASMOL_HTML.write("<input type=\"hidden\" name=\"working_dir\" value=%s id=\"working_dir\">" %working_dir)
        RASMOL_HTML.write("<input type=\"hidden\" name=\"chain\" value=%s id=\"chain\">" %chain)
        RASMOL_HTML.write("<input type=\"hidden\" name=\"ATOMS_with_ConSurf_Scores\" value=%s id=\"ATOMS_with_ConSurf_Scores\">" %ATOMS_with_ConSurf_Scores)
        RASMOL_HTML.write("<input type=\"hidden\" name=\"ATOMS_with_ConSurf_Scores_isd\" value=%s id=\"ATOMS_with_ConSurf_Scores_isd\">" %ATOMS_with_ConSurf_Scores_isd)
        RASMOL_HTML.write("<input type=\"hidden\" name=\"chimerax_script_for_figure\" value=%s id=\"chimerax_script_for_figure\">" %chimerax_script_for_figure)
        RASMOL_HTML.write("<input type=\"hidden\" name=\"chimerax_script_for_figure_isd\" value=%s id=\"chimerax_script_for_figure_isd\">" %chimerax_script_for_figure_isd)
        RASMOL_HTML.write("<input type=\"hidden\" name=\"pdbFileName\" value=%s id=\"pdbFileName\">" %pdb_file_name)
        RASMOL_HTML.write("&nbsp&nbsp&nbsp")
        RASMOL_HTML.write("<input type=\"submit\" value=\"Calculate & Download\">")
        RASMOL_HTML.write("</form>")

    RASMOL_HTML.write("<br>")
    RASMOL_HTML.write("""<?php
    include(\"%stemplates/footer.tpl\");
?>
</body></html>""" %server_html_dir)

    RASMOL_HTML.close()

    return("OK")

def color_with_chimera(query_name, msa_file, rate4site_filename, chimera_scf, headers, chimera_scf_isd, headers_isd):

    # For chimera scripts, in order to color the positions in the MSA only according to query sequence and ignore positions where there are gaps, the MSA and the R4S files should be read simultaneously.
    # the rorutine will only create 'insufficient data' scripts, if insufficient data was calculated from r4s grades.
    # input: MSA - CLUSTALW format only!!

    r4s_position_grade = []
    aln_position_grade = {}
    aln_position_grade_isd = {}
    isd = False

    ans = assign_colors_according_to_r4s_layers(r4s_position_grade, rate4site_filename)
    if ans != "OK":

        return(ans)

    # open MSA file, read only the query sequence lines. assign each position with increasing number. if it is not a '-' : read the grade from the array r4s_position_grade and put it in the hash aln_position_grade.
    # if there is insufficient data: save it to a saparate hash

    position_in_msa = 0
    position_in_r4s = 0
    FOUND_QUERY = False
    r4s_position_grade_ended = False # True if we reached the end of the array r4s_position_grade

    try:

        MSA = open(msa_file, 'r')

    except:

        return("can't' open the file " + msa_file + " for reading.")

    line = MSA.readline()
    while line != "" and not r4s_position_grade_ended:

        line = line.strip()
        match = re.match(r'^(\S+)\s+([A-Za-z-]+)$', line)
        if match:

            seq_name = match.group(1)
            seq = match.group(2)
            if query_name != seq_name:

                line = MSA.readline()
                continue

            # we are in the query name
            FOUND_QUERY = True
            for char in seq:

                if position_in_r4s >= len(r4s_position_grade):

                    r4s_position_grade_ended = True
                    break

                if char != '-':

                    if r4s_position_grade[position_in_r4s]['ISD'] == 1:

                        # we reached the end of the array r4s_position_grade
                        aln_position_grade_isd[position_in_msa] = 0
                        isd = True

                    else:

                        aln_position_grade_isd[position_in_msa] = r4s_position_grade[position_in_r4s]['COLOR']

                aln_position_grade[position_in_msa] = r4s_position_grade[position_in_r4s]['COLOR']
                position_in_r4s += 1

            position_in_msa += 1

        line = MSA.readline()

    MSA.close()

    if not FOUND_QUERY:

        # DONT FIND QUERY EXACT NAME THUS LOOP AGAIN AND LOOK FOR LINE CONTAIN THE QUERY NAME (FOR TRUNCATED NAMES)

        try:

            MSA = open(msa_file, 'r')

        except:

            return("can't' open the file " + msa_file + " for reading.")

        line = MSA.readline()
        while line != "" and not r4s_position_grade_ended:

            line = line.strip()
            match = re.match(r'^(\S+)\s+([A-Za-z-]+)$', line)
            if match:

                seq_name = match.group(1)
                seq = match.group(2)
                if query_name in seq_name:

                    line = MSA.readline()
                    continue

                # we are in the query name
                for char in seq:

                    if position_in_r4s >= len(r4s_position_grade):

                        # we reached the end of the array r4s_position_grade
                        r4s_position_grade_ended = True
                        break

                    if char != '-':

                        if r4s_position_grade[position_in_r4s]['ISD'] == 1:

                            aln_position_grade_isd[position_in_msa] = 0
                            isd = True

                        else: aln_position_grade_isd[position_in_msa] = r4s_position_grade[position_in_r4s]['COLOR']

                    aln_position_grade[position_in_msa] = r4s_position_grade[position_in_r4s]['COLOR']
                    position_in_r4s += 1

                position_in_msa += 1

            line = MSA.readline()

        MSA.close()

    print_chimera_scf_and_histogram(chimera_scf, headers, aln_position_grade, position_in_msa)
    if isd:

        print_chimera_scf_and_histogram(chimera_scf_isd, headers_isd, aln_position_grade_isd, position_in_msa)

    return("OK", isd)

def print_chimera_scf_and_histogram(chimera_scf, chimera_headers, ref_aln_position_grade, position_in_msa):

    consurf_colors_chimera = [" 255 255 150", " 16 200 209", " 140 255 255", " 215 255 255", " 234 255 255", " 255 255 255", " 252 237 244", " 250 201 222", " 240 125 171", " 160 37 96"]

    try:

        SCF = open(chimera_scf, 'w')

    except:

        return("Can't open the file " + chimera_scf + " for writing.")

    try:

        HEADER = open(chimera_headers, 'w')

    except:

        return("Can't open the file " + chimera_headers + " for writing.")

    HEADER.write("name: Conservation scores\nstyle: character\n")

    i = 0
    while i < position_in_msa:

        if i in ref_aln_position_grade:

            SCF.write("%d %d 0 0 %s\n" %(i, i, consurf_colors_chimera[ref_aln_position_grade[i]]))
            HEADER.write("\t%d\t%s\tblack\n" %(i + 1, consurf_colors_chimera[ref_aln_position_grade[i]]))

        i += 1

    SCF.close()

    HEADER.write("#  ConSurf data shown as histogram\nname: ConSurf histogram\nstyle: numeric\n")

    i = 0
    while i < position_in_msa:

        if i in ref_aln_position_grade:

            histo = ""
            if ref_aln_position_grade[i] == 9:

                histo = "1.0"

            else:

                histo = "." + str(ref_aln_position_grade[i] + 1)

            HEADER.write("\t%d\t%s\tblack\n" %(i + 1, histo))

        i += 1

    HEADER.close()

def create_chimera_image_script(script_file, pdb_file_name, pdb_path):

    try:

        CHIMERAX = open(script_file, 'w')

    except:

        return("Failed to open " + script_file + " for writing.")

    CHIMERAX.write("""<?xml version="1.0"?>
  <ChimeraPuppet type="std_webdata">
<web_files>
<file name="%s" format="text" loc="%s"/>
</web_files>
<commands>
  <mid_cmd>colordef CONS10 1.00 1.00 0.59</mid_cmd>
  <mid_cmd>colordef CONS9 0.63 0.15 0.38</mid_cmd>
  <mid_cmd>colordef CONS8 0.94 0.49 0.67</mid_cmd>
  <mid_cmd>colordef CONS7 0.98 0.79 0.87</mid_cmd>
  <mid_cmd>colordef CONS6 0.99 0.93 0.96</mid_cmd>
  <mid_cmd>colordef CONS5 1.00 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS4 0.92 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS3 0.84 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS2 0.55 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS1 0.06 0.78 0.82</mid_cmd>
  <mid_cmd>color CONS10 @/bfactor=0</mid_cmd>
  <mid_cmd>color CONS9 @/bfactor=9 </mid_cmd>
  <mid_cmd>color CONS8 @/bfactor=8</mid_cmd>
  <mid_cmd>color CONS7 @/bfactor=7</mid_cmd>
  <mid_cmd>color CONS6 @/bfactor=6</mid_cmd>
  <mid_cmd>color CONS5 @/bfactor=5</mid_cmd>
  <mid_cmd>color CONS4 @/bfactor=4</mid_cmd>
  <mid_cmd>color CONS3 @/bfactor=3</mid_cmd>
  <mid_cmd>color CONS2 @/bfactor=2</mid_cmd>
  <mid_cmd>color CONS1 @/bfactor=1</mid_cmd>
  <mid_cmd>preset apply pub 3;repr cpk;show;focus;color red ligand</mid_cmd>
</commands>
</ChimeraPuppet>""" %(pdb_file_name, pdb_path + pdb_file_name))

    CHIMERAX.close()
    os.chmod(script_file, 0o755)

    return("OK")

def create_chimera_script_align_tree(www_dir, chimerax_file, msa_file, tree_file, scf_file, hdr_file):

    try:

        CHIMERAX = open(chimerax_file, 'w')

    except:

        return("Failed to open " + chimerax_file + " for writing.")

    CHIMERAX.write("""<?xml version="1.0"?>
  <ChimeraPuppet type="std_webdata">
<web_files>""")

    if msa_file != "":

        CHIMERAX.write("""<file  name=\"%s\" format=\"text\" loc=\"%s\"/>
</web_files>
<commands>

<!-- the following 3 lines locate the Multalign Viewer instance
    that was created by opening the alignment file, and stores a reference
    to the instance as the variable 'mav' -->
  <py_cmd>from MultAlignViewer.MAViewer import MAViewer</py_cmd>
  <py_cmd>from chimera.extension import manager</py_cmd>
  <py_cmd>mav = [inst for inst in manager.instances if isinstance(inst, MAViewer)][-1]</py_cmd>
<!-- read in/show Consurf tree -->

<!-- hide initial headers, show phylogeny tree, load the coloring file,
    and make residue letters black.
    This uses two possible sets of calls:  one for the 1.2540 release
    and one for later versions that uses a better API.
    The 'if' condition guarantees that the code will work no
    matter what version the user has -->
<py_cmd>
if hasattr(mav, 'loadScfFile'):
    mav.hideHeaders(mav.headers(shownOnly=True))
    mav.usePhylogenyFile("%s", askReorder=False)
    mav.loadScfFile("%s")
    mav.useColoringFile(None)
else:
    mav.hideHeaders(mav.seqCanvas.headerDisplayOrder())
    mav.usePhylogenyFile("%s")
    mav.regionBrowser.loadScfFile("%s")
    from MultAlignViewer.prefs import RC_BLACK
    mav.seqCanvas.setColorFunc(RC_BLACK)
</py_cmd>
<!-- read in/show Consurf headers -->
    <py_cmd>mav.readHeaderFile(\"%s\")</py_cmd>""" %(msa_file, www_dir + msa_file, www_dir + tree_file, www_dir + scf_file, www_dir + tree_file, www_dir + scf_file, www_dir + hdr_file))

    else:

        CHIMERAX.write("""</web_files>
<commands>
  <mid_cmd>colordef CONS10 1.00 1.00 0.59</mid_cmd>
  <mid_cmd>colordef CONS9 0.63 0.15 0.38</mid_cmd>
  <mid_cmd>colordef CONS8 0.94 0.49 0.67</mid_cmd>
  <mid_cmd>colordef CONS7 0.98 0.79 0.87</mid_cmd>
  <mid_cmd>colordef CONS6 0.99 0.93 0.96</mid_cmd>
  <mid_cmd>colordef CONS5 1.00 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS4 0.92 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS3 0.84 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS2 0.55 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS1 0.06 0.78 0.82</mid_cmd>
  <mid_cmd>color CONS10 @/bfactor=0</mid_cmd>
  <mid_cmd>color CONS9 @/bfactor=9 </mid_cmd>
  <mid_cmd>color CONS8 @/bfactor=8</mid_cmd>
  <mid_cmd>color CONS7 @/bfactor=7</mid_cmd>
  <mid_cmd>color CONS6 @/bfactor=6</mid_cmd>
  <mid_cmd>color CONS5 @/bfactor=5</mid_cmd>
  <mid_cmd>color CONS4 @/bfactor=4</mid_cmd>
  <mid_cmd>color CONS3 @/bfactor=3</mid_cmd>
  <mid_cmd>color CONS2 @/bfactor=2</mid_cmd>
  <mid_cmd>color CONS1 @/bfactor=1</mid_cmd>
  <mid_cmd>preset apply pub 3;repr cpk;show;focus;color red ligand</mid_cmd>""")

    CHIMERAX.write("\n  </commands>\n</ChimeraPuppet>")
    CHIMERAX.close()
    os.chmod(chimerax_file, 0o755)

    return("OK")

def create_chimera_script(pdb_file_name, www_dir, chimerax_file, msa_file, tree_file, scf_file, hdr_file):

    try:

        CHIMERAX = open(chimerax_file, 'w')

    except:

        return("Failed to open " + chimerax_file + " for writing.")

    CHIMERAX.write("""<?xml version="1.0"?>
  <ChimeraPuppet type="std_webdata">
<web_files>
<file  name="%s" format="text" loc="%s"/>""" %(pdb_file_name, www_dir + pdb_file_name))

    if msa_file != "":

        CHIMERAX.write("""<file  name=\"%s\" format=\"text\" loc=\"%s\"/>
</web_files>
<commands>
  <mid_cmd>preset apply pub 3;repr cpk;show;focus;color red ligand</mid_cmd>

<!-- the following 3 lines locate the Multalign Viewer instance
    that was created by opening the alignment file, and stores a reference
    to the instance as the variable 'mav' -->
  <py_cmd>from MultAlignViewer.MAViewer import MAViewer</py_cmd>
  <py_cmd>from chimera.extension import manager</py_cmd>
  <py_cmd>mav = [inst for inst in manager.instances if isinstance(inst, MAViewer)][-1]</py_cmd>
<!-- read in/show Consurf tree -->

<!-- hide initial headers, show phylogeny tree, load the coloring file,
    and make residue letters black.
    This uses two possible sets of calls:  one for the 1.2540 release
    and one for later versions that uses a better API.
    The 'if' condition guarantees that the code will work no
    matter what version the user has -->
<py_cmd>
if hasattr(mav, 'loadScfFile'):
    mav.hideHeaders(mav.headers(shownOnly=True))
    mav.usePhylogenyFile("%s", askReorder=False)
    mav.loadScfFile("%s")
    mav.useColoringFile(None)
else:
    mav.hideHeaders(mav.seqCanvas.headerDisplayOrder())
    mav.usePhylogenyFile("%s")
    mav.regionBrowser.loadScfFile("%s")
    from MultAlignViewer.prefs import RC_BLACK
    mav.seqCanvas.setColorFunc(RC_BLACK)
</py_cmd>
<!-- read in/show Consurf headers -->
    <py_cmd>mav.readHeaderFile(\"%s\")</py_cmd>
<!-- show chains other than the one in the alignment as gray ribbon -->
<py_cmd>
import chimera
m = chimera.openModels.list(modelTypes=[chimera.Molecule])[0]
for seq in m.sequences():
    if hasattr(seq.residues[0], "mavConSurfHistogram"):
        continue
    from chimera.colorTable import getColorByName
    gray = getColorByName("gray")
    for r in seq.residues:
        for a in r.atoms:
            a.display = False
        r.ribbonDisplay = True
        r.ribbonDrawMode = chimera.Residue.Ribbon_Round
        r.ribbonColor = gray
</py_cmd>""" %(msa_file, www_dir + msa_file, www_dir + tree_file, www_dir + scf_file, www_dir + tree_file, www_dir + scf_file, www_dir + hdr_file))

    else:

        CHIMERAX.write("""</web_files>
<commands>
  <mid_cmd>colordef CONS10 1.00 1.00 0.59</mid_cmd>
  <mid_cmd>colordef CONS9 0.63 0.15 0.38</mid_cmd>
  <mid_cmd>colordef CONS8 0.94 0.49 0.67</mid_cmd>
  <mid_cmd>colordef CONS7 0.98 0.79 0.87</mid_cmd>
  <mid_cmd>colordef CONS6 0.99 0.93 0.96</mid_cmd>
  <mid_cmd>colordef CONS5 1.00 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS4 0.92 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS3 0.84 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS2 0.55 1.00 1.00</mid_cmd>
  <mid_cmd>colordef CONS1 0.06 0.78 0.82</mid_cmd>
  <mid_cmd>color CONS10 @/bfactor=0</mid_cmd>
  <mid_cmd>color CONS9 @/bfactor=9 </mid_cmd>
  <mid_cmd>color CONS8 @/bfactor=8</mid_cmd>
  <mid_cmd>color CONS7 @/bfactor=7</mid_cmd>
  <mid_cmd>color CONS6 @/bfactor=6</mid_cmd>
  <mid_cmd>color CONS5 @/bfactor=5</mid_cmd>
  <mid_cmd>color CONS4 @/bfactor=4</mid_cmd>
  <mid_cmd>color CONS3 @/bfactor=3</mid_cmd>
  <mid_cmd>color CONS2 @/bfactor=2</mid_cmd>
  <mid_cmd>color CONS1 @/bfactor=1</mid_cmd>
  <mid_cmd>preset apply pub 3;repr cpk;show;focus;color red ligand</mid_cmd>""")

    CHIMERAX.write("\n  </commands>\n</ChimeraPuppet>")
    CHIMERAX.close()
    os.chmod(chimerax_file, 0o755)

    return("OK")

def create_pymol_page(PyMol_Instruction_Page, PDB_ID, Chain_ID, Conf_Link, PyMol_ConSurf_Script, PDB_With_ConSurf_Scores, PDB_With_ConSurf_Scores_Showing_ISD, identical_chains, run_number, working_dir, chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, chimerax_script_for_figure, chimerax_script_for_figure_isd, pdb_file_name_NoPath, PyMol_ConSurf_Script_CBS):

    server_html_dir = GENERAL_CONSTANTS.CONSURF_HTML_DIR
    consurf2016_url = GENERAL_CONSTANTS.CONSURF_URL + "2016/"

    try:

        PYMOL_PAGE = open(PyMol_Instruction_Page, 'w')

    except:

        return("cp_rasmol_gradesPE_and_pipe::create_pymol_page: Can't open PyMol Output '" + PyMol_Instruction_Page + "' for writing.\n")

    PYMOL_PAGE.write("""<?php
\$TITLE_LINE = "ConSurf Image using PyMOL for %s Chain: %s";

include (\"%stemplates/definitions.tpl\");
include (\"%stemplates/header.tpl\");
?>
<center> <h2>Create a high resolution PyMol figure for %s Chain: %s</h2></center>
It is recommended to create a new directory where the attached files will be saved.<br>
Please note: For all PyMOL commands, you should type absolute paths for the files you are about to download, from their new location (on your computer).<br><br>""" %(PDB_ID, Chain_ID, server_html_dir, server_html_dir, PDB_ID, Chain_ID))

    if PDB_With_ConSurf_Scores_Showing_ISD != "":

        PYMOL_PAGE.write("<br> You may view your protein colored according to conservation scores with a unique color for positions of low confidence cutoff. If this is your preferred option, please choose option 1(a). Otherwise - choose option 1(b).<br>\n")
        PYMOL_PAGE.write("Click <a href = \"%s\">HERE</a> for more information regarding confidence cutoff.\n" %Conf_Link)
        PYMOL_PAGE.write("<br><br><ol><li>Download the <b>PDB_FILE</b> updated with ConSurf\'s colors.<br>\n")
        PYMOL_PAGE.write("(a) <a href= \"%s\">PDB_FILE</a> showing Insufficient Data<br>\n" %PDB_With_ConSurf_Scores_Showing_ISD)
        PYMOL_PAGE.write("(b) <a href= \"%s\">PDB_FILE</a> hiding Insufficient Data\n" %PDB_With_ConSurf_Scores)
        PYMOL_PAGE.write("</li>\n")

    else:

        PYMOL_PAGE.write("<br><br><ol><li>Download the <b>PDB_FILE</b> updated with ConSurf\'s colors.<br>\n")
        PYMOL_PAGE.write("<a href= \"%s\">PDB_FILE</a> hiding Insufficient Data\n" %PDB_With_ConSurf_Scores)

    PYMOL_PAGE.write("""<li>Download the file <a href= \"%s\" target =\"PyMOL_py\"><b>consurf_new.py</b></a> [or the color-blind freindly version: <a href= \"%s\" target =\"PyMOL_py\"><b>consurf_new_CBS.py</b></a>].</li>
<li>Start PyMOL.</li>
<li>A. Load the <font color=red> Consurf's modified pdb file (not the original pdb file that used as an input)</font><br>In the PyMOL viewer window type:<br>
<font face=Courier New>PyMOL>\"load <b> PDB_FILE</b>\"
(no quotes) and hit return.<br></font>
</li>
B. Run the script to define ConSurf\'s color; Type:<br>
<font face=Courier New>PyMOL>\"run <b> consurf_new.py</b>\" (no quotes) and hit return.<br><font color=red><b>or</b></font><br><font face=Courier New>PyMOL>\"run <b> consurf_new_CBS.py</b>\"
      (no quotes) and hit return.<br> </font>
     <br><font color=red>
     Note: Please don't forget to write the path of the script file unless the script file is located at the pymol's root directory.
     <br></font>

</font>
</ol>""" %(PyMol_ConSurf_Script, PyMol_ConSurf_Script_CBS))

    if identical_chains != chain and chain != "NONE":

        PYMOL_PAGE.write("<form action=\"%sproject_on_identical_chains.php\" method=\"post\">" %consurf2016_url)
        PYMOL_PAGE.write("<br>Act on identical chains: ")

        for identical_chain in identical_chains.split():

            if identical_chain == chain:

                PYMOL_PAGE.write("&nbsp<input type=checkbox name=\"identicalChainToProjectOn[]\" value=%s id=\"identicalChainToProjectOn%s\" checked=\"checked\">%s  " %(identical_chain, identical_chain, identical_chain))

            else:

                PYMOL_PAGE.write("&nbsp<input type=checkbox name=\"identicalChainToProjectOn[]\" value=%s id=\"identicalChainToProjectOn%s\">%s  " %(identical_chain, identical_chain, identical_chain))

        # hidden params
        PYMOL_PAGE.write("<input type=\"hidden\" name=\"run_number\" value=%s id=\"run_number\">" %run_number)
        PYMOL_PAGE.write("<input type=\"hidden\" name=\"working_dir\" value=%s id=\"working_dir\">" %working_dir)
        PYMOL_PAGE.write("<input type=\"hidden\" name=\"chain\" value=%s id=\"chain\">" %chain)
        PYMOL_PAGE.write("<input type=\"hidden\" name=\"ATOMS_with_ConSurf_Scores\" value=%s id=\"ATOMS_with_ConSurf_Scores\">" %ATOMS_with_ConSurf_Scores)
        PYMOL_PAGE.write("<input type=\"hidden\" name=\"ATOMS_with_ConSurf_Scores_isd\" value=%s id=\"ATOMS_with_ConSurf_Scores_isd\">" %ATOMS_with_ConSurf_Scores_isd)
        PYMOL_PAGE.write("<input type=\"hidden\" name=\"chimerax_script_for_figure\" value=%s id=\"chimerax_script_for_figure\">" %chimerax_script_for_figure)
        PYMOL_PAGE.write("<input type=\"hidden\" name=\"chimerax_script_for_figure_isd\" value=%s id=\"chimerax_script_for_figure_isd\">" %chimerax_script_for_figure_isd)
        PYMOL_PAGE.write("<input type=\"hidden\" name=\"pdbFileName\" value=%s id=\"pdbFileName\">" %pdb_file_name_NoPath)
        PYMOL_PAGE.write("&nbsp&nbsp&nbsp")
        PYMOL_PAGE.write("<input type=\"submit\" value=\"Calculate & Download\">")
        PYMOL_PAGE.write("</form>")

    PYMOL_PAGE.write("<br>")

    PYMOL_PAGE.write("""<?php
    include(\"%stemplates/footer.tpl\");
?>
</body></html>""" %server_html_dir)

    return("OK")

def create_chimera_page(chimera_instructions_file, conf_link, chimera_consurf_commands, chimerax_script, PDB_atoms_ConSurf, isd_chimerax_script, isd_PDB_atoms_ConSurf, url, identical_chains, run_number, working_dir, chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, chimerax_script_for_figure, chimerax_script_for_figure_isd, pdb_file_name_NoPath, chimera_consurf_commands_CBS):

    server_html_dir = GENERAL_CONSTANTS.CONSURF_HTML_DIR
    consurf2016_url = GENERAL_CONSTANTS.CONSURF_URL + "2016/"
    CHIMERA_SAVING_FIGURE_LINK = "http://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/print.html"

    # Create the html describing how to create Chimera Hige resolution images

    try:

        CHIMERA_HTML = open(chimera_instructions_file, 'w')

    except:

        return("Cann't open the file " + chimera_instructions_file + " for writing.\n")

    CHIMERA_HTML.write("""<?php
\$TITLE_LINE = "ConSurf Image using Chimera";

include (\"%stemplates/definitions.tpl\");
include (\"%stemplates/header.tpl\");
?>
    <center><h2>Create a high resolution Chimera figure</h2></center>""" %(server_html_dir, server_html_dir))

    if isd_PDB_atoms_ConSurf != "":

        CHIMERA_HTML.write("""	You may view your protein colored according to conservation scores with a unique color for positions of low confidence cutoff. If this is your preferred option, please choose option 1(a). Otherwise - choose option 1(b).
    <br>
    Click <a href = "%s">HERE</a> for more information regarding confidence cutoff.
    <br /><br />

        <b><font color="green">You can open the molecule with Chimera on a click!</font></b>
    <ol>
        <li>(a) <a href = "%s" type="application/x-chimerax">Click here</a> for PDB showing Insufficient Data<br />
        (b) <a href = "%s" type="application/x-chimerax">Click here</a> for PDB hiding Insufficient Data</li>
    </ol>

""" %(conf_link, url + isd_chimerax_script, url + chimerax_script))

    else:

        CHIMERA_HTML.write("<b><font color=\"green\">You can open the molecule with Chimera <a href=\"" + chimerax_script + "\" type=\"application/x-chimerax\">on a click!</a></font></b><br /><br />")

    CHIMERA_HTML.write("    <b><font color=\"green\">Alternatively, you can download the files and view it in a later stage.</font></b><br/>\n\n\tIt is recommended to create a new directory where the attached files will be saved.<br>\n\t<ol>")

    if isd_PDB_atoms_ConSurf != "":

        CHIMERA_HTML.write("""
        <li>Download the  <b>PDB_FILE</b> updated with ConSurf's colors.<br>\n
        (a) <a href=%s>PDB_FILE</a> showing Insufficient Data<br>\n
        (b) <a href=%s>PDB_FILE</a> hiding Insufficient Data
        </li>\n""" %(isd_PDB_atoms_ConSurf, PDB_atoms_ConSurf))

    else:

        CHIMERA_HTML.write("<li>Download the <a href=" + PDB_atoms_ConSurf + "><b>PDB_FILE</b></a> updated with ConSurf's colors.</li>\n")

    CHIMERA_HTML.write("""		<li>Download the <a href="%s"><b>coloring script</b></a> or <a href="%s"><B>the coloring script with color-blind friendly scale</B></a> for chimera.</li>
        <li>Start Chimera program.</li>
        <li>Load the PDB_FILE from the File/Open menu.</li>
        <li>Load the coloring script from the File/Open menu.</li>
    </ol>
    <br />
    Once the molecule is loaded and coloured, you may save a high quality figure following <a href=%s>these instructions</a>.<br />

""" %(chimera_consurf_commands, chimera_consurf_commands_CBS, CHIMERA_SAVING_FIGURE_LINK))

    CHIMERA_HTML.write("<br>")

    if identical_chains != chain and chain != "NONE":

        CHIMERA_HTML.write("<form action=\"%sproject_on_identical_chains.php\" method=\"post\">" %consurf2016_url)
        CHIMERA_HTML.write("<br>Act on identical chains: ")

        for identical_chain in identical_chains.split():

            if identical_chain == chain:

                CHIMERA_HTML.write("&nbsp<input type=checkbox name=\"identicalChainToProjectOn[]\" value=%s id=\"identicalChainToProjectOn%s\" checked=\"checked\">%s  " %(identical_chain, identical_chain, identical_chain))


            else:

                CHIMERA_HTML.write("&nbsp<input type=checkbox name=\"identicalChainToProjectOn[]\" value=%s id=\"identicalChainToProjectOn%s\">%s  " %(identical_chain, identical_chain, identical_chain))

        # hidden params
        CHIMERA_HTML.write("<input type=\"hidden\" name=\"run_number\" value=%s id=\"run_number\">" %run_number)
        CHIMERA_HTML.write("<input type=\"hidden\" name=\"working_dir\" value=%s id=\"working_dir\">" %working_dir)
        CHIMERA_HTML.write("<input type=\"hidden\" name=\"chain\" value=%s id=\"chain\">" %chain)
        CHIMERA_HTML.write("<input type=\"hidden\" name=\"ATOMS_with_ConSurf_Scores\" value=%s id=\"ATOMS_with_ConSurf_Scores\">" %ATOMS_with_ConSurf_Scores)
        CHIMERA_HTML.write("<input type=\"hidden\" name=\"ATOMS_with_ConSurf_Scores_isd\" value=%s id=\"ATOMS_with_ConSurf_Scores_isd\">" %ATOMS_with_ConSurf_Scores_isd)
        CHIMERA_HTML.write("<input type=\"hidden\" name=\"chimerax_script_for_figure\" value=%s id=\"chimerax_script_for_figure\">" %chimerax_script_for_figure)
        CHIMERA_HTML.write("<input type=\"hidden\" name=\"chimerax_script_for_figure_isd\" value=%s id=\"chimerax_script_for_figure_isd\">" %chimerax_script_for_figure_isd)
        CHIMERA_HTML.write("<input type=\"hidden\" name=\"pdbFileName\" value=%s id=\"pdbFileName\">" %pdb_file_name_NoPath)
        CHIMERA_HTML.write("&nbsp&nbsp&nbsp")
        CHIMERA_HTML.write("<input type=\"submit\" value=\"Calculate & Download\">")
        CHIMERA_HTML.write("</form>")

    CHIMERA_HTML.write("<br>")

    CHIMERA_HTML.write("""<?php
    include(\"%stemplates/footer.tpl\");
?>
</body></html>

""" %server_html_dir)

    CHIMERA_HTML.close()

    return("OK")

def create_gradesPE_ConSeq(Output, ref_residue_freq, ref_Solv_Acc_Pred, gradesPE_file, type, ALGORITHM, layers_array):

    # printing the the ConSurf gradesPE file

    try:

        PE = open(gradesPE_file, 'w')

    except:

        return("create_gradesPE_Conseq : can't open " + gradesPE_file + " for writing.")


    PE.write("\t %s Acid Conservation Scores\n" %type)
    PE.write("\t===============================\n\n")
    PE.write("The layers for assigning grades are as follows.\n")
    for i in range(1, len(layers_array)):
        
        if layers_array[i - 1] < 0:
        
            left_end = "%.3f" %layers_array[i - 1]
            
        else:
            
            left_end = " %.3f" %layers_array[i - 1]
    
            
        if layers_array[i] < 0:
        
            right_end = "%.3f" %layers_array[i]
            
        else:
            
            right_end = " %.3f" %layers_array[i]
            
        PE.write("from %s to %s the grade is %d\n" %(left_end, right_end, 10 - i))
        
    PE.write("\nIf the difference between the colors of the CONFIDENCE INTERVAL COLORS is more than 3 or the msa number (under the column titled MSA) is less than 6, there is insufficient data and an * appears in the COLOR column.\n")

    PE.write("\n- POS: The position of the %s Acid in the SEQRES derived sequence.\n" %type)
    PE.write("- SEQ: The SEQRES derived sequence in one letter code.\n")
    PE.write("- SCORE: The normalized conservation scores.\n")
    PE.write("- COLOR: The color scale representing the conservation scores (9 - conserved, 1 - variable).\n")
    PE.write("- CONFIDENCE INTERVAL: When using the bayesian method for calculating rates, a confidence interval is assigned to each of the inferred evolutionary conservation scores.\n")
    PE.write("- CONFIDENCE INTERVAL COLORS: When using the bayesian method for calculating rates. The color scale representing the lower and upper bounds of the confidence interval.\n")
    PE.write("- B/E: Burried (b) or Exposed (e) residue.\n")
    PE.write("- FUNCTION: functional (f) or structural (s) residue (f - highly conserved and exposed, s - highly conserved and burried).\n")
    PE.write("- MSA DATA: The number of aligned sequences having an %s acid (non-gapped) from the overall number of sequences at each position.\n" %type.lower())
    PE.write("- RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.\n\n")
    
    if ALGORITHM == "Bayes":
        
        CONFIDENCE_INTERVAL = "\tCONFIDENCE INTERVAL\tCONFIDENCE INTERVAL COLORS"
        
    else:
        
        CONFIDENCE_INTERVAL = ""
        
    PE.write(" POS\t SEQ\tSCORE\t\tCOLOR" + CONFIDENCE_INTERVAL)
    if type == "Amino":
        
        PE.write("\tB/E\tFUNCTION")
            
    PE.write("\t   MSA DATA  RESIDUE VARIETY\n    \t    \t(normalized)\t        \t               \n")

    for elem in Output:

        pos = elem['POS']
        var = ""
        found = False
        """
        if str(pos) in ref_Solv_Acc_Pred:

            Solv_Acc_Pred = ref_Solv_Acc_Pred[str(pos)]
        """

        score = elem['COLOR']

        PE.write("%4d" %pos)
        PE.write("\t%4s" %elem['SEQ'])
        PE.write("\t%6.3f" %elem['GRADE'])

        if elem['ISD'] == 1:

            PE.write("\t\t%3d" %elem['COLOR'])
            PE.write("*")

        else:

            PE.write("\t\t%3d" %elem['COLOR'])
            
        if ALGORITHM == "Bayes":

            PE.write("\t%6.3f" %elem['INTERVAL_LOW'])
            PE.write(", ")
            PE.write("%6.3f" %elem['INTERVAL_HIGH'])
            PE.write("\t\t\t%5d" %ColorScale[elem['INTERVAL_LOW_COLOR']])
            PE.write(",")
            PE.write("%1d\t\t" %ColorScale[elem['INTERVAL_HIGH_COLOR']])
            
        if type == "Amino":
            
            Solv_Acc_Pred = "" #ref_Solv_Acc_Pred[pos]

            PE.write("\t%3s" %Solv_Acc_Pred)

            # FUNCT/STRUCT COL
            if Solv_Acc_Pred == "e":

                if elem['COLOR'] == 9 or elem['COLOR'] == 8:

                    PE.write("\tf\t")

                else:

                    PE.write("\t\t ")

            elif elem['COLOR'] == 9:

                PE.write("\ts\t")

            else:

                PE.write("\t\t ")

        PE.write("\t%8d/%s" %(elem['MSA_NUM'], elem['MSA_DENUM']))

        keys = list((ref_residue_freq[pos].keys()))
        keys.sort()
        for aa in keys:

            if var == "":

                var += aa

            else:

                var += ", " + aa

            if (elem['SEQ']).upper() == aa.upper():

                found = True

        PE.write("\t" + var + "\n")

        # the amino-acid in that position, must be part of the residue variety in this column
        if not found:

            PE.close()
            return ("create_gradesPE_ConSeq : in position %d, the %s-acid %s does not match the residue variety: %s." %(pos, type.lower(), elem['SEQ'], var))

        # printing the residue to the rasmol script
        # assigning grades to seq3d strings

    PE.write("\n\n*Below the confidence cut-off - The calculations for this site were performed on less than 6 non-gaped homologue sequences,\nor the confidence interval for the estimated score is equal to- or larger than- 4 color grades.\n")
    PE.close()
    return("OK")

