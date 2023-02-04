
import re
import os
from Bio import SeqIO
from Bio import AlignIO


def determine_msa_format(working_dir):

    ret = []
    first_char = ""
    rest = ""
    found_pir = 0

    try:

        MSA = open(working_dir, 'r')

    except:

        return ["err", "cannot open the file " + working_dir]

    line = MSA.readline()
    while line != "":

        if re.search(r'\S', line):

            match = re.search(r'^\s*(\S)(.*)', line)
            first_char = match.group(1)
            rest = match.group(2)

            if first_char == 'M':

                if rest and re.search(r'^SF:' ,rest):

                    ret = ["format", "gcg"]

            elif first_char == '>':

                if rest and re.search(r'^P1;', rest):

                    line = MSA.readline()

                    while line != "":

                        if re.search(r'^>', line):

                            ret = ["format", "fasta"]
                            break

                        if re.search(r'\*\s*$', line):

                            found_pir = 1
                            ret = ["format", "pir"]
                            break

                else:

                    ret = ["format", "fasta"]
                    break

            elif re.search(r'\s*MSF:\s+', line):

                ret = ["format", "gcg"]
                break

            break

    MSA.close()
    if len(ret) == 0:

        if first_char == '#':

            ret = ["format", "nexus"]

        elif first_char == 'C':

            ret = ["format", "clustal"]

        elif first_char == 'C':

            ret = ["format", "gcg"]

        else:

            ret = ["err", "unknown_format"]

    return ret

def check_msa_licit_and_size(msa_file, msa_format):

    # if the format is "gcg" pass "msf" to AlignIO
    if msa_format == "gcg":

        msa_format = "msf"

    # if the format is "clustalw" pass "clustal" to AlignIO
    if msa_format == "clustalw":

        msa_format = "clustal"

    ans = ""
    msa_fh = ""
    seq = ""
    seq_name = ""
    ret = []
    read_seq = 0

    if msa_format == "fasta" or msa_format == "nexus" or msa_format == "pir":

        try:

            MSA = open(msa_file, 'r')

        except:

            return ['user_error', "could not read msa"]

        try:

            for sequence in SeqIO.parse(MSA, msa_format):

                read_seq += 1
                ans = check_sequence_licit(sequence.seq)

                if ans != "":

                    break

        except:

            MSA.close()
            return ['user_error', "exception"]

        MSA.close()

    elif msa_format == "clustal" or msa_format == "msf":

        try:

            MSA = open(msa_file, 'r')

        except:

            return ['user_error', "could not read msa"]

        try:

            for alignment in AlignIO.parse(MSA, msa_format):

                for record in alignment:

                    read_seq += 1
                    ans = check_sequence_licit(str(record.seq))
                    if ans != "":

                        break

                if ans != "":

                    break

        except:

            MSA.close()
            return ['user_error', "exception"]

        MSA.close()

    else:

        return ['user_error', "could not read msa"]

    if ans != "":

        ret = ['user_error', "SEQ_NAME: " + seq_name, "IRR_CHAR: " + ans]

    elif read_seq == 0:

        ret = ['user_error', "could not read msa"]

    else:

        ret = ["OK", read_seq]

    return ret





def get_info_from_msa(msa_file, msa_format, ref_hash_seq_names):    #, msa_out_file = "", file_header_msa = ""):

    # if the format is "gcg" pass "msf" to AlignIO
    if msa_format == "gcg":

        msa_format = "msf"

    # if the format is "clustalw" pass "clustal" to AlignIO
    if msa_format == "clustalw":

        msa_format = "clustal"

    ans = ""
    msa_fh = ""
    seq = ""
    seq_name = ""
    ret = ["OK"]


    if msa_format == "fasta" or msa_format == "nexus" or msa_format == "pir":

        try:

            MSA = open(msa_file, 'r')

        except:

            return ['user_error', "could not read msa"]

        try:

            for sequence in SeqIO.parse(MSA, msa_format):

                if re.search(r'^\s*$', sequence.id):

                    ret = ['user_error', "no seq id"]
                    break

                elif sequence.id in ref_hash_seq_names:

                    ret = ['user_error', "duplicity " + sequence.id]
                    break

                else:

                    ref_hash_seq_names[sequence.id] = str(sequence.seq)

        except:

            MSA.close()
            return ['user_error', "exception"]

        MSA.close()

    elif msa_format == "clustal" or msa_format == "msf":

        try:

            MSA = open(msa_file, 'r')

        except:

            return ['user_error', "could not read msa"]

        try:

            for alignment in AlignIO.parse(MSA, msa_format):

                for record in alignment:

                    if re.search(r'^\s*$', record.id):

                        ret = ['user_error', "no seq id"]
                        break

                    elif record.id in ref_hash_seq_names:

                        ret = ['user_error', "duplicity " + record.id]
                        break

                    else:

                        ref_hash_seq_names[record.id] = str(record.seq)

                if ret[0] != "OK":

                    break

        except:

            MSA.close()
            return ['user_error', "exception"]

        MSA.close()

    else:

        return ['user_error', "could not read msa"]
    """
    if msa_out_file != "":

        if file_header_msa == "yes":

            try:

                MSA_OUT = open(msa_out_file, 'w')

            except:

                return ['user_error', "could not write msa out file"]

            for key in ref_hash_seq_names:

                MSA_OUT.write(key + "\n")

            MSA_OUT.close()

        else:

            try:

                MSA_OUT = open(msa_out_file, 'w')

            except:

                return ['user_error', "could not write msa out file"]

            for key in ref_hash_seq_names:

                MSA_OUT.write(key + "\n" + ref_hash_seq_names[key] + "\n")

            MSA_OUT.close()
    """
    return ret








# look for iregular chars in AA sequence
def check_sequence_licit(sequence):

    no_regular_format_char = ""
    no_regular_note = ""
    found = []

    for char in sequence:

        if not re.search(r'[ACDEFGHIKLMNPQRSTUVWXY\-]', char, re.IGNORECASE):

            if not char in found:

                if re.search(r'[BJOUZ]', char, re.IGNORECASE):

                    no_regular_note += "\"" + char + "\", "

                else:

                    no_regular_format_char += "\"" + char + "\", "

                found.append(char)

    no_regular_format_char += no_regular_note

    if no_regular_format_char:

        no_regular_format_char = re.sub(r', $', '', no_regular_format_char)

    return no_regular_format_char


def convert_msa_format(infile, infileformat, outfile, outfileformat):

    """
    try:

        MSA_INPUT = AlignIO.parse(infile, infileformat)

    except:

        return "could not open " + infile + " for reading."

    #try:

    MSA_OUTPUT = open(outfile, 'w')
    AlignIO.write(MSA_INPUT, MSA_OUTPUT, outfileformat)
    MSA_OUTPUT.close()

    #except:

    #    return "could not open " + outfile + " for writing."
    """
    """
    try:

        for aln in AlignIO.parse(MSA_INPUT, "clustal"):

            AlignIO.write(aln, MSA_OUTPUT, "fasta")

    except:

        ret = "exception"

    MSA_INPUT.close()
    MSA_OUTPUT.close()

    """
    try:

        AlignIO.convert(infile, infileformat, outfile, outfileformat)

    except:

        return "exception"

    if not os.path.exists(outfile) or os.path.getsize(outfile) == 0:

        return "MSA not created."

    return "OK"

def read_residue_variety(msa_file, msa_ref_sequence, msa_format, ref_residue_frequency, ref_position_totalAA):

    # the routine creates 2 hashes, according to information it reads from a MSA
    # for each position in the query sequence :
    # 1. it collects all the residues that are aligned to this positions. it also counts the number of time each residue appeared in that position in the MSA
    # 2. it counts the total number of residues which aligned to this position

    num_of_seqs = 0
    ret_seq = ""
    elements_to_remove = []

    # open msa file
    try:

        MSA = open(msa_file, 'r')

    except:

        return(["MSA_parser.read_residue_variety : could not read msa " + msa_format])

    try:

        for seq_record in SeqIO.parse(MSA, msa_format):

            num_of_seqs += 1
            ret_seq = set_position(seq_record.seq, seq_record.id, msa_ref_sequence, elements_to_remove, ref_residue_frequency, ref_position_totalAA)

    except:

        MSA.close()
        return(["MSA_parser.read_residue_variety : exception"])

    MSA.close()

    # remove all positions where the query seq was skipped
    while len(elements_to_remove) > 0:

        index = elements_to_remove.pop()
        residues_shift_left(index, ref_residue_frequency)
        residues_shift_left(index, ref_position_totalAA)

    # If there was an unexpected character in the query sequence in the MSA, report it and exit
    if ret_seq != "":

        return(ret_seq)

    else:

        return("OK", num_of_seqs)

def residues_shift_left(start_position, ref_residue_frequency):

    index = start_position
    while index + 1 in ref_residue_frequency:

        ref_residue_frequency[index] = ref_residue_frequency[index + 1]
        index += 1

    if index in ref_residue_frequency:

        del ref_residue_frequency[index]

def set_position(seq_data, seq_name, msa_ref_sequence, ref_elements_to_remove, ref_residue_frequency, ref_position_totalAA):

    # reads the given sequence, seq_data. For each position in the sequence, which is not a gap (-) put in the hash ref_residue_frequency the other residues from the MSA which aligns to this position

    ret_seq = ""
    position_in_MSA = 0
    for char in seq_data:

        position_in_MSA += 1

        # if this is the query seq and no aa is found we need to save
        # this position in order to erase it in the end
        if seq_name == msa_ref_sequence and not char.upper() in "ABCDEFGHIKLMNPQRSTVWY":

            ref_elements_to_remove.append(position_in_MSA)
            if not char in "-Xx":

                ret_seq += "in seq %s ignored %s in position %s " %(seq_name, char, position_in_seq)

        elif char != "-":

            # save residue value
            if not position_in_MSA in ref_residue_frequency:

                ref_residue_frequency[position_in_MSA] = {}
                ref_residue_frequency[position_in_MSA][char] = 1

            elif not char in ref_residue_frequency[position_in_MSA]:

                ref_residue_frequency[position_in_MSA][char] = 1

            else:

                ref_residue_frequency[position_in_MSA][char] += 1

            if not position_in_MSA in ref_position_totalAA:

                ref_position_totalAA[position_in_MSA] = 1

            else:

                ref_position_totalAA[position_in_MSA] += 1

    return ret_seq











