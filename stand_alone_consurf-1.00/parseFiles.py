
import re
import os
import subprocess
from Bio import SearchIO
from Bio import SeqIO

def choose_homologoues_from_search_with_lower_identity_cutoff(searchType, query, redundancyRate, frag_overlap, min_length_percent, min_id_percent, min_num_of_homologues, search_output, fasta_output, rejected_seqs, ref_search_hash, Nuc_or_AA):

    query_seq_length = len(query)
    sequences = {}

    # Defining the minimum length a homologue should have
    # 60% the query's length
    min_length = query_seq_length * min_length_percent

    # Reading blast/hmmer output and collect the homologues
    # Printing the selected homologues to a file and insert the e-value info to hash
    try:

        OUT_REJECT = open(rejected_seqs, 'w')

    except:

        return ("sys", "parseFiles.choose_homologoues_from_search_with_lower_identity_cutoff : can't open file " + rejected_seqs + " for writing\n")

    try:

        OUT = open(fasta_output, 'w')

    except:

        return("sys", "parseFiles.choose_homologoues_from_search_with_lower_identity_cutoff : can't open file " + fasta_output + " for writing\n")

    try:

        RAW_OUT = open(search_output, 'r')

    except:

        return("sys", "parseFiles.choose_homologoues_from_search_with_lower_identity_cutoff : can't open file " + search_output + " for reading\n")

    num_homologoues = 0
    num_rejected = 1
    final_num_homologues = 0
    if searchType == "HMMER":

        # Moving the start position to the last round of the output
        position = 0
        line = RAW_OUT.readline()
        while line != "":

            if line.rstrip() == "[ok]":

                RAW_OUT.seek(position)
                break

            else:

                match =re.match(r'^@@ Round:', line)
                if match:

                    position = RAW_OUT.tell()

            line = RAW_OUT.readline()

        # we skip the beginning
        line = RAW_OUT.readline()
        while line != "":

            if re.match(r'^>>\s*(\S+)\s+(\S.*)', line):

                break

            else:

                 line = RAW_OUT.readline()

        # we process the output
        while line != "":

            # each seq may have more than one domain

            match1 = re.match(r'^>>\s*(\S+)\s+(\S.*)', line)
            if match1:

                # new seq found

                num_homologoues += 1
                seq_name = match1.group(1)
                regex = r"^\s*" + re.escape(seq_name) + r"\s*(\S+)\s+(\S+)\s+(\S+)"
                #seq_name = change_name(seq_name)
                seq_description = match1.group(2)
                sequences = [] # The array will hold the beginning and ending of the domains
                line = RAW_OUT.readline()

                # we extract data from the domains
                while line != "":
                    
                    if Nuc_or_AA == "AA":
                        
                        match2 = re.match(r'.*E-value:\s+(\S+)\s*', line)
                        
                    else:
                        
                        match2 = re.match(r'^\s+!\s+\S+\s+\S+\s+(\S+)', line)
                        
                    if match2:

                        # new domain found

                        seq_eval = match2.group(1)
                        prev_line = line # the line between the db seq and the query seq contains the residues that appear in both
                        line = RAW_OUT.readline()

                        while line != "":

                            match3 = re.match(regex, line)
                            if match3:

                                seq_beg = int(match3.group(1))
                                seq = match3.group(2) # the seq in fasta with gaps
                                seq_end = int(match3.group(3))
                                seq_ident = (float(count_letters(prev_line)) / float(len(seq))) * 100
                                seq = re.sub(r'-', "", seq) # delete gaps

                                # deciding if we take the fragment
                                # in case there is already a fragemnt with the same name, we do another validity test for overlapping

                                if len(sequences) > 0:

                                    # there is more than one domain. check if there is overlap
                                    ans = check_if_seq_valid_with_min_id(redundancyRate, min_length, min_id_percent, seq_ident, seq, seq_name, Nuc_or_AA)
                                    if ans == "yes":

                                        ans = check_if_no_overlap(frag_overlap, sequences, seq_beg, seq_end)

                                else:

                                    ans = check_if_seq_valid_with_min_id(redundancyRate, min_length, min_id_percent, seq_ident, seq, seq_name, Nuc_or_AA)

                                # after taking the info, check if the currecnt sequence is valid. If so - insert it to the hash
                                if ans == "yes":

                                    # in case there is more than one fragment for this seq_name: add another hash to the details array
                                    sequences.append([seq_beg, seq_end])

                                    # printing the selected homologues to a file and insert the e-value info to hash
                                    seq_frag_name = ""
                                    match = re.match(r'(\|)$', seq_name)
                                    if not match:

                                        seq_frag_name = "%s_start_%d_end_%d_Evalue_%s" %(seq_name, seq_beg, seq_end, seq_eval)

                                    else:

                                        seq_frag_name = "%sFRAMENT_%d_%d" %(seq_name, seq_beg, seq_end)

                                    final_num_homologues += 1
                                    OUT.write(">%s | %s\n%s\n" %(seq_frag_name, seq_description, seq))
                                    ref_search_hash[seq_frag_name] = seq_eval

                                else:

                                    OUT_REJECT.write("%d Fragment %s_start_%d_end_%d rejected: %s\n" %(num_rejected, seq_name, seq_beg, seq_end, ans))
                                    num_rejected += 1

                                # end of domain, move to next domain
                                break

                            prev_line = line
                            line = RAW_OUT.readline()

                    elif line[:2] == ">>":

                        # end of seq, move to next seq
                        break

                    else:

                        line = RAW_OUT.readline()

    else:

        searchio = SearchIO.parse(search_output, "blast-xml")
        last_blast_round = ""
        for qresult in searchio:

            last_blast_round = qresult

        for hit in last_blast_round:

            s_description = hit.description
            #s_name = change_name(hit.id)
            s_name = hit.id
            seq_name_exists = "no"
            for hsp in hit:

                # hsp is the next available High Scoring Pair, Bio::Search::HSP::HSPI object or null if finished
                num_homologoues += 1
                # extracting relevant details from the fragment
                [s_beg, s_end] = hsp.hit_range
                s_beg += 1
                AAseq = str(hsp.hit.seq)
                AAseq = re.sub(r'-', "", AAseq)
                s_ident = (float(hsp.ident_num) / float(hsp.aln_span)) * 100

                """
                local_min_id_percent = min_id_percent
                if isinstance(min_id_percent, str) and (min_id_percent).upper() == "ROST":

                    align_length = hsp.aln_span
                    local_min_id_percent = get_min_identity_cutoff(align_length)
                """

                s_eval = str(hsp.evalue)
                s_eval = re.sub(r',', "", s_eval)
                match = re.match(r'^e', s_eval)
                if match:

                    s_eval = "1" + s_eval

                # deciding if we take the fragment
                # in case there is already a fragemnt with the same name, we do another validity test for overlapping

                if s_name in sequences:

                    seq_name_exists = "yes"
                    seq_details = sequences[s_name]
                    ans = check_if_seq_valid_with_min_id(redundancyRate, min_length, min_id_percent, s_ident, AAseq, s_name, Nuc_or_AA)
                    if ans == "yes":

                        ans = check_if_no_overlap(frag_overlap, seq_details, s_beg, s_end)

                else:

                    ans = check_if_seq_valid_with_min_id(redundancyRate, min_length, min_id_percent, s_ident, AAseq, s_name, Nuc_or_AA)

                # after taking the info, check if the currecnt sequence is valid. If so - insert it to the hash
                if ans == "yes":

                    # in case there is more than one fragment for this seq_name: add another hash to the details array
                    if seq_name_exists == "yes":

                        (sequences[s_name]).append([s_beg, s_end])

                    # in case it is the first fragment for this seq_name: insert a details array as a value for this seq_name key
                    else:

                        sequences[s_name] = [[s_beg, s_end]]

                    # Printing the selected homologues to a file and insert the e-value info to hash
                    seq_frag_name = ""
                    match = re.match(r'(\|)$', s_name)
                    if not match:

                        seq_frag_name = "%s_start_%s_end_%s_Evalue_%s" %(s_name, s_beg, s_end, s_eval)

                    else:

                        seq_frag_name = "%sFRAMENT_%s_%s" %(s_name, s_beg, s_end)

                    final_num_homologues += 1
                    OUT.write(">%s | %s\n%s\n" %(seq_frag_name, s_description, AAseq))
                    ref_search_hash[seq_frag_name] = s_eval

                else:

                    OUT_REJECT.write("Fragment %s_%s_%s rejected: %s\n" %(s_name, s_beg, s_end, ans))

    OUT_REJECT.close()
    OUT.close()
    # Checking that the number of homologues found is legal

    ret = ""
    if final_num_homologues == 1:

        ret = "only <a href=\"<?=$orig_path?>/" + fasta_output + "\" style=\"color: #400080; text-decoration:underline;\">one unique sequence</a> "

    elif final_num_homologues < min_num_of_homologues:

        ret = "only <a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\">%d unique sequences</a> " %(fasta_output, final_num_homologues)

    else:

        ret ="ok"

    if ret != "ok":

        if searchType != "HMMER":

            ret += "were chosen from <a href=\"<?=$orig_path?>/" + search_output + "\" style=\"color: #400080; text-decoration:underline;\" download>BLAST output</a>." #CSI-BLAST???

        else:

            ret += "were chosen from <a href=\"<?=$orig_path?>/" + search_output + "\" style=\"color: #400080; text-decoration:underline;\" download>HMMER output</a>."

        if os.path.getsize(rejected_seqs) > 0:

            if searchType == "HMMER":

                ret += "(<a href=\"<?=$orig_path?>/" + rejected_seqs + "\" style=\"color: #400080; text-decoration:underline;\" download>Click here</a> if you wish to view the list of sequences which produced significant alignments in HMMER, but were not chosen as hits.)."

            else:

                ret += "(<a href=\"<?=$orig_path?>/" + rejected_seqs + "\" style=\"color: #400080; text-decoration:underline;\" download>Click here</a> if you wish to view the list of sequences which produced significant alignments in blast, but were not chosen as hits.)."

        ret += "<br>The minimal number of sequences required for the calculation is " + str(min_num_of_homologues) + ".<br>"
        return ["user", ret]

    return [ret, num_homologoues, final_num_homologues]



def change_name(s_name):


    match1 = re.match(r'.+\|(\S+\|\S+_\S+)', s_name)
    if match1:

        s_name = match1.group(1)

    match2 = re.match(r'.+\|(\S+_\S+)\|(\S+)', s_name)
    if match2:

        s_name = match2.group(2) + "|" + match2.group(1)

    else:

        match3 = re.match(r'(\S+\|\S+_\S+)', s_name)
        if match3:

            s_name = match3.group(1)

        else:

            match4 = re.match(r'(\S+_\S+)\|(\S+)', s_name)
            if match4:

                s_name = match4.group(2) + "|" + match4.group(1)

            else:

                match5 = re.match(r'(\S+_\S+)\|(\S+)', s_name)
                if match5:

                    s_name = match5.group(2) + "|" + match5.group(1)

                else:

                    match6 = re.match(r'(.*)\|([A-Za-z0-9._]+)\|', s_name)
                    if match6:

                        s_name = match6.group(2)

    return s_name




def check_if_seq_valid_with_min_id(redundancyRate, min_length, min_id, ident_percent, aaSeq, seqName, Nuc_or_AA):

    seq_length = len(aaSeq)
    ans = "yes"

    if ident_percent >= redundancyRate:

        # the sequence identity is not too high
        ans = "identity percent " + str(ident_percent) + " is too big"

    if ident_percent < min_id:

        # the sequence identity is higher than the minium idnentity percent that was defined for homologus
        ans = "identity percent %f is too low (below %d)" %(ident_percent, min_id)

    elif seq_length < min_length:

        # the sequnece length is greater than the minimum sequence length
        ans = "the sequence length " + str(seq_length) + " is too short. The minimum is " + str(min_length)

    # the sequnece letters should be legal to rate4site
    if Nuc_or_AA == "AA":

        # AA seq
        if not re.match(r'^[ACDEFGHIKLMNPQRSTVWYBZXacdefghiklmnpqrstvwybzx]+$', aaSeq):

            ans = "illegal character was found in sequence: " + seqName

    else:

        # Nuc seq
        if not re.match(r'^[ACGTUINacgtuin]+$', aaSeq):

            ans = "illegal character was found in sequence: " + seqName

    return ans

def check_if_no_overlap(max_overlap, ref_seq_details, s_bgn, s_end):

    ans = "check_if_no_overlap : no ans was picked"

    i = 0
    while i < len(ref_seq_details):

        fragment_beg = ref_seq_details[i][0]
        fragment_end = ref_seq_details[i][1]
        fragment_length = int(fragment_end) - int(fragment_beg) + 1

        if s_bgn <= fragment_beg and s_end >= fragment_end:

            # fragment is inside subjct
            return "previous fragment found %s_%s is fully inside new fragment" %(fragment_beg, fragment_end)

        elif s_bgn >= fragment_beg and s_end <= fragment_end:

            # subjct is inside fragment
            return "new fragment is fully inside previous fragment found " + str(fragment_beg + fragment_end)

        elif fragment_end < s_end and fragment_end > s_bgn:

            # fragment begins before subjct
            overlap_length = fragment_end - s_bgn + 1
            if overlap_length > fragment_length * max_overlap:

                return "overlap length of fragment is %d which is greater than maximum overlap: %d" %(overlap_length, fragment_length * max_overlap)

            else:

                # when the fragment might be a good match, we can only insert it if it did not match to all the fragments
                if i == len(ref_seq_details) - 1:

                    ans = "yes"

        elif fragment_beg > s_bgn and  fragment_beg < s_end:

            # fragment begins after subjct
            overlap_length = s_end - fragment_beg + 1
            if overlap_length > fragment_length * max_overlap:

                return "overlap length of fragment is %d which is greater than maximum overlap: %d" %(overlap_length, fragment_length * max_overlap)

            else:

                # when the fragment might be a good match, we can only insert it if it did not match to all the fragments
                if i == len(ref_seq_details) - 1:

                    ans = "yes"

        elif fragment_beg >= s_end or fragment_end <= s_bgn:

            # no overlap
            if i == len(ref_seq_details) - 1:

                ans = "yes"

        i += 1

    return ans


def create_cd_hit_output(input_file, output_file, cutoff, cd_hit_dir, ref_cd_hit_hash, type):


    seq = ""
    seq_name = ""
    cmd = []

    # running cd-hit

    if type == "AA":

        cmd.extend([cd_hit_dir + "cd-hit", "-i", input_file, "-o", output_file])

        if cutoff > 0.7 and cutoff < 1:

            cmd.extend(["-c", str(cutoff), "-n", "5"])

        elif cutoff > 0.6 and cutoff <= 0.7:

            cmd.extend(["-c", str(cutoff), "-n", "4"])

        elif cutoff > 0.5 and cutoff <= 0.6:

            cmd.extend(["-c", str(cutoff), "-n", "3"])

        elif cutoff > 0.4 and cutoff <= 0.5:

            cmd.extend(["-c", str(cutoff), "-n", "2"])

    else:

        # DNA
        cmd.extend([cd_hit_dir + "cd-hit-est", "-i", input_file, "-o", output_file])

        if cutoff > 0.9 and cutoff < 1:

            cmd.extend(["-c", str(cutoff), "-n", '8'])

        elif cutoff > 0.88 and cutoff <= 0.9:

            cmd.extend(["-c", str(cutoff), "-n", '7'])

        elif cutoff > 0.85 and cutoff <= 0.88:

            cmd.extend(["-c", str(cutoff), "-n", "6"])

        elif cutoff > 0.8 and cutoff <= 0.85:

            cmd.extend(["-c", str(cutoff), "-n", "5"])

        elif cutoff > 0.75 and cutoff <= 0.8:

            cmd.extend(["-c", str(cutoff), "-n", "4"])

    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    p.communicate()

    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:

        return("sys", "parseFiles.create_cd_hit_output : " + str(cmd) + ": CD-HIT produced no output!\n")

    num_cd_hits = 0
    # inserting chosen homologues to a hash
    for seq_record in SeqIO.parse(output_file, "fasta"):

        seq = seq_record.seq
        seq_name = seq_record.id
        description = seq_record.description
        if not seq_name in ref_cd_hit_hash:

            num_cd_hits += 1
            ref_cd_hit_hash[seq_name] =seq
        """
        ref_cd_hit_hash[seq_name]['SEQ'] = seq
        ref_cd_hit_hash[seq_name]['DESCRIPTION'] = description
        """
    return ["ok", num_cd_hits]

def add_sequences_removed_by_cd_hit_to_rejected_report(cd_hit_clusters_file, rejected_fragments_file):

    try:

        REJECTED = open(rejected_fragments_file, 'a')

    except:

        return "extract_sequences_removed_by_cd_hit: Can't open '" + rejected_fragments_file + "' for writing."

    try:

        CDHIT = open(cd_hit_clusters_file, 'r')

    except:

        return "extract_sequences_removed_by_cd_hit: Can't open '" + cd_hit_clusters_file + "' for reading.\n"

    cluster_members = {}
    cluster_head = ""

    line = CDHIT.readline()
    while line != "":

        match = re.match(r'^>Cluster', line)
        if match:

            # New Cluster
            for cluster_member in cluster_members.keys():

                REJECTED.write("Fragment %s rejected: the sequence shares %s identity with %s (%s was preserved)\n" %(cluster_member, cluster_members[cluster_member], cluster_head, cluster_head))

            cluster_members = {}
            cluster_head = ""

        else:

            # Clusters Members
            lines = line.split()
            if len(lines) > 2:

                x = re.sub(r'[>]', "", lines[2])
                if lines[3] == "*":

                    cluster_head  = x

                elif len(lines) > 3:

                    cluster_members[x] = lines[4]

        line = CDHIT.readline()

    return "ok"

def count_letters(s):

    l = 0
    for c in s:

        if c.isalpha():

            l += 1

    return l


def sort_sequences_from_eval(ref_to_seqs_hash, ref_to_cd_hash, max_num_homologs, witch_unifrom, output_file, input_file = ""):

    # input: reference to hash which holds sequences and their evalue
    # optional: query sequence input file (to extract the query sequence and add it)
    # output: a file with the top-evalued sequences out of the hash

    query_name = ""
    query_AAseq = ""
    counter = 1

    # open original seq fasta file
    if input_file != "" and os.path.exists(input_file):

        try:

            QUERY = open(input_file, 'r')

        except:

            return("err", "Can't open '" + input_file + "' for reading.")

        # get query name & amino acid seq from file
        line = QUERY.readline()
        while line != "":

            line = line.rstrip()
            match = re.match(r'^>(.+)', line)
            if match:

                query_name = match.group(1)

            else:

                query_AAseq += match.group()

            line = QUERY.readline()

        QUERY.close()

    # choose best homologs and create final file final_homologues_filename

    try:

        FINAL = open(output_file, 'a')

    except:

        return("err", "Can't open '" + output_file + "' for writing.")

    # write query details
    if query_AAseq != "":

        FINAL.write(">%s\n%s\n" %(query_name, query_AAseq))

    size_cd_hit_hash = len(ref_to_cd_hash)
    uniform = 1
    jump = 1

    if not witch_unifrom == "best":

        uniform = int(size_cd_hit_hash / max_num_homologs)
        if uniform == 0:

            uniform = 1

    final_number_of_homologoues = 1
    # write homologs
    for s_name in sorted(ref_to_seqs_hash.keys(), key = ref_to_seqs_hash.get):

        # quit if reached max number of homologs
        if counter > max_num_homologs * uniform:

            break

        # write next homolog
        if s_name in ref_to_cd_hash: # and 'SEQ' in ref_to_cd_hash[s_name]:

            if counter != jump:

                counter += 1
                continue
            """
            s_aa_sq = ref_to_cd_hash[s_name]['SEQ']
            s_description = ref_to_cd_hash[s_name]['DESCRIPTION']
            s_eval = ref_to_seqs_hash[s_name]
            FINAL.write(">%s | E_val=%s\n%s\n" %(s_description, s_eval, s_aa_sq))
            """
            final_number_of_homologoues += 1  
            FINAL.write(">%s\n%s\n" %(s_name, ref_to_cd_hash[s_name]))
            counter += 1
            jump += uniform

    FINAL.close()
    return("ok", "uniform is: " + str(uniform), final_number_of_homologoues)

def fasta_max_median_E_value(blst_f):

    # input fasta file with e-values in description.
    # output results[0] = min E-val (or -1 if empty), results[1] = median E-val

    e_vals = []
    for seq_record in SeqIO.parse(blst_f, "fasta"):

        match = re.match(r'E_val=(\S+)', seq_record.description)
        if match:

            e_vals.append(match.group(1))

    return_array = []
    if len(e_vals) > 0:

        sorted_vals = e_vals.sort(reverse = True)
        return_array.append(sorted_vals[0])
        return_array.append(median(e_vals))

    else:

        return_array.append(-1)
        return_array.append(-1)

    return return_array

def median(array):

    # median of an array of values
    vals = array.sort()
    len = len(vals)
    if len % 2:

        # The length is odd
        return vals[int(len / 2)]

    else:

        # the length is even
        return (vals[int(len / 2) - 1] + vals[int(len / 2)]) / 2

def strip_html(HTML_File, PlainText_File):

    # strip HTML tags - Not perfect - not working for multilines tag
    try:

        HTML = open(HTML_File, 'r')

    except:

        return("parseFiles.strip_html: Can not Open the HTML File: " + HTML_File + " for reading.")

    try:

        PLAIN_TEXT = open(PlainText_File, 'w')

    except:

        return("parseFiles.strip_html: Can not open " + PlainText_File + " for writing.")

    line = HTML.readline()
    while line != "":

        line = re.sub(r"<(?:[^>'\"]*|(['\"]).*?\1)*>", "", line)
        PLAIN_TEXT.write(line)
        line = HTML.readline()

    HTML.close()
    PLAIN_TEXT.close()
    return "ok"

def print_blast_according_to_round(blast_output_filename, round_number, output_filename, round_found):

    # prints the content of the file from blast, starting from round_number
    # Returns:
    #	If error : ("err", error_description)
    #	Else:  ("no_err", 0 if round not found; 1 otherwise)

    # open files
    try:

        BLAST = open(BLAST, 'r')

    except:

        return("err", "Can't open file '" + blast_output_filename + "' for reading.")

    try:

        NEW_FILE = open(output_filename, 'w')

    except:

        return("err", "Can't open file '" + output_filename + "' for writing.")

    # read blast file
    line = BLAST.readline()
    while line != "":

        match = re.match(r'Results from round (\d+)', line)
        if match:

            if round_number == match.group(1):

                output_data = 1
                round_found = 1

            else:

                output_data = 0

        if output == 1:

            NEW_FILE.write(line)

    # close files
    NEW_FILE.close()
    BLAST.close()

    return("no_err", round_found)





