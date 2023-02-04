
import GENERAL_CONSTANTS
import os
import re

def check_if_rate4site_failed(res_flag, r4s_log):

    # There are some tests to see if rate4site failed.
    # Since we can't trust only one of them, we do all of them. If one of them is tested to be true - than a flag will get TRUE value
    # 1. the .res file might be empty.
    # 2. if the run failed, it might be written to the log file of r4s.
    # 3. in a normal the r4s.log file there will lines that describe the grades. if it fail - we won't see them
    # In one of these cases we try to run the slower version of rate4site.
    # We output this as a message to the user.

    return_ans = "no"
    print_to_html = ""
    print_to_log = ""
    r4s_process_id = ""
    error_found = "no"

    if not os.path.exists(res_flag):

        print_to_log = "rate4site_routines.check_if_rate4site_failed : the file " + res_flag + " does not exsits. \n"
        error_found = "yes"

    elif os.path.exists(res_flag) and os.path.getsize(res_flag) == 0: # 1

        print_to_log = "rate4site_routines.check_if_rate4site_failed : the file " + res_flag + " was found to be of size 0. \n"
        error_found = "yes"

    if os.path.exists(res_flag):

        try:

            R4S_RES = open(res_flag, 'r')

        except:

            print_to_log = "rate4site_routines.check_if_rate4site_failed : can not open file: " + res_flag + ". aborting.\n"
            error_found = "yes"

        line = R4S_RES.readline()
        while line != "":

            match = re.match(r'In the tree file there is the name: (.+) that is not found in the sequence file', line)
            if match:

                print_to_html += "The sequence name %s was found in the tree file, but was not found in your MSA.<br>\nPlease correct your tree file, so it will include the same names as they appear in the MSA and re-run your query.<br>\n" %match.group(1)
                error_found = "yes"
                print_to_log = "rate4site_routines.check_if_rate4site_failed : sequence name %s was found in the tree file, was not found in the MSA\n" %match.group(1)
                break

            line = R4S_RES.readline()

        R4S_RES.close()

    if os.path.exists(r4s_log) and os.path.getsize(r4s_log) == 0: # 2, 3

        try:

            R4S_RES = open(res_flag, 'r')

        except:

            print_to_log = "rate4site_routines.check_if_rate4site_failed : can not open file: " + res_flag + ". aborting.\n"
            error_found = "yes"

        line = R4S_RES.readline()
        while line != "":

            match1 = re.march(r'/^.Process_id= (\d+)', line)
            if match1:

                r4s_process_id = match1.group(1)

            match2 = rem.match(r'likelihood of pos was zero', line)
            if match2:

                print_to_log = "rate4site_routines.check_if_rate4site_failed : the line: \"likelihood of pos was zero\" was found in %s.\n" %r4s_log
                error_found = "yes"
                break

            match3 = re.match(r'rate of pos\:\s\d\s=', line)
            if match3:

                # if we see this line, we change the flag
                return_ans = "no"
                break

            match4 = re.match(r'The amino-acid sequences contained the character: (.+)', line)
            if match4:

                print_to_html += "The illegal character %s was found in your MSA. Please make sure your MSA contains only the following characters:<br />\nA, B, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, X, Z, -" %match4.group(1)
                print_to_log = "rate4site_routines.check_if_rate4site_failed : illegal character %s was found in the MSA\n" %match4.group(1)
                error_found = "yes"
                break

            match5 = re.match(r'Could not find a sequence that matches the sequence name', line)
            if match5:

                seq_name = R4S_RES.readline()
                print_to_log = "rate4site_routines.check_if_rate4site_failed : the submitted query sequence name %s was not found in MSA" %seq_name
                print_to_html += "The query sequence name %s you have submitted was not found in the uploaded MSA. Please note that the sequence name should be exactly as it is in the MSA file<br>" %seq_name
                error_found = "yes"
                break

            match6 = re.match(r'The sequence name: (.+)was found in the tree file but not found in the sequence file', line)
            if match6:

                seq_name = match6.group(1)
                print_to_log = "rate4site_routines.check_if_rate4site_failed : the sequence name %s was found in the tree file, but not in the MSA" %seq_name
                print_to_html += " The tree file is inconsistant with the uploaded MSA. The sequence: \'%s\' was found in the tree file, but was not found in the MSA.<br>" %seq_name
                error_found = "yes"
                break

            match7 = re.match(r'Bad format in tree file', line)
            if match7:

                print_to_log = "rate4site_routines.check_if_rate4site_failed : "
                print_to_html += " There is an error in the tree file format. Please check that your tree is in the <a href = \"" + GENERAL_CONSTANTS.CONSURF_TREE_FAQ + "\">requested format</a> and reupload it to the server.<br>"
                error_found = "yes"
                break

            match8 = re.match(r'not all sequences are of the same lengths', line)
            if match8:

                print_to_log = "rate4site_routines.check_if_rate4site_failed : problem with the MSA : not all sequences are of the same lengths"
                error_found = "yes"
                break

            line = R4S_RES.readline()

        R4S_RES.close()

    if error_found == "yes":

        return_ans = "yes"

    return(return_ans, print_to_log, r4s_process_id, print_to_html)

def extract_diversity_matrix_info(r4s_log_file):

    # extracting diversity matrix info

    matrix_disINFO = "\"\""
    matrix_lowINFO = "\"\""
    matrix_upINFO = "\"\""

    try:

        RES_LOG = open(r4s_log_file, 'r')

    except:

        return("rate4site_routines.extract_diversity_matrix_info: Can't open '" + r4s_log_file + "' for reading\n")

    line = RES_LOG.readline()
    while line != "":

        line = line.rstrip()
        match1 = re.match(r'\#Average pairwise distance\s*=\s+(.+)', line)
        if match1:

            matrix_disINFO = match1.group(1)

        else:

            match2 = re.match(r'\#lower bound\s*=\s+(.+)', line)
            if match2:

                matrix_lowINFO = match2.group(1)

            else:

                match3 = re.match(r'\#upper bound\s*=\s+(.+)', line)
                if match3:

                    matrix_upINFO = match3.group(1)
                    break

        line = RES_LOG.readline()

    RES_LOG.close()

    return("OK", matrix_disINFO, matrix_lowINFO, matrix_upINFO)