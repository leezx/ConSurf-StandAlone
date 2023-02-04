#!/usr/bin/env /powerapps/share/python3.8/bin/python

#import cgi
#import cgitb
import sys

#print("Content-type:text/html\n")

#sys.path.append("/bioseq/consurf/python/")
import GENERAL_CONSTANTS
import parseFiles
import BlastResultsForm
import pdbParser
import MSA_parser
import TREE_parser
import wasabiUtils
import rate4site_routines
import cp_rasmol_gradesPE_and_pipe
import ConSeq_gradesPE_and_Outputs



import socket
import getpass
import re
import os
import json
import shutil
import gzip
import requests
import Bio
import subprocess
import time
import tempfile
import urllib
import argparse

from datetime import date
from datetime import datetime
from Bio import AlignIO
from zipfile import ZipFile
from Bio import SeqIO



def get_form():
    
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("--msa", default = None, help = "MSA file name.")
    arg_parser.add_argument("--query", default = "", help = "MSA query sequence name.")
    arg_parser.add_argument("--tree", default = None, help = "Tree file name.")
    arg_parser.add_argument("--pdb", default = None, help = "PDB file name.")
    arg_parser.add_argument("--chain", help = "PDB chain.")
    arg_parser.add_argument("--iterations", default = 1, help = "Number of iterations.")
    arg_parser.add_argument("--cutoff", default = 0.0001, help = "E-value cutoff.")
    arg_parser.add_argument("--DB", default = "UNIREF90", help = "Choose protein data base: SWISS-PROT, CLEAN_UNIPROT, UNIREF90, UniProt or NR_PROT_DB. Nucleotide data base is NT.")
    arg_parser.add_argument("--MAX_HOMOLOGS", default = 150, help = "Maximum number of homologs.")
    arg_parser.add_argument("--closest", action='store_true', help = "Select homologs closest to the query (otherwise sample the list of homologs).")
    arg_parser.add_argument("--MAX_ID", default = 95, help = "Maximal %%ID between sequences. Filter out redundant sequences. Sequences are clustered according to the given sequence identity cutoff, one representative of each cluster is reserved.")
    arg_parser.add_argument("--MIN_ID", default = 35, help = "Minimal %%ID for homologs. Minimal sequence identity with the query sequence. Hits that share less than the given identity cutoff are ignored.")
    arg_parser.add_argument("--align", default = "MAFFT", help = "Choose alignment method: MAFFT, PRANK, MUSCLE  or CLUSTALW.")
    arg_parser.add_argument("--Maximum_Likelihood", action='store_true', help = "The calculation method for the rate of evolution at each site in the MSA. The default is Bayesian. This flag changes it to the Maximum Likelihood method.")
    arg_parser.add_argument("--model", default = "BEST", help = "Choose Evolutionary substitution model: JTT, LG, mtREV, cpREV, WAG or Dayhoff.")
    arg_parser.add_argument("--seq", default = "protein_seq.fas", help = "Sequence file name in fasta format.")
    arg_parser.add_argument("--Nuc", action='store_true', help = "Nucleic acid.")
    arg_parser.add_argument("--algorithm", default = "BLAST", help = "Choose algorithm: HMMER, BLAST or CSI-BLAST.")
    arg_parser.add_argument("--dir", required = True, help = "working directory")

    args = arg_parser.parse_args()
    form['msa_SEQNAME'] = args.query
    form['pdb_ID'] = None
    form['FASTA_text'] = None
    form['modeller_key'] = None
    form['user_select_seq'] = "no"
    form['pdb_FILE'] = args.pdb
    vars['user_msa_file_name'] = args.msa
    form['tree_name'] = args.tree
    form['ITERATIONS'] = args.iterations
    vars['protein_seq'] = args.seq
    vars['working_dir'] = args.dir + "/"
    vars['pdb_file_name'] = args.pdb
    form['proteins_DB'] = args.DB
    form['E_VALUE'] = args.cutoff
    form['Homolog_search_algorithm'] = args.algorithm
    form['MAX_NUM_HOMOL'] = args.MAX_HOMOLOGS
    form['MAX_REDUNDANCY'] = args.MAX_ID
    form['MIN_IDENTITY'] = args.MIN_ID
    form['MSAprogram'] = args.align
    form['PDB_chain'] = args.chain
    form['SUB_MATRIX'] = args.model
    vars['tree_file'] = args.tree
    if args.Nuc:
        
        form['DNA_AA'] = "Nuc"
        form['proteins_DB'] = "NT"
        
    else:
        
        form['DNA_AA'] = "AA"
        
    if form['DNA_AA'] == "Nuc" and form['Homolog_search_algorithm'] != "BLAST":
        
        form['Homolog_search_algorithm'] = "HMMER"
        
    if args.closest:
        
        form['best_uniform_sequences'] = "best"
        
    else:
        
        
        form['best_uniform_sequences'] = "uniform"
        
    if args.Maximum_Likelihood:
        
        form['ALGORITHM'] = "Bayes"
        
    else:
        
        form['ALGORITHM'] = "LikelihoodML"
        
    if vars['pdb_file_name'] is None:
        
        vars['pdb_file_name'] = "pdb_file.ent"
        
    if vars['tree_file'] is None:
        
        vars['tree_file'] = "TheTree.txt"




        
            

        
        


def create_cd_hit_output(input_file, output_file, cutoff, cd_hit_dir, ref_cd_hit_hash, type):


    seq = ""
    seq_name = ""
    cmd = "" 
    n = 0

    # running cd-hit

    if type == "AA":

        cmd += "%scd-hit -i %s -o %s " %(cd_hit_dir, input_file, output_file)
        if cutoff > 0.7 and cutoff < 1:

            n = 5

        elif cutoff > 0.6 and cutoff <= 0.7:

            n = 4

        elif cutoff > 0.5 and cutoff <= 0.6:

            n = 3

        elif cutoff > 0.4 and cutoff <= 0.5:

            n = 2
            
    else:

        # DNA
        cmd += "%scd-hit-est -i %s -o %s " %(cd_hit_dir, input_file, output_file)
        if cutoff > 0.9 and cutoff < 1:

            n = 8

        elif cutoff > 0.88 and cutoff <= 0.9:

            n = 7

        elif cutoff > 0.85 and cutoff <= 0.88:

            n = 6

        elif cutoff > 0.8 and cutoff <= 0.85:

            n = 5

        elif cutoff > 0.75 and cutoff <= 0.8:

            n = 4
            
    cmd += "-c %f -n %d" %(cutoff, n)

    submit_job_to_Q("CD-HIT", cmd)

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

	
	
def send_mail_to_select_sequences():
	
    email_subject = "Your ConSurf run titled " + form['JOB_TITLE'] + " is waiting for you to select sequences."
    email_message = "Hello,\\n\\nYour ConSurf run (number " + form['Run_Number'] + ") is waiting for you to select sequences.\\nPlease have a look at " + vars['run_url'] + "  to select your sequences\\n\\nConSurf Team\\n"
    send_email(email_subject, email_message)
	
	

def make_sequences_file_HTML(plain_txt_sequences, HTML_sequences):

    try:

        HTML_SEQUENCES = open(HTML_sequences, 'w')

    except:

        exit_on_error('sys_error', "make_sequences_file_HTML : cannot open the file " + HTML_sequences + " for writing.")

    try:

        TXT_SEQUENCES = open(plain_txt_sequences, 'r')

    except:

        exit_on_error('sys_error', "make_sequences_file_HTML : cannot open the file " + plain_txt_sequences + " for reading.")
	
    counter = 1	
    line = TXT_SEQUENCES.readline()
    while line != "":

        line = line.rstrip()
        match = re.match(r'^>(Input_pdb.*)', line)
        if match and form['pdb_ID'] is not None:

            line = line.replace(">", "")
            HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3><A HREF=\"https://www.rcsb.org/pdb/explore/explore.do?structureId=%s\">>%d_%s</A></FONT><BR>\n" %(form['pdb_ID'], counter, line.strip()))

        else:

            match = re.match(r'^>(UniRef90_[A-Za-z0-9]+)', line)
            if match:

                line = line.replace(">", "")
                HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3><A HREF=\"https://www.uniprot.org/uniref/%s\">>%d_%s</A></FONT><BR>\n" %(match.group(1), counter, line.strip()))

            else:

                match = re.match(r'>sp\|([A-Za-z0-9]+)', line)
                if match:

                    line = line.replace(">", "")
                    HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3><A HREF=\"https://www.uniprot.org/uniprot/%s\">>%d_%s</A></FONT><BR>\n" %(match.group(1), counter, line.strip()))

                else:
				
                    match = re.match(r'>tr\|([A-Za-z0-9]+)', line)
                    if match:

                        line = line.replace(">", "")
                        HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3><A HREF=\"https://www.uniprot.org/uniprot/%s\">>%d_%s</A></FONT><BR>\n" %(match.group(1), counter, line.strip()))

                    else:
			
                        match = re.match(r'>gi\|([A-Za-z0-9]+)', line)
                        if match:

                            line = line.replace(">", "")
                            HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3><A HREF=\"https://www.ncbi.nlm.nih.gov/nuccore/%s\">>%d_%s</A></FONT><BR>\n" %(match.group(1), counter, line.strip()))

                        else:
		
                            match = re.match(r'>\S+\|([A-Za-z0-9]+)', line)
                            if match:

                                line = line.replace(">", "")
                                HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3><A HREF=\"https://www.uniprot.org/uniprot/%s\">>%d_%s</A></FONT><BR>\n" %(match.group(1), counter, line.strip()))
					
                            else:
				
                                counter += 1		
                                HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3>" + line.strip() + "</FONT><BR>\n")


        line = TXT_SEQUENCES.readline()			
		
    """
    line = TXT_SEQUENCES.readline()
    while line != "":

        line = line.rstrip()
        match1 = re.match(r'^>([A-Z0-9]+)\|(.*)', line)
        if match1:

            HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3>> <A HREF=\"https://www.uniprot.org/uniprot/%s\">%s</A> \| %s </FONT><BR>\n" %(match1.group(1), match1.group(1), match1.group(2)))

        else:

            match2 = re.match(r'^>(Input_pdb.*)', line)
            if match2 and form['pdb_ID'] != "":

                HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3>> <A HREF=\"https://www.rcsb.org/pdb/explore/explore.do?structureId=%s\">%s</A> \| %s </FONT><BR>\n" %(form['pdb_ID'], form['pdb_ID'], match2.group(1)))

            else:

                match3 = re.match(r'^>(UniRef90_[A-Za-z0-9]+)(.*)', line)
                if match3:

                    HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3>> <A HREF=\"https://www.uniprot.org/uniref/%s\">%s</A>%s </FONT><BR>\n" %(match3.group(1), match3.group(1), match3.group(2)))

                else:

                    match4 = re.match(r'^>([A-Za-z0-9]+_[A-Za-z0-9]+)(_[0-9]+_[0-9]+)(.*)', line)
                    if match4:

                        HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3>> <A HREF=\"https://www.uniprot.org/uniprot/%s\">%s</A>%s %s </FONT><BR>\n" %(match4.group(1), match4.group(1), match4.group(2), match4.group(3)))

                    else:

                        re.sub(r'(.{100,100})', r'\1\<BR\>\n', line)
                        line.rstrip() # we make sure there exactly one new line char
                        line += "\n"
                        HTML_SEQUENCES.write("<FONT FACE=\"courier new\" SIZE=3>" + line + "</FONT><BR>\n")

        line = TXT_SEQUENCES.readline()	
    """

    HTML_SEQUENCES.close()
    TXT_SEQUENCES.close()


	
def final_homologoues_choosen_manually(alignments):

    try:

        FINAL = open(vars['FINAL_sequences'], 'w')

    except:

        exit_on_error('sys_error',"final_homologoues_choosen_manually : cannot open the file %s for writing" %vars['FINAL_sequences'])

    try:

        CHOOSEN_SEQS = open("selected_seqs.txt", 'r')

    except:

        exit_on_error('sys_error',"final_homologoues_choosen_manually : cannot open the file selected_seqs.txt for reading")

    choosen = (CHOOSEN_SEQS.read()).split()
    CHOOSEN_SEQS.close()

    vars['final_number_of_homologoues'] = len(choosen) + 1
    vars['unique_seqs'] = vars['final_number_of_homologoues']

    FINAL.write(">%s\n%s\n" %(vars['query_string'], vars['protein_seq_string']))
	
    for number in choosen:
	
        FINAL.write(">%s\n%s\n" %(alignments[int(number)]['name'], alignments[int(number)]['seq']))

    FINAL.close()



def Create_HTML_file(choose_file, parsed_output, raw_output, query, job_number):

    if form['DNA_AA'] == "AA":

        database = GENERAL_CONSTANTS.UNIPOROT_SITE

    else:

        database = GENERAL_CONSTANTS.NCBI_SITE

    if form['Homolog_search_algorithm'] == "HMMER":
	
        LOG.write("BlastResultsForm.Create_HTML_file_for_hmmer(%s, %s, %s, %s, %s, %s)\n" %(choose_file, parsed_output, raw_output, query, job_number, database))
        ans = BlastResultsForm.Create_HTML_file_for_hmmer(choose_file, parsed_output, raw_output, query, job_number, database)
		
    else:
	
        LOG.write("BlastResultsForm.Create_HTML_file_for_blast(%s, %s, %s, %s, %s, %s)\n" %(choose_file, parsed_output, raw_output, query, job_number, database))
        ans = BlastResultsForm.Create_HTML_file_for_blast(choose_file, parsed_output, raw_output, query, job_number, database)	
	  
    vars['number_of_homologoues'] = ans[2]
    if ans[0] != "OK":
    
        exit_on_error('sys_error', ans[1])

    return ans[1]
    
   

def downloadMmtf():

    mmtfFileName = (form['pdb_ID']).lower() + ".mmtf.gz"
    http = "http://mmtf.rcsb.org/v1.0/full/"

    try:

        response = urllib.request.urlopen(http + mmtfFileName, timeout=10).read()

    except:

        LOG.write("downloadMmtf : the file " + mmtfFileName + " was not downloaded.\n")
        return None

    try:

        GET_Mmtf = open(mmtfFileName, 'wb')

    except:

        exit_on_error('system_error', "downloadMmtf : can't open the file " + mmtfFileName + " for writing.")

    GET_Mmtf.write(response)
    GET_Mmtf.close()
	
    return mmtfFileName



def Update_ConSurf_Users_Log():

    user_ip = socket.gethostbyname(socket.gethostname())
    curr_time = time.asctime(time.localtime(time.time()))
    pid = os.getpid()
    try:

        LIST = open(GENERAL_CONSTANTS.CONSURF_LOG, 'a')

    except:

        exit_on_error('sys_error', "Update_ConSurf_Users_Log : Can't open the file " + GENERAL_CONSTANTS.CONSURF_LOG + " for writing.")

    LIST.write(curr_time + " " + form['Run_Number'] + " " + vars['running_mode'] + " " + user_ip + "  PID: " + str(pid) + " email: " + form['user_email'] + "\n")
    LIST.close()

    LOG.write("Program: " + sys.argv[0] + " time: " + curr_time + " " + form['Run_Number'] + " " + vars['running_mode'] + " " + user_ip + "  PID: " + str(pid) + " email: " + form['user_email'] + "\n")


def Failed():

    PROGRESS = open(vars['PROGRESS_REPORT'], 'r')
    lines = PROGRESS.readlines()
    PROGRESS.close()

    # delete the old file
    PROGRESS = open(vars['PROGRESS_REPORT'], 'w')
    PROGRESS.close()

    for i in range(len(lines)):

        if i == vars['progress_stage']:

            print_progress(lines[i], "x")

        else:

            print_progress(lines[i], "")

    try:

        DONE = open("FAILED", 'w')

    except:

        exit_on_error('sys_error', "Failed : can't open the file FAILED for writing.")




def Update_Progress(failed = ""):

    PROGRESS = open(vars['PROGRESS_REPORT'], 'r')
    lines = PROGRESS.readlines()
    PROGRESS.close()

    # delete the old file
    PROGRESS = open(vars['PROGRESS_REPORT'], 'w')
    PROGRESS.close()

    for i in range(len(lines)):

        if i == vars['progress_stage']:

            if failed == "":

                print_progress(lines[i], "v")

            else:

                print_progress(lines[i], "x")
          
        else:

            print_progress(lines[i], "")


    TIME_FILE = open(vars['time_table'][vars['progress_stage']], 'a')
    new_current_time = time.time()
    TIME_FILE.write("%s,%.3f\n" %(vars['Used_PDB_Name'], new_current_time - vars['current_time']))
    vars['current_time'] = new_current_time
    TIME_FILE.close()

    vars['progress_stage'] += 1

def Print_Initial_Running_Progress():

    vars['progress_stage'] = 0
    vars['stages'] = "1"
    if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_pdb_msa_tree":

        # PDB_mode
        vars['stages'] += "1"
        print_progress("Extract sequence from PDB file", "dot")

    else:

        vars['stages'] += "0"

    if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_no_pdb_no_msa":

        # MSA was not provided
        vars['stages'] += "1"
        print_progress("Find sequence homologs", "dot")
        if form['user_select_seq'] is not None:

            vars['stages'] += "1"
            print_progress("Remove redundant and unrelated sequences", "dot")

        else:

            vars['stages'] += "0"

        print_progress("Align sequences", "dot")

    else:

        vars['stages'] += "0"


    if form['SUB_MATRIX'] == "BEST":

        vars['stages'] += "1"
        print_progress("Select best evolutionary model", "dot")

    else:

        vars['stages'] += "0"

    print_progress("Calculate conservation scores", "dot")

    if vars['running_mode'] == "_mode_msa" or vars['running_mode'] == "_mode_no_pdb_no_msa" or vars['running_mode'] == "_mode_msa_tree":

        # ConSeq Mode
        vars['stages'] += "1"
        if form['DNA_AA'] == "AA":

            vars['stages'] += "1"
            print_progress("Search for 3D structure for the protein sequence", "dot")
            if form['modeller_key'] is not None:

                vars['stages'] += "1"
                print_progress("Predict 3D structure using HHPred and MODELLER", "dot")

            else:

                vars['stages'] += "0"

        else:

            vars['stages'] += "0"

    else:

        vars['stages'] += "0"

    print_progress("Project conservation scores onto the molecule", "dot")

def print_progress(text, symbol):

    vars['PROGRESS_REPORT'] = "ProgressReport.php"
    #LOG.write("starting  :%s\n" %vars['PROGRESS_REPORT'])
    try:

        PROGRESS = open(vars['PROGRESS_REPORT'], 'a')

    except:

        exit_on_error('sys_error', "print_progress : could not open the file %s for writing" %vars['PROGRESS_REPORT'])

    if symbol == "v":

        PROGRESS.write(text.replace("dotsymbol", "vsymbol"))

    elif symbol == "x":

        PROGRESS.write(text.replace("dotsymbol", "xsymbol"))

    elif symbol == "dot":

        PROGRESS.write("<div class=\"box\" style=\"margin-bottom: 35px;float: left;\"><?include(\"dotsymbol.php\");?><div class=\"box\" style=\"margin-left: 35px;float: left;\">%s</div></div>\n" %text)

    else:

        PROGRESS.write(text)

    PROGRESS.close()

    if symbol == "dot":
	
        vars['time_table'].append("/bioseq/consurf/time-tables/" + text +".txt")

def create_output_php():

    try:
	
        JOB_DETAILS = open("job_details.txt", 'w')
		
    except:
	
        exit_on_error('sys_error', "could not open the file job_details.txt for writing.")

    try:

        STRUCTURE = open("structure.php", 'w')

    except:

        exit_on_error('sys_error', "could not open the file structure.php for writing.")


    try:

        PHYLOGENETIC = open("phylogenetic.php", 'w')

    except:

        exit_on_error('sys_error', "could not open the file phylogenetic.php for writing.")


    try:

        SCORE = open("score.php", 'w')

    except:

        exit_on_error('sys_error', "could not open the file score.php for writing.")


    try:

        ALIGNMENT = open("alignment.php", 'w')

    except:

        exit_on_error('sys_error', "could not open the file alignment.php for writing.")


    SCORE.write("Method of Calculation: ")
    if form['ALGORITHM'] == "Bayes":

        SCORE.write("Bayesian")

    elif form['ALGORITHM'] == "LikelihoodML":

        SCORE.write("Max. Likelihood")

    if vars['running_mode'] == "_mode_pdb_msa_tree" or vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_pdb_no_msa":

        if form['pdb_FILE'] is not None:

            STRUCTURE.write("PDB File: " + form['pdb_FILE'])

        else:

            STRUCTURE.write("PDB ID: " + form['pdb_ID'])

        STRUCTURE.write("<br><br>Chain identifier: " + form['PDB_chain'])

    elif vars['running_mode'] == '_mode_msa' or vars['running_mode'] == '_mode_msa_tree':

        ALIGNMENT.write("MSA File: " + form['msa_input'])

        if form['AA_or_NUC'] == "AA":

            STRUCTURE.write("Not Given. Running in 'ConSeq' mode")

        else:

            STRUCTURE.write("Nucleic Acids Sequnce")

    elif vars['running_mode'] == '_mode_no_pdb_no_msa':

        STRUCTURE.write("Not Given. Running in 'ConSeq' mode")


    if vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_pdb_msa_tree" or vars['running_mode'] == "_mode_msa" or vars['running_mode'] == "_mode_msa_tree":

        # include MSA
        ALIGNMENT.write("MSA File: %s<br><br>Query sequence name in MSA file: %s\"" %(form['msa_input'], form['msa_SEQNAME']))

    else:

        # no MSA
        ALIGNMENT.write("Multiple Sequence Alignment is built using %s<br><br>The Homologues are collected from %s<br><br>Homolog search algorithm: %s<br><br>%s E-value: %s<br><br>No. of %s Iterations: %s" %(form['MSAprogram'], form['proteins_DB'], form['Homolog_search_algorithm'], vars['blast_algorithm'], form['E_VALUE'], vars['blast_algorithm'], form['ITERATIONS']))
        if form['user_select_seq'] != "yes":

            if form['best_uniform_sequences'] == "best":

                ALIGNMENT.write("<br><br>Maximal %%ID Between Sequences : %s<br><br>Minimal %%ID For Homologs : %s<br><br>%s sequences closest to the query." %(form['MAX_REDUNDANCY'], form['MIN_IDENTITY'], form['MAX_NUM_HOMOL']))

            else:

                ALIGNMENT.write("<br><br>Maximal %%ID Between Sequences : %s<br><br>Minimal %%ID For Homologs : %s<br><br>%s sequences that sample the list of homologues to the query." %(form['MAX_REDUNDANCY'], form['MIN_IDENTITY'], form['MAX_NUM_HOMOL']))

        else:

            ALIGNMENT.write("<br><br>Sequences are selected by user out of %s resaults" %vars['blast_algorithm'])

        matrix = form['SUB_MATRIX']
        if matrix == "BEST":

            matrix = "Best fit"

        SCORE.write("<br><br>Model of substitution for %s: %s" %(vars['protein_or_nucleotide'], matrix))

    if vars['running_mode'] == "_mode_pdb_msa_tree" or vars['running_mode'] == "_mode_msa_tree":

        # include tree
        PHYLOGENETIC.write("TREE file: " + form['tree_name'][1])

    else:

        # no tree
        PHYLOGENETIC.write("Neighbor Joining with ML distance")

    STRUCTURE.close()
    PHYLOGENETIC.close()
    SCORE.close()
    ALIGNMENT.close()

    Print_Initial_Running_Progress()

    #today = date.today()
    JOB_DETAILS.write("%s\n%s\n%s\n%s" %(vars['Used_PDB_Name'], form['JOB_TITLE'], vars['stages'], vars['date']))
    JOB_DETAILS.close()







def print_nucleotide_precentage():

    # print a file that details percentage of each nucleotide in the MSA
    LOG.write("cp_rasmol_gradesPE_and_pipe.print_precentage_nuc(%s, %s, %s)\n" %("residue_freq", "position_totalAA", vars['Msa_percentageFILE']))
    ans = cp_rasmol_gradesPE_and_pipe.print_precentage_nuc(residue_freq, position_totalAA, vars['Msa_percentageFILE'])
    if ans != "OK":

        exit_on_error('sys_error',str(ans))

    elif not os.path.exists(vars['Msa_percentageFILE']) or os.path.getsize(vars['Msa_percentageFILE']) == 0:

        exit_on_error('sys_error', "The output " + vars['Msa_percentageFILE'] + " was not found or empty")



def convert_rna_to_dna(Seqs, Seqs_dna):

    # replace the u with t and return the sequences names replaced

    try:

        OUT = open(Seqs_dna, 'w')

    except:

        return("convert_rna_to_dna: Can't open file " + Seqs_dna + " for writing.")

    try:

        SEQS = open(Seqs, 'r')

    except:

        return("convert_rna_to_dna: Can't open file " + Seqs + " for reading.")

    Seqs_Names = []
    seq_name = ""

    line = SEQS.readline()
    while line != "":

        line = line.rstrip()
        match1 = re.match(r'^>(.*)', line)
        if match1:

            seq_name = match1.group(1)


        elif 'u' in line or 'U' in line:

            Seqs_Names.append(seq_name)
            line = line.replace('u', 't')
            line = line.replace('U', 'T')

        OUT.write(line + "\n")
        line = SEQS.readline()

    OUT.close()
    SEQS.close()

    return("OK", Seqs_Names)




def pisa_html():

    try:

        PISA_HTML = open("pisa.html", 'w')

    except:

        exit_on_error('sys_error', "pisa_html : cannot open the file pisa.html for writing.")

    PISA_HTML.write("""
                                <div class="collapse">
                                    <a href="javascript:void(0)" class="accordion">
                                        <div class="plus-minus-toggle"></div><h3>PISA PDB file colored by ConSurf</h3>
                                    </a>
                                    <div class="panel"><br>
                                        <ul>
                                           <li><a href="<?=$ngl_viewer_pisa?>" target="_blank">View ConSurf results projected on PISA first model</a> with NGL viewer</li>
                                           <li><a href="<?=$orig_path?>/<?=$name?>_pisa_consurf_grades.txt"  class="download" onclick="gtag('event', 'download', { 'event_category' : 'download_pymol', 'event_action' : 'download_pymol'});" download>Amino Acid Conservation Scores, Confidence Intervals and Conservation Colors</a></li>
                                           <li><a href="<?=$orig_path?>/<?=$name?>_pisa_With_Conservation_Scores.pdb" class="download" onclick="gtag('event', 'download', { 'event_category' : 'download_pymol', 'event_action' : 'download_pymol'});" download>PDB File with Conservation Scores in the tempFactor field</a></li>
                                           <li><a href="<?="figure_instructions.php?number=$number&name=$name&cbs=$cbs&pisa=yes"?>" target="_blank">Follow the instructions to produce a high resolution figure</a></li>
                                        </ul>
                                    </div>
                                </div>
                    """)
    PISA_HTML.close()

def more_parameters_html():

    try:

        FILES = open("files.html", 'a')

    except:

        exit_on_error('sys_error', "more_parameters_html : cannot open the file files.html for writing.")

    LOG.write("parameters_html: rate4site_routines.extract_diversity_matrix_info(%s)\n" %(vars['r4s_log']))
    ans = rate4site_routines.extract_diversity_matrix_info(vars['r4s_log'])
    if ans[0] == "OK":   
		
        FILES.write("""
					<li><b>Alignment details</b></li>
					<ul>
					   <li>The average number of replacements between any two sequences in the alignment;<br>A distance of 0.01 means that on average, the expected replacement for every 100 positions is 1.</li>
					   <li>Average pairwise distance : %s</li>
                                           <li>Lower bound : %s</li>
					   <li>Upper bound : %s</li>
                                           <li><a href="<?=$orig_path?>/%s" onclick="gtag('event', 'download', { 'event_category' : 'download_jackhmmer', 'event_action' : 'download_jackhmmer'});" download>Residue variety per position in the MSA</a> (The table is best viewed with an editor that respects Comma-Separated Values)</li>											   
				       </ul>
	""" %(ans[1], ans[2], ans[3], vars['Msa_percentageFILE']))


    FILES.write("""
                                         <ul>
                                            <li><a href="%swasabi/?url=%sresults/<?=$number?>/%s" target="WASABI"> View MSA and phylogenetic tree using WASABI</a></li>
                                            <li><a href="<?=$orig_path?>/%s" onclick="gtag('event', 'download', { 'event_category' : 'download_tree', 'event_action' : 'download_tree'});" download>Download Phylogenetic Tree </a> (Newick format)</li>
    """ %(GENERAL_CONSTANTS.CONSURF_URL, GENERAL_CONSTANTS.CONSURF_URL, vars['WASABI_XML'], vars['tree_file']))
	
    if vars['best_fit']:
	
        if form['SUB_MATRIX'] == "JC_Nuc":

            best_matrix = "JC"

        else:

            best_matrix = form['SUB_MATRIX']
			
        FILES.write("""                                            <li>The best evolutionary model was selected to be: %s. See details <a href="<?=$orig_path?>/model_selection.txt" target="_blank">here</a></li>\n""" %best_matrix)
	
    FILES.write("""
	                                    </ul>
                                         </ul>
                                    </div>
                                </div>
    """)
    FILES.close()


def parameters_html(msa, tree):

    try:

        FILES = open("files.html", 'w')

    except:

        exit_on_error('sys_error', "parameters_html : cannot open the file files.html for writing.")

    if not msa:

        if form['best_uniform_sequences'] == "best":

            sample_or_best = "the unique hits colsest to the query were chosen"

        else:

            sample_or_best = "sampled from the unique hits"

        FILES.write("""
                                <div class="collapse">
                                    <a href="javascript:void(0)" class="accordion">
                                        <div class="plus-minus-toggle"></div><h3>Homologues, Alignment and Phylogeny</h3>
                                    </a>
                                    <div class="panel"><br> 
                                      <ul>
                                        <ul style="width: 90%%;">
                                            <li><a href="<?=$orig_path?>/%s" onclick="gtag('event', 'download', { 'event_category' : 'download_jackhmmer', 'event_action' : 'download_jackhmmer'});" download>%d homologues</a> were collected from the %s database using %s.</li>
        """ %(vars['BLAST_out_file'], vars['number_of_homologoues'], form['proteins_DB'], form['Homolog_search_algorithm']))

        if form['user_select_seq'] != "yes":

            FILES.write("""
                                            <li>Of these, <a href="<?=$orig_path?>/%s" onclick="gtag('event', 'download', { 'event_category' : 'download_before_cdhit', 'event_action' : 'download_before_cdhit'});">%d homologues</a> passed the thresholds (min/max similarity, coverage, etc), <a href="<?=$orig_path?>/%s" onclick="gtag('event', 'download', { 'event_category' : 'download_cdhit', 'event_action' : 'download_cdhit'});">%d of them</a> are CD-HIT unique.</li>
                                            <li>The calculations were conducted on <a href="<?=$orig_path?>/%s" target="_blank">%d hits</a> (query included), %s. Click <a href="<?=$orig_path?>/%s" onclick="gtag('event', 'download', { 'event_category' : 'download_rejected', 'event_action' : 'download_rejected'});">here</a> if you wish to view the list of sequences which produced significant alignments, but were not chosen as hits.</li>
                                        </ul>
            """ %(vars['HITS_fasta_file'], vars['number_of_homologoues_before_cd-hit'], vars['cd_hit_out_file'], vars['unique_seqs'], vars['FINAL_sequences_html'], vars['final_number_of_homologoues'], sample_or_best, vars['HITS_rejected_file']))

        else:

            FILES.write("""
                                            <li>The calculations were conducted on <a href="<?=$orig_path?>/%s" target="_blank">%d hits</a> (query included), selected manually.</li>
                                        </ul>
            """ %(vars['FINAL_sequences_html'], vars['final_number_of_homologoues']))


    else:

        FILES.write("""
                                <div class="collapse">
                                    <a href="javascript:void(0)" class="accordion">
                                        <div class="plus-minus-toggle"></div><h3>Alignment and Phylogeny</h3>
                                    </a>
                                    <div class="panel"><br>
                                        <ul>
        """)

    FILES.close()

    try:

        PARAMETERS = open("parameters.html", 'w')

    except:

        exit_on_error('sys_error', "parameters_HTML : cannot open the file parameters.html for writing.")

    if form['best_uniform_sequences'] == "best":

        how_homologues_are_chosen = " These are sampled from the list of unique homologues."

    else:

        how_homologues_are_chosen = ""

    if not msa:

        PARAMETERS.write("""
                                <div class="collapse">
                                    <a href="javascript:void(0)" class="accordion">
                                        <div class="plus-minus-toggle"></div><h3>Running Parameters</h3>
                                    </a>
                                    <div class="panel"><br>
                                        <ul>
                                            <li><b>Homologues Search:</b></li>
                                            <ul>
                                                <li>Homologues were collected from %s database.</li>
                                                <li>Homologues search algorithm is %s.</li>
                                                <li>E-value cutoff is %s.</li>
                                                <li>Number of Iterations is %s.</li>
                                            </ul>
        """ %(form['proteins_DB'], form['Homolog_search_algorithm'], form['E_VALUE'], form['ITERATIONS']))

        if form['user_select_seq'] != "yes":

            PARAMETERS.write("""
                                            <li><b>Homologues Thresholds:</b></li>
                                            <ul style="width: 90%%;">
                                                <li>CD-HIT cutoff is %d%% (This is the maximal sequence identity between homologues).</li>
                                                <li>Maximal number of final homologues is %s.%s</li>
                                                <li>Maximal overlap between homologues is 10%% (If overlap between two homologues exceeds 10%%, the highest scoring homologue is chosen).</li>
                                                <li>Coverage is 60%% (This is the minimal percentage of the query sequence covered by the homologue).</li>
                                                <li>Minimal sequence identity with the query sequence is %s%%.</li>
                                            </ul>
            """ %(vars['hit_redundancy'], form['MAX_NUM_HOMOL'], how_homologues_are_chosen, form['MIN_IDENTITY']))

        PARAMETERS.write("""
                                            <li><b>Alignment, Phylogeny and Conservation Scores:</b></li>
                                            <ul>
                                                <li>Multiple Sequence Alignment was built using %s.</li>
        """ %(form['MSAprogram']))

    else:

        PARAMETERS.write("""
                                <div class="collapse">
                                    <a href="javascript:void(0)" class="accordion">
                                        <div class="plus-minus-toggle"></div><h3>Running Parameters</h3>
                                    </a>
                                    <div class="panel"><br>
                                        <ul>
        """)



    if not tree:

        PARAMETERS.write("""                                                <li>Phylogenetic tree was built using Neighbor Joining with ML distance.</li>\n""")

    if form['ALGORITHM'] == "Bayes":
	
        algorithm = "Bayesian"
		
    else:
	
        algorithm = "Max. Likelihood (ML)"
	
    if form['SUB_MATRIX'] == "BEST":
	
        matrix = "chosen by best fit"
		
    else:
	
        matrix = form['SUB_MATRIX']


    PARAMETERS.write("""
                                            <li>Conservation Scores were calculated with the %s method.</li>
                                            <li>Amino acid substitution model was %s.</li>
                                        </ul>
                                    </ul>
                                </div>
                            </div>
    """ %(algorithm, matrix))



    PARAMETERS.close()



def create_pdf(color_column, B_E_column, prediction_method):

    url = "https://consurf.tau.ac.il/barak/colored_sequence.php?number=" + form['Run_Number'] + "&cbs=0&color_column=" + color_column + "&B_E_column=" + B_E_column + "&prediction_method=" + prediction_method + "&date=" + vars['date'] + "&title=" + form['JOB_TITLE']
    url_CBS = "https://consurf.tau.ac.il/barak/colored_sequence.php?number=" + form['Run_Number'] + "&cbs=1&color_column=" + color_column + "&B_E_column=" + B_E_column + "&prediction_method=" + prediction_method + "&date=" + vars['date'] + "&title=" + form['JOB_TITLE']
    
    try:
    
        response = urllib.request.urlopen(url, timeout=10).read()
        
    except:
    
        exit_on_error('sys_error', "create_pdf : the url " + url + " can't be reached.")

    try:

        PDF = open(vars['Colored_Seq_PDF'], 'wb')

    except:

        exit_on_error('sys_error', "create_pdf : cannot open the file " + vars['Colored_Seq_PDF'] + " for writing!")

    PDF.write(response)
    PDF.close()

    try:
    
        response_CBS = urllib.request.urlopen(url_CBS, timeout=10).read()
        
    except:
    
        exit_on_error('sys_error', "create_pdf : the url " + url_CBS + " can't be reached.")

    try:

        PDF_CBS = open(vars['Colored_Seq_CBS_PDF'], 'wb')

    except:

        exit_on_error('sys_error', "create_pdf : cannot open the file " + vars['Colored_Seq_CBS_PDF'] + " for writing!")

    PDF_CBS.write(response_CBS)
    PDF_CBS.close()






def Create_Colored_MSA():

    msa_html_title = "ConSurf Color-Coded MSA for Job:%s Date:%s" %(form['JOB_TITLE'], vars['date'])
    msa_colored_css = "https://consurfdb.tau.ac.il/scripts/css/MSA_Colored.css"
    LOG.write("ConSeq_gradesPE_and_Outputs.print_msa_colors_FASTA_clustalwLike(%s, %s, %s, %s, %s, %s, %s)\n" %("gradesPE_Output", vars['msa_fasta'], vars['msa_SEQNAME'], vars['Colored_MSA_HTML'], msa_colored_css, msa_html_title, form['DNA_AA']))
    ConSeq_gradesPE_and_Outputs.print_msa_colors_FASTA_clustalwLike(gradesPE_Output, vars['msa_fasta'], vars['msa_SEQNAME'], vars['Colored_MSA_HTML'], msa_colored_css, msa_html_title, form['DNA_AA'])

    msa_colored_css_cbs = "https://consurfdb.tau.ac.il/scripts/css/MSA_Colored_CBS.css"
    LOG.write("ConSeq_gradesPE_and_Outputs.print_msa_colors_FASTA_clustalwLike(%s, %s, %s, %s, %s, %s, %s)\n" %("gradesPE_Output", vars['msa_fasta'], vars['msa_SEQNAME'], vars['Colored_MSA_HTML_CBS'], msa_colored_css_cbs, msa_html_title, form['DNA_AA']))
    ConSeq_gradesPE_and_Outputs.print_msa_colors_FASTA_clustalwLike(gradesPE_Output, vars['msa_fasta'], vars['msa_SEQNAME'], vars['Colored_MSA_HTML_CBS'], msa_colored_css_cbs, msa_html_title, form['DNA_AA'])



def create_pymol(input, prefix):

    cmd = "module load miniconda3-4.5.12-AND-perl-5.28.1-AND-chimera.1.13\n"
    cmd += "pymol -qc " + input + " -d \"run " + vars['pymol_color_script_isd'] + "\"\n"
    cmd += "pymol -qc " + input + " -d \"run " + vars['pymol_color_script_CBS_isd'] + "\"\n"

    LOG.write("create_pymol : %s\n" %cmd)
    submit_job_to_Q("PYMOL", cmd)

    pymol_session = "consurf_pymol_session.pse"
    pymol_session_CBS = "consurf_CBS_pymol_session.pse"

    while not os.path.exists(pymol_session) or os.path.getsize(pymol_session) == 0:

        time.sleep(1)

    while not os.path.exists(pymol_session_CBS) or os.path.getsize(pymol_session_CBS) == 0:

        time.sleep(1)

    os.chmod(pymol_session, 0o664)
    os.chmod(pymol_session_CBS, 0o664)
    os.rename(pymol_session, prefix + pymol_session)
    os.rename(pymol_session_CBS, prefix + pymol_session_CBS)
    vars['zip_list'].append(prefix + pymol_session)
    vars['zip_list'].append(prefix + pymol_session_CBS)



def redirect():

    try:

        OUTPUT = open(vars['output_page'], 'w')

    except:

        exit_on_error("sys_error", "redirect: Could not open the file" + vars['output_page'] + " for writing.")

    OUTPUT.write("""
    <!DOCTYPE html>
        <html>
            <head>
                <title>HTML Meta Tag</title>
                <meta http-equiv = \"refresh\" content = \"3; url = https://consurf.tau.ac.il/barak/final_output.php?number=%s&name=%s&job_name=%s\" />
            </head>
            <body>
                <p>Redirecting to another URL</p>
            </body>
        </html>""" %(form['Run_Number'], vars['Used_PDB_Name'], form['JOB_TITLE']))

    OUTPUT.close()


def create_chimera(input, prefix):

    cmd = "module load miniconda3-4.5.12-AND-perl-5.28.1-AND-chimera.1.13\nchimera --nogui --script '%s %s'\n" %(vars['chimera_color_script'], input)
    cmd += "chimera --nogui --script '%s %s'" %(vars['chimera_color_script_CBS'], input)
    LOG.write("create_chimera : %s\n" %cmd)
    submit_job_to_Q("CHIMERA", cmd)

    cimera_python = "consurf_chimera_session.py"
    cimera_python_CBS = "consurf_CBS_chimera_session.py"

    while not os.path.exists(cimera_python) or os.path.getsize(cimera_python) == 0:

        time.sleep(1)

    while not os.path.exists(cimera_python_CBS) or os.path.getsize(cimera_python_CBS) == 0:

        time.sleep(1)

    os.chmod(cimera_python, 0o664)
    os.chmod(cimera_python_CBS, 0o664)
    os.rename(cimera_python, prefix + cimera_python)
    os.rename(cimera_python_CBS, prefix + cimera_python_CBS)
    vars['zip_list'].append(prefix + cimera_python)
    vars['zip_list'].append(prefix + cimera_python_CBS)





def extract_round_from_hmmer():

    last_round_number = 0
    LOG.write("parseFiles.parse_hmmer(\"%s%s%s\", %s, \"%s%s%s\" );\n" %(vars['working_dir'], dir, vars['BLAST_out_file'], form['ITERATIONS'], vars['working_dir'], dir, vars['BLAST_last_round']))
    ans = parseFiles.parse_hmmer(vars['BLAST_out_file'], form['ITERATIONS'], vars['BLAST_last_round'])

    if ans[0] == "err":

        exit_on_error('sys_error', "parse_hmmer :" + ans[1])

    elif ans[0] == "no_hits":

        err = "No HMMER hits were found.  You may try to:<OL>"

        if form['database'] == "SWISS-PROT":

            err += "<LI>Run your query using UniProt or NR database."

        err += " <li> Increase the Evalue.</li>"
        err += " <li> Decrease the Minimal %ID For Homologs.</li>"

        if form['ITERATIONS']:

            err += " <li> Increase the number of " + vars['blast_algorithm'] + " iterations.</li>"

        err += "</OL>"
        exit_on_error('user_error',err)

    else:

        last_round_number = ans[1]
        LOG.write("last_round_number: " + last_round_number + "\n")

    LOG.write("cp " + vars['BLAST_out_file'] + " " + vars['HMMER_out_file'])
    shutil.copyfile(vars['BLAST_out_file'],  vars['HMMER_out_file'])

    if os.path.exists(vars['BLAST_last_round']) and os.path.size(vars['BLAST_last_round']) != 0:

        shutil.move(vars['BLAST_last_round'], vars['BLAST_out_file'])

    else:

        LOG.write("extract_round_from_hmmer : the file " + vars['BLAST_last_round'] + " was not created. The hits will be collected from the original blast file")
        exit_on_error('sys_error', "extract_round_from_hmmer : the file " + vars['BLAST_last_round'] + " was not created. The hits will be collected from the original blast file")

    return last_round_number




def Get_PISA_Complex(PDB_ID, out_file):

    cmd = "wget https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimer.pdb?%s:1,1 -O %s"  %(PDB_ID.lower(), out_file)
    LOG.write("Get the first model of PISA: %s\n" %str(cmd))
    submit_job_to_Q("PISA", cmd)

    # see if there is a astructure
    if not os.path.exists(out_file):

        return("NA")

    try:

        PISA = open(out_file, 'r')

    except:

        exit_on_error ("sys_error", "Get_PISA_Complex: Could not open PISA: " + out_file)

    line = PISA.readline()
    while line != "":

        if re.match(r'Assembly structure .+ unavailable', line):

            LOG.write(line)
            PISA.close()
            return("NA")

        if re.match(r'\*\*\* No symmetry operations', line):

            LOG.write(line)
            PISA.close()
            return("NA")

        line = PISA.readline()

    PISA.close()
    return("OK")




def read_Solv_Acc_Pred(ref_Solv_Acc_Pred):

    # read the prediction file of PACC

    try:

        PRED = open(vars['Solv_ACC_Pred'], 'r')

    except:

        exit_on_error('sys_error', "read_Solv_Acc_Pred: Can\'t open the file " + vars['Solv_ACC_Pred'])

    trantab = str.maketrans("EB", "eb")
    i = 1
    line = PRED.readline()
    while line != "":

        match = re.match(r'(E|B)\s+\S+\s+\S+', line)
        if match:

            ref_Solv_Acc_Pred[i] = (match.group(1)).translate(trantab)
            i += 1

        line = PRED.readline()

    PRED.close()



def fill_r4s2pdb(ref_r4s2pdb):

    try:

        POSITIONS = open(vars['atom_positionFILE'], 'r')

    except:

        exit_on_error('sys_error', "fill_r4s2pdb: Can't open " + vars['atom_positionFILE'])

    LENGTH_OF_ATOM = 0

    line = POSITIONS.readline()
    while line != "":

        frags = line.split()
        Residue_Name = frags[0]
        FAS_Pos = int(frags[1])
        ATOM_pos = frags[2]
        ref_r4s2pdb[FAS_Pos] ="%s:%s:%s" %(Residue_Name, ATOM_pos ,form['PDB_chain'])
        LENGTH_OF_ATOM += 1
        line = POSITIONS.readline()

    POSITIONS.close()



def check_msa_tree_match(ref_msa_seqs, ref_tree_nodes):

    err_msg = "<br />Note that the search is case-sensitive!<br />Please correct your files and re-upload your query to the server.<br />\n"
    LOG.write("check_msa_tree_match : check if all the nodes in the tree are also in the MSA\n")

    for node in ref_tree_nodes:

        if not node in ref_msa_seqs:

            exit_on_error('user_error', "The uploaded <a href=\"" + vars['tree_file'] + "\">TREE file</a> is inconsistant with the uploaded <a href=\"" + vars['user_msa_file_name'] + "\">MSA file</a>.<br />The node '" + node + "' is found in the TREE file, but there is no sequence in the MSA file with that exact name. " + err_msg)

    LOG.write("check_msa_tree_match : check if all the sequences in the MSA are also in the tree\n")

    for seq_name in ref_msa_seqs: #check that all the msa nodes are in the tree

        if not seq_name in ref_tree_nodes:

            exit_on_error('user_error', "The uploaded <a href=\"" + vars['user_msa_file_name'] + "\">MSA file</a> is inconsistant with the uploaded <a href=\"" + vars['tree_file'] + "\">TREE file</a>.<br />The Sequence name '" + seq_name + "' is found in the MSA file, but there is no node with that exact name in the TREE file. " + err_msg)

    vars['unique_seqs'] = len(ref_msa_seqs)
    LOG.write("There are " + str(vars['unique_seqs']) + " On the MSA\n")


def extract_nodes_from_tree(ref_to_tree_nodes):

    try:

        TREEFILE = open(vars['tree_file'], 'r')

    except:

        exit_on_error('sys_error', "extract_nodes_from_tree : could not open" + vars['tree_file'])

    tree = ""

    line = TREEFILE.readline()

    while line != "":

        line.rstrip()
        if re.search(r'\S', line):

            tree += line

        line = TREEFILE.readline()

    TREEFILE.close()

    LOG.write("extract_nodes_from_tree : calling TREE_parser.extract_nodes_from_tree()\n")
    ans = TREE_parser.extract_nodes_from_tree(tree, ref_to_tree_nodes)
    if ans[0] != "OK":

        match = re.search(r'^duplicity: (.*)$', ans[1])
        if match:

            exit_on_error('user_error', "The uploaded <a href=\"" + vars['tree_file'] + "\">TREE file</a>, which appears to be in Newick format, contains the same node identifier: '" + match.group(1) + "' more than once.<br />Note that each node in the tree should have a uniuqe identifier. Please correct your file and re-upload your query to the server.<br />\n")



def check_validity_tree_file():

    ERR = {}
    LOG.write("check_validity_tree_file : calling TREE_parser.check_validity_tree_file(" + vars['tree_file']  + ")\n")
    TREE_parser.check_validity_tree_file(vars['tree_file'], ERR)
    if 'left_right' in ERR or 'noRegularFormatChar' in ERR:

        msg = "The uploaded <a href=\"" + form['tree_name'][1] + "\">TREE file</a>, which appears to be in Newick format, "
        if 'left_right' in ERR:

            msg += "is missing parentheses."

        elif 'noRegularFormatChar' in ERR:

            msg += "contains the following non-standard characters: %s" %ERR['noRegularFormatChar']

        msg += "<br />\nPlease fix the file and re-run your query.<br />"
        exit_on_error('user_error', msg)

    LOG.write("check_validity_tree_file : tree is valid\n")


def upload_tree():

    try:

        TREE = open(vars['tree_file'], 'bw')

    except:

        return("upload_tree: Can't open file " + vars['tree_file'] + " for writing.")

    TREE.write(form['tree_FILE'])
    TREE.close()


def remove_core(r4s_process_id):

    if os.path.exists("core." + r4s_process_id):

        LOG.write("remove core file : core.%s\n" %r4s_process_id)
        os.unlink("core." + r4s_process_id)


def determine_msa_format():

    # the routine trys to read the MSA format. If there are errors - reports it to the user.

    # alternative messages to the user
    clustal_msg = "As an alternative you can also re-run the ConSurf session using a different alignment format, preferably ClustAlW format.<br />\n"
    #conversion_progam = "It is recommend to use <a href=\"" + GENERAL_CONSTANTS.MSA_CONVERT + "\">EBI's biosequence conversion tool</a>.<br />\n"
    #msa_info_msg = "<a href =\"" + GENERAL_CONSTANTS.MSA_FORMATS + "\">Read more on MSA formats</a><br />\n"

    LOG.write("determine_msa_format : calling MSA_parser.determine_msa_format(" + vars['user_msa_file_name'] + ")\n")
    msa_format = MSA_parser.determine_msa_format(vars['user_msa_file_name'])

    if msa_format[0] == "err":

        exit_on_error('user_error',"The uploaded " + vars['user_msa_file_name'] + " MSA file is not in one of the formats supported by ConSurf: NBRF/PIR, Pearson (Fasta), Nexus, Clustal, GCG/MSF." )

    LOG.write("determine_msa_format : MSA format is : " + msa_format[1] + "\n")
    msa_format = msa_format[1]

    # format is known. Now check if the sequences are readable and contain legal characters
    LOG.write("determine_msa_format : calling MSA_parser.check_msa_licit_and_size(%s, %s)\n" %(vars['user_msa_file_name'], msa_format))
    # check each sequence in the MSA for illegal characters
    ans = MSA_parser.check_msa_licit_and_size(vars['user_msa_file_name'], msa_format)

    # there was an error
    if ans[0] != "OK":

        LOG.write("determine_msa_format : an error was found while check_msa_licit_and_size\n")
        msg = "The uploaded <a href=\"%s\">MSA file</a>, which appears to be in %s format, " %(vars['user_msa_file_name'], msa_format)

        # if there were illegal characters - report them
        match = re.search(r'^SEQ_NAME: (.+)?', ans[1])
        if match:

            # report in which sequence were the characters found
            try:

                seq_name = match.group(1)

            except:

                seq_name = ""

            try:

                irr_chars = re.search(r'^IRR_CHAR: (.+)?', ans[2]).group(1)

            except:

                irr_chars = ""

            LOG.write("determine_msa_format : found irregular chars: %s in sequence: %s\n" %(irr_chars, seq_name))
            msg += "contains non-standard characters"

            if re.search(r'\S+', irr_chars):

                msg += ": '%s'" %irr_chars

            if re.search(r'\S+', seq_name):

                msg += " in the sequence named '%s'" %seq_name

            msg += ". To fix the format please replace all non-standard characters with standard characters (gaps : \"-\" Amino Acids : \"A\" , \"C\" , \"D\" .. \"Y\") and resubmit your query.<br />\n"

        elif ans[1] == "exception": # the MSA file was not opened correctly

            LOG.write(ans[1] + "\n")
            msg = "could not be read by the server. Please convert it, preferably to ClustAlW format, then resubmit your query.\n" 

        elif ans[1] == "exception":

            LOG.write(ans[1] + "\n")
            msg += "could not be read by the server. Please note that the sequences should contain only standard characters (gaps : \"-\" Amino Acids : \"A\" , \"C\" , \"D\" .. \"Y\"). Please replace all non-standard characters with standard characters and resubmit your query.<br />\n"

            if msa_format != "clustal":

                msg += clustal_msg 

        exit_on_error('user_error', msg)

    else:

        LOG.write("determine_msa_format : MSA is legal!\n")

    vars['unique_seqs'] = ans[1]
    LOG.write("There are " + str(vars['unique_seqs']) + " sequences in the MSA\n")

    return msa_format





def get_info_from_msa(msa_ref):

    # extract all the sequence identifiers and the sequences themselves to the hash MSA_sequences (here - transmitted by reference)
    # create a fasta format file for this msa, since it is more convenient to work with
    # NOTE!! on some formats (I saw it in "pir") the bioPerl modoule which reads the MSA shortens the sequences!
    # also create file containing the reference seq
    # return the reference seq as a string

    msa_format = determine_msa_format() # determine the msa format
    LOG.write("get_info_from_msa : MSA_parser.get_info_from_msa(%s, %s, %s)\n" %(vars['user_msa_file_name'], msa_format, msa_ref))
    ans = MSA_parser.get_info_from_msa(vars['user_msa_file_name'], msa_format, msa_ref)
    LOG.write("get_info_from_msa : answer is: " + ans[0] + "\n")

    if ans[0] != "OK": # an error was found

        LOG.write("get_info_from_msa : Error was found: " + ans[1] + "\n")
        # general error message
        msg = "The uploaded <a href=\"%s\">MSA file</a>, which appears to be in %s format, " %(vars['user_msa_file_name'], msa_format)
        if ans[1] == "exception":

            msg += "could not be read by the server. Please note that the sequences should contain only standard characters (gaps : \"-\" Amino Acids : \"A\" , \"C\" , \"D\" .. \"Y\"). Please replace all non-standard characters with standard characters and resubmit your query.<br />\n"

        elif ans[1] == "could not read msa":

            msg += "could not be read by the server. Please convert it, preferably to ClustAlW format, then resubmit your query.<br />\n"

        elif ans[1] == "no seq id" or re.match(r'^duplicity', ans[1]): # the MSA was read correctly , but there was a problem with the sequences ids

            msg = "An error was found in the uploaded <a href=\"" + vars['user_msa_file_name'] + "\">MSA file</a>:<br />\n"
            if ans[1] == "no seq id":

                msg += "At lease one of the sequences did not have an identifier (a sequence id).<br />\n"

            else:

                match = re.search(r'duplicity (.*)$', ans[1])
                if match:

                    seq_id = match.group(1)
                    msg += "The same sequence identifier: '" + seq_id + "' appeared more than once. "

        exit_on_error('user_error', msg)

    # we save to copies of the msa, one in fasta format and another in clustal format.
    if msa_format == "fasta":

        #os.rename(vars['user_msa_file_name'], vars['msa_fasta'])
        vars['msa_fasta'] = vars['user_msa_file_name']
        convert_msa_format(vars['msa_fasta'], "fasta", vars['msa_clustal'], "clustal")

    elif msa_format == "clustal":

        #os.rename(vars['user_msa_file_name'], vars['msa_clustal'])
        vars['msa_clustal'] = vars['user_msa_file_name']
        convert_msa_format(vars['msa_clustal'], "clustal", vars['msa_fasta'], "fasta")

    else:

        exit_on_error('user_errot', "The msa should be either in format clustal or fasta.")


    num_of_seq = len(msa_ref)
    LOG.write("MSA contains " + str(num_of_seq) + " sequences\n")
    if num_of_seq < 5:

        exit_on_error('user_error',"The MSA file contains only " + str(num_of_seq) + " sequences. The minimal number of homologues required for the calculation is 5.")

    if not form['msa_SEQNAME'] in msa_ref:

        msg = "The query sequence name '%s' is not found in the <a href = \"%s\">MSA file</a>. (It should be written exactly as it appears in the MSA file)" %(form['msa_SEQNAME'], vars['user_msa_file_name'])
        if re.search(r'^\s+', form['msa_SEQNAME']) or re.search(r'\s+$', form['msa_SEQNAME']):

            msg += ".<br />Looks like there are extra spaces. Please check."

        if re.search(r'[\[\]\/\\;]', form['msa_SEQNAME']):

            msg += ".<br />Please note that signs such as: '[', ']', '/', may be problematic as sequence identifier."

        exit_on_error('user_error', msg)

    query_seq = msa_ref[form['msa_SEQNAME']]
    query_seq = query_seq.replace("-", "")

    return query_seq





def create_pipe_file(pipeFile, pipeFile_CBS, seq3d_grades, seq3d_grades_isd, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, length_of_seqres, length_of_atom, pdb_file_name, user_chain, IN_pdb_id_capital, identical_chains):

    # CREATE PART of PIPE
    partOfPipe = "partOfPipe"
    partOfPipe_CBS = "partOfPipe_CBS"

    LOG.write("Calling cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n" %(partOfPipe, vars['unique_seqs'], form['proteins_DB'], seq3d_grades_isd, seq3d_grades, length_of_seqres, length_of_atom, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, form['E_VALUE'], form['ITERATIONS'], form['MAX_NUM_HOMOL'], form['MSAprogram'], form['ALGORITHM'], form['SUB_MATRIX'], "legacy"))
    ans = cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new(partOfPipe, vars['unique_seqs'], form['proteins_DB'], seq3d_grades_isd, seq3d_grades, length_of_seqres, length_of_atom, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, form['E_VALUE'], form['ITERATIONS'], form['MAX_NUM_HOMOL'], form['MSAprogram'], form['ALGORITHM'], form['SUB_MATRIX'], "legacy")
    if ans[0] != "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new FAILED: " + str(ans))

    elif not os.path.exists(partOfPipe) or os.path.getsize(partOfPipe) == 0:

        exit_on_error('sys_error', "create_pipe_file: The file '" + partOfPipe + "' was not found or empty")

    # create the color blind friendly version
    ans = cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new(partOfPipe_CBS, vars['unique_seqs'], form['proteins_DB'], seq3d_grades_isd, seq3d_grades, length_of_seqres, length_of_atom, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, form['E_VALUE'], form['ITERATIONS'], form['MAX_NUM_HOMOL'], form['MSAprogram'], form['ALGORITHM'], form['SUB_MATRIX'], "cb")
    if ans[0] != "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new FAILED: " + str(ans))

    elif not os.path.exists(partOfPipe_CBS) or os.path.getsize(partOfPipe_CBS) == 0:

        exit_on_error('sys_error', "create_pipe_file: The file '" + partOfPipe_CBS + "' was not found or empty")

    LOG.write("going to extract data from the pdb, calling: cp_rasmol_gradesPE_and_pipe.extract_data_from_pdb(%s)\n" %pdb_file_name)
    header_pipe = cp_rasmol_gradesPE_and_pipe.extract_data_from_pdb(pdb_file_name)
    if header_pipe[0] != "OK":

        LOG.write(str(header_pipe))

    # GET THE FILE NAMES
    msa_filename = ""
    msa_query_seq_name = ""
    if form['msa_input'] is not None:

        msa_filename = form['msa_input']
        msa_query_seq_name = form['msa_SEQNAME']

    tree_filename = ""
    if form['tree_name'] is not None:

        tree_filename = vars['tree_file']

    # GET THE CURRENT TIME
    completion_time = str(datetime.now().time())
    run_date = str(datetime.now().date())

    # USE THE CREATED PART of PIPE to CREATE ALL THE PIPE TILL THE PDB ATOMS (DELETE THE PART PIPE)
    LOG.write("Calling cp_rasmol_gradesPE_and_pipe.create_consurf_pipe_new(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n" %(vars['working_dir'], IN_pdb_id_capital, user_chain, header_pipe, pipeFile, identical_chains, partOfPipe, vars['working_dir'], form['Run_Number'], msa_filename, msa_query_seq_name, tree_filename, vars['submission_time'], completion_time, run_date))
    ans = cp_rasmol_gradesPE_and_pipe.create_consurf_pipe_new(vars['working_dir'], IN_pdb_id_capital, user_chain, header_pipe, pipeFile, identical_chains, partOfPipe, vars['working_dir'], form['Run_Number'], msa_filename, msa_query_seq_name, tree_filename, vars['submission_time'], completion_time, run_date)
    if ans != "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_consurf_pipe_new FAILED: " + ans)

    # USE THE CREATED PART of PIPE to CREATE ALL THE PIPE TILL THE PDB ATOMS (DELETE THE PART PIPE) - Color friendly version
    LOG.write("Calling cp_rasmol_gradesPE_and_pipe.create_consurf_pipe_new(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n" %(vars['working_dir'], IN_pdb_id_capital, user_chain, header_pipe, pipeFile_CBS, identical_chains, partOfPipe_CBS, vars['working_dir'], form['Run_Number'], msa_filename, msa_query_seq_name, tree_filename, vars['submission_time'], completion_time, run_date))
    ans = cp_rasmol_gradesPE_and_pipe.create_consurf_pipe_new(vars['working_dir'], IN_pdb_id_capital, user_chain, header_pipe, pipeFile_CBS, identical_chains, partOfPipe_CBS, vars['working_dir'], form['Run_Number'], msa_filename, msa_query_seq_name, tree_filename, vars['submission_time'], completion_time, run_date)
    if ans != "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_consurf_pipe_new FAILED: " + ans)

    # Add the PDB data to the pipe
    LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.add_pdb_data_to_pipe(%s, %s)\n" %(pdb_file_name, pipeFile))
    ans = cp_rasmol_gradesPE_and_pipe.add_pdb_data_to_pipe(pdb_file_name, pipeFile)
    if ans != "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.add_pdb_data_to_pipe FAILED: " + ans)

    if not os.path.exists(pipeFile) and os.path.getsize(pipeFile) == 0:

        exit_on_error('sys_error', "create_pipe_file : Did not create the FGiJ output")

    # Add the PDB data to the pipe - color blind version
    LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.add_pdb_data_to_pipe(%s, %s)\n" %(pdb_file_name, pipeFile_CBS))
    ans = cp_rasmol_gradesPE_and_pipe.add_pdb_data_to_pipe(pdb_file_name, pipeFile_CBS)
    if ans != "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.add_pdb_data_to_pipe FAILED: " + ans)

    if not os.path.exists(pipeFile_CBS) and os.path.getsize(pipeFile_CBS) == 0:

        exit_on_error('sys_error', "create_pipe_file : Did not create the FGiJ output")




def create_gradesPE_ConSeq(ref_Solv_Acc_Pred, type):

    # Create ConSeq gradesPE file

    ans = cp_rasmol_gradesPE_and_pipe.create_gradesPE_ConSeq(gradesPE_Output, residue_freq, ref_Solv_Acc_Pred, vars['gradesPE'], type, form['ALGORITHM'], vars['layers_array'])
    if not ans == "OK":

        exit_on_error('sys_error', "create_gradesPE_ConSeq : rasmol_gradesPE_and_pipe." + ans)

    elif not os.path.exists(vars['gradesPE']) or os.path.getsize(vars['gradesPE']) == 0:

        exit_on_error('sys_error', "create_gradesPE_ConSeq : the file " + vars['gradesPE'] + " was not found or empty")






def compare_atom_to_query(Query_seq, ATOM_seq, pairwise_aln, PDB_Name):

    # in case there are both seqres and atom fields, checks the similarity between the 2 sequences.

    clustalw = GENERAL_CONSTANTS.CLUSTALW
    two_fastas = PDB_Name + ".fasta2"
    clustalw_out = PDB_Name + ".out"

    [first_seq, second_seq] = compare_two_fastas(two_fastas, Query_seq, ATOM_seq, clustalw_out, pairwise_aln, "QUERY")
    return(first_seq, second_seq)



def Get_NACSES_buried_Exposed(structure_file, chain = "", atom_positionFILE = ""):
    
    return("failed", "")

    buried_cutoff = 16
    exposed_cutoff = 16

    ATOM_pos_To_SEQRES = {} # will have values if atom_positionFILE provided
    Buried_Exposed = {} # key: poistion in model, value: b|e according to the cutoff
    NACSSESS_DIR = "" #GENERAL_CONSTANTS.NACSSESS_DIR
    cmd = NACSSESS_DIR + "naccess " + structure_file
    submit_job_to_Q("NACSSESS", cmd)

    basename = structure_file[:-4] # delete ending
    out_Files = [basename + ".asa", basename + ".rsa", basename + ".log"]
    for file in out_Files:

        if not os.path.exists(file):

            LOG.write("ERROR, command: '%s' FAILED to create %s\n" %(str(cmd), file))

    """
    if atom_positionFILE != "":

        try:

            ATOM_TO_SEQRES = open(atom_positionFILE, 'r')

        except:

            return("failed", "")

        line = ATOM_TO_SEQRES.readline()
        while line != "":

            [Res, SEQRES_POS, ATOM_POS] = line.rstrip().split()

            if not ATOM_POS in ATOM_pos_To_SEQRES:

                ATOM_pos_To_SEQRES[ATOM_POS] = {}

            ATOM_pos_To_SEQRES[ATOM_POS]['RES'] = Res
            ATOM_pos_To_SEQRES[ATOM_POS]['SEQRES'] = SEQRES_POS

            line = ATOM_TO_SEQRES.readline()

        ATOM_TO_SEQRES.close()
    """
    # parse rsa file for b/e
    try:

        RSA = open(basename + ".rsa", 'r')

    except:

        return("failed", "")

    """
    try:

        OUT = open(basename + ".NACSSES.rsa", 'w')

    except:

        return("failed", "")

    OUT.write("ATOM_POS\tAll_atoms_REL\tb/e")
    if atom_positionFILE != "":

        OUT.write("\tSEQRES_POS\tRES\n")

    else:

        OUT.write("\n")
    """
    line = RSA.readline()
    while line != "":

        words = line.split()
        if len(words) >= 13: # informative line

            if len(words) == 13: # pdb file with no chain id

                pos = words[2]
                """
                if atom_positionFILE != "" and pos in ATOM_pos_To_SEQRES and 'SEQRES' in ATOM_pos_To_SEQRES[pos]:

                    pos = ATOM_pos_To_SEQRES[pos]['SEQRES']
                """
                if float(words[4]) >= exposed_cutoff: # exposed according to All atoms REL

                    Buried_Exposed[pos] = "e"

                elif float(words[4]) < buried_cutoff: # buried according to All atoms REL

                    Buried_Exposed[pos] = "b"

                """
                OUT.write("%s\t%s\t%s" %(words[2], words[4], Buried_Exposed[pos]))

                if atom_positionFILE != "" and words[2] in ATOM_pos_To_SEQRES and 'SEQRES' in ATOM_pos_To_SEQRES[words[2]]:

                    OUT.write("\t%s\t%s\n" %(ATOM_pos_To_SEQRES[words[2]]['SEQRES'], ATOM_pos_To_SEQRES[words[2]]['RES']))

                else:

                    OUT.write("\n")
                """

            elif chain == words[2]: # the relevant chain

                pos = words[3]
                """
                if atom_positionFILE != "" and pos in ATOM_pos_To_SEQRES and 'SEQRES' in ATOM_pos_To_SEQRES[pos]:

                    pos = ATOM_pos_To_SEQRES[pos]['SEQRES']
                """

                if float(words[5]) >= exposed_cutoff: # exposed according to All atoms REL

                    Buried_Exposed[pos] = "e"

                elif float(words[5]) < buried_cutoff: # buried according to All atoms REL

                    Buried_Exposed[pos] = "b"

                """
                OUT.write("%s:%s\t%s\t%s" %(words[2], words[3], words[5], Buried_Exposed[pos]))

                if atom_positionFILE != "" and words[3] in ATOM_pos_To_SEQRES and 'SEQRES' in ATOM_pos_To_SEQRES[words[3]]:

                    OUT.write("\t%s\t%s\n" %(ATOM_pos_To_SEQRES[words[3]]['SEQRES'], ATOM_pos_To_SEQRES[words[3]]['RES']))

                else:

                    OUT.write("\n")
                """

        line = RSA.readline()

    RSA.close()
    #OUT.close()

    return("OK", Buried_Exposed)




def Project_ConSurf_On_Model(PDB_File, chain, OutDir, Query_Seq, r4s_out, protein_MSA_clustalw, pisa, Type):

    if pisa:
	
        prefix = form['pdb_ID'] + "_" + chain + "_pisa" # this begins many of the variable names

    else:
	
        prefix = vars['Used_PDB_Name']

    ATOM_Pos_File = prefix + ".ATOM_Pos_File.txt" 
    PDB_Object = pdbParser.pdbParser()
    PDB_Object.read(PDB_File, chain, Type, ATOM_Pos_File)

    #D_Nuc = PDB_Object.IS_D_NUC(chain)

    # FIND IDENTICAL CHAINS
    identical_chains_local = find_identical_chains_in_PDB_file(PDB_Object, chain, 2)

    [SEQRES_seq_local, ATOM_seq_loacl, ATOM_seq_without_X_local] = get_seqres_atom_seq(PDB_Object, chain, PDB_File, True)
    if SEQRES_seq_local == "user_error":

        return("")
	
		
    # an array to hold all the information that should be printed to gradesPE
    # in each array's cell there is a hash for each line from r4s.res.
    # POS: position of that aa in the sequence ; SEQ : aa in one letter ;
    # GRADE : the given grade from r4s output ; COLOR : grade according to consurf's scale
    gradesPE_Output_local = []

    residue_freq_local = {} # for each position in the MSA, detail the residues
    position_totalAA_local = {} # for each position in the MSA, details the total number of residues

    # these arrays will hold for each grade, the residues which corresponds to it.
    # there are 2 arrays: in the isd_residue_color, a grade with insufficient data, *, will classify to grade 10
    # in the no_isd_residue_color, the grade will be given regardless of the * mark
    # PLEASE NOTE : the [0] position in those arrays is empty, because each position corresponds a color on a 1-10 scale
    no_isd_residue_color_local = [[],[],[],[],[],[],[],[],[],[],[]]
    isd_residue_color_loacl = [[],[],[],[],[],[],[],[],[],[],[]]


    pairwise_aln = prefix + ".aln" # created by compare_atom_to_query()

    gradesPE = prefix + "_consurf_grades.txt"
    scoresFile = prefix + ".consurf.scores"

    Server_Results_Path = "results/" + form['Run_Number'] + "/" # in the global dir

    pymol_color_script = GENERAL_CONSTANTS.CONSURF_URL + "pyMOL/consurf_new.py"


    # Atoms Section with consurf grades instead TempFactor Field
    ATOMS_with_ConSurf_Scores = prefix + "_ATOMS_section_With_ConSurf.pdb"
    ATOMS_with_ConSurf_Scores_isd =  prefix + "_ATOMS_section_With_ConSurf_isd.pdb"

    pdb_file_with_score_at_TempFactor = prefix + "_With_Conservation_Scores.pdb"

    vars['zip_list'].append(gradesPE)
    vars['zip_list'].append(ATOMS_with_ConSurf_Scores)
    #vars['zip_list'].append(ATOMS_with_ConSurf_Scores_isd)
    vars['zip_list'].append(pdb_file_with_score_at_TempFactor)

    assign_colors_according_to_r4s_layers(gradesPE_Output_local, r4s_out)
    [Query_Seq_with_gaps, ATOM_seq_loacl_with_gaps] = compare_atom_to_query(Query_Seq, ATOM_seq_loacl, pairwise_aln, prefix) # Compare the PDB Atom With Reference Seq
    read_residue_variety(residue_freq_local, position_totalAA_local)
    #res = create_atom_position_file(chain, PDB_File, ATOM_Pos_File, Type, pisa) # this file will be used later to create the output which aligns rate4site sequence with the ATOM records

    r4s2pdb_local = {} # key: poistion in SEQRES/MSA, value: residue name with position in atom (i.e: ALA22:A)
    [length_of_seqres_loacl, length_of_atom_loacl] = match_pdb_to_seq(r4s2pdb_local, chain, Query_Seq_with_gaps, ATOM_seq_loacl_with_gaps, ATOM_Pos_File, Type)

    # The following routine, apart from creating the file "consurf.grades" also collects the information in order to create
    # the RasMol scripts and a variable which will be used in the "pipe" file, that holds a string with the grades.
    # In the pipe file this var is called: seq3d_grades_isd and seq3d_grades
    #[seq3d_grades_isd_local, seq3d_grades_isd_local] = create_gradesPE_ConSurf(gradesPE, r4s2pdb_local, gradesPE_Output_local, residue_freq_local, no_isd_residue_color_local, isd_residue_color_loacl, Type, D_Nuc)
    [seq3d_grades_isd_local, seq3d_grades_local] = create_gradesPE_ConSurf(gradesPE, r4s2pdb_local, gradesPE_Output_local, residue_freq_local, no_isd_residue_color_local ,isd_residue_color_loacl, Type, "")
    cp_rasmol_gradesPE_and_pipe.create_score_file(gradesPE_Output_local, scoresFile, r4s2pdb_local)

    # This will create the 2 rasmol scripts (one with isd and one without)
    #create_rasmol(chain, rasmolFILE, rasmol_isdFILE, rasmolFILE_CBS, rasmol_isdFILE_CBS, no_isd_residue_color_local, isd_residue_color_loacl)

    # This will create the pipe file for FGiJ
    pipeFile_local = prefix + "_consurf_firstglance.pdb"
    pipeFile_CBS_local = prefix + "_consurf_firstglance_CBS.pdb"
    create_pipe_file(pipeFile_local, pipeFile_CBS_local, seq3d_grades_local, seq3d_grades_isd_local, isd_residue_color_loacl, no_isd_residue_color_local, length_of_seqres_loacl, length_of_atom_loacl, PDB_File, chain, PDB_File.upper(), identical_chains_local)
    """
    if pisa:

        saveJobInfoJson(jobInfoFile, scoresFile, PDB_File, "", form['JOB_TITLE'], "PISA", "", chain, identical_chains)

    else:

        saveJobInfoJson(jobInfoFile, scoresFile, PDB_File, "", form['JOB_TITLE'], "", "", "A", "A")
    """

    if pisa:

        jobInfoFile = "PISA1.job_info.json"
        Hash_Json = {'scoresFile' : scoresFile, 'pdbFile' : PDB_File, 'jobTitle' : form['JOB_TITLE'], 'subTitle' : "PISA", 'chainId' : chain}
        saveJobInfoJson(jobInfoFile, Hash_Json)

    else:

        # json for NGL
        jobInfoFile = "job_info.json"
        Hash_Json = {'scoresFile' : scoresFile, 'pdbFile' : PDB_File, 'jobTitle' : form['JOB_TITLE'], 'chainId' : "A"}
        saveJobInfoJson(jobInfoFile, Hash_Json)

        # json for embedded NGL
        relative_path = "../results/" + form['Run_Number'] + "/"
        Hash_Json = {'scoresFile' : relative_path + scoresFile, 'pdbFile' : relative_path + PDB_File, 'chainId' : "A"}
        saveJobInfoJson("job_info_small.json", Hash_Json)



    # Create ATOMS section and replace the TempFactor Column with the ConSurf Grades (will create also isd file if relevant)
    chimera_pymol_input = replace_TmpFactor_Consurf_Scores(chain, PDB_File, gradesPE, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd)
    create_chimera(chimera_pymol_input, prefix + "_")
    create_pymol(chimera_pymol_input, prefix + "_")

    # Create PDB file and replace the TempFactor Column with the Rate4Site grades
    replace_TmpFactor_Rate4Site_Scores(chain, PDB_File, gradesPE, pdb_file_with_score_at_TempFactor)

    if pisa:

        pisa_html()

def Run_HHPred():

    # Get a file with sequence and run the HHPred algorithm (to find template and model structure using modeller)

    HHPRED = GENERAL_CONSTANTS.HH_PRED
    #OutDir = "HH_Pred/"
    SeqFile = vars['working_dir'] + vars['protein_seq']
    #OutDir = vars['working_dir'] + "HH_Pred/"
    #OutModel = vars['working_dir'] + vars['HHPred_Model']

    #os.mkdir(OutDir)
    
    #cmd = "cd %s\n" %OutDir
    cmd = "/powerapps/share/centos7/hh-suite3/hh-suite/build/bin/hhblits -i %s -o query.hhr -oa3m query.a3m -n 3 -cov 20 -d /bioseq/hhpred-databases/db/uniclus30/UniRef30_2020_06 -cpu 16\n" %SeqFile
    cmd += "/powerapps/share/centos7/hh-suite3/hh-suite/build/bin/hhsearch -i query.a3m -d /bioseq/hhpred-databases/db/pdb70 -o results.hhr -cpu 16\n"
    cmd += "module load miniconda/miniconda3-4.7.12-environmentally\n"
    cmd += "conda activate /powerapps/share/centos7/miniconda/miniconda3-4.7.12-environmentally/envs/PDBX\n"
    cmd += "pip list installed | grep pdbx\n"
    cmd += "python -c \"import pdbx\"\n"
    cmd += "python /bioseq/hhpred-databases/db/hhpred/hhmakemodel.py results.hhr /bioseq/hhpred-databases/db/pdb70/ PIR.pir ./ -m 1\n"
    cmd += "conda deactivate\n"
    cmd += "cp /bioseq/hhpred-databases/db/hhpred/create_model.py create_model.py\n"
    cmd += "module load modeller-9.22\n"
    cmd += "mod9.22 create_model.py\n"

    submit_job_to_Q("HHPRED", cmd, "NO")
    if os.path.exists("UKNP.B99990001.pdb"):
        
        os.rename("UKNP.B99990001.pdb", vars['HHPred_Model'])
        return("OK")
  
    else:
	 

        try:
		
            HHPRED_FAILED = open("HHPRED_FAILED", 'w')
			
        except:
		
            exit_on_error('sys_error', "Run_HHPred : Cannot open the file HHPRED_FAILED for writing!")
			
        HHPRED_FAILED.close()
        return("NOT_OK")

"""
def Run_HHPred():

    # Get a file with sequence and run the HHPred algorithm (to find template and model structure using modeller)

    HHPRED = GENERAL_CONSTANTS.HH_PRED
    OutDir = "HH_Pred/"
    SeqFile = vars['working_dir'] + vars['protein_seq']
    OutDir = vars['working_dir'] + "HH_Pred/"
    OutModel = vars['working_dir'] + vars['HHPred_Model']

    os.mkdir(OutDir)

    cmd = ["ssh", "bioseq@power", "perl", HHPRED, "-i", SeqFile, "-d", OutDir, "-o", OutModel, ">", OutDir + "HHPred.std"]
    LOG.write("Running HHPred: %s\n" %str(cmd))
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    p.communicate()

    while not os.path.exists(OutModel):

        time.sleep(1)
"""

def predict_solvent_accesibility(aln, hssp, pred):

    # runs the PACC algorithm to calculate burried/exposed

    query_name = "'" + vars['msa_SEQNAME'] + "'"
    MSA_to_HSSP = GENERAL_CONSTANTS.MSA_to_HSSP
    PACC_path = GENERAL_CONSTANTS.PACC_path

    HSSP = open(hssp, "w")
    HSSP.close()
    PRED = open(pred, "w")
    PRED.close()

    # run the script that turns the MSA to hssp file
    cmd = "perl %s %s %s > %s" %(MSA_to_HSSP, query_name, aln, hssp)
    LOG.write("predict_solvent_accesibility:\n%s\n" %cmd)
    submit_job_to_Q("HSSP", cmd)

    # check if the hssp file has been created and is not empty
    while not os.path.exists(hssp) or os.path.getsize(hssp) == 0:

        time.sleep(1)
        #exit_on_error('sys_error', "predict_solvent_accesibility: The file " + hssp + " does not exist or contains no data\n")

    LOG.write("predict_solvent_accesibility: The file " + hssp + " was created and now contains data\n")

    # run the script that produce the prediction file
    cmd = "cd %s\n./run.sh %s %s\ncd %s" %(PACC_path, vars['working_dir'] + hssp, vars['working_dir'] + pred, vars['working_dir'])
    LOG.write("predict_solvent_accesibility: %s\n" %cmd)
    submit_job_to_Q("SOLVENT", cmd)

    # check if the pred file has been created and is not empty
    while not os.path.exists(pred) or os.path.getsize(pred) == 0:

        time.sleep(1)
        #exit_on_error('sys_error', "predict_solvent_accesibility: The file " + pred + " does not exist or contains no data\n")



def parse_blast_alignment(BLAST_PDB_out_file):

    PDB = ""
    BLAST_OUTPUT = open(BLAST_PDB_out_file, 'r+')
    line = BLAST_OUTPUT.readline()
    while line != "":

        match1 = re.match(r'^\s*', line)
        if match1:

            break

        else:

            print(line)

        line = BLAST_OUTPUT.readline()

    Alignment = line

    BLAST_OUTPUT.seek(0, 0)
    line = BLAST_OUTPUT.readline()
    while line != "":

        match2 = re.match(r'>([A-Za-z0-9]{5})', line)
        if match2:

            if PDB != "":

                vars['Templates_Alignments']['PDB'] = Alignment

            PDB = match2.group(1)
            line = re.sub(r' ', "&nbsp;", line)
            line = re.sub(r'Sbjct:', PDB + ":", line)
            Alignment = line
            Alignment += "<BR>"
            line = BLAST_OUTPUT.readline()

        else:

            line = re.sub(r' ', "&nbsp;", line)
            line = re.sub(r'Sbjct:', PDB + ":", line)
            match3 = re.match(r'Database', line)
            if match3:

                break

            Alignment += line
            Alignment += "<BR>"

        line = BLAST_OUTPUT.readline()

    BLAST_OUTPUT.close()
    vars['Templates_Alignments']['PDB'] = Alignment



def blast_vs_PDB(query, Blast_PDB_Out_File):

    LOG.write("blast_vs_PDB(%s, %s)\n" %(query, Blast_PDB_Out_File))
    cmd = [pgp, "-i", query, "-e", form['E_VALUE'], "-d", vars['pdb_db'], "-j", str(form['ITERATIONS']), "-v", str(vars['max_PDB_homologues_to_display']), "-b", str(vars['max_PDB_homologues_to_display']), "-o", Blast_PDB_Out_File, "-F", "F", "-T", "F"]
    LOG.write("run_blast : running: %s\n" %cmd)
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    p.communicate()



def upload_protein_sequence():

    if form['FASTA_text'] is not None:

        if form['FASTA_text'].count('>') > 1:

            exit_on_error('user_error', "The protein input <a href = \"%s\">sequence</a> contains more than one FASTA sequence. If you wish to upload MSA, please upload it as a file." %form['FASTA_text'])

        # delete sequence name and white spaces
        protein_seq_string = re.sub(r'>.*\n', "", form['FASTA_text'])
        protein_seq_string = re.sub(r'\s', "", protein_seq_string)

        try:

            UPLOADED = open(vars['protein_seq'], 'w')

        except:

            exit_on_error('sys_error', "upload_protein_sequence : Cannot open the file " + vars['protein_seq'] + "for writing!")

        LOG.write("upload_protein_sequence : add '>' to the protein sequence\n")
        UPLOADED.write(">SEQ\n" + protein_seq_string)
        UPLOADED.close()

    elif os.path.exists(vars['protein_seq']) and os.path.getsize(vars['protein_seq']) != 0:

        try:

            UPLOADED = open(vars['protein_seq'], 'r')

        except:

            exit_on_error('sys_error', "upload_protein_sequence : Cannot open the file " + vars['protein_seq'] + "for writing!")

        protein_seq_string = UPLOADED.read()
        UPLOADED.close()

        if protein_seq_string.count('>') > 1:

            exit_on_error('user_error', "The protein input <a href = \"%s\">sequence</a> contains more than one FASTA sequence. If you wish to upload MSA, please upload it as a file." %protein_seq_string)

        # delete sequence name and white spaces
        protein_seq_string = re.sub(r'>.*\n', "", protein_seq_string)
        protein_seq_string = re.sub(r'\s', "", protein_seq_string)

        try:

            UPLOADED = open(vars['protein_seq'], 'w')

        except:

            exit_on_error('sys_error', "upload_protein_sequence : Cannot open the file " + vars['protein_seq'] + "for writing!")

        LOG.write("upload_protein_sequence : add '>' to the protein sequence\n")
        UPLOADED.write(">SEQ\n" + protein_seq_string)
        UPLOADED.close()

    else:

        exit_on_error('sys_error', 'upload_protein_sequence : no user sequence.')



    if form['DNA_AA'] == "AA":

        if re.match(r'^[actguACTGUNn]+$', protein_seq_string):

            exit_on_error('user_error',"It seems that the protein input <a href = \"%s\">sequence</a> is only composed of Nucleotides (i.e. :A,T,C,G). Please note that you chose to run the server based on amino acids sequnce and not DNA / RNA sequence.<br />You may translate your sequence to amino acids and resubmit your query, or alternatively choose to analyze nucleotides.<br />" %protein_seq_string)

    else:

        if not re.match(r'^[actguACTGUNn]+$', protein_seq_string):

            exit_on_error('user_error',"It seems that the protein input <a href = \"%s\">sequence</a> is only composed of Amino Acids. Please note that you chose to run the server based on nucleotides sequnce and not protein sequence.<br />You may resubmit your query and choose to analyze Amino Acids.<br />" %protein_seq_string)

    vars['protein_seq_string'] = protein_seq_string
    vars['query_string'] = "Input_protein_seq"



def send_finish_email_to_user():
	
    email_subject = "'The results of your ConSurf run titled " + form['JOB_TITLE'] + " are ready'"
    email_message = "'Hello,\\n\\nThe results for your ConSurf run are ready at:\\n" + vars['run_url'] + "\\n\\nRunning Parameters:\\n"
	
    if form['pdb_ID'] is not None: #ConSurf Mode

        email_message += "PDB: %s\\nCHAIN: %s\\n" %(form['pdb_ID'], form['PDB_chain'])

    if 'user_msa_fasta' in vars:

        email_message += "Alignment: %s\\n" %vars['user_msa_fasta']

    else:

        email_message += "Alignment: Created by ConSurf using %s and %s\\n" %(form['proteins_DB'], form['MSAprogram'])

    email_message += "\\nPlease note: the results will be kept on the server for three months.'"
    send_email(email_subject, email_message)
	

def send_email(email_subject, email_message):

    msg = ["perl", "-w", vars['send_email_cmd'], "-f", GENERAL_CONSTANTS.ADMIN_EMAIL, "-t", form['user_email'], "-u", email_subject, "-xu", vars['userName'], "-xp", vars['userPass'], "-s", vars['smtp_server'], "-m", email_message]
    LOG.write("MESSAGE:%s\nCOMMAND:%s\n" %(email_message, str(msg)))
    #email_dir = vars['working_dir'] + "mail"
    p = subprocess.Popen(msg, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    out, err = p.communicate()
    if not "successfully" in out:

        LOG.write("send_mail: The message was not sent successfully. system returned: %s\n" %out)

def send_administrator_mail_on_error(email_message):

    email_subject = "'System ERROR has occurred on ConSurf: %s'" %vars['run_url']
    msg = ["perl", "-w", vars['send_email_cmd'], "-f", "'bioSequence@tauex.tau.ac.il'", "-t", "'bioSequence@tauex.tau.ac.il'", "-u", email_subject, "-xu", vars['userName'], "-xp", vars['userPass'], "-s", vars['smtp_server'], "-m", email_message]
    LOG.write("MESSAGE:%s\nCOMMAND:%s\n" %(email_message, str(msg)))
    #email_dir = vars['working_dir'] + "mail"
    p = subprocess.Popen(msg, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    out, err = p.communicate()
    if not "successfully" in out:

        LOG.write("send_mail: The message was not sent successfully. system returned: %s\n" %out)








def extract_and_print_diversity_matrix():

    LOG.write("extract_and_print_diversity_matrix Calling: rate4site_routines.extract_diversity_matrix_info(%s)\n" %(vars['r4s_log']))
    [ans, diversity_matrix_av_dist, diversity_matrix_low_bound, diversity_matrix_up_bound] = rate4site_routines.extract_diversity_matrix_info(vars['r4s_log'])
    if ans == "OK":

        if diversity_matrix_av_dist != "" and diversity_matrix_low_bound != "" and diversity_matrix_up_bound != "":

            OUTPUT.write("<b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Alignment details</b><br>\n")
            OUTPUT.write("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The average number of replacements between any two sequences in the alignment;<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A distance of 0.01 means that on average, the expected replacement for every 100 positions is 1.<br>\n")
            OUTPUT.write("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<i>Average pairwise distance</i> : %s<br>\n" %diversity_matrix_av_dist)
            OUTPUT.write("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<i>Lower bound</i> : %s<br>\n" %diversity_matrix_low_bound)
            OUTPUT.write("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<i>Upper bound</i> : %s<br>\n" %diversity_matrix_up_bound)

        else:

            LOG.write("diversity_matrix : av_dist is: '%s', low_bound is: '%s', up_bound is: '%s'\n" %(diversity_matrix_av_dist, diversity_matrix_low_bound, diversity_matrix_up_bound))

    else:

        LOG.write(ans)


def zip_all_outputs():

    zipObj = ZipFile(vars['All_Outputs_Zip'], 'w')
	
    for file in vars['zip_list']:
	    
        zipObj.write(file)
	
    zipObj.close()




def saveJobInfoJson(jobInfoFile, hash):


    try:

        JSON = open(jobInfoFile ,'w')

    except:

        exit_on_error('sys_error', "Could not open the file " + jobInfoFile + " for writing.")

    json.dump(hash, JSON)

    JSON.close()


def buildHtml_OutMsg(txt):

    return("\n<ul><li>" + txt + "</li></ul>\n")

def buildHtml_3dViewerLinks(runNumber, nglFile, pipeUrlNoSuffix, description = ""):

    # Generates HTML containing links to 3D model viewers (NGL, FGiJ)

    nglLink = GENERAL_CONSTANTS.CONSURF_URL + "ngl/viewer.php?job=" + str(runNumber)
    if nglFile != "":

        nglLink += "&file=" + nglFile

    nglLinkCbs = nglLink + "&cbs=1"
    fgijLink = "/fgij/fg.htm?mol=" + pipeUrlNoSuffix + "_pipe.pdb"
    fgijLinkCbs = "/fgij/fg.htm?mol=" + pipeUrlNoSuffix + "_pipe_CBS.pdb" # pipe for color blind friendly

    #  projected on PDB_Name_For_Out
    res = buildHtml_OutMsg("<A HREF='" + nglLink + "' TARGET='_blank'><b>View ConSurf results" + description + "</b></A> with NGL viewer [load results using color-blind friendly scale <a href='" + nglLinkCbs + "' TARGET='_blank'>here</a>]<br>")
    res += buildHtml_OutMsg("<A HREF='" + fgijLink + "' TARGET='_blank'>View ConSurf results" + description + "</A> with FirstGlance in Jmol [load results using color-blind friendly scale <a href='" + fgijLinkCbs + "' TARGET='_blank'>here</a>]<br><font size='-1'><span style='color:red;'>Note:</span> FirstGlance version was updated. In case you can't see the molecule and get warnings, press OK and afterwards press Ctrl+F5 keys to reload the new version. You can read more about the issue <a href='https://en.wikipedia.org/wiki/Wikipedia:Bypass_your_cache' TARGET='_blank'>here.</a></font><br>")

    return res




def Prepare_tree_View():

    LOG.write("Prepare_tree_View()\n")
    try:

        TREE_VIEWER = open(vars['tree_viewer_file'], 'w')

    except:

        LOG.write("Can't Open the TreeViewer File For Writing '%s'\n" %vars['tree_viewer_file'])

    ConSurf_CGI_Dir = GENERAL_CONSTANTS.CONSURF_CGI_DIR
    ConSurf_URL = GENERAL_CONSTANTS.CONSURF_URL
    WASABI_page = "%s/wasabi/?url=%s" %(GENERAL_CONSTANTS.CONSURF_URL, vars['WASABI_XML'])

    TREE_VIEWER.write("""<HTML>
    <p>
    <!-- Applet blocked by Java Security? enable it: Control Panel &#8594; Programs &#8594; Java &#8594; Security &#8594; Edit Site List &#8594; add https://ploidb.tau.ac.il/
    <br>
    "This plugin is not supported" (Google Chrome)? follow these <a href="https://support.google.com/chrome/answer/6213033" target="_blank"><i>Instructions</i></a> -->
    If you do not see the phylogenetic tree below, your web browser block Java.<br>
    You can: <ol><li>See the phylogenetic tree together with the MSA using WASABI (no need of Java) <a href="%s">here</a>
    <li>Follow the instructions <a href="https://proteopedia.org/wiki/index.php/Installing_and_enabling_Java" target="_blank">here</a> to enable Java<br>
    <br>
    </p>
    <APPLET CODEBASE="%s"
    CODE="treeViewer.Main.class"
        ARCHIVE="treeViewer/freehep-base.jar,
            treeViewer/freehep-graphics2d.jar,
            treeViewer/freehep-graphicsio-pdf.jar,
            treeViewer/freehep-graphicsio-ps.jar,
            treeViewer/freehep-graphicsio.jar"
        WIDTH=100%% HEIGHT=100%%>
    <PARAM NAME=TREEFILE VALUE="results/%s">
    <PARAM NAME=USERDIR VALUE="%s">
    <PARAM NAME=CGIPATH VALUE="%s">
    <PARAM NAME=HELPURL VALUE="%s/treeHelp/index.html">
    <PARAM NAME=SCRIPTS VALUE="%s">
    </APPLET>
    <!-- both ARCHIVE and TREEFILE are using the relative path given in CODEBASE -->
    <!-- the HTLPURL is opened from the applet in the "help" menu on a new window -->""" %(WASABI_page, ConSurf_URL, vars['tree_file'], vars['working_dir'], ConSurf_CGI_Dir, ConSurf_URL, ConSurf_CGI_Dir))

    TREE_VIEWER.close()


def find_identical_chains_in_PDB_file(pdb_Object, query_chain, number):

    # Looking in the PDB for chains identical to the original chain
    ATOMS = pdb_Object.get_ATOM_withoutX()

    # string with identical chains
    identicalChains = query_chain

    # looking for chains identical to the original chain
    for chain in ATOMS:

        if query_chain != chain:

            chain_length = len(ATOMS[chain])
            OrgChain_length = len(ATOMS[query_chain])

            # if length not similar, skip
            if min(OrgChain_length, chain_length)/max(OrgChain_length, chain_length) <= 0.9:

                continue

            # compare the two chains with clustalw

            # file for two fastas
            fastaFileName = query_chain + "_" + chain + "_twoFastas"
            # file with clustalw output
            clustalwOutputFile = query_chain + "_" + chain + "_twoFastas.aln"

            try:

                FASTAS = open(fastaFileName, 'w')

            except:

                exit_on_error('sys_error', "find_identical_chains_in_PDB_file: can't open the file " + fastaFileName + " for writing.")

            FASTAS.write(">%s\n%s\n>%s\n%s\n" %(query_chain, ATOMS[query_chain], chain, ATOMS[chain]))
            FASTAS.close()

            """
            cmd = ["/bioseq/Programs/ClustalW_2.0.10/clustalw-2.0.10-linux-i386-libcppstatic/clustalw2", "-INFILE=" + fastaFileName, "-gapopen=1", "-OUTPUT=" + clustalwOutputFile]
            # run clustalw
            p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
            out, err = p.communicate()
            """

            # run clustalw
            cmd = "/bioseq/Programs/ClustalW_2.0.10/clustalw-2.0.10-linux-i386-libcppstatic/clustalw2 -INFILE=" + fastaFileName + " -gapopen=1 -OUTPUT=" + clustalwOutputFile
            submit_job_to_Q("compare_chain_%s_to_chain_%s_%d" %(chain, query_chain, number), cmd)


            # read cluatal output
            try:

                ALN = open(clustalwOutputFile, 'r')

            except:

                exit_on_error('sys_error', "find_identical_chains_in_PDB_file: can't open the file " + clustalwOutputFile + " for reading.")

            # count number of asterisks
            countAsterisk = (ALN.read()).count("*")
            ALN.close()

            # check if to add chain
            if countAsterisk/max(OrgChain_length, chain_length) > 0.95:

                identicalChains += " " + chain

    return identicalChains






def replace_TmpFactor_Rate4Site_Scores(chain, pdb_file, gradesPE, pdb_file_with_score_at_TempFactor):

    # This will create a PDB file that contains the Rate4Site scores instead of the TempFactor Column

    Rate4Site_Grades = {}
    LOG.write("Calling:cp_rasmol_gradesPE_and_pipe.read_Rate4Site_gradesPE(%s, %s)\n" %(gradesPE, str(Rate4Site_Grades)))
    ans = cp_rasmol_gradesPE_and_pipe.read_Rate4Site_gradesPE(gradesPE, Rate4Site_Grades)
    if ans != "OK":

        exit_on_error("sys_error", ans)

    LOG.write("Calling:cp_rasmol_gradesPE_and_pipe.replace_tempFactor(%s, %s, %s, %s)\n" %(pdb_file, chain, "Rate4Site_Grades", pdb_file_with_score_at_TempFactor))
    ans = cp_rasmol_gradesPE_and_pipe.replace_tempFactor(pdb_file, chain, Rate4Site_Grades, pdb_file_with_score_at_TempFactor)
    if ans != "OK":

        exit_on_error("sys_error", ans)



def replace_TmpFactor_Consurf_Scores(chain, pdb_file, gradesPE, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd):

    # This Will create a File containing the ATOMS records with the ConSurf grades instead of the TempFactor column

    LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.ReplaceTempFactConSurfScores(%s, %s, %s, %s, %s);\n" %(chain, pdb_file, gradesPE, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd))
    ans = cp_rasmol_gradesPE_and_pipe.ReplaceTempFactConSurfScore(chain, pdb_file, gradesPE, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd)
    if ans[0] != "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.ReplaceTempFactConSurf FAILED: " + ans)

    elif ans[1] == "yes":

        # insufficient data
        #vars['zip_list'].append(ATOMS_with_ConSurf_Scores_isd)
        return(ATOMS_with_ConSurf_Scores_isd)

    else:

        # no insufficient
        return(ATOMS_with_ConSurf_Scores)






def create_gradesPE_ConSurf(gradesPE, ref_r4s2pdb, gradesPE_Output_ArrRef, residue_freq_HashRef, no_isd_residue_color_ArrRef ,isd_residue_color_ArrRef, Type, Pos_Solv_Acc_Pred_NACSES_hashRef):

    ans = cp_rasmol_gradesPE_and_pipe.create_gradesPE_ConSurf(gradesPE_Output_ArrRef, ref_r4s2pdb, residue_freq_HashRef, no_isd_residue_color_ArrRef, isd_residue_color_ArrRef, gradesPE, Type, Pos_Solv_Acc_Pred_NACSES_hashRef, form['ALGORITHM'], vars['layers_array'])

    if not ans[0] == "OK":

        exit_on_error('sys_error', "create_gradesPE_ConSurf : cp_rasmol_gradesPE_and_pipe." + str(ans))

    elif not os.path.exists(gradesPE) or os.path.getsize(gradesPE) == 0:

        exit_on_error('sys_error', "create_gradesPE_ConSurf : the file '" + gradesPE + "' was not found or empty")

    if ans[1] == "" or ans[2] == "":

        exit_on_error('sys_error', "create_gradesPE_ConSurf : there is no data in the returned values seq3d_grades_isd or seq3d_grades from the routine")

    return(ans[1], ans[2])



def match_pdb_to_seq(ref_r4s2pdb, chain, query_seq, pdbseq, atom_positionFILE, Type):

    LOG.write("match_pdb_to_seq : calling cp_rasmol_gradesPE_and_pipe.match_seqres_pdb(%s ,%s, %s, %s, %s, %s)\n" %(pdbseq, query_seq, atom_positionFILE, chain, ref_r4s2pdb, Type))
    ans = cp_rasmol_gradesPE_and_pipe.match_seqres_pdb(pdbseq, query_seq, atom_positionFILE, chain, ref_r4s2pdb, Type)
    if not ans[0] == "OK":

        exit_on_error('sys_error', "match_pdb_to_seq : rasmol_gradesPE_and_pipe.%s" %ans[0])

    elif len(ref_r4s2pdb) < 1:

        LOG.write("match_pdb_to_seq : Total residues in the msa sequence: %s. Total residues in the ATOM : %s\n" %(ans[1], ans[2]))

    return(ans[1], ans[2])



def create_atom_position_file(chain, pdb_file, atom_positionFILE, Type, isPisa):

    if chain == "NONE":

        chain = " "

    output = {}
    LOG.write("create_atom_position_file : calling cp_rasmol_gradesPE_and_pipe.create_atom_position_file(%s, %s, %s, %s, %s)\n" %(pdb_file, atom_positionFILE, chain, str(output), Type))
    cp_rasmol_gradesPE_and_pipe.create_atom_position_file(pdb_file, atom_positionFILE, chain, output, Type)

    if 'INFO' in output:

        LOG.write("create_atom_position_file : " + output['INFO'])

    if 'WARNING' in output:

        LOG.write("create_atom_position_file : " + output['WARNING'])

    if 'ERROR' in output:

        if isPisa == 0:

            exit_on_error('sys_error', output['ERROR'])

        else:

            LOG.write("create_atom_position_file : %s\n" %output['ERROR'])

    if not os.path.exists(atom_positionFILE) or os.path.getsize(atom_positionFILE) == 0:

        if isPisa == 0:

            exit_on_error('sys_error', "create_atom_position_file : The file '" + atom_positionFILE + "' does not exist or of size 0\n")


        else:

            LOG.write("create_atom_position_file : The file '" + atom_positionFILE + "' does not exist or of size 0\n")

    if 'ERROR' in output:

        return 0

    else:

        return 1




def print_residue_precentage():

    # print a file that details percentage of each AA in the MSA
    LOG.write("cp_rasmol_gradesPE_and_pipe.print_precentage(%s, %s, %s, %s)\n" %("residue_freq", "position_totalAA", vars['Msa_percentageFILE'], "gradesPE_Output"))
    ans = cp_rasmol_gradesPE_and_pipe.print_precentage(residue_freq, position_totalAA, vars['Msa_percentageFILE'], gradesPE_Output)
    if ans[0] != "OK":

        exit_on_error('sys_error',str(ans))

    elif not os.path.exists(vars['Msa_percentageFILE']) or os.path.getsize(vars['Msa_percentageFILE']) == 0:

        exit_on_error('sys_error', "The output " + vars['Msa_percentageFILE'] + " was not found or empty")




def read_residue_variety(ref_to_res_freq, ref_to_positionAA):

    protein_MSA = vars['msa_fasta']
    query_string = vars['query_string']
    msa_format = "fasta"

    LOG.write("read_residue_variety : Calling: MSA_parser.read_residue_variety(%s, %s, %s, %s, %s)\n" %(protein_MSA, query_string, msa_format, ref_to_res_freq, ref_to_positionAA))
    ans = MSA_parser.read_residue_variety(protein_MSA, query_string, msa_format, ref_to_res_freq, ref_to_positionAA)
    if ans[0] != "OK":

        exit_on_error('sys_error', ans[0])

    if len(ref_to_res_freq) < 1 or len(ref_to_positionAA) < 1:

        exit_on_error('sys_error', "could not extract information from MSA '" + protein_MSA + "' in routine MSA_parser.read_residue_variety")

    return ans[1]




def assign_colors_according_to_r4s_layers(ref_to_gradesPE, r4s_out):

    LOG.write("assign_colors_according_to_r4s_layers : %s\n" %r4s_out)
    ans = cp_rasmol_gradesPE_and_pipe.assign_colors_according_to_r4s_layers(ref_to_gradesPE, r4s_out, form['ALGORITHM'])
    if ans[0] != "OK":

        exit_on_error('sys_error', ans)

    else:

        vars['layers_array'] = ans[1]
        LOG.write("assign_colors_according_to_r4s_layers : color layers are %s\n" %str(vars['layers_array']))



def run_rate4site():

    rate4s = GENERAL_CONSTANTS.RATE4SITE
    rate4s_ML = GENERAL_CONSTANTS.RATE4SITE_ML
    rate4s_slow = GENERAL_CONSTANTS.RATE4SITE_SLOW

    algorithm = ""
    tree_file_r4s = ""
    msa = ""
    did_r4s_fail = ""

    MatrixHash = {'JTT' : '-Mj', 'MTREV' : '-Mr', 'CPREV' : '-Mc', 'WAG' : '-Mw', 'DAYHOFF' : '-Md', 'T92' : '-Mt', 'HKY' : '-Mh', 'GTR' : '-Mg', 'JC_NUC' : '-Mn', 'JC_AA' : '-Ma', 'LG' : '-Ml'}

    r4s_comm = ""
    comm_end = ""
    if form['ALGORITHM'] == "Bayes":

        algorithm = "-ib" # Save the algorithm, it maybe used later
        r4s_comm += "%s -ib -a %s -s %s -zn %s " %(rate4s, vars['query_string'], vars['msa_fasta'], MatrixHash[(form['SUB_MATRIX']).upper()]) 
        comm_end += "-bn -l %s -o %s -n 32 -v 9 " %( vars['r4s_log'], vars['r4s_out'])

    else:

        algorithm = "-im" # Save the algorithm, it maybe used later
        r4s_comm += "%s %s -a %s -s %s -zn %s " %(rate4s_ML, algorithm, vars['query_string'], vars['msa_fasta'], MatrixHash[(form['SUB_MATRIX']).upper()])
        comm_end += "-bn -l %s -o %s -v 9 " %(vars['r4s_log'], vars['r4s_out']) 

    if vars['running_mode'] == "_mode_pdb_msa_tree" or vars['running_mode'] == "_mode_msa_tree":

        r4s_comm += "-t %s " %vars['tree_file']

    r4s_comm += comm_end

    LOG.write("run_rate4site : running command: %s\n" %r4s_comm)
    submit_job_to_Q("rate4site", r4s_comm)

    # if the run failed - we rerun using the slow verion
    if check_if_rate4site_failed(vars['r4s_log']):

        LOG.write("run_rate4site : The run of rate4site failed. Sending warning message to output.\nThe same run will be done using the SLOW version of rate4site.\n")
        #print_message_to_output("<font color='red'><b>Warning:</b></font> The given MSA is very large, therefore it will take longer for ConSurf calculation to finish. The results will be sent to the e-mail address provided.<br>The calculation continues nevertheless.")
        form['send_user_mail'] = "yes"
        r4s_comm = "%s %s -a %s -s %s -zn %s %s -bn -l %s -o %s -n 32 -v 9" %(rate4s_slow, algorithm, vars['query_string'], vars['msa_fasta'],  MatrixHash[(form['SUB_MATRIX']).upper()], tree_file_r4s, vars['r4s_slow_log'], vars['r4s_out'])
        LOG.write("run_rate4site : running command: %s\n" %str(r4s_comm))
        submit_job_to_Q("rate4siteSlow", r4s_comm)

        if check_if_rate4site_failed(vars['r4s_slow_log']):

            exit_on_error('user_error', "The calculation could not be completed due to memory problem, since the MSA is too large. Please run ConSurf again with fewer sequences.")

def check_if_rate4site_failed(r4s_log):

    # There are some tests to see if rate4site failed.
    # Since I can't trust only one of them, I do all of them. If onw of them is tested to be true - than a flag will get TRUE value
    # 1. the .res file might be empty.
    # 2. if the run failed, it might be written to the log file of r4s.
    # 3. in a normal the r4s.log file there will lines that describe the grades. if it fail - we won't see them
    # In one of these cases we try to run the slower version of rate4site.
    # We output this as a message to the user.
    ret = False
    did_r4s_failed = rate4site_routines.check_if_rate4site_failed(vars['r4s_out'], r4s_log)
    if did_r4s_failed[0] == "yes":

        ret = True
        LOG.write("check_if_rate4site_failed : " + did_r4s_failed[1])
        remove_core(did_r4s_failed[2])
        # if there was an error in user input which rate4site reported: we output a message to the user and exit
        if len(did_r4s_failed) > 3 and did_r4s_failed[3] != "":

            exit_on_error('user_error', did_r4s_failed[3])

    return ret



def find_best_substitution_model():

    # convert fasta to phylip
    msa_phy_filepath = "input_msa.phy"
    convert_msa_format(vars['msa_fasta'], "fasta", msa_phy_filepath, "phylip-relaxed")
    model = ""

    if form['DNA_AA'] == "Nuc":

        model = run_jmt(msa_phy_filepath)

    else:

        model = run_protest(msa_phy_filepath)

    model = model.strip('()')
    if model == "JC" and form['DNA_AA'] == "Nuc":

        model = "JC_Nuc"

    return model

def run_protest(msa_file_path):

    PRT_JAR_FILE = "/bioseq/Programs/ModelTest/prottest-3.4.1/prottest-3.4.1.jar"
    output_file_path = "model_selection.txt"
    cmd = "java -jar %s -log disabled -i %s -AICC -o %s -S 1 -JTT -LG -MtREV -Dayhoff -WAG -CpREV -threads 1" %(PRT_JAR_FILE, msa_file_path, output_file_path)
    submit_job_to_Q("protest", cmd)
    LOG.write("run_protest: %s\n" %cmd)

    try:

        f = open(output_file_path, 'r')

    except:

        exit_on_error('sys_error', "Cannot open the file " + output_file_path + " for reading.")

    line = f.read()
    model = re.search(r"(?<=Best model according to AICc: ).*", line).group()
    f.close()

    return model

def run_jmt(msa_file_path):

    JMT_JAR_FILE = "/bioseq/Programs/ModelTest/jmodeltest-2.1.7/jModelTest.jar"
    output_file_path = "model_selection.txt"
    cmd = "java -jar %s -d %s -t BIONJ -AICc -f -o %s" %(JMT_JAR_FILE, msa_file_path, output_file_path)
    submit_job_to_Q("jmt", cmd)
    LOG.write("run_jmt: %s\n" %cmd)

    try:

        f = open(output_file_path, 'r')

    except:

        exit_on_error('sys_error', "Cannot open the file " + output_file_path + " for reading.")

    start_reading = 0
    test_table_headers = "Model\s+-lnL\s+K\s+AICc\s+delta\s+weight\s+cumWeight(\s+uDelta)?"
    JMT_VALID_MODELS = ["JC","HKY","GTR"]

    # extract best model from table
    line = f.readline()
    while line != "":

        if start_reading:

            split_row = line.strip().split()
            if split_row[0] in JMT_VALID_MODELS:

                f.close()
                return split_row[0]

        elif re.search(test_table_headers, line, re.M):

            start_reading = 1

        line = f.readline()

    f.close()
    return ""

def convert_msa_format(infile, infileformat, outfile, outfileformat):

    try:

        AlignIO.convert(infile, infileformat, outfile, outfileformat)

    except:

        exit_on_error('sys_error', "convert_msa_format : exception")

    if not os.path.exists(outfile) or os.path.getsize(outfile) == 0:

        exit_on_error('sys_error', "convert_msa_format : MSA not created.")


def create_MSA():

    if form['MSAprogram'] == "CLUSTALW":

        cmd = "%s -infile=%s -outfile=%s" %(GENERAL_CONSTANTS.CLUSTALW, vars['FINAL_sequences'], vars['msa_clustal'])
        LOG.write("create_MSA : run %s\n" %cmd)
        submit_job_to_Q("clustalw", cmd)
        convert_msa_format(vars['msa_clustal'], "clustal", vars['msa_fasta'], "fasta")

    elif form['MSAprogram'] == "MAFFT":

        cmd = "%s --localpair --maxiterate 1000 --quiet %s > %s" %(GENERAL_CONSTANTS.MAFFT_LINSI_GUIDANCE, vars['FINAL_sequences'], vars['msa_fasta'])
        LOG.write("create_MSA : run %s\n" %cmd)
        submit_job_to_Q("MAFFT", cmd)
        convert_msa_format(vars['msa_fasta'], "fasta", vars['msa_clustal'], "clustal")

    elif form['MSAprogram'] == "PRANK":

        cmd = "%s -d=%s -o=%s -F" %(GENERAL_CONSTANTS.PRANK, vars['FINAL_sequences'], vars['msa_fasta'])
        #print_message_to_output("<font color='red'><b>Warning:</b></font> PRANK is accurate but slow MSA program, please be patient.")
        LOG.write("create_MSA : run %s\n" %cmd)
        submit_job_to_Q("PRANK", cmd)

        if os.path.exists(vars['msa_fasta'] + ".2.fas"):

            vars['msa_fasta'] += ".2.fas"

        elif os.path.exists(vars['msa_fasta'] + ".1.fas"):

            vars['msa_fasta'] += ".1.fas"

        elif os.path.exists(vars['msa_fasta'] + ".best.fas"):

            vars['msa_fasta'] +=  ".best.fas"

        convert_msa_format(vars['msa_fasta'], "fasta", vars['msa_clustal'], "clustal")

    else:

        cmd = "%s -in %s -out %s -clwstrict -quiet" %(GENERAL_CONSTANTS.MUSCLE, vars['FINAL_sequences'], vars['msa_clustal'])
        LOG.write("create_MSA : run %s\n" %cmd)
        submit_job_to_Q("MUSCLE", cmd)
        convert_msa_format(vars['msa_clustal'], "clustal", vars['msa_fasta'], "fasta")


def choose_final_homologoues(ref_cd_hit_hash, ref_blast_hash):

    size_cd_hit_hash = len(ref_cd_hit_hash)
    size_blast_hash = len(ref_blast_hash)

    LOG.write("In sub choose_final_homologoues: size cd_hit hash: " + str(size_cd_hit_hash) + ", size blast hash: " + str(size_blast_hash) + "\n")

    try:

        FINAL = open(vars['FINAL_sequences'], 'w')

    except:

        exit_on_error('sys_error',"choose_final_homologoues : cannot open the file %s for writing" %vars['FINAL_sequences'])

    FINAL.write(">%s\n%s\n" %(vars['query_string'], vars['protein_seq_string']))
    FINAL.close()

    final_file_size = os.path.getsize(vars['FINAL_sequences']) # take the size of the file before we add more sequences to it
    LOG.write("parseFiles.sort_sequences_from_eval(%s ,%s , %f, %s, %s)\n" %("ref_blast_hash" ,"ref_cd_hit_hash", float(form['MAX_NUM_HOMOL']) -1, form['best_uniform_sequences'], vars['FINAL_sequences']))
    ans = parseFiles.sort_sequences_from_eval(ref_blast_hash ,ref_cd_hit_hash, float(form['MAX_NUM_HOMOL']) -1, form['best_uniform_sequences'], vars['FINAL_sequences'])

    if ans[0] == "err":

        exit_on_error('sys_error', ans[1])

    LOG.write(ans[1] + "\n")
    vars['final_number_of_homologoues'] = ans[2]

    # check that more sequences were added to the file
    if not final_file_size < os.path.getsize(vars['FINAL_sequences']):

        exit_on_error('sys_error', "choose_final_homologoues : the file " + vars['FINAL_sequences'] + " doesn't contain sequences")


def add_sequences_removed_by_cd_hit_to_rejected_report():

    LOG.write("add_sequences_removed_by_cd_hit_to_rejected_report : running parseFiles.add_sequences_removed_by_cd_hit_to_rejected_report(%s, %s)\n" %(vars['cd_hit_out_file'] + ".clstr", vars['HITS_rejected_file']))

    ans = parseFiles.add_sequences_removed_by_cd_hit_to_rejected_report(vars['cd_hit_out_file'] + ".clstr", vars['HITS_rejected_file'])

    if ans != "ok":

        exit_on_error('sys_error', ans)

def print_message_to_output(msg):

    try:

        OUTPUT = open(vars['output_page'], 'a')

    except:

        exit_on_error('sys_error', "print_message_to_output : could not open the file " + vars['output_page'] + " for writing.")

    OUTPUT.write("<br><br>" + msg + "\n")
    OUTPUT.close()

def cluster_homologoues(ref_cd_hit_hash):

    msg = ""
    #LOG.write("cluster_homologoues : running parseFiles.create_cd_hit_output(%s, %s, %f, %s, %s, %s);\n" %(vars['HITS_fasta_file'], vars['cd_hit_out_file'], vars['hit_redundancy']/100, GENERAL_CONSTANTS.CD_HIT_DIR, ref_cd_hit_hash, form['DNA_AA']))
    LOG.write("cluster_homologoues : create_cd_hit_output(%s, %s, %f, %s, %s, %s);\n" %(vars['HITS_fasta_file'], vars['cd_hit_out_file'], vars['hit_redundancy']/100, GENERAL_CONSTANTS.CD_HIT_DIR, ref_cd_hit_hash, form['DNA_AA']))
    #ans = parseFiles.create_cd_hit_output(vars['HITS_fasta_file'], vars['cd_hit_out_file'], vars['hit_redundancy']/100, GENERAL_CONSTANTS.CD_HIT_DIR, ref_cd_hit_hash, form['DNA_AA'])
    ans = create_cd_hit_output(vars['HITS_fasta_file'], vars['cd_hit_out_file'], vars['hit_redundancy']/100, GENERAL_CONSTANTS.CD_HIT_DIR, ref_cd_hit_hash, form['DNA_AA'])

    if ans[0] == "sys":

        exit_on_error('sys_error', ans[1])

    total_num_of_hits = ans[1]

    if form['MAX_NUM_HOMOL'] == 'all' or form['MAX_NUM_HOMOL'] == 'ALL':

        form['MAX_NUM_HOMOL'] == total_num_of_hits

    if total_num_of_hits < vars['min_num_of_hits']: # less seqs than the minimum: exit

        if total_num_of_hits <= 1:

            msg = "There is only 1 "

        else:

            msg = "There are only %d " %total_num_of_hits

        msg += "<a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\">unique %s</a> hits. The minimal number of sequences required for the calculation is %d. You may try to:<ol>" %(vars['cd_hit_out_file'], vars['blast_algorithm'], vars['min_num_of_hits'])

        if form['proteins_DB'] == "SWISS-PROT":

            msg += "<li> Run your query using UniProt, UniRef90 or NR database. "

        msg += "<li>Re-run the server and manually select the homologous sequences.<li>Re-run the server with a multiple sequence alignment file of your own.</li><li>Increase the Evalue.</li></li><li>Decrease the Minimal %ID For Homologs"

        if int(form['ITERATIONS']) < 5:

            msg += " <li> Increase the number of " + vars['blast_algorithm'] + " iterations.</li>"

        msg += "</OL>\n"
        exit_on_error('user_error',msg)

    elif total_num_of_hits + 1 < vars['low_num_of_hits']: # less seqs than 10 : output a warning.

        msg = "<font color='red'><b>Warning:</font></b> There are "

        if total_num_of_hits + 1 < vars['number_of_homologoues_before_cd-hit']: # because we will add the query sequence itself to all the unique sequences.

            msg += "%d <a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\" download>%s</a> hits, only %d of them are" %(vars['number_of_homologoues_before_cd-hit'], vars['BLAST_out_file'], vars['blast_algorithm'], total_num_of_hits+1)

        else:

            msg += str(total_num_of_hits + 1)

        msg += " unique sequences. The calculation is performed on the %d <a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\">unique sequences</a>, but it is recommended to run the server with a multiple sequence alignment file containing at least %s sequences." %(total_num_of_hits + 1, vars['FINAL_sequences_html'], vars['low_num_of_hits'])

    else:

        msg = "There are <a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\" download>%d %s hits</a>. %d of them are unique, including the query.<br/>The calculation is performed on " %(vars['BLAST_out_file'], vars['number_of_homologoues_before_cd-hit'], vars['blast_algorithm'], total_num_of_hits + 1)
        #msg = "There are %d <a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\" download>%s</a> hits. %d of them are unique, including the query.<br/>The calculation is performed on " %(vars['number_of_homologoues_before_cd-hit'], vars['BLAST_out_file'], vars['blast_algorithm'], total_num_of_hits + 1)

        if total_num_of_hits <= int(form['MAX_NUM_HOMOL']):

            msg += "%d <a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\">unique sequences</a>." %(total_num_of_hits + 1, vars['FINAL_sequences_html'])

        elif form['best_uniform_sequences'] == "best":

            msg += "the %s <a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\">sequences</a> closest to the query (with the lowest E-value)." %(form['MAX_NUM_HOMOL'], vars['FINAL_sequences_html'])

        else:

            msg += "a sample of <a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\">%s sequences</a> that represent the list of homologues to the query." %(vars['FINAL_sequences_html'], form['MAX_NUM_HOMOL'])
            #msg += "a sample of %s <a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\">sequences</a> that represent the list of homologues to the query." %(form['MAX_NUM_HOMOL'], vars['FINAL_sequences_html'])

    """
    print_message_to_output(msg)

    if os.path.exists(vars['HITS_rejected_file']) and os.path.getsize(vars['HITS_rejected_file']) != 0:

        print_message_to_output("Here is the <a href=\"<?=$orig_path?>/" + vars['HITS_rejected_file'] + "\" TARGET=Rejected_Seqs style=\"color: #400080; text-decoration:underline;\">list of sequences</a> that produced significant alignments, but were not chosen as hits.")
        #print_message_to_output("<a href=\"<?=$orig_path?>/" + vars['HITS_rejected_file'] + "\" TARGET=Rejected_Seqs style=\"color: #400080; text-decoration:underline;\">Click here</a> if you wish to view the list of sequences which produced significant alignments, but were not chosen as hits.")
    """
    return(total_num_of_hits + 1)


def choose_homologoues_from_search_with_lower_identity_cutoff(ref_blast_hash):

    LOG.write("Gershon_updated CGI:  choose_homologoues_from_search : running parseFiles::choose_homologoues_from_search_with_lower_identity_cutoff\n")
    LOG.write("parseFiles::choose_homologoues_from_search_with_lower_identity_cutoff(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);\n" %(form['Homolog_search_algorithm'], vars['protein_seq_string'], vars['hit_redundancy'], vars['hit_overlap'], vars['hit_min_length'], form['MIN_IDENTITY'], vars['min_num_of_hits'], vars['BLAST_out_file'], vars['HITS_fasta_file'], vars['HITS_rejected_file'], str(ref_blast_hash), form['DNA_AA']))
    LOG.write("before\n")
    ans = parseFiles.choose_homologoues_from_search_with_lower_identity_cutoff(form['Homolog_search_algorithm'], vars['protein_seq_string'], vars['hit_redundancy'], vars['hit_overlap'], vars['hit_min_length'], float(form['MIN_IDENTITY']), vars['min_num_of_hits'], vars['BLAST_out_file'], vars['HITS_fasta_file'], vars['HITS_rejected_file'], ref_blast_hash, form['DNA_AA'])
    LOG.write("after\n")
    LOG.write("ans[0]: %s, ans[1]: %s\n" %(ans[0], ans[1]))

    if ans[0] == "sys":

        exit_on_error('sys_error', ans[1])

    elif ans[0] == "user":

        message = "According to the parameters of this run, " + ans[1] + " You can try to:<OL>"

        if form['proteins_DB'] == "SWISS-PROT":

            message += "<li> Run your query using UniProt, UniRef90 or NR database. "

        message += "<li>Re-run the server and manually select the homologous sequences.<li>Re-run the server with a multiple sequence alignment file of your own.</li><li>Increase the Evalue.</li></li><li>Decrease the Minimal %ID For Homologs"

        if int(form['ITERATIONS']) < 5:

            message += " <li> Increase the number of " + vars['blast_algorithm'] + " iterations.</li>"
            if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_pdb_msa_tree":
			
                message += " <li> Try the web server <a href=\"https://evorator.tau.ac.il/\" style=\"color: #400080; text-decoration:underline;\">EvoRator</a>.</li>"

        message += "</OL>\n"
        exit_on_error('user_error', message)

    if not os.path.exists(vars['HITS_fasta_file']) or os.path.getsize(vars['HITS_fasta_file']) == 0:

        exit_on_error('sys_error',"choose_homologoues_from_search_with_lower_identity_cutoff : the file " + vars['HITS_fasta_file'] + " was not created or contains no data")

    vars['number_of_homologoues'] = ans[1]
    vars['number_of_homologoues_before_cd-hit'] = ans[2]

def run_search(Search_Out_File, Output_Type = "PlainText"):

    # Search_Out_File - The output File
    # Output_Type - The output can be in HTML format or xml

    cmd = ""

    if form['DNA_AA'] == "AA":

        if form['Homolog_search_algorithm'] == "BLAST":

            cmd = "%s -m 7 -i %s -e %s -d %s -j %s -J T -v %s -b %s -o %s -F F" %(GENERAL_CONSTANTS.BLASTPGP, vars['protein_seq'], form['E_VALUE'], vars['protein_db'], form['ITERATIONS'], vars['max_homologues_to_display'], vars['max_homologues_to_display'], Search_Out_File)

        elif form['Homolog_search_algorithm'] == "HMMER":

            cmd = "jackhmmer --notextw -N %s --domE %s -E %s --incE %s --cpu 1 %s  %s > %s" %(form['ITERATIONS'], form['E_VALUE'], form['E_VALUE'], form['E_VALUE'],  vars['protein_seq'], vars['protein_db'], Search_Out_File)

        else: # if form['Homolog_search_algorithm'] == "CS_BLAST"

            cmd = "%s -i %s -e %s -d %s -D %s/K4000.lib -j %s -v %s -b %s -o %s -m 7 -F F --blast-path %s" %(GENERAL_CONSTANTS.CS_BLAST, vars['protein_seq'], form['E_VALUE'], vars['protein_db'], GENERAL_CONSTANTS.CS_BLAST_DATA, form['ITERATIONS'], vars['max_homologues_to_display'], vars['max_homologues_to_display'], Search_Out_File, GENERAL_CONSTANTS.BLAST_PATH)

    else:

        if form['Homolog_search_algorithm'] == "HMMER":

            cmd = "nhmmer --notextw --incE %s -E %s --cpu 6 %s %s > %s" %(form['E_VALUE'], form['E_VALUE'], vars['protein_seq'], GENERAL_CONSTANTS.NR_NUC_DB_FASTA, Search_Out_File)

        else:

            cmd = "%s -p blastn -m 7 -i %s -e %s -d %s -v %s -b %s -o %s -F F" %(GENERAL_CONSTANTS.BLASTALL, vars['protein_seq'], form['E_VALUE'], GENERAL_CONSTANTS.NR_NUC_DB, vars['max_homologues_to_display'], vars['max_homologues_to_display'], Search_Out_File)

    if form['Homolog_search_algorithm'] != "HMMER":

        if Output_Type == "HTML":

            cmd += " -T T\n" # Blast Result will be formated to HTML

        else:

            cmd += " -T F\n" # Plain Text Format

    else:

        cmd += "\n"

    LOG.write("run_search : running: " + cmd + "\n")
    submit_job_to_Q("BLAST", cmd)

    if not os.path.exists(Search_Out_File) or os.path.getsize(Search_Out_File) == 0:

        exit_on_error('sys_error',"run_search : run of search fail. " + Search_Out_File + " is zero or not exists")

def submit_job_to_Q(job_name_prefix, cmd):
    
    #q_cmd = ["module", "load", "python/python-3.8", "hmmr/hmmr-3.1b2", "clustalw/2.1;"]
    #q_cmd = cmd.split()
    LOG.write("submit_job_to_Q : command: %s\n" %cmd)
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8", shell=True)
    out, err = p.communicate()
    LOG.write("waiting for the file %s\n" %err)
    
"""
def submit_job_to_Q(job_name_prefix, cmd, load_module = ""):

    qsub_script = "qsub_" + job_name_prefix + ".sh"; # script to run in the queue
    #job_end_flag = job_name_prefix + "_ENDS"
    #cmd += "\necho > %s\n" %job_end_flag
    err_file = job_name_prefix + ".err"
    out_file = job_name_prefix + ".out"

    try:

        QSUB_SH = open(qsub_script, 'w')

    except:

        exit_on_error('sys_error', "submit_job_to_Q : cannot open the file " + qsub_script + " for writing!")

    QSUB_SH.write("#!/bin/bash\n")
    QSUB_SH.write("#PBS -N %s_%s\n"  %(form['Run_Number'], job_name_prefix))
    QSUB_SH.write("#PBS -l nodes=1:ppn=3\n")
    QSUB_SH.write("#PBS -r y\n")
    QSUB_SH.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n")
    QSUB_SH.write("#PBS -e %s\n"  %(vars['working_dir'] + err_file))
    QSUB_SH.write("#PBS -o %s\n"  %(vars['working_dir'] + out_file))
    QSUB_SH.write("cd %s\n"  %vars['working_dir'])
    if load_module == "":

        QSUB_SH.write("module load python/python-3.8 hmmr/hmmr-3.1b2 clustalw/2.1\n")

    QSUB_SH.write(cmd)
    QSUB_SH.close()

    q_cmd = ['ssh', 'bioseq@power', 'qsub', '-q', 'bentalweb', '-l', 'nodes=1:ppn=2', vars['working_dir'] + qsub_script] # HA comment
    LOG.write("\nsubmit_job_to_Q :\n['ssh', 'bioseq@powerweb1', 'qsub', '-q', 'bentalweb', '-l', 'nodes=1:ppn=2', '" + qsub_script + "']\n")
    p = subprocess.Popen(q_cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    p.communicate()
    #LOG.write("waiting for the file %s\n" %job_end_flag)
    #start_job = time.time()
    LOG.write("waiting for the file %s\n" %err_file)
    while True:

        #if os.path.exists(job_end_flag):
        if os.path.exists(err_file):

            break

        time.sleep(1)

    #end_job = time.time()
    #JOB_TIME = open("/bioseq/consurf-rel/" + job_name_prefix + ".txt", 'a')
    #JOB_TIME.write(str(end_job - start_job) + "\n")
    #JOB_TIME.close()

    try:
	
        ERROR_FILE = open(err_file, 'r')
		
    except:
	
        exit_on_error('sys_error', "submit_job_to_Q : cannot open the file " + err_file + " for reading")
		
    error_message = ERROR_FILE.read()
    ERROR_FILE.close()
    if job_name_prefix == "CD-HIT" and "not enough memory, please set -M option greater than " in error_message:
        
        match = re.search(r'not enough memory, please set -M option greater than (\S+)', error_message)
        if match:
            
            limit = str(int(match.group(1)) + 100)
            submit_job_to_Q(job_name_prefix + "_increase_M", cmd + " -M " + limit)
        
    elif "exceeded limit" in error_message and "Terminated" in error_message:

        print_message_to_output("Some processes are being moved to a low priority queue. The calculation may take a longer time.")
        submit_job_to_Q_low_priority(job_name_prefix + "_LOW_PRIORITY", cmd, load_module)

	

def submit_job_to_Q_low_priority(job_name_prefix, cmd, load_module):

    qsub_script = "qsub_" + job_name_prefix + ".sh"; # script to run in the queue
    #job_end_flag = job_name_prefix + "_ENDS"
    #cmd += "\necho > %s\n" %job_end_flag
    err_file = job_name_prefix + ".err"
    out_file = job_name_prefix + ".out"
	
    try:

        QSUB_SH = open(qsub_script, 'w')

    except:

        exit_on_error('sys_error', "submit_job_to_Q_low_priority : cannot open the file " + qsub_script + " for writing!")

    QSUB_SH.write("#!/bin/bash\n")
    QSUB_SH.write("#PBS -N %s_%s\n"  %(form['Run_Number'], job_name_prefix))
    QSUB_SH.write("#PBS -l nodes=1:ppn=3:ncpus=4\n")
    QSUB_SH.write("#PBS -r y\n")
    QSUB_SH.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n")
    QSUB_SH.write("#PBS -e %s\n"  %(vars['working_dir'] + err_file))
    QSUB_SH.write("#PBS -o %s\n"  %(vars['working_dir'] + out_file))
    QSUB_SH.write("cd %s\n"  %vars['working_dir'])
    if load_module == "":

        QSUB_SH.write("module load python/python-3.8 hmmr/hmmr-3.1b2 clustalw/2.1\n")

    QSUB_SH.write(cmd)
    QSUB_SH.close()

    q_cmd = ['ssh', 'bioseq@power', 'qsub', '-q', 'bentalweb', vars['working_dir'] + qsub_script] # HA comment
    LOG.write("\nsubmit_job_to_Q_low_priority :\n['ssh', 'bioseq@powerweb1', 'qsub', '-q', 'bentalweb', '-l', 'nodes=1:ppn=2', '" + qsub_script + "']\n")
    p = subprocess.Popen(q_cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    p.communicate()
    #start_job = time.time()
    LOG.write("waiting for the file %s\n" %err_file)
    while True:

        if os.path.exists(err_file):

            break

        time.sleep(1)

    #end_job = time.time()
    #JOB_TIME = open("/bioseq/consurf-rel/" + job_name_prefix + ".txt", 'a')
    #JOB_TIME.write(str(end_job - start_job) + "\n")
    #JOB_TIME.close()
"""


def compare_atom_seqres_or_msa(what_to_compare):

    # in case they are both seqres and atom fields, checks the similarity between the 2 sequences.

    two_fastas = "PDB_" + what_to_compare + ".fasta2"
    clustalw_out = "PDB_" + what_to_compare + ".out"
    pairwise_aln = "PDB_" + what_to_compare + ".aln"
    atom_length = len(vars['ATOM_without_X_seq'])
    alignment_score = 0
    other_query_length = 0
    seqres_or_msa_seq = ""
    query_line = {}
    atom_line = "sequence extracted from the ATOM field of the PDB file"
    query_line['SEQRES'] = "sequence extracted from the SEQRES field of the PDB file"
    query_line['MSA'] = "sequence extracted from the MSA file"

    if what_to_compare == "SEQRES":

        other_query_length = len(vars['SEQRES_seq'])
        seqres_or_msa_seq = vars['SEQRES_seq']

    else:

        other_query_length = len(vars['MSA_query_seq'])
        seqres_or_msa_seq = vars['MSA_query_seq']

    """
    # compare the length of ATOM and SEQRES. output a message accordingly
    if other_query_length != 0 and other_query_length < atom_length:

        print_message_to_output("The %s is shorter than the %s.<br>The %s sequence has %f residues and the ATOM sequence has %f residues. The calculation continues nevertheless." %(query_line[what_to_compare],atom_line ,what_to_compare, other_query_length, atom_length))

    if atom_length < other_query_length:

        if atom_length < other_query_length * 0.2:

            print_message_to_output("<font color=\'red\'><b>Warning:</font></b> The %f is significantly shorter than the %s. The %s sequence has %f residues and the ATOM sequence has only %f residues. The calculation continues nevertheless." %(atom_line, query_line[what_to_compare], what_to_compare, other_query_length, atom_length))

        else:

            print_message_to_output("The %s is shorter than the %s. The %s sequence has %d residues and the ATOM sequence has %d residues. The calculation continues nevertheless." %(atom_line, query_line[what_to_compare], what_to_compare, other_query_length, atom_length))
    """
    # run clustalw to see the match between ATOM and SEQRES sequences
    LOG.write("compare_atom_seqres_or_msa : run clustalw to see the match between ATOM and " + what_to_compare + " sequences\n")
    [vars['seqres_or_msa_seq_with_gaps'], vars['ATOM_seq_with_gaps']] = compare_two_fastas(two_fastas, seqres_or_msa_seq, vars['ATOM_seq'], clustalw_out, pairwise_aln, what_to_compare)

    try:

        OUT = open(clustalw_out, 'r')

    except:

        exit_on_error('sys_error', "compare_atom_seqres_or_msa : Cannot open the file %s for reading." %clustalw_out)

    line = OUT.readline()
    while line != "":

        match = re.match(r'Sequences.+Aligned.+Score:\s+(\d+)', line)
        if match:

            alignment_score = int(match.group(1))
            break

        line = OUT.readline()

    OUT.close()

    if alignment_score < 100:

        if alignment_score < 30:

            exit_on_error('user_error',"The Score of the alignment between the %s and the %s is ONLY %d%% identity.<br>See <a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\" TARGET=PairWise_Align>pairwise alignment</a>." %(query_line[what_to_compare], atom_line, alignment_score, pairwise_aln))

        """
        else:

            print_message_to_output("The Score of the alignment between the %s and the %s is %d%% identity.<br>See <a href=\"<?=$orig_path?>/%s\" style=\"color: #400080; text-decoration:underline;\" TARGET=PairWise_Align>pairwise alignment</a>. The calculation continues nevertheless." %(query_line[what_to_compare], atom_line, alignment_score, pairwise_aln))
        """

def exit_on_error(which_error, error_msg):

    error_definition = "<font size=+2 color='red'>ERROR! ConSurf session has been terminated:</font><br />\n"
    syserror = "<font size=+1 color='red'>A SYSTEM ERROR OCCURRED!</font><br />Please try to run ConSurf again in a few minutes.<br />We apologize for the inconvenience.<br />\n"

    LOG.write("\n\t EXIT on error:\n" + error_msg + "\n")
    """
    if which_error == 'user_error':

        print_message_to_output(error_definition + error_msg)

    elif which_error == 'sys_error':

        print_message_to_output(syserror)
        send_administrator_mail_on_error(error_msg)
        
    email_subject = "'Your ConSurf run titled " + form['JOB_TITLE'] + " FAILED'"
    email_message = "'Hello,\\n\\nUnfortunately your ConSurf run (number %s) has failed.\\nPlease have a look at %s for further details\\n\\nSorry for the inconvenience\\nConSurf Team'" %(form['Run_Number'], vars['run_url'])
    send_email(email_subject, email_message)

    Failed()
    """
    LOG.write("\nExit Time: %s \n" %datetime.now())
    LOG.close()
    print("%s\n%s"  %(which_error, error_msg))
    exit(1)





def compare_two_fastas(two_fastas, first_seq ,second_seq ,clustalw_out, clustalw_aln, seq_type):

    vars['compare_two_fastas_count'] += 1
    clustalw = GENERAL_CONSTANTS.CLUSTALW

    try:

        FAS = open(two_fastas, 'w')

    except:

        exit_on_error('sys_error', "compare_two_fastas : Cannot open the file " + two_fastas + " for writing.")

    FAS.write(">%s_SEQ\n%s\n>ATOM_SEQ\n%s\n" %(seq_type, first_seq, second_seq))
    FAS.close()

    """
    command = ["ssh", "bioseq@power", "cd", vars['working_dir'] + ";", clustalw, two_fastas, "-gapopen=1", ">", clustalw_out]
    LOG.write("compare_two_fastas : run %s\n" %command)
    p = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    out, err = p.communicate()
    """

    command = "cd %s; %s %s -gapopen=1 > %s" %(vars['working_dir'], clustalw, two_fastas, clustalw_out)
    submit_job_to_Q("compare_two_fastas_%d" %vars['compare_two_fastas_count'], command)


    if not os.path.exists(clustalw_out) or os.path.getsize(clustalw_out) == 0 or not os.path.exists(clustalw_aln) or os.path.getsize(clustalw_aln) == 0:

        exit_on_error('sys_error', "compare_two_fastas : one of clustalw outputs were not created; %s or %s" %(clustalw_out, clustalw_aln))

    # return the two sequences with gaps

    try:

        CLUSTALW_ALN = open(clustalw_aln, 'r')

    except:

        exit_on_error('sys_error', "compare_two_fastas : could not open " + clustalw_out + " for reading.")

    first_seq_with_gaps = ""
    second_seq_with_gaps = ""

    line = CLUSTALW_ALN.readline()
    while line != "":

        match1 = re.match(r'^' + seq_type + r'_SEQ\s+(\S*)', line)
        if match1:

            first_seq_with_gaps += match1.group(1)

        else:

            match2 = re.match(r'^ATOM_SEQ\s+(\S*)', line)
            if match2:

                second_seq_with_gaps += match2.group(1)

        line = CLUSTALW_ALN.readline()


    CLUSTALW_ALN.close()

    return(first_seq_with_gaps, second_seq_with_gaps)




def analyse_seqres_atom():

    # there is no ATOM field in the PDB

    if vars['ATOM_without_X_seq'] == "":

        exit_on_error('user_error', "There is no ATOM derived information in the PDB file.<br>Please refer to the OVERVIEW for detailed information about the PDB format.")

    # there is no SEQRES field in the PDB

    if vars['SEQRES_seq'] == "":

        msg = "<font color='red'>Warning:</font> There is no SEQRES derived information in the PDB file. The calculation will be based on the ATOM derived sequence. "

        if vars['running_mode'] == "_mode_pdb_no_msa":

            msg += "If this sequence is incomplete, we recommend to re-run the server using an external multiple sequence alignment file, which is based on the complete protein sequence."

        LOG.write("analyse_seqres_atom : There is no SEQRES derived information in the PDB file.\n")
        #print_message_to_output(msg)

    if form['DNA_AA'] == "AA":

        # check if seqres contains nucleic acid
        #type_SEQRES = vars['pdb_object'].get_type_SEQRES()
        #if form['PDB_chain'] in type_SEQRES and type_SEQRES[form['PDB_chain']] == "Nuc":
        if vars['pdb_object'].get_type() == "Nuc":

            exit_on_error('user_error', "The selected chain: " + form['PDB_chain'] + " contains nucleic acid, and you have selected amino acid")

    else:

        # check if seqres contains amino acid
        #type_SEQRES = vars['pdb_object'].get_type_SEQRES()
        #if form['PDB_chain'] in type_SEQRES and type_SEQRES[form['PDB_chain']] == "AA":
        if vars['pdb_object'].get_type() == "AA":

            exit_on_error('user_error', "The selected chain: " + form['PDB_chain'] + " contains amino acid, and you have selected nucleic acid")

    # if modified residues exists, print them to the screen

    MODIFIED_COUNT = vars['pdb_object'].get_MODIFIED_COUNT()
    if MODIFIED_COUNT > 0:
        
        if form['DNA_AA'] == "AA":

            if len(vars['SEQRES_seq']) > 0 and MODIFIED_COUNT / len(vars['SEQRES_seq']) > GENERAL_CONSTANTS.MAXIMUM_MODIFIED_PERCENT:

                LOG.write("MODIFIED_COUNT %d\nSEQRES_seq %s\n" %(MODIFIED_COUNT, vars['SEQRES_seq']))
                exit_on_error('user_error', "Too many modified residues were found in SEQRES field; %0.3f%% of the residues are modified, the maximum is %0.3f%%. Please read the <a href=\"%s\">overview</a> for more information." %(MODIFIED_COUNT / len(vars['SEQRES_seq']) ,GENERAL_CONSTANTS.MAXIMUM_MODIFIED_PERCENT ,GENERAL_CONSTANTS.OVERVIEW_PAGE))

            LOG.write("analyse_seqres_atom : modified residues found\n")
            #print_message_to_output("<font color=\"red\">Please note:</font> Before the analysis took place, modified residues read from SEQRES field in the <a href=\"<?=$orig_path?>/" + vars['pdb_file_name'] + "\" style=\"color: #400080; text-decoration:underline;\">PDB</a> were converted back to the original residues:<br>" + vars['pdb_object'].get_MODIFIED_LIST() + ".<br>For a description of modified residues in PDB files <a href = \"" + GENERAL_CONSTANTS.PDB_FILE_MODRES_GUIDE + "\" style=\"color: #400080; text-decoration:underline;\" TARGET=Modified_Res>click here</a>.")

        else:

            LOG.write("analyse_seqres_atom : modified residues found\n")
            #print_message_to_output("<font color=\"red\">Please note:</font> Before the analysis took place, modified nucleotides read from SEQRES field in the <a href=\"<?=$orig_path?>/" + vars['pdb_file_name'] + "\" style=\"color: #400080; text-decoration:underline;\">PDB</a> were converted back to the original nucleotides:<br>" + vars['pdb_object'].get_MODIFIED_LIST() + ".<br>For a description of modified nucleotides in PDB files <a href = \"" + GENERAL_CONSTANTS.PDB_FILE_MODRES_GUIDE + "\" style=\"color: #400080; text-decoration:underline;\" TARGET=Modified_Res>click here</a>.")




def get_seqres_atom_seq(PDB_Obj, query_chain, pdb_file_name, model = False):

    # extract the sequences from the pdbParser

    seqres = ""
    atom = PDB_Obj.get_ATOM()
    atom_without_X = PDB_Obj.get_ATOM_withoutX()
    if query_chain in atom_without_X:
	
	    atom_without_X = atom_without_X[query_chain]

    all_seqres = PDB_Obj.get_SEQRES()
    for chain in all_seqres:

        if chain == query_chain:

            seqres = all_seqres[chain]
            break

    """
    all_atoms = PDB_Obj.get_ATOM()
    for chain in all_atoms:

        if chain == query_chain:

            atom = all_atoms[chain]
            atom_without_X = PDB_Obj.get_ATOM_withoutX(chain)
            break
    """

    if seqres == "" and atom == "":

        if model:

            return('user_error', "The protein sequence for chain '%s' was not found in SEQRES nor ATOM fields in the <a href=\"%s\">PDB file</a>." %(query_chain, pdb_file_name), "1")

        else:

            exit_on_error('user_error', "The protein sequence for chain '%s' was not found in SEQRES nor ATOM fields of the PDB file %s." %(query_chain, pdb_file_name))

    # output a message in case there is no seqres relevant sequence, but there are other chains.

    if seqres == "" and len(all_seqres) > 0:

        all_chains = ""

        for chain in all_seqres:

            all_chains += chain + ", "

        all_chains = all_chains[:-2]

        if model:

            return('user_error', "Chain \"%s\" does not exist in the SEQRES field of the <a href = \"%s\">PDB file</a>." %(chain, pdb_file_name), "1")

        else:

            exit_on_error('user_error', "Chain \"%s\" does not exist in the SEQRES field of the <a href = \"%s\">PDB file</a>.<br>Please check your file or run ConSurf again with a different chain." %(chain, pdb_file_name))

    return [seqres, atom, atom_without_X]





def determine_mode():

    mode = ""

    # determine the running mode and the pdb file mode
    #_mode_pdb_no_msa : (pdb_ID OR pdb_FILE) AND _form_pdb_chain AND build MSA parameters
    #_mode_pdb_msa : (pdb_ID OR pdb_FILE) AND _form_pdb_chain AND msa_FILE AND msa_SEQ_name
    #_mode_pdb_msa_tree: (pdb_ID OR pdb_FILE) AND _form_pdb_chain AND msa_FILE AND msa_SEQ_name AND tree_FILE
    #_mode_msa : msa_FILE AND msa_SEQ_name
    #_mode_msa_tree : msa_FILE AND msa_SEQ_name AND tree_FILE
    #_mode_no_pdb_no_msa : protein_SEQ or uploaded_SEQ AND build MSA parameters

    if form['pdb_FILE'] is not None or form['pdb_ID'] is not None:

        if vars['user_msa_file_name'] is not None:

            if form['tree_name'] is not None:

                mode = "_mode_pdb_msa_tree"

            else:

                mode = "_mode_pdb_msa"

        else:

            mode = "_mode_pdb_no_msa"

    else:

        if vars['user_msa_file_name'] is not None:

            if form['tree_name'] is not None:

                mode = "_mode_msa_tree"

            else:

                mode = "_mode_msa"

        else:

            mode = "_mode_no_pdb_no_msa"

    LOG.write("determine_mode : Mode is: " + mode + "\n")
    return mode





def prepare_pdb():

    LOG.write("form['pdb_ID']:%s*\n" %form['pdb_ID'])
    LOG.write("form['pdb_FILE']:%s*\n" %form['pdb_FILE'])
    RCSB_WGET = GENERAL_CONSTANTS.RCSB # web site to get pdb from
    local_PDB_path = GENERAL_CONSTANTS.PDB_DIVIDED # folder to get pdb from

    if form['pdb_ID'] is not None:

        pdb_ID = form['pdb_ID'].lower()
        LOG.write("PDB CASE is: ID, %s\n" %pdb_ID)
        PdbFileName = "pdb" + pdb_ID + ".ent.gz"  # pdb file name to check
        PdbFilePath = local_PDB_path + pdb_ID[1:3] + "/" + PdbFileName # pdb file absolute path to check
        LOG.write("PdbFilePath:%s\n" %PdbFilePath)

        if not os.path.exists(PdbFilePath):

            # if the file doesn't exist on local database, we try to get it from rcsb
            LOG.write("PDB id: %s was not found on our database. Trying to get it using command:\n\t%s%s\n" %(pdb_ID, RCSB_WGET, PdbFileName))

            try:

                response = request.urlopen(RCSB_WGET + PdbFileName, timeout=10).read()

            except:

                exit_on_error('user_error', "prepare_pdb : <b>The PDB file with the ID " + pdb_ID + " is not found on our database, or not available right now.</b><br />You may try to download it from <a href=\"RCSB\" target=\"RCSB\">RCSB</a> website, and upload it to ConSurf as FILE input.\n")

            try:

                GET_PDB = open(PdbFileName, 'wb')

            except:

                exit_on_error('system_error', "prepare_pdb : can't open the file " + PdbFileName + " for writing.")

            GET_PDB.write(response)
            GET_PDB.close()

        else:

            LOG.write("prepare_pdb : shutil.copyfile(%s, %s)\n" %(PdbFilePath, PdbFileName))
            shutil.copyfile(PdbFilePath, PdbFileName)

        # unzip the file

        try:

            PDB_UNZIP = open(vars['pdb_file_name'], 'wb')

        except:

            exit_on_error('system_error', "prepare_pdb : can't open the file " + vars['pdb_file_name'] + " for writing.")

        try:

            PDB_ZIP = gzip.open(PdbFileName, 'rb')

        except:

            exit_on_error('system_error', "prepare_pdb : can't open the file " + PdbFileName + " for reading.")

        PDB_UNZIP.write(PDB_ZIP.read())

        PDB_UNZIP.close()
        PDB_ZIP.close()

    elif(form['pdb_FILE'] is not None):

        LOG.write("PDB CASE is: FILE, %s\n" %form['pdb_FILE'])
        #upload_file(form['pdb_FILE'], vars['pdb_file_name'])

    else:

        exit_on_error('user_error',"There is no pdb input!!")






vars = {} # Hash with variables information
form = {} # Hash with form information
get_form()
    
"""
fields = cgi.FieldStorage()
for field in fields:

    if fields.getvalue(field): # field is not empty

        form[field] = fields.getvalue(field)

    else: # field is empty

        form[field] =  None
"""
if form['Homolog_search_algorithm'] == "BLAST":

    vars['BLAST_out_file'] = "protein_query.xml"

    if form['DNA_AA'] == "Nuc":

        vars['blast_algorithm'] = "BLAST"

    else:

        vars['blast_algorithm'] = "PSI-BLAST"

elif form['Homolog_search_algorithm'] == "HMMER":

    vars['blast_algorithm'] = "HMMER"
    vars['BLAST_out_file'] = "sequences_found_hmmer.txt"

else:

    vars['blast_algorithm'] = "CSI-BLAST"
    vars['BLAST_out_file'] = "protein_query.xml"

if form['DNA_AA'] == "AA": # proteins

    vars['protein_or_nucleotide'] = "proteins"
    vars['Msa_percentageFILE'] = "msa_aa_variety_percentage.csv"

else: # nucleotides

    vars['protein_or_nucleotide'] = "nucleotides"
    vars['Msa_percentageFILE'] = "msa_nucleic_acids_variety_percentage.csv"
"""
if form['Run_Number'] is None:

    vars['run_number'] = str(int(time.time()))  # creates a new number for the query
    vars['working_dir'] = GENERAL_CONSTANTS.SERVERS_RESULTS_DIR + vars['run_number'] + "/"

    # In case a directory with this name already exists, sleeps for 1 second and change the directory name
    while os.path.isdir(vars['working_dir']):

        time.sleep(1)
        vars['run_number'] = str(int(time.time()))
        vars['working_dir'] = GENERAL_CONSTANTS.SERVERS_RESULTS_DIR + vars['run_number'] + "/"

    os.mkdir(vars['working_dir'])


else:
"""

#vars['working_dir'] = GENERAL_CONSTANTS.SERVERS_RESULTS_DIR + form['Run_Number'] + "/"
#vars['server_html_dir'] = GENERAL_CONSTANTS.CONSURF_HTML_DIR
vars['run_log'] = vars['working_dir'] + "log.txt"
#vars['Confidence_link'] = GENERAL_CONSTANTS.CONSURF_URL + "overview.php#CONFIDENCE"
#vars['output_page'] = "output.php"
#vars['user_msa_file_name'] = "msa_file.aln"
#vars['tree_viewer_file'] = "treeView.html"
#vars['protein_seq'] = "protein_seq.fas" # a fasta file with the protein sequence from PDB or from protein seq input
#vars['All_Outputs_Zip'] = "Consurf_Outputs_" + form['Run_Number'] + ".zip"
"""
vars['send_email_cmd'] = "/bioseq/consurf/perl/share/sendEmail.pl"
vars['userName'] = GENERAL_CONSTANTS.ADMIN_USER_NAME
vars['userPass'] = GENERAL_CONSTANTS.ADMIN_PASSWORD
vars['smtp_server'] = GENERAL_CONSTANTS.SMTP_SERVER
vars['Server_Results_Path'] = "results/" + form['Run_Number'] + "/"
"""
if form['pdb_FILE'] is not None: # User PDB

    match = re.search(r'([^\/]+)$', form['pdb_FILE'])
    if match:

        vars['Used_PDB_Name'] = match.group(1)
        # Ofer bug fix: replace spaces with underscores in the used pdb name:
        vars['Used_PDB_Name'] = re.sub(r" ", r"_", vars['Used_PDB_Name'])
        vars['Used_PDB_Name'] = re.sub(r".pdb", r"", vars['Used_PDB_Name'])

elif form['pdb_ID'] is not None: # Given PDB_ID

    vars['Used_PDB_Name'] = form['pdb_ID']
    if form['PDB_chain'] != "NONE":

        vars['Used_PDB_Name'] += "_" + form['PDB_chain']

elif form['DNA_AA'] == "AA" and form['modeller_key'] is not None: # ConSeq Mode 

    vars['Used_PDB_Name'] = "hhpred_model"

else: # ConSeq Mode, No Model

    vars['Used_PDB_Name'] = "no_model"

vars['compare_two_fastas_count'] = 0
vars['FINAL_sequences'] = "query_final_homolougs.fasta"
vars['FINAL_sequences_html'] = "query_final_homolougs.html"
vars['submission_time'] = str(datetime.now())
vars['date'] = date.today().strftime("%d/%m/%Y")
vars['time_table'] = []
vars['current_time'] = time.time()
vars['r4s_log'] = "r4s.log"
vars['r4s_out'] = "r4s.res"
vars['r4s_slow_log'] = "r4s_slow.log"
vars['atom_positionFILE'] = "atom_pos.txt"
vars['gradesPE'] = "consurf_grades.txt"
vars['insufficient_data'] = "no"
vars['protein_MSA_for_usr'] = ""
#vars['zip_list'] = []

"""
vars['scf_for_chimera'] = "ConSurf_%s_RUN_%s.scf" %(vars['Used_PDB_Name'], vars['run_number'])
vars['header_for_chimera'] = "ConSurf_%s_RUN_%s.hdr" %(vars['Used_PDB_Name'], vars['run_number'])
vars['chimerax_file'] = "ConSurf_%s_RUN_%s.chimerax" %(vars['Used_PDB_Name'], vars['run_number'])
vars['isd_scf_for_chimera'] = "ConSurf_%s_RUN_%s_isd.scf" %(vars['Used_PDB_Name'], vars['run_number'])
vars['isd_header_for_chimera'] = "ConSurf_%s_RUN_%s_isd.hdr" %(vars['Used_PDB_Name'], vars['run_number'])
vars['isd_chimerax_file'] = "ConSurf_%s_RUN_%s_isd.chimerax" %(vars['Used_PDB_Name'], vars['run_number'])
"""


vars['WASABI_XML'] = "MSA_and_Tree.xml"
vars['msa_fasta'] = "msa_fasta.aln" # if the file is not in fasta format, we create a fasa copy of it
vars['msa_clustal'] = "msa_clustal.aln" # if the file is not in clustal format, we create a clustal copy of it

HTML_TXT_FOR_PISA_MODEL = ""

residue_freq = {} # for each position in the MSA, detail the residues
position_totalAA = {} # for each position in the MSA, details the total number of residues
gradesPE_Output = [] # an array to hold all the information that should be printed to gradesPE
# in each array's cell there is a hash for each line from r4s.res.
# POS: position of that aa in the sequence ; SEQ : aa in one letter ;
# GRADE : the given grade from r4s output ; COLOR : grade according to consurf's scale

os.chdir(vars['working_dir'])

LOG = open(vars['run_log'], 'w')
vars['running_mode'] = determine_mode()
"""
create_output_php()
vars['run_url'] = "%sbarak/progress.php?number=%s" %(GENERAL_CONSTANTS.CONSURF_URL, form['Run_Number'])
print("<meta http-equiv=\"refresh\" content=\"0; URL=%s\" />" %vars['run_url'])

sys.stdin.close()
sys.stdout.close()
sys.stderr.close()
os.close(0)
os.close(1)
#os.close(2)
sys.stderr = open("STD_ERROR", 'w')
Update_ConSurf_Users_Log()

vars['zip_list'].append(vars['tree_file'])
vars['zip_list'].append(vars['gradesPE'])
vars['zip_list'].append(vars['Msa_percentageFILE'])
vars['zip_list'].append(vars['Colored_Seq_PDF'])
vars['zip_list'].append(vars['Colored_Seq_CBS_PDF'])
vars['zip_list'].append(vars['msa_fasta'])
"""
## mode : include pdb

# create a pdbParser, to get various info from the pdb file
if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_pdb_msa_tree":

    """
    if form['PDB_chain'] != "NONE": 

        vars['Used_PDB_Name'] += "_" + form['PDB_chain']
    """

    #prepare_pdb()
    vars['pdb_object'] = pdbParser.pdbParser()
    vars['pdb_object'].read(vars['pdb_file_name'], form['PDB_chain'], form['DNA_AA'], vars['atom_positionFILE'])

    #vars['D_Nuc'] = vars['pdb_object'].IS_D_NUC(form['PDB_chain'])

    # check if there is no seqres
    [vars['SEQRES_seq'], vars['ATOM_seq'], vars['ATOM_without_X_seq']] = get_seqres_atom_seq(vars['pdb_object'], form['PDB_chain'], vars['pdb_file_name'])
    analyse_seqres_atom()

    try:

        FAS = open(vars['protein_seq'], 'w')

    except:

        exit_on_error('sys_error',"cannot open the file " + vars['protein_seq'] + " for writing!")

    if vars['SEQRES_seq'] == "":

        vars['query_string'] = "Input_pdb_ATOM_" + form['PDB_chain']
        vars['protein_seq_string'] = vars['ATOM_without_X_seq']
        #FAS.write(">PDB_ATOM\n" + vars['ATOM_without_X_seq'])
        FAS.write(">" + vars['query_string'] + "\n" + vars['ATOM_without_X_seq'])

    else:

        vars['query_string'] = "Input_pdb_SEQRES_" + form['PDB_chain']
        vars['protein_seq_string'] = vars['SEQRES_seq']
        #FAS.write(">PDB_SEQRES\n" + vars['SEQRES_seq'])
        FAS.write(">" + vars['query_string'] + "\n" + vars['SEQRES_seq'])

    FAS.close()

    #Update_Progress()

## mode : only protein sequence

# if there is only protein sequence: we upload it.
elif vars['running_mode'] == "_mode_no_pdb_no_msa":

    upload_protein_sequence()

## mode : no msa - with PDB or without PDB

if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_no_pdb_no_msa":

    # if there is pdb : we compare the atom and seqres
    if vars['running_mode'] ==  "_mode_pdb_no_msa" and ('SEQRES_seq' in vars and len(vars['SEQRES_seq']) > 0):

        # align seqres and pdb sequences
        compare_atom_seqres_or_msa("SEQRES")

    """
    vars['BLAST_out_file'] = "protein_query.xml" # file to hold blast output
    vars['HMMER_out_file'] = "protein_query.hmmer"
    vars['BLAST_last_round'] = "last_round.blast" # file to hold blast output, last round
    """

    vars['max_homologues_to_display'] = GENERAL_CONSTANTS.BLAST_MAX_HOMOLOGUES_TO_DISPLAY

    if form['Homolog_search_algorithm'] != "HMMER":

        if form['proteins_DB'] == "SWISS-PROT":

            vars['protein_db'] = GENERAL_CONSTANTS.SWISSPROT_DB

        elif form['proteins_DB'] == "UNIREF90":

            vars['protein_db'] = GENERAL_CONSTANTS.UNIREF90_DB

        elif form['proteins_DB'] == "CLEAN_UNIPROT":

            vars['protein_db'] = GENERAL_CONSTANTS.CLEAN_UNIPROT_DB

        elif form['proteins_DB'] == "NR_PROT_DB":

            vars['protein_db'] = GENERAL_CONSTANTS.NR_PROT_DB

        else:

            vars['protein_db'] = GENERAL_CONSTANTS.UNIPROT_DB

    else: # form['Homolog_search_algorithm'] == "HMMER"

        if form['proteins_DB'] == "SWISS-PROT":

            vars['protein_db'] = GENERAL_CONSTANTS.SWISSPROT_DB_FASTA

        elif form['proteins_DB'] == "UNIREF90":

            vars['protein_db'] = GENERAL_CONSTANTS.UNIREF90_DB_FASTA

        elif form['proteins_DB'] == "CLEAN_UNIPROT":

            vars['protein_db'] = GENERAL_CONSTANTS.CLEAN_UNIPROT_DB_FASTA

        elif form['proteins_DB'] == "NR_PROT_DB":

            vars['protein_db'] = GENERAL_CONSTANTS.NR_PROT_DB_FASTA

        else:

            vars['protein_db'] = GENERAL_CONSTANTS.UNIPROT_DB_FASTA


    blast_hash = {}

    # Manual Selection of Sequences Mode
    # create the seqres fasta file, run blast
    if form['user_select_seq'] == "yes":

        form['MAX_NUM_HOMOL'] = "2000"
        last_round_number = 0
        found_converged = 0

        # create the seqres fasta file, run blast and create form so user can select sequences
        vars['BLAST_out_file_html'] = "hmmer_or_blast.htm"
        vars['Form_Selected_Seq_ID_File'] = vars['BLAST_out_file'] + ".User_Selection"
        vars['HITS_rejected_file'] = "query_rejected_homolougs.fas"

        run_search(vars['BLAST_out_file'], "PlainText")

        """
        # if AA use HMMER or BLAST, if NUC use BLAST

        if form['Homolog_search_algorithm'] == "HMMER":

            # HMMER
            run_search(vars['BLAST_out_file'], "PlainText")
            #last_round_number = extract_round_from_hmmer()
            #alignments = Create_HTML_file_for_hmmer(vars['BLAST_out_file_html'], "hmmer_seqs.htm", vars['BLAST_out_file'], vars['query_string'], form['Run_Number'])

        else:

            # AA BLAST
            run_search(vars['BLAST_out_file_html'], "HTML")
            #last_round_number = get_blast_round(vars['BLAST_out_file_html'])
        """

        alignments = Create_HTML_file(vars['BLAST_out_file_html'], "hmmer_seqs.htm", vars['BLAST_out_file'], vars['query_string'], form['Run_Number'])
        try:

            choose_flag = open("CHOOSE", 'w')

        except:

            exit_on_error('sys_error', "can't open the file CHOOSE for writing.")

        """
        LOG.write("last_round_number: " + str(last_round_number) + "\n")
        vars['BLAST_Form'] = vars['BLAST_out_file'] + "_form.tpl"
        create_blast_form(last_round_number)
        OUTPUT.close()
        if vars['output_page'] != vars['output_page_without_form']:

            shutil.copyfile(vars['output_page'], vars['output_page_without_form'])

        add_form_to_output()
        """
        #Update_Progress()
        send_mail_to_select_sequences()
        while os.path.exists("CHOOSE"):

            time.sleep(1)

        final_homologoues_choosen_manually(alignments)
        """
        while not os.path.exists(vars['Form_Selected_Seq_ID_File']): # User did not choose his sequences

            time.sleep(1)

        remove_form()
        vars['output_page'] = vars['output_page_without_form'] # return to the output without the form
        try:

            OUTPUT = open(vars['output_page'], 'a')

        except:

            exit_on_error('sys_error',"create_output_php : could not open the file " + vars['output_page'] + " for writing")

        if form['Homolog_search_algorithm'] != "HMMER":

            remove_html_tags(vars['BLAST_out_file_html'], vars['BLAST_out_file'])
            LOG.write("parseFiles.print_blast_according_to_round(\"%s\", %s, \"%s\");\n" %(vars['BLAST_out_file'], last_round_number, vars['BLAST_last_round']))
            ans = parseFiles.print_blast_according_to_round(vars['BLAST_out_file'], last_round_number, vars['BLAST_last_round'])
            LOG.write("ans: %s, %s\n" %(ans[0], ans[1]))
            if ans[0] == "err":

                exit_on_error('sys_error',"extract_round_from_blast : " + ans[1])

            if not os.path.exists(vars['BLAST_last_round']) and os.path.getsize(vars['BLAST_last_round']) == 0:

                LOG.write("the file " + vars['BLAST_last_round'] + " exist and not empty!\n")
                shutil.copyfile(vars['BLAST_last_round'], vars['BLAST_out_file'])

            else:

                LOG.write("extract_round_from_blast : the file " + vars['BLAST_last_round'] + " was not created. The hits will be collected from the original blast file")

            vars['unique_seqs'] = choose_final_homologoues_according_user(blast_hash)
            print_message_to_output(vars['unique_seqs'] + " sequences were selected for the analysis. You can observe the selected sequences <A HREF=\"" + vars['FINAL_sequences_html'] + "\" TARGET=Seq_for_Analysis>here</A>.")

            # exit if there are not enough seqs
            if vars['unique_seqs'] < vars['min_num_of_hits']:

                msg = "Only " + vars['unique_seqs']
                msg += "  <a href=\"%s\">seqences</a> were selected. The minimal number of sequences required for the calculation is %d. You may try to:<ol>" %(vars['FINAL_sequences'], vars['min_num_of_hits'])
                msg += "<li>Re-run the server and select additional homologous sequences.<li>Re-run the server with a multiple sequence alignment file of your own.</li>"
                msg += "</OL>\n"
                exit_on_error('user_error', msg)

            elif vars['unique_seqs'] < vars['low_num_of_hits']:

                msg = "<font color='red'><b>Warning:</font></b> only %s <a href=\"%s\">seqences</a> were selected.   The calculation is performed on the <a href=\"%s\">selected sequences</a>, but it is recommended to run the server with a multiple sequence alignment file containing at least %d sequences." %(vars['unique_seqs'], vars['FINAL_sequences'], vars['FINAL_sequences_html'], vars['low_num_of_hits'])
                print_message_to_output(msg)

            # end if user select
        """

    # Automatic Mode - Without manual homologues selection
    # create the seqres fasta file, run blast
    else:

        run_search(vars['BLAST_out_file'], "PlainText")
        #Update_Progress()

        """
        if form['Homolog_search_algorithm'] == "HMMER":

            extract_round_from_hmmer()

        else:

            extract_round_from_blast()
        """
        # choosing homologs, create fasta file for all legal homologs
        cd_hit_hash = {}
        vars['hit_redundancy'] = float(form['MAX_REDUNDANCY']) # Now taken as argument from user #OLD: #CONSURF_CONSTANTS.FRAGMENT_REDUNDANCY_RATE
        vars['hit_overlap'] = GENERAL_CONSTANTS.FRAGMENT_OVERLAP
        vars['hit_min_length'] = GENERAL_CONSTANTS.FRAGMENT_MINIMUM_LENGTH
        vars['min_num_of_hits'] = GENERAL_CONSTANTS.MINIMUM_FRAGMENTS_FOR_MSA
        vars['low_num_of_hits'] = GENERAL_CONSTANTS.LOW_NUM_FRAGMENTS_FOR_MSA
        vars['HITS_fasta_file'] = "query_homolougs.txt"
        vars['HITS_rejected_file'] = "query_rejected_homolougs.txt"

        choose_homologoues_from_search_with_lower_identity_cutoff(blast_hash)

        vars['cd_hit_out_file'] = "query_cdhit.out"
        vars['unique_seqs'] = cluster_homologoues(cd_hit_hash)
        LOG.write("num_of_unique_seq: %d\n" %vars['unique_seqs'])
        add_sequences_removed_by_cd_hit_to_rejected_report()
        choose_final_homologoues(cd_hit_hash, blast_hash)
    if form['DNA_AA'] == "Nuc":

        # convert rna to dna

        LOG.write("convert_rna_to_dna(%s, %s)\n" %(vars['FINAL_sequences'], vars['FINAL_sequences'] + ".dna"))
        ans = convert_rna_to_dna(vars['FINAL_sequences'], vars['FINAL_sequences'] + ".dna")
        if ans[0] == "OK":

            vars['FINAL_sequences'] += ".dna"
            LOG.write("Seqs with u or U: " + str(ans[1]))
            """
            for seq in ans[1]:

                print_message_to_output ("<b><font color='red'>Warnning: </font></b>The seqeunce '" + seq + "' contains a 'U' replaced by 'T'")
            """
        else:

            exit_on_error('sys_error', ans)

    """

    else:

        # calculate max and media e-value of used sequences
        LOG.write("parseFiles.fasta_max_median_E_value %s\n" %(vars['FINAL_sequences']))
        results = parseFiles.fasta_max_median_E_value(vars['FINAL_sequences'])

        if results[0] != -1:

            vars['MaxFinalSequenceEvalue'] = results[0]
            vars['MedianFinalSequenceEvalue'] = results[1]

        results[0] = -1
        results[1] = -1

        LOG.write("max Eval = " + str(results[0]) + ", median Eval = " + str(results[1]) + "\n")
    """

    LOG.write("make_sequences_file_HTML(%s, %s)\n" %(vars['FINAL_sequences'], vars['FINAL_sequences_html']))
    #make_sequences_file_HTML(vars['FINAL_sequences'], vars['FINAL_sequences_html'])
    #Update_Progress()

    # we save to copies of the msa, one in fasta format and another in clustal format.
    #vars['msa_fasta'] = "msa_fasta.aln"
    #vars['msa_clustal'] = "msa_clustal.aln"
    create_MSA()
    #Update_Progress()
    vars['msa_SEQNAME'] = vars['query_string']
    #parameters_html(False, False)

## mode :  include msa

elif vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_msa" or vars['running_mode'] == "_mode_pdb_msa_tree" or vars['running_mode'] == "_mode_msa_tree":

    # check that there are at least 5 sequecnes in the MSA
    # extract the query sequence from the MSA
    # change the MSA sequences names to numbers
    # change the MSA format to clustalw
    # align the seqres/atom sequence with that of the query
    MSA_sequences = {} # a hash to hold all the MSA sequences, : key - sequence id, value - sequence
    #vars['msa_fasta'] = "msa_fasta.aln" # if the file is not in fasta format, we create a fasa copy of it
    #vars['msa_clustal'] = "msa_clustal.aln" # if the file is not in clustal format, we create a clustal copy of it
    vars['msa_SEQNAME'] = form['msa_SEQNAME']
    vars['query_string'] = form['msa_SEQNAME']
    vars['MSA_query_seq'] = get_info_from_msa(MSA_sequences)

    if vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_pdb_msa_tree":

        compare_atom_seqres_or_msa("MSA")

    else:

        # there is no input seq, use msa seq instead
        vars['protein_seq_string'] = vars['MSA_query_seq']
        try:

            QUERY_FROM_MSA = open(vars['protein_seq'], 'w')

        except:

            exit_on_error('sys_error', "could not open " + vars['protein_seq'] )

        QUERY_FROM_MSA.write(">" + form['msa_SEQNAME'] + "\n")
        QUERY_FROM_MSA.write(vars['MSA_query_seq'] + "\n")
        QUERY_FROM_MSA.close()

    ## mode : include tree

    if vars['running_mode'] == "_mode_pdb_msa_tree" or vars['running_mode'] == "_mode_msa_tree":

        #parameters_html(True, True)
        #upload_tree()
        check_validity_tree_file()
        tree_nodes = [] # an array to hold all the nodes in the tree (as keys)
        extract_nodes_from_tree(tree_nodes)
        check_msa_tree_match(MSA_sequences, tree_nodes)
        TREE_parser.removeBPvalues(vars['tree_file'], vars['tree_file'] + "USER_ORIG")

    """
    else:

        parameters_html(True, False)
        
    """
if form['SUB_MATRIX'] == "BEST":

    vars['best_fit'] = True
    best_matrix = find_best_substitution_model()
    form['SUB_MATRIX'] = best_matrix

    if best_matrix == "JC_Nuc":

        best_matrix = "JC"

    #print_message_to_output ("The best evolutionary model was selected to be: " + best_matrix + ". (<A href=\"<?=$orig_path?>/model_selection.txt\" style=\"color: #400080; text-decoration:underline;\">details</A>).<br>")
    #print_message_to_output ("The best evolutionary model was selected to be: " + best_matrix + ". See details <A href=\"<?=$orig_path?>/model_selection.txt\" style=\"color: #400080; text-decoration:underline;\">here</A><br>")
    #Update_Progress()

else:

    vars['best_fit'] = False

run_rate4site()
#more_parameters_html()
#Update_Progress()
assign_colors_according_to_r4s_layers(gradesPE_Output, vars['r4s_out'])
read_residue_variety(residue_freq, position_totalAA)

if form['DNA_AA'] == "AA":

    print_residue_precentage()

else:

    print_nucleotide_precentage()


## mode : include pdb

# in order to create 3D outputs, we need to compare the ATOM to the sequence from rate4site
if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_pdb_msa_tree":

    #create_atom_position_file(form['PDB_chain'], vars['pdb_file_name'], vars['atom_positionFILE'], form['DNA_AA'], 0) # this file will be used later to create the output which aligns rate4site sequence with the ATOM records
    r4s2pdb = {} # key: poistion in SEQRES/MSA, value: residue name with position in atom (i.e: ALA22:A)

    if vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_pdb_msa_tree" or vars['SEQRES_seq'] != "":

        [LENGTH_OF_SEQRES, LENGTH_OF_ATOM] = match_pdb_to_seq(r4s2pdb, form['PDB_chain'], vars['seqres_or_msa_seq_with_gaps'], vars['ATOM_seq_with_gaps'], vars['atom_positionFILE'], form['DNA_AA'])

    else: # no seqres and no msa

        LENGTH_OF_ATOM = fill_r4s2pdb(r4s2pdb) # fill length_of_atom
        LENGTH_OF_SEQRES = 0

    scoresFile = "consurf.scores"
    no_isd_residue_color = [[],[],[],[],[],[],[],[],[],[],[]]
    isd_residue_color = [[],[],[],[],[],[],[],[],[],[],[]]
    # these arrays will hold for each grade, the residues which corresponds to it.
    # there are 2 arrays: in the @isd_residue_color, a grade with insufficient data, *, will classify to grade 10
    # in the @no_isd_residue_color, the grade will be given regardless of the * mark
    # PLEASE NOTE : the [0] position in those arrays is empty, because each position corresponds a color on a 1-10 scale

    # The following routine, apart from creating the file "consurf.grades" also collects the information in order to create
    # the RasMol scripts and a variable which will be used in the "pipe" file, that holds a string with the grades.
    # In the pipe file this var is called: seq3d_grades_isd and seq3d_grades
    #[SEQ3D_GRADES_ISD, SEQ3D_GRADES] = create_gradesPE_ConSurf(vars['gradesPE'], r4s2pdb, gradesPE_Output, residue_freq, no_isd_residue_color, isd_residue_color, form['DNA_AA'], vars['D_Nuc'])

    if form['DNA_AA'] == "AA":

        #LOG.write("Get_NACSES_buried_Exposed(%s, %s, %s)\n" %(vars['pdb_file_name'], form['PDB_chain'], vars['atom_positionFILE']))
        [NACSES_RESULTS, Pos_Solv_Acc_Pred_NACSES_hashRef] = Get_NACSES_buried_Exposed(vars['pdb_file_name'], form['PDB_chain'], vars['atom_positionFILE'])
        #LOG.write("Get_NACSES_buried_Exposed:%s\n" %NACSES_RESULTS)
        #LOG.write("Get_NACSES_buried_Exposed:%s\n" %str(Pos_Solv_Acc_Pred_NACSES_hashRef))
        #B_E_column = "9" # the number of the column showing buried or exposed in grades file.
        #prediction_method = "NACSES-algorithm"

    else:

        Pos_Solv_Acc_Pred_NACSES_hashRef = ""
        #B_E_column = "none" # the number of the column showing buried or exposed in grades file.
        #prediction_method = "no-prediction"

    [SEQ3D_GRADES_ISD, SEQ3D_GRADES] = create_gradesPE_ConSurf(vars['gradesPE'], r4s2pdb, gradesPE_Output, residue_freq, no_isd_residue_color, isd_residue_color, form['DNA_AA'], Pos_Solv_Acc_Pred_NACSES_hashRef)
    #cp_rasmol_gradesPE_and_pipe.create_score_file(gradesPE_Output, scoresFile, r4s2pdb)
    """
    # FIND IDENTICAL CHAINS
    identical_chains = find_identical_chains_in_PDB_file(vars['pdb_object'], form['PDB_chain'], 1)

    # This will create the pipe file for FGiJ
    pipeFile = vars['Used_PDB_Name'] + "_consurf_firstglance.pdb"
    pipeFile_CBS = vars['Used_PDB_Name'] + "_consurf_firstglance_CBS.pdb" # pipe for color blind friendly
    create_pipe_file(pipeFile, pipeFile_CBS, SEQ3D_GRADES, SEQ3D_GRADES_ISD, isd_residue_color, no_isd_residue_color, LENGTH_OF_SEQRES, LENGTH_OF_ATOM, vars['pdb_file_name'], form['PDB_chain'], (vars['Used_PDB_Name']).upper(), identical_chains)

    if form['pdb_ID'] is not None:

        mmtf = downloadMmtf()

    else:

        mmtf = None

    # json for embedded NGL
    relative_path = "../results/" + form['Run_Number'] + "/"
    Hash_Json = {'scoresFile' : relative_path + scoresFile, 'pdbFile' : relative_path + vars['pdb_file_name'], 'chainId' : form['PDB_chain']}
    if mmtf is not None:

        Hash_Json['mmtfFile'] = relative_path + mmtf

    saveJobInfoJson("job_info_small.json", Hash_Json)

    # json for NGL
    Hash_Json = {'scoresFile' : scoresFile, 'pdbFile' : vars['pdb_file_name'], 'chainId' : form['PDB_chain']}
    if mmtf is not None:

        Hash_Json['mmtfFile'] = mmtf

    saveJobInfoJson("job_info.json", Hash_Json)
    """
    vars['ATOMS_with_ConSurf_Scores'] = vars['Used_PDB_Name'] + "_ATOMS_section_With_ConSurf.pdb"
    vars['ATOMS_with_ConSurf_Scores_isd'] = vars['Used_PDB_Name'] + "_ATOMS_section_With_ConSurf_isd.pdb"
    vars['pdb_file_with_score_at_TempFactor'] = vars['Used_PDB_Name'] + "_With_Conservation_Scores.pdb"

    #vars['zip_list'].append(vars['ATOMS_with_ConSurf_Scores'])
    #vars['zip_list'].append(vars['ATOMS_with_ConSurf_Scores_isd'])
    #vars['zip_list'].append(vars['pdb_file_with_score_at_TempFactor'])

    input_chimera_pymol = replace_TmpFactor_Consurf_Scores(form['PDB_chain'], vars['pdb_file_name'], vars['gradesPE'], vars['ATOMS_with_ConSurf_Scores'], vars['ATOMS_with_ConSurf_Scores_isd']) # Create ATOMS section and replace the TempFactor Column with the ConSurf Grades (will create also isd file if relevant)
    #create_chimera(input_chimera_pymol, vars['Used_PDB_Name'] + "_")
    #create_pymol(input_chimera_pymol, vars['Used_PDB_Name'] + "_")
    replace_TmpFactor_Rate4Site_Scores(form['PDB_chain'], vars['pdb_file_name'], vars['gradesPE'], vars['pdb_file_with_score_at_TempFactor']) # Create PDB file and replace the TempFactor Column with the Rate4Site grades
    
    """
    create_pdf("5", B_E_column, prediction_method)

    if form['pdb_ID'] is not None and form['DNA_AA'] == "AA":
        # get complex from PISA
        vars['PISA_PDB'] = (form['pdb_ID']).upper() + "_PISA1.pdb"
        ans2 = Get_PISA_Complex(form['pdb_ID'], vars['PISA_PDB'])
        if ans2 != "OK":

            LOG.write("No PISA complex for %s\n" %form['pdb_ID'])

        else:

            LOG.write("PISA complex model was saved to %s\n" %(vars['PISA_PDB']))
            LOG.write("Project_ConSurf_On_Model(%s, %s, %s, %s, %s, %s, %s, %s)\n" %(vars['PISA_PDB'], form['PDB_chain'], vars['working_dir'], vars['protein_seq_string'], "r4s.res", vars['msa_clustal'], "True", form['DNA_AA']))
            Project_ConSurf_On_Model(vars['PISA_PDB'], form['PDB_chain'], vars['working_dir'], vars['protein_seq_string'], "r4s.res", vars['msa_clustal'], True, form['DNA_AA'])
	"""
    
## mode : ConSeq - NO PDB

if vars['running_mode'] == "_mode_msa" or vars['running_mode'] == "_mode_no_pdb_no_msa" or vars['running_mode'] == "_mode_msa_tree":

    Pos_Solv_Acc_Pred = {} # key: poistion in SEQRES/MSA, value: PACC Solv Acc Pred - relevant only for AA case

    if form['DNA_AA'] == "AA":

        vars['protein_MSA_hssp'] = vars['msa_fasta'] + ".hssp"
        vars['Solv_ACC_Pred'] = "Solv_Acc_Pred.PACC"
        vars['HHPred_Model'] = "HHPred_Model.ent"
        #Update_Progress()
        form['modeller_key'] = None
        if form['modeller_key'] is not None:

            # run HHPred
            LOG.write("Run_HHPred()\n")
            ans = Run_HHPred()
            if ans == "OK":

                #Update_Progress()
                Project_ConSurf_On_Model(vars['HHPred_Model'], "NONE", vars['working_dir'], vars['protein_seq_string'], "r4s.res", vars['msa_clustal'], False, form['DNA_AA'])

            # predict b/e by neural-network algorithm
            predict_solvent_accesibility(vars['msa_clustal'], vars['protein_MSA_hssp'], vars['Solv_ACC_Pred']) #runs the PACC algorithm to calculate burried/exposed
            read_Solv_Acc_Pred(Pos_Solv_Acc_Pred)
            create_gradesPE_ConSeq(Pos_Solv_Acc_Pred, "Amino")
            #B_E_column = "8" # the number of the column showing buried or exposed in grades file.
            #prediction_method = "neural-network-algorithm"

            """
            ConSeq_gradesPE_and_Outputs.ConSeq_HTML_Output(gradesPE_Output, Pos_Solv_Acc_Pred, vars['Colored_Seq_HTML'])
            Protein_Length = len(gradesPE_Output) + 1
            ConSeq_gradesPE_and_Outputs.ConSeq_PDF_Output(vars['gradesPE'], vars['Colored_Seq_PDF'], Protein_Length, "neural-network algorithm")
            """
        else: # v/e based on prediction

            #predict_solvent_accesibility(vars['msa_clustal'], vars['protein_MSA_hssp'], vars['Solv_ACC_Pred']) #runs the PACC algorithm to calculate burried/exposed
            #read_Solv_Acc_Pred(Pos_Solv_Acc_Pred)
            create_gradesPE_ConSeq(Pos_Solv_Acc_Pred, "Amino")
            #B_E_column = "8" # the number of the column showing buried or exposed in grades file.
            #prediction_method = "neural-network-algorithm"

            """
            ConSeq_gradesPE_and_Outputs.ConSeq_HTML_Output(gradesPE_Output, Pos_Solv_Acc_Pred, vars['Colored_Seq_HTML'])
            Protein_Length = len(gradesPE_Output) + 1
            ConSeq_gradesPE_and_Outputs.ConSeq_PDF_Output(vars['gradesPE'], vars['Colored_Seq_PDF'], Protein_Length, "neural-network algorithm")
            """
    else: # Nucleic Acid sequence

        #vars['D_Nuc'] = ""
        create_gradesPE_ConSeq("", "Nucleic")
        B_E_column = "none" # the number of the column showing buried or exposed in grades file.
        #prediction_method = "no-prediction"

        """
        ConSeq_gradesPE_and_Outputs.ConSeq_HTML_Output(gradesPE_Output, Pos_Solv_Acc_Pred, vars['Colored_Seq_HTML'], "neural-network algorithm", True) # The last argument ("yes") tell us that only the sequence should be printed
        Protein_Length = len(gradesPE_Output) + 1
        ConSeq_gradesPE_and_Outputs.ConSeq_NUC_PDF_Output(vars['gradesPE'], vars['Colored_Seq_PDF'], Protein_Length)
        # predict RNA SS
        [vars['RNA_ss_ps'], vars['RNA_dp_ps'], vars['RNA_ss_pdf']] = Predict_RNA_SS(vars['protein_seq'], vars['working_dir'])
        if os.path.exists(vars['working_dir'] + vars['RNA_ss_ps']):

            [vars['RNA_ss_colored_by_ConSurf_ps'], vars['RNA_ss_colored_by_ConSurf_pdf']] = Color_RNA_SS_by_ConSurf(vars['RNA_ss_ps'], vars['gradesPE'], vars['working_dir'])
            [vars['RNA_ss_colored_by_ConSurf_ps_CBS'], vars['RNA_ss_colored_by_ConSurf_pdf_CBS']] = Color_RNA_SS_by_ConSurf(vars['RNA_ss_ps'], vars['gradesPE'], vars['working_dir'], "CBS")
        """

    #create_pdf("4", B_E_column, prediction_method)


# create wasabi xml
#LOG.write("Create file for wasabi: wasabiUtils.createWasabiXml(%s, %s, %s)" %(vars['msa_fasta'], vars['tree_file'], vars['WASABI_XML']))
#wasabiUtils.createWasabiXml(vars['msa_fasta'], vars['tree_file'], vars['WASABI_XML'])

#Prepare_tree_View()
#Create_Colored_MSA()

#ConSeq_gradesPE_and_Outputs.consurf_HTML_Output(gradesPE_Output, vars['Colored_Seq_HTML'], False)
#ConSeq_gradesPE_and_Outputs.consurf_HTML_Output(gradesPE_Output, vars['Colored_Seq_HTML_CBS'], True)

#zip_all_outputs()
#Update_Progress()

## Arrange The HTML Output File

#send_finish_email_to_user()
LOG.write("\nDONE")
LOG.close()
"""
try:

    DONE = open("DONE", 'w')

except:

    exit_on_error('sys_error', "can't open the file DONE for writing.")

DONE.close()
"""

