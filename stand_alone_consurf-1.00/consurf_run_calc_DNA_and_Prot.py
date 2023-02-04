

import shutil
import json
import re
import sys
import os
import GENERAL_CONSTANTS
import CONSURF_CONSTANTS
import MSA_parser
import subprocess
import rate4site_routines
import cp_rasmol_gradesPE_and_pipe
import time
import datetime
import pdbParser
import ConSeq_gradesPE_and_Outputs
import wasabiUtils
from Bio import AlignIO

dir = "/"
consurf_perl_dir = ""
server_html_dir = GENERAL_CONSTANTS.CONSURF_HTML_DIR


stored_data_file = sys.argv[1]
stored_form_data = sys.argv[2]



# loading vars hash
try:

    VARS_LOAD = open(stored_data_file, 'r')

except:

    exit_on_error("sys_error", "store_data : could not open %s for storing." %stored_data_file)

vars = json.load(VARS_LOAD)
VARS_LOAD.close()

# loading form hash
try:

    FORM_LOAD = open(stored_form_data, 'r')

except:

    exit_on_error("sys_error", "store_data : could not open %s for storing." %stored_form_data)

form = json.load(FORM_LOAD)
FORM_LOAD.close()

os.chdir(vars['working_dir'])

gershon = ""
HTML_TXT_FOR_MODEL_SECTION = "" # will be filled if a model is created by HHPred [ConSeq mode]
HTML_TXT_FOR_MODEL_SECTION_FGIJ = "" # will be filled if a model is created by HHPred [ConSeq mode] - Only the FGIJ link
HTML_TXT_FOR_PISA_MODEL = "" # will be filed  if PISA complex is available
BURRIED_CUTOFF = 16
EXPOSED_CUTOFF = 16
Pos_Solv_Acc_Pred_NACSES_hashRef = "" # will hold ref to hash for solv acc by NACSES key: pos in PDB; value: b/e
# here we handle rerun of ConSurf with sub tree data called directly by: /bioseq/wasabi/scripts/runConsurfOnSubtree.php
if len(sys.argv) > 3:

    # update VARS and FORM variables
    vars['old_run_number'] = vars['run_number']
    vars['run_number'] = sys.argv[3] # the new run number for the sub tree

    vars['old_WorkingDir'] = vars['working_dir']
    vars['working_dir'] = CONSURF_CONSTANTS.CONSURF_RESULTS_DIR + dir + vars['run_number'] + dir

    vars['old_run_url'] = vars['run_url']
    vars['run_url'] = CONSURF_CONSTANTS.CONSURF_RESULTS_URL + dir + vars['run_number']

    vats['old_run_log'] = vars['run_log']
    vars['run_log'] = CONSURF_CONSTANTS.CONSURF_LOGS_DIR + dir + vars['run_number'] + ".log"

    vars['run_log_Q'] = CONSURF_CONSTANTS.CONSURF_LOGS_DIR + dir + vars['run_number'] + "Q.log"
    vars['run_url_old'] = GENERAL_CONSTANTS.CONSURF_URL + "results" + dir + vars['run_number'] # To remove from all places... leave it here till then...

    # copy query seq name and extract query seq
    Query_Seq = ""
    Query_SeqName = ""
    print("QUERY STRING: " + vars['query_string'] + "\n")
    if os.path.exists(vars['old_WorkingDir'] + dir + vars['protein_seq']):

        shutil.copyfile(vars['old_WorkingDir'] + dir + vars['protein_seq'], vars['working_dir'] + vars['protein_seq'])
        print("shutil.copyfile(%s, %s)\n" %(vars['old_WorkingDir'] + dir + vars['protein_seq'], vars['working_dir'] + vars['protein_seq']))
        (Query_Seq, Query_SeqName) = Extract_Fasta_Single_Seq(vars['working_dir'] + vars['protein_seq'])

    elif os.path.exists(vars['old_WorkingDir'] + dir + vars['query_seq_file_from_MSA']):

        shutil.copyfile(vars['old_WorkingDir'] + dir + vars['query_seq_file_from_MSA'], vars['working_dir'] + vars['query_seq_file_from_MSA'])
        print("shutil.copyfile(%s, %s)\n" %(vars['old_WorkingDir'] + dir + vars['query_seq_file_from_MSA'], vars['working_dir'] + dir + vars['query_seq_file_from_MSA']))
        (Query_Seq, Query_SeqName) = Extract_Fasta_Single_Seq(vars['working_dir'] + dir + vars['query_seq_file_from_MSA'])

    if not query_string in vars:

        vars['query_string'] = Query_SeqName

    # cp HHPred model from previous run
    if not HHPred_Model in vars and os.path.exists(vars['old_WorkingDir'] + dir + vars['HHPred_Model']):

        shutil.copyfile(vars['old_WorkingDir'] + dir + vars['HHPred_Model'], vars['working_dir'] + dir + vars['HHPred_Model'])

    vars['old_running_mode'] = vars['running_mode']
    vars['SubTree_sequences'] = "subtree.fasta"
    # run number and working dir were updated so we can open outputs
    open_log_file()
    try:

        OUTPUT = open(vars['working_dir'] + dir + vars['output_page'], 'a')

    except:

        exit_on_error('sys_error',"create_output_php : could not open the file " + vars['working_dir'] + dir + vars['output_page'] + " for writing.")

    # define running mode!
    if (form['pdb_FILE'] != "" or form[pdb_ID] != "") and form['chain'] != "": # PDB file provided

        vars['running_mode'] = '_mode_pdb_no_msa'

    else: # NO PDB file

        vars['running_mode'] = '_mode_no_pdb_no_msa'

    # add the query seqouence to the sub tree seq
    # read subtree seqs
    try:

        SUB_TREE_SEQS = open(vars['working_dir'] + dir + vars['SubTree_sequences'], 'r')

    except:

        exit_on_error("sys_error","Can't open SUB_TREE_SEQS for reading " + vars['working_dir'] + dir + vars['SubTree_sequences'])

    # check if query seq was included in selected sub tree
    QueryIncluded = "NO"
    NumOfSeq = 0

    line = SUB_TREE_SEQS.readline()
    lines = []
    while line != "":

        match = re.search(r'^>(.*)', line)
        if match:

            NumOfSeq += 1
            tmp = re.sub(r'\s+', "", match.group(1))
            if tmp == vars['query_string'] and QueryIncluded == "NO":

                QueryIncluded = "YES"

        lines.append(line)
        line = SUB_TREE_SEQS.readline()

    SUB_TREE_SEQS.close()

    # new file
    vars['FINAL_sequences'] = "Final_Sequences.fasta"
    try:

        SEQS = open(vars['working_dir'] + dir + vars['FINAL_sequences'], 'w')

    except:

        exit_on_error("sys_error","Can't open SEQS for writing " + vars['working_dir'] + dir + vars['FINAL_sequences'])

    if QueryIncluded == "NO":

        LOG.write("Adding the query sequence to the file of selected sub tree:\n%s\n%s\n" %(vars['guery_string'], Query_Seq))
        SEQS.write(">%s\n%s\n" %(vars['guery_string'], Query_Seq))
        NumOfSeq += 1 # to include the query seq

    for l in lines:

        SEQS.write(l + "\n")

    SEQS.close()

    # create query_final_homolougs.html
    make_sequences_file_HTML(vars['working_dir'] + dir + vars['FINAL_sequences'], vars['working_dir'] + dir + vars['FINAL_sequences_html'], form['AA_or_NUC'])

    print_message_to_output("ConSurf analysis started for <A href=\"%s\" target='_blank'>these %d selected sequences</A> based on the subtree of ConSurf <A href=\"%s/output.php\" target='_balnk'>run number %d</A></font><br>" %(vars['FINAL_sequences_html'], NumOfSeq, vars['old_run_url'], vars['old_run_number']))
    # build the MSA
    vars['protein_MSA'] = "query_msa.aln"
    vars['protein_MSA_clustalw'] = "query_msa_clustalw.aln"
    vars['protein_MSA_FASTA'] = "query_msa_fasta.aln"

    if vars['old_running_mode'] == "_mode_pdb_msa" or vars['old_running_mode'] == "_mode_pdb_msa_tree" or vars['old_running_mode'] == "_mode_msa" or vars['old_running_mode'] == "_mode_msa_tree":

        print_message_to_output("The <A HREF='%s' TARGET=sequences>sequences</A> are automatically aligned using MAFFT-L-INSi procedure</font><br>" %vars['FINAL_sequences_html'])
        form['MSAprogram'] = "MAFFT"

    # create_MSA; NO NEED IT WILL BE CREATED AFTERWARD...

# end sub tree run handeling, from here regular run

else: # regular run, no need to change the run number and out file so we can open them

    vars['run_log_Q'] = CONSURF_CONSTANTS.CONSURF_LOGS_DIR + dir + vars['run_number'] + "_Q.log"
    open_log_file()

    try:

        OUTPUT = open(vars['output_page'], 'a')

    except:

        exit_on_error('sys_error',"create_output_php : could not open the file " + vars['output_page'] + " for writing")

vars['protein_MSA_for_usr'] = ""
vars['run_url_old'] = GENERAL_CONSTANTS.CONSURF_URL + "results" + dir + vars['run_number']

# for nucleotides sequence PDB determine if the nucleotides are D or not
if (form['AA_or_NUC']).upper() == "AA":

    vars['D_Nuc'] = ""

elif (form['AA_or_NUC']).upper() == "NUC": # Nucleic Acid sequence

    if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_pdb_msa_tree":

        vars['D_Nuc'] = is_D_Nuc(vars['pdb_file_name'])
        print("\n\n=========== IS D NUC: '%s' =======\n\n" %vars['D_Nuc'])

    else:

        vars['D_Nuc'] = ""

gradesPE_Output = [] # an array to hold all the information that should be printed to gradesPE
# in each array's cell there is a hash for each line from r4s.res.
# POS: position of that aa in the sequence ; SEQ : aa in one letter ;
# GRADE : the given grade from r4s output ; COLOR : grade according to consurf's scale
residue_freq = {} # for each position in the MSA, detail the residues
position_totalAA = {} # for each position in the MSA, details the total number of residues

# these arrays will hold for each grade, the residues which corresponds to it.
# there are 2 arrays: in the @isd_residue_color, a grade with insufficient data, *, will classify to grade 10
# in the @no_isd_residue_color, the grade will be given regardless of the * mark
# PLEASE NOTE : the [0] position in those arrays is empty, because each position corresponds a color on a 1-10 scale
no_isd_residue_color = [[],[],[],[],[],[],[],[],[],[],[]]
isd_residue_color = [[],[],[],[],[],[],[],[],[],[],[]]
# these variables will be used in the pipe block, for view with FGiJ.
SEQ3D_GRADES_ISD = "" # a string. each position in the string corresponds to that ATOM (from the PDB) ConSurf grade. For Atoms with insufficient data - the grade will be 0
SEQ3D_GRADES = "" # same, only regardeless of insufficient data

#These variables will hold the length of pdb ATOMS and the length of SEQRES/MSA_REFERENCE seq
# The data is filled by cp_rasmol_gradesPE_and_pipe::match_seqres_pdb
LENGTH_OF_SEQRES = 0
LENGTH_OF_ATOM = 0

# In case running on ConSeq Mode, the procedure that finds templates will find templates and fill the vars
Templates_PDB = [] # PDB IDs
Templates_PDB_Descriptions = {}
Templates_Alignments = {} # Hold alignments
# programs
pgp = GENERAL_CONSTANTS.BLASTPGP
rate4s = GENERAL_CONSTANTS.RATE4SITE
rate4s_slow = GENERAL_CONSTANTS.RATE4SITE_SLOW
rate4s_ML = GENERAL_CONSTANTS.RATE4SITE_ML

ModelTest = "python " + consurf_perl_dir + "inferModel.py"

vars['send_email_cmd'] = "perl -w " + consurf_perl_dir + "share/sendEmail.pl"
#vars['send_email_dir'] = tempdir(CLEANUP => 1)

# files
if form['uploaded_TREE'] != "":

    vars['tree_file'] = vars['user_tree_file_name']

else:

    vars['tree_file'] = "TheTree.txt"

vars['tree_viewer_file'] = "treeView.html"

if 'pdb_FILE' in vars and vars['pdb_FILE'] != "": # User PDB

    match = re.search(r'([^\/]+)$', form['pdb_FILE'])
    if match:

        vars['Used_PDB_Name'] = match.group(1)
        # Ofer bug fix: replace spaces with underscores in the used pdb name:
        vars['Used_PDB_Name'] = re.sub(" ", "_", vars['Used_PDB_Name'])

elif 'pdb_ID' in form and form['pdb_ID'] != "":

    vars['Used_PDB_Name'] = form['pdb_ID'] # Given PDB_ID

else:

    vars['Used_PDB_Name'] = "" #ConSeq Mode

vars['r4s_log'] = "r4s.log"
vars['r4s_out'] = "r4s.res"
vars['r4s_slow_log'] = "r4s_slow.log"
vars['atom_positionFILE'] = "atom_pos.txt"
vars['gradesPE'] = "consurf.grades"
vars['scoresFile'] = "consurf.scores"
vars['rasmolFILE'] = "rasmol.scr"
vars['rasmol_isdFILE'] = "isd_rasmol.scr"

vars['rasmolFILE_CBS'] = "rasmol_CBS.scr"
vars['rasmol_isdFILE_CBS'] = "isd_rasmol_CBS.scr"

vars['pipeFile'] = vars['Used_PDB_Name'] + "_consurf" + vars['run_number'] + "_pipe.pdb"
vars['pipeFile_CBS'] = vars['Used_PDB_Name'] + "_consurf" + vars['run_number'] + "_pipe_CBS.pdb" # pipe for color blind friendly
vars['Server_Results_Path'] = vars['run_number']

vars['Confidence_link'] = GENERAL_CONSTANTS.CONSURF_URL + "overview.php#CONFIDENCE"
if (form['AA_or_NUC']).upper() == "NUC":

    vars['Msa_percentageFILE'] = "msa_nucleic_acids_variety_percentage.csv"

else:

    vars['Msa_percentageFILE'] = "msa_aa_variety_percentage.csv"

#chimera files
vars['chimerax_script_for_figure'] = vars['Used_PDB_Name'] + '_consurf_' + vars['run_number'] + '_Figure.chimerax'
vars['chimerax_script_for_figure_isd'] = vars['Used_PDB_Name'] + '_consurf_' + vars['run_number'] + '_Figure_isd.chimerax'

vars['chimera_color_script'] = GENERAL_CONSTANTS.CONSURF_URL + "chimera/chimera_consurf.cmd"
vars['chimera_color_script_CBS'] = GENERAL_CONSTANTS.CONSURF_URL + "chimera/chimera_consurf_CBS.cmd"
vars['chimera_instructions'] = "chimera_instructions.php"

vars['pymol_color_script'] = GENERAL_CONSTANTS.CONSURF_URL + "pyMOL/consurf_new.py"
vars['pymol_color_script_CBS'] = GENERAL_CONSTANTS.CONSURF_URL + "pyMOL/consurf_new_CBS.py"
vars['pymol_instructions'] = "PyMol_instructions.php"
vars['rasmol_instructions'] = "rasmol_instructions.php"

vars['scf_for_chimera'] = "ConSurf_%s_RUN_%s.scf" %(vars['Used_PDB_Name'], vars['run_number'])
vars['header_for_chimera'] = "ConSurf_%s_RUN_%s.hdr" %(vars['Used_PDB_Name'], vars['run_number'])
vars['chimerax_file'] = "ConSurf_%s_RUN_%s.chimerax" %(vars['Used_PDB_Name'], vars['run_number'])

vars['isd_scf_for_chimera'] = "ConSurf_%s_RUN_%s_isd.scf" %(vars['Used_PDB_Name'], vars['run_number'])
vars['isd_header_for_chimera'] = "ConSurf_%s_RUN_%s_isd.hdr" %(vars['Used_PDB_Name'], vars['run_number'])
vars['isd_chimerax_file'] = "ConSurf_%s_RUN_%s_isd.chimerax" %(vars['Used_PDB_Name'], vars['run_number'])


# ConSeq Outputs
vars['Colored_Seq_HTML'] = "Multicolored_Seq.html"
vars['Colored_Seq_PDF'] = "Multicolored_Seq.php"
vars['Colored_MSA_HTML'] = "Multicolored_MSA.html"

vars['insufficient_data_pdb'] = ""
vars['BLAST_PDB_out_file'] = "Blast_vs_PDB"
vars['max_PDB_homologues_to_display'] = GENERAL_CONSTANTS.BLAST_PDB_MAX_HOMOLOGUES_TO_DISPLAY
# Atoms Section with consurf grades instead TempFactor Field
vars['ATOMS_with_ConSurf_Scores'] = vars['Used_PDB_Name'] + "_ATOMS_section_With_ConSurf.pdb"
vars['ATOMS_with_ConSurf_Scores_isd'] = vars['Used_PDB_Name'] + "_ATOMS_section_With_ConSurf_isd.pdb"

vars['pdb_file_with_score_at_TempFactor'] = vars['Used_PDB_Name'] + "_With_Conservation_Scores.pdb"
vars['insufficient_data'] = "no"
vars['All_Outputs_Zip'] = "Consurf_Outputs_" + vars['run_number'] + ".zip"
vars['pdb_db'] = GENERAL_CONSTANTS.CULLED_PDB

## mode : no msa - with PDB or without PDB

if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_no_pdb_no_msa":

    # we save to copies of the msa, one in fasta format and another in clustal format.
    vars['msa_fasta'] = "msa_fasta.aln"
    vars['msa_clustal'] = "msa_clustal.aln"
    create_MSA()
    Update_Progress(vars['PROGRESS_REPORT'], "Align sequences")
    vars['msa_SEQNAME'] = vars['query_string']

## mode : include MSA

elif vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_msa" or vars['running_mode'] == "_mode_pdb_msa_tree" or vars['running_mode'] == "_mode_msa_tree":

    vars['msa_SEQNAME'] = form['msa_SEQNAME']

if form['matrix'] == "BEST":

    best_matrix = find_best_substitution_model(vars['msa_fasta'],  form['AA_or_NUC'])
    form['matrix'] = best_matrix

    if best_matrix == "JC_Nuc":

        best_matrix = "JC"

    print_message_to_output ("The best evolutionary model was selected to be: " + best_matrix + ". See details <A href=\"ModelTest/model_selection.out\">here</A><br>")
    Update_Progress(vars['PROGRESS_REPORT'], "Select best evolutionary model")

run_rate4site(vars['query_string'])
Update_Progress(vars['PROGRESS_REPORT'], "Calculate conservation scores")
assign_colors_according_to_r4s_layers(gradesPE_Output, vars['r4s_out'])
vars['num_of_seqs_in_MSA'] = read_residue_variety(residue_freq, position_totalAA)

if form['AA_or_NUC'] == "AA":

    print_residue_precentage()

elif form['AA_or_NUC'] == "NUC":

    print_nucleotide_precentage()

## mode : include pdb

# in order to create 3D outputs, we need to compare the ATOM to the sequence from rate4site
if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_pdb_msa" or vars['running_mode'] == "_mode_pdb_msa_tree":

    create_atom_position_file(form['chain'], vars['pdb_file_name'], vars['atom_positionFILE'], form['AA_or_NUC'], 0) # this file will be used later to create the output which aligns rate4site sequence with the ATOM records
    r4s2pdb = {} # key: poistion in SEQRES/MSA, value: residue name with position in atom (i.e: ALA22:A)

    if 'seqres_or_msa_seq_with_gaps' in vars and len(vars['seqres_or_msa_seq_with_gaps']) > 0:

        [LENGTH_OF_SEQRES, LENGTH_OF_ATOM] = match_pdb_to_seq(r4s2pdb, form['chain'], vars['seqres_or_msa_seq_with_gaps'], vars['ATOM_seq_with_gaps'], vars['atom_positionFILE'], form['AA_or_NUC'])

    else: # NO SEQRES

        fill_r4s2pdb(r4s2pdb) # fill length_of_atom
        LENGTH_OF_SEQRES = 0

    # The following routine, apart from creating the file "consurf.grades" also collects the information in order to create
    # the RasMol scripts and a variable which will be used in the "pipe" file, that holds a string with the grades.
    # In the pipe file this var is called: seq3d_grades_isd and seq3d_grades
    [SEQ3D_GRADES_ISD, SEQ3D_GRADES] = create_gradesPE_ConSurf(vars['gradesPE'], vars['scoresFile'], r4s2pdb, gradesPE_Output, residue_freq, no_isd_residue_color, isd_residue_color, form['AA_or_NUC'])

    create_rasmol(form['chain'], vars['rasmolFILE'], vars['rasmol_isdFILE'], vars['rasmolFILE_CBS'], vars['rasmol_isdFILE_CBS'], no_isd_residue_color, isd_residue_color) #This will create the 2 rasmol scripts (one with isd and one without)

    # This will create the pipe file for FGiJ
    vars['identical_chains'] = create_pipe_file(vars['pipeFile'], vars['pipeFile_CBS'], SEQ3D_GRADES, SEQ3D_GRADES_ISD, isd_residue_color, no_isd_residue_color, LENGTH_OF_SEQRES, LENGTH_OF_ATOM, vars['pdb_file_name'], form['chain'], (vars['Used_PDB_Name']).upper())
    Hash_Json = {'scoresFile' : vars['scoresFile'], 'pdbFile' : vars['pdb_file_name'], 'mmtfFile' : vars['mmtfFile'], 'jobTitle' : form['job_title'], 'pdbId' : form['pdb_ID'], 'chainId' : form['chain'], 'identicalChains' : vars['identical_chains']}
    saveJobInfoJson("job_info.json", Hash_Json)

    replace_TmpFactor_Consurf_Scores(form['chain'], vars['pdb_file_name'], vars['gradesPE'], vars['ATOMS_with_ConSurf_Scores'], vars['ATOMS_with_ConSurf_Scores_isd']) #Create ATOMS section and replace the TempFactor Column with the ConSurf Grades (will create also isd file if relevant)
    replace_TmpFactor_Rate4Site_Scores(form['chain'], vars['pdb_file_name'], vars['gradesPE'], vars['pdb_file_with_score_at_TempFactor']) # Create PDB file and replace the TempFactor Column with the Rate4Site grades
    create_chimera(vars['msa_clustal'], vars['working_dir'], vars['r4s_out'], vars['scf_for_chimera'], vars['header_for_chimera'], vars['isd_scf_for_chimera'], vars['isd_header_for_chimera'], vars['ATOMS_with_ConSurf_Scores'], vars['chimerax_file'], vars['ATOMS_with_ConSurf_Scores_isd'], vars['isd_chimerax_file'], vars['chimerax_script_for_figure'], vars['chimerax_script_for_figure_isd'], vars['chimera_instructions'], vars['identical_chains'], form['chain'], vars['pdb_file_name']) #This Will create the script for chimera coloring
    create_pymol_instructions(vars['insufficient_data'], vars['pymol_instructions'], vars['Used_PDB_Name'], form['chain'], vars['ATOMS_with_ConSurf_Scores'], vars['ATOMS_with_ConSurf_Scores_isd'], vars['identical_chains'], vars['working_dir'], vars['chimerax_script_for_figure'], vars['chimerax_script_for_figure_isd'], vars['pdb_file_name']) # This will create the PyMol instruction page

    # create 'ConSeq' output for protein structure burried exposed
    LOG.write("Get_NACSES_buried_Exposed(%s, %s, %s, %s)\n" %(vars['pdb_file_name'], form['chain'], BURRIED_CUTOFF, EXPOSED_CUTOFF))
    [NACSES_RESULTS, Pos_Solv_Acc_Pred_NACSES_hashRef] = Get_NACSES_buried_Exposed(vars['pdb_file_name'], form['chain'], BURRIED_CUTOFF, EXPOSED_CUTOFF)
    LOG.write("Get_NACSES_buried_Exposed:%s\n" %NACSES_RESULTS)
    ans = cp_rasmol_gradesPE_and_pipe.create_rasmol_page(vars['rasmol_instructions'], vars['identical_chains'], vars['run_number'], vars['working_dir'], form['chain'], vars['ATOMS_with_ConSurf_Scores'], vars['ATOMS_with_ConSurf_Scores_isd'], vars['chimerax_script_for_figure'], vars['chimerax_script_for_figure_isd'], vars['rasmolFILE'], vars['rasmol_isdFILE'], vars['pdb_file_name'], vars['Used_PDB_Name'], vars['rasmolFILE_CBS'], vars['rasmol_isdFILE_CBS']) # Used_PDB_Name is just for file download name purposes
    if ans != "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_rasmol_page FAILED: " + ans)

    if form['pdb_ID'] != "" and form['AA_or_NUC'] == "AA":

        # get complex from PISA
        vars['PISA_PDB'] = (form['pdb_ID']).upper() + "_PISA1.pdb"
        ans2 = Get_PISA_Complex(form['pdb_ID'], vars['PISA_PDB'])
        if ans2 != "OK":

            LOG.write("No PISA complex for %s\n" %form['pdb_ID'])
            del vars['PISA_PDB']

        else:

            LOG.write("PISA complex model was saved to %s\n" %(vars['PISA_PDB']))
            LOG.write("Project_ConSurf_On_Model(%s, %s, %s, %s, %s, %s, %s, %s)\n" %(vars['PISA_PDB'], form['chain'], vars['working_dir'], vars['protein_seq_string'], "r4s.res", vars['msa_clustal'], "PISA first model", form['AA_or_NUC']))
            [HTML_TXT_FOR_PISA_FILES, HTML_TXT_FOR_PISA_FGIJ] = Project_ConSurf_On_Model(vars['PISA_PDB'], form['chain'], vars['working_dir'], vars['protein_seq_string'], "r4s.res", vars['msa_clustal'], "PISA first model", form['AA_or_NUC'])
            if HTML_TXT_FOR_PISA_FILES != "":

                HTML_TXT_FOR_PISA_MODEL += HTML_TXT_FOR_PISA_FGIJ + HTML_TXT_FOR_PISA_FILES

## mode : ConSeq - NO PDB

if vars['running_mode'] == "_mode_msa" or vars['running_mode'] == "_mode_no_pdb_no_msa" or vars['running_mode'] == "_mode_msa_tree": #ConSeq Mode AA

    Pos_Solv_Acc_Pred = {} # key: poistion in SEQRES/MSA, value: PACC Solv Acc Pred - relevant only for AA case
    if form['AA_or_NUC'] == "AA":

        vars['protein_MSA_hssp'] = vars['msa_fasta'] + ".hssp"
        vars['Solv_ACC_Pred'] = "Solv_Acc_Pred.PACC"

        """
        msa_format = MSA_parser.determine_msa_format(vars['protein_MSA'])
        if msa_format[0] == "err":

            msa_info_msg = "<a href =\"" + GENERAL_CONSTANTS.MSA_FORMATS + "\">Read more on MSA formats</a><br />\n"
            exit_on_error('user_error', "The uploaded <a href=\"" + vars['user_msa_file_name'] + "\">MSA file</a> is not in one of the formats supported by ConSurf: NBRF/PIR, Pearson (Fasta), Nexus, Clustal, GCG/MSF.<br /><br />\nPlease check the following items and try to run ConSurf again:<br />\n1. The file should be saved as plain text (e.g. file type 'txt' in windows or 'MS-Dos' from Word in Mac).<br />\n2. The file should not contain unnecessary characters (You can check it with 'Notepad' editor).<br />\n3. The same sequence name must not be repeated more then once.<br />\n" + msa_info_msg)

        else:

            LOG.write("determine_msa_format : MSA format is : " + msa_format[1] + "\n")
            vars['msa_format'] = msa_format[1]
        """

        vars['HHPred_Model'] = "HHPred_Model.ent"
        blast_vs_PDB(vars['protein_seq'], vars['BLAST_PDB_out_file'])
        #Find_Good_Templates(vars['protein_seq_string'], vars['BLAST_PDB_out_file'])
        parse_blast_alignment()
        Update_Progress(vars['PROGRESS_REPORT'], "Search for 3D structure for the protein sequence")
        # Store data for projection on PDB structures
        store_data()
        if 'Modeller_Test' in vars and vars['Modeller_Test'] == "OK":

            # run HHPred
            LOG.write("Run_HHPred(%s, %s, %s)\n" %(vars['protein_seq'], "HH_Pred/", vars['HHPred_Model']))
            if not os.path.exists(vars['HHPred_Model']): # JUST TO SAVE TIME...

                Run_HHPred(vars['protein_seq'], "HH_Pred/", vars['HHPred_Model'])
                vars['base_HH_PRED_ModelName'] = fileparse(vars['protein_seq'], r'\.[^.]*')

            Update_Progress(vars['PROGRESS_REPORT'], "Predict 3D structure using HHPred and MODELLER")
            if os.path.exists(vars['HHPred_Model']): # HHPred OK, project scores on model

                [HTML_TXT_FOR_MODEL_SECTION, HTML_TXT_FOR_MODEL_SECTION_FGIJ] = Project_ConSurf_On_Model(vars['HHPred_Model'], "NONE", vars['working_dir'], vars['protein_seq'], "r4s.res", vars['msa_clustal'], "HHPred model", form['AA_or_NUC'])

        if os.path.exists(vars['working_dir'] + vars['HHPred_Model']): # b/e/ based on HHPred_model

            LOG.write("Get_NACSES_buried_Exposed(%s, \"\", %s, %s)" %(vars['HHPred_Model'], BURRIED_CUTOFF, EXPOSED_CUTOFF))
            [NACSES_RESULTS, Pos_Solv_Acc_Pred_NACSES_hashRef] = Get_NACSES_buried_Exposed(vars['HHPred_Model'], "", BURRIED_CUTOFF, EXPOSED_CUTOFF)
            LOG.write("Get_NACSES_buried_Exposed:" + NACSES_RESULTS + "\n")
            if NACSES_RESULTS == "OK":

                create_gradesPE_ConSeq(Pos_Solv_Acc_Pred_NACSES_hashRef)
                ConSeq_gradesPE_and_Outputs.ConSeq_HTML_Output(gradesPE_Output, Pos_Solv_Acc_Pred_NACSES_hashRef, vars['Colored_Seq_HTML'], "HHPred 3D model")
                Protein_Length = len(gradesPE_Output) + 1
                ConSeq_gradesPE_and_Outputs.ConSeq_PDF_Output(vars['gradesPE'], vars['Colored_Seq_PDF'], Protein_Length, "HHPred 3D model")

            else: # NACSESS failed predict b/e by neural-network algorithm

                predict_solvent_accesibility(vars['msa_clustal'], vars['protein_MSA_hssp'], vars['Solv_ACC_Pred']) #runs the PACC algorithm to calculate burried/exposed
                read_Solv_Acc_Pred(Pos_Solv_Acc_Pred)
                create_gradesPE_ConSeq(Pos_Solv_Acc_Pred)
                ConSeq_gradesPE_and_Outputs.ConSeq_HTML_Output(gradesPE_Output, Pos_Solv_Acc_Pred, vars['Colored_Seq_HTML'], "neural-network algorithm")
                Protein_Length = len(gradesPE_Output) + 1
                ConSeq_gradesPE_and_Outputs.ConSeq_PDF_Output(vars['gradesPE'], vars['Colored_Seq_PDF'], Protein_Length, "neural-network algorithm")

        else: # v/e based on prediction

            predict_solvent_accesibility(vars['msa_clustal'], vars['protein_MSA_hssp'], vars['Solv_ACC_Pred']) #runs the PACC algorithm to calculate burried/exposed
            read_Solv_Acc_Pred(Pos_Solv_Acc_Pred)
            create_gradesPE_ConSeq(Pos_Solv_Acc_Pred)
            ConSeq_gradesPE_and_Outputs.ConSeq_HTML_Output(gradesPE_Output, Pos_Solv_Acc_Pred, vars['Colored_Seq_HTML'], "neural-network algorithm")
            Protein_Length = len(gradesPE_Output) + 1
            ConSeq_gradesPE_and_Outputs.ConSeq_PDF_Output(vars['gradesPE'], vars['Colored_Seq_PDF'], Protein_Length, "neural-network algorithm")

    elif (form['AA_or_NUC']).upper() == "NUC": # Nucleic Acid sequence

        vars['D_Nuc'] = ""
        create_gradesPE_ConSeq_Nuc()
        ConSeq_gradesPE_and_Outputs.ConSeq_HTML_Output(gradesPE_Output, Pos_Solv_Acc_Pred, vars['Colored_Seq_HTML'], True) # The last argument ("yes") tell us that only the sequence should be printed
        Protein_Length = len(gradesPE_Output) + 1
        ConSeq_gradesPE_and_Outputs.ConSeq_NUC_PDF_Output(vars['gradesPE'], vars['Colored_Seq_PDF'], Protein_Length)
        # predict RNA SS
        [vars['RNA_ss_ps'], vars['RNA_dp_ps'], vars['RNA_ss_pdf']] = Predict_RNA_SS(vars['protein_seq'], vars['working_dir'])
        if os.path.exists(vars['working_dir'] + vars['RNA_ss_ps']):

            [vars['RNA_ss_colored_by_ConSurf_ps'], vars['RNA_ss_colored_by_ConSurf_pdf']] = Color_RNA_SS_by_ConSurf(vars['RNA_ss_ps'], vars['gradesPE'], vars['working_dir'])
            [vars['RNA_ss_colored_by_ConSurf_ps_CBS'], vars['RNA_ss_colored_by_ConSurf_pdf_CBS']] = Color_RNA_SS_by_ConSurf(vars['RNA_ss_ps'], vars['gradesPE'], vars['working_dir'], "CBS")

    vars['scf_for_chimera'] = "consurf_" + vars['run_number'] + ".scf"
    vars['header_for_chimera'] = "consurf_" + vars['run_number'] + ".hdr"
    vars['chimerax_file'] = "consurf_" + vars['run_number'] + ".chimerax"
    vars['isd_scf_for_chimera'] = "consurf_" + vars['run_number'] + "_isd.scf"
    vars['isd_header_for_chimera'] = "consurf_" + vars['run_number'] + "_isd.hdr"
    vars['isd_chimerax_file'] = "consurf_" + vars['run_number'] + "_isd.chimerax"

    create_chimera_align_tree()

# create wasabi xml
vars['WASABI_XML'] = "MSA_and_Tree.xml"
LOG.write("Create file for wasabi: wasabiUtils.createWasabiXml(%s, %s, %s)" %(vars['msa_fasta'], vars['tree_file'], vars['WASABI_XML']))
wasabiUtils.createWasabiXml(vars['msa_fasta'], vars['tree_file'], vars['WASABI_XML'])

Prepare_tree_View()
Create_Colored_MSA()

if (form['AA_or_NUC']).upper() == "NUC": # Nucleic Acid sequence

    # predict RNA SS
    [vars['RNA_ss_ps'], vars['RNA_dp_ps'], vars['RNA_ss_pdf']] = Predict_RNA_SS(vars['protein_seq'], vars['working_dir'])

    if os.path.exists(vars['RNA_ss_ps']):

        [vars['RNA_ss_colored_by_ConSurf_ps'], vars['RNA_ss_colored_by_ConSurf_pdf']] = Color_RNA_SS_by_ConSurf(vars['RNA_ss_ps'], vars['gradesPE'], vars['working_dir'], "legacy")
        [vars['RNA_ss_colored_by_ConSurf_ps_CBS'], vars['RNA_ss_colored_by_ConSurf_pdf_CBS']] = Color_RNA_SS_by_ConSurf(vars['RNA_ss_ps'], vars['gradesPE'], vars['working_dir'], "CBS")

zip_all_outputs()
Update_Progress(vars['PROGRESS_REPORT'], "Project conservation scores onto the molecule")

## Arrange The HTML Output File

OUTPUT.write("\n<H1><center><a name=finish class='unlink'>ConSurf calculation is finished:</a></center></H1>\n")

try:

    IN = open(vars['gradesPE'], 'r')

except:

    exit_on_error("sys_error", "Can't open IN " + vars['gradesPE'])

Pos = 0
ISD_Pos = 0
line = IN.readline()
while line != "":

    line = line.rstrip()
    line = re.sub(r'^\s+|\s+$', "", line)
    phrases = re.split(r'\s+', line) # remember to check if input could be float or negative
    if (phrases[0]).isdigit():

        Pos += 1
        match = re.match(r'\*', line)
        if match:

            ISD_Pos += 1

    line = IN.readline()

IN.close()
if ISD_Pos > 30 or ISD_Pos / Pos > 0.1:

    type = ""
    if (form['AA_or_NUC']).upper() == "AA":

        type = "residues"

    elif (form['AA_or_NUC']).upper() == "NUC":

        type = "nucleic acids"

    OUTPUT.write("\n<H3><center><font color=red><a name=finish class='unlink'><i>Warning: %s of %s %s have unreliable conservation scores due to <A href='%s'>insufficient data</A> in the multiple sequence alignment</center></H3> </i></font>" %(ISD_Pos, Pos, type, vars['Confidence_link']))
"""
try:

    OUTPUT.write(vars['working_dir'] + dir + vars['output_page'], 'w')

except:

    exit_on_error('sys_error', "create_output_php : could not open the file " + vars['working_dir'] + dir + vars['output_page'] + " for writing.")
"""
OUTPUT.write("<h4 class=output_title>Final Results</h4>\n")
if form['pdb_FILE'] == "" and form['pdb_ID'] == "" and form['chain'] == "":

    # No PDB is available - Uniq ConSeq Outputs
    print_message_to_output ("The query sequence colored according to the conservation scores <A HREF='%s' TARGET=results>(HTML)</A><A HREF='%s' TARGET=results_PDF>(PDF)</A></font><br>" %(vars['Colored_Seq_HTML'], vars['Colored_Seq_PDF']))
    print_message_to_output ("<A HREF='%s' TARGET=Conservation_MSA>Multiple Sequence Alignment Color-Coded by Conservation</A>" %vars['Colored_MSA_HTML'])
    print_message_to_output ("<A HREF='%swasabi/?url=%s' TARGET='WASABI'>View MSA and phylogenetic tree using WASABI and run ConSurf on sub-tree</a>" %(GENERAL_CONSTANTS.CONSURF_URL, vars['run_url'] + dir + vars['WASABI_XML']))
    # CHIMERA
    OUTPUT.write("\n<ul><li><a href = \"")
    if vars['insufficient_data'] == "yes" and os.path.exists(vars['isd_chimerax_file']) and os.path.getsize(vars['isd_chimerax_file']) != 0 and os.path.exists(vars['isd_header_for_chimera']) and os.path.getsize(vars['isd_header_for_chimera']) != 0:

        OUTPUT.write(vars['run_url_old'] + dir + vars['isd_chimerax_file'])
        OUTPUT.write("\"  type=\"application/x-chimerax\">Multiple Sequence Alignment Color-Coded by Conservation and the (neighbor-joining) ConSurf tree</a></b> (<a href=\"" + GENERAL_CONSTANTS.CHIMERA_DOWNLOAD + "\">Download Chimera</a>; required)\n")
        OUTPUT.write("<a href=\"javascript:toggle_help('chimerax');\" title=\"click for help\"><img src=\"<?= \$consurf2016_url; ?>images/i_quest.gif\" border=\"0\"></a><div id=\"chimerax\">To display the file:<ol><li>Save the chimerax file in your computer<li>Open chimera program<li>Load (File->Open...) the chimerax file from your computer in chimera</div><br>")
        OUTPUT.write(" If you wish to avoid seeing the <a href='%s'>insufficient data</a>, use <a href = \"%s\" type=\"application/x-chimerax\">this link</a> please.<br>" %(vars['Confidence_link'], vars['run_url_old'] + dir + vars['chimerax_file']))

    else:

        OUTPUT.write(vars['chimerax_file'])
        OUTPUT.write("\"  type=\"application/x-chimerax\">Multiple Sequence Alignment Color-Coded by Conservation and the (neighbor-joining) ConSurf tree </a> (<a href=\"" + GENERAL_CONSTANTS.CHIMERA_DOWNLOAD + "\">Download Chimera</a>; required)\n")
        OUTPUT.write("<a href=\"javascript:toggle_help('chimerax');\" title=\"click for help\"><img src=\"<?= \$consurf2016_url; ?>images/i_quest.gif\" border=\"0\"></a><div id=\"chimerax\">To display the file:<ol><li>Save the chimerax file in your computer<li>Open chimera program<li>Load (File->Open...) the chimerax file from your computer in chimera</div>\n")

    OUTPUT.write("</li></ul>\n")
    if form['AA_or_NUC'] == "AA":

        print_message_to_output("<A HREF='" + vars['gradesPE'] + "' TARGET=Conservation_window>Amino Acid Conservation Scores, Confidence Intervals and Conservation Colors</A>")

    elif form['AA_or_NUC'] == "NUC":

        print_message_to_output("<A HREF='" + vars['gradesPE'] + "' TARGET=Conservation_window>Nucleic Acid Conservation Scores, Confidence Intervals and Conservation Colors</A>")
        if RNA_ss_colored_by_ConSurf_pdf in vars and RNA_ss_colored_by_ConSurf_ps in vars and os.path.exists(vars['working_dir'] + vars['RNA_ss_colored_by_ConSurf_pdf']):

            print_message_to_output ("Secondary structure prediction for the RNA molecule (by RNAFold) colored by ConSurf scores (<A href=\"%s\" target=\"_blank\">ps file</A> <A href=\"%s\" target=\"_blank\">pdf file</A>)" %(vars['RNA_ss_colored_by_ConSurf_ps'], vars['RNA_ss_colored_by_ConSurf_pdf']))
            print_message_to_output ("Secondary structure prediction for the RNA molecule (by RNAFold) colored by ConSurf scores [color blind friendly scale] (<A href=\"%s\" target=\"_blank\">ps file</A> <A href=\"%s\" target=\"_blank\">pdf file</A>)" %(vars['RNA_ss_colored_by_ConSurf_ps_CBS'], vars['RNA_ss_colored_by_ConSurf_pdf_CBS']))

    print_message_to_output("<B>Download all ConSurf outputs in a <A HREF='" + vars['All_Outputs_Zip'] + "'>click!</A><br></B>")

else:

    # UNIQ ConSurf Results
    pipeUrlNoSuffix = dir + vars['Server_Results_Path'] + dir + vars['Used_PDB_Name'] + "_consurf" + vars['run_number']
    OUTPUT.write(buildHtml_3dViewerLinks(vars['run_number'], "", pipeUrlNoSuffix))
    OUTPUT.write("\n<ul><li><a href=\"")

    # if there was insufficient data, we add also the chimerax file showing insufficient data
    if vars['insufficient_data'] == "yes" and os.path.exists(vars['working_dir'] + vars['isd_chimerax_file']) and os.path.getsize(vars['working_dir'] + vars['isd_chimerax_file']) != 0 and os.path.exists(vars['working_dir'] + vars['isd_header_for_chimera']) and os.path.getsize(vars['isd_header_for_chimera']) != 0:

        OUTPUT.write(vars['run_url_old'] + dir + vars['isd_chimerax_file'])
        OUTPUT.write("\" type=\"application/x-chimerax\">View ConSurf results</a> with Chimera (<a href=\"" + GENERAL_CONSTANTS.CHIMERA_DOWNLOAD + "\">Download Chimera</a>)<br>\n")
        OUTPUT.write(" If you wish to view the molecule and avoid the <a href='%s'>insufficient data</a>, use <a href = \"%s\" type=\"application/x-chimerax\" TARGET=\"_blank\">this link</a> please.</li></ul>" %(vars['Confidence_link'], vars['run_url_old'] + dir + vars['chimerax_file']))

    else:

        OUTPUT.write(vars['chimerax_file'])
        OUTPUT.write("\"  type=\"application/x-chimerax\">View ConSurf results</a></b> with Chimera (<a href=\"" + GENERAL_CONSTANTS.CHIMERA_DOWNLOAD + "\">Download Chimera</a>)</li></ul>\n")

    print_message_to_output("<A HREF='%swasabi/?url=%s' TARGET='WASABI'>View MSA and phylogenetic tree using WASABI and run ConSurf on sub-tree</a>" %(GENERAL_CONSTANTS.CONSURF_URL, vars['run_url'] + dir + vars['WASABI_XML']))
    print_message_to_output("<A HREF='%s' TARGET=Conservation_MSA>Multiple Sequence Alignment Color-Coded by Conservation</A>" %vars['Colored_MSA_HTML'])
    if form['AA_or_NUC'] == "NUC":

        print_message_to_output("<A HREF='%s' TARGET=Conservation_window>Nucleic Acid Conservation Scores, Confidence Intervals and Conservation Colors</A>" %vars['gradesPE'])
        if RNA_ss_colored_by_ConSurf_pdf in vars and RNA_ss_colored_by_ConSurf_ps in vars and os.path.exists(vars['RNA_ss_colored_by_ConSurf_pdf']):

            print_message_to_output("Secondary structure prediction for the RNA molecule (by RNAFold) colored by ConSurf scores (<A href=\"%s\" target=\"_blank\">ps file</A> <A href=\"%s\" target=\"_blank\">pdf file</A>)" %(vars['RNA_ss_colored_by_ConSurf_ps'], vars['RNA_ss_colored_by_ConSurf_pdf']))
            print_message_to_output("Secondary structure prediction for the RNA molecule (by RNAFold) colored by ConSurf scores [color blind friendly scale] (<A href=\"%s\" target=\"_blank\">ps file</A> <A href=\"%s\" target=\"_blank\">pdf file</A>)" %(vars['RNA_ss_colored_by_ConSurf_ps_CBS'], vars['RNA_ss_colored_by_ConSurf_pdf_CBS']))

    else:

        print_message_to_output("<A HREF='%s' TARGET=Conservation_window>Amino Acid Conservation Scores, Confidence Intervals and Conservation Colors</A>" %vars['gradesPE'])

    print_message_to_output("<B>Download all ConSurf outputs in a <A HREF='%s'>click!</A><br></B>" %vars['All_Outputs_Zip'])

    OUTPUT.write("<h4 class=output_title>PDB Files</h4>")
    print_message_to_output("<A HREF='%s' TARGET=PDB_window> PDB File with Conservation Scores in the tempFactor field</A>" %vars['pdb_file_with_score_at_TempFactor'])
    print_message_to_output("<A HREF='%s' TARGET=PDB_window> PDB File with ConSurf Results in its Header, for FirstGlance in Jmol</A>; [<A href='%s'>with color-blind friendly scale</A>]<br>" %(vars['Used_PDB_Name'] + "_consurf" + vars['run_number'] + "_pipe.pdb", vars['Used_PDB_Name'] + "_consurf" + vars['run_number'] + "_pipe_CBS.pdb"))

    OUTPUT.write("<h4 class=output_title>Create high resolution figures</h4>\n")
    print_message_to_output("<A HREF='%s' TARGET=Chimera_HighRes_Instruct_window>Follow the instructions to produce a Chimera figure</A><font size=\"-1\"> (For users of Chimera)</font><br>" %vars['chimera_instructions'])
    print_message_to_output("<A HREF='%s' TARGET=PyMol_HighRes_Instruct_window>Follow the instructions to produce a PyMol figure</A><font size=\"-1\"> (For users of PyMol)</font><br>" %vars['pymol_instructions'])
    print_message_to_output("<A HREF='%s' TARGET=rasmol_HighRes_Instruct_window>Follow the instructions to produce a RasMol figure</A><font size=\"-1\"> (For users of RasMol)</font><br>" %vars['rasmol_instructions'])

# Create a high resolution figure
#    Produce a figure of your protein, colored by ConSurf's colors:
#
#     * Follow the instructions to produce a PyMOL figure (For users of PyMOL)
#     * Follow the instructions to produce a Chimera figure (For users of Chimera)

# COMMON OUTPUTS
if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_no_pdb_no_msa":

    # The Server create the MSA
    blast_algorithm = ""
    if form['HomologSearchAlg'] == "BLAST":

        if form['AA_or_NUC'] == "NUC":

            blast_algorithm = "BLAST"

        else:

            blast_algorithm = "PSI-BLAST"

    elif form['HomologSearchAlg'] == "HMMER":

        blast_algorithm = "HMMER"

    else:

        blast_algorithm = "CSI-BLAST"

    OUTPUT.write("<h4 class=output_title>Sequence Data</h4>")

    if blast_algorithm != "HMMER":

        print_message_to_output("<A HREF='%s' TARGET=Blast_window>%s output</A> (%s hits with E-values and pairwise alignments)<br>" %(vars['BLAST_out_file'], blast_algorithm, blast_algorithm))

    else:

        print_message_to_output("<A HREF='%s' TARGET=Blast_window>%s output</A> (%s hits with E-values and pairwise alignments)<br>" %(vars['HMMER_out_file'], blast_algorithm, blast_algorithm))

    evalue_max_str = ''
    if 'MaxFinalSequenceEvalue' in vars:

        evalue_max_str = ", max E-value = %s, median E-value = %s" %(vars['MaxFinalSequenceEvalue'], vars['MedianFinalSequenceEvalue'])

    print_message_to_output("<A HREF= '%s' TARGET=Homologues_window>Sequences Used</A> (displayed in FASTA format, linked to sequence data-base)%s<br>" %(vars['FINAL_sequences_html'], evalue_max_str))

OUTPUT.write("<h4 class=output_title>Alignment</h4>\n")
if vars['protein_MSA_for_usr'] != "":

    print_message_to_output("<A HREF= '%s' TARGET=Alignment_window>Multiple Sequence Alignment </A><br>" %vars['protein_MSA_for_usr'])

else:

    print_message_to_output("<A HREF= '%s' TARGET=Alignment_window>Multiple Sequence Alignment </A><br>" %vars['msa_fasta'])

extract_and_print_diversity_matrix()
if form['AA_or_NUC'] == "NUC":

    print_message_to_output("<A HREF= '%s' TARGET=Msa_percentage_window>Nucleic acids variety per position in the MSA</A></font><font size=-1> (The table is best viewed with an editor that respects Comma-Separated Values)</font><br>" %vars['Msa_percentageFILE'])

else:

    print_message_to_output("<A HREF= '%s' TARGET=Msa_percentage_window>Residue variety per position in the MSA</A><font size=-1> (The table is best viewed with an editor that respects Comma-Separated Values)</font><br>" %vars['Msa_percentageFILE'])

OUTPUT.write("<h4 class=output_title>Phylogenetic Tree</h4>\n")
print_message_to_output("<A HREF= '%s' TARGET=Tree_window>Phylogenetic Tree in Newick format</A><br>" %vars['tree_file'])
print_message_to_output("<A HREF= '%s' TARGET=TreeView_window>View Phylogenetic Tree </A><br>" %vars['tree_viewer_file'])

if vars['running_mode'] == "_mode_msa" or vars['running_mode'] == "_mode_no_pdb_no_msa" or vars['running_mode'] == "_mode_msa_tree": #ConSeq Mode

    # Projection on PDB for ConSeq AA
    if form['AA_or_NUC'] == "AA":

        Num_Of_Templates = len(Templates_PDB_Descriptions)
        if Num_Of_Templates > 0:

            OUTPUT.write("<br><br><br>\n<span class=\"PrintOutH\"><div id=\"PrintOutH\">Project ConSurf results on known structures similar to your sequence</div></span>\n")
            OUTPUT.write("<div class=\"mH\"></div>\n")
            OUTPUT.write("When looking for known structures that share high similarity with your sequence, we have found %d relevant structures.<br> Press <span class=\"mH1\" onclick=\"toggleMenu('TEMPLATES_SECTION')\">here</span> to show the list.<br>" %Num_Of_Templates)
            OUTPUT.write("<div id='TEMPLATES_SECTION'>\n")
            print_templates_projection_links()

        # print Model info
        if os.path.exists(vars['HHPred_Model']): # HHPred OK

            OUTPUT.write("<br><br><br>\n<span class=\"PrintOutH\"><div id=\"PrintOutH\">Template-based 3D model for the query sequence</div></span>\n<div class=\"mH\"></div>")
            OUTPUT.write(HTML_TXT_FOR_MODEL_SECTION_FGIJ)
            print_message_to_output("<A HREF='%s' TARGET=3Dmodel>3D model predicted by HHPred</A><br>" %vars['HHPred_Model']);
            print_message_to_output("<A HREF='HH_Pred/%s.pir' TARGET=AlnHHPred>The alignment between the protein query and available templates (pir format)</A><br>" %vars['base_HH_PRED_ModelName']);
            print_message_to_output("<A HREF='HH_Pred/%s.hhr' TARGET=hhr_HHPred>Templates report (hhr)</A><br>" %vars['base_HH_PRED_ModelName']);
            print_message_to_output("<A HREF='HH_Pred/%s.py' TARGET=py_HHPred>Modeller python script</A><br>" %vars['base_HH_PRED_ModelName']);
            print_message_to_output("<A HREF='HH_Pred/%s.new.rsr' TARGET=rsr_HHPred>Modeller HHPred restrains</A><br>" %vars['base_HH_PRED_ModelName'])
            OUTPUT.write(HTML_TXT_FOR_MODEL_SECTION)

if 'PISA_PDB' in vars: # PISA file was downloaded

    if HTML_TXT_FOR_PISA_MODEL != "":

        OUTPUT.write("<br><br><br>\n<span class=\"PrintOutH\"><div id=\"PrintOutH\">Project ConSurf scores on the protein most probable assembly</div></span>\n<div class=\"mH\"></div>")
        OUTPUT.write(HTML_TXT_FOR_PISA_MODEL + "\n")

    else:

        OUTPUT.write("<br><br><br>\n<span class=\"PrintOutH\"><div id=\"PrintOutH\">Cannot calculate PISA model results</div></span>\n<div class=\"mH\"></div>")

OUTPUT.write("""<br>
<?php
    include(\"%stemplates/footer.tpl\");
?>""" %server_html_dir)
OUTPUT.close()
update_output_that_run_finished(vars['output_page'])
send_finish_email_to_user()
OUTPUT.close()
LOG.close()
system("cp %s %s" %(vars['output_page'], "output_with_form.php")) # For those bookmark the page with the form
END_FLAG = open(vars['CONSURF_qsub_name'] + ".END_OK", 'w')
END_FLAG.close()

## functions

def zip_all_outputs():

    # -q option means it will create no output, also > /dev/null
    cmd = "zip -q %s %s %s %s %s > /dev/null" %(vars['All_Outputs_Zip'], vars['Colored_MSA_HTML'], vars['tree_file'], vars['gradesPE'], vars['Msa_percentageFILE'])

    if form['pdb_FILE'] == "" and form['pdb_ID'] == "" and form['chain'] == "": # UNIQ TO CONSEQ

        cmd += " " + vars['Colored_Seq_HTML']

    else: # UNIQ TO CONSURF

        cmd += " %s_consurf%s_pipe.pdb" %(vars['Used_PDB_Name'], vars['run_number'])
        cmd += " " + vars['isd_chimerax_file']
        cmd += " " + vars['chimerax_file']
        cmd += " " + vars['isd_header_for_chimera']
        cmd += " " + vars['header_for_chimera']
        cmd += " " + vars['rasmol_isdFILE']
        cmd += " " + vars['rasmolFILE']
        cmd += " " + vars['pdb_file_with_score_at_TempFactor']
        cmd += " " + vars['ATOMS_with_ConSurf_Scores']
        cmd += " " + vars['ATOMS_with_ConSurf_Scores_isd']
        cmd += " " + vars['scf_for_chimera']
        cmd += " " + vars['isd_scf_for_chimera']

    if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_no_pdb_no_msa": # The Server create the MSA

        cmd += " %s %s" %(vars['BLAST_out_file'], vars['FINAL_sequences_html'])

        if vars['protein_MSA_for_usr'] != "":

            cmd += " " + vars['protein_MSA_for_usr']

        else:

            cmd += " " + vars['protein_MSA']

        if os.path.exists(vars['protein_MSA_clustalw']):

            cmd += vars['protein_MSA_clustalw']

    LOG.write("zip_all_outputs cmd: %s\n" %cmd)
    exit_status = system(cmd)

    # check if zip file was created
    filesize = os.path.getsize(vars['All_Outputs_Zip'])
    if filesize > 0:

        LOG.write("zip file %s - created successfully - Size: %d , cmd exit_status: %s\n" %(vars['All_Outputs_Zip'], filesize, exit_status))

    else:

        LOG.write("!!! zip file %s - failed to create - not exist, or size is zero !!! cmd exit_status: %s \n" %(vars['All_Outputs_Zip'], exit_status))

def print_templates_projection_links():

    Num_Of_Templates = len(Templates_PDB_Descriptions)
    OUTPUT.write("Press on the plus sign ('+') near each structure id to show the relevant links for that structure.<BR>\n")
    OUTPUT.write("<div class=\"mC\">\n")
    if Num_Of_Templates > 0:

        PDB_Counter = 1
        for PDB in Templates_PDB_Descriptions.keys():

            PDB_ID = PDB[0:4]
            PDB_CHAIN = PDB[4:5]
            OUTPUT.write("<div class=\"mH\"></div>\n")
            OUTPUT.write("<span class=\"mH1\"   onclick=\"toggleMenu_and_Project_PDB(\'Prot_%s\',\'%s\',\'%s\','%s\',\'%s\',\'Prot_%s_Results\')\"><FONT FACE= \"Courier New\">+ PDB:%s CHAIN: %s </span><span class=\"mSpacer\"></span> | </FONT>\n" %(PDB, vars['run_number'], PDB_ID, PDB_CHAIN, vars['protein_seq'], PDB, PDB_ID, PDB_CHAIN))
            OUTPUT.write("\t<span class=\"mH2\" onclick=\"toggleMenu(\'Prot_" + PDB + "_Align\')\"><FONT FACE= \"Courier New\">Show Alignment</span> | </FONT>\n")
            OUTPUT.write("\t<span class=\"mH2\" onclick=\"toggleMenu(\'Prot_" + PDB + "_Details\')\"><FONT FACE= \"Courier New\">Show PDB Details</span></FONT>\n")
            OUTPUT.write("\t\t<div id=\'Prot_" + PDB + "\' class=\"mL\"><span id=\'Prot_" + PDB + "_Results\'></span></p></a></div>\n")
            OUTPUT.write("\t\t<div id=\'Prot_" + PDB + "_Align\' class=\"mL\"><FONT FACE= \"Courier New\">" + Templates_Alignments(PDB) + "</FONT></div>\n")
            OUTPUT.write("\t\t<div id=\'Prot_" + PDB + "_Details\' class=\"mL\">" + print_PDB_details(PDB) + "</div>\n")

            PDB_Counter += 1

    OUTPUT.write("</div></div>\n")

def parse_blast_alignment():

    PDB = ""
    BLAST_OUTPUT = open(vars['BLAST_PDB_out_file'], 'r+')
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

                Templates_Alignments['PDB'] = Alignment

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
    Templates_Alignments['PDB'] = Alignment

def print_PDB_details(PDB):

    PDB_ID = PDB[0:4]
    PDB_CHAIN = PDB[4:5]
    PDB_URL = GENERAL_CONSTANTS.PDB_DB
    PROTEOPEDIA_URL = GENERAL_CONSTANTS.PROTEOPEDIA
    Details = "&nbsp;&nbsp;&nbsp;<b>PDB ID:</b> %s Chain: %s<br>" %(PDB_ID, PDB_CHAIN)
    Details += "&nbsp;&nbsp;&nbsp;<b>\%Identity between query and PDB:</b> " + Templates_PDB_Descriptions[PDB]['IDENTITY'] + "<BR>"
    Details += "&nbsp;&nbsp;&nbsp;<b>\%Similarity between query and PDB:</b> " + Templates_PDB_Descriptions[PDB]['SIMILARITY'] + "<BR>"
    Details += "&nbsp;&nbsp;&nbsp;<b>PSI-blast E-Value:</b> " + Templates_PDB_Descriptions[PDB]['E_VAL'] + "<BR>"
    Details += "&nbsp;&nbsp;&nbsp;<b>Length of alignment overlap with query:</b> " + Templates_PDB_Descriptions[PDB]['Q_AlignmentLength'] + "<BR>"
    Details += "&nbsp;&nbsp;&nbsp;<b>Length of alignment overlap with PDB:</b> " + Templates_PDB_Descriptions[PDB]['S_AlignmentLength'] + "<BR>"
    Details += "&nbsp;&nbsp;&nbsp;<b>Descriptions:</b> " + Templates_PDB_Descriptions[PDB]['DESCR'] + "<BR>"
    Details += "&nbsp;&nbsp;&nbsp;Learn more about the protein on <A HREF= \'" + PROTEOPEDIA_URL + PDB_ID.lower() + "\' TARGET=Proteopedia_window>Proteopedia</A> | <A HREF= \'" + PDB_URL + PDB_ID.lower() + "\' TARGET=PDB_window>PDB</A>"

    return(Details)

def store_data():

    if stored_data_file != "":

        shutil.copyfile(stored_data_file, stored_data_file + ".AfterCGI")
        os.chmod(stored_data_file + ".AfterCGI", 0o600)

    if stored_form_data != "":

        shutil.copyfile(stored_form_data, stored_form_data + ".AfterCGI")
        os.chmod(stored_form_data + ".AfterCGI", 0o600)

    LOG.write("store_data : storing hashes to files %s and %s\n" %(stored_data_file, stored_form_data))

    # storing vars hash
    try:

        VARS_STORE = open(stored_data_file, 'w')

    except:

        exit_on_error("sys_error", "store_data : could not open %s for storing." %stored_data_file)

    json.dump(vars, VARS_STORE)
    VARS_STORE.close()

    # storing form hash
    try:

        FORM_STORE = open(stored_form_data, 'w')

    except:

       exit_on_error("sys_error", "store_data : could not open %s for storing." %stored_form_data)

    match1 = re.match(r'([^\/]+)$', form['pdb_FILE'])
    if match1:

        form['pdb_FILE'] = match1.group(1)

    match2 = re.match(r'([^\/]+)$', form['uploaded_MSA'])
    if match2:

        form['uploaded_MSA'] = match2.group(1)

    match3 = re.match(r'([^\/]+)$', form['uploaded_TREE'])
    if match3:

        form['uploaded_TREE'] = match3.group(1)

    json.dump(form, FORM_STORE)
    FORM_STORE.close()

    os.chmod(stored_data_file, 0o600)
    os.chmod(stored_form_data, 0o600)

def Find_Good_Templates(Query_Seq_File, Blast):

    Min_Overlap_Percent = CONSURF_CONSTANTS.PDB_FRAGMENT_MINIMUM_LENGTH
    Min_ID = CONSURF_CONSTANTS.PDB_FRAGMENT_MINIMUM_IDENTITY

    LOG.write("Going to find the best templates: ConSeq_gradesPE_and_Outputs.Find_Good_Templates(%s, %s, %s, %s, Templates_PDB, Templates_PDB_Descriptions);\n" %(Query_Seq_File, Blast, Min_Overlap_Percent, Min_ID))
    ans = ConSeq_gradesPE_and_Outputs.Find_Good_Templates(Query_Seq_File, Blast, Min_Overlap_Percent, Min_ID, Templates_PDB, Templates_PDB_Descriptions)
    LOG.write("%s: %s\n" %(ans[0], nas[1]))
    if ans[0] == "OK":

        retrun(len(Templates_PDB_Descriptions) + 1)

    else:

        exit_on_error('sys_error', "ConSeq_gradesPE_and_Outputs.Find_Good_Templates FAILED: " + str(ans))

def send_finish_email_to_user():

    email_subject = ""
    HttpPath = vars['run_url'] + dir + vars['output_page']
    if form['job_title'] != "":

        email_subject = "'The results of your ConSurf run titled " + form['job_title'] + " are ready'"

    elif vars['Used_PDB_Name'] != "":

        email_subject = "'Your ConSurf results for PDB: %s chain: %s are ready'" %(form['pdb_ID'], form['chain'])

    else:

        email_subject = "'Your ConSurf results for the seq " + vars['msa_SEQNAME'] + " are ready'"

    email_message = "'Hello,\\n\\nThe results for your ConSurf run are ready at:\\n" + HttpPath + "\\n\\nRunning Parameters:\\n"
    if vars['Used_PDB_Name'] != "": #ConSurf Mode

        email_message += "PDB: %s\\nCHAIN: %s\\n" %(vars['Used_PDB_Name'], form['chain'])

    if 'user_msa_fasta' in vars:

        email_message += "Alignment: %s\\n" %vars['user_msa_fasta']

    else:

        email_message += "Alignment: Created by ConSurf using %s and %s\\n" %(form['database'], form['MSAprogram'])

    email_message += "\\nPlease note: the results will be kept on the server for three months.'"
    msg = "%s -f \'%s\' -t \'%s\' -u %s -xu %s -xp %s -s %s -m %s" %(vars['send_email_cmd'], GENERAL_CONSTANTS.ADMIN_EMAIL, form['userEmail'], email_subject, vars['userName'], vars['userPass'], vars['smtp_server'], email_message)
    LOG.write("MESSAGE:%s\nCOMMAND:%s\n" %(email_message, msg))
    os.mkdir(vars['send_email_dir'])
    email_system_return = msg
    match = re.match(r'successfully', email_system_return)
    if not match:

        LOG.write("send_mail: The message was not sent successfully. system returned: %s\n" %email_system_return)

def update_output_that_run_finished(OutHtmlFile):

    # write to the output that the job has finished and stop reload
    OUTPUT = open(OutHtmlFile, 'r')
    lines = OUTPUT.readlines()
    OUTPUT.close()

    OUTPUT = open(OutHtmlFile, 'w')
    for line in lines:

        match1 = re.match(r'.*consurf_output_text_running.tpl.*', line)
        if match1:

            OUTPUT.write("include (\"" + server_html_dir + "/templates/consurf_output_text_finished.tpl\");\n")

        else:

            match2 = re.match(r'.*output_header.tpl.*', line)
            if match2:

                OUTPUT.write("include (\"" + server_html_dir + "/templates/output_header_no_refresh.tpl\");\n")

            else:

                OUTPUT.write(line)

    OUTPUT.close()

def create_pymol_instructions(insufficient_data, pymol_instructions, Used_PDB_Name, chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, identical_chains, working_dir, chimerax_script_for_figure, chimerax_script_for_figure_isd, pdb_file_name_NoPath):

    LOG.write("Calling\n")
    if insufficient_data == "yes":

        ans = cp_rasmol_gradesPE_and_pipe.create_pymol_page(pymol_instructions, Used_PDB_Name, chain, vars['Confidence_link'], vars['pymol_color_script'], ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, identical_chains, vars['run_number'], working_dir, chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, chimerax_script_for_figure, chimerax_script_for_figure_isd, pdb_file_name_NoPath, vars['pymol_color_script_CBS'])
        if ans != "OK":

            exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_pymol_page FAILED: " + ans)

    else:

        ans = cp_rasmol_gradesPE_and_pipe.create_pymol_page(pymol_instructions, Used_PDB_Name, chain, vars['Confidence_link'], vars['pymol_color_script'], ATOMS_with_ConSurf_Scores, "", identical_chains, vars['run_number'], working_dir, chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, chimerax_script_for_figure, chimerax_script_for_figure_isd, pdb_file_name_NoPath, vars['pymol_color_script_CBS'])
        if ans != "OK":

            exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_pymol_page FAILED: " + ans)

def blast_vs_PDB(query, Blast_PDB_Out_File):

    LOG.write("blast_vs_PDB(%s, %s)\n" %(query, Blast_PDB_Out_File))
    cmd = [pgp, "-i", query, "-e", str(form['ESCORE']), "-d", vars['pdb_db'], "-j", str(form['iterations']), "-v", str(vars['max_PDB_homologues_to_display']), "-b", str(vars['max_PDB_homologues_to_display']), "-o", Blast_PDB_Out_File, "-F", "F", "-T", "F"]
    LOG.write("run_blast : running: %s\n" %cmd)
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    p.communicate()

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

def create_MSA():

    if form['MSAprogram'] == "CLUSTALW":

        cmd = [GENERAL_CONSTANTS.CLUSTALW, "-infile=" + vars['FINAL_sequences'], "-outfile=" + vars['msa_clustal']]
        p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
        p.communicate()
        ans = MSA_parser.convert_msa_format(vars['msa_clustal'], "clustal", vars['msa_fasta'], "fasta")
        if ans != "OK":

            exit_on_error('sys_error', "create_MSA: MSA_parser.convert_msa_format - %s\n" %ans)

    elif form['MSAprogram'] == "T_COFFEE":

        print_message_to_output("<font color='red'><b>Warning:</b></font> T-COFFEE (EXPRESSO) is accurate but slow MSA program, please be patient.")
        cmd = "%s %s -mode expresso -blast_server=LOCAL -pdb_db %s -cache=ignore -output=fasta 1> t_coffe.std 2>&1" %(GENERAL_CONSTANTS.T_COFFEE, vars['FINAL_sequences'], GENERAL_CONSTANTS.PDBAA_NCBI)
        LOG.write("create_MSA:%s\n" %cmd)
        cmd # should execute command
        if os.path.exists("query_final_homolougs.fasta_aln"):

            shutil.copyfile("query_final_homolougs.fasta_aln", vars['protein_MSA'])

        shutil.copyfile(vars['protein_MSA'], vars['protein_MSA_FASTA'])
        ans = MSA_parser.convert_msa_format(vars['protein_MSA'], "fasta", vars['protein_MSA_clustalw'], "clustalw")
        if ans[0] != "OK":

            exit_on_error('sys_error', "create_MSA: MSA_parser.convert_msa_format - %s\n" %str(ans))

        vars['protein_MSA_for_usr'] = vars['protein_MSA']
        vars['protein_MSA'] = vars['protein_MSA_clustalw']
        vars['msa_format'] = "clustalw"

    elif form['MSAprogram'] == "MAFFT":

        cmd = "%s --quiet %s >  %s" %(GENERAL_CONSTANTS.MAFFT_LINSI_GUIDANCE, vars['FINAL_sequences'], vars['protein_MSA'])
        cmd # should execute command
        shutil.copyfile(vars['protein_MSA'], vars['protein_MSA_FASTA'])
        ans = MSA_parser.convert_msa_format(vars['protein_MSA'], "fasta", vars['protein_MSA_clustalw'], "clustalw")
        if ans[0] != "OK":

            exit_on_error('sys_error', "create_MSA: MSA_parser.convert_msa_format - %s\n" %str(ans))

        vars['protein_MSA_for_usr'] = vars['protein_MSA']
        vars['protein_MSA'] = vars['protein_MSA_clustalw']
        vars['msa_format'] = "clustalw"

    elif form['MSAprogram'] == "PRANK":

        cmd = "%s -d=%s -o=%s -F" %(GENERAL_CONSTANTS.PRANK, vars['FINAL_sequences'], vars['protein_MSA'])
        print_message_to_output("<font color='red'><b>Warning:</b></font> PRANK is accurate but slow MSA program, please be patient.")
        LOG.write("create_MSA : run %s\n" %cmd)
        cmd # should execute command
        if os.path.exists(vars['protein_MSA'] + ".2.fas"):

            shutil.copyfile(vars['protein_MSA'] + ".2.fas", vars['protein_MSA'])

        elif os.path.exists(vars['protein_MSA'] + ".1.fas"):

            shutil.copyfile(vars['protein_MSA'] + ".1.fas", vars['protein_MSA'])

        elif os.path.exists(vars['protein_MSA'] + ".best.fas"):

            shutil.copyfile(vars['protein_MSA'] + ".best.fas", vars['protein_MSA'])

        shutil.copyfile(vars['protein_MSA'], vars['protein_MSA_FASTA'])
        vars['protein_MSA_for_usr'] = vars['protein_MSA']
        ans = MSA_parser.convert_msa_format(vars['protein_MSA'], "fasta", vars['protein_MSA_clustalw'], "clustalw")
        if ans[0] != "OK":

            exit_on_error('sys_error', "create_MSA: MSA_parser.convert_msa_format - %s\n" %str(ans))

        vars['protein_MSA'] = vars['protein_MSA_clustalw']
        vars['msa_format'] = "clustalw"
        vars['protein_MSA_FASTA'] = vars['protein_MSA']

    else:

        cmd = "%s -in %s -out %s -clwstrict -quiet" %(GENERAL_CONSTANTS.MUSCLE, vars['FINAL_sequences'], vars['protein_MSA'])
        vars['msa_format'] = "clustalw"
        cmd # should execute command
        ans = MSA_parser.convert_msa_format(vars['protein_MSA'], "clustalw", vars['protein_MSA_FASTA'], "fasta")
        if ans[0] != "OK":

            exit_on_error('sys_error', "create_MSA: MSA_parser.convert_msa_format - %s\n" %str(ans))
    """
    LOG.write("create_MSA : run %s\n" %cmd)
    if not os.path.exists(vars['working_dir'] + dir + vars['protein_MSA']) or os.path.getsize(vars['working_dir'] + dir + vars['protein_MSA']) == 0:

        exit_on_error('sys_error', "create_MSA : the file " + vars['working_dir'] + dir + vars['protein_MSA'] + " was not created or of size zero" )
    """
def open_log_file():

    global LOG
    if not os.path.exists(vars['run_log_Q']) or os.path.getsize(vars['run_log_Q']) == 0:

        try:

            LOG = open(vars['run_log_Q'], 'w')

        except:

            exit_on_error('sys_error', "Cannot open the log file '" + vars['run_log_Q'] + "' for writing.")

        LOG.write("--------- ConSurf Log " + vars['run_number'] + "_Q.log -------------\n")
        LOG.write("Begin Time: %d\n" %time.time())

    else:

        try:

            LOG = open(vars['run_log_Q'], 'a')

        except:

            exit_on_error('sys_error', "Cannot open the log file '" + vars['run_log_Q'] + "' for writing.")

def exit_on_error(which_error, error_msg):

    error_definition = "<font size=+1 color='red'>ERROR! ConSurf session has been terminated:</font><br />\n"
    syserror = "<font size=+1 color='red'>A SYSTEM ERROR OCCURRED!</font><br />Please try to run ConSurf again in a few minutes.<br />We apologize for the inconvenience.<br />\n"

    if which_error == "user_error":

        LOG.write("\nEXIT on error:\n%s\n" %error_msg)
        OUTPUT.write(error_definition + error_msg)
        USER_ERROR = open("USER_ERROR", 'w')
        USER_ERROR.write(error_msg + "\n")
        USER_ERROR.close()
        # print error_msg to the screen

    elif which_error == "sys_error":

        send_administrator_mail_on_error(error_msg)
        LOG.write("\n%s\n" %error_msg)
        OUTPUT.write(syserror)
        #print error_msg to the log file
        SYS_ERROR = open("SYS_ERROR", 'w')
        SYS_ERROR.write(error_msg + "\n")
        SYS_ERROR.close()

    OUTPUT.write(""""<br>
    <?php
        include(\"%stemplates/footer.tpl\");
    ?>""" %server_html_dir)
    OUTPUT.close()
    # finish the output page
    time.sleep(10)
    OUTPUT = open(vars['output_page'], 'r')
    lines = OUTPUT.readlines()
    OUTPUT.close()
    # remove the refresh commands from the output page
    OUTPUT = open(vars['output_page'], 'w')
    for line in lines:

        match = re.match(r'consurf_output_text_running.tpl', line)
        if match:

            OUTPUT.write("include (\"" + server_html_dir + "templates/consurf_output_text_failed.tpl\");\n")

        elif line == "include (\"" + server_html_dir + "templates/output_header.tpl\");\n":

            OUTPUT.write("include (\"" + server_html_dir + "templates/output_header_no_refresh.tpl\");\n")

        else:

            OUTPUT.write(line)

    OUTPUT.close()
    if form['send_user_mail']:

        send_mail_on_error()

    LOG.write("\nExit Time: " + BIOSEQUENCE_FUNCTIONS.printTime + "\n")
    LOG.close()
    os.chmod(vars['working_dir'], 0o755)
    exit(1)

def run_rate4site(query_name):

    algorithm = ""
    tree_file_r4s = ""
    msa = ""
    did_r4s_fail = ""
    MatrixHash = {'JTT' : '-Mj', 'MTREV' : '-Mr', 'CPREV' : '-Mc', 'WAG' : '-Mw', 'DAYHOFF' : '-Md', 'T92' : '-Mt', 'HKY' : '-Mh', 'GTR' : '-Mg', 'JC_Nuc' : '-Mn', 'JC_AA' : '-Ma', 'LG' : '-Ml'}

    if vars['running_mode'] == "_mode_pdb_msa_tree" or vars['running_mode'] == "_mode_msa_tree":

        tree_file_r4s = "-t " + vars['user_tree_file_name']
    """
    if vars['running_mode'] == "_mode_pdb_no_msa" or vars['running_mode'] == "_mode_no_pdb_no_msa":


        msa = vars['protein_MSA']
        if form['MSAprogram'] == 'CLUSTALW' or form['MSAprogram'] == 'MUSCLE' or form['MSAprogram'] == 'MAFFT' or form['MSAprogram'] == 'PRANK' or form['MSAprogram'] == 'T_COFFEE':

            vars['rate4site_msa_format'] = "clustal"

        else:

            vars['rate4site_msa_format'] = "fasta"

    else:

        LOG.write("run_rate4site : Please note: MSA for rate4site run is '%s' (and not original user file : '%s')\n" %(vars['user_msa_fasta'], vars['user_msa_file_name']))
        msa = vars['user_msa_fasta']
        vars['rate4site_msa_format'] = "fasta"
    """
    r4s_comm = []
    if form['algorithm'] == "Bayes":

        algorithm = "-ib" # Save the algorithm, it maybe used later
        r4s_comm += [rate4s, "-ib", "-a", query_name, "-s", vars['msa_fasta'], "-zn", MatrixHash[(form['matrix']).upper()], tree_file_r4s, "-bn", "-l", vars['r4s_log'], "-o", vars['r4s_out'], "-n", "32", "-v", "9"]

    else:

        algorithm = "-im" # Save the algorithm, it maybe used later
        r4s_comm += [rate4s_ML, "-im", "-a", query_name, "-s", vars['msa_fasta'], "-zn", MatrixHash[(form['matrix']).upper()], tree_file_r4s, "-bn", "-l", vars['r4s_log'], "-o", vars['r4s_out'], "-v", "9"]

    LOG.write("run_rate4site : running command: %s\n" %r4s_comm)
    p = subprocess.Popen(r4s_comm, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    out, err = p.communicate()

    did_r4s_fail = check_if_rate4site_failed(vars['r4s_out'], vars['r4s_log'])
    # if the run failed - we rerun using the slow verion
    if did_r4s_fail == "yes":

        LOG.write("run_rate4site : The run of rate4site failed. Sending warning message to output.\nThe same run will be done using the SLOW version of rate4site.\n")
        print_message_to_output("<font color='red'><b>Warning:</b></font> The given MSA is very large, therefore it will take longer for ConSurf calculation to finish. The results will be sent to the e-mail address provided.<br>The calculation continues nevertheless.")
        form['send_user_mail'] = "yes"
        r4s_comm = [rate4s_slow, algorithm, "-a", query_name, "-s", vars['msa_fasta'], "-zn", MatrixHash[(form['matrix']).upper()], tree_file_r4s, "-bn", "-l", vars['r4s_slow_log'], "-o", vars['r4s_out'], "-n", "32", "-v", "9"]
        LOG.write("run_rate4site : running command: %s\n" %str(r4s_comm))
        p = subprocess.Popen(r4s_comm, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
        out, err = p.communicate()

        did_r4s_fail = check_if_rate4site_failed(vars['r4s_out'], vars['r4s_slow_log'])
        if did_r4s_fail == "yes":

            err = "The run %s for %s failed.\nThe MSA %s was too large.\n" %(r4s_comm, vars['run_number'], vars['msa_fasta'])
            exit_on_error('user_error', "The calculation could not be completed due to memory problem, since the MSA is too large. Please run ConSurf again with fewer sequences.")

def check_if_rate4site_failed(res_flag, r4s_log):

    # There are some tests to see if rate4site failed.
    # Since I can't trust only one of them, I do all of them. If onw of them is tested to be true - than a flag will get TRUE value
    # 1. the .res file might be empty.
    # 2. if the run failed, it might be written to the log file of r4s.
    # 3. in a normal the r4s.log file there will lines that describe the grades. if it fail - we won't see them
    # In one of these cases we try to run the slower version of rate4site.
    # We output this as a message to the user.
    ret = "no"
    did_r4s_failed = rate4site_routines.check_if_rate4site_failed(res_flag, r4s_log)
    if did_r4s_failed[0] == "yes":

        ret = "yes"
        LOG.write("check_if_rate4site_failed : " + did_r4s_failed[1])
        vars['r4s_process_id'] = did_r4s_failed[2]
        remove_core()
        # if there was an error in user input which rate4site reported: we output a message to the user and exit
        if 3 in did_r4s_failed and did_r4s_failed[3] != "":

            exit_on_error('user_error', did_r4s_failed[3])

    return ret

def print_message_to_output(msg):

    txt = "\n<ul><li>%s</li></ul>\n" %msg
    OUTPUT.write(txt)
    return txt

def remove_core():

    if os.path.exists(vars['working_dir'] + "/core." + vars['r4s_process_id']):

        LOG.write("remove core file : core.%s\n" %vars['r4s_process_id'])
        os.unlink(vars['working_dir'] + "/core." + vars['r4s_process_id'])

def assign_colors_according_to_r4s_layers(ref_to_gradesPE, r4s_out):

    LOG.write("assign_colors_according_to_r4s_layers : %s\n" %r4s_out)
    ans = cp_rasmol_gradesPE_and_pipe.assign_colors_according_to_r4s_layers(ref_to_gradesPE, r4s_out)
    if ans != "OK":

        exit_on_error('sys_error', ans)

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

def print_residue_precentage():

    # print a file that details percentage of each AA in the MSA
    LOG.write("cp_rasmol_gradesPE_and_pipe.print_precentage(%s, %s, %s, %s)\n" %(str(residue_freq), str(position_totalAA), vars['Msa_percentageFILE'], str(gradesPE_Output)))
    ans = cp_rasmol_gradesPE_and_pipe.print_precentage(residue_freq, position_totalAA, vars['Msa_percentageFILE'], gradesPE_Output)
    if ans[0] != "OK":

        exit_on_error('sys_error',str(ans))

    elif not os.path.exists(vars['Msa_percentageFILE']) or os.path.getsize(vars['Msa_percentageFILE']) == 0:

        exit_on_error('sys_error', "The output " + vars['Msa_percentageFILE'] + " was not found or empty")

def print_nucleotide_precentage():

    # print a file that details percentage of each nucleotide in the MSA
    LOG.write("cp_rasmol_gradesPE_and_pipe.print_precentage_nuc(%s, %s, %s)\n" %(str(residue_freq), str(position_totalAA), vars['Msa_percentageFILE']))
    ans = cp_rasmol_gradesPE_and_pipe.print_precentage_nuc(residue_freq, position_totalAA, vars['Msa_percentageFILE'])
    if ans[0] != "OK":

        exit_on_error('sys_error',str(ans))

    elif not os.path.exists(vars['Msa_percentageFILE']) or os.path.getsize(vars['Msa_percentageFILE']) == 0:

        exit_on_error('sys_error', "The output " + vars['Msa_percentageFILE'] + " was not found or empty")

def create_atom_position_file(user_chain, pdb_file, atom_positionFILE, Type, isPisa):

    chain = " "
    if not re.match(r'none', user_chain, re.IGNORECASE):

        chain = user_chain

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

def fill_r4s2pdb(ref_r4s2pdb):

    chain = trim(form['chain'])
    LENGTH_OF_ATOM = 0
    match = re.match(r'none', form['chain'], re.IGNORECASE)
    if match:

        chain = " "

    try:

        POSITIONS = open(vars['atom_positionFILE'], 'r')

    except:

        exit_on_error('sys_error', "fill_r4s2pdb: Can't open " + vars['atom_positionFILE'])

    line = POSITIONS.readline()
    while line != "":

        line = line.rstrip()
        frags = line.split()
        Residue_Name = trim(frags[0])
        FAS_Pos = trim(frags[1])
        ATOM_pos = trim(frags[2])
        ref_r4s2pdb[FAS_Pos] ="%d%s:%s" %(Residue_Name, ATOM_pos ,chain)
        LENGTH_OF_ATOM += 1
        line = POSITIONS.readline()

    POSITIONS.close()

def match_pdb_to_seq(ref_r4s2pdb, user_chain, pdbseq, query_seq, atom_positionFILE, Type):

    chain = " "
    match = re.match(r'none', form['chain'], re.IGNORECASE)
    if not match:

        chain = user_chain

    LOG.write("match_pdb_to_seq : calling cp_rasmol_gradesPE_and_pipe.match_seqres_pdb(%s ,%s, %s, %s, %s, %s)\n" %(pdbseq, query_seq, atom_positionFILE, chain, ref_r4s2pdb, Type))
    ans = cp_rasmol_gradesPE_and_pipe.match_seqres_pdb(pdbseq, query_seq, atom_positionFILE, chain, ref_r4s2pdb, Type)
    if not ans[0] == "OK":

        exit_on_error('sys_error', "match_pdb_to_seq : rasmol_gradesPE_and_pipe.%s" %ans[0])

    elif len(ref_r4s2pdb) < 1:

        LOG.write("match_pdb_to_seq : Total residues in the msa sequence: %s. Total residues in the ATOM : %s\n" %(ans[1], ans[2]))

    return(ans[1], ans[2])

def replace_TmpFactor_Consurf_Scores(user_chain, pdb_file, gradesPE, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd):

    # This Will create a File containing the ATOMS records with the ConSurf grades instead of the TempFactor column

    chain = " "
    match = re.match(r'none', form['chain'], re.IGNORECASE)
    if not match:

        chain = user_chain

    LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.ReplaceTempFactConSurfScores(%s, %s, %s, %s, %s);\n" %(chain, pdb_file, gradesPE, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd))
    ans = cp_rasmol_gradesPE_and_pipe.ReplaceTempFactConSurfScore(chain, pdb_file, gradesPE, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd)
    if ans != "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.ReplaceTempFactConSurf FAILED: " + ans)

def replace_TmpFactor_Rate4Site_Scores(user_chain, pdb_file, gradesPE, pdb_file_with_score_at_TempFactor):

    # This will create a PDB file that contains the Rate4Site scores instead of the TempFactor Column

    chain = " "
    match = re.match(r'none', form['chain'], re.IGNORECASE)
    if not match:

        chain = user_chain

    Rate4Site_Grades = {}
    LOG.write("Calling:cp_rasmol_gradesPE_and_pipe.read_Rate4Site_gradesPE(%s, %s)\n" %(gradesPE, str(Rate4Site_Grades)))
    ans = cp_rasmol_gradesPE_and_pipe.read_Rate4Site_gradesPE(gradesPE, Rate4Site_Grades)
    if ans != "OK":

        exit_on_error("sys_error", ans)

    LOG.write("Calling:cp_rasmol_gradesPE_and_pipe.replace_tempFactor(%s, %s, %s, %s)\n" %(pdb_file, chain, str(Rate4Site_Grades), pdb_file_with_score_at_TempFactor))
    ans = cp_rasmol_gradesPE_and_pipe.replace_tempFactor(pdb_file, chain, Rate4Site_Grades, pdb_file_with_score_at_TempFactor)
    if ans != "OK":

        exit_on_error("sys_error", ans)

def create_chimera(aln, working_dir, r4s_out, scf_for_chimera, header_for_chimera, isd_scf_for_chimera, isd_header_for_chimera, ATOMS_with_ConSurf_Scores, chimerax_file, ATOMS_with_ConSurf_Scores_isd, isd_chimerax_file, chimerax_script_for_figure, chimerax_script_for_figure_isd, chimera_instructions, identical_chains, chain, pdb_file_name_NoPath):

    # Chimera Output includes several files:
    # a. header file *.hdr
    # b. scf file *.scf
    # c. script to show the colored MSA, Tree, and colored 3D structure. (.chimerax)
    # d. the Script for Chimera Image.
    # e. The Html with the istructions how to create Chimera High resolution Image
    #
    # A.+B. creating the *.hdr and *.scf files
    # creating view page for the chimera alingment requires the query name for the input sequence in the MSA. In case MSA was uploaded -
    # this name might not be exact, so the option of viewing the alignment only applies for cases where the user did not supply MSA

    insufficient_data = "no"

    """
    # the script must get MSA in ClustalW format
    aln = protein_MSA
    if msa_format != "clustalw":

        LOG.write("create_chimera: Calls 'MSA_parser.convert_msa_format(%s, %s, %s, %s);'\n" %(protein_MSA, msa_format, vars['protein_MSA_clustalw'], "clustalw"))
        ans = MSA_parser.convert_msa_format(protein_MSA, msa_format, vars['protein_MSA_clustalw'], "clustal")
        if ans != "OK":

            exit_on_error('user_error', "create_chimera: MSA_parser.convert_msa_format - %s\n" %ans)

        aln = vars['protein_MSA_clustalw']

    elif not os.path.exists(vars['protein_MSA_clustalw']):

        # we are already in MSA created by clustalw, just the name is incorrect. for consistency we copy it to be with the desired name
        try:

            shutil.copyfile(vars['protein_MSA'], vars['protein_MSA_clustalw'])

        except:

            exit_on_error('sys_error', "create_chimera: Copy failed.\n")

        aln = vars['protein_MSA_clustalw']


    LOG.write("MSA: %s\tFORMAT:%s\n" %(aln, msa_format))
    """

    LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.color_with_chimera(%s, %s, %s, %s, %s, %s, %s)\n" %(vars['msa_SEQNAME'], aln, r4s_out, scf_for_chimera, header_for_chimera, isd_scf_for_chimera, isd_header_for_chimera))
    ans = cp_rasmol_gradesPE_and_pipe.color_with_chimera(vars['msa_SEQNAME'], aln, r4s_out, scf_for_chimera, header_for_chimera, isd_scf_for_chimera, isd_header_for_chimera)
    if ans[0] != "OK":

        exit_on_error('sys_error', "create_chimera: %s\n" %str(ans))

    else:

        insufficient_data = ans[1]

    # C. creating the script that shows the MSA, Tree and colored 3D structure
    LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.create_chimera_script(%s, %s, %s, %s, %s, %s, %s, %s)\n" %(ATOMS_with_ConSurf_Scores, vars['run_url'], working_dir, chimerax_file, aln, vars['tree_file'], scf_for_chimera, header_for_chimera))
    cp_rasmol_gradesPE_and_pipe.create_chimera_script(ATOMS_with_ConSurf_Scores, vars['run_url'], chimerax_file, aln, vars['tree_file'], scf_for_chimera, header_for_chimera)

    # create also ignoring insufficient data file, only in case the selected algorithm was bayes
    if insufficient_data == "yes" and form['algorithm'] == "Bayes":

        LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.create_chimera_script(%s, %s, %s, %s, %s, %s, %s, %s)\n" %(ATOMS_with_ConSurf_Scores_isd, vars['run_url'], working_dir, isd_chimerax_file, aln, vars['tree_file'], isd_scf_for_chimera, isd_header_for_chimera))
        cp_rasmol_gradesPE_and_pipe.create_chimera_script(ATOMS_with_ConSurf_Scores_isd, vars['run_url'], working_dir, isd_chimerax_file, aln, vars['tree_file'], isd_scf_for_chimera, isd_header_for_chimera)

    # D. The Script For Chimera Image
    LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.create_chimera_image_script(%s, %s, %s)\n" %(chimerax_script_for_figure, ATOMS_with_ConSurf_Scores, vars['run_url']))
    cp_rasmol_gradesPE_and_pipe.create_chimera_image_script(chimerax_script_for_figure, ATOMS_with_ConSurf_Scores, vars['run_url'])
    if insufficient_data == "yes":

        LOG.wrrite("Calling: cp_rasmol_gradesPE_and_pipe.create_chimera_image_script(%s, %s, %s)\n" %(chimerax_script_for_figure_isd, ATOMS_with_ConSurf_Scores_isd, vars['run_url']))
        cp_rasmol_gradesPE_and_pipe.create_chimera_image_script(chimerax_script_for_figure_isd, ATOMS_with_ConSurf_Scores_isd, vars['run_url'])

    # E. Create Chimera HTML Page
    if insufficient_data == "yes":

        LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.create_chimera_page(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n" %(chimera_instructions, vars['Confidence_link'], vars['chimera_color_script'], chimerax_script_for_figure, ATOMS_with_ConSurf_Scores, chimerax_script_for_figure_isd, ATOMS_with_ConSurf_Scores_isd, vars['run_url_old'], identical_chains, vars['run_number'], working_dir, chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, chimerax_script_for_figure, chimerax_script_for_figure_isd, pdb_file_name_NoPath, vars['chimera_color_script_CBS']))
        ans = cp_rasmol_gradesPE_and_pipe.create_chimera_page(chimera_instructions, vars['Confidence_link'], vars['chimera_color_script'], chimerax_script_for_figure, ATOMS_with_ConSurf_Scores, chimerax_script_for_figure_isd, ATOMS_with_ConSurf_Scores_isd, vars['run_url_old'], identical_chains, vars['run_number'], working_dir, chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, chimerax_script_for_figure, chimerax_script_for_figure_isd, pdb_file_name_NoPath, vars['chimera_color_script_CBS'])
        if ans != "OK":

            exit_on_error('sys_error',"cp_rasmol_gradesPE_and_pipe.create_chimera_page FAILED: %s" %ans)

    else:

        LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.create_chimera_page(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n" %(chimera_instructions, vars['Confidence_link'], vars['chimera_color_script'], chimerax_script_for_figure, ATOMS_with_ConSurf_Scores, "", "", "", identical_chains, vars['run_number'], working_dir, chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, chimerax_script_for_figure, chimerax_script_for_figure_isd, pdb_file_name_NoPath, vars['chimera_color_script_CBS']))
        ans = cp_rasmol_gradesPE_and_pipe.create_chimera_page(chimera_instructions, vars['Confidence_link'], vars['chimera_color_script'], chimerax_script_for_figure, ATOMS_with_ConSurf_Scores, "", "", "", identical_chains, vars['run_number'], working_dir, chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, chimerax_script_for_figure, chimerax_script_for_figure_isd, pdb_file_name_NoPath, vars['chimera_color_script_CBS'])
        if ans != "OK":

            exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_chimera_page FAILED: %s" %ans)

    return(insufficient_data)

def send_administrator_mail_on_error(message):

    email_subject = "'System ERROR has occurred on ConSurf:" + vars['run_url'] + "'"
    email_message = "'Hello,\\n\\nUnfortunately a system System ERROR has occurred on ConSurf: " + vars['run_url'] + ".\\nERROR: " + message + ".'"
    Admin_Email = GENERAL_CONSTANTS.ADMIN_EMAIL
    msg = "%s -f 'bioSequence@tauex.tau.ac.il' -t 'bioSequence\@tauex.tau.ac.il' -u  %s -xu %s -xp %s -s %s -m %s" %(vars['send_email_cmd'], email_subject, vars['userName'], vars['userPass'], vars['smtp_server'], email_message)
    os.mkdir(vars['send_email_dir'])
    email_system_return = msg

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

def create_gradesPE_ConSeq(ref_Solv_Acc_Pred):

    # Create ConSeq gradesPE file

    ans = cp_rasmol_gradesPE_and_pipe.create_gradesPE_ConSeq(gradesPE_Output, residue_freq, ref_Solv_Acc_Pred, vars['gradesPE'])
    if not ans == "OK":

        exit_on_error('sys_error', "create_gradesPE_ConSeq : rasmol_gradesPE_and_pipe." + ans)

    elif not os.path.exists(vars['gradesPE']) or os.path.getsize(vars['gradesPE']) == 0:

        exit_on_error('sys_error', "create_gradesPE_ConSeq : the file " + vars['gradesPE'] + " was not found or empty")



def create_gradesPE_ConSeq_Nuc(ref_Solv_Acc_Pred):

    # Create ConSeq gradesPE file for Nucleotides

    ans = cp_rasmol_gradesPE_and_pipe.create_gradesPE_ConSeq_Nuc(gradesPE_Output, residue_freq, vars['gradesPE'])
    if not ans[0] == "OK":

        exit_on_error('sys_error', "create_gradesPE_ConSeq_Nuc : rasmol_gradesPE_and_pipe_Nuc." + str(ans[0]))

    elif not os.path.exists(vars['gradesPE']) or os.path.getsize(vars['gradesPE']) == 0:

        exit_on_error('sys_error', "create_gradesPE_ConSeq_Nuc : the file " + vars['gradesPE'] + " was not found or empty")

def create_rasmol(user_chain, rasmolFILE, rasmol_isdFILE, rasmolFILE_CBS, rasmol_isdFILE_CBS, no_isd_residue_color_ArrRef, isd_residue_color_ArrRef):

    # print 2 rasmol files, one showing insufficient data, one hiding it. for color blind scale and legacy scale

    chain = " "
    match = re.match(r'none', form['chain'], re.IGNORECASE)
    if not match:

        chain = user_chain

    LOG.write("Calling cp_rasmol_gradesPE_and_pipe.print_rasmol for files %s and %s\n" %(rasmolFILE, rasmol_isdFILE))
    LOG.write("cp_rasmol_gradesPE_and_pipe.print_rasmol(%s, %s, %s, %s, %s, %s)\n" %(rasmolFILE, "no", no_isd_residue_color_ArrRef, chain, "no", "legacy"))
    ans = cp_rasmol_gradesPE_and_pipe.print_rasmol(rasmolFILE, "no", no_isd_residue_color_ArrRef, chain, "no", "legacy") # Without isd residue Color
    if not ans[0] == "OK":

        exit_on_error('sys_error', "create_rasmol : cp_rasmol_gradesPE_and_pipe." + str(ans[0]))

    LOG.write("cp_rasmol_gradesPE_and_pipe.print_rasmol(%s, %s, %s, %s, %s, %s)\n" %(rasmol_isdFILE, "yes", no_isd_residue_color_ArrRef, chain, "no", "legacy"))
    ans = cp_rasmol_gradesPE_and_pipe.print_rasmol(rasmol_isdFILE, "yes", no_isd_residue_color_ArrRef, chain, "no", "legacy") # With isd Residue Color
    if not ans[0] == "OK":

        exit_on_error('sys_error', "create_rasmol : cp_rasmol_gradesPE_and_pipe." + str(ans[0]))

    if not os.path.exists(rasmolFILE) or os.path.getsize(rasmolFILE) == 0 or not os.path.exists(rasmol_isdFILE) or os.path.getsize(rasmol_isdFILE) == 0:

        exit_on_error('sys_error', "create_rasmol : Did not create one of rasmol outputs")

    # color blind friendly
    LOG.write("Calling cp_rasmol_gradesPE_and_pipe.print_rasmol for files %s and %s\n" %(rasmolFILE_CBS, rasmol_isdFILE_CBS))
    LOG.write("cp_rasmol_gradesPE_and_pipe.print_rasmol(%s, %s, %s, %s, %s, %s)\n" %(rasmolFILE_CBS, "no", no_isd_residue_color_ArrRef, chain, "no", "CBS"))
    ans = cp_rasmol_gradesPE_and_pipe.print_rasmol(rasmolFILE_CBS, "no", no_isd_residue_color_ArrRef, chain, "no", "CBS") # Without isd residue Color
    if not ans[0] == "OK":

        exit_on_error('sys_error', "create_rasmol : cp_rasmol_gradesPE_and_pipe." + str(ans[0]))

    LOG.write("cp_rasmol_gradesPE_and_pipe.print_rasmol(%s, %s, %s, %s, %s, %s)\n" %(rasmol_isdFILE_CBS, "yes", no_isd_residue_color_ArrRef, chain, "no", "CBS"))
    ans = cp_rasmol_gradesPE_and_pipe.print_rasmol(rasmol_isdFILE_CBS, "yes", no_isd_residue_color_ArrRef, chain, "no", "CBS") # With isd Residue Color
    if not ans[0] == "OK":

        exit_on_error('sys_error', "create_rasmol : cp_rasmol_gradesPE_and_pipe." + str(ans[0]))

    if not os.path.exists(rasmolFILE_CBS) or os.path.getsize(rasmolFILE_CBS) == 0 or not os.path.exists(rasmol_isdFILE_CBS) or os.path.getsize(rasmol_isdFILE_CBS) == 0:

        exit_on_error('sys_error', "create_rasmol : Did not create one of rasmol outputs")

def create_chimera_align_tree():

    # Chimera Output includes several files:
    # a. header file *.hdr
    # b. scf file *.scf
    # c. script to show the colored MSA and Tree. (.chimerax)
    #

    # A.+B. creating the *.hdr and *.scf files
    # creating view page for the chimera alingment requires the query name for the input sequence in the MSA. In case MSA was uploaded - this name might not be exact, so the option of viewing the alignment only applies for cases where the user did not supply MSA

    # the script must get MSA in ClustalW format
    aln = vars['msa_clustal']
    """
    if vars['msa_format'] != "clustalw":

        LOG.write("create_chimera: Calls 'MSA_parser.convert_msa_format(%s, %s, %s, %s)'\n" %(vars['working_dir'] + dir + vars['protein_MSA'], vars['msa_format'], vars['working_dir'] + dir + vars['protein_MSA_clustalw'], "clustalw"))
        ans = MSA_parser.convert_msa_format(vars['working_dir'] + dir + vars['protein_MSA'], vars['msa_format'], vars['working_dir'] + dir + vars['protein_MSA_clustalw'], "clustalw")
        if ans[0] == "OK":

            exit_on_error('user_error', "create_chimera: MSA_parser.convert_msa_format - %s\n" %str(ans))

        aln = vars['protein_MSA_clustalw']

    elif protein_MSA_clustalw in vars and os.path.exists(vars['working_dir'] + dir + vars['protein_MSA_clustalw']):

        # we are already in MSA created by clustalw, just the name is incorrect. for consistency we copy it to be with the desired name
        try:

            os.rename(vars['working_dir'] + dir + vars['protein_MSA'], vars['working_dir'] + dir + vars['protein_MSA_clustalw'])

        except:

            exit_on_error('sys_error', "create_chimera: change name failed\n")

        aln = vars['protein_MSA_clustalw']

    LOG.write("MSA: %s\tFORMAT:%s\n" %(aln, vars['msa_format']))
    """
    LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.color_with_chimera(%s, %s, %s, %s, %s, %s, %s)\n" %(vars['msa_SEQNAME'], aln, vars['r4s_out'], vars['scf_for_chimera'], vars['header_for_chimera'], vars['isd_scf_for_chimera'], vars['isd_header_for_chimera']))
    ans = cp_rasmol_gradesPE_and_pipe.color_with_chimera(vars['msa_SEQNAME'], aln, vars['r4s_out'], vars['scf_for_chimera'], vars['header_for_chimera'], vars['isd_scf_for_chimera'], vars['isd_header_for_chimera'])
    if ans[0] != "OK":

        exit_on_error('sys_error', "create_chimera: %s\n" %str(ans))

    else:

        vars['insufficient_data'] = ans[1]

    # C. creating the script that shows the MSA, Tree and colored 3D structure
    LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.create_chimera_script_align_tree(%s, %s, %s, %s, %s, %s)\n" %(vars['run_url'], vars['chimerax_file'], aln, vars['tree_file'], vars['scf_for_chimera'], vars['header_for_chimera']))
    cp_rasmol_gradesPE_and_pipe.create_chimera_script_align_tree(vars['run_url'], vars['chimerax_file'], aln, vars['tree_file'], vars['scf_for_chimera'], vars['header_for_chimera'])

    # create also ignoring insufficient data file, only in case the selected algorithm was bayes
    if vars['insufficient_data'] == "yes" and form['algorithm'] == "Bayes":

        LOG.write("Calling: cp_rasmol_gradesPE_and_pipe.create_chimera_script_align_tree(%s, %s, %s, %s, %s, %s)\n" %(vars['run_url'], vars['isd_chimerax_file'], aln, vars['tree_file'], vars['isd_scf_for_chimera'], vars['isd_header_for_chimera']))
        cp_rasmol_gradesPE_and_pipe.create_chimera_script_align_tree(vars['run_url'], vars['isd_chimerax_file'], aln, vars['tree_file'], vars['isd_scf_for_chimera'], vars['isd_header_for_chimera'])

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

def Create_Colored_MSA():

    """
    msa_format = MSA_parser.determine_msa_format(vars['protein_MSA']) #TO CHECK
    if msa_format[0] == "err":

        msa_info_msg = "<a href =\"" + GENERAL_CONSTANTS.MSA_FORMATS + "\">Read more on MSA formats</a><br />\n"
        exit_on_error('user_error', "The uploaded <a href=\"" + vars['user_msa_file_name'] + "\">MSA file</a> is not in one of the formats supported by ConSurf: NBRF/PIR, Pearson (Fasta), Nexus, Clustal, GCG/MSF.<br /><br />\nPlease check the following items and try to run ConSurf again:<br />\n1. The file should be saved as plain text (e.g. file type 'txt' in windows or 'MS-Dos' from Word in Mac).<br />\n2. The file should not contain unnecessary characters (You can check it with 'Notepad' editor).<br />\n3. The same sequence name must not be repeated more then once.<br />\n" + msa_info_msg)

    else:

        LOG.write("determine_msa_format : MSA format is : %s\n" %msa_format[1])
        vars['msa_format'] = msa_format[1]

    # NEW CORRECTED ONE
    """
    msa_colored_css = GENERAL_CONSTANTS.CONSURF_URL + "MSA_Colored.NEW.EM.css"
    """
    if vars['msa_format'] != "fasta":

        LOG.write("Create_Colored_MSA Calls 'MSA_parser.convert_msa_format(%s, %s, %s, %s)'\n" %(vars['protein_MSA'], vars['msa_format'], vars['protein_MSA'] + ".fs", "fasta"))
        ans = MSA_parser.convert_msa_format(vars['protein_MSA'], vars['msa_format'], vars['protein_MSA'] + ".fs", "fasta")
        if ans != "OK":

            exit_on_error('user_error', "Create_Colored_MSA: MSA_parser.convert_msa_format - %s\n" %ans)

        LOG.write("ConSeq_gradesPE_and_Outputs.print_msa_colors_FASTA_clustalwLike(%s, %s, %s, %s, %s, %s)\n" %(gradesPE_Output, vars['protein_MSA'] + ".fs", vars['msa_SEQNAME'], vars['Colored_MSA_HTML'], msa_colored_css, "ConSurf Color-Coded MSA"))
        ConSeq_gradesPE_and_Outputs.print_msa_colors_FASTA_clustalwLike(gradesPE_Output, vars['protein_MSA'] + ".fs", vars['msa_SEQNAME'], vars['Colored_MSA_HTML'], msa_colored_css, "ConSurf Color-Coded MSA")

    elif vars['msa_format'] == "fasta":
    """
    LOG.write("ConSeq_gradesPE_and_Outputs.print_msa_colors_FASTA_clustalwLike(%s, %s, %s, %s, %s, %s)\n" %(gradesPE_Output, vars['msa_fasta'] , vars['msa_SEQNAME'], vars['Colored_MSA_HTML'], msa_colored_css, "ConSurf Color-Coded MSA"))
    ConSeq_gradesPE_and_Outputs.print_msa_colors_FASTA_clustalwLike(gradesPE_Output, vars['msa_fasta'], vars['msa_SEQNAME'], vars['Colored_MSA_HTML'], msa_colored_css, "ConSurf Color-Coded MSA")

def predict_solvent_accesibility(aln, hssp, pred):

    # runs the PACC algorithm to calculate burried/exposed

    query_name = "\'" + vars['msa_SEQNAME'] + "\'"
    MSA_to_HSSP = GENERAL_CONSTANTS.MSA_to_HSSP
    PACC_path = GENERAL_CONSTANTS.PACC_path
    Predict_PACC = GENERAL_CONSTANTS.PREDICT_PACC

    HSSP = open(hssp, "w")
    HSSP.close()
    PRED = open(pred, "w")
    PRED.close()

    # run the script that turns the MSA to hssp file
    cmd = ["ssh", "barakyariv@power", "(", MSA_to_HSSP, query_name, vars['working_dir'] + aln, ")", ">", vars['working_dir'] + hssp]
    LOG.write("predict_solvent_accesibility:\n%s\n" %cmd)
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    p.wait()

    # check if the hssp file has been created and is not empty
    while not os.path.exists(hssp) or os.path.getsize(hssp) == 0:

        continue
        #exit_on_error('sys_error', "predict_solvent_accesibility: The file " + hssp + " does not exist or contains no data\n")

    LOG.write("predict_solvent_accesibility: The file " + hssp + " was created and now contains data\n")

    # run the script that produce the prediction file
    cmd = ["ssh", "barakyariv@power", "cd", PACC_path + ";", Predict_PACC, vars['working_dir'] + hssp, vars['working_dir'] + pred]
    LOG.write("predict_solvent_accesibility: %s\n" %cmd)
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    p.communicate()

    # check if the pred file has been created and is not empty
    if not os.path.exists(pred) or os.path.getsize(pred) == 0:

        exit_on_error('sys_error', "predict_solvent_accesibility: The file " + pred + " does not exist or contains no data\n")


def Update_Progress(ProgressFile, message):

    PROGRESS = open(ProgressFile, 'r')
    lines = PROGRESS.readlines()
    PROGRESS.close()

    PROGRESS = open(ProgressFile, 'w')

    for line in lines:

        if re.search(message, line):

            PROGRESS.write(re.sub(r'in_progress', r'finished', line))

        else:

            PROGRESS.write(line)


    PROGRESS.close()


def find_best_substitution_model(msa_file_path, AA_or_NUC):

    # convert fasta to phylip
    msa_phy_filepath = "input_msa.phy"
    AlignIO.convert(msa_file_path, "fasta", msa_phy_filepath, "phylip-relaxed")
    model = ""

    if AA_or_NUC == "NUC":

        model = run_jmt(msa_phy_filepath)

    else:

        model = run_protest(msa_phy_filepath)

    model = model.strip('()')
    if model == "JC" and AA_or_NUC == "NUC":

        model = "JC_Nuc"

    return model

def run_protest(msa_file_path):

    PRT_JAR_FILE = "/bioseq/Programs/ModelTest/prottest-3.4.1/prottest-3.4.1.jar"
    output_file_path = "model_selection.out"
    cmd = ["java", "-jar", PRT_JAR_FILE, "-log", "disabled", "-i", msa_file_path, "-AICC", "-o", output_file_path, "-S", "1", "-JTT", "-LG", "-MtREV", "-Dayhoff", "-WAG", "-CpREV", "-threads", "1"]
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    out, err = p.communicate()
    LOG.write("find_best_substitution_model: %s\n" %str(cmd))
    LOG.write(out + "\n")

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
    output_file_path = "model_selection.out"
    cmd = ["java", "-jar", JMT_JAR_FILE, "-d", msa_file_path, "-t", "BIONJ", "-AICC", "-f", "-o", output_file_path]
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    out, err = p.communicate()
    LOG.write("find_best_substitution_model: %s\n" %str(cmd))
    LOG.write(out + "\n")

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


def create_gradesPE_ConSurf(gradesPE, scoresFile, ref_r4s2pdb, gradesPE_Output_ArrRef, residue_freq_HashRef, no_isd_residue_color_ArrRef, isd_residue_color_ArrRef, Type):

    ans = cp_rasmol_gradesPE_and_pipe.create_gradesPE_ConSurf(gradesPE_Output_ArrRef, ref_r4s2pdb, residue_freq_HashRef, no_isd_residue_color_ArrRef, isd_residue_color_ArrRef, gradesPE, scoresFile, Type, vars['D_Nuc'])

    if not ans[0] == "OK":

        exit_on_error('sys_error', "create_gradesPE_ConSurf : cp_rasmol_gradesPE_and_pipe." + str(ans))

    elif not os.path.exists(gradesPE) or os.path.getsize(gradesPE) == 0:

        exit_on_error('sys_error', "create_gradesPE_ConSurf : the file '" + gradesPE + "' was not found or empty")

    if ans[1] == "" or ans[2] == "":

        exit_on_error('sys_error', "create_gradesPE_ConSurf : there is no data in the returned values seq3d_grades_isd or seq3d_grades from the routine")

    return(ans[1], ans[2])
"""
def create_part_of_pipe_new(pipeFile, pipeFile_CBS, seq3d_grades, seq3d_grades_isd, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, length_of_seqres, length_of_atom, pdb_file_name, user_chain, IN_pdb_id_capital):

    # CREATE PART of PIPE

    partOfPipe = vars['working_dir'] + dir + "partOfPipe"
    partOfPipe_CBS = vars['working_dir'] + dir + "partOfPipe_CBS"
    LOG.write("Calling cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)" %(partOfPipe, vars['unique_seqs'], form['database'], seq3d_grades_isd, seq3d_grades, length_of_seqres, length_of_atom, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, form['ESCORE'], form['iterations'], form['MAX_NUM_HOMOL'], form['MSAprogram'], form['algorithm'], form['matrix'], "legacy"))
    ans = cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new(partOfPipe, vars['unique_seqs'], form['database'], seq3d_grades_isd, seq3d_grades, length_of_seqres, length_of_atom, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, form['ESCORE'], form['iterations'], form['MAX_NUM_HOMOL'], form['MSAprogram'], form['algorithm'], form['matrix'], "legacy")
    if not ans[0] == "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new FAILED: " + str(ans))

    elif not os.path.exists(partOfPipe) or os.path.getsize(partOfPipe) == 0:

        exit_on_error('sys_error', "create_pipe_file: The file '" + partOfPipe + "' was not found or empty")

    # create the color blind friendly version
    ans = cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new(partOfPipe_CBS, vars['unique_seqs'], form['database'], seq3d_grades_isd, seq3d_grades, length_of_seqres, length_of_atom, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, form['ESCORE'], form['iterations'], form['MAX_NUM_HOMOL'], form['MSAprogram'], form['algorithm'], form['matrix'], "cb")
    if not ans[0] == "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new FAILED: " + str(ans))

    elif not os.path.exists(partOfPipe_CBS) or os.path.getsize(partOfPipe_CBS) == 0:

        exit_on_error('sys_error', "create_pipe_file: The file '" + partOfPipe_CBS + "' was not found or empty")
"""
def create_pipe_file(pipeFile, pipeFile_CBS, seq3d_grades, seq3d_grades_isd, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, length_of_seqres, length_of_atom, pdb_file_name, user_chain, IN_pdb_id_capital):

    # CREATE PART of PIPE
    partOfPipe = "partOfPipe"
    partOfPipe_CBS = "partOfPipe_CBS"

    LOG.write("Calling cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n" %(partOfPipe, vars['unique_seqs'], form['database'], seq3d_grades_isd, seq3d_grades, length_of_seqres, length_of_atom, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, form['ESCORE'], form['iterations'], form['MAX_NUM_HOMOL'], form['MSAprogram'], form['algorithm'], form['matrix'], "legacy"))
    ans = cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new(partOfPipe, vars['unique_seqs'], form['database'], seq3d_grades_isd, seq3d_grades, length_of_seqres, length_of_atom, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, form['ESCORE'], form['iterations'], form['MAX_NUM_HOMOL'], form['MSAprogram'], form['algorithm'], form['matrix'], "legacy")
    if ans[0] != "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new FAILED: " + str(ans))

    elif not os.path.exists(partOfPipe) or os.path.getsize(partOfPipe) == 0:

        exit_on_error('sys_error', "create_pipe_file: The file '" + partOfPipe + "' was not found or empty")

    # create the color blind friendly version
    ans = cp_rasmol_gradesPE_and_pipe.create_part_of_pipe_new(partOfPipe_CBS, vars['unique_seqs'], form['database'], seq3d_grades_isd, seq3d_grades, length_of_seqres, length_of_atom, isd_residue_color_ArrRef, no_isd_residue_color_ArrRef, form['ESCORE'], form['iterations'], form['MAX_NUM_HOMOL'], form['MSAprogram'], form['algorithm'], form['matrix'], "cb")
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
    if form['uploaded_MSA'] != "":

        msa_filename = form['uploaded_MSA']

    tree_filename = ""
    if form['uploaded_TREE'] != "":

        tree_filename = form['uploaded_TREE']

    msa_query_seq_name = ""
    if form['msa_SEQNAME'] != "":

        msa_query_seq_name = form['msa_SEQNAME']

    # GET THE CURRENT TIME
    completion_time = str(datetime.datetime.now().time())
    run_date = str(datetime.datetime.now().date())

    # FIND IDENTICAL CHAINS
    identical_chains = find_identical_chains_in_PDB_file(pdb_file_name, user_chain)

    # USE THE CREATED PART of PIPE to CREATE ALL THE PIPE TILL THE PDB ATOMS (DELETE THE PART PIPE)
    LOG.write("Calling cp_rasmol_gradesPE_and_pipe.create_consurf_pipe_new(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n" %(vars['working_dir'], IN_pdb_id_capital, user_chain, header_pipe, pipeFile, identical_chains, partOfPipe, vars['working_dir'], vars['run_number'], msa_filename, msa_query_seq_name, tree_filename, vars['submission_time'], completion_time, run_date))
    ans = cp_rasmol_gradesPE_and_pipe.create_consurf_pipe_new(vars['working_dir'], IN_pdb_id_capital, user_chain, header_pipe, pipeFile, identical_chains, partOfPipe, vars['working_dir'], vars['run_number'], msa_filename, msa_query_seq_name, tree_filename, vars['submission_time'], completion_time, run_date)
    if ans != "OK":

        exit_on_error('sys_error', "cp_rasmol_gradesPE_and_pipe.create_consurf_pipe_new FAILED: " + ans)

    # USE THE CREATED PART of PIPE to CREATE ALL THE PIPE TILL THE PDB ATOMS (DELETE THE PART PIPE) - Color friendly version
    LOG.write("Calling cp_rasmol_gradesPE_and_pipe.create_consurf_pipe_new(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n" %(vars['working_dir'], IN_pdb_id_capital, user_chain, header_pipe, pipeFile_CBS, identical_chains, partOfPipe_CBS, vars['working_dir'], vars['run_number'], msa_filename, msa_query_seq_name, tree_filename, vars['submission_time'], completion_time, run_date))
    ans = cp_rasmol_gradesPE_and_pipe.create_consurf_pipe_new(vars['working_dir'], IN_pdb_id_capital, user_chain, header_pipe, pipeFile_CBS, identical_chains, partOfPipe_CBS, vars['working_dir'], vars['run_number'], msa_filename, msa_query_seq_name, tree_filename, vars['submission_time'], completion_time, run_date)
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

    return(identical_chains)
"""
def find_identical_chains_in_PDB_file(pdb_file, chain):

    # find identical chains by extracting chains from ATOMS file and then pairwise align
    LOG.write("find_identical_chains_in_PDB_file\n")

    if not os.path.exists(pdb_file):

        print_message_to_output( "find_identical_chains_in_PDB_file - File " + pdb_file + " does not Exists!")
        LOG.write("find_identical_chains_in_PDB_file - File " + pdb_file + " does not Exists!")
"""

def find_identical_chains_in_PDB_file(filename, OrgChain):

    # Looking in the PDB for chains identical to the original chain

    # this array shows how to convert three letters to one letter
    arrayOneLetter = {'VAL' : 'V', 'ILE' : 'I', 'LEU' : 'L', 'GLU' : 'E', 'GLN' : 'Q', 'ASP' : 'D', 'ASN' : 'N', 'HIS' : 'H', 'TRP' : 'W', 'PHE' : 'F', 'TYR' : 'Y', 'ARG' : 'R', 'LYS' : 'K', 'SER' : 'S', 'THR' : 'T', 'MET' : 'M', 'ALA' : 'A', 'GLY' : 'G', 'PRO' : 'P', 'CYS' : 'C'}

    # open PDB to get chains
    try:

        PDB = open(filename, 'r')

    except:

        exit_on_error('sys_error', "find_identical_chains_in_PDB_file: can't open the file " + filename + " for reading.")

    # fill hash with chains
    HashChainsSeq = {}
    line = PDB.readline()
    while line != "":

        # skip non ATOM recs
        if not re.match(r'^ATOM', line):

            line = PDB.readline()
            continue

        # skip non CA recs
        if not line[13:15] == "CA":

            line = PDB.readline()
            continue

        chain = line[21:22]
        if chain in HashChainsSeq:

            HashChainsSeq[chain] += arrayOneLetter[line[17:20]]

        else:

            HashChainsSeq[chain] = arrayOneLetter[line[17:20]]

        line = PDB.readline()

    PDB.close()

    # string with identical chains
    identicalChains = OrgChain

    # looking for chains identical to the original chain
    for chain in HashChainsSeq:

        if OrgChain != chain:

            chain_length = len(HashChainsSeq[chain])
            OrgChain_length = len(HashChainsSeq[OrgChain])

            # if length not similar, skip
            if min(OrgChain_length, chain_length)/max(OrgChain_length, chain_length) <= 0.9:

                continue

            # compare the two chains with clustalw

            # file for two fastas
            fastaFileName = OrgChain + "_" + chain + "_twoFastas"
            # file with clustalw output
            clustalwOutputFile = OrgChain + "_" + chain + "_twoFastas.aln"

            try:

                FASTAS = open(fastaFileName, 'w')

            except:

                exit_on_error('sys_error', "find_identical_chains_in_PDB_file: can't open the file " + fastaFileName + " for writing.")

            FASTAS.write(">%s\n%s\n>%s\n%s\n" %(OrgChain, HashChainsSeq[OrgChain], chain, HashChainsSeq[chain]))
            FASTAS.close()

            cmd = ["/bioseq/Programs/ClustalW_2.0.10/clustalw-2.0.10-linux-i386-libcppstatic/clustalw2", "-INFILE=" + fastaFileName, "-gapopen=1", "-OUTPUT=" + clustalwOutputFile]
            print(cmd)
            # run clustalw
            p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
            out, err = p.communicate()
            print(err)

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

                identicalChains += chain

    return identicalChains

def Get_NACSES_buried_Exposed(structure_file, chain, buried_cutoff, exposed_cutoff, atom_positionFILE = ""):

    ATOM_pos_To_SEQRES = {} # will have values if atom_positionFILE provided
    Buried_Exposed = {} # key: poistion in model, value: b|e according to the cutoff
    NACSSESS_DIR = GENERAL_CONSTANTS.NACSSESS_DIR
    cmd = [NACSSESS_DIR + "naccess", structure_file]
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    p.communicate()

    basename = structure_file[:-4] # delete ending
    out_Files = [basename + ".asa", basename + ".rsa", basename + ".log"]
    for file in out_Files:

        if not os.path.exists(file):

            LOG.write("ERROR, command: '%s' FAILED to create %s\n" %(str(cmd), file))

    if atom_positionFILE != "":

        try:

            ATOM_TO_SEQRES = open(atom_positionFILE, 'r')

        except:

            return("failed", "Get_NACSES_buried_Exposed: could not open: " + atom_positionFILE + " for writing.")

        line = ATOM_TO_SEQRES.readline()
        while line != "":

            [Res, SEQRES_POS, ATOM_POS] = line.rstrip().split()

            if not ATOM_POS in ATOM_pos_To_SEQRES:

                ATOM_pos_To_SEQRES[ATOM_POS] = {}

            ATOM_pos_To_SEQRES[ATOM_POS]['RES'] = Res
            ATOM_pos_To_SEQRES[ATOM_POS]['SEQRES'] = SEQRES_POS

            line = ATOM_TO_SEQRES.readline()

        ATOM_TO_SEQRES.close()

    # parse rsa file for b/e
    try:

        RSA = open(basename + ".rsa", 'r')

    except:

        return("failed", "Get_NACSES_buried_Exposed: could not open: " + basename + ".rsa for reading.")

    try:

        OUT = open(basename + ".NACSSES.rsa", 'w')

    except:

        return("failed", "Get_NACSES_buried_Exposed: could not open: " + basename + ".NACSSES.rsa for writing.")

    OUT.write("ATOM_POS\tAll_atoms_REL\tb/e")
    if atom_positionFILE != "":

        OUT.write("\tSEQRES_POS\tRES\n")

    else:

        OUT.write("\n")

    line = RSA.readline()
    while line != "":

        words = line.split()
        if len(words) >= 13: # informative line

            if len(words) == 13: # pdb file with no chain id

                pos = words[2]
                if atom_positionFILE != "" and pos in ATOM_pos_To_SEQRES and 'SEQRES' in ATOM_pos_To_SEQRES[pos]:

                    pos = ATOM_pos_To_SEQRES[pos]['SEQRES']

                if float(words[4]) >= exposed_cutoff: # exposed according to All atoms REL

                    Buried_Exposed[pos] = "e"

                elif float(words[4]) < buried_cutoff: # buried according to All atoms REL

                    Buried_Exposed[pos] = "b"

                OUT.write("%s\t%s\t%s" %(words[2], words[4], Buried_Exposed[pos]))

                if atom_positionFILE != "" and words[2] in ATOM_pos_To_SEQRES and 'SEQRES' in ATOM_pos_To_SEQRES[words[2]]:

                    OUT.write("\t%s\t%s\n" %(ATOM_pos_To_SEQRES[words[2]]['SEQRES'], ATOM_pos_To_SEQRES[words[2]]['RES']))

                else:

                    OUT.write("\n")

            elif chain == words[2]: # the relevant chain

                pos = words[3]
                if atom_positionFILE != "" and pos in ATOM_pos_To_SEQRES and 'SEQRES' in ATOM_pos_To_SEQRES[pos]:

                    pos = ATOM_pos_To_SEQRES[pos]['SEQRES']

                if float(words[5]) >= exposed_cutoff: # exposed according to All atoms REL

                    Buried_Exposed[pos] = "e"

                elif float(words[5]) < buried_cutoff: # buried according to All atoms REL

                    Buried_Exposed[pos] = "b"

                OUT.write("%s:%s\t%s\t%s" %(words[2], words[3], words[5], Buried_Exposed[pos]))

                if atom_positionFILE != "" and words[3] in ATOM_pos_To_SEQRES and 'SEQRES' in ATOM_pos_To_SEQRES[words[3]]:

                    OUT.write("\t%s\t%s\n" %(ATOM_pos_To_SEQRES[words[3]]['SEQRES'], ATOM_pos_To_SEQRES[words[3]]['RES']))

                else:

                    OUT.write("\n")

        line = RSA.readline()

    RSA.close()
    OUT.close()

    return("OK", Buried_Exposed)



def Project_ConSurf_On_Model(PDB_File, chain, OutDir, Query_Seq, r4s_out, protein_MSA_clustalw, PDB_Name_For_Out, Type):

    PDB_Object = pdbParser.pdbParser()
    PDB_Object.read(PDB_File)

    # This function is used either for PISA or for HHPred
    pisa = False
    if 'PISA1' in PDB_File:

        pisa = True

    prefix = PDB_File + "_" + chain # this begins many of the variable names

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

    [SEQRES_seq_local, ATOM_seq_loacl] = get_seqres_atom_seq(PDB_Object, chain, PDB_File)
    if SEQRES_seq_local == "user_error":

        return("<font size=+1 color='red'>ERROR! " + ATOM_seq_loacl + "</font><br>")

    ATOM_Pos_File = prefix + ".ATOM_Pos_File.txt" # created by create_atom_position_file()
    pairwise_aln = prefix + "_QUERY.aln" # created by compare_atom_to_query()

    rasmolFILE = prefix + ".rasmol.scr"
    rasmol_isdFILE = prefix + ".isd_rasmol.scr"

    rasmolFILE_CBS = prefix + ".rasmol_CBS.scr"
    rasmol_isdFILE_CBS = prefix + ".isd_rasmol_CBS.scr"
    gradesPE = prefix + ".consurf.grades"
    scoresFile = prefix + ".consurf.scores"
    if pisa:

        jobInfoFile = "PISA1.job_info.json"

    else:

        jobInfoFile = "job_info.json"

    prefix = PDB_File + "_consurf_" + vars['run_number']

    pipeFile = prefix  + "_pipe.pdb"
    pipeFile_CBS = prefix  + "_pipe_CBS.pdb"
    Server_Results_Path = "results" + dir + vars['run_number'] + dir # in the global dir

    # chimera files
    chimerax_script_for_figure = prefix + "_Figure.chimerax"
    chimerax_script_for_figure_isd = prefix + "_Figure_isd.chimerax"

    chimera_color_script = GENERAL_CONSTANTS.CONSURF_URL + "chimera/chimera_consurf.cmd"
    chimera_instructions = PDB_File + "_chimera_instructions.php"

    pymol_color_script = GENERAL_CONSTANTS.CONSURF_URL + "pyMOL/consurf_new.py"
    pymol_instructions = PDB_File + "_PyMol_instructions.php"

    scf_for_chimera = prefix + ".scf"
    header_for_chimera = prefix + ".hdr"
    chimerax_file = prefix + ".chimerax"

    isd_scf_for_chimera = prefix + "_isd.scf"
    isd_header_for_chimera = prefix + "_isd.hdr"
    isd_chimerax_file = prefix + "_isd.chimerax"

    # Atoms Section with consurf grades instead TempFactor Field
    ATOMS_with_ConSurf_Scores = PDB_File + "_ATOMS_section_With_ConSurf.pdb"
    ATOMS_with_ConSurf_Scores_isd =  PDB_File + "_ATOMS_section_With_ConSurf_isd.pdb"

    pdb_file_with_score_at_TempFactor= PDB_File + "_With_Conservation_Scores.pdb"

    assign_colors_according_to_r4s_layers(gradesPE_Output_local, r4s_out)
    [Query_Seq_with_gaps, ATOM_seq_loacl_with_gaps] = compare_atom_to_query(Query_Seq, ATOM_seq_loacl, pairwise_aln, PDB_File + "_" + chain) # Compare the PDB Atom With Reference Seq
    read_residue_variety(residue_freq_local, position_totalAA_local)
    res = create_atom_position_file(chain, PDB_File, ATOM_Pos_File, Type, pisa) # this file will be used later to create the output which aligns rate4site sequence with the ATOM records
    if res == 0:

        return("", "")

    r4s2pdb_local = {} # key: poistion in SEQRES/MSA, value: residue name with position in atom (i.e: ALA22:A)
    [length_of_seqres_loacl, length_of_atom_loacl] = match_pdb_to_seq(r4s2pdb_local, chain, Query_Seq_with_gaps, ATOM_seq_loacl_with_gaps, ATOM_Pos_File, Type)

    # The following routine, apart from creating the file "consurf.grades" also collects the information in order to create
    # the RasMol scripts and a variable which will be used in the "pipe" file, that holds a string with the grades.
    # In the pipe file this var is called: seq3d_grades_isd and seq3d_grades
    [seq3d_grades_isd_local, seq3d_grades_isd_local] = create_gradesPE_ConSurf(gradesPE, scoresFile, r4s2pdb_local, gradesPE_Output_local, residue_freq_local, no_isd_residue_color_local, isd_residue_color_loacl, Type)

    # This will create the 2 rasmol scripts (one with isd and one without)
    create_rasmol(chain, rasmolFILE, rasmol_isdFILE, rasmolFILE_CBS, rasmol_isdFILE_CBS, no_isd_residue_color_local, isd_residue_color_loacl)

    # This will create the pipe file for FGiJ
    identical_chains = create_pipe_file(pipeFile, pipeFile_CBS, seq3d_grades_isd_local, seq3d_grades_isd_local, isd_residue_color_loacl, no_isd_residue_color_local, length_of_seqres_loacl, length_of_atom_loacl, PDB_File, chain, PDB_File.upper())
    """
    if pisa:

        saveJobInfoJson(jobInfoFile, scoresFile, PDB_File, "", form['job_title'], "PISA", "", chain, identical_chains)

    else:

        saveJobInfoJson(jobInfoFile, scoresFile, PDB_File, "", form['job_title'], "", "", "A", "A")
    """

    if pisa:

        Hash_Json = {'scoresFile' : scoresFile, 'pdbFile' : PDB_File, 'jobTitle' : form['job_title'], 'subTitle' : "PISA", 'chainId' : chain, 'identicalChains' : identical_chains}
        saveJobInfoJson(jobInfoFile, Hash_Json)

    else:

        Hash_Json = {'scoresFile' : scoresFile, 'pdbFile' : PDB_File, 'jobTitle' : form['job_title'], 'chainId' : "A", 'identicalChains' : "A"}
        saveJobInfoJson(jobInfoFile, Hash_Json)

    # Create ATOMS section and replace the TempFactor Column with the ConSurf Grades (will create also isd file if relevant)
    replace_TmpFactor_Consurf_Scores(chain, PDB_File, gradesPE, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd)

    # Create PDB file and replace the TempFactor Column with the Rate4Site grades
    replace_TmpFactor_Rate4Site_Scores(chain, PDB_File, gradesPE, pdb_file_with_score_at_TempFactor)

    # This Will create the script for chimera coloring
    insufficient_data = create_chimera(protein_MSA_clustalw, OutDir, r4s_out, scf_for_chimera, header_for_chimera, isd_scf_for_chimera, isd_header_for_chimera, ATOMS_with_ConSurf_Scores, chimerax_file, ATOMS_with_ConSurf_Scores_isd, isd_chimerax_file, chimerax_script_for_figure, chimerax_script_for_figure_isd, chimera_instructions, identical_chains, chain, PDB_File)

    # This will create the PyMol instruction page
    create_pymol_instructions(insufficient_data, pymol_instructions, PDB_File.upper(), chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, identical_chains, OutDir, chimerax_script_for_figure, chimerax_script_for_figure_isd, PDB_File)

    HTML_TXT = ""
    HTML_TXT_FGIJ = ""

    pipeUrlNoSuffix = Server_Results_Path + dir + PDB_File + "_consurf_" + vars['run_number']
    if pisa:

        HTML_TXT_FGIJ += buildHtml_3dViewerLinks(vars['run_number'], "PISA1", pipeUrlNoSuffix, " projected on " + PDB_Name_For_Out)

    else:

        HTML_TXT_FGIJ += buildHtml_3dViewerLinks(vars['run_number'], "", pipeUrlNoSuffix, " projected on " + PDB_Name_For_Out)

    prefix = "<A HREF='" + dir + vars['Server_Results_Path'] + dir

    HTML_TXT += buildHtml_OutMsg(prefix + gradesPE + "' TARGET=Conservation_window>Amino Acid Conservation Scores, Confidence Intervals and Conservation Colors</A>")
    HTML_TXT += "<h4 class=output_title>RasMol Coloring Scripts</h4>"
    HTML_TXT += buildHtml_OutMsg(prefix + rasmol_isdFILE + "' TARGET=spt_window> RasMol Coloring Script Showing Insufficient Data</A>")
    HTML_TXT += buildHtml_OutMsg(prefix + rasmolFILE + "' TARGET=spt_window> RasMol Coloring Script Hiding Insufficient Data</A><br>")
    HTML_TXT += "<h4 class=output_title>PDB Files</h4>"
    HTML_TXT += buildHtml_OutMsg(prefix + pdb_file_with_score_at_TempFactor + "' TARGET=PDB_window> PDB File with Conservation Scores in the tempFactor field</A>")
    HTML_TXT += buildHtml_OutMsg(prefix + PDB_File + "_consurf_" + vars['run_number'] + "_pipe.pdb' TARGET=PDB_window> PDB File with ConSurf Results in its Header, for FirstGlance in Jmol or Protein Explorer</A>; [<A href='" + vars['Server_Results_Path'] + dir + PDB_File + "_consurf_" + vars['run_number'] + "_pipe_CSB.pdb'>with color-blind friendly scale</A>]<br>")
    HTML_TXT += "<h4 class=output_title>Create a high resolution figures</h4>\n"
    HTML_TXT += buildHtml_OutMsg(prefix + chimera_instructions + "' TARGET=Chimera_HighRes_Instruct_window>Follow the instructions to produce a Chimera figure</A><font size=\"-1\"> (For users of Chimera)</font><br>")
    HTML_TXT += buildHtml_OutMsg(prefix + pymol_instructions + "' TARGET=PyMol_HighRes_Instruct_window>Follow the instructions to produce a PyMol figure</A><font size=\"-1\"> (For users of PyMol)</font><br>")

    return(HTML_TXT, HTML_TXT_FGIJ)


"""
def Project_ConSurf_On_Model(PDB_File, chain, OutDir, Query_Seq_File_User, r4s_out, protein_MSA_clustalw, PDB_Name_For_Out, Type):

    PDB_Object = pdbParser.pdbParser()
    PDB_Object.read(OutDir + PDB_File)

    # This function is used either for PISA or for HHPred
    pisa = False
    if 'PISA1' in PDB_File:

        pisa = True

    prefix = OutDir + PDB_File + "_" + chain # this begins many of the variable names

    # Constant file names
    ATOM_seq_File = prefix + ".ATOMS_SEQ.fs"
    Query_seq_File = prefix + ".QUERY_SEQ.fs"

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
    no_isd_residue_color_local = []
    isd_residue_color_loacl = []

    [SEQRES_seq_local, ATOM_seq_loacl] = get_seqres_atom_seq(PDB_Object, chain, PDB_File)
    if SEQRES_seq_local == "user_error":

        return("<font size=+1 color='red'>ERROR! " + ATOM_seq_loacl + "</font><br>")

    # create file with ATOM seq
    try:

        ATOM_SEQ = open(ATOM_seq_File, 'w')

    except:

        exit_on_error('sys_error', "Project_ConSurf_On_Model: could not open the file '" + ATOM_seq_File + "' for writing.")

    if ATOM_seq_loacl != "":

        ATOM_SEQ.write(">ATOM_%s_%s\n%s" %(PDB_File, chain, ATOM_seq_loacl))

    else:

        exit_on_error('sys_error', "Project_ConSurf_On_Model: The Atom Seq is Empty: %s_%s" %(PDB_File, chain))

    ATOM_SEQ.close()

    # create one line FASTA for QUERY
    try:

        QUERY = open(Query_Seq_File_User, 'r')

    except:

        exit_on_error('sys_error', "Project_ConSurf_On_Model: could not open the file '" + Query_Seq_File_User + "' for reading.")

    try:

        QUERY_SEQ = open(Query_seq_File, 'w')

    except:

        exit_on_error('sys_error', "Project_ConSurf_on_Model: could not open the file '" + Query_seq_File + "' for writing.")

    line = QUERY.readline()
    while line != "":

        if re.match(r'^>', line):

            QUERY_SEQ.write(">QUERY\n")

        else:

            QUERY_SEQ.write(line.rstrip())

        line = QUERY.readline()

    QUERY_SEQ.write("\n")

    QUERY.close()
    QUERY_SEQ.close()

    ATOM_Pos_File = prefix + ".ATOM_Pos_File.txt" # created by create_atom_position_file()
    pairwise_aln = prefix + "vs_QUERY.aln" # created by compare_atom_to_query()

    rasmolFILE = prefix + ".rasmol.scr"
    rasmol_isdFILE = prefix + ".isd_rasmol.scr"

    rasmolFILE_CBS = prefix + ".rasmol_CBS.scr"
    rasmol_isdFILE_CBS = prefix + ".isd_rasmol_CBS.scr"
    gradesPE = prefix + ".consurf.grades"
    scoresFile = prefix + ".consurf.scores"

    jobInfoFile = vars['working_dir'] + dir
    if pisa:

        jobInfoFile += "PISA1.job_info.json"

    else:

        jobInfoFile += "job_info.json"

    prefix = PDB_File + "_consurf_" + str(vars['run_number'])

    pipeFile = OutDir + dir + prefix  + "_pipe.pdb"
    pipeFile_CBS = OutDir + dir + prefix  + "_pipe_CBS.pdb"

    # chimera files
    chimerax_script_for_figure = prefix + "_Figure.chimerax"
    chimerax_script_for_figure_isd = prefix + "_Figure_isd.chimerax"

    chimera_color_script = GENERAL_CONSTANTS.CONSURF_URL + "chimera/chimera_consurf.cmd"
    chimera_instructions = PDB_File + "_chimera_instructions.php"

    pymol_color_script = GENERAL_CONSTANTS.CONSURF_URL + "pyMOL/consurf_new.py"
    pymol_instructions = PDB_File + "_PyMol_instructions.php"

    scf_for_chimera = prefix + ".scf"
    header_for_chimera = prefix + ".hdr"
    chimerax_file = prefix + ".chimerax"

    isd_scf_for_chimera = prefix + "_isd.scf"
    isd_header_for_chimera = prefix + "_isd.hdr"
    isd_chimerax_file = prefix + "_isd.chimerax"

    # Atoms Section with consurf grades instead TempFactor Field
    ATOMS_with_ConSurf_Scores = PDB_File + "_ATOMS_section_With_ConSurf.pdb"
    ATOMS_with_ConSurf_Scores_isd =  PDB_File + "_ATOMS_section_With_ConSurf_isd.pdb"

    pdb_file_with_score_at_TempFactor= PDB_File + "_With_Conservation_Scores.pdb"

    assign_colors_according_to_r4s_layers(gradesPE_Output_local, vars['working_dir'] + dir + r4s_out)
    compare_atom_to_query(OutDir + dir + Query_seq_File, OutDir + dir + ATOM_seq_File, pairwise_aln, OutDir + dir + PDB_File + "_" + chain) # Compare the PDB Atom With Reference Seq
    read_residue_variety(residue_freq_local, position_totalAA_local)
    res = create_atom_position_file(chain, PDB_File, ATOM_Pos_File, Type, pisa) # this file will be used later to create the output which aligns rate4site sequence with the ATOM records
    if res == 0:

        return("", "")

    r4s2pdb_local = {} # key: poistion in SEQRES/MSA, value: residue name with position in atom (i.e: ALA22:A)
    [length_of_seqres_loacl, length_of_atom_loacl] = match_pdb_to_seq(r4s2pdb_local, chain, pairwise_aln, ATOM_Pos_File, Type)

    # The following routine, apart from creating the file "consurf.grades" also collects the information in order to create
    # the RasMol scripts and a variable which will be used in the "pipe" file, that holds a string with the grades.
    # In the pipe file this var is called: seq3d_grades_isd and seq3d_grades
    [seq3d_grades_isd_local, seq3d_grades_isd_local] = create_gradesPE_ConSurf(OutDir + dir + gradesPE, OutDir + dir + scoresFile, r4s2pdb_local, gradesPE_Output_local, residue_freq_local, no_isd_residue_color_local, isd_residue_color_loacl, Type)

    # This will create the 2 rasmol scripts (one with isd and one without)
    create_rasmol(chain, OutDir + dir + rasmolFILE, OutDir + dir + rasmol_isdFILE, OutDir + dir + rasmolFILE_CBS, OutDir + dir + rasmol_isdFILE_CBS, no_isd_residue_color_local, isd_residue_color_loacl)

    # This will create the pipe file for FGiJ
    identical_chains = create_pipe_file(pipeFile, pipeFile_CBS, seq3d_grades_isd_local, seq3d_grades_isd_local, isd_residue_color_loacl, no_isd_residue_color_local, length_of_seqres_loacl, length_of_atom_loacl, isd_residue_color_loacl, no_isd_residue_color_local, length_of_seqres_loacl, length_of_atom_loacl, PDB_File, chain, PDB_File_Name.upper(), OutDir)

    if pisa:

        Hash_Json = {'scoresFile' : scoresFile, 'pdbFile' : PDB_File_Name, 'jobTitle' : form['job_title'], 'subTitle' : "PISA", 'chainId' : chain, 'identicalChains' : identical_chains}
        saveJobInfoJson(jobInfoFile, Hash_Json)

    else:

        Hash_Json = {'scoresFile' : scoresFile, 'pdbFile' : PDB_File_Name, 'jobTitle' : form['job_title'],'chainId' : "A", 'identicalChains' : "A"}
        saveJobInfoJson(jobInfoFile, Hash_Json)

    # Create ATOMS section and replace the TempFactor Column with the ConSurf Grades (will create also isd file if relevant)
    ace_TmpFactor_Consurf_Scores(chain, PDB_File, OutDir + dir + gradesPE, OutDir + dir + ATOMS_with_ConSurf_Scores, OutDir + dir + ATOMS_with_ConSurf_Scores_isd)

    # Create PDB file and replace the TempFactor Column with the Rate4Site grades
    replace_TmpFactor_Rate4Site_Scores(chain, PDB_File, OutDir + dir + gradesPE, OutDir + dir + pdb_file_with_score_at_TempFactor)

    # This Will create the script for chimera coloring
    insufficient_data = create_chimera(protein_MSA_clustalw, "clustalw", OutDir, r4s_out, scf_for_chimera, header_for_chimera, isd_scf_for_chimera, isd_header_for_chimera, ATOMS_with_ConSurf_Scores, chimerax_file, ATOMS_with_ConSurf_Scores_isd, isd_chimerax_file, chimerax_script_for_figure, chimerax_script_for_figure_isd, OutDir + dir + chimera_instructions, identical_chains, chain, PDB_File_Name)

    # This will create the PyMol instruction page
    create_pymol_instructions(insufficient_data, OutDir + dir + pymol_instructions, PDB_File_Name.upper(), chain, ATOMS_with_ConSurf_Scores, ATOMS_with_ConSurf_Scores_isd, identical_chains, OutDir, chimerax_script_for_figure, chimerax_script_for_figure_isd, PDB_File_Name)

    HTML_TXT = ""
    HTML_TXT_FGIJ = ""

    pipeUrlNoSuffix = Server_Results_Path + dir + PDB_File_Name + "_consurf_" + vars['run_number']
    if pisa:

        HTML_TXT_FGIJ += buildHtml_3dViewerLinks(vars['run_number'], "PISA1", pipeUrlNoSuffix, " projected on " + PDB_Name_For_Out)

    else:

        HTML_TXT_FGIJ += buildHtml_3dViewerLinks(vars['run_number'], "", pipeUrlNoSuffix, " projected on " + PDB_Name_For_Out)

    prefix = "<A HREF='" + dir + vars['Server_Results_Path'] + dir

    HTML_TXT += buildHtml_OutMsg(prefix + gradesPE + "' TARGET=Conservation_window>Amino Acid Conservation Scores, Confidence Intervals and Conservation Colors</A>")
    HTML_TXT += "<h4 class=output_title>RasMol Coloring Scripts</h4>"
    HTML_TXT += buildHtml_OutMsg(prefix + rasmol_isdFILE + "' TARGET=spt_window> RasMol Coloring Script Showing Insufficient Data</A>")
    HTML_TXT += buildHtml_OutMsg(prefix + rasmolFILE + "' TARGET=spt_window> RasMol Coloring Script Hiding Insufficient Data</A><br>")
    HTML_TXT += "<h4 class=output_title>PDB Files</h4>"
    HTML_TXT += buildHtml_OutMsg(prefix + pdb_file_with_score_at_TempFactor + "' TARGET=PDB_window> PDB File with Conservation Scores in the tempFactor field</A>")
    HTML_TXT += buildHtml_OutMsg(prefix + PDB_File_Name + "_consurf_" + vars['run_number'] + "_pipe.pdb' TARGET=PDB_window> PDB File with ConSurf Results in its Header, for FirstGlance in Jmol or Protein Explorer</A>; [<A href='" + vars['Server_Results_Path'] + dir + PDB_File_Name + "_consurf_" + vars['run_number'] + "_pipe_CSB.pdb'>with color-blind friendly scale</A>]<br>")
    HTML_TXT += "<h4 class=output_title>Create a high resolution figures</h4>\n"
    HTML_TXT += buildHtml_OutMsg(prefix + chimera_instructions + "' TARGET=Chimera_HighRes_Instruct_window>Follow the instructions to produce a Chimera figure</A><font size=\"-1\"> (For users of Chimera)</font><br>")
    HTML_TXT += buildHtml_OutMsg(prefix + pymol_instructions + "' TARGET=PyMol_HighRes_Instruct_window>Follow the instructions to produce a PyMol figure</A><font size=\"-1\"> (For users of PyMol)</font><br>")

    return(HTML_TXT, HTML_TXT_FGIJ)
"""
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

def Get_PISA_Complex(PDB_ID, out_file):

    cmd = ["wget", "https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimer.pdb?" + PDB_ID.lower() +":1,1", "-O", out_file]
    LOG.write("Get the first model of PISA: %s\n" %str(cmd))
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    p.communicate()

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

def saveJobInfoJson(jobInfoFile, hash):


    try:

        JSON = open(jobInfoFile ,'w')

    except:

        exit_on_error('sys_error', "Could not open the file " + jobInfoFile + " for writing.")

    json.dump(hash, JSON)

    JSON.close()


def compare_atom_to_query(Query_seq, ATOM_seq, pairwise_aln, PDB_Name):

    # in case there are both seqres and atom fields, checks the similarity between the 2 sequences.

    clustalw = GENERAL_CONSTANTS.CLUSTALW_2_0_10
    two_fastas = PDB_Name + "_QUERY.fasta2"
    clustalw_out = PDB_Name + "_QUERY.out"

    [first_seq, second_seq] = compare_two_fastas(two_fastas, Query_seq, ATOM_seq, clustalw_out, pairwise_aln)
    """
    LOG.write("compare_atom_to_query : run clustalw to see the match between ATOM and QUERY sequences\n")
    cmd = [clustalw, "-INFILE=" + two_fastas, "-gapopen=1", "-OUTFILE=" + pairwise_aln, ">", clustalw_out]
    LOG.write("compare_atom_seqres_or_msa : run %s\n" %str(cmd))
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    out, err = p.communicate()
    print(err)
    if not os.path.exists(pairwise_aln) or os.path.getsize(pairwise_aln) == 0 or not os.path.exists(clustalw_out) or os.path.getsize(clustalw_out) == 0:

        exit_on_error('sys_error', "compare_atom_seqres_or_msa : one of clustalw outputs were not create; %s or %s\n" %(pairwise_aln, clustalw_out))
    """
    return(first_seq, second_seq)

def compare_two_fastas(two_fastas, first_seq ,second_seq ,clustalw_out, clustalw_aln):

    clustalw = GENERAL_CONSTANTS.CLUSTALW_2_0_10

    try:

        FAS = open(two_fastas, 'w')

    except:

        exit_on_error('sys_error', "compare_two_fastas : Cannot open the file " + two_fastas + " for writing.")

    FAS.write(">first\n%s\n>second\n%s\n" %(first_seq ,second_seq))
    FAS.close()
    command = ["ssh", "barakyariv@power", "(", clustalw, vars['working_dir'] + two_fastas, "-gapopen=1", ")", ">", vars['working_dir'] + clustalw_out]
    LOG.write("compare_two_fastas : run %s\n" %command)
    p = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    out, err = p.communicate()
    print(err)

    if not os.path.exists(clustalw_out) or os.path.getsize(clustalw_out) == 0 or not os.path.exists(clustalw_aln) or os.path.getsize(clustalw_aln) == 0:

        exit_on_error('sys_error', "compare_two_fastas : one of clustalw outputs were not created; %s or %s" %(clustalw_out, clustalw_aln))

    # return the two sequences with gaps

    try:

        CLUSTALW_ALN = open(clustalw_aln, 'r')

    except:

        exit_on_error('sys_error', "compare_two_fastas : could not open " + clustalw_out + " for reading.")

    """
    CLUSTALW_ALN.readline() # skip first line, the name of the first sequence
    first_seq_with_gaps = CLUSTALW_ALN.readline()
    CLUSTALW_ALN.readline() # skip third line, the name of the second sequence
    second_seq_with_gaps = CLUSTALW_ALN.readline()

    """
    first_seq_with_gaps = ""
    second_seq_with_gaps = ""

    line = CLUSTALW_ALN.readline()
    while line != "":

        match1 = re.match(r'^first\s+(\S*)', line)
        if match1:

            first_seq_with_gaps += match1.group(1)

        else:

            match2 = re.match(r'^second\s+(\S*)', line)
            if match2:

                second_seq_with_gaps += match2.group(1)

        line = CLUSTALW_ALN.readline()

    CLUSTALW_ALN.close()

    return(first_seq_with_gaps, second_seq_with_gaps)


"""
def compare_atom_to_query(Query_seq_File, PDB_seq_File, pairwise_aln, PDB_Name):

    # in case there are both seqres and atom fields, checks the similarity between the 2 sequences.
    two_fastas = PDB_Name + "_QUERY.fasta2"
    clustalw_out = PDB_Name + "_QUERY.out"
    vars['pairwise_aln'] = PDB_Name + "_QUERY.aln"

    # run clustalw to see the match between ATOM and SEQRES sequences
    if not os.path.exists(Query_seq_File) or os.path.getsize(Query_seq_File) == 0 or not os.path.exists(PDB_seq_File) or .os.path.getsize(PDB_seq_File) == 0:

        cmd = ["cat", Query_seq_File, PDB_seq_File, ">", two_fastas]
        subprocess.call(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")

    else:

        exit_on_error('sys_error', "compare_atom_to_query : one of the inputs are missing; %s or %s\n" %(PDB_seq_File, Query_seq_File))

    LOG.write("compare_atom_to_query : run clustalw to see the match between ATOM and QUERY sequences\n")
    cmd = ["clustalw", "-INFILE=" + two_fastas, "-gapopen=1", "-OUTFILE=" + pairwise_aln, ">", clustalw_out]
    LOG.write("compare_atom_seqres_or_msa : run %s\n" %str(cmd))
    subprocess.call(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding = "utf-8")
    if not os.path.exists(pairwise_aln) or os.path.getsize(pairwise_aln) == 0 or not os.path.exists(clustalw_out) and os.path.getsize(clustalw_out) == 0:

        exit_on_error('sys_error', "compare_atom_seqres_or_msa : one of clustalw outputs were not create; %s or %s\n" %(pairwise_aln, clustalw_out))
"""



def get_seqres_atom_seq(PDB_Obj, chain, pdb_file_name):

    # extract the sequences from the pdbParser

    seqres = ""
    atom = ""

    for chainid in PDB_Obj.get_SEQRES_chains():

        if (chainid == " " and chain == "NONE") or chain == chainid:

            seqres = PDB_Obj.get_SEQRES(chainid)
            break

    for chainid in PDB_Obj.get_ATOM_chains():

        if (chainid == " " and chain == "NONE") or chain == chainid:

            atom = PDB_Obj.get_ATOM(chainid)
            break

    if seqres == "" and atom =="":

        exit_on_error('user_error', "The protein sequence for chain '%s' was not found in SEQRES nor ATOM fields in the <a href=\"%s\">PDB file</a>. For more details please see <a href = \"%s\"\>PDB File Format Documentation</a>" %(chain, pdb_file_name, CONSURF_CONSTANTS.PDB_FILE_FORMAT_CONTENTS_GUIDE))

    # output a message in case there is no seqres relevant sequence, but there are other chains.

    if seqres == "" and PDB_Obj.get_SEQRES_chains():

        all_chains = ""

        for chainid in PDB_Obj.get_SEQRES_chains():

            if chainid == " ":

                chainid = "NONE"

            all_chains += chainid + ", "

        all_chains = all_chains[:-2]

        if chain == "NONE":

            exit_on_error('user_error', "The chain column in SEQRES field is not empty, but contains the chains: %s. Please check the <a href = \"%s\">PDB file</a>, or run ConSurf again with a specific chain identifier." %(all_chains, pdb_file_name))

        else:

            exit_on_error('user_error', "Chain \"%s\" does not exist in the SEQRES field of the <a href = \"%s\">PDB file</a>.<br>Please check your file or run ConSurf again with a different chain. If there is no chain identifier, choose \"NONE\" as your Chain Identifier." %(chain, pdb_file_name))

    return [seqres, atom]





