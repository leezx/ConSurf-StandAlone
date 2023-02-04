import re
import GENERAL_CONSTANTS

def print_msa_colors_FASTA_clustalwLike(grades, MSA, SeqName, Out_file, CSS_File, Header, AA_DNA):

    # Print the results: the colored sequence of the msa acording to the query sequence

    MSA_To_RefSeq_Pos = [] # array that holds for each position on the MSA the Refernce seq position
    blockSize = 50

    if AA_DNA == "AA":

        unknownChar = 'X'

    else:

        unknownChar = 'N'

    try:

        OUT = open(Out_file, 'w')

    except:

        return("could not open the file " + Out_file + " for writing.")

    OUT.write("<!DOCTYPE html>\n")
    OUT.write("<html lang=\"en\">\n")
    OUT.write("<head>\n")
    OUT.write("<meta charset=\"utf-8\">\n")
    OUT.write("<meta http-equiv=\"X-UA-Compatible\" content=\"IE=edge\">\n")
    OUT.write("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n")
    OUT.write("<link rel=\"stylesheet\" type=\"text/css\" href=\"%s\">\n" %CSS_File)
    OUT.write("<title>%s</title>\n" %Header)
    OUT.write("</head>\n")
    OUT.write("<H1 align=center><u>%s</u></H1>\n\n" %Header) # MSA color-coded by GAIN probability

    [MSA_Hash, Seq_Names_In_Order] = ReadMSA(MSA)
    NumOfBlocks = int(len(MSA_Hash[SeqName]) / blockSize) + 1

    seq_pos = 0
    ind = 0
    for block in range(0, NumOfBlocks):

        MSA_Pos = blockSize * block
        OUT.write("<table>\n")
        seqNum = 0
        seq_pos += ind
        for name in Seq_Names_In_Order:

            seqNum += 1
            OUT.write("<tr>\n")
            if name == SeqName:

                OUT.write("<td class=\"Seq_Name\"><u><b>%d %s</u></b></td>\n" %(seqNum, name))

            else:

                OUT.write("<td class=\"Seq_Name\">%d %s</td>\n" %(seqNum, name))

            ind = 0
            for (pos, char) in list(zip(MSA_Hash[SeqName], MSA_Hash[name]))[MSA_Pos : (block + 1) * blockSize]:

                if pos != "-" and pos != unknownChar:

                    ScoreClass = "Score" + str(grades[seq_pos + ind]['COLOR'])
                    if grades[seq_pos + ind]['ISD'] == 1:

                        ScoreClass = "Score_ISD"

                    OUT.write("<td class=\"%s\">%s</td>" %(ScoreClass, char))
                    ind += 1

                else:

                    OUT.write("<td class=\"white\">%s</td>" %char)

            OUT.write("</tr>\n")

        OUT.write("</table><br><br>\n")

    # print the color scale
    OUT.write("<table style = 'table-layout: auto;margin-left: 0em;margin-right: 0em;padding:1px 1px 1px 1px; margin:1px 1px 1px 1px; border-collapse: collapse;' border=0 cols=1 width=310>\n<tr><td align=center>\n<font face='Courier New' color='black' size=+1><center>\n<tr>")
    for i in range(1,10):

        OUT.write("<td class=\"%s\">%d</td>\n" %("Score" + str(i), i))

    OUT.write("</tr></font></center>\n<center><table style = 'table-layout: auto;margin-left:0em;margin-right: 0em;padding:1px 1px 1px 1px; margin:1px 1px 1px 1px; border-collapse: collapse;' border=0 cols=3 width=310>\n<tr>\n<td align=left><td align=left><b>Variable</b></td><td></td><td align=center><b>Average</b></td><td></td>\n<td align=right><b>Conserved</b></td>\n</tr><tr></tr><tr></tr>\n</table></center>\n")
    OUT.write("<table><tr><b><td class=\"Score_ISD\">X</td><td class=\"white\"> - Insufficient data - the calculation for this site was performed on less than 10% of the sequences.</b><br></td></tr></table>\n")
    OUT.write("</body>\n</table>\n")
    OUT.close()

def ReadMSA(msa):

    Seq_Names_In_Order = [] # array to hold sequences names in order
    MSA_Hash = {} # hash to hold sequnces
    Seq = ""
    Seq_Name = ""

    try:

        MSA = open(msa, 'r')

    except:

        return("ReadMSA: Can't read the MSA: " + msa, "")

    line = MSA.readline()
    while line != "":

        line = line.rstrip()
        match = re.match(r'^>(.*)', line)
        if match:

            if Seq != "":

                MSA_Hash[Seq_Name] = Seq
                Seq = ""
                Seq_Name = ""

            Seq_Name = match.group(1)
            Seq_Names_In_Order.append(Seq_Name)

        else:

            Seq += line

        line = MSA.readline()

    MSA.close()
    MSA_Hash[Seq_Name] = Seq # last sequence

    return(MSA_Hash, Seq_Names_In_Order)

def Find_Good_Templates(query, Blast_vs_PDB, Min_Overlap_Percent, Min_ID, Selected_Templates_Nems_ref, Selected_Templates_Details_ref):

    Min_ID = Min_ID * 100
    templates = {} # hash of templates names. each sequence name (unique) points to array of hashes.
    query_seq_length = len(query)
    query_min_length = query_seq_length * Min_Overlap_Percent # min length of overlap with the query





def ConSeq_HTML_Output(output, ref_Solv_Acc_Pred, Out, b_e_prediction_method = "neural-network algorithm", Just_Seq_Flag = False):

    # Print the results: the colored sequence and the B/E information

    consurf_html_colors = {1 : "#0A7D82", 2 : "#4BAFBE", 3 : "#A5DCE6", 4 : "#D7F0F0", 5 : "#FFFFFF", 6 : "#FAEBF5", 7 : "#FAC8DC", 8 : "#F07DAA", 9 : "#A0285F", 'ISD' : "#FFFF96"}

    try:

        COLORS = open(Out, 'w')

    except:

        return("ConSeq_gradesPE_and_Outputs::ConSeq_HTML_Output: Can\'t open the file " + Out + "for writing.")

    COLORS.write("<html>\n<title>ConSeq Results</title>\n")
    COLORS.write("<body bgcolor='white'>\n")
    COLORS.write("<H1 align=center><u>ConSeq Results</u></H1>\n\n")
    COLORS.write("\n<table border=0 width=750>\n")
    COLORS.write("<tr><td>\n")

    # print the colored sequence

    count = 1
    letter_str = ""
    pred_str = ""
    func_str = ""

    for elem in output:

        # print the counter above the beginning of each 10 characters
        if count % 50 == 1:

            count_num = count
            while count_num < count + 50:

                if count_num <= len(output):

                    space_num = 11 - len(str(count_num))
                    spaces = ""
                    for i in range(0, space_num):

                        spaces += "&nbsp;"

                    COLORS.write("<font face='Courier New' color='black' size=+1>" + str(count_num) + spaces + "</font>")

                count_num += 10

            COLORS.write("<br>\n")

        # print the colored letters and 'e' for the exposed residues

        # after 50 characters - print newline
        if count % 50 == 0 or count == len(output):

            if elem['ISD'] == 1: # INSUFFICIENT DATA

                letter_str += "<b><font face='Courier New' color='%s' size=+1><span style='background: %s;'>%s</span></font></b><br>" %(consurf_html_colors['ISD'], consurf_html_colors[elem['COLOR']], elem['SEQ'])

            elif elem['COLOR'] == 9 or elem['COLOR'] == 1: # MOST OR LEAST CONSERVED

                letter_str += "<b><font face='Courier New' color='white' size=+1><span style='background: %s;'>%s</span></font></b><br>\n" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

            else:

                letter_str += "<b><font face='Courier New' color='black' size=+1><span style='background: %s;'>%s</span></font></b><br>\n" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

            if not Just_Seq_Flag:

                if ref_Solv_Acc_Pred[elem['POS']] == "e": # exposed

                    pred_str += "<b><font face='Courier New' color='orange' size=+1>e</font></b><br>\n"

                    if elem['COLOR'] == 9 or elem['COLOR'] == 8: # exposed & conserved = functional

                        func_str += "<b><font face='Courier New' color='red' size=+1>f</font></b>\n"

                    else: # not functional

                        func_str += "<font face='Courier New' size=+1>&nbsp;</font>\n"

                else: # burried

                    pred_str += "<b><font face='Courier New' color='#00cc00' size=+1>b</font></b><br>\n"

                    if elem['COLOR'] == 9: # burried & conserved = structural

                        func_str += "<b><font face='Courier New' color='#000099' size=+1>s</font></b>\n"

                    else: # not structural

                        func_str += "<font face='Courier New' size=+1>&nbsp;</font>\n"

            COLORS.write(letter_str)
            if not Just_Seq_Flag:

                COLORS.write(pred_str)
                COLORS.write(func_str)

            COLORS.write("</td></tr>\n")
            COLORS.write("<tr><td>\n")

            letter_str = ""
            pred_str = ""
            func_str = ""

        elif count % 10 == 0: # after 10 characters - print a space ('&nbsp;')

            if elem['ISD'] == 1:

                letter_str += "<b><font face='Courier New' color='%s' size=+1><span style='background: %s;'>%s</span> </font></b>\n" %(consurf_html_colors['ISD'], consurf_html_colors[elem['COLOR']], elem['SEQ'])

            elif elem['COLOR'] == 9 or elem['COLOR'] == 1: # MOST OR LEAST CONSERVED

                letter_str += "<b><font face='Courier New' color='white' size=+1><span style='background: %s;'>%s</span> </font></b>\n" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

            else:

                letter_str += "<b><font face='Courier New' color='black' size=+1><span style='background: %s;'>%s</span> </font></b>\n" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

            if not Just_Seq_Flag: # exposed

                if ref_Solv_Acc_Pred[elem['POS']] == "e":

                    pred_str += "<b><font face='Courier New' color='orange' size=+1>e&nbsp;</font></b>"

                    if elem['COLOR'] == 9 or elem['COLOR'] == 8: # exposed & conserved = functional

                        func_str += "<b><font face='Courier New' color='red' size=+1>f&nbsp;</font></b>"

                    else: # not functional

                        func_str += "<font face='Courier New' size=+1>&nbsp;&nbsp;</font>"

                else: # burried

                    pred_str += "<b><font face='Courier New' color='#00cc00' size=+1>b&nbsp;</font></b>"

                    if elem['COLOR'] == 9: # burried & conserved = structural

                        func_str += "<b><font face='Courier New' color='#000099' size=+1>s&nbsp;</font></b>"

                    else: # not structural

                        func_str += "<font face='Courier New' size=+1>&nbsp;&nbsp;</font>"

        else:

            if elem['ISD'] == 1:

                letter_str += "<b><font face='Courier New' color='%s' size=+1><span style='background: %s;'>%s</span> </font></b>" %(consurf_html_colors['ISD'], consurf_html_colors[elem['COLOR']], elem['SEQ'])

            elif elem['COLOR'] == 9 or elem['COLOR'] == 1: # MOST OR LEAST CONSERVED

                letter_str += "<b><font face='Courier New' color='white' size=+1><span style='background: %s;'>%s</span> </font></b>" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

            else:

                letter_str += "<b><font face='Courier New' color='black' size=+1><span style='background: %s;'>%s</span> </font></b>" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

            if not Just_Seq_Flag: # exposed

                if ref_Solv_Acc_Pred[elem['POS']] == "e":

                    pred_str += "<b><font face='Courier New' color='orange' size=+1>e</font></b>"

                    if elem['COLOR'] == 9 or elem['COLOR'] == 8: # exposed & conserved = functional

                        func_str += "<b><font face='Courier New' color='red' size=+1>f</font></b>"

                    else: # not functional

                        func_str += "<font face='Courier New' size=+1>&nbsp;</font>"

                else: # burried

                    pred_str += "<b><font face='Courier New' color='#00cc00' size=+1>b</font></b>"

                    if elem['COLOR'] == 9: # burried & conserved = structural

                        func_str += "<b><font face='Courier New' color='#000099' size=+1>s</font></b>"

                    else: # not structural

                        func_str += "<font face='Courier New' size=+1>&nbsp;</font>"

        count += 1

    COLORS.write("</td></tr>\n</table><br>\n")
    # print the color scale
    COLORS.write("\n<br><b><u>Legend:</u><br><br>\nThe conservation scale:</b><br>\n<table border=0 cols=1 width=310>\n<tr><td align=center>\n<font face='Courier New' color='black' size=+1><center>\n")
    for i in range(1, 10):

        if i == 9:

            COLORS.write("<font color='white'><span style='background: %s;'>&nbsp;%d&nbsp;</span></font>" %(consurf_html_colors[i], i))

        else:

            COLORS.write("<span style='background: %s;'>&nbsp;%d&nbsp;</span>" %(consurf_html_colors[i], i))

    COLORS.write("</font></center>\n<center><table border=0 cols=3 width=310>\n<tr>\n<td align=left><b>Variable</b></td>\n<td align=center><b>Average</b></td>\n<td align=right><b>Conserved</b></td>\n</tr>\n</table></center>\n</td>\n</tr>\n</table>\n")
    if not Just_Seq_Flag:

        COLORS.write("<b><font face='Courier New' color='orange' size=+1>e</font> - An exposed residue according to the %s.</b><br>\n" %b_e_prediction_method)
        COLORS.write("<b><font face='Courier New' color='#00cc00' size=+1>b</font> - A buried residue according to the %s.</b><br>\n" %b_e_prediction_method)
        COLORS.write("<b><font face='Courier New' color='red' size=+1>f</font> - A predicted functional residue (highly conserved and exposed).</b><br>\n")
        COLORS.write("<b><font face='Courier New' color='#000099' size=+1>s</font> - A predicted structural residue (highly conserved and buried).</b><br>\n")

    COLORS.write("<b><font face='Courier New' color='%s' size=+1><span style='background: %s;'>X</span></font> - Insufficient data - the calculation for this site was performed on less than 10%% of the sequences.</font></b><br>\n" %(consurf_html_colors['ISD'], consurf_html_colors[2]))
    COLORS.write("</body>\n</html>\n")
    COLORS.close()

    return("ok")

def ConSeq_PDF_Output(Grades_PE, Out, Protein_Length, b_e_prediction_method):

    ConSeq_PDF_Lib = GENERAL_CONSTANTS.CONSURF_HTML_DIR + "ConSurf_PDF/ConSurf_PDF_Class.php"

    try:

        PDF = open(Out, 'w')

    except:

        return("ConSeq_gradesPE_and_Outputs::ConSeq_PDF_Output: Can\'t open the file " + Out + " for writing.")

    ConSurf_Grades = {}
    ans = read_ConSeq_Grades_PE(Grades_PE, ConSurf_Grades)
    if ans[0] != "ok":

        return("ConSeq_gradesPE_and_Outputs.ConSeq_PDF_Output: \'" + ans + "\'")

    IS_THERE_FUNCT_RES = ans[1]
    IS_THERE_STRUCT_RES = ans[2]
    IS_THERE_INSUFFICIENT_DATA = ans[3]

    PDF.write("""<?php
    require("%s");
    \$pdf=new ConSurf_PDF();
    \$pdf->setCBS('no');
    \$pdf->AddPage();
    \$pdf->SetFont('Times','BU',20);
    \$pdf->Cell(0,0,'ConSurf Results',0,0,C);
    \$pdf->SetY(\$pdf->GetY()+10);""" %ConSeq_PDF_Lib)

    maxPosPerPage = 600 # 12 rows
    for pos in range(1, Protein_Length):

        if ((pos - 1)  % maxPosPerPage == 0) and pos != 1:

            PDF.write("\$pdf->AddPage();")

        if (pos - 1) % 50 == 0:

            PDF.write("\$pdf->Ln();\n\$pdf->Ln();\n\$pdf->Ln();\n\$pdf->Ln();\n")
            PDF.write("\$Rows_Pos=\$pdf->GetY();\n")

        elif (pos - 1) % 10 == 0:

            PDF.write("\$pdf->Print_ForegroundColor\('','B',10,0,4\);\n")

        if (pos - 1) % 10 == 0:

            PDF.write("\$pdf->Print_4_Lines_Element\(\$Rows_Pos,%s,%s,%d,%s,10,%d,'%s'\);\n" %(ConSurf_Grades[pos]['COLOR'], ConSurf_Grades[pos]['AA'], pos, ConSurf_Grades[pos]['B_E'], ConSurf_Grades[pos]['Insufficient_Data'], ConSurf_Grades[pos]['Struct_Funct']))

        else:

            PDF.write("\$pdf->Print_4_Lines_Element\(\$Rows_Pos,%s,%s,'',%s,10,%d,'%s'\);\n" %(ConSurf_Grades[pos]['COLOR'], ConSurf_Grades[pos]['AA'], ConSurf_Grades[pos]['B_E'], ConSurf_Grades[pos]['Insufficient_Data'], ConSurf_Grades[pos]['Struct_Funct']))

    PDF.write("\$pdf->Ln();\n")
    PDF.write("\$pdf->Ln();\n")
    PDF.write("\$pdf->Ln();\n")
    PDF.write("\$pdf->Ln();\n")
    PDF.write("\$pdf->Print_NEW_Legend(%d, %d, %d, '%s');\n" %(IS_THERE_FUNCT_RES, IS_THERE_STRUCT_RES, IS_THERE_INSUFFICIENT_DATA, b_e_prediction_method))
    PDF.write("\$pdf->Output();\n")
    PDF.write("?>\n")
    PDF.write(" </span><font color='white'><span style='background: #A02560;'>&nbsp;&nbsp;</span></font></font></center>\n")
    PDF.close()

def read_ConSeq_Grades_PE(Grades_PE, ref_GradesPE_Hash):

    IS_THERE_FUNCT_RES = 0
    IS_THERE_STRUCT_RES = 0
    IS_THERE_INSUFFICIENT_DATA = 0

    try:

        GRADES = open(Grades_PE, 'r')

    except:

        return("ConSeq_gradesPE_and_Outputs.ConSeq_PDF_Output: Can\'t open the GradesPE file \'" + Grades_PE + "\' for reading.")

    line = GRADES.readline()
    while line != "":

        line = line.rstrip()
        if len(line) >= 59:

            pos = (line[:4]).strip()
            AA = (line[5:9]).strip()
            Score = (line[9:16]).strip()
            Color = (line[20:22]).strip()
            B_E = (line[50:52]).strip()
            S_F = (line[59:60]).strip()

            if pos.isdigit():

                pos = int(pos)
                ref_GradesPE_Hash[pos] = {}
                ref_GradesPE_Hash[pos]['AA'] = AA
                ref_GradesPE_Hash[pos]['SCORE'] = float(Score)
                ref_GradesPE_Hash[pos]['COLOR'] = Color
                ref_GradesPE_Hash[pos]['B_E'] = B_E
                ref_GradesPE_Hash[pos]['Struct_Funct'] = S_F

                if S_F == "s":

                    IS_THERE_STRUCT_RES = 1

                elif S_F == "f":

                    IS_THERE_FUNCT_RES = 1

                if len(Color) > 1 and (Color[0]).isdigit() and Color[1] == "*":

                    ref_GradesPE_Hash[pos]['Insufficient_Data'] = 1
                    ref_GradesPE_Hash[pos]['COLOR'] = Color[0]
                    IS_THERE_INSUFFICIENT_DATA = 1

                else:

                    ref_GradesPE_Hash[pos]['Insufficient_Data'] = 0

        line = GRADES.readline()

    return("ok", IS_THERE_FUNCT_RES, IS_THERE_STRUCT_RES, IS_THERE_INSUFFICIENT_DATA)






def consurf_HTML_Output(output, Out, cbs):

    # Print the results: the colored sequence and the B/E information
    if cbs:

        consurf_html_colors = {1 : "#0F5A23", 2 : "#5AAF5F", 3 : "#A5DCA0", 4 : "#D7F0D2", 5 : "#FFFFFF", 6 : "#E6D2E6", 7 : "#C3A5CD", 8 : "#9B6EAA", 9 : "#782882", 'ISD' : "#FFFF96"}

    else:

        consurf_html_colors = {1 : "#0A7D82", 2 : "#4BAFBE", 3 : "#A5DCE6", 4 : "#D7F0F0", 5 : "#FFFFFF", 6 : "#FAEBF5", 7 : "#FAC8DC", 8 : "#F07DAA", 9 : "#A0285F", 'ISD' : "#FFFF96"}

    try:

        COLORS = open(Out, 'w')

    except:

        return("ConSeq_gradesPE_and_Outputs::consurf_HTML_Output: Can\'t open the file " + Out + "for writing.")

    COLORS.write("<html>\n<title>ConSurf Results</title>\n")
    COLORS.write("<head>\n<style>\nb { float: left;}\n</style>\n</head>\n")
    COLORS.write("<body bgcolor='white'>\n")
    COLORS.write("\n<table border=0 width=750>\n")
    COLORS.write("<tr><td>\n")

    # print the colored sequence

    count = 1
    letter_str = ""
    pred_str = ""
    func_str = ""

    for elem in output:

        # print the counter above the beginning of each 10 characters
        if count % 50 == 1:

            count_num = count
            while count_num < count + 50:

                if count_num <= len(output):

                    space_num = 11 - len(str(count_num))
                    spaces = ""
                    for i in range(0, space_num):

                        spaces += "&nbsp;"

                    COLORS.write("<font face='Courier New' color='black' size=+1>" + str(count_num) + spaces + "</font>")

                count_num += 10

            COLORS.write("<br>\n")

        # print the colored letters and 'e' for the exposed residues

        # after 50 characters - print newline
        if count % 50 == 0 or count == len(output):

            if elem['ISD'] == 1: # INSUFFICIENT DATA

                letter_str += "<b><font face='Courier New' color='black' size=+1><span style='background: %s;'>%s</span></font></b><br>" %(consurf_html_colors['ISD'], elem['SEQ'])

            elif elem['COLOR'] == 9 or elem['COLOR'] == 1: # MOST OR LEAST CONSERVED

                letter_str += "<b><font face='Courier New' color='white' size=+1><span style='background: %s;'>%s</span></font></b><br>\n" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

            else:

                letter_str += "<b><font face='Courier New' color='black' size=+1><span style='background: %s;'>%s</span></font></b><br>\n" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

            COLORS.write(letter_str)
            COLORS.write("</td></tr>\n")
            COLORS.write("<tr><td>\n")

            letter_str = ""
            pred_str = ""
            func_str = ""

        elif count % 10 == 1 and count % 50 != 1: # after 10 characters - print a space ('&nbsp;')

            letter_str += "<b><font face='Courier New' color='black' size=+1>&nbsp;</font></b>"

            if elem['ISD'] == 1:

                letter_str += "<b><font face='Courier New' color='black' size=+1><span style='background: %s;'>%s</span> </font></b>\n" %(consurf_html_colors['ISD'], elem['SEQ'])

            elif elem['COLOR'] == 9 or elem['COLOR'] == 1: # MOST OR LEAST CONSERVED

                letter_str += "<b><font face='Courier New' color='white' size=+1><span style='background: %s;'>%s</span> </font></b>\n" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

            else:

                letter_str += "<b><font face='Courier New' color='black' size=+1><span style='background: %s;'>%s</span> </font></b>\n" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

        else:

            if elem['ISD'] == 1:

                letter_str += "<b><font face='Courier New' color='black' size=+1><span style='background: %s;'>%s</span> </font></b>\n" %(consurf_html_colors['ISD'], elem['SEQ'])

            elif elem['COLOR'] == 9 or elem['COLOR'] == 1: # MOST OR LEAST CONSERVED

                letter_str += "<b><font face='Courier New' color='white' size=+1><span style='background: %s;'>%s</span> </font></b>\n" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

            else:

                letter_str += "<b><font face='Courier New' color='black' size=+1><span style='background: %s;'>%s</span> </font></b>\n" %(consurf_html_colors[elem['COLOR']], elem['SEQ'])

        count += 1

    COLORS.write("</td></tr>\n</table><br>\n")
    COLORS.write("</body>\n</html>\n")
    COLORS.close()








