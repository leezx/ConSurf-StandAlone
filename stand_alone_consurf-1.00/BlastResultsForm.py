import re
from Bio import SearchIO

def Create_HTML_file_for_blast(choose_file, parsed_output, raw_output, query, job_number, database):
    
    max_seqs = 2000

    try:

        OUT = open(choose_file, 'w')

    except:

        return("error", "BlastResultsForm.Create_HTML_file_for_hmmer : can't open file " + choose_file + " for writing\n")
    
    try:

        PARSED_OUT = open(parsed_output, 'w')

    except:

        return("error", "BlastResultsForm.Create_HTML_file_for_hmmer : can't open file " + parsed_output + " for writing\n")

    try:

        RAW_OUT = open(raw_output, 'r')

    except:

        return("error", "BlastResultsForm.Create_HTML_file_for_hmmer : can't open file " + raw_output + " for reading\n")

    alignments_details = []
 
    searchio = SearchIO.parse(RAW_OUT, "blast-xml")
    last_blast_round = ""
    for qresult in searchio:

        last_blast_round = qresult

    number_of_sequnces = len(last_blast_round)

    seq_count = 0
    for hit in last_blast_round[:max_seqs]:

        if seq_count >= max_seqs:
            
            break
       
        description = hit.description
        name = hit.id
        hsp_count = 1
        for hsp in hit:
            
            if seq_count >= max_seqs:
                
                break
                
            [s_beg, s_end] = hsp.hit_range
            s_beg += 1

            if len(hit) > 1:
                
                name_domain = "%s domain %d" %(name, hsp_count)
                
            else:
                
                name_domain = name
                
            seq = str(hsp.hit.seq)
            length = len(seq)
            
            E_val = str(hsp.evalue)
            E_val = re.sub(r',', "", E_val)
            match = re.match(r'^e', E_val)
            if match:

                E_val = "1" + E_val
                
            score = hsp.bitscore

            
            new_domain = "\nseqgment start: %d\nsegment end: %d\nquery seq = %s\ndatabase seq = %s\n\n   query seq: %s\n              %s\ndatabase seq: %s\n" %(hsp.hit_start, hsp.hit_end, query, name, hsp.query.seq, hsp.aln_annotation['similarity'], hsp.hit.seq)
            title = "\tScore = %.3f, Expect = %s" %(score, E_val)  
            if getattr(hsp, 'gap_num', False):
                
                title += ", gaps = %d/%d (%d%%)" %(hsp.gap_num, length, int((float(hsp.gap_num) / length)
                                                                          * 100))
            else:
                
                title += ", gaps = 0" 
                
            if getattr(hsp, 'ident_num', False):
                
                title += ", Identities = %d/%d (%d%%)" %(hsp.ident_num, length, int((float(hsp.ident_num)
                                                                                   / length) * 100))
            else:
                  
                title += ", Identities = 0"  
                  
            if getattr(hsp, 'pos_num', False):
                
                title += ", Positives = %d/%d (%d%%)" %(hsp.pos_num, length, int((float(hsp.pos_num) / length) * 100))
                
            else:
                    
                title += ", Positives = 0" 
                    
            new_domain = title + "\n" + new_domain + "\n\n\n"
            
            PARSED_OUT.write("<pre <?=$hidden[%d];?>>\n" %seq_count)
            PARSED_OUT.write(new_domain)
            PARSED_OUT.write("\n</pre>\n\n")
            
            OUT.write(write_title(name_domain, description, seq_count, database_address(name, database), E_val, job_number))
            
            details = {}
            details['name'] = "%s_%d_%d" %(name, s_beg, s_end)
            details['seq'] = re.sub(r'-', "", seq)
            alignments_details.append(details)

            seq_count += 1
            hsp_count +=1

    RAW_OUT.close()
    OUT.close()
    PARSED_OUT.close()
    return("OK", alignments_details, number_of_sequnces)
   
    
    
    
    
def Create_HTML_file_for_hmmer(choose_file, parsed_output, raw_output, query, job_number, database):

    max_seqs = 2000

    try:

        OUT = open(choose_file, 'w')

    except:

        return("error", "BlastResultsForm.Create_HTML_file_for_hmmer : can't open file " + choose_file + " for writing\n")
    
    try:

        PARSED_OUT = open(parsed_output, 'w')

    except:

        return("error", "BlastResultsForm.Create_HTML_file_for_hmmer : can't open file " + parsed_output + " for writing\n")

    try:

        RAW_OUT = open(raw_output, 'r')

    except:

        return("error", "BlastResultsForm.Create_HTML_file_for_hmmer : can't open file " + raw_output + " for reading\n")

    last_round = ((RAW_OUT.read()).split("@@ Round:"))[-1]
    RAW_OUT.close()
    number_of_sequnces = (re.search(r'@@ New targets included:\s+(\S+)', last_round)).group(1)
    last_round = (last_round.split("\n\nInternal pipeline statistics summary:"))[0] # we delete the summary:
    alignments = (last_round.split(">>"))[1 : max_seqs + 1] # we only take the first 2000 sequences.
    alignments_details = []
    
    i = 0
    for alignment in alignments:

        if i >= max_seqs:
            
            break
        
        domains = alignment.split("  == domain")
        match = re.match(r'^\s*(\S+)(\s.+)', domains[0])
        #name = match.group(1)
        description = match.group(2)
        #RepID = (re.search(r'RepID=(\S+)', domains[0])).group(1)
        for domain in domains[1:]:
            
            if i >= max_seqs:
                
                break
            

            lines = domain.split("\n")
            query_line = 0
            for line in lines:
                
                if query in line:
                    
                    break
                
                else:
                    
                    query_line += 1
                    
            [name, seq_start, seq, seq_end] = lines[query_line + 2].split()
            if len(domains[1:]) > 1:
                
                number = re.search(r'^\s+(\S+)', lines[0]).group(1)
                #name_domain = name + "_" + seq_start + "_" + seq_end # no space
                name_domain_space = name + " domain " + number # space
                
            else:
                
                #name_domain = name
                name_domain_space = name

            score = re.search(r'score:\s+(\S+\s+bits)', lines[0]).group(1)
            Evalue = re.search(r'E-value:\s+(\S+)', lines[0]).group(1)
            #seq = (lines[3].split())[2]
            
            new_domain = "\n" + lines[query_line] + "\n" + lines[query_line + 1] + "\n" + lines[query_line + 2]
            title = "\tScore = " + score + ", Expect = " + Evalue + get_details(lines[query_line], lines[query_line + 1], lines[query_line + 2])
            new_domain = title + "\n" + new_domain + "\n\n\n"
            
            PARSED_OUT.write("<pre <?=$hidden[%d];?>>\n" %i)
            PARSED_OUT.write(new_domain)
            PARSED_OUT.write("\n</pre>\n\n")
            
            OUT.write(write_title(name_domain_space, description, i, database_address(name, database), Evalue, job_number))
            
            details = {}
            details['name'] = name + "_" + seq_start + "_" + seq_end
            details['seq'] = re.sub(r'-', "", seq)
            alignments_details.append(details)

            i += 1

    """
    
    i = 1
    for alignment in alignments_details:
        
        OUT.write(write_title(alignment['name'], alignment['description'], i, "", 1, job_number))
        i += 1
        
    
    i = 1
    for alignment in alignments_details:
        
        PARSED_OUT.write("<pre <?=$hidden[%d];?>>\n" %(i -1))
        #PARSED_OUT.write(str(i) + " " + write_name(alignment['name'], alignment['description'], "") + "\n")
        for domain in alignment['domains']:
            
            PARSED_OUT.write(domain + "\n")
            
        PARSED_OUT.write("</pre>\n\n")
        i += 1
    """
    OUT.close()
    PARSED_OUT.close()
    return("OK", alignments_details, int(number_of_sequnces))

def get_details(first_line, second_line, third_line):

    query_seq = (first_line.split())[2]
    db_seq = (third_line.split())[2]

    gaps = 0
    positives = 0
    identities = 0
    length = len(query_seq)

    for char in query_seq:

        if char == "-":

            gaps += 1

    for char in db_seq:

        if char == "-":

            gaps += 1

    for char in second_line:

        if char == "+":

            positives += 1

        if char.isalpha():

            identities += 1
            
    details = ""
    
    if gaps != 0:
        
        details += ", Gaps = %d/%d (%d%%)" %(gaps, length, int((float(gaps) / length) * 100))
        
    else:
        
        details += ", Gaps = 0"
        
    if positives != 0:
        
        details += ", Positives = %d/%d (%d%%)" %(positives, length, int((float(positives) / length) * 100))
        
    else:
        
        details += ", Positives = 0"
        
    if identities != 0:
        
        details += ", Identities = %d/%d (%d%%)" %(identities, length, int((float(identities) / length) * 100))
        
    else:
        
        details += ", Identities = 0"

    return details


def write_title(name, description, number, database, evalue, job_number):
    
    title = "<br><br><div class=\"radio_btns\"><label class=\"container\">Add Sequence Number %d<input type = \"checkbox\" value = \"%d\" name = \"MarkedSeq[]\" id = \"%d.evalue\" onchange = \"document.getElementById('%d.align').checked=document.getElementById('%d.evalue').checked\"><span class=\"checkmark\"></span></label></div><br>" %(number + 1, number, number, number, number)
    #title += "Sequence Name: %s<br><br>" %name
    title += "<a class=\"ChooseSubmit ChooseSubmit1\" href=\"%s\"target=\"_blank\">%s</a><br><br>" %(database, name)
    title += "<a class=\"ChooseSubmit ChooseSubmit1\" href=\"hide_show.php?name=%d&number=%s\"target=\"_blank\">Show Alignment</a><br><br>\n" %(number, job_number)
    title += "<ul>Description: %s</ul><br>" %description
    title += "<ul>E-value: %s</ul><br><br>" %evalue
    
    return title
    
def database_address(name, database): 
    
    match = re.match(r'(UniRef90_[A-Za-z0-9]+)', name)
    if match:

        return(database + "uniref/" + match.group(1))

    match = re.match(r'sp\|([A-Za-z0-9]+)', name)
    if match:

        return(database + "uniprot/" + match.group(1))
	
    match = re.match(r'tr\|([A-Za-z0-9]+)', name)
    if match:

        return(database + "uniprot/" + match.group(1))
                
    match = re.match(r'gi\|([A-Za-z0-9]+)', name)
    if match:

        return(database + match.group(1))	
				
    match = re.match(r'\S+\|([A-Za-z0-9]+)', name)
    if match:

        return(database + "uniprot/" + match.group(1))
					
def write_name(name, description, database):
    
    return "<input type = \"checkbox\" id = \"%s.align\" onchange=\"document.getElementById('%s.evalue').checked=document.getElementById('%s.align').checked\"> ><b><a href=\"%s=%s\">%s</a></b>%s</br>" %(name, name, name, database, name, name, description)






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







def MakeForm(Blast_HTML_In, Round, Blast_Form_Out, Selected_Seq_File, DB = "", Multi_Select = "yes", Redirect_Output = "", search_type = "", num_of_hits = "NA"):

    max_hits_for_display=2000
    try:

        BLAST_FORM = open(Blast_Form_Out, 'w')

    except:

        return("BlastResultsForm.MakeForm : Could Not Open The Blast Form Output " + Blast_Form_Out + " For writing\n")

    try:

        BLAST_HTML_RESULTS = open(Blast_HTML_In, 'r')

    except:

        return("BlastResultsForm.MakeForm : Could Not Open The Blast HTML Output " + Blast_HTML_In + " For reading\n")

    BLAST_FORM.write("""<\?php
    \$selected_seq = \$_POST[\"MarkedSeq\"];
    if (!isset(\$_POST['submit'])) { // if page is not submitted to itself echo the form
    \?>

    <div id="Select_MSA_Seq">
    <form name="Select_Blast_Hit" id="Select_Blast_Hit" method="post">
    <H1 align=center><font color='red'> Please choose which sequences you want to use for ConSurf calculation</font></h1>
    <H3 align=center><font color='black'> Please noticed that the query sequence is automatically included in the analysis without choosing it</font></h3>""")

    if num_of_hits != "NA" and num_of_hits > max_hits_for_display:

        BLAST_FORM.write("<H3 align=center><font color='red'>Please note: " + num_of_hits + " hits were found for your query sequence, ConSurf display here only the first " + max_hits_for_display + " hits</font></h3>\n")

    BLAST_FORM.write("<FONT FACE= \"Courier New\">\n")

    if Multi_Select == "yes":

        BLAST_FORM.write("""Select the first <INPUT TYPE="text" NAME="FirstHomologs" VALUE="50" id = "FirstHomologs" size=6> sequences (select <input type="radio" name="SelectJump" value="1" checked>all or only every <input type="radio" name="SelectJump" value="2">2nd, <input type="radio" name="SelectJump" value="3">3rd, <input type="radio" name="SelectJump" value="4">4th sequence)<INPUT TYPE="button" NAME="button1" Value="Update Selection" onClick="clearAll();SelectFirstSeqWithJump(document.Select_Blast_Hit.FirstHomologs.value,getRadioVal('SelectJump'))">


        <input type="button" name="CheckAll" value="Check All(MAX:500)" onClick="javascript:checkAll()"> <input type="button" name="UnCheckAll" value="Uncheck All" onClick="javascript:clearAll()">
        <a href="#submit">Go to submit</a>""")

    if search_type == "HMMER":

        counter_table = 1
        counter_align = 1
        SeqID = ""
        line = BLAST_HTML_RESULTS.readline()
        while line != "":

            match1 = re.match(r'^(\<tr\>\<td\>)(\<a href\=\"https:\/\/www\.ncbi\.nlm\.nih\.gov\/)(nucleotide|protein)(\?term\=([A-Za-z0-9_\-\.]+)\")(>.*)(\<td\>0\<\/td\>)(.*)', line)
            if not match1:

                match2 = re.match(r'Reformatted HTML of HMMSEARCH Search Report', line)
                if match2: # junk...

                    line =  BLAST_HTML_RESULTS.readline()
                    continue

                match3 = re.match(r'\<table border\=0\>', line)
                if match3:

                    BLAST_FORM.write("<table border=0>")
                    line =  BLAST_HTML_RESULTS.readline()
                    continue

                match4 = re.match(r'Sequences producing significant alignments', line)
                if match4: # table

                    # add column to tha table
                    BLAST_FORM.write("<tr><th></th><th></th><th>Sequences producing significant alignments:</th>")
                    line =  BLAST_HTML_RESULTS.readline()
                    continue

                #first seq (table):
                match5 = re.match(r'^(.*)(\<tr\>\<td\>)(\<a href\=\"https:\/\/www\.ncbi\.nlm\.nih\.gov\/)(nucleotide|protein)(\?term\=([A-Za-z0-9_\-\.]+)\")(>.*)(\<td\>0\<\/td\>)(.*)', line)
                if match5:

                    SeqID = match5.group(6)
                    BLAST_FORM.write("<th>E<br>value</th></tr>\n")
                    BLAST_FORM.write("<tr><td>%d</td><td> <input type = \"checkbox\" value = \"%s.HmmerHit_%d\" name = \"MarkedSeq[]\" id = \"%s.evalue\" onchange=\"document.getElementById('%s.align').checked=document.getElementById('%s.evalue').checked\"></td><td>" %(counter_table, SeqID, counter_table, SeqID, SeqID, SeqID))
                    BLAST_FORM.write(match5.group(3) + match5.group(4) + match5.group(5) + "target=\"_blank\"" + match5.group(7) + match5.group(9) + "\n")
                    counter_table += 1
                    line =  BLAST_HTML_RESULTS.readline()
                    continue

                # end of table and start align
                match6 = re.match(r'^(\>\<b\>\<a href\=\"https\:\/\/www\.ncbi\.nlm\.nih\.gov\/)(nucleotide|protein)(\?term\=([A-Za-z0-9_\-\.]+)\">)(.*)', line)
                if match6:

                    if counter_align > max_hits_for_display:

                        break

                    else:

                        SeqID = match6.group(4)
                        BLAST_FORM.write("%d <input type = \"checkbox\" id = \"%s.align\" onchange=\"document.getElementById('%s.evalue').checked=document.getElementById('%s.align').checked\"> %s" %(counter_align, SeqID, SeqID, SeqID, line))
                        counter_align += 1
                        line =  BLAST_HTML_RESULTS.readline()
                        continue

                # it is seq without align
                match7 = re.match(r'^(\<tr\>\<td\>)(\<a href\=\"https:\/\/www\.ncbi\.nlm\.nih\.gov\/)(nucleotide|protein)(\?term\=([A-Za-z0-9_\-\.]+)\">)(.*)', line)
                if match7:

                    SeqID = match7.group(5)
                    line =  BLAST_HTML_RESULTS.readline()
                    continue

                match8 = re.match(r'Length \=', line)
                if match8:

                    line =  BLAST_HTML_RESULTS.readline()
                    continue

                match9 = re.match(r'Produced by Bioperl module Bio', line)
                if match9: # end of file

                    break

                BLAST_FORM.write(line)

            else: # table

                if counter_table > max_hits_for_display:

                    line =  BLAST_HTML_RESULTS.readline()
                    continue

                else:

                    SeqID = match1.group(5)
                    BLAST_FORM.write("<tr><td>%d</td><td> <input type = \"checkbox\" value = \"%s.HmmerHit_%d\" name = \"MarkedSeq[]\" id = \"%s.evalue\" onchange=\"document.getElementById('%s.align').checked=document.getElementById('%s.evalue').checked\"></td><td>" %(counter_table, SeqID, counter_table, SeqID, SeqID, SeqID))
                    BLAST_FORM.write(match1.group(2) + match1.group(3) + match1.group(4) + "target=\"_blank\"" + match1.group(6) + match1.group(8) + "\n")
                    counter_table += 1

            line =  BLAST_HTML_RESULTS.readline()

    else: # blast

        Counter_E_Value_Section = 1
        Counter_Align_Section = 1
        Desired_Round = 0
        Header = 1 # Flag to indicate blast header
        if Round > 1: # More than one Round - Go Till the results of the right run

            line =  BLAST_HTML_RESULTS.readline()
            while line != "" and Header == 1:

                match1 = re.match(r'Results from round', line)
                if not match1:

                    BLAST_FORM.write(line)

                else:

                    match2 = re.match(r'Results from round ([0-9]+)', line)
                    if match2:

                        if match2.group(1) == Round:

                            Desired_Round = 1
                            BLAST_FORM.write(line)

                        Header = 0

                line =  BLAST_HTML_RESULTS.readline()

            BLAST_HTML_RESULTS.seek(0, 0)
            line =  BLAST_HTML_RESULTS.readline()
            while line != "" and Desired_Round == 0: # Go till the desired round

                match = re.match(r'Results from round (\d+)', line)
                if match and match.group(1) == Round:

                    Desired_Round = 1
                    BLAST_FORM.write(line)

                line =  BLAST_HTML_RESULTS.readline()

        BLAST_HTML_RESULTS.seek(0, 0)
        line =  BLAST_HTML_RESULTS.readline()
        while line != "":

            printed = 0 # Was this line processed
            SeqID = ""
            match1 = re.match(r'^([A-Za-z0-9_]+)\|(.*)', line)
            match2 = re.match(r'(^UniRef.*)', line)
            if match1 or match2:

                match3 = re.match(r'^([A-Za-z0-9]+)\|(.*)\|', line)
                if match3:

                    SeqID = match3.group(2)

                if DB == "UniProt" or DB == "SWISS-PROT":

                    match4 = re.match(r'([A-Za-z0-9]+)\|([A-Za-z0-9]+_[A-Za-z0-9]+)', line)
                    if match4:

                        SeqID = match4.group(2)

                elif DB == "CLEAN_UNIPROT":

                    match5 = re.match(r'(\S+_\S+)\|(\S+)', line)
                    if match5:

                        SeqID = match5.group(1)

                elif DB == "UNIREF90":

                    match6 = re.match(r'(\S+_\S+)', line)
                    if match6:

                        SeqID = match6.group(1)

                if Counter_E_Value_Section < 10:

                    BLAST_FORM.write(str(Counter_E_Value_Section) + "&nbsp;&nbsp;&nbsp;")

                elif Counter_E_Value_Section >= 10 and Counter_E_Value_Section < 100:

                    BLAST_FORM.write(str(Counter_E_Value_Section) + "&nbsp;&nbsp;")

                elif Counter_E_Value_Section >= 100:

                    BLAST_FORM.write(str(Counter_E_Value_Section) + " ")

                if Multi_Select == "no":

                    BLAST_FORM.write("<input type=\"radio\" ")

                else:

                    BLAST_FORM.write("<input type=\"checkbox\" ")

                BLAST_FORM.write("value=\"%s.BlastHit_%d\" name=\"MarkedSeq[]\" id=\"%s.evalue\" onchange=\"document.getElementById('%s.align').checked=document.getElementById('%s.evalue').checked\">%s" %(SeqID, Counter_E_Value_Section, SeqID, SeqID, SeqID, line))
                Counter_E_Value_Section += 1
                printed = 1

            elif DB == "CULLED_PDB" or DB == "PDB":

                match3 = re.match(r'^([A-Z0-9]{5}) [0-9]+', line)
                if match3:

                    SeqID = match3.group(1)
                    if Counter_E_Value_Section < 10:

                        BLAST_FORM.write(str(Counter_E_Value_Section) + "&nbsp;&nbsp;&nbsp;")

                    elif Counter_E_Value_Section >= 10 and Counter_E_Value_Section < 100:

                        BLAST_FORM.write(str(Counter_E_Value_Section) + "&nbsp;&nbsp;")

                    elif Counter_E_Value_Section >= 100:

                        BLAST_FORM.write(str(Counter_E_Value_Section) + " ")

                    if Multi_Select == "no":

                        BLAST_FORM.write("<input type=\"radio\"  value=\"%s\" name=\"MarkedSeq\" id=\"%s.evalue\" >%s" %(SeqID, SeqID, line))

                    else:

                        BLAST_FORM.write("<input type=\"checkbox\" value=\"%s.BlastHit_%d\" name=\"MarkedSeq[]\" id=\"%s.evalue\" onchange=\"document.getElementById('%s.align').checked=document.getElementById('%s.evalue').checked\">%s" %(SeqID, Counter_E_Value_Section, SeqID, SeqID, SeqID, line))

                    Counter_E_Value_Section += 1
                    printed = 1

            elif DB == "NT": # for blast agains nt or nr

                match3 = re.match(r'^<a href=(.*)>([A-Za-z0-9_]+)\|(.*)\|(.*)?(\|)?(.*)?<\/a>', line)
                if match3:

                    SeqID = match3.group(3)
                    if match3.group(2) == "pdb":

                        SeqID = match3.group(3) + match3.group(4)

                    if Counter_E_Value_Section < 10:

                        BLAST_FORM.write(str(Counter_E_Value_Section) + "&nbsp;&nbsp;&nbsp;")

                    elif Counter_E_Value_Section >= 10 and Counter_E_Value_Section < 100:

                        BLAST_FORM.write(str(Counter_E_Value_Section) + "&nbsp;&nbsp;")

                    elif Counter_E_Value_Section >= 100:

                        BLAST_FORM.write(str(Counter_E_Value_Section) + " ")

                    if Multi_Select == "no":

                        BLAST_FORM.write("<input type=\"radio\" ")

                    else:

                        BLAST_FORM.write("<input type=\"checkbox\" ")

                    BLAST_FORM.write("value=\"%s.BlastHit_%d\" name=\"MarkedSeq[]\" id=\"%s.evalue\" onchange=\"document.getElementById('%s.align').checked=document.getElementById('%s.evalue').checked\">%s" %(SeqID, Counter_E_Value_Section, SeqID, SeqID, SeqID, line))
                    Counter_E_Value_Section += 1
                    printed = 1

            elif DB == "NR_PROT_DB":

                match3 = re.match(r'^<a href=(.*)>([A-Za-z0-9_]+)\|(.*)\|(.*)?(\|)?(.*)?<\/a>', line)
                if match3:

                    SeqID = match3.group(3)
                    if match3.group(2) == "pdb":

                        SeqID = match3.group(3) + match3.group(4)

                    if Counter_E_Value_Section < 10:

                        BLAST_FORM.write(str(Counter_E_Value_Section) + "&nbsp;&nbsp;&nbsp;")

                    elif Counter_E_Value_Section >= 10 and Counter_E_Value_Section < 100:

                        BLAST_FORM.write(str(Counter_E_Value_Section) + "&nbsp;&nbsp;")

                    elif Counter_E_Value_Section >= 100:

                        BLAST_FORM.write(str(Counter_E_Value_Section) + " ")

                    if Multi_Select == "no":

                        BLAST_FORM.write("<input type=\"radio\" ")

                    else:

                        BLAST_FORM.write("<input type=\"checkbox\" ")

                    BLAST_FORM.write("value=\"%s.BlastHit_%d\" name=\"MarkedSeq[]\" id=\"%s.evalue\" onchange=\"document.getElementById('%s.align').checked=document.getElementById('%s.evalue').checked\">%s" %(SeqID, Counter_E_Value_Section, SeqID, SeqID, SeqID, line))
                    Counter_E_Value_Section += 1
                    printed = 1

            match3 = re.match(r'^>(<a name = [0-9]+><\/a>)(.*)', line)
            if match3:

                SeqID = match3.group(2)
                if DB == "NT" or DB == "NR_PROT_DB":

                    match4 = re.match(r'<a href=(.*)>(.*)\|([A-Za-z0-9._]+)\|(.*)?(\|)?(.*)?', match3.group(2))
                    if match4:

                        SeqID = match4.group(3)
                        if match4.group(2) == "pdb":

                            SeqID = match4.group(3) + match4.group(4)

                if DB == "UniProt" or DB == "SWISS-PROT":

                    match4 = re.match(r'([A-Za-z0-9]+)\|([A-Za-z0-9]+_[A-Za-z0-9]+)', match3.group(2))
                    if match4:

                        SeqID = match4.group(2)

                elif DB == "CLEAN_UNIPROT":

                    match4 = re.match(r'(\S+_\S+)\|(\S+)', match3.group(2))
                    if match4:

                        SeqID = match4.group(1)

                elif DB == "UNIREF90":

                    match4 = re.match(r'(\S+_\S+)', match3.group(2))
                    if match4:

                        SeqID = match4.group(1)

                elif DB == "CULLED_PDB" or DB == "PDB":

                    match4 = re.match(r'^([A-Z0-9]{5}) [0-9]+', match3.group(2))
                    if match4:

                        SeqID = match4.group(1)

                if Multi_Select == "no":

                    BLAST_FORM.write("%d<input type=\"radio\" value=\"%s\" name=\"MarkedSeq\" id=\"%s.align\"\">%s" %(Counter_Align_Section, SeqID, SeqID, line))

                else:

                    BLAST_FORM.write("%d<input type=\"checkbox\" id=\"%s.align\" onchange=\"document.getElementById('%s.evalue').checked=document.getElementById('%s.align').checked\">%s" %(Counter_Align_Section, SeqID, SeqID, SeqID, line))

                Counter_Align_Section += 1
                printed = 1

            if printed == 0:

                BLAST_FORM.write(line)

            line =  BLAST_HTML_RESULTS.readline()

    BLAST_FORM.write("</FONT>\n")

    if Multi_Select == "yes":

        BLAST_FORM.write("""<a href="#Select_Blast_Hit">Go to top</a>
        <input type="submit" value="submit" id="submit" name="submit" onclick="javascript:return VarifySelect_Box();">
        </form>
        </div>
        <!--Select MSA Seq ends here-->
        <SCRIPT LANGUAGE="JavaScript">
        <!-- Begin


        function getRadioVal(RadioName)
        {
            var radios = document.getElementsByName(RadioName);
            var RadVal;
            for (var i = 0, length = radios.length; i < length; i++) {
                if (radios[i].checked) {
                    RadVal=radios[i].value;
                    }
        }
        return RadVal;
        }
        function checkAll() {
            SelectFirstSeqWithJump (500,'1');
        //     var form=document.Select_Blast_Hit;
        //     var boxes = form.getElementsByTagName('input');

        //     for (var i = 0; ((i < boxes.length) && (i < 508)); i++) {
        //          if (boxes[i].type == 'checkbox'){
        //                boxes[i].checked = true;
        //          }
        //     }
        }

        function clearAll() {
            if(confirm("Uncheck selected sequences?"))
            {
                var form=document.Select_Blast_Hit;
                var boxes = form.getElementsByTagName('input');

                for (var i = 0; i < boxes.length; i++) {
                    if (boxes[i].type == 'checkbox'){
                        boxes[i].checked = false;
                    }
                }
            }
        }

        function VarifySelect_Box() {
            var form=document.Select_Blast_Hit;
            var boxes = form.getElementsByTagName('input');
            var NoneChecked=1;
            var Num_of_Seq=0;

            for (var i=0; i<boxes.length;i++){
                if (boxes[i].type == 'checkbox'){
                    if (boxes[i].checked == true){
                        NoneCheacked=0;
                        Num_of_Seq++;
                    }
                else{
                    NoneCheacked=1;
                    }
                }
            }
            Num_of_Seq=Num_of_Seq/2;
            if(Num_of_Seq < 1){
                alert("No sequences were selected, please choose sequences and continue;");
                return false;
            }
            if (Num_of_Seq<5)
            {
                alert("Only "+Num_of_Seq+" sequences were selected, please choose at least 5 sequences and continue;");
                return false;
            }
            if (Num_of_Seq<10)
            {
                confirm("Only "+Num_of_Seq+" sequences were selected. It is recomended to select at least 10 sequences for ConSurf analysis. Do you wnat to continue?");
            }
        }
        //function VarifySelect_Box() {
        //	var form=document.Select_Blast_Hit;
        //	var boxes = form.getElementsByTagName('input');
        //	var NoneChecked=1;
        //
        //	for (var i=0; i<boxes.length;i++){
        //	   if (boxes[i].type == 'checkbox'){
        //	   	if (boxes[i].checked == true){
        //		    NoneCheacked=0;
        //		    break;
        //		    }
        //		else{
        //		    NoneCheacked=1;
        //		    }
        //	   }
        //	}
        //	if(NoneCheacked == 1){
        //	alert("No sequences were selected, please choose sequences and continue;");
        //	return false;
        //	}
        //}


        function VarifySelect_Radio() {
            var boxes = document.getElementById('Select_Blast_Hit').getElementsByTagName('input');
            var NoneChecked=1;

            for (var i=0; i<boxes.length;i++){
                if (boxes[i].type == 'radio'){
                    if (boxes[i].checked == true){
                        NoneCheacked=0;
                        break;
                    }
                    else{
                        NoneCheacked=1;
                    }
                }
            }
            if(NoneCheacked == 1){
                alert("No sequences were selected, please choose sequences and continue;");
                return false;
            }
        }
        function SelectFirstSeqWithJump (First_Homologs,Jump) {
            var First_Homologs_Int=parseInt(First_Homologs);
            var JumpInt=parseInt(Jump);
            if (First_Homologs_Int > 500 )
            {
                alert ("The maximum value allowed for homologues is 500");
            }
            else
            {

                if (Jump===JumpInt)
                {
                    alert ("Conversion failed...");
                }

                var form=document.getElementById('Select_Blast_Hit');
                var boxes = form.getElementsByTagName('input');
                var TotalNum_of_Seq=(boxes.length-9)/2;
                for (var i = 8; i < ((First_Homologs_Int*JumpInt)+8) && (i < TotalNum_of_Seq+8); i=i+JumpInt) {
                    if (boxes[i].type == 'checkbox'){ /select the e-value section checkbox/
                        boxes[i].checked = true;
                    }
                    if (boxes[(i+TotalNum_of_Seq)].type == 'checkbox'){ /select the alignment section checkbox/
                        boxes[(i+TotalNum_of_Seq)].checked = true;
                    }
                }
            }
        }
        function SelectFirstSeq (First_Homologs) {
            var form=document.Select_Blast_Hit;
            var boxes = form.getElementsByTagName('input');
            var TotalNum_of_Seq=(boxes.length-5)/2;
            for (var i = 0; i < First_Homologs; i++) {
                if (boxes[i].type == 'checkbox'){ /select the e-value section checkbox/
                    boxes[i].checked = true;
                }
                if (boxes[(i+TotalNum_of_Seq)].type == 'checkbox'){ /select the alignment section checkbox/
                    boxes[(i+TotalNum_of_Seq)].checked = true;
                }
            }
        }

        //  End -->
        </script>

        <\?
        }
        else
        {
            \$fh = fopen("%s", 'w') or die("Can't open file");
            fwrite(\$fh, "The Selected Seq Are:\n");
            foreach (\$selected_seq as \$f) {
                fwrite (\$fh,"\$f\n");
            }
            fclose(\$fh);
            header('Location: %s');
        }
        \?>""" %(Selected_Seq_File, Redirect_Output))

    else: # Only One seq is selceted

        BLAST_FORM.write("""<input type=\"submit\" value=\"submit\" name=\"submit\">
        <\?
        }
           else
        {
            \$fh = fopen("%s", 'w') or die("Can't open file");
            fwrite (\$fh,"\$selected_seq");
            fclose(\$fh);
        }
        \?>""" %Selected_Seq_File)

    BLAST_FORM.close()
    BLAST_HTML_RESULTS.close()
    return("OK")








