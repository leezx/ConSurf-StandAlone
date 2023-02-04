

import re

def createWasabiXml(fname_msa, fname_tree, fname_out):

    try:

        OUT = open(fname_out, 'w')

    except:

        return("Could not open file " + fname_out + " for writing.")

    OUT.write("<ms_alignment>\n")

    printNewickTree(fname_tree, OUT)
    printNodes(fname_msa, OUT)

    OUT.write("</ms_alignment>\n")
    OUT.close()

    return("OK")

def printNewickTree(fname, OUT):

    OUT.write("\t<newick>")

    try:

        TREE = open(fname, 'r')

    except:

        return("Could not open " + fname + " for reading.")

    line = TREE.readline()
    while line != "":

        line = line.strip()
        line = line[:-1] # remove the character ;
        OUT.write(line)
        line = TREE.readline()

    OUT.write("</newick>\n")
    TREE.close()

    return("OK")

def printNodes(fname, OUT):

    OUT.write("<nodes>\n")

    try:

        NODES = open(fname, 'r')

    except:

        return("Could not open file " + fname + " for reading.")

    Seq = ""

    line = NODES.readline()
    while line != "":

        line = line.strip()
        if re.match(r'^>.*', line):

            if Seq != "":

                OUT.write("  <sequence>\n")
                OUT.write("    %s\n" %Seq)
                OUT.write("  </sequence>\n")
                OUT.write("</leaf>\n")

                Seq = ""

            SeqName = line[1:]
            OUT.write("<leaf id=\"%s\" name=\"%s\">\n" %(SeqName, SeqName))

        else:

            Seq += line

        line = NODES.readline()

    OUT.write("  <sequence>\n")
    OUT.write("    %s\n" %Seq)
    OUT.write("  </sequence>\n")
    OUT.write("</leaf>\n")
    OUT.write("</nodes>\n")
    NODES.close()

    return("OK")


















