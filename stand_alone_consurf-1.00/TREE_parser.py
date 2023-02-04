import re
import os

def check_validity_tree_file(treeFile, ref_to_err):

    # check the validity of the newick format of the uploaded tree
    # input: path to a tree file, reference to error hash
    #
    # if find error, mark it through the error hash:
    # the marked items are: 'bootstrap' - if found bootstrap value,
    #   'internal_nodes' - if found internal nodes
    #   'left_right' - if the brackets are not equal in number
    #   'noRegularFormatChar' - a list of all the non standard characters

    lineCounter = 0
    rightBrackets = 0
    leftBrackets = 0
    tree = ""
    noRegularFormatChar  = ""
    tree = ""
    errorBool = 0
    internal_nodes = 0
    bootstrap = 0
    missingSemicolon="no"
    read_right_bracket = "no"

    # in case the uploaded tree file contains more than one line -
    # read the tree and rewrite it
    try:

        TREEFILE = open(treeFile, 'r')

    except:

        ref_to_err['error'] = "could not read " + treeFile

    line = TREEFILE.readline()
    while line != "":

        line.rstrip()
        tree += line
        lineCounter += 1
        line = TREEFILE.readline()

    TREEFILE.close()

    # add a semi-colon if missing
    if tree[-1] != ';':

        tree += ";"
        missingSemicolon = "yes"

    if lineCounter > 1 or missingSemicolon == "yes":

        try:

            TREEFILE = open(treeFile, 'w')

        except:

            ref_to_err['error'] = "could not write to file " + treeFile

        TREEFILE.write(tree)
        TREEFILE.close()

    # legal tree: same number of left and right brackets, no irregular chars
    lineArr = tree.split()
    for chain in lineArr:

        if chain == '(':

            leftBrackets += 1
            read_right_bracket = "no"

        elif chain == ')':

            rightBrackets += 1
            read_right_bracket = "yes"

        else:

            match1 = re.match(r'([\!\@\#\$\^\&\*\~\`\{\}\'\?\\\/\<\>])', chain)
            if match1:

                char = match1.group(1)
                if not char in noRegularFormatChar:

                    noRegularFormatChar += " '" + char + "', "
                    read_right_bracket = "no"

            else:

                # if right after a right Bracket we read a character which is not legal (ie: , : ;)
                # we output a message to the user, since we don't handle bootstrap values or internal node names
                if read_right_bracket == "yes":

                    match2 = re.match(r'\d', chain)
                    if match2:

                        ref_to_err['bootstrap'] = 1

                    else:

                        match3 = re.match(r'[,|:|;]', chain)
                        if not match3:

                            ref_to_err['internal_nodes'] = 1

                read_right_bracket = "no"

    if leftBrackets != rightBrackets:

        ref_to_err['left_right'] = 1

    if noRegularFormatChar != "":

        noRegularFormatChar = re.sub(r'\,\s$', "", noRegularFormatChar)
        ref_to_err['noRegularFormatChar'] = noRegularFormatChar

def extract_nodes_from_tree(tree, ref_tree_nodes):

    tree_arr = tree.split('(')
    sub_tree = {}
    temp_arr = []
    sub_counter = 0

    # building the array sub_tree, so that each cell will hold maximum one sequence name
    i = 0
    while i < len(tree_arr):

        if tree_arr[i] != "":

            tree_arr[i] = "(" + tree_arr[i]

        match1 = re.match(r'.*,.+', tree_arr[i])
        if match1:

            temp_arr = (tree_arr[i]).split(',')
            for a in temp_arr:

                if a != "":

                    sub_tree[sub_counter] = a + ","
                    sub_counter += 1

        else:

            sub_tree[sub_counter] = tree_arr[i]
            sub_counter += 1

        i += 1

    # extract the nodes
    exp = ""
    new_rest_exp = ""
    seq_found = "no"
    k = 1
    while k < len(sub_tree):

        #in this part we wish to split the expression to 2 parts; left part : (?seq_name ; right part: all the rest
        if sub_tree[k] != "":

            match2 = re.match(r'(.+)(:.+)', sub_tree[k])
            if match2:

                while match2:

                    exp = match2.group(1)
                    match2 = re.match(r'(.+)(:.+)', exp)

            else:

                # in case the expression is of format:  seq_name:distance,
                match3 = re.match(r'(.+)(\);.+)', sub_tree[k])
                if match3:

                    exp = match3.group(1)
                    match4 = re.match(r'(.+)(\))', exp)
                    while match4:

                        exp = match4.group(1)
                        match4 = re.match(r'(.+)(\))', exp)

                else:

                    #  in case the expression is of format:  seq_name)*,
                    match5 = re.match(r'(.+)(\)?.+)', sub_tree[k])
                    if match5:

                        exp = match5.group(1)
                        match6 = re.match(r'(.+)(\))', exp)
                        while match6:

                            exp = match6.group(1)
                            match6 = re.match(r'(.+)(\))', exp)

            match5 = re.match(r'(\(?)(.+)', exp)
            if match5:

                if match5.group(2) in ref_tree_nodes:

                    return('user_error', "duplicity: " + match5.group(2))

                else:

                    ref_tree_nodes.append(match5.group(2))

        k += 1

    return ["OK"]

def removeBPvalues(IN_treeFile, OLD_treeFile):

    TREEFILE = open(IN_treeFile, 'r')
    treeFileOneLine = re.sub(r'\n', "", TREEFILE.read())
    TREEFILE.close()

    # replace bootstrap values  which look like this: ((A:0.02,B:0.03)40:0.3);
    treeFileOneLine = re.sub(r'\)\d*\.?\d+\:', "\)\:", treeFileOneLine)

    # replace bootstrap values  which look like this:(A:0.4,(B:0.1,C:0.1):0.3[40]);
    treeFileOneLine = re.sub(r'(\d*\.?\d+)\[\d*\.?\d+\]', r'\1', treeFileOneLine)

    os.rename(IN_treeFile, OLD_treeFile)

    TREEFILE = open(IN_treeFile, 'w')
    TREEFILE.write(treeFileOneLine + "\n")
    TREEFILE.close()













