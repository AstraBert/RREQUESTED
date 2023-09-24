import os
import subprocess as sp
import edlib
from datetime import datetime
from argparse import ArgumentParser
import sys

argparse = ArgumentParser()
argparse.add_argument("-i", "--infile", help="Path to the input file", required=True)
# argparse.add_argument("-dn", "--demultiplex_number", help="Number of demulexing units in raw reads file", required=True, type=int)

args = argparse.parse_args()

inf = args.infile
# dem_num = args.demultiplex_number
############################################################################
                #THE FOLLOWING BLOCKS, WHOSE TITLE LINES ARE
                #STARTING WITH #@, ENCLOSE GENERAL FUNCTIONS
                #THAT SHOULD BE IGNORED IN THE MAINTENANCE OF 
                #THE CODE

#@==========================================================================
#@GENERAL OPERATIONS ON PATH AND FOLDERS FUNCTIONS

def makedir_orchange(path):
    try:
        os.mkdir(path)
        return path
    except FileExistsError:
        return path

def get_base_dir(path):
    base, ext= os.path.splitext(path)
    asedir = base.split("/")
    b = "/"
    basedir = b.join(asedir[:len(asedir)-1])
    base_basedir = b.join(asedir[:len(asedir)-2])
    basename = asedir[len(asedir)-1]
    return basedir, basename, base_basedir

#@==========================================================================
#OPERATIONS ON LISTS
def find_n_remove(l1, l2):
    for j in range(len(l1)):
        for k in range(len(l2)):
            if k<len(l2) and l1[j]==l2[k]:
                l2.remove(l2[k])
            else:
                continue
    return l2

def reverse_complement(seq):
    rev = {'\n':'\n','A': 'T','T': 'A','C': 'G','G': 'C','U': 'A','R': 'Y','Y': 'R','M': 'K','K': 'M','S': 'S','W': 'W','H': 'D','B': 'V','V': 'B','D': 'H','N': 'N'}
    rev_comp=[]
    for base in seq:
        rev_comp.append(base)
    rev_comp.reverse()
    for ind in range(len(rev_comp)):
        rev_comp[ind] = rev[rev_comp[ind]]
    sequence = ""
    rev_comp_sequence = sequence.join(rev_comp)
    return rev_comp_sequence
  
#@==========================================================================

###############################################################################

                    #THE FOLLOWING FUNCTIONS, WHOSE 
                    #TITLE LINES START WHIT #>
                    #SHOULD BE INSTEAD CONSIDERED,
                    #AS THEY DO MOST OF THE WORK :)
###############################################################################
                    
#>==========================================================================
#FIND REFERENCE SEQUENCES

def find_the_num(infile):
    print("Started the search for higly diverging sequences (>50 perc apart from any other)", file=sys.stderr)
    the_num=[]
    try:
        with open(infile, "r") as f:
            lines = f.readlines()
        f.close()
        if os.path.splitext(infile)[1]==".fastq":
            seqs=0
            for i in range(len(lines)):
                if i==0 or (i>0 and lines[i].startswith("@") and lines[i-1].startswith("+")==False):
                    seqs+=1
                    if len(the_num)==0:
                        the_num.append(lines[i+1])
                    else:
                        c=0
                        for j in the_num:
                            alignment1=edlib.align(lines[i+1], j, k=len(j)*0.5, task="distance", mode="NW")
                            if alignment1["editDistance"]==-1:
                                alignment2=edlib.align(lines[i+1], reverse_complement(j), k=len(j)*0.5, task="distance", mode="NW")
                                if alignment2["editDistance"]==-1:
                                    c+=1
                                else:
                                    break
                            else:
                                break
                        if c==len(the_num):
                            the_num.append(lines[i+1])
                        else:
                            continue
        if os.path.splitext(infile)[1]==".fasta":
            seqs=0
            for i in range(len(lines)):
                if i==0 or (i>0 and lines[i].startswith(">")):
                    seqs+=1
                    if len(the_num)==0:
                        the_num.append(lines[i+1])
                    else:
                        c=0
                        for j in the_num:
                            alignment1=edlib.align(lines[i+1], j, k=len(j)*0.5, task="distance", mode="NW")
                            if alignment1["editDistance"]==-1:
                                alignment2=edlib.align(lines[i+1], reverse_complement(j), k=len(j)*0.5, task="distance", mode="NW")
                                if alignment2["editDistance"]==-1:
                                    c+=1
                                else:
                                    break
                            else:
                                break
                        if c==len(the_num):
                            the_num.append(lines[i+1])
                        else:
                            continue
        print("Search finished:\nNumber of analysed sequences: %d\nNumber of higly divergent sequences: %d" %(seqs, len(the_num)), file=sys.stderr)
        return the_num
    except KeyboardInterrupt:
        sys.exit()

def find_the_num_list(lines):
    print("Started the search for higly diverging sequences (>50 perc apart from any other)", file=sys.stderr)
    the_num=[]
    try:
        seqs=0
        for i in range(len(lines)):
            seqs+=1
            if len(the_num)==0:
                the_num.append(lines[i])
            else:
                c=0
                for j in the_num:
                    alignment1=edlib.align(lines[i], j, k=len(j)*0.5, task="distance", mode="NW")
                    if alignment1["editDistance"]==-1:
                        alignment2=edlib.align(lines[i], reverse_complement(j), k=len(j)*0.5, task="distance", mode="NW")
                        if alignment2["editDistance"]==-1:
                            c+=1
                        else:
                            break
                    else:
                        break
                    if c==len(the_num):
                        the_num.append(lines[i])
                    else:
                        continue
        print("Search finished:\nNumber of analysed sequences: %d\nNumber of higly divergent sequences: %d" %(seqs, len(the_num)), file=sys.stderr)
        return the_num
    except KeyboardInterrupt:
        sys.exit()

#>==========================================================================
#DEMULTIPLEX

def demultiplex(infile):
    try:
        value=False
        refseqs=find_the_num(infile)
        groups=[]
        ind=-1
        with open(infile, "r") as f:
            lines = f.readlines()
        f.close()
        seq=[]
        for i in range(len(lines)):
            if i==0 or (i>0 and (lines[i].startswith("@") or lines[i].startswith(">")) and lines[i-1].startswith("+")==False):
                seq.append(lines[i+1])
            else:
                continue
        tot=len(seq)
        print("Pairwise alignment assigning of raw reads to the higly divergent sequences started", file=sys.stderr)
        processed=[]
        for j in refseqs:
            groups.append([])
            ind+=1
            if os.path.splitext(infile)[1]==".fastq":
                seqs=0
                for i in range(len(lines)):
                    if i==0 or (i>0 and lines[i].startswith("@") and lines[i-1].startswith("+")==False):
                        seqs+=1                        
                        alignment1=edlib.align(lines[i+1], j, k=len(j)*0.3, task="distance", mode="NW")
                        if alignment1["editDistance"]==-1:
                            alignment2=edlib.align(lines[i+1], reverse_complement(j), k=len(j)*0.3, task="distance", mode="NW")
                            if alignment2["editDistance"]==-1:
                                continue
                            else:
                                groups[ind].append(lines[i+1])
                                processed.append(lines[i+1])
                        else:
                            groups[ind].append(lines[i+1])
                            processed.append(lines[i+1])
                print("Assigned %d sequences to group %d (%g perc)" %(len(groups[ind]), ind, round(((len(groups[ind])/seqs)*100),3)), file=sys.stderr)
            if os.path.splitext(infile)[1]==".fasta":
                seqs=0
                for i in range(len(lines)):
                    if lines[i].startswith(">"):
                        seqs+=1                        
                        alignment1=edlib.align(lines[i+1], j, k=len(j)*0.3, task="distance", mode="NW")
                        if alignment1["editDistance"]==-1:
                            alignment2=edlib.align(lines[i+1], reverse_complement(j), k=len(j)*0.3, task="distance", mode="NW")
                            if alignment2["editDistance"]==-1:
                                continue
                            else:
                                groups[ind].append(lines[i+1])
                                processed.append(lines[i+1])
                        else:
                            groups[ind].append(lines[i+1])
                            processed.append(lines[i+1])
                print("Assigned %d sequences to group %d (%g perc)" %(len(groups[ind]), ind, round(((len(groups[ind])/seqs)*100),3)), file=sys.stderr)
        no_group=find_n_remove(processed, seq)
        run=0
        if len(no_group)>0:
            print("Non-grouped sequences are %d, starting reassignment..." %(len(no_group)), file=sys.stderr)
            while len(no_group)>0 and run<5:
                print("Started round %d" %(run+1), file=sys.stderr)
                refs=find_the_num_list(no_group)
                ranks=[]
                proc=[]
                indx=-1
                for j in refs:
                    ranks.append([])
                    indx+=1
                    for i in range(len(no_group)):
                        alignment1=edlib.align(no_group[i], j, k=len(j)*0.3, task="distance", mode="NW")
                        if alignment1["editDistance"]==-1:
                            alignment2=edlib.align(no_group[i], reverse_complement(j), k=len(j)*0.3, task="distance", mode="NW")
                            if alignment2["editDistance"]==-1:
                                continue
                            else:
                                ranks[indx].append(no_group[i])
                                proc.append(no_group[i])
                        else:
                            ranks[indx].append(no_group[i])
                            proc.append(no_group[i])
                    ranks[indx].append(j)
                    print("Assigned %d sequences to no_group %d (%g perc)" %(len(ranks[indx]), indx, round(((len(ranks[indx])/len(no_group))*100),3)), file=sys.stderr)
                for r in ranks:
                    groups.append(r)
                no_group=find_n_remove(proc, no_group)
                run+=1
            if run>=4 and len(no_group)>0:
                print("Non-grouped sequences are still %d (%g perc), they will be grouped together" %(len(no_group), round((len(no_group)/tot)*100, 3)), file=sys.stderr)
                groups.append(no_group)
                value=True
        else:
            print("All sequences have been demultiplexed at first attempt, proceeding with the analysis")
        path = makedir_orchange(os.path.join(get_base_dir(infile)[0], get_base_dir(infile)[1]+"-"+os.path.splitext(infile)[1][1]+os.path.splitext(infile)[1][len(os.path.splitext(infile)[1])-1]+"-demultiplexed"))
        print("Creating a new folder at " + path + " to write demultiplexed files in there", file=sys.stderr)
        print("Writing demultiplexed sequences in separated files", file=sys.stderr)
        if os.path.splitext(infile)[1]==".fastq":
            count=0
            for g in groups:
                if len(g)>=(tot)*0.01 and groups.index(g)!=len(groups)-1:
                    count+=1
                    with open(os.path.join(path, str(count)+".fasta"), "w") as fp:
                        c=1
                        for s in g:
                            fp.write(">"+str(c)+"\n")
                            fp.write(s)
                            c+=1
                if len(g)>=(tot)*0.01 and groups.index(g)==len(groups)-1 and value:
                    with open(os.path.join(path, "nogroup.fasta"), "w") as fp:
                        c=1
                        for s in g:
                            fp.write(">"+str(c)+"\n")
                            fp.write(s)
                            c+=1
                if len(g)>=(tot)*0.01 and groups.index(g)==len(groups)-1 and value==False:
                    count+=1
                    with open(os.path.join(path, str(count)+".fasta"), "w") as fp:
                        c=1
                        for s in g:
                            fp.write(">"+str(c)+"\n")
                            fp.write(s)
                            c+=1
                else:
                    continue
        if os.path.splitext(infile)[1]==".fasta":
            count=0
            for g in groups:
                if len(g)>=(tot)*0.01 and groups.index(g)!=len(groups)-1:
                    count+=1
                    with open(os.path.join(path, str(count)+".fasta"), "w") as fp:
                        c=1
                        for s in g:
                            fp.write(">"+str(c)+"\n")
                            fp.write(s)
                            c+=1
                if len(g)>=(tot)*0.01 and groups.index(g)==len(groups)-1 and value:
                    with open(os.path.join(path, "nogroup.fasta"), "w") as fp:
                        c=1
                        for s in g:
                            fp.write(">"+str(c)+"\n")
                            fp.write(s)
                            c+=1
                if len(g)>=(tot)*0.01 and groups.index(g)==len(groups)-1 and value==False:
                    count+=1
                    with open(os.path.join(path, str(count)+".fasta"), "w") as fp:
                        c=1
                        for s in g:
                            fp.write(">"+str(c)+"\n")
                            fp.write(s)
                            c+=1
                else:
                    continue
        print("The program ended its execution successfully", file=sys.stderr)            
    except KeyboardInterrupt:
        sys.exit()

#>==========================================================================

if __name__=="__main__":
    start = datetime.now()
    demultiplex(inf)
    end = datetime.now()
    print('Duration: {}'.format(end - start), file=sys.stderr)
