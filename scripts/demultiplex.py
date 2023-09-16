import os
import subprocess as sp
import pandas as pd
from argparse import ArgumentParser
from datetime import datetime
import sys
import edlib

argparse = ArgumentParser()
argparse.add_argument("-i", "--infile", help="Path to the file with raw reads", required=True)
argparse.add_argument("-p", "--primers", help="Path to the csv file containing primers", required=True)

args = argparse.parse_args()

inf = args.infile
prim = args.primers


############################################################################
                #THE FOLLOWING BLOCKS, WHOSE TITLE LINES ARE
                #STARTING WITH #@, ENCLOSE GENERAL FUNCTIONS
                #THAT SHOULD BE IGNORED IN THE MAINTENANCE OF 
                #THE CODE

#@==========================================================================
#@GENERAL OPERATIONS ON PATH AND FOLDERS FUNCTIONS

def reverse_str(string):
    a = []
    b = ""
    for i in string:
        a.append(i)
    a.reverse()
    return b.join(a)

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
#@GENERAL OPERATION ON LISTS FUNCTIONS:

def whole_length(ls_of_ls):
    a = 0
    for ls in ls_of_ls:
        a+=len(ls)
    return a

def merge_lists_wisely(l1, l2, target):
    merged_list = target
    for i in l1:
        if i not in merged_list:
            merged_list.append(i)
    for j in l2:
        if j not in merged_list:
            merged_list.append(j)
    return merged_list

def double_ls_of_ls(lis):
    overall = []
    doubles = []
    for i in lis:
        if type(i)==list:
            for el in i:
                overall.append(el)
        else:
            overall.append(i)
    for j in overall:
        if overall.count(j)>1 and j not in doubles:
            doubles.append(j)
    return len(doubles), doubles


def find_n_remove(query,target_r):
    a = 0
    for i in target_r:
        for el in range(len(i)):
            if el<len(i):
                if i[el]==query:
                    i.remove(i[el])
                    a+=1
                else:
                    continue
            else:
                pass
    return a
      
def listcompare(l1, l2):
    if len(l1)!=len(l2):
        return False
    else:
        c=0
        for i in range(len(l1)):
            if l1[i]==l2[i]:
                c+=1
        if c==len(l1):
            return True
        return False

##############################################################################
def reverse_complement(seq):
    rev = {'\n': '\n', 'A': 'T','T': 'A','C': 'G','G': 'C','U': 'A','R': 'Y','Y': 'R','M': 'K','K': 'M','S': 'S','W': 'W','H': 'D','B': 'V','V': 'B','D': 'H','N': 'N'}
    rev_comp=[]
    for base in seq:
        rev_comp.append(base)
    rev_comp.reverse()
    for ind in range(len(rev_comp)):
        rev_comp[ind] = rev[rev_comp[ind]]
    sequence = ""
    rev_comp_sequence = sequence.join(rev_comp)
    return rev_comp_sequence

def get_aln_pos(aln, query, seq):
    c = edlib.getNiceAlignment(aln, seq, query)
    d=[]
    for i in query:
        d.append(i)
    index = []
    for i in range(len(c['target_aligned'])):
        if c['target_aligned'][i] in d:
            index.append(i)
        else:
            pass
    return index[len(index)-1]


def unzip_if_gzipped(gzfile):
    if os.path.splitext(gzfile)[1]==".gz":
        sp.run("pigz -d " + gzfile, shell=True, check=True)
        return gzfile.split(".")[0]+"."+gzfile.split(".")[1]
    return gzfile

def read_raw_fastaq(infile):
    fastaq = unzip_if_gzipped(infile)
    if os.path.splitext(fastaq)[1]==".fasta":
        lines_dict = {}
        with open(fastaq, "r") as fq:
            lines = fq.readlines()
        fq.close()
        for i in range(len(lines)):
            if lines[i].startswith(">"):
                lines_dict.update({lines[i]:lines[i+1]})
            else:
                continue
        return lines_dict
    if os.path.splitext(fastaq)[1]==".fastq":
        lines_dict = {}
        with open(fastaq, "r") as fa:
            lines = fa.readlines()
        fa.close()
        for i in range(len(lines)):
            if i==0 or (i>0 and lines[i].startswith("@") and lines[i-1].startswith("+")==False):
                lines_dict.update({lines[i]:lines[i+1]})
            else:
                continue
        return lines_dict
    else:
        print("Invalid file format, file should be provided as fasta or fastq")
        raise ValueError

def read_primer_table(csv):
    primers = pd.read_csv(csv)
    fwd_prim = primers["FWD"]
    rev_prim = primers["REV"]
    ids_prim = primers["ID"]
    fp = fwd_prim
    rp = rev_prim
    ids = ids_prim
    fwd_dict={}
    rev_dict= {}
    for i in range(len(ids)):
        fwd_dict.update({ids[i]:fp[i]})
        rev_dict.update({ids[i]:rp[i]})
    return fwd_dict, rev_dict

def demultiplex(fastaq, csv):
    print("Demultiplexing process started")
    correspondences = [("R", "A"),("R", "G"),("Y", "C"),("Y", "T"),("M", "A"),("M", "C"),("K", "G"),("K", "T"),("S", "G"),("S", "C"),("W", "A"),("W", "T"),("B", "C"),("B", "G"),("B", "T"),("D", "A"),("D", "G"),("D", "T"),("H", "A"),("H", "C"),("H", "T"),("V", "A"),("V", "C"),("V", "G"),("N", "A"),("N", "C"),("N", "G"),("N", "T")]
    read_dict = read_raw_fastaq(fastaq)
    fwd_dict, rev_dict = read_primer_table(csv)
    nums = []
    reads = list(read_dict.values())
    headers = list(read_dict.keys())
    path = makedir_orchange(os.path.join(get_base_dir(fastaq)[0], get_base_dir(fastaq)[1]+os.path.splitext(fastaq)[1][1]+os.path.splitext(fastaq)[1][len(os.path.splitext(fastaq)[1])-1]+"_demultiplexed"))
    demultiplexed = []
    for fwd in list(fwd_dict.keys()):
        ass = 0
        with open(os.path.join(path, "id_"+fwd+".fasta"), "w") as id_fa:
            for i in range(len(reads)):
                if i<len(reads):
                    alignment1 = edlib.align(reads[i], fwd_dict[fwd],mode="SHW", task="path", additionalEqualities=correspondences)
                    if get_aln_pos(alignment1,fwd_dict[fwd], reads[i]) < 100:
                        alignment2=edlib.align(reverse_complement(reads[i]), rev_dict[fwd], mode="SHW", task="path", additionalEqualities=correspondences)
                        if get_aln_pos(alignment2,rev_dict[fwd], reads[i])<100:
                            ass += 1
                            id_fa.write(">"+headers[i][1:])
                            demultiplexed.append(">"+headers[i][1:])
                            id_fa.write(reads[i])
                    else:
                        alignment1 = edlib.align(reads[i],rev_dict[fwd],mode="SHW", task="path", additionalEqualities=correspondences)
                        if get_aln_pos(alignment1,rev_dict[fwd], reads[i]) < 100:
                            alignment2=edlib.align(reverse_complement(reads[i]), fwd_dict[fwd], mode="SHW", task="path", additionalEqualities=correspondences)
                            if get_aln_pos(alignment2,fwd_dict[fwd], reads[i])< 100:
                                ass += 1
                                id_fa.write(">"+headers[i][1:])
                                id_fa.write(reads[i])
                        else:
                            continue
                else:
                    pass  
        nums.append(ass)         
        print(str(ass) + " sequences were assigned to id " + fwd)
        id_fa.close()
    print("Demultiplexing finished, %d reads were assigned" %(sum(nums)))

if __name__=="__main__":
    start = datetime.now()
    demultiplex(inf, prim)
    end = datetime.now()
    print('Duration: {}'.format(end - start), file=sys.stderr)
                        