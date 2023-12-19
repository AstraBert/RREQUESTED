import os 
from argparse import ArgumentParser
import sys
import gzip
from Bio import SeqIO

arg_parser = ArgumentParser()
arg_parser.add_argument("-c", "--consensus", required=False, help="Path to the file with consensus sequences")
arg_parser.add_argument("-min", "--min_len", required=False, type=int, default=500, help="Minimum PCR product length in bp (default: 500bp)")
arg_parser.add_argument("-max", "--max_len", required=False, type=int, default=10000, help="Maximum PCR product length in bp (default: 10000bp)")

args = arg_parser.parse_args()
cons = args.consensus
minimum = args.min_len
maximum = args.max_len


#==============================================
#GENERAL OPERATIONS ON PATH FUNCTION

def get_base_dir(path):
    base, ext= os.path.splitext(path)
    asedir = base.split("/")
    b = "/"
    basedir = b.join(asedir[:len(asedir)-1])
    basename = asedir[len(asedir)-1]
    return basedir, basename

def load_data(infile):
    """Load data from infile if it is in fastq format (after having unzipped it, if it is zipped)"""
    print("Reading data from input file...", file=sys.stderr)
    if infile.endswith(".gz"):  # If file is gzipped, unzip it
        y = gzip.open(infile, "rt", encoding="latin-1")
        if infile.endswith(".fastq.gz"):  # Read file as fastq if it is fastq
            records = SeqIO.parse(y, "fastq")
            seq_dict = {}  # Create a dictionary to store everything from the file
            for record in records:
                # Update dictionary with header as key, sequence as [0] element of the value list, and base quality string as [1] element of the value list
                seq_dict.update(
                    {record.id: [str(record.seq), record.format("fastq").split("\n")[3]]})
            y.close()
            return seq_dict
        elif infile.endswith(".fasta.gz") or infile.endswith(".fna.gz") or infile.endswith(".fas.gz") or infile.endswith(".fa.gz"):
                records = SeqIO.parse(y, "fasta")
                sequences = {}
                for record in records:
                    sequences.update({str(record.id): str(record.seq)})
                y.close()
                return sequences
    # Read file directly as fastq if it is a not zipped fastq
    elif infile.endswith(".fastq"):
        with open(infile, "r") as y:
            records = SeqIO.parse(y, "fastq")
            seq_dict = {}  # Create a dictionary to store everything from the file
            for record in records:
                # Update dictionary with header as key, sequence as [0] element of the value list, and base quality string as [1] element of the value list
                seq_dict.update(
                    {record.id: [str(record.seq), record.format("fastq").split("\n")[3]]})
            y.close()
            return seq_dict
    elif infile.endswith(".fasta") or infile.endswith(".fna") or infile.endswith(".fas") or infile.endswith(".fa"):
        with open(infile, "r") as y:
            records = SeqIO.parse(y, "fasta")
            sequences = {}
            for record in records:
                sequences.update({str(record.id): str(record.seq)})
            y.close()
            return sequences
    else:
        raise ValueError("File is the wrong format")
    print("Done", file=sys.stderr)

#==============================================
#FINAL SIZE FILTERING FUNCTION
def size_filter(infile, minlen, maxlen):
    print("Size selection process of sequences with lenght over " + str(round(minlen*0.9,0)) + "bp and under " + str(round(maxlen*1.1,0)) + " has been started...")
    with open(infile, "r") as ifp:
        ls = ifp.readlines()
    ifp.close()
    ext = os.path.splitext(infile)[1]
    if ext==".fasta":
        outfile = os.path.join(get_base_dir(infile)[0], get_base_dir(infile)[1].split(".")[0]+"_size_selected.fasta")
        readdict = load_data(infile)
        with open(outfile, "w") as ofp:
            print("Writing size selected sequences from %s to %s" %(infile, outfile))
            total=int(len(list(readdict.keys())))
            ass = 0
            for el in list(readdict.keys()):
                if minlen*0.9<=len(readdict[el])<=maxlen*1.1:
                    ofp.write(">"+el+"\n")
                    ofp.write(readdict[el]+"\n")
                else:
                    ass+=1
                    pass    
        print("Finished writing")
        ofp.close()
        print("Size selection process ended: %d sequences were discared out of %d total sequences (%g percent)" %(ass, total, 100-(round(1-ass/total, 4))*100))
    if ext==".fastq":
        outfile = os.path.join(get_base_dir(infile)[0], get_base_dir(infile)[1].split(".")[0]+"_size_selected.fastq")
        readdict = load_data(infile)
        total=int(len(list(readdict.keys())))
        ass = 0
        with open(outfile, "w") as ofp:
            print("Writing size selected sequences from %s to %s" %(infile, outfile))
            for el in list(readdict.keys()):
                if minlen*0.9<=len(readdict[el][0])<=maxlen*1.1:
                    ofp.write("@"+el+"\n")
                    ofp.write(readdict[el][0]+"\n")
                    ofp.write("+"+el+"\n")
                    ofp.write(readdict[el][1]+"\n")
                else:
                    ass+=1
                    pass    
        print("Finished writing")
        ofp.close()
        print("Size selection process ended: %d sequences were discared out of %d total sequences (%g percent)" %(ass, total, 100-(round(1-ass/total, 4))*100))
#==============================================

if __name__ == "__main__":
    size_filter(cons, minimum, maximum)
