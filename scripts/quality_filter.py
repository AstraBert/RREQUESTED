from argparse import ArgumentParser
from datetime import datetime
import sys
from Bio import SeqIO
import gzip

argparse = ArgumentParser()
argparse.add_argument("-i", "--infile", help="Path to the file with raw reads", required=True)
argparse.add_argument("-q", "--quality", help="Minimum average quality score for filtering", required=True, type=float)

args = argparse.parse_args()

inf = args.infile
qlt = args.quality

#FILTERING FUNCTIONS

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
    else:
        raise ValueError("File is the wrong format")
    print("Done", file=sys.stderr)

def ascii_conv_and_mean(line):
    phred_quality_dict = {
    '!' : 0, '"' : 1, '#' : 2, '$' : 3, '%' : 4, '&' : 5, "'" : 6, '(' : 7, ')' : 8, '*' : 9,
    '+' : 10, ',' : 11, '-' : 12, '.' : 13, '/' : 14, '0' : 15, '1' : 16, '2' : 17, '3' : 18, '4' : 19,
    '5' : 20, '6' : 21, '7' : 22, '8' : 23, '9' : 24, ':' : 25, ';' : 26, '<' : 27, '=' : 28, '>' : 29,
    '?' : 30, '@' : 31, 'A' : 32, 'B' : 33, 'C' : 34, 'D' : 35, 'E' : 36, 'F' : 37, 'G' : 38, 'H' : 39,
    'I' : 40, 'J' : 41, 'K' : 42, 'L' : 43, 'M' : 44, 'N' : 45, 'O' : 46, 'P' : 47, 'Q' : 48, 'R' : 49,
    'S' : 50, 'T' : 51, 'U' : 52, 'V' : 53, 'W' : 54, 'X' : 55, 'Y' : 56, 'Z' : 57, '[' : 58, '\\' : 59,
    ']' : 60, '^' : 61, '_' : 62, '`' : 63, 'a' : 64, 'b' : 65, 'c' : 66, 'd' : 67, 'e' : 68, 'f' : 69,
    'g' : 70, 'h' : 71, 'i' : 72, 'j' : 73, 'k' : 74, 'l' : 75, 'm' : 76, 'n' : 77, 'o' : 78, 'p' : 79,
    'q' : 80, 'r' : 81, 's' : 82, 't' : 83, 'u' : 84, 'v' : 85, 'w' : 86, 'x' : 87, 'y' : 88, 'z' : 89,
    '{' : 90, '|' : 91, '}' : 92, '~' : 93
    }
    mean_list = []
    for i in line:
        if i != '\n':
            mean_list.append(phred_quality_dict[i])
    return round(sum(mean_list)/len(mean_list), 3)

def filter(fastq, q):
    print("Filtering of reads over quality threshold " + str(q) + " has been started")
    readdict=load_data(fastq)
    with open(fastq, "r+") as f:
        f.truncate()
    f.close()
    disc = 0
    total = int(len(list(readdict.keys())))
    with open(fastq, "w") as fp:
        for i in list(readdict.keys()):
            if ascii_conv_and_mean(readdict[i][1])<q:
                pass
            else:
                fp.write("@"+i+"\n")
                fp.write(readdict[i][0]+"\n")
                fp.write("+"+i+"\n")
                fp.write(readdict[i][1]+"\n")
    fp.close()
    print("Filtering finished: %d reads were discarded out of %d (%g percent), now re-compiling the file with filtered reads..." %(disc, total, 100*(round(disc/total, 4))))


if __name__=="__main__":
    filter(inf, qlt)
