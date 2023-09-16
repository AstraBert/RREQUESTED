from argparse import ArgumentParser
from datetime import datetime
import sys

argparse = ArgumentParser()
argparse.add_argument("-i", "--infile", help="Path to the file with raw reads", required=True)
argparse.add_argument("-q", "--quality", help="Minimum average quality score for filtering", required=True, type=float)

args = argparse.parse_args()

inf = args.infile
qlt = args.quality

#FILTERING FUNCTIONS
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
    with open(fastq, "r+") as f:
        a = f.readlines()
        f.seek(0)
        f.truncate()
    f.close()
    disc = 0
    total = int(len(a)/4)
    for i in range(len(a)):
        if i > 2 and str(a[i-1]).startswith("+")==True and str(a[i-3]).startswith("@")==True:
            mean_q = ascii_conv_and_mean(a[i]) 
            if mean_q < q:
                a[i-3] = 0
                a[i-2] = 0
                a[i-1] = 0
                a[i] = 0
                disc+=1
    print("Filtering finished: %d reads were discarded out of %d (%g percent), now re-compiling the file with filtered reads..." %(disc, total, 100*(round(disc/total, 4))))
    with open(fastq, "w") as fp:
        for el in a:
            if el != 0:
                fp.write(el)
    fp.close()

if __name__=="__main__":
    filter(inf, qlt)