import os 
from argparse import ArgumentParser


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

#==============================================
#FINAL SIZE FILTERING FUNCTION
def size_filter(infile, minlen, maxlen):
    print("Size selection process of consensus sequences with lenght over " + str(minlen*0.9) + "bp and under " + str(maxlen*1.1) + " has been started...")
    with open(infile, "r") as ifp:
        ls = ifp.readlines()
    ifp.close()
    ext = os.path.splitext(infile)[1]
    if ext==".fasta":
        outfile = os.path.join(get_base_dir(infile)[0], get_base_dir(infile)[1].split(".")[0]+"_size_selected.fasta")
        with open(outfile, "w") as ofp:
            print("Writing size selected consensus sequences from %s to %s" %(infile, outfile))
            total=int(len(ls))/2
            ass = 0
            for el in range(len(ls)):
                if ls[el].startswith(">")==True and maxlen*1.1 > len(ls[el+1])-1 > minlen*0.9:
                    ofp.write(ls[el])
                if ls[el].startswith(">")==False and maxlen*1.1 > len(ls[el])-1 > minlen*0.9:
                    ofp.write(ls[el])
                    ass+=1
                else:
                    continue
            print("Finished writing")
        ofp.close()
        print("Size selection process ended: %d sequences were discared out of %d total sequences (%g percent)" %(total-ass, total, 100-(round(ass/total, 4))*100))
    if ext==".fastq":
        outfile = os.path.join(get_base_dir(infile)[0], get_base_dir(infile)[1].split(".")[0]+"_size_selected.fastq")
        total=int(len(ls)/4)
        ass=0
        with open(outfile, "w") as ofp:
            print("Writing size selected consensus sequences from %s to %s" %(infile, outfile))
            for el in range(len(ls)):
                if ls[el].startswith("@")==True and ls[el-1].startswith("+")==False and maxlen*1.1 > len(ls[el+1])-1 > minlen*0.9:
                    ofp.write(ls[el])
                    ass+=1
                if ls[el].startswith("@")==False and ls[el-1].startswith("@")==True and maxlen*1.1 > len(ls[el])-1 > minlen*0.9:
                    ofp.write(ls[el])
                else:
                    continue
            print("Finished writing")
        ofp.close()
        print("Size selection process ended: %d sequences were discarded out of %d total sequences (%g percent)" %(total-ass, total, 100-(round(ass/total, 4))*100))
#==============================================

if __name__ == "__main__":
    size_filter(cons, minimum, maximum)
