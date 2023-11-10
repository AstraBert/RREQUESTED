import os
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from datetime import datetime
import sys

argparse = ArgumentParser()
argparse.add_argument("-d", "--directory", help="Path to the input files directory", required=True)

args = argparse.parse_args()

indr = args.directory


def get_base_dir(path):
    base, ext= os.path.splitext(path)
    asedir = base.split("/")
    b = "/"
    basedir = b.join(asedir[:len(asedir)-1])
    basename = asedir[len(asedir)-1]
    return basedir, basename

def check_clusters(files):
    dict_of_dicts={}
    qcovs=[]
    pidents=[]
    qlens=[]
    for fp in files:
        print("Reading BLAST results file to get out annotation stats...", file=sys.stderr)
        with open(fp, "r") as f:
            lines = f.readlines()
        f.close()
        diction = {}
        qcov=0
        piden=0
        qlen=0
        for i in range(len(lines)):
            if lines[i].split("\t")[1] not in diction.keys():
                a = lines[i].split("\t")[1]
                diction.update({a: 1})
                qcov+=float(lines[i].split("\t")[5])
                piden+=float(lines[i].split("\t")[4])
                qlen+=float(lines[i].split("\t")[3])
            else:
                diction[lines[i].split("\t")[1]]+=1
                qcov+=float(lines[i].split("\t")[5])
                piden+=float(lines[i].split("\t")[4])
                qlen+=float(lines[i].split("\t")[3])
        if len(lines)>0:
            dict_of_dicts.update({get_base_dir(fp)[1]: diction})
            qcovs.append(round((qcov/len(lines)),3))
            pidents.append(round((piden/len(lines)),3))
            qlens.append(round((qlen/len(lines)),3))
        else:
            continue
    keys=list(dict_of_dicts.keys())
    rate_of_contamination=[]
    qkeys=[]
    for j in range(len(keys)):
        qkeys.append(keys[j])
        print(keys[j]+":")
        print("Average query coverage = %g\nAverage percentage identity = %g\nAverage query length = %g" %(qcovs[j], pidents[j], qlens[j]))
        for item in list(dict_of_dicts[keys[j]].items()):
            print(str(item[0])+"\t"+str(item[1]))
        rate_of_contamination.append(100-round((max(list(dict_of_dicts[keys[j]].values()))/sum(list(dict_of_dicts[keys[j]].values())))*100, 3))
    plt.style.use("ggplot")
    plt.bar(range(len(rate_of_contamination)), rate_of_contamination, color="#b134eb")
    plt.xticks(range(len(rate_of_contamination)), qkeys, rotation=60)
    plt.title("Contamination rate in demultiplexing")
    plt.xlabel("Demultiplexed file")
    plt.ylabel("Contamination(%)")
    plt.axhline(100, color="white")
    plt.savefig(os.path.join(get_base_dir(files[0])[0],"demultiplexing_contamination.png"))
    plt.show()

if __name__=="__main__":
    files=[os.path.join(indr, f) for f in os.listdir(indr) if os.path.isfile(os.path.join(indr, f)) and f.endswith(".blast")]
    check_clusters(files)
