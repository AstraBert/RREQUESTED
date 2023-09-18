#  QUASARd MANUAL #

### -General purpose and applicability

QUASARd (QUAlity Size and All Reads demultiplexing) is a modular shellscript tool that combines three different python scripts and allows preprocessing (quality filtering, size selection and demultiplexing) of raw basecalled reads, especially produced by Third Generation Sequences technologies (it was specifically designed for ONT).

It is easy-to-use, it requires only two python packages to be installed for the analysis and it is fairly quick, given the computational amount that it has to parse. 

As said, it was specifically taylored on ONT reads and was tested on them, but nothing prevents you to use it (cautiously) also with NGS or PacBio reads.

### -Installation 
Follow instructions in [install.sh](./scripts/install.sh) to view a comprehensive guide and set up everything.

To run the program, you will need `pandas` and `edlib` installed, at least in your working environment; if you haven't got them yet, please download the required packages at the following links:
- [pandas](https://pandas.pydata.org/)
- [edlib](https://pypi.org/project/edlib/)

After having fulfilled these basic requirements, just download the QUASARd folder and place it within your working environment with:

```bash
cd /path/to/working/environment
git init
git clone https://github.com/AstraBert/QUASARd/
```

Now run:

```bash
echo -n 'alias QUASARd="bash /absolute/path/to/QUASARd.sh"' | cat >> ~/.bash_aliases
source ~/.bash_aliases
```

to create a permanent alias that will make running the tool less verbose, and test the installation with:

```bash
QUASARd -h
```

### -Command-line arguments
QUASARd takes one required argument (the directory in which you placed all your files with raw reads, in fasta/fastq format, gzipped or not) and four optional ones (minimum average read quality, minimum and maximum length, demultiplexing table with primes): 

```
Usage: QUASARd -d, --directory RAW_READS_DIRECTORY [-dmx, --demultipex DEMULTIPLEXING_CSV] [-q,--quality QUALITY] [-mi,--min MINIMUM_LENGTH] [-ma,--max MAXIMUM_LENGTH]

  REQUIRED ARGUMENTS:
  -RAW_READS_DIRECTORY: Provide the path to the directory where raw reads are stored in fasta/fasta.gz/fastq/fastq.gz format (can be also a mixture of these extensions)

  OPTIONAL ARGUMENTS:
  -DEMULTIPLEXING_CSV: Provide a csv file (only comma-separated) whose fields must be "ID,FWD,REV", standing for sample id, forward primer and reverse primer. If not given, no demultiplexing will take place.
  -QUALITY: Provide the value of minimum average read quality for filtering (default: 7)
  -MINIMUM_LENGTH: Provide the minimum length of PCR product for size selection (default is 500)
  -MAXIMUM_LENGTH: Provide the maximum length of PCR product for size selection (default is 10000)

  
  Input QUASARd -h to show this message
```

As said, the demultiplexing table should be a csv file, and has **mandatorily** to respect the following format:

```
ID,FWD,REV
id1,TGGCCTCTTCCTTGGCCGTC,AACTCGTCTGGCTTTTCGCC
id2,AACTACTGCTCRCCAGAAAARC,GGAAATGCATAGTTGTCTGCAA
id3,TACGTGCCTATGTCCAAYGC,TGCTTGTTCATGCAGATGTAGA
id4,AGTTTTCCATCGACTCSCAGTA,AGGTGGATTTTGGTGTGTCTYTT
id5,AATATCATGGACYTGGGNATGG,GGYTTCTTGTCCTTCTGTTTSAG
id6,AGCTGTAGTCAGTAYCACAARATG,GTGTAGAGCCAGTGRTGYTT
id7,CTGTTTCCCATCAACGAGGA,CCGCTACTGAGGGAATCCTT
```

Obviously, primer sequences and demultiplexing IDs can be modified, but **not** the names of the fields.

Here's an example run:

Example run: 	

```bash
QUASARd -d /path/to/reads/folder -dmx demultiplexing_table.csv -q 10 -mi 550 -ma 8000
```

### -How does it work? ###
There are four steps in QUASARd analysis:
1. QUASARd will check that the provided directory exists *and* is actually a directory, then it will search for fasta, fastq, fasta.gz and fastq.gz files: the working directory could potentially contain also a mixture of these file types, as the program is set to recognize and treat them differently based on their format.
2. After having found the files, QUASARd creates a results directory where it will copy all the raw reads files, and then the analysis can begin.
3. The quality filtering step is based on the easiest implementation one could think of: for every read, the filtering algorithms takes the mean quality and discards the reads that are under a given value (default is 7, so this step will take place nevertheless if the file is fastq/fastq.gz)
4. The size filtering step is also based on the easiest implementation one could think of: for every read, the filtering algorithms takes the length and, it this is below the minimum or above the maximum allowed, the read gets discarded. Minimum and maximum are set by default at 500 and 10000bp, so this step will take place nevertheless with all formats of files.
5. The demultiplexing step is the last portion of the program, and takes place only if a demultiplexing table is given in the form of a csv file delimited by commas with demultiplexing ID, forward and reverse primers. It is based on super-fast local alignment of the first 100bp of every read with the forward and reverse primer (for reverse primer, the sequence gets reverse-complemented): this is accomplished thanks to edlib library, in particular through the "SHW" mode of the function edlib.align(). If the forward primer gets completely locally aligned within the first 100 bp of the read and the reverse primer gets aligned within the first 100 bp of the reverse-complemented one, the read gets written in a file starting with the ID of the primers. This is done for every ID, for every individual, and the results are stored in folders named "individualf?_demultiplexed", where thew ? could be replaced by an a if the source file is a fasta or by a q if the source file is a fastq.  

### -Benchmarking 
Benchmarking was conducted on a folder with 4 subsampled files from the first two raw reads files by Genchen et al. (2022), with SRA label SRR17577143 and SRR175177144 (can be found on NCBI). The subsample were a fasta file, a fasta.gz file, a fastq file and a fastq.gz file, each with a total of 2000 reads. So, in total, we benchmarked the program on 8000 raw reads.

QUASARd, on a normal PC, took 4min12.771s of real time and 4min07.499s of CPU time, successfully demultiplexing all of the 28 primers for all of the four files. 

Nevertheless, the program will be much slower if the input number of reads is high: for example, tried with the entire SRR17577143 as fastq (99820 reads) it took almost half an hour to complete the task. 

### -Final considerations ###
Please note that QUASARd is still experimental and may contain errors, may fail/take really long while performing large analyses with limited computational power (e.g. on a normal laptop) and may output not-100%-reliable results, so always check them and pull issues whenever you feel it to be the case, we'll be on your back as soon as possible to fix/implement/enhance whatever you suggest!


### -License and rights of usage ###
The code is protected by the GNU v.3 license. As the license provider reports: "Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications, which include larger works using a licensed work, under the same license. Copyright and license notices must be preserved. Contributors provide an express grant of patent rights".
