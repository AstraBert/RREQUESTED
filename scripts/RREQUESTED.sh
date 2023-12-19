#!/bin/bash

# Usage function
usage() {
  echo "Usage: RREQ -d, --directory RAW_READS_DIRECTORY [-q,--quality QUALITY] [-mi,--min MINIMUM_LENGTH] [-ma,--max MAXIMUM_LENGTH] [-v, --version]

  REQUIRED ARGUMENTS:
  -RAW_READS_DIRECTORY: Provide the path to the directory where raw reads are stored in fasta/fasta.gz/fastq/fastq.gz format (can accept also format mixture, but would be better to have all the file in the same format)

  OPTIONAL ARGUMENTS:
  -QUALITY: Provide the value of minimum average read quality for filtering (default: 7)
  -MINIMUM_LENGTH: Provide the minimum length of PCR product for size selection (default is 500)
  -MAXIMUM_LENGTH: Provide the maximum length of PCR product for size selection (default is 10000)
  -VERSION: Print the version of the code
  
  Input RREQ -h to show this message again"
  exit 1
}

ConPath=$(which conda)
tmp=${ConPath#* }
Conda=${tmp%%/bin/co*}


# Initialize variables with default values
directory=""
quality=7
min=""
max=""
version="1.1.0"
sourcedir=$(dirname $0)
wd=$(realpath "$sourcedir")


# Loop through the arguments
while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -v|--version)
            echo Version: $version
            exit 0
            ;;
        -d|--directory)
            directory="$2"
            shift 2
            ;;
        -q|--quality)
            quality="$2"
            shift 2
            ;;
        -mi|--min)
            min="$2"
            shift 2
            ;;
        -ma|--max)
            max="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    usage
fi

if [ -z "$directory" ] && [ ! "$1" = "-h" ] && [ ! "$1" = "--help" ] && [ ! "$1" = "-v" ] && [ ! "$1" = "--version" ]; then
    echo "Missing required argument: RAW_READS_DIRECTORY"
    usage
fi



if [ -z "$quality" ]; then
  quality=7
fi

if [ -z "$min" ]; then
  min=500
fi

if [ -z "$max" ]; then
  max=10000
fi


if [ ! -d "$directory" ] && [ ! "$1" = "-h" ] && [ ! "$1" = "--help" ] && [ ! "$1" = "-v" ] && [ ! "$1" = "--version" ]
then
    echo "Either $directory does not exist or it is not a directory"
    usage
else
    source ${Conda}/etc/profile.d/conda.sh
    conda activate ${wd}/environments/rrequested
    now=$(date +%Y-%m-%d' '%H:%M:%S)
    echo "Program started its execution at ${now}"
    count_fast=$(find $directory -type f -name "*.fast*" | wc -l)
    if [ $count_fast -eq 0 ]
    then
        echo "No fasta/fastq/fasta.gz/fastq.gz files, quitting"
        conda deactivate
        exit 1
    else
        mkdir -p ${directory}/results
        echo "Copying files..."
        for file in $directory/*.fast*
        do
            flname=$(basename "$file")
            drname=$(dirname "$file")
            bsname_no_ext="${flname%.*}"
            nogz="$drname/$bsname_no_ext"
            echo $nogz
            if [[ $file == *.fast?.gz ]]
            then
                pigz -d $file 
                cp -n ${nogz} ${directory}/results
            else
                cp -n $file ${directory}/results 
            fi
        done
        echo "Files copied, starting with the analysis"
        for file in $directory/results/*.fast*
        do
            if [[ ${file} == *_size_selected.fast? ]]
            then
                continue
            else
                if [[ ${file} == *.fasta* ]]
                then
                    echo "Program started current run on file $file"
                    python3 ${wd}/size_filter.py -c ${file} -min ${min} -max ${max}
                    filename=$(basename "$file")
                    dirname=$(dirname "$file")
                    basename_no_ext="${filename%.*}"
                    combined_path="$dirname/$basename_no_ext"
                    python3 "${wd}/unref_demult.py" -i ${combined_path}_size_selected.fasta
                    echo "Demultiplexing finished, the program ended current run"
                    now=$(date +%Y-%m-%d' '%H:%M:%S)
                    echo "Program ended current run at ${now}"
                else
                    echo "Program started current run on file $file"
                    python3 "${wd}/quality_filter.py" -i ${file} -q ${quality}
                    python3 "${wd}/size_filter.py" -c ${file} -min ${min} -max ${max}
                    filename=$(basename "$file")
                    dirname=$(dirname "$file")
                    basename_no_ext="${filename%.*}"
                    combined_path="$dirname/$basename_no_ext"
                    python3 "${wd}/unref_demult.py" -i ${combined_path}_size_selected.fastq
                    echo "Demultiplexing finished, the program ended current run"
                    now=$(date +%Y-%m-%d' '%H:%M:%S)
                    echo "Program ended current run at ${now}"
                fi
            fi
        done
    fi
fi

now=$(date +%Y-%m-%d' '%H:%M:%S)
echo "Program ended its execution at ${now}"
conda deactivate
