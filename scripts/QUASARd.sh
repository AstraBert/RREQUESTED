#!/bin/bash

#QUASARd stands for QUAlity Size and All Reads demultiplexing
# Usage function
usage() {
  echo "Usage: QUASARd -d, --directory RAW_READS_DIRECTORY [-dmx, --demultipex DEMULTIPLEXING_CSV] [-q,--quality QUALITY] [-mi,--min MINIMUM_LENGTH] [-ma,--max MAXIMUM_LENGTH]

  REQUIRED ARGUMENTS:
  -RAW_READS_DIRECTORY: Provide the path to the directory where raw reads are stored in fasta/fasta.gz/fastq/fastq.gz format (can accept also format mixture, but would be better to have all the file in the same format)

  OPTIONAL ARGUMENTS:
  -DEMULTIPLEXING_CSV: Provide a csv file (only comma-separated) whose fields must be "ID,FWD,REV", standing for sample id, forward primer and reverse primer. If not given, no demultiplexing will take place.
  -QUALITY: Provide the value of minimum average read quality for filtering (default: 7)
  -MINIMUM_LENGTH: Provide the minimum length of PCR product for size selection (default is 500)
  -MAXIMUM_LENGTH: Provide the maximum length of PCR product for size selection (default is 10000)

  
  Input QUASARd -h to show this message again"
  exit 1
}

#!/bin/bash

# Initialize variables with default values
directory=""
demultiplex=""
quality=7
min=""
max=""
source ~/.bash_aliases
alias_string=$(alias QUASARd)
command_part=$(echo "$alias_string" | awk -F"'" '{print $2}')
command_part=${command_part#bash }
absolute_path=$(realpath $(eval echo "$command_part"))
wd=$(dirname "$absolute_path")
echo $wd

# Loop through the arguments
while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -d|--directory)
            directory="$2"
            shift 2
            ;;
        -dmx|--demultiplex)
            demultiplex="$2"
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

if [ -z "$directory" ] && [ ! "$1" = "-h" ] && [ ! "$1" = "--help" ]; then
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


if [ ! -d "$directory" ] && [ ! "$1" = "-h" ] && [ ! "$1" = "--help" ]
then
    echo "Either $directory does not exist or it is not a directory"
    usage
else
    now=$(date +%Y-%m-%d' '%H:%M:%S)
    echo "Program started its execution at ${now}"
    count_fast=$(find $directory -type f -name "*.fast*" | wc -l)

    if [ $count_fast -eq 0 ]
    then
        echo "No fasta/fastq/fasta.gz/fastq.gz files, proceeding with the search of fastq"
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
                    python3 "${wd}/size_filter.py" \
                        -c ${file} \
                        -min ${min} \
                        -max ${max}
                    filename=$(basename "$file")
                    dirname=$(dirname "$file")
                    basename_no_ext="${filename%.*}"
                    combined_path="$dirname/$basename_no_ext"
                    if [ -z "$demultiplex" ]; then
                        echo "No demultiplexing table provided, the program ended current run"
                        now=$(date +%Y-%m-%d' '%H:%M:%S)
                        echo "Program ended current run at ${now}"
                    else
                        python3 "${wd}/demultiplex.py" \
                            -i ${combined_path}_size_selected.fasta \
                            -p ${demultiplex}
                        echo "Demultiplexing finished, the program ended current run"
                        now=$(date +%Y-%m-%d' '%H:%M:%S)
                        echo "Program ended current run at ${now}"
                    fi
                else
                    echo "Program started current run on file $file"
                    python3 "${wd}/quality_filter.py" \
                        -i ${file} \
                        -q ${quality}
                    python3 "${wd}/size_filter.py" \
                        -c ${file} \
                        -min ${min} \
                        -max ${max}
                    filename=$(basename "$file")
                    dirname=$(dirname "$file")
                    basename_no_ext="${filename%.*}"
                    combined_path="$dirname/$basename_no_ext"
                    if [ -z "$demultiplex" ]; then
                        echo "No demultiplexing table provided, the program ended current run"
                        now=$(date +%Y-%m-%d' '%H:%M:%S)
                        echo "Program ended current run at ${now}"
                    else
                        python3 "${wd}/demultiplex.py" \
                            -i ${combined_path}_size_selected.fastq \
                            -p ${demultiplex}
                        echo "Demultiplexing finished, the program ended current run"
                        now=$(date +%Y-%m-%d' '%H:%M:%S)
                        echo "Program ended current run at ${now}"
                    fi
                fi
            fi
        done
    fi
fi

now=$(date +%Y-%m-%d' '%H:%M:%S)
echo "Program ended its execution at ${now}"
