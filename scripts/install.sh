#PLEASE NOTE THAT THIS INSTALLATION PAGE ASSUMES YOU ARE RUNNING ON A LINUX OR LINUX-LIKE OPERATING SYSTEM
#AND THAT YOU ALREADY HAVE A PYTHON (VERSION 3 OR FOLLOWING) ENVIRONMENT INSTALLED. 
#IF THIS IS NOT THE CASE, PLEASE CONSIDER TO ADJUST YOUR WORKING SPACE SO THAT IT COULD FIT THE REQUIREMENTS

#DOWNLOAD THIS FOLDER IN YOUR WORKING DIRECTORY
cd /path/to/working/directory
git init
git clone https://github.com/AstraBert/RREQUESTED/

#INSTALL DEPENDENCIES
## -EDLIB
python3 -m pip install edlib
## -PANDAS
python3 -m pip install pandas

#SET UP THE ALIAS MAKING SURE TO MODIFY THE PATH TO QUASARd.sh SO THAT IT MATCHES THE ONE IN WHICH YOU PUT IT
echo -n 'alias RREQ="bash /absolute/path/to/RREQUESTED.sh"' | cat >> ~/.bash_aliases
source ~/.bash_aliases

#TEST INSTALLATION
RREQ -h
