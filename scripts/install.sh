##INSTALL DEPENDENCIES
python3 -m pip install edlib
python3 -m pip install pandas

##SET UP THE ALIAS MAKING SURE TO MODIFY THE PATH TO QUASARd.sh SO THAT IT MATCHES THE ONE IN WHICH YOU PUT IT

echo -n 'alias QUASARd="bash /absolute/path/to/QUASARd.sh"' | cat >> ~/.bash_aliases
source ~/.bash_aliases

##TEST INSTALLATION

QUASARd -h