import os,sys

f = open(sys.argv[1],'r')
for line in f:
    os.system('git submodule add --force ' + line)
