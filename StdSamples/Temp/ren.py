import os,sys
import glob
dirs = glob.glob(sys.argv[1])
for f in dirs:
	s = f.replace(sys.argv[2],sys.argv[3])
	print(sys.argv[1])
	os.system("mv " + f + " " + s)
