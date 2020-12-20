# ----
#
#  merge console output files according to time stamps
#
# usage:
#	python mergelog.py  INFILENAME N OUTFILENAME
#
#   merges INFILENAME.0 INFILENAME.1 ... INFILENAME.N-1 into OUTFILENAME 
#
#------

import sys
from heapq import *

infilename=sys.argv[1]
nfiles=int(sys.argv[2])
outfilename=sys.argv[3]

iflist=list()
msg_map=list()

# open all files
for i in range(nfiles):
	iflist.append(open(infilename+"."+str(i)))

outfile=open(outfilename,'w')

# initialize time -> message map by reading the 1st message from each file
for i in range(nfiles):
	msg=iflist[i].readline()
	t_msg=float(msg.split()[1])
	heappush(msg_map,(t_msg,msg,i))
	
# merge files
while (len(msg_map)>0):
	# pop msg with lowest time
	t,msg,idx=heappop(msg_map)
	# write it to output file
	outfile.write("<"+str(idx)+">"+msg)
	# replace with message from same file
	msg=iflist[idx].readline()
	if(len(msg)>0):
		t_msg=float(msg.split()[1])
		heappush(msg_map,(t_msg,msg,idx))
	
# close all input files
for i in range(nfiles):
	iflist[i].close()
	
outfile.close()
