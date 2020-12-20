#!/usr/bin/env python
# A Python script used to write or append a line to a file given
# the file name.
# Author: J. Rahardjo 2014

import sys
import os

def ensure_file_exists(filename = "output.txt"):
	""" Checks if the file to be appended to exists. If not,
		Create the file in the root location
	"""
	if not os.path.exists(filename):
		f = open(filename, 'w')
		f.write("# Test results file\n")
		f.close()

def append_to_file(string, filename = "output.txt"):
	ensure_file_exists(filename)
	f = open(filename, 'a')
	f.write(string + "\n")
	f.close()

def write_to_file(string, filename = "output.txt"):
	ensure_file_exists(filename)
	f = open(filename, 'w')
	f.write(string + "\n")
	f.close()	

if __name__ == '__main__':
	command = str(sys.argv[1])
	filename = str(sys.argv[2])
	string = str(sys.argv[3])
	if(command == 'a'):
		append_to_file(string, filename)
	elif(command == 'w'):
		write_to_file(string, filename)
	else:
		sys.exit("Command not recognized:\n" + \
			"--- Please use 'a' to append or 'w' to write")