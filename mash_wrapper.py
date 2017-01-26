#!/usr/bin/env python

## Last update: 26/1/2017
## Author: T.F. Jesus
## This script runs MASH in plasmid databases.
## Note: each header in fasta is considered a reference

import argparse
import os
from subprocess import Popen, PIPE


def folderexist(directory):
	if not directory.endswith("/"):
		directory = directory + "/"
	if not os.path.exists(os.path.join(directory)):		
		os.makedirs(os.path.join(directory))
		print os.path.join(directory) + " does not exist. One will be created..."
	else:
		print os.path.join(directory) + " exists!"	

def header_fix(input_header):
	problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/","[","]",":","{","}"]
	for char in problematic_characters:
		input_header=input_header.replace(char, '_')
	return input_header

def master_fasta(fastas, output_tag):
	master_fasta = open("master_fasta_" + output_tag + ".fas", "w")
	for filename in fastas:
		fasta = open(filename,"r")
		for line in fasta:
			if line.startswith(">"):
				line = header_fix(line)
			master_fasta.write(line)

	master_fasta.close()
	return "master_fasta_" + output_tag + ".fas"

def sketch_references(inputfile, output_tag, threads):
	out_folder = os.path.join(os.path.dirname(os.path.abspath(inputfile)), output_tag)
	out_file = os.path.join(out_folder, output_tag +"_reference")
	folderexist(out_folder)
	sketcher_command = "mash sketch -o " + out_file + " -p "+ threads + " -i " + inputfile
	print
	print sketcher_command
	print
	p=Popen(sketcher_command, stdout = PIPE, stderr = PIPE, shell=True)
	p.wait()
	return out_file + ".msh"

## According to mash turorial it is useful to provide the -m 2 option in order to remove single-copy k-mers
def sketch_reads(read, inputfile, output_tag, threads):
	out_folder = os.path.join(os.path.dirname(os.path.abspath(inputfile)), output_tag)
	out_file = os.path.join(out_folder, output_tag + "_" + os.path.basename(read).split(".")[0]) 
	folderexist(out_folder)
	sketcher_command = "mash sketch -o " + out_file + " -p "+ threads +" -m 2 -r " + read
	print
	print sketcher_command
	print
	p=Popen(sketcher_command, stdout = PIPE, stderr = PIPE, shell=True)
	p.wait()
	return out_file + ".msh"

def masher(ref_sketch, read_sketch, output_tag, threads):
	out_folder = os.path.join(os.path.dirname(os.path.abspath(ref_sketch)))
	out_file = os.path.join(out_folder , output_tag)
	mash_command = "mash dist -p "+ threads+ " " + ref_sketch +" "+ read_sketch +" > " + os.path.basename(read_sketch).split(".")[0]+"_distances.txt"
	print
	print mash_command
	print	
	p=Popen(mash_command, stdout = PIPE, stderr = PIPE, shell=True)
	p.wait()

def main():
	parser = argparse.ArgumentParser(description="Retrieves all gb files given an input fasta")
	parser.add_argument('-i','--input_references', dest='inputfile', nargs='+', required=True, help='Provide the input fasta files to parse.')
	parser.add_argument('-r','--reads', dest='reads', nargs='+', required=True, help='Provide the input read files to parse.')	
	parser.add_argument('-o','--output', dest='output_tag', required=True, help='Provide an output tag')
	parser.add_argument('-t', '--threads', dest='threads', help='Provide the number of threads to be used')
	args = parser.parse_args()
	if args.threads is None:
		threads = "1"
	else:
		threads = args.threads
	fastas = []
	for filename in args.inputfile:
		if any (x in filename for x in [".fas",".fasta",".fna",".fsa", ".fa"]):
			fastas.append(filename)

	main_fasta = master_fasta(fastas, args.output_tag)
	ref_sketch=sketch_references(main_fasta,args.output_tag,threads)
	for read in args.reads:
		read_sketch = sketch_reads(read, main_fasta, args.output_tag, threads)
		mash_output = masher(ref_sketch, read_sketch, args.output_tag, threads)

	## remove master_fasta

if __name__ == "__main__":
	main()