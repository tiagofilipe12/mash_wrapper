#!/usr/bin/env python

## Last update: 1/2/2017
## Author: T.F. Jesus
## This script runs MASH in plasmid databases.
## Note: each header in fasta is considered a reference

import argparse
import os
import operator
from subprocess import Popen, PIPE
import plotly
import plotly.graph_objs as go


## Checks if a directory exists and if not creates one.
def folderexist(directory):
	if not directory.endswith("/"):
		directory = directory + "/"
	if not os.path.exists(os.path.join(directory)):		
		os.makedirs(os.path.join(directory))
		print os.path.join(directory) + " does not exist. One will be created..."
	else:
		print os.path.join(directory) + " exists!"	

## Function to fix several issues that fasta header names can have with some programs 
def header_fix(input_header):
	problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/","[","]",":","{","}"]
	for char in problematic_characters:
		input_header=input_header.replace(char, '_')
	return input_header

## Function to create a master fasta file from several fasta databases. One fasta is enought though
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

## Makes the sketch command of mash for the reference
def sketch_references(inputfile, output_tag, threads,kmer_size):
	out_folder = os.path.join(os.path.dirname(os.path.abspath(inputfile)), output_tag)
	out_file = os.path.join(out_folder, output_tag +"_reference")
	folderexist(out_folder)
	sketcher_command = "mash sketch -o " + out_file +" -k " + kmer_size+ " -p "+ threads + " -i " + inputfile
	print
	print sketcher_command
	print
	p=Popen(sketcher_command, stdout = PIPE, stderr = PIPE, shell=True)
	p.wait()
	stdout,stderr= p.communicate()
	print stderr
	return out_file + ".msh"

## Makes the sketch command of mash for the reads to be compare to the reference.
## According to mash turorial it is useful to provide the -m 2 option in order to remove single-copy k-mers
def sketch_reads(read, inputfile, output_tag, threads,kmer_size):
	out_folder = os.path.join(os.path.dirname(os.path.abspath(inputfile)), output_tag)
	out_file = os.path.join(out_folder, output_tag + "_" + os.path.basename(read).split(".")[0]) 
	folderexist(out_folder)
	sketcher_command = "mash sketch -o " + out_file +" -k " + kmer_size + " -p "+ threads +" -m 2 -r " + read
	print
	print sketcher_command
	print
	p=Popen(sketcher_command, stdout = PIPE, stderr = PIPE, shell=True)
	p.wait()
	stdout,stderr= p.communicate()
	print stderr
	return out_file + ".msh"

## Executes mash dist
def masher(ref_sketch, read_sketch, output_tag, threads):
	out_folder = os.path.join(os.path.dirname(os.path.abspath(ref_sketch)))
	out_file = os.path.basename(read_sketch).split(".")[0]+"_distances.txt"
	out_file_path = os.path.join(out_folder + "/" + out_file)
	mash_command = "mash dist -p "+ threads+ " " + ref_sketch +" "+ read_sketch +" > " + out_file_path
	print
	print mash_command
	print	
	p=Popen(mash_command, stdout = PIPE, stderr = PIPE, shell=True)
	p.wait()
	stdout,stderr= p.communicate()
	print stderr
	return out_file_path

## Reads the output of mash dist and performes a barplot for each reads
def mashdist2graph(list_mash_files, tag):
	dist_dict = {}
	trace_list=[]
	## First reads mash output files and creates a sorted dictionary for each file by mash distance (from highest to lowest values)
	for in_file in list_mash_files:
		dist_file = open(in_file, 'r')
		for line in dist_file:
			tab_split = line.split("\t")
			reference_id = tab_split[0]
			query_id = tab_split[1].split("/")[-1].split(".")[0]
			mash_distance = tab_split[2]
			p_value = tab_split[3]			
			if float(p_value) < 0.05:		## filters distances with no significant p-values
				dist_dict[reference_id] = mash_distance
		## Orders the dictionary from the minor distances to the major
		sorted_dist_dict = sorted(dist_dict.items(), key=operator.itemgetter(1))
		## creates two lists to be given as the x and y for the bar plot
		global_reference = []
		global_values = []		
		for k,v in sorted_dist_dict:
			if k not in global_reference:
				global_reference.append(k)
			global_values.append(1-float(v))		#converts mash distances to mash "similarities", i.e., inverts the scale.
		trace = go.Bar(x=global_reference, y=global_values, name=query_id)
		trace_list.append(trace)

	## Creates the plot itself		
	layout = go.Layout(barmode='group', yaxis= dict(title='Mash distances (mutation rate between pairwise comparisons)'))
	fig = go.Figure(data=trace_list, layout=layout)
	plot_url = plotly.offline.plot(fig, filename= tag + '.html',auto_open=False)	

def main():
	parser = argparse.ArgumentParser(description="Retrieves all gb files given an input fasta")
	parser.add_argument('-i','--input_references', dest='inputfile', nargs='+', required=True, help='Provide the input fasta files to parse.')
	parser.add_argument('-r','--reads', dest='reads', nargs='+', required=True, help='Provide the input read files to parse.')	## should implement a parser for a given directory with reads or a list file with all full path to each read library
	parser.add_argument('-o','--output', dest='output_tag', required=True, help='Provide an output tag')
	parser.add_argument('-t', '--threads', dest='threads', help='Provide the number of threads to be used. Default: 1')
	parser.add_argument('-k', '--kmers', dest='kmer_size', help='Provide the number of k-mers to be provided to mash sketch. Default: 21')
	parser.add_argument('-no_rm', '--no-remove', dest='no_remove', action='store_true', help='Specify if you do not want to remove the output concatenated fasta.')
	args = parser.parse_args()
	if args.threads is None:
		threads = "1"
	else:
		threads = args.threads
	if args.kmer_size is None:
		kmer_size = "21"
	else:
		kmer_size = args.kmer_size
	fastas = []
	for filename in args.inputfile:
		if any (x in filename for x in [".fas",".fasta",".fna",".fsa", ".fa"]):
			fastas.append(filename)

	main_fasta = master_fasta(fastas, args.output_tag)
	ref_sketch=sketch_references(main_fasta,args.output_tag,threads,kmer_size)
	list_mash_files=[]
	for read in args.reads:
		read_sketch = sketch_reads(read, main_fasta, args.output_tag, threads,kmer_size)
		mash_output = masher(ref_sketch, read_sketch, args.output_tag, threads)
		list_mash_files.append(mash_output)
	mashdist2graph(list_mash_files, args.output_tag)

	## remove master_fasta
	if not args.no_remove:
		os.remove(main_fasta)

if __name__ == "__main__":
	main()