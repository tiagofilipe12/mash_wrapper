# mash sketch -p 7 test_run.fasta -o mash_screen_test

# gunzip *fastq.gz ... if gz is found

# mash screen ../../../plasmid_db_07202017/fasta/mash_screen_test.msh 5474_EP42_15_1_trimmed.fastq 5474_EP42_15_2_trimmed.fastq -p 7 > screen.tab

import os
#from subprocess import Popen, PIPE
import numpy as np
import json
from .shell_functions import shell_stdout_write

def mashscreen(ref_sketch, output_tag, threads, pvalue,
               winner, min_identity, reads):
    '''function to execute mash screen

    This function executes mash screen given a database in .msh generated by
    mash sketch for the reference sequences and takes reads in fastq (not in
    fastq.gz) to execute mash screen.

    Parameters
    ----------
    ref_sketch: str
    output_tag: str
    threads: str
    pvalue: str
    winner: bool
    min_identity: str
    reads: list or str

    Returns
    -------
    out_file: String
    '''

    out_file = os.path.join(output_tag + ".tab")

    # checks if there is two or one argument in reads
    # establishes default command
    if len(reads) > 1:
        command = ["mash",
                   "screen",
                   "-i",
                   min_identity,
                   "-v",
                   pvalue,
                   "-p",
                   threads,
                   ref_sketch,
                   reads[0],
                   reads[1]]
    # if single read file is provided...
    elif len(reads) == 1:
        command = ["mash",
                   "screen",
                   "-i",
                   min_identity,
                   "-v",
                   pvalue,
                   "-p",
                   threads,
                   ref_sketch,
                   reads[0]]
    else:
        print("no reads were provided. Please provide one to two read "
              "sequences in fastq.")

    # checks if winner takes it all is defined
    if winner:
        command = command.append("-w")

    # Runs mash screen
    print("\nRunning mash screen for reads...\n{}\n".format(" ".join(command)))
    shell_stdout_write(command, out_file)
    return out_file

def sort_mash_screen(out_file):
    '''sort mash screen output by identity

    Parameters
    ----------
    out_file: str
        this string stores th path to the file to be sorted

    Returns
    -------
    sorted_out_file: str

    '''
    sorted_out_file = " ".join(out_file.split(".")[:-1]) + "_sorted.tab"
    sorter_command = ["sort",
                      "-gr",
                      out_file]
    print("\nSorting {}...\n{}\n".format(out_file, " ".join(sorter_command)))
    shell_stdout_write(sorter_command, sorted_out_file)
    return sorted_out_file

def screen2json(mash_output):
    '''converts top results to json

    Parameters
    ----------
    mash_output: str
        this is a string that stores the path to this file

    Returns
    -------

    '''

    read_mash_output = open(mash_output)

    dic = {}
    median_list = []

    for line in read_mash_output:
        #print(line)
        tab_split = line.split("\t")
        identity = tab_split[0]
        #shared_hashes = tab_split[1]
        median_multiplicity = tab_split[2]
        #p_value = tab_split[3]
        query_id = tab_split[4]
        #query-comment should not exist here and it is irrelevant

        # here identity is what in fact interests to report to json but
        # median_multiplicity also is important since it gives an rough
        # estimation of the coverage depth for each plasmid.
        # Plasmids should have higher coverage depth due to their increased
        # copy number in relation to the chromosome.
        # TODO this only is valid for relative estimate of copy number (
        # TODO relative to all plasmids in results)
        # TODO basically a median of the coverage of all plasmids
        dic[query_id] = [identity, median_multiplicity]
        median_list.append(float(median_multiplicity))

    # median cutoff is twice the median of all median_multiplicity values
    # reported by mash screen. In the case of plasmids, since the database
    # has 9k entries and reads shouldn't have that many sequences it seems ok...
    # TODO but needs further testing
    # TODO maybe the cutoff can be ignored for filtering
    median_cutoff = np.median(median_list)*2

    output_json = open(" ".join(mash_output.split(".")[:-1]) + ".json", "w")

    filtered_dic = {}
    for k,v in dic.items():
        # estimated copy number
        copy_number = int(float(v[1])/median_cutoff)
        # assure that plasmid as at least twice the median coverage depth
        if float(v[1]) > median_cutoff: filtered_dic[k] = [v[0],
                                                           str(copy_number)]
    print("Exported dictionary has {} entries".format(len(filtered_dic)))
    output_json.write(json.dumps(filtered_dic))
    output_json.close()
