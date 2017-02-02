# mash_wrapper.py

This script runs MASH in plasmid databases.

Note: each header in fasta is considered a reference ('-i' option of mash).

###Options:

**'-i'**,**'--input_references'**, dest='inputfile', nargs='+', required=True, help='Provide the input fasta files to parse.'

**'-r'**,**'--reads'**, dest='reads', nargs='+', required=True, help='Provide the input read files to parse.'	

**'-o'**,**'--output'**, dest='output_tag', required=True, help='Provide an output tag'

**'-t'**, **'--threads'**, dest='threads', help='Provide the number of threads to be used'

**'-k'**,**'--kmers'**, dest='kmer_size', help='Provide the number of k-mers to be provided to mash sketch. Default: 21'

**'-no_rm'**, **'--no-remove'**, dest='no_remove', action='store_true', help='Specify if you do not want to remove the output concatenated fasta.'


---

####Outputs:

* A **text** file with MASH distance outputs as well as .msh files for each input read file.

* A **html** file with a bar plot graphical visualization, with each input read file as a series of data.

