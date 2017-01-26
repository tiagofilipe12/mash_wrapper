# mash_wrapper.py

This script runs MASH in plasmid databases.

Note: each header in fasta is considered a reference

###Options:

**'-i'**,**'--input_references'**, dest='inputfile', nargs='+', required=True, help='Provide the input fasta files to parse.'

**'-r'**,**'--reads'**, dest='reads', nargs='+', required=True, help='Provide the input read files to parse.'	

**'-o'**,**'--output'**, dest='output_tag', required=True, help='Provide an output tag'

**'-t'**, **'--threads'**, dest='threads', help='Provide the number of threads to be used'
