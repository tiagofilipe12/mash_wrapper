# mash_wrapper.py

This script runs MASH in plasmid databases, using mash dist

Note: each header in fasta is considered a reference ('-i' option of mash).

### Options:

```
optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE [INPUTFILE ...], --input_references INPUTFILE [INPUTFILE ...]
                        Provide the input fasta files to parse. Not required
                        for multiple comparisons between assemblies.
  -rs REF_SKETCH, --reference_sketch REF_SKETCH
                        If you have a reference sketch for references provide
                        it with this option.
  -r READS [READS ...], --reads READS [READS ...]
                        Provide the input read files to parse. Usually fastq
                        files. This option is mutually with "-f".
  -f SEQUENCES [SEQUENCES ...], --sequences SEQUENCES [SEQUENCES ...]
                        Provide the input sequence files to parse. Usually
                        fasta files. This option is mutually exclusive with
                        "-r".
  -o OUTPUT_TAG, --output OUTPUT_TAG
                        Provide an output tag
  -t THREADS, --threads THREADS
                        Provide the number of threads to be used. Default: 1
  -no_rm, --no-remove   Specify if you do not want to remove the output
                        concatenated fasta.
  -j, --json            If you desire to export a json file with all
                        significant entries use this options.
  -ms, --mashix         Perform a matrix of all mash distance, taking all
                        files.

MASH related options:
  -k KMER_SIZE, --kmers KMER_SIZE
                        Provide the number of k-mers to be provided to mash
                        sketch. Default: 21
  -p PVALUE, --pvalue PVALUE
                        Provide the p-value to consider a distance
                        significant. Default: 0.05.
  -md MASHDISTANCE, --mashdist MASHDISTANCE
                        Provide the maximum mash distance to be parsed to the
                        matrix. Default: 0.1.

```

---

#### Outputs:


* A **text** file with MASH distance outputs as well as .msh files for each input read file.

* A **html** file with a bar plot graphical visualization, with each input read file as a series of data.

