## a function to convert mash screen results to json for importing with
#pATLAS

# mash sketch -p 7 test_run.fasta -o mash_screen_test

# gunzip *fastq.gz ... if gz is found

# mash screen ../../../plasmid_db_07202017/fasta/mash_screen_test.msh 5474_EP42_15_1_trimmed.fastq 5474_EP42_15_2_trimmed.fastq -p 7 > screen.tab