Phylogenetic Intersection Analysis
==================================
Metagenomic phylogenetic assignment of mixed environmental assemblages
Allaby lab, University of Warwick
Version 5.6
2021-05-14

Phylogenetic intersection analysis (PIA) takes standard-format BLAST output and a corresponding FASTA file. It assigns reads to phylogenetic intersections based on their BLAST hits, assuming that the true taxon will be inside that phylogenetic intersection. It is designed to be robust to the uneven representation of taxa in databases.

A very early version was published in Smith et al. ("Sedimentary DNA from a submerged site reveals wheat in the British Isles 8000 years ago.", Science, 2015) and revised version in Cribdon et al. ("PIA: More Accurate Taxonomic Assignment of Metagenomic Data Demonstrated on sedaDNA From the North Sea", Frontiers in Ecology and Evolution, 2020).

Test FASTAs, their corresponding BLAST files, and expected PIA output directories (based on the BLAST nt database and taxdump files from 2021-05-03) are included.

Please email r.g.allaby@warwick.ac.uk with any questions, problems, or suggestions. Thanks for using PIA!


Prerequisites
-------------
-   Perl 5
-   Perl module List::MoreUtils
-   Perl module DB_File


BLASTing your input FASTA
-------------------------
-   PIA interprets BLAST output from your query FASTA.
-   This is the basic recommended BLAST command:
    blastn -db [database] -num_threads [n] -query [input FASTA] -out [output]  -max_target_seqs 500 -outfmt "6 std staxids"
-   -max_target_seqs 500 returns up to the first ~500 hits for each read. PIA only considers a finite amount that is probably some way below 500. 500 is a safe bet for now.
-   -outfmt "6 std staxids" is necessary. However, extra optional fields can be included between the standard ones and staxids, e.g. "6 std staxid ssciname staxids". PIA will ignore these fields.


Setting up PIA
--------------
-   PIA requires names.dmp and nodes.dmp from the NCBI taxdump. These should be from a similar date as your BLAST database to maximise compatibility (e.g. taxa are added and taxonomic IDs can change over time).
-   To download the files, go to https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/ and download taxdump.tar.gz. Uncompress. If you need older versions, look here: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/
-   Move your names.dmp and nodes.dmp to Reference_files/.
-   Have all PIA files and directories in the same directory (see tree below). Input FASTAs and BLAST files can be in other directories.
-   Note that you should be able to run PIA on multiple FASTA files in the same directory simultaneously, but this has not been thoroughly tested.

Reference_files/
    names.dmp
    nodes.dmp
PIA.pl
PIA_index_maker.pl
PIA_inner.pl


Before running PIA, use PIA_index_maker.pl to build DBM index files from names.dmp and nodes.dmp:
>perl PIA_index_maker.pl
This will take a few minutes. You will need to make new index files if you update the .dmp files (do this when you update the BLAST database) or if you want to run the PIA on a new machine. The index files do not transfer between machines (something to do with the DB_File module).

If you get an error like this:
"Cannot open Reference_files/nodes.dmp.dbm_test0.fasta.read_infoaaa: Inappropriate file type or format"
Try making DBM files on your current machine.


Running PIA
-----------
PIA_inner.pl does the analysis. PIA.pl is a wrapper that allows threading. To run PIA.pl:
>perl PIA.pl -f [input FASTA] -b [input BLAST file] -t [optional number of threads] [optional options for advanced users]

While PIA.pl is running, each thread will generate a PIA_inner.log.txt. Look at that file to see how far the thread has progressed through its list of sequences.


Outputs
-------
Stored in an output directory derived from the input FASTA: [input].PIA_output/.
Note that output files will be produced even if the input FASTA was empty.

-   [input].Full.txt: a tab-delimited table listing information for each read in the FASTA. It should be used for investigating how PIA worked, rather than as a final output. Besides 'Taxonomic diversity', all hit information refers to the final list of processed hits rather than the raw BLAST output.
   -   Read_name: From the FASTA and BLAST input files.
   -   Coverage: Match length / read length. Reads with insufficient coverage are not processed further.
   -   Top_hit_name, Top_hit_ID: Taxon name and ID of the top BLAST hit. Used to calculate phylogenetic intersection and phylogenetic range.
   -   Expect_value: Expect (E) value of the top hit. For interest only.
   -   Identities_value: Identities value of the top hit; the percentage of identical matches. For interest only.
   -   Next_hit_name, Next_hit_ID: Taxon name and ID of the second-best hit. Used to calculate phylogenetic intersection.
   -   Last_hit_name, Last_hit_ID: Taxon name and ID of the last hit. Used to calculate phylogenetic range.
   -   Phylogenetic_range_name, Phylogenetic_range_ID: Taxon name and ID of the phylogenetic intersection of the top and last hits. For interest only.
   -   Raw_hit_count: How many hits were in the BLAST file. For interest only.
   -   Taxonomic_diversity_(up_to_cap_if_met): Number of taxa in the raw BLAST hits, or cap if met. Used to calculate taxonomic diversity score.
   -   Taxonomic_diversity_score: (taxonomic diversity / cap) - (1/cap). A read is added to the Summary Basic/Reads files if this is above a threshold (default 0.1).
   -   Phylogenetic_intersection_name, Phylogenetic_intersection_ID: Taxon name and ID of the phylogenetic intersection of the top and next hits; what a read is assigned to.

-   [input].PIA_inner_logs.txt: Collected logs from PIA_inner.pl. Notes BLAST hits that had trouble being identified. Lots of these suggest that you might want to update/synchronise your BLAST database and NCBI taxdump files.

-   [input].Post-PIA.fasta: A subset of the input FASTA containing only reads that passed PIA. Each read name has "_IDx" appended, where x is the taxonomic ID the read was assigned to.

-   [input].Summary_Basic.txt: The main PIA output. Lists the IDs, names, and counts of taxa that reads passing PIA were assigned to, excluding 'root'. Taxa are listed alphabetically. The header states run parameters and total read count.

-   [input].Summary_Reads.txt: Lists reads that passed PIA and which taxa they were assigned to.

-   [input].Summary_Reads_MEGAN.csv: Same information as the Summary Reads file, but in a CSV format that MEGAN can read. See "How to open Summary_Reads_MEGAN.csv files in MEGAN" below.


How to open Summary_Reads_MEGAN.csv files in MEGAN
--------------------------------------------------
The CSVs contain each read name, its taxonomic ID, and a stand-in bitscore of 50 (the default minimum pass score for MEGAN's LCA).

To open a CSV in MEGAN, go to go to "File -> Import -> Text (CSV) Format...". Keep the default format (Read Class) and separator (Comma). Tick the classification "Taxonomy" and, if available, tick "Parse accessions ids" and accept the default taxonomy analysis settings.

MEGAN will then try to run the LCA. To turn the LCA off, keep the min score at 50, change top percent to 0.001, change min support percent to 0, and keep min support at 1.

Remember to uncollapse the tree to view it fully.


Known issues
------------
-   If you run PIA.pl with multiple threads, each thread makes a copy of the DBM index files to avoid conflict. These files can be big, so watch your memory usage.
-   BLAST hits are identified using their taxonomic ID. These IDs apparently change sometimes. If your BLAST database and NCBI reference files are not well synchronised (e.g. one dates from years ago and the other from yesterday), the PIA may not be able to calculate intersections as accurately and hits may be unfairly excluded.
