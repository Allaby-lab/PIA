PROVISIONAL - THIS IS A BETA VERSION

The Phylogenetic Intersection Analysis
======================================
Metagenomic phylogenetic assignment of mixed environmental assemblages
Allaby lab, University of Warwick
Version 5.0
2019-10-11

Phylogenetic intersection analysis (PIA) takes standard-format BLAST output and a corresponding FASTA file. It assigns reads to phylogenetic intersections based on their BLAST hits, assuming that the true taxon will be inside that phylogenetic intersection. It is designed to be robust to the uneven representation of taxa in databases.

A very early version was published in Smith et al. ("Sedimentary DNA from a submerged site reveals wheat in the British Isles 8000 years ago.", Science, 2015) and a more current version will be in a forthcoming paper this year. For more information, email r.cribdon@warwick.ac.uk.


Prerequisites
-------------
-   Perl 5
-   Perl module List::MoreUtils
-   Perl module DB_File


BLASTing your input FASTA
-------------------------
-   PIA interprets BLAST output from your query FASTA.
-   Read names must be identical in the FASTA and BLAST files. If the read names in your FASTA contains whitespace, the BLAST file will truncate them and they will no longer be identical.
-   This is the basic recommended BLAST command:
    blastn -db [database] -num_threads [n] -query [input FASTA] -out [output]  -max_target_seqs 500 -outfmt "6 std staxids"
-   -max_target_seqs 500 returns up to the first ~500 hits for each read. PIA only considers a finite amount that is probably some way below 500. 500 is a safe bet for now.
-   -outfmt "6 std staxids" is necessary. However, extra optional fields can be included between the standard ones and staxids, e.g. "6 std staxid ssciname staxids". PIA will ignore these fields.


Setting up PIA
--------------
-   Go to https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/ and download taxdump.tar.gz. Uncompress.
-   Move names.dmp and nodes.dmp to Reference_files/.
-   For Collate_summary_basics.pl only: move fullnamelineage.dmp to Reference_files/.
-   Have all PIA files and directories in the same directory (see tree below). Note that you should be able to run PIA on multiple FASTA files in the same directory simultaneously, but this has not been thoroughly tested.

Reference_files/
    names.dmp
    nodes.dmp
    fullnamelineage.dmp (for Collate_summary_basics.pl only)
(Collate_summary_basics.pl)
PIA.pl
PIA_index_maker.pl
PIA_inner.pl
[input FASTA]
[input BLAST]


Usage
------
-   Run PIA_index_maker.pl:
    >perl PIA_index_maker.pl
    This generates DBM index files for names.dmp and nodes.dmp. It will take a few minutes. You will need to make new index files if you change the .dmp files or if you want to run the PIA on a new machine. The index files do not transfer between machines (something to do with DB_File).
-   PIA_inner.pl does the analysis. PIA.pl is a wrapper that allows threading. To run PIA.pl:
    >perl PIA.pl -f [input FASTA] -b [input BLAST file] -t [number of threads] [optional options for advanced users]
    While PIA.pl is running, each thread will generate a PIA_inner.log.txt. Look at that file to see how far the thread has progressed through its list of sequences.
-   OPTIONAL: collate Summary_Basic.txt output from multiple FASTAs using Collate_summary_basics.pl as follows:
    >perl Collate_summary_basics.pl [optional options] [summary basics]

Please report any problems to r.cribdon@warwick.ac.uk.


Outputs
-------
Will be in [input FASTA].header_out/.
-   [FASTA].header_out.intersects.txt_Summary_Basic.txt: the main PIA output. Summarises any reads that passed both the quality and taxonomic diversity filters, excluding reads assigned to 'none' or 'root'. Reads are grouped by taxon for easy interpretation. The header states run parameters.
-   [FASTA].header_PIA_inner_logs.txt: collected logs from PIA_inner.pl. Notes BLAST hits that had trouble being identified. Lots of these suggest that you might want to update/synchronise your BLAST database and NCBI reference files.
-   [FASTA].header_out.intersects.txt: lists information for each read that passed the initial quality filter. Apart from 'raw hit count' and 'taxonomic diversity', all hit information refers to the final list of processed hits.
   -   'Query': read name. From the FASTA and BLAST input files.
   -   'top hit': taxon name and ID of the top hit. Used to calculate classification intersect and phylogenetic range.
   -   'expect': expect (E) value of the top hit. For interest only.
   -   'identities': Identities value of the top hit; the percentage of identical matches. For interest only.
   -   'next hit': taxon name and ID of the second-best hit. Used to calculate classification intersect.
   -   'last hit': taxon name and ID of the last hit in the final list. Used to calculate phylogenetic range.
   -   'taxon count': how many taxa were in the final list. Also how many hits were in the final list because it's one hit per taxon. For interest only.
   -   'phylogenetic range': phylogenetic intersection of the top and last hits. For interest only.
   -   'raw hit count': how many hits were in the BLAST file. For interest only.
   -   'taxonomic diversity (up to cap if met)': number of taxa in the BLAST file, or cap if met. Equal to the number of hits analysed because it's one hit per taxon. Used to calculate taxonomic diversity score.
   -   'taxonomic diversity score': (taxonomic diversity / cap) - (1/cap). Must be above a threshold (default 0.1) for a read to make it to the summary basic.
   -   'classification intersect': phylogenetic intersection of the top and next hits. This is what the read is assigned to.



Known issues
------------
-   If you run PIA.pl with multiple threads, each thread makes a copy of the DBM index files to avoid conflict. These files can be big, so watch your memory usage.
-   BLAST hits are identified using their taxonomic ID. These IDs apparently change sometimes. If your BLAST database and NCBI reference files are not well synchronised (e.g. one dates from years ago and the other from yesterday), the PIA may not be able to calculate intersections as accurately and hits may be unfairly excluded.
-   If you have to interrupt a run, delete any output. It may confuse future runs.
