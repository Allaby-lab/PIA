PROVISIONAL - THIS IS A BETA VERSION

The Phylogenetic Intersection Analysis
======================================
Metagenomic phylogenetic assignment of mixed environmental assemblages
Allaby lab, University of Warwick
Version 4.8
2019-09-10

The phylogenetic intersection analysis (PIA) takes standard-format BLAST output and a corresponding FASTA file. It assigns reads to phylogenetic intersections based on their BLAST hits, assuming that the true taxon will be inside that phylogenetic intersection. It is designed to be robust to the uneven representation of taxa in databases.

A very early version was published in Smith et al. ("Sedimentary DNA from a submerged site reveals wheat in the British Isles 8000 years ago.", Science, 2015) and a more current version will be in a forthcoming paper this year. For more information, email r.cribdon@warwick.ac.uk.


Prerequisites
-------------
-   Perl 5
-   Perl module List::MoreUtils
-   Perl module DB_File
-   Add to Reference_files/ names.dmp and nodes.dmp from the NCBI: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
-   For Collate_summary_basics.pl only: add to Reference_files/ fullnamelineage.dmp from the NCBI


Files and directories
---------------------
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
-   The input FASTA and corresponding input BLAST must have the same headers (sequence names). If the FASTA headers are very long, BLAST may crop them. Change header names if necessary to prevent this.
-   BLAST your data with this format option: -outfmt "6 std staxid".
-   Have all PIA files and directories in the same directory. Note that you should be able to run the PIA on multiple FASTA files in the same directory simultaneously, but this has not been thoroughly tested.
-   Run PIA_index_maker.pl. This generates DBM index files for names.dmp and nodes.dmp. You will need to make new index files every time you change these files and when you want to run the PIA on a new machine. The index files do not reliably transfer between machines (something to do with DB_File).
-   PIA_inner.pl does the analysis. PIA.pl is a wrapper that allows threading. To run PIA.pl:
    >perl PIA.pl -f [input FASTA] -b [input BLAST file] -t [number of threads] [other options]
-   To run PIA_inner.pl alone, see notes at the top of PIA_inner.pl.
-   OPTIONAL: collate Summary_Basic.txt output from multiple FASTAs using Collate_summary_basics.pl as follows:
    >perl Collate_summary_basics.pl -f Reference_files/fullnamelineage.dmp [other options] [summary basics]

Please report any problems to r.cribdon@warwick.ac.uk.


Outputs
-------
Will be in [input FASTA]_out/.
-   [FASTA].header_out.intersects.txt: lists information for each read that passed the initial quality filter.
-   [FASTA].header_out.intersects.txt_Summary_Basic.txt: the main PIA output. Summarises any reads that passed both the quality and taxonomic diversity filters. Reads are grouped by taxon for easy interpretation. Also contains a header section that states run parameters.
-   [FASTA].header_PIA_inner_logs.txt: collected logs from PIA_inner.pl. Notes BLAST hits that had trouble being identified. Lots of these suggest that you might want to update/synchronise your BLAST database and NCBI reference files.


Known issues
------------
-   If you run PIA.pl with multiple threads, each thread makes a copy of the DBM index files to avoid conflict. These files can be big, so watch your memory usage.
-   BLAST hits are identified using their taxonomic ID. These IDs apparently change sometimes. If your BLAST database and NCBI reference files are not well synchronised (e.g. one dates from years ago and the other from yesterday), the PIA may not be able to calculate intersections as accurately and hits may be unfairly excluded.
-   If you have to interrupt a run, delete any output. It will confuse future runs.
