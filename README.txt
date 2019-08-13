PROVISIONAL - THIS IS A BETA VERSION

The Phylogenetic Intersection Analysis
======================================
Metagenomic phylogenetic assignment of mixed environmental assemblages
Allaby lab, University of Warwick
Version 4.4
2019-08-13

The phylogenetic intersection analysis (PIA) takes standard-format BLAST output and a corresponding FASTA file. It assigns reads to phylogenetic intersections based on their BLAST hits, assuming that the true taxon will be inside that phylogenetic intersection. It is designed to be robust to the uneven representation of taxa in databases.

A very early version was published in Smith et al. ("Sedimentary DNA from a submerged site reveals wheat in the British Isles 8000 years ago.", Science, 2015) and a more current version will be in a forthcoming paper this year. For more information, email r.cribdon@warwick.ac.uk.


Prerequisites
-------------
-   Perl 5
-   Perl module List::MoreUtils
-   Perl module DB_File
-   Directory called Reference_files/ containing names.dmp and nodes.dmp from the NCBI: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/


Files and directories
---------------------
Modules/
    FileChecks.pm
    FileManipulations.pm
    FileMerge.pm
    TreeOfLife.pm
Reference_files/
    names.dmp
    nodes.dmp
Helpfile_PIA.txt
PIA.pl
PIA_inner.pl
[input FASTA]
[input BLAST]


Usage
------
-   Have all PIA scripts and directories in the same directory. Note that you should be able to run the PIA on multiple FASTA files in the same directory simultaneously, but this has not been thoroughly tested.
-   PIA_inner.pl does the analysis. PIA.pl is a wrapper that allows threading. To run PIA.pl:
    >perl PIA.pl -f [input FASTA] -b [input BLAST file] -p [taxonomic ID of expected phylogenetic range] -t [number of threads] [other options]
-   The input FASTA and corresponding input BLAST must have the same headers (sequence names). If the FASTA headers are very long, BLAST may crop them. Change header names if necessary to prevent this.
-   To run PIA_inner.pl alone, see notes at the top of PIA_inner.pl.
-   Outputs will be in [input FASTA]_out/.
    -   The intersects file is the full output.
    -   The summary basic file is a simple list of how many reads were assigned to which phylogenetic intersections.
    -   The log tracks which reads (headers) have been processed. Use this to pick up where you left off if the PIA is interrupted. It also states when a BLAST hit could not be matched to names.dmp.
    -   The timer states how long PIA.pl took to run.

Please report any problems to r.cribdon@warwick.ac.uk.


Outputs
-------
-   [FASTA].header_out.intersects.txt: the main PIA output. Lists metrics for each read that passed the initial quality filter and its BLAST hits.
-   [FASTA].header_out.intersects.txt_Summary_Basic.txt: summarises any reads that passed every filter. Reads are grouped by taxon for easy interpretation.
-   [FASTA].header_PIA_inner_logs.txt: collected logs from PIA_inner.pl. Notes any BLAST hits that could not be matched to the names file.
-   timer_[FASTA].txt: brief log written by PIA.pl, if using. Times the run.


Known issues
------------
-   The extended summary function is untested.
-   BLAST hits are assigned to taxa using names.dmp. If there are matches to synonymous taxa, such as Iris the angiosperm and Iris the mantid, the PIA will exclude matches that do not fall in the expected phylogenetic range and then choose the match with the highest rank. If multiple matches have the highest rank, it will choose between them randomly.
-   BLAST hits that could not be found in names.dmp will be noted in the log file. These hits will still count towards the taxonomic diversity score, but cannot be used to generate an intersection. This may also cause reads to be unnecessarily discarded if it affects the top or second-best hit.
-   A few BLAST hits (often to human sequences) do not follow the normal naming convention. I am slowly improving the name search function to pick these up.
-   Also, more recent names.dmp files contain more names so may allow more matches, so results can change with names.dmp versions.
