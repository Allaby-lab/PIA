PROVISIONAL - THIS IS A BETA VERSION

The Phylogenetic Intersection Analysis
======================================
Metagenomic phylogenetic assignment of mixed environmental assemblages
Allaby lab, University of Warwick
Version 4.0
2019-07-20

The phylogenetic intersection analysis (PIA) takes standard-format BLAST output and a corresponding FASTA file. It assigns reads to phylogenetic intersections based on their BLAST hits, assuming that the true taxon will be inside that phylogenetic intersection. It is designed to be robust to the uneven representation of taxa in databases.

A very early version was published in Smith et al. ("Sedimentary DNA from a submerged site reveals wheat in the British Isles 8000 years ago.", Science, 2015) and a more current version will be in a forthcoming paper this year. For more information, email r.cribdon@warwick.ac.uk.


Prerequisites
-------------
-   Perl 5
-   Perl module List::MoreUtils


Files and directories
---------------------
Modules/
    FileChecks.pm
    FileManipulations.pm
    FileMerge.pm
    TreeOfLife.pm
Helpfile_PIA.txt
PIA.pl
PIA_inner.pl

You need to make a Reference_files/ directory containing names.dmp and nodes.dmp, which can be downloaded from the NCBI here: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/
Reference_files/
    names.dmp
    nodes.dmp


Usage
------
-   Have all PIA scripts and directories in the same directory.
-   PIA_inner.pl does the analysis. PIA.pl is a wrapper that allows threading. To run PIA.pl:
    >perl PIA.pl -f [input FASTA] -b [input BLAST file] -t [x number of threads] [other options]
-   To run PIA_inner.pl alone, see notes at the top of PIA_inner.pl.
-   Outputs will be in [input FASTA]_out/.
    -   The intersects file is the full output.
    -   The summary basic file is a simple list of how many reads were assigned to which phylogenetic intersections.
    -   The log tracks which reads (headers) have been processed. Use this to pick up where you left off if the PIA is interrupted. It also states when a BLAST hit could not be found in names.dmp.
    -   The timer states how long PIA.pl took to run.

Please report any problems to r.cribdon@warwick.ac.uk.


Known issues
------------
-   The extended summary function is untested.
-   BLAST hits are assigned to taxa using names.dmp. The PIA currently cannot distinguish between synonymous taxa, as can happen across kingdoms. If a BLAST hit is to Iris, it will be    randomly assigned to either Iris the angiosperm or Iris the mantid. This may cause reads to be unnecessarily discarded if it affects the top or second-best hit.
-   Results will change depending on the names.dmp file used. BLAST hits that could not be found in names.dmp will be noted in the log file. These hits will still count towards the taxonomic diversity score, but cannot be used to generate an intersection. This may also cause reads to be unnecessarily discarded if it affects the top or second-best hit.
