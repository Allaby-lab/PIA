#!/usr/bin/perl 

######################################################
####################### PIA.pl #######################
######################################################
## split and run Phylogenetic Intersection Analysis ##
############## Roselyn ware, UoW 2018 ################
######################################################
############## Version 5.6, 2021-05-06 ###############
######################################################

my $PIA_version = '5.6'; # String instead of numeric because otherwise perl can 'helpfully' remove '.0'.

# Edited by Roselyn Ware, UoW 2018
# Further edited by Becky Cribdon, UoW 2021
		
# A method of metagenomic phylogenetic assignation, designed to be robust to partial representation of organisms in the database		
# Please report any problems to r.cribdon@warwick.ac.uk.

	use strict;
	use warnings;
	use Getopt::Std;
	use File::Copy qw(copy);
    use File::Path;
	use Data::Dumper qw(Dumper); # Only used for testing

# Run as follows:
# perl PIA.pl -f [FASTA file] -b [BLAST file] -t [number of threads]

# Summary:
# - Check arguments and inputs.
# - Generate a command file that will run PIA_inner.pl on n threads.
# - Collate the n sets of output.


######################################################
########### Check arguments and Input Data ###########	
######################################################

##### Get arguments from command line #####
	my %options=();
	getopts('hf:b:c:C:s:t:', \%options);

    
##### Display help text and exit if -h flag called #####
	if ($options{h}){
        print "Usage: perl PIA.pl -H [header file] -b [BLAST file] [other options]

Main Arguments
	Option	Description	Input		Explanation
	-f	FASTA filename	Y		FASTA of reads to analyse.
	-b	BLAST filename	Y		BLAST filename containing entries for all reads in the FASTA (can contain other entries too).


Optional
	Option	Description		Input		Explanation
 	-c 	cap			Optional	Maximum unique BLAST taxa examined. Impacts taxonomic diversity score. Default is 100.
	-C	min % coverage		Optional	Minimum percentage coverage a top BLAST hit must have for a read to be taken forward. Default is 95.
	-h	help			N		Print this help text.
	-s	min diversity score	Optional	Minimum taxonomic diversity score for a read to make it to Summary_Basic.txt and Summary_Reads.txt. Depends on cap. Default is 0.1.
	-t	threads			Optional	PIA.pl only. Split the header file into x subfiles and run PIA_inner.pl on each one. Default is 1.
";
        exit;
	}
	
##### Check FASTA file path #####
	# The header file is the read names (headers) extracted from a FASTA file. PIA.pl makes the header file.
	my $fasta_filepath = '';
    if ($options{f}) {
        $fasta_filepath = $options{f};
    } else {
        print "\nSpecify a FASTA file with -f.\n\n";
        exit;
    }
	
##### Check BLAST file path #####
 	my $blast_filepath = '';
    if ($options{b}) {
        $blast_filepath = $options{b};
    } else {
        print "\nSpecify a BLAST file with -b.\n\n";
        exit;
    }

##### See if cap input #####
	# BLAST hits are all assigned to taxa. $cap is the maximum number of taxa to be looked at. If $cap is hit, no more BLAST hits will be considered. $cap is used to calculate the taxonomic diversity score, so affects the minimum score threshold (option s). Default is 100.	
	my $cap = 100;
	if ($options{c}) { # If there is an option, overwrite the default.
        $cap = $options{c};
    }
    print "\nMax taxa to consider per read: $cap\n";

##### See if min % coverage input #####
	# The minimum percentage coverage a top BLAST hit must have for a read to be taken forward. Default is 95.
	my $min_coverage_perc = 95;
	if ($options{C}) { # If there is an option, overwrite the default.
        $min_coverage_perc = $options{C};
    }
    print "Min coverage %: $min_coverage_perc\n";

##### See if min taxonomic diversity score input #####
	# The minimum taxonomic diversity score a read must have to make it to Summary_Basic.txt. Depends on $cap. Defaults to 0.1.
	my $min_taxdiv_score = 0.1;
	if ($options{s}) { # If there is an option, overwrite the default.
        $min_taxdiv_score = $options{s};
    }
    print "Min taxonomic diversity score: $min_taxdiv_score\n\n";

##### Check see if thread number input #####
	my $threads = 1; # Set to 1 by default.
	if ($options{t}) { $threads = $options{t};}	


# Start a timer
#--------------
my $timestamp_start = localtime();


# Generate a read info file
#--------------------------
# For every read in the FASTA, this will contain its name and sequence length.
my $read_info_filepath = $fasta_filepath . '.read_info';
open (my $read_info_filehandle, '>', $read_info_filepath) or die "Cannot open $read_info_filepath for writing: $!\n";
open (my $fasta_filehandle, $fasta_filepath) or die "Cannot open FASTA $fasta_filepath: $!\n";

$/ = '>'; # Set the input record separator to '>', the first character of a FASTA record.
while (1) { # Run this loop until "last" is called.
    my $record = readline($fasta_filehandle); # Read the next line from the names file.
    if (! defined $record) { last }; # If there is no next line, exit the loop. You've processed the whole file. 
    
    my @record = split("\n", $record); # Split the record into its two lines.
    if ($record[1]) { # Is there ever a record without a sequence? Whoever first wrote PIA thought so.
        my @sequence_name = split (' ', $record[0]); # The 0th record is the sequence name (qseqid in BLAST lingo). If there's a space BLAST only takes what comes before the space, so we must do the same.
        my $sequence_length = length $record[1]; # The sequence is the 1st line.
        print $read_info_filehandle $sequence_name[0] . "\t$sequence_length\n";
    }
}
close $read_info_filehandle;
close $fasta_filehandle;
$/ = "\n"; # Set it back to default.


# Generate a command file to run PIA_inner.pl on x threads
#---------------------------------------------------------
# How many reads will each thread have to process?
my $line_count =`wc -l $read_info_filepath`; # wc returns the number of lines in the read_info file followed by the name of the read_info file.
$line_count =~ s/^\s+//; # Remove any leading whitespace.
my @line_count = split (' ', $line_count); # Then split on whitespace. So, there will be two elements: the number of lines and the filename.
$line_count = "$line_count[0]"; # Redefine $line_count to be just the number of lines. This is equivalent to the number of reads.

my $line_count_per_split_read_info = $line_count/$threads; # Redefine $line_count to be the number of reads to process per thread.
$line_count_per_split_read_info = $line_count_per_split_read_info + 1; # Add 1. Does this ensure it's never 0-point-something?
$line_count_per_split_read_info =~ s/\.\d+$//; # It's probably now a decimal. Remove any "." characters followed by digits. So, $line_count_per_split_read_info gains 1 but loses whatever fraction it had. It's rounded up to the nearest integer.

# Split the read_info file accordingly:
system("split -a 3 -l $line_count_per_split_read_info $read_info_filepath $read_info_filepath"); # Split the read_info file such that each new file has $line_count lines. The new files will be named [read_info file].aaa, [read_info file].aab, [read_info file].aac...

# Collect the new read_info file names.
my $read_info_filename = $read_info_filepath; # Default the name to the path; assume there isn't a path.
my $read_info_parentpath = ''; # Default the path to empty.
if (index($read_info_filepath, '/') != -1) { # If there was in fact a path:
    my @read_info_filepath_elements = split('/', $read_info_filepath); # Split $read_info_filepath on "/" symbols to separate elements of the path.
    $read_info_filename = pop @read_info_filepath_elements; # Separate the read_info file name.
    $read_info_parentpath = join('/', @read_info_filepath_elements); # Put the rest of the path back together.
}

my $grep_search_expression = '^' . $read_info_filename . '[a-z]{3}$';
my $list_of_split_read_infos = `ls $read_info_parentpath . | grep -E $grep_search_expression`; # "ls ." lists files and directories in the current directory. The grep -E narrows this list to files and directories starting with the read_info filename, then three lower-case letters. It ends up being a list of the split read_info files.
my @split_read_infos = split("\n", $list_of_split_read_infos); # Save that list.

# Make the command file.
my $shellscript_filename = 'shellscript_' . $read_info_filename . '.txt';
open (my $shellscript_filehandle, '>'."$shellscript_filename") or die "\nUnable to create $shellscript_filename: $!\n";

foreach my $split_read_info_filename (@split_read_infos) { # For each split file, print to the shell script the command to run the PIA_inner.pl script with the relevant options we checked earlier.
    my $split_read_info_filepath = $split_read_info_filename; # Default to no parent path.
    if ($read_info_parentpath ne '') { # If there was a parent path, add it on now.
        $split_read_info_filepath = $read_info_parentpath . '/' . $split_read_info_filename;
    }
	print $shellscript_filehandle "perl PIA_inner.pl -r $split_read_info_filepath -b $blast_filepath -c $cap -C $min_coverage_perc -s $min_taxdiv_score &\n";
} # Ending in "&" means that the command is told to run in the background, so new commands can be started on top. The split files will be processed alongside each other. That's threading.

close $shellscript_filehandle; # Close the shell script file. We're finished editing its contents.
`bash $shellscript_filename`; # Run it.


# Set up to collate the PIA_inner.pl outputs from each thread
#------------------------------------------------------------
$grep_search_expression = '^' . $read_info_filename . '[a-z]{3}.PIA_output$'; # Add '_out' to pick up only the output directories from each split read_info file.
my @split_directories = `ls $read_info_parentpath . | grep -E $grep_search_expression `;

my @outFulls;
my @outSBs;
my @outSRs;
my @logs;

foreach my $split_directory (@split_directories){ # For every directory that PIA_inner.pl just made, look at its files:
	chomp $split_directory;
    
    my $base_output_filepath = $split_directory . '/' . $split_directory; # Default to no parent path.
    if ($read_info_parentpath ne '') { # If there was a parent path, add it on now.
        $base_output_filepath = $read_info_parentpath . '/' . $base_output_filepath;
    }
    $base_output_filepath = substr($base_output_filepath, 0, -11); # Remove the '.PIA_output' from the very end. The files inside the output directory don't have this.
    push (@outFulls, $base_output_filepath . '.Full.txt');
    push (@outSBs, $base_output_filepath . '.Summary_Basic.txt');
    push (@outSRs, $base_output_filepath . '.Summary_Reads.txt');
    push (@logs, $base_output_filepath . '.PIA_inner_log.txt');
}


# Make a master output folder:
my @read_info_filename_elements = split('\.', $read_info_filename);
my $sample_name = join('.', @read_info_filename_elements[0..(scalar(@read_info_filename_elements)-3)]); # The sample name is $read_info_filename minus 'read_info' and whatever FASTA extension there was: "test_sample.fa.read_infoaaa" becomes "test_sample".
my $collated_output_directory_path = $sample_name . '.PIA_output'; # Default to no parent path.
if ($read_info_parentpath ne '') { # If there was a parent path, add it on now.
    $collated_output_directory_path = $read_info_parentpath . '/' . $collated_output_directory_path;
}
mkdir $collated_output_directory_path;


# Make a variable for the filename of each of the PIA_inner.pl output files.
my $output_full_filepath = $collated_output_directory_path . '/' . $sample_name . '.Full.txt';
my $output_SB_filepath = $collated_output_directory_path . '/' . $sample_name . '.Summary_Basic.txt';
my $output_SR_filepath = $collated_output_directory_path . '/' . $sample_name . '.Summary_Reads.txt';
my $output_SRM_filepath = $collated_output_directory_path . '/' . $sample_name . '.Summary_Reads_MEGAN.txt'; # To be filled in later.
my $output_fasta_filepath = $collated_output_directory_path . '/' . $sample_name . '.Post-PIA.fasta'; # To be filled in later.
my $log_filepath = $collated_output_directory_path . '/' . $sample_name . '.PIA_inner_logs.txt';


# Read and merge the Full output files:
open(my $output_full_filehandle, '>', $output_full_filepath) or die "Could not open final Full output file $output_full_filepath for writing: $!\n";
# Print the header again:
print $output_full_filehandle "Read_name\tCoverage\tTop_hit_name\tTop_hit_ID\tExpect_value\tIdentities_value\tNext_hit_name\tNext_hit_ID\tLast_hit_name\tLast_hit_ID\tPhylogenetic_range_name\tPhylogenetic_range_ID\tRaw_hit_count\tTaxonomic_diversity_(up_to_cap_if_met)\tTaxonomic_diversity_score\tPhylogenetic_intersection_name\tPhylogenetic_intersection_ID\n";

foreach my $split_full_filepath (@outFulls) {
    open(my $split_full_filehandle, $split_full_filepath) or die "Could not open split Full output file $split_full_filepath: $!\n";
    readline($split_full_filehandle); # Skip the header.
    while(1) {
        my $line = readline($split_full_filehandle);
        if (! defined $line) { last }; # If there is no next line, exit the loop. You've processed the whole file.
        print $output_full_filehandle $line; # Print other lines to the output file.
    }
    close $split_full_filehandle;
}
close $output_full_filehandle;


# Simply concatenate the logs, but overwrite any old file first by writing a simple header:
my $timestamp_end = localtime();
`printf "# Start: $timestamp_start\tEnd: $timestamp_end\n# PIA version:\t$PIA_version\n# Input FASTA:\t$fasta_filepath\n# Input BLAST:\t$blast_filepath\n" > $log_filepath`;
foreach my $split_log_filepath (@logs) {
	chomp($split_log_filepath);
	system("cat $split_log_filepath >> $log_filepath");
}


# Collate the Summary Basics
#---------------------------
my %taxa_and_hits = ();
foreach my $split_SB_filepath (@outSBs) {
    open(my $split_SB_filehandle, $split_SB_filepath) or die "Could not open split Summary Basic $split_SB_filepath: $!\n";
    readline($split_SB_filehandle); # Skip the header.
    
    while(1) {
        my $line = readline($split_SB_filehandle);
        if (! defined $line) { last; } # If you've reached the end of the file, exit the loop.   

        chomp $line;
        my @line = split("\t", $line);
        
        if ($line[0] == 1) { next; } # Skip any lines for taxon ID 1 (root).
        my $name_and_ID = $line[1] . "\t" . $line[0]; # The first two columns are ID and name.
        my $hit_count = $line[2]; # The final column is count.
        if (exists $taxa_and_hits{$name_and_ID}) {
            $taxa_and_hits{$name_and_ID} = $taxa_and_hits{$name_and_ID} + $hit_count;
        } else {
            $taxa_and_hits{$name_and_ID} = $hit_count;
        }
    }
    close $split_SB_filehandle;
}

my $total_reads = 0;
foreach my $name_and_ID (keys %taxa_and_hits) {
    $total_reads = $total_reads + $taxa_and_hits{$name_and_ID};
}

open(my $output_SB_filehandle, '>', $output_SB_filepath) or die "Could not open final Summary Basic $output_SB_filepath for writing: $!\n";

# Print a new header section including the end time. Preface every new line with '#' to make ignoring them easier.
print $output_SB_filehandle "# Start: $timestamp_start\tEnd: $timestamp_end\n# PIA version:\t$PIA_version\n# Input FASTA:\t$fasta_filepath\n# Input BLAST:\t$blast_filepath\n# Minimum coverage for top BLAST hit:\t$min_coverage_perc %\n# Cap of BLAST taxa to examine:\t\t$cap\n# Minimum taxonomic diversity score:\t$min_taxdiv_score\n# Number of threads:\t$threads\n#\n# Total passed reads: $total_reads\n#\n# ID\tName\tReads\n";

foreach my $taxon_information (sort {lc $a cmp lc $b} keys %taxa_and_hits) { # Print the hits in alphabetical order of taxon name (not just ASCII order).
    my @taxon_information = split("\t", $taxon_information);
    print $output_SB_filehandle "$taxon_information[1]\t$taxon_information[0]\t$taxa_and_hits{$taxon_information}\n"; # If there weren't any hits in any Summary Basic, this hash will be empty. But that shouldn't throw an error.
}
close $output_SB_filehandle;


# Read and merge the Summary Reads files, and generate from them a CSV that MEGAN can read
#-----------------------------------------------------------------------------------------
open(my $output_SR_filehandle, '>', $output_SR_filepath) or die "Could not open final Summary Reads $output_SR_filepath for writing: $!\n"; # Open the collated Summary Reads file.
open (my $output_SRM_filehandle, '>', $output_SRM_filepath) or die "Could not open $output_SRM_filepath: $!\n"; # Open the MEGAN CSV.
my %SR_data = (); # Also store the data in a hash with which to filter the input FASTA. Keys are read names, values are taxonomic IDs.

# Print the same new header section as for the new Summary Basic file.
print $output_SR_filehandle "# Start: $timestamp_start\tEnd: $timestamp_end\n# PIA version:\t$PIA_version\n# Input FASTA:\t$fasta_filepath\n# Input BLAST:\t$blast_filepath\n# Minimum coverage for top BLAST hit:\t$min_coverage_perc %\n# Cap of BLAST taxa to examine:\t\t$cap\n# Minimum taxonomic diversity score:\t$min_taxdiv_score\n# Number of threads:\t$threads\n#\n# Total passed reads: $total_reads\n#\n# Read\tID\tName\n";

foreach my $split_SR_filepath (@outSRs) {
    open(my $split_SR_filehandle, $split_SR_filepath) or die "Could not open split Summary Reads $split_SR_filepath: $!\n";
    readline($split_SR_filehandle); # Skip the header.
    while(1) {
        my $line = readline($split_SR_filehandle);
        if (! defined $line) { last }; # If there is no next line, exit the loop. You've processed the whole file.
        
        my @line = split("\t", $line); # Split the line on tabs.
        if ($line[1] == 1) { next; } # Skip any lines for taxon ID 1 (root).
        print $output_SR_filehandle $line; # Print straight to the collated Summary Reads file.
        print $output_SRM_filehandle "$line[0],$line[1],50\n"; # Print to the MEGAN CSV the read name, taxonomic ID, and a stand-in bitscore of 50 (the default minimum pass score for the LCA).
        $SR_data{$line[0]} = $line[1]; # Store the read name and taxonomic ID to use with the FASTA.
    }
}

close $output_SR_filehandle;
close $output_SRM_filehandle;


# Use the Summary Reads data to produce a post-PIA FASTA
#-------------------------------------------------------
# Filter the input FASTA.
open($fasta_filehandle, $fasta_filepath) or die "Could not open original FASTA file $fasta_filepath: $!\n";
open(my $output_fasta_filehandle, '>', $output_fasta_filepath) or die "Could not open post-PIA FASTA file $output_fasta_filepath: $!\n";

$/ = '>'; # Set the record separator to '>', which separates FASTA records.

while (1) { # Run this loop until "last" is called.
    my $record = readline($fasta_filehandle);
    if (! defined $record) { last }; # If there is no next record, exit the loop. You've processed the whole file.
    
    my @record = split ("\n", $record); # Split the record into lines.
    my $read_name = $record[0];
    if (exists $SR_data{$read_name} ) {
        $read_name = $read_name . '_ID' . $SR_data{$read_name}; # Add the taxonomic ID on the end of the read name.
        $record = join("\n", $read_name, $record[1]); # Slot the new name back into the record.
        $record =~ s/>//g; # Remove the '>' character off the end so we can add a new one on the front.
        print $output_fasta_filehandle ">$record\n";
    }
}

close $fasta_filehandle;
close $output_fasta_filehandle;
$/ = "\n"; # Set the record separator back to default.


# Tidy up
#--------
foreach my $split_directory_filepath (@split_directories) { # Delete the split output directories.
	chomp $split_directory_filepath;
    if ($read_info_parentpath ne '') { # If there was a parent path, add it on now.
        $split_directory_filepath = $read_info_parentpath . '/' . $split_directory_filepath;
    }
	rmtree $split_directory_filepath or warn "Could not rmtree $split_directory_filepath: $!";
}

foreach my $split_read_info_filepath (@split_read_infos) { # Delete the split read_info files.
	chomp $split_read_info_filepath;
    if ($read_info_parentpath ne '') { # If there was a parent path, add it on now.
        $split_read_info_filepath = $read_info_parentpath . '/' . $split_read_info_filepath;
    }
	unlink $split_read_info_filepath or warn "Could not unlink $split_read_info_filepath: $!";
}

unlink $shellscript_filename; # Delete the list of PIA_inner.pl commands.
unlink $read_info_filepath; # Delete the original read_info file.

exit;