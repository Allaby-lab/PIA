#!/usr/bin/perl 

######################################################
####################### PIA.pl #######################
######################################################
## split and run Phylogenetic Intersection Analysis ##
############## Roselyn ware, UoW 2018 ################
######################################################
############## Version 5.5, 2021-05-03 ###############
######################################################

my $PIA_version = '5.5'; # String instead of numeric because otherwise perl can 'helpfully' remove '.0'.

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
# - Generate a command file that will run PIA_inner.pl on x threads.
# - Collate the x sets of output.
# - Collapse the Summary_Basic.txts (if there are any).


######################################################										#####################
########### Check arguments and Input Data ###########										###### Modules ######
######################################################										#####################

##### Get arguments from command line #####
	my %options=();
	getopts('hf:b:c:C:s:t:', \%options); 														#Getopt::Std
    
	# If other text found on command line, do:
	print "Other things found on the command line:\n" if $ARGV[0];
	foreach (@ARGV)	{
        print "$_\n";
	}
    
##### Display help file and exit if -h flag called #####
	my $helpfile="Helpfile_PIA.txt";
	if ($options{h}){
        print "Usage: perl PIA.pl -f <file> -b <blast.txt> [options]

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
	
##### Check FASTA filename and open file #####
	# The header file is the read names (headers) extracted from a FASTA file. PIA.pl makes the header file.
	my $fasta_filename = ($options{f});
	
##### Check BLAST filename and open file #####
 	my $blast_filename = ($options{b});

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


# Generate a header file
#-----------------------
# This will contain every header from the FASTA file plus the corresponding sequence length.
my $header_filename = $fasta_filename . '.header';
open (my $header_filehandle, '>', $header_filename) or die "Cannot open $header_filename for writing: $!\n";
open (my $fasta_filehandle, $fasta_filename) or die "Cannot open FASTA $fasta_filename: $!\n";

$/ = '>'; # Set the input record separator to '>', the first character of a header line and of a FASTA record.
while (1) { # Run this loop until "last" is called.
    my $record = readline($fasta_filehandle); # Read the next line from the names file.
    if (! defined $record) { last }; # If there is no next line, exit the loop. You've processed the whole file. 
    
    my @record = split("\n", $record); # Split the record into its two lines.
    if ($record[1]) { # Is there ever a record without a sequence? Whoever first wrote PIA thought so.
        my @sequence_name = split (' ', $record[0]); # The 0th record is the sequence name (qseqid in BLAST lingo). If there's a space BLAST only takes what comes before the space, so we must do the same.
        my $sequence_length = length $record[1]; # The sequence is the 1st line.
        print $header_filehandle $sequence_name[0] . "\t$sequence_length\n";
    }
}
close $fasta_filehandle;
$/ = "\n"; # Set it back to default.


# Generate a command file to run PIA_inner.pl on x threads
#---------------------------------------------------------
# How many headers will each thread have to process?
my $lines =`wc -l $header_filename`; # wc returns the number of lines in the header file followed by the name of the header file.
$lines =~ s/^\s+//; # Remove any leading whitespace.
my @lines = split (/ /, $lines); # Then split on whitespace. So, there will be two elements: the number of lines and the filename.
$lines = "$lines[0]"; # Redefine $lines to be just the number of lines. This is equivalent to the number of headers.

$lines = $lines/$threads; # Redefine $lines to be the number of headers to process per thread.
$lines = $lines + 1; # Add 1. Does this ensure it's never 0-point-something?
$lines =~ s/\.\d+$//; # $lines is probably now a decimal. Remove any "." characters followed by digits. So, $lines gains 1 but loses whatever fraction it had. It's rounded up to the nearest integer.

# Split the header file accordingly:
system("split -a 3 -l $lines $header_filename $header_filename"); # Split the header file such that each new file has $lines lines. The new files will be named [header file].aaa, [header file].aab, [header file].aac...
my @header_filename_elements = split/\//, $header_filename; # Split $header_filename on "/" symbols in case it's not actually a filename, but a path.
$header_filename = pop @header_filename_elements; # Redefine $header_filename as the final element, so the actual filename if $header_filename is a path. If it isn't a path, $header_filename stays the same.

my $grep_search_expression = '^' . $header_filename . '[a-z]{3}$';
my $splitfiles = `ls . | grep -E $grep_search_expression`; # "ls ." lists files and directories in the current directory. The grep -E narrows this list to files and directories starting with the header filename, then three lower-case letters. It ends up being a list of the split header files.
my @splitfiles= split /\n/, $splitfiles; # Save that list in @splitfiles.

# Make the command file.
my $shellscript_filename = 'shellscript_' . $header_filename . '.txt';
open (my $shellscript_filehandle, '>'."$shellscript_filename") or die "\nUnable to create $shellscript_filename: $!\n";

foreach my $splitfile (@splitfiles){ # For each split header file, print to $shellscript the command to run the PIA_inner.pl file with the relevant options we checked earlier.
	print $shellscript_filehandle "perl PIA_inner.pl -H $splitfile -b $blast_filename -c $cap -C $min_coverage_perc -s $min_taxdiv_score &\n";
} # Ending in "&" means that the command is told to run in the background, so new commands can be started on top. The split header files will be processed alongside each other. That's threading.

close $shellscript_filehandle; # Close the shellscript file. We're finished editing its contents.
`bash $shellscript_filename`; # Run it.


# Concatenate the PIA_inner.pl outputs from each thread
#------------------------------------------------------
$grep_search_expression = '^' . $header_filename . '[a-z]{3}_out$'; # Add '_out' to pick up only the output directories from each split header file.
my @splitfolders = `ls . | grep -E $grep_search_expression `;

my @outIntersects;
my @outBasics;
my @outReads;
my @logs;

foreach my $splitfolder (@splitfolders){ # For every directory that PIA_inner.pl just made, look at its files:
		chomp $splitfolder;

        my $intersects_filename = $splitfolder . '/' . $splitfolder . '.intersects.txt';
        push (@outIntersects, $intersects_filename);
        push (@outBasics, $intersects_filename . '_Summary_Basic.txt');
        push (@outReads, $intersects_filename . '_Summary_Reads.txt');
        
        my $log_filename = $splitfolder . '/' . $splitfolder . '_PIA_inner_log.txt';
        push (@logs, $log_filename);
}


mkdir "$header_filename"."_out"; # Make a master output folder to collect all output in.

# Now make a variable for the filename of each of the PIA_inner.pl output files.
my $intersects_filename = $header_filename."_out/".$header_filename."_out.intersects.txt";
my $SB_filename = $header_filename."_out/".$header_filename."_out.intersects.txt_Summary_Basic.txt";
my $SR_filename = $header_filename."_out/".$header_filename."_out.intersects.txt_Summary_Reads.txt";
my $SRM_filename = $header_filename."_out/".$header_filename."_out.intersects.txt_Summary_Reads_MEGAN.csv"; # To be filled in later.
my $output_fasta_filename = $header_filename."_out/".$header_filename."_out.intersects.fasta"; # To be filled in later.
my $log_filename = $header_filename."_out/".$header_filename."_PIA_inner_logs.txt";

foreach my $data_file (@outIntersects) {
	chomp $data_file;
	`cat $data_file >> $intersects_filename`;
}
foreach my $data_file (@outBasics) {
	chomp($data_file);
	system("cat $data_file >>$SB_filename");
}
foreach my $data_file (@outReads) {
	chomp($data_file);
	system("cat $data_file >>$SR_filename");
}
foreach my $data_file (@logs) {
	chomp($data_file);
	system("cat $data_file >>$log_filename");
}


# Collapse duplicates in the Summary Basic file
#----------------------------------------------
if (-e $SB_filename) { # If there were actually any Summary Basics to start with.
    open (my $SB_filehandle, $SB_filename) or die "Could not open $SB_filename for collapsing: $!\n"; # Read in the combined file.
    my %taxa_and_hits = ();
    
    while (1) { # Run this loop until last is called.
        my $line = readline($SB_filehandle);
        if (! defined $line) { last; } # If you've reached the end of the file, exit the loop.   
        if (index ($line, '#') != -1) { next; } # If the line contains a hash symbol, which indicates the header line, skip it.

        chomp $line;
        my @line = split("\t", $line);
    
        my $ID_and_name = $line[0] . "\t" . $line[1]; # The first two columns are ID and name.
        my $hit_count = $line[2]; # The final column is count.
        if (exists $taxa_and_hits{$ID_and_name}) {
            $taxa_and_hits{$ID_and_name} = $taxa_and_hits{$ID_and_name} + $hit_count;
        } else {
            $taxa_and_hits{$ID_and_name} = $hit_count;
        }
    }
    close $SB_filehandle;
    
    my $total_reads = 0;
    foreach my $ID_and_name (keys %taxa_and_hits) {
        $total_reads = $total_reads + $taxa_and_hits{$ID_and_name};
    }
    
    my @intersects_filename = split ('/', $intersects_filename); # Take the intersects file name. If $intersects_filenamew is a path, find the file name on the end.
    my $sample_filename = $intersects_filename[-1];
    
    open ($SB_filehandle, '>', $SB_filename) or die "Could not open $SB_filename: $!\n"; # Open the file again, this time to overwrite.
    
    # Print a new header section including the end time. Preface every new line with '#' to make ignoring them easier.
    my $timestamp_end = localtime();
    print $SB_filehandle "# Start: $timestamp_start\tEnd: $timestamp_end\n# PIA version:\t$PIA_version\n# Input FASTA:\t$fasta_filename\n# Input BLAST:\t$blast_filename\n# Minimum coverage for top BLAST hit:\t$min_coverage_perc %\n# Cap of BLAST taxa to examine:\t\t$cap\n# Minimum taxonomic diversity score:\t$min_taxdiv_score\n# Number of threads:\t$threads\n#\n# Total passed reads: $total_reads\n#\n# ID\tName\tReads\n";
    
    foreach my $taxon (keys %taxa_and_hits) { # If there weren't any hits in any Summary Basic, this hash will be empty. But that shouldn't throw an error.
        print $SB_filehandle "$taxon\t$taxa_and_hits{$taxon}\n";
    }
    close $SB_filehandle;
        
} else {
    print "\nPIA finished, but no Summary Basic files found.\n\n";
}


# Re-format the Summary Reads file and generate from it a CSV that MEGAN can read
#--------------------------------------------------------------------------------
my @SR_reformatted = ();
if (-e $SR_filename) { # If there were any Summary Reads files,
    open (my $SR_filehandle, $SR_filename) or die "Could not open $SR_filename for collapsing: $!\n"; # Read in the combined file.
    
    while(1) {
        my $line = readline($SR_filehandle);
        if (! defined $line) { last; }
            
        if (index ($line, '#') != -1) {
             next; # If the line contains a hash symbol, which indicates the header line, skip it.
        } else {
            push (@SR_reformatted, $line)
        }
    }
    close $SR_filehandle;
    my $total_reads = scalar (@SR_reformatted);
    
    my @intersects_filename = split ('/', $intersects_filename); # Take the intersects file name. If $intersects_filenamew is a path, find the file name on the end.
    my $sample_filename = $intersects_filename[-1];
    
    open ($SR_filehandle, '>', $SR_filename) or die "Could not open $SR_filename: $!\n"; # Open the file again, this time to overwrite.
        
    # Print the same new header section as for the new Summary Basic file.
    my $timestamp_end = localtime();
    print $SR_filehandle "# Start: $timestamp_start\tEnd: $timestamp_end\n# PIA version:\t$PIA_version\n# Input FASTA:\t$fasta_filename\n# Input BLAST:\t$blast_filename\n# Minimum coverage for top BLAST hit:\t$min_coverage_perc %\n# Cap of BLAST taxa to examine:\t\t$cap\n# Minimum taxonomic diversity score:\t$min_taxdiv_score\n# Number of threads:\t$threads\n#\n# Total passed reads: $total_reads\n#\n# Read\tID\tName\n";
    
    # Then print the rest of the data.
    foreach my $line (@SR_reformatted) {
        print $SR_filehandle $line;
    }
    close $SR_filehandle;
    
    # Print in a different format to the MEGAN CSV.
    open (my $SRM_filehandle, '>', $SRM_filename) or die "Could not open $SRM_filename: $!\n"; # Open the MEGAN CSV.
    foreach my $line (@SR_reformatted) {
        my @line = split("\t", $line); # Split the line on tabs.
        print $SRM_filehandle "$line[0],$line[1],50\n" # Export the read name, ID, and a stand-in bitscore of 50 (the default minimum pass score for the LCA).
    }
    close $SRM_filehandle;    
}


# Use the Summary Reads data to produce a post-PIA FASTA
#-------------------------------------------------------
# Convert the Summary Reads data into a hash for easier lookup.
my %SR_data = (); # Keys are read names, values are taxonomic IDs.
foreach my $SR_line (@SR_reformatted) {
    my @SR_line = split("\t", $SR_line);
    $SR_data{$SR_line[0]} = $SR_line[1];
}

# Filter the input FASTA.
open($fasta_filehandle, $fasta_filename) or die "Could not open original FASTA file $fasta_filename: $!\n";
open(my $output_fasta_filehandle, '>', $output_fasta_filename) or die "Could not open post-PIA FASTA file $output_fasta_filename: $!\n";

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

$/ = "\n"; # Set the record separator back to default.


# Tidy up
#--------
foreach my $directory (@splitfolders){ # Delete the split output directories.
	chomp $directory;
	rmtree $directory or warn "Could not rmtree $directory: $!";
}

foreach my $file (@splitfiles){
	chomp $file;
	unlink $file or warn "Could not unlink $file: $!"; # Delete the split header files.
}

if (-z $output_fasta_filename) { unlink $output_fasta_filename }; # If the input FASTA was empty, the only output file made will be an empty output FASTA. Delete it if it is empty.

unlink $shellscript_filename; # Delete the list of PIA_inner.pl commands.
unlink $header_filename; # Delete the original header file.

exit;
