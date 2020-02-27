#!/usr/bin/perl 

######################################################
####################### PIA.pl #######################
######################################################
## split and run Phylogenetic Intersection Analysis ##
############## Roselyn ware, UoW 2018 ################
######################################################
############## Version 5.1, 2020-02-27 ###############
######################################################

my $PIA_version = '5.1'; # String instead of numeric because otherwise it can 'helpfully' remove '.0'.

# Edited by Roselyn Ware, UoW 2018
# Further edited by Becky Cribdon, UoW 2019
		
# A method of metagenomic phylogenetic assignation, designed to be robust to partial representation of organisms in the database		
# Please report any problems to r.cribdon@warwick.ac.uk.

	use strict;
	use warnings;
	use Getopt::Std;
	use File::Copy qw(copy);
    use File::Path;
	#use Data::Dumper qw(Dumper); # Only used for testing

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
	Option	Description			Input		Explanation
	-f	FASTA or header filename	Y		FASTA of reads to PIA. Sequence names must be short enough for BLAST to not crop them.
	-b	BLAST filename			Y		BLAST filename containing entries for all reads in the FASTA (can contain other entries too).


Optional
	Option	Description			Input		Explanation
 	-c 	cap				Optional	Maximum unique BLAST taxa examined. Impacts taxonomic diversity score. Default is 100.
	-C	min % coverage			Optional	Minimum percentage coverage a top BLAST hit must have for a read to be taken forward. Default is 95.
	-h	help				N		Print this help text.
	-s	min tax diversity score		Optional	Minimum taxonomic diversity score for a read to make it to Summary_Basic.txt. Depends on cap. Default is 0.1.
	-t	threads				Optional	PIA.pl only. Split the header file into x subfiles and run PIA_inner.pl on each one. Default is 1.
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
    print "\nMax taxa to consider: $cap\n";

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
    my $record = <$fasta_filehandle>; # Read the next line from the names file.
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
$lines=~ s/\.\d+$//; # $lines is probably now a decimal. Remove any "." characters followed by digits. So, $lines gains 1 but loses whatever fraction it had. It's rounded up to the nearest integer.

# Split the header file accordingly:
system("split -a 3 -l $lines $header_filename $header_filename"); # Split the header file such that each new file has $lines lines. The new files will be named [header file].aaa, [header file].aab, [header file].aac...
my @header_filename_elements = split/\//, $header_filename; # Split $header_filename on "/" symbols in case it's not actually a filename, but a path.
$header_filename = pop @header_filename_elements; # Redefine $header_filename as the final element, so the actual filename if $header_filename is a path. If it isn't a path, $header_filename stays the same.

my $splitfiles =`ls . | grep $header_filename.`; # "ls ." lists files and directories in the current directory. The grep narrows this list to files and directories starting with the header filename. It ends up being a list of the split header files.
my @splitfiles= split /\n/, $splitfiles; # Save that list in @splitfiles.

# Make the command file.
my $shellscript = 'shellscript_' . $header_filename . '.txt';
unless(open FILE, '>'."$shellscript") { die "\nUnable to create $shellscript\n"; }

foreach my $file (@splitfiles){ # For each split header file, print to $shellscript the command to run the PIA_inner.pl file with the relevant options we checked earlier.
		print FILE "perl PIA_inner.pl -f $file -b $blast_filename -c $cap -C $min_coverage_perc -s $min_taxdiv_score &\n";
} # Ending in "&" means that the command is told to run in the background, so new commands can be started on top. The split header files will be processed alongside each other. That's threading.

close FILE; # Close $shellscript. We're finished editing its contents.
chmod 0755, "$shellscript"; # But still need to change the permissions.
`./$shellscript `; # Finally, run it.


# Collate the PIA_inner.pl outputs from each thread
#--------------------------------------------------
my @splitfiles2 = `ls . | grep $header_filename... | grep -v "out" `; # -v means "invert". This line searches for files or directories starting with $tophitfile but that don't contain "out". The items are separated by newlines.

my @splitfolders = `ls . | grep $header_filename... | grep "out" `; # This searches for files or directories starting with $tophitfile that do contain "out". These will be directories that PIA_inner.pl has just made, such as test.headeraaa_out/. The items are separated by newlines.

my @outIntersects;
my @Summary;
my @logs;

foreach my $file (@splitfolders){ # For every directory that PIA_inner.pl just made, look at its files:
		chomp $file; # Not sure why, but this is necessary.
	
		if (`ls $file | grep "_out.intersects.txt"`){
				push @outIntersects, ($file."/".`ls $file | grep "_out.intersects.txt"| grep -v "_out.intersects.txt_Summary_Basic.txt"`);
		}
		if (`ls $file | grep "_out.intersects.txt_Summary_Basic.txt" `){
				push @Summary, ($file."/".`ls $file | grep "_out.intersects.txt_Summary_Basic.txt"`);
		}
		if (`ls $file | grep "_PIA_inner_log.txt"`){
				push @logs, ($file."/".`ls $file | grep "_PIA_inner_log.txt"`); # List all log files in @logs.
		}
}


mkdir "$header_filename"."_out"; # Make a master output folder to collect all output in.

# Now make a variable for the filename of each of these PIA_inner.pl output files we've just looked at.
my $OI = $header_filename."_out/".$header_filename."_out.intersects.txt"; # So, $OI represents a collated, master out.intersects file inside the output folder. Note that this and $S get "out." because they are more of an output than the log.
my $S = $header_filename."_out/".$header_filename."_out.intersects.txt_Summary_Basic.txt";
my $L = $header_filename."_out/".$header_filename."_PIA_inner_logs.txt";

foreach my $data_file (@outIntersects) {
	chomp $data_file;
	`cat $data_file >>$OI`;
}
foreach my $data_file (@Summary) {
	chomp($data_file);
	system("cat $data_file >>$S");
}
foreach my $data_file (@logs) {
	chomp($data_file);
	system("cat $data_file >>$L");
}


# Collapse duplicates in the summary basic file
#----------------------------------------------
if (-e $S) { # If there were actually any summary basics to start with.
        open (my $S_filehandle, $S) or die "Could not open summary basic file $S for collapsing: $!\n"; # Read in the uncollapsed file.
        my %taxa_and_hits = ();
        foreach my $line (readline($S_filehandle)) {
            
            if (index ($line, '#') != -1) {
                next; # If the line contains a hash symbol, which indicates the header line, skip it.
            }
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
        close $S_filehandle;
        
        my @OI = split ('/', $OI); # Take the intersects file name. If $OI is a path, find the file name on the end.
        my $sample_filename = $OI[-1];
        
        open ($S_filehandle,'>', $S) or die "Could not open summary basic file $S: $!\n"; # Open the file again, this time to overwrite.
        
        # Print a new header section including the end time. Preface every new line with '#' to make ignoring them easier.
        my $timestamp_end = localtime();
        print $S_filehandle "# Start: $timestamp_start\tEnd: $timestamp_end\n# PIA version:\t$PIA_version\n# Input FASTA:\t$fasta_filename\n# Input BLAST:\t$blast_filename\n# Minimum coverage for top BLAST hit:\t$min_coverage_perc %\n# Cap of BLAST taxa to examine:\t\t$cap\n# Minimum taxonomic diversity score:\t$min_taxdiv_score\n# Number of threads:\t$threads\n#\n# ID\tName\tReads\n";
        
        foreach my $taxon (keys %taxa_and_hits) { # If there weren't any hits in any summary basic, this hash will be empty. But that shouldn't throw an error.
            print $S_filehandle "$taxon\t$taxa_and_hits{$taxon}\n";
        }
} else {
    print "\nPIA finished, but no summary basic files found.\n\n";
}


# Tidy up
#--------
foreach my $file (@splitfolders){ # For every directory that PIA_inner.pl made,
	chomp $file;
	rmtree $file or warn "Could not rmtree $file: $!"; # Delete it and any subdirectories.
}

foreach my $file (@splitfiles2){
	chomp $file;
	unlink $file or warn "Could not unlink $file: $!"; # Delete the other temporary outputs too.
}

unlink $shellscript; # Delete the list of PIA_inner.pl commands.
unlink $header_filename; # Delete the original header file.

exit;
