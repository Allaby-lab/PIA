#!/usr/bin/perl 

######################################################
####################### PIA.pl #######################
######################################################
## split and run Phylogenetic Intersection Analysis ##
############## Roselyn ware, UoW 2018 ################
######################################################
############## Version 4.4, 2019-08-13 ###############
######################################################

# Edited by Roselyn Ware, UoW 2018

# Further edited by Becky Cribdon, UoW 2019
# - Deals with log files from PIA_inner.pl and makes its own timer log.
# - Collapses hits in the final Summary_Basic.txt.
# - Does not expect re-PIAing or re-BLASTing or whatnot.
# - Converts the non-BLAST input file into a header file. This allows it to accept FASTAs as well as ready-made header files.
# - Should be able to run simultanously on multiple FASTAs in the same directory.

	use strict;
	use warnings;
	use lib './Modules'; # Add the /Modules directory to @INC: the list of places Perl will look for modules. 
	use FileChecks; # A bespoke module in /Modules.
	use TreeOfLife; # A bespoke module in /Modules.
	use Getopt::Std;
	use File::Find;
	use Data::Dumper qw(Dumper);
	use FileMerge; # A bespoke module in /Modules.
	use FileManipulations; # A bespoke module in /Modules.
	use File::Copy qw(copy);
	use File::Path;

# Run as follows:
# perl PIA.pl -f [FASTA or corresponding header file] -b [BLAST file] -p [taxonomic ID of expected phylogenetic range] -t [number of threads]

# Summary:
# - Check arguments and inputs.
# - Generate a command file that will run PIA_inner.pl on x threads.
# - Collate the x sets of output.
# - Collapse the Summary_Basic.txts.


######################################################										#####################
########### Check arguments and Input Data ###########										###### Modules ######
######################################################										#####################

##### Get arguments from command line #####
	my %options=();
	getopts('f:b:c:C:ehn:N:p:t:', \%options); 														#Getopt::Std
    
	# If other text found on command line, do:
	print "Other things found on the command line:\n" if $ARGV[0];
	foreach (@ARGV)	{
        print "$_\n";
	}			
	
##### Check FASTA filename and open file #####
	# The header file is simply the read names (headers) extracted from a FASTA file. PIA.pl makes the header file.
	# Ensure that a filename has been given as an argument. -f flag
	my $fasta_filename = FileChecks::process_filename($options{f});									#FileChecks.pm
	chomp $fasta_filename ;
	
##### Check BLAST filename and open file #####
	# Ensure that a filename has been given as an argument. -b flag
 	my $blastfile = FileChecks::process_filename($options{b});									#FileChecks.pm
	chomp $blastfile ;

##### See if cap input #####
	# The cap impacts the taxon diversity score. Default is 100.	
	# Ensure that a number has been given as an argument. -c flag
	my $capopt = ($options{c});
	my $cap = 100;
	if ($capopt) { $cap = $capopt;} # If there is an option, overwrite the default.

##### See if min % coverage input #####
	# The minimum percentage coverage a top BLAST hit must have for a read to be taken forward. Default is 95.
	# Ensure that a number has been given as an argument. -C flag
	my $min_coverage_perc_opt = ($options{C});
	my $min_coverage_perc = 95;
	if ($min_coverage_perc_opt) { $min_coverage_perc = $min_coverage_perc_opt;} # If there is an option, overwrite the default.

##### See if extended summary required #####
	# Create extended summary if -e flag called #		
	my $summary;
	if ($options{e}){ $summary = 1;}	

##### Display help file and exit if -h flag called #####
	my $helpfile="Helpfile_PIA.txt";
	if ($options{h}){
		FileChecks::process_help($helpfile);														#FileChecks.pm	
	}

##### Locate nodes.dmp and names.dmp files #####
	# If nodes.dmp and names.dmp aren't in their default location (Reference_files/), use -n and -N to specify where they are.	
	my $nodesfile="Reference_files/nodes.dmp";
	if ($options{n}){
		$nodesfile=($options{n});
	}
	my $namesfile="Reference_files/names.dmp";   
	if ($options{N}){
		$namesfile=($options{N});
	}

##### See if expected phylogenetic range input #####	
	# Ensure that a number has been given as an argument. -p flag
	my $expected_phylogenetic_range_opt = ($options{p});
	my $expected_phylogenetic_range = 1;
	if ($expected_phylogenetic_range_opt) {$expected_phylogenetic_range = $expected_phylogenetic_range_opt;}

	##### Check see if thread number input #####
	# Ensure that a number has been given as an argument. -t flag
	my $threadopt = ($options{t});																	#FileChecks.pm
	my $thread = 2; # Set to 2 by default.
	if ($threadopt) { $thread = $threadopt;}	


# Start timer file
#-----------------
my $timer_filename = 'timer_' . $fasta_filename . '.txt';
open( my $timer_filehandle, '>', $timer_filename) or die "Cannot open $timer_filename for writing.\n$!\n"; # Note: this will overwrite old timer. I didn't see any point in appending old ones when you need to remove the other PIA outputs before running it again anyway.
use IO::Handle; # Enable autoflush.
$timer_filehandle -> autoflush(1); # Set autoflush to 1 for the log filehandle. This means that Perl won't buffer its output to the log, so the log should be updated in real time.
my $timestamp = localtime();
print $timer_filehandle "Input FASTA: $fasta_filename\nInput BLAST: $blastfile\n\nPIA started at $timestamp\n\n";


# Generate a header file
#-----------------------
my $header_filename = $fasta_filename . '.header';
open (my $header_filehandle, '>', $header_filename) or die "Cannot open $header_filename for writing.\n$!\n";
open (my $fasta_filehandle, $fasta_filename) or die "Cannot open $fasta_filename.\n$!\n";

while (1) { # Run this loop until "last" is called.
            my $line = <$fasta_filehandle>; # Read the next line from the names file.
            if (! defined $line) { last }; # If there is no next line, exit the loop. You've processed the whole file.
            
            if (index ($line, '>') != -1) { # If this line starts with a ">", it's a header line, so print it to the header file.
                print $header_filehandle $line;
            }
}
close $fasta_filehandle;


# Generate a command file to run PIA_inner.pl on x threads
#---------------------------------------------------------
# How many headers will each thread have to process?
my $lines=`wc -l $header_filename`; # wc returns the number of lines in the header file followed by the name of the header file.
$lines=~ s/^\s+//; # Remove any leading whitespace.
my @lines= split / /,$lines; # Then split on whitespace. So, there will be two elements: the number of lines and the filename.
$lines="$lines[0]"; # Redefine $lines to be just the number of lines. This is equivalent to the number of headers.

$lines=$lines/$thread; # Redefine $lines to be the number of headers to process per thread.
$lines=$lines+1; # Add 1. Does this ensure it's never 0. something?
chomp $lines; # Not sure what the chomp is for. If $lines was a number and a whitespace character, the whitespace would have been removed when we treated it as a numeric and added 1 to it. If it couldn't be treated as a numeric, there would have been an error.
$lines=~ s/\.\d+$//; # $lines is probably now a decimal. Remove any "." characters followed by digits. So, $lines gains 1 but loses whatever fraction it had. It's rounded up to the nearest integer.
#print "Number of headers each thread will process: $lines\n";

# Split the header file accordingly:
system("split -a 3 -l $lines $header_filename $header_filename"); # Split the header file such that each new file has $lines lines. The new files will be named [header file].aaa, [header file].aab, [header file].aac...
my @header_filename_elements = split/\//, $header_filename; # Split $header_filename on "/" symbols in case it's not actually a filename, but a path.
$header_filename = pop @header_filename_elements; # Redefine $header_filename as the final element, so the actual filename if $header_filename is a path. If it isn't a path, $header_filename stays the same.

my $splitfiles =`ls . | grep $header_filename.`; # "ls ." lists the files and directories in the current directory. The grep narrows this list to files and directories starting with the header filename. It ends up being a list of the split header files.
my @splitfiles= split /\n/, $splitfiles; # Save that list in @splitfiles.

# Make the command file.
my $shellscript = 'shellscript_' . $header_filename . '.txt';
unless(open FILE, '>'."$shellscript") { die "\nUnable to create $shellscript\n"; } # If the $shellscript file can't be created, die with an error.

foreach my $file (@splitfiles){ # For each split header file, print to $shellscript the command to run the PIA_inner.pl file with the relevant options we checked earlier.
		print FILE "perl PIA_inner.pl -f $file -b $blastfile -c $cap -C $min_coverage_perc -n $nodesfile -N $namesfile -p $expected_phylogenetic_range &\n";
} # So, if we run these commands simultaneously, we'll be processing each of the split header files simultaneously, and that's threading.
# The fact that the command ends in "&" means that the command is told to run in the background, so new commands can be started on top. The commands should run alongside each other.

close FILE; # Close $shellscript. We're finished editing its contents.

chmod 0755, "$shellscript"; # But still need to change the permissions.

`./$shellscript `; # Finally, run it.


# Collate the PIA_inner.pl outputs from each thread
#---------------------------------------------
my @splitfiles2= `ls . | grep $header_filename... | grep -v "out" `; # -v means "invert". This line searches for files or directories starting with $tophitfile but that don't contain "out".
my @splitfolders= `ls . | grep $header_filename... | grep "out" `; # This searches for files or directories starting with $tophitfile that do contain "out". These will be directories that PIA_inner.pl has just made, such as test.headeraaa_out/.

my @outIntersects;
my @Summary;
my @logs;

foreach my $file (@splitfolders){ # For every directory that PIA_inner.pl just made:
		chomp $file;
	
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
	chomp($data_file);
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
open (my $S_filehandle, $S) or die "Could not open summary basic file $S for collapsing.\n$!\n"; # Read in the uncollapsed file.
my %taxa_and_hits = ();
foreach my $line (readline($S_filehandle)) {
    
    if (index ($line, '#') != -1) {
        next; # If the line contains a hash symbol, which indicates the header line, skip it.
    }
    chomp $line;
	my @line = split("\t", $line);

    my $taxon_and_ID = $line[0] . "\t" . $line[1]; # The first two columns are name and ID.
	my $hit_count = $line[2]; # The final column is count.
	if (exists $taxa_and_hits{$taxon_and_ID}) {
		$taxa_and_hits{$taxon_and_ID} = $taxa_and_hits{$taxon_and_ID} + $hit_count;
	} else {
		$taxa_and_hits{$taxon_and_ID} = $hit_count;
    }
}
close $S_filehandle;

my @OI = split ('/', $OI); # Take the path to the intersects file.
my $summary_basic_header = $OI[-1]; # If it contains directories as well as the file name, take just the file name.

open ($S_filehandle,'>', $S) or die "Could not open summary basic file $S.\n$!\n"; # Open the file again, this time to overwrite.
print $S_filehandle "#Series:\t$summary_basic_header\n"; # Add a new header.
foreach my $taxon (keys %taxa_and_hits) {
	print $S_filehandle "$taxon\t$taxa_and_hits{$taxon}\n";
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

$timestamp = localtime();
print $timer_filehandle "PIA finished at $timestamp\n";

# Move the timer into the final output folder.
my $timer_final_destination = $header_filename ."_out";
system ("mv $timer_filename $timer_final_destination");

exit;
