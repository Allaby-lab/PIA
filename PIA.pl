#!/usr/bin/perl 

######################################################
####################### PIA.pl #######################
######################################################
## split and run Phylogenetic Intersection Analysis ##
############## Roselyn ware, UoW 2018 ################
######################################################
############## Version 4.0, 2019-07-20 ###############
######################################################

# Edited by Roselyn Ware, UoW 2018

# Further edited by Becky Cribdon, UoW 2019
# - Deals with log files from PIA_inner.pl and makes its own timer log.
# - Collapses hits in the final Summary_Basic.txt.
# - Does not expect re-PIAing or re-BLASTing or whatnot.
		
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
# perl PIA.pl -f [FASTA or corresponding header file] -b [BLAST file] -t [number of threads]

# Summary:
# - Check arguments and inputs.
# - Generate a command file that will run PIA_inner.pl on x threads.
# - Collate the x sets of output.
# - Collate Summary_Basic.txt.


######################################################										#####################
########### Check arguments and Input Data ###########										###### Modules ######
######################################################										#####################

##### Get arguments from command line #####
	my %options=();
	getopts('hf:b:c:egn:N:G:t:', \%options); 														#Getopt::Std

	# If other text found on command line, do:
	print "Other things found on the command line:\n" if $ARGV[0];
	foreach (@ARGV)	{
							print "$_\n";
	}


##### Check other commandline arguments and store #####
	#Display help file and exit if -h flag called #		
	my $helpfile="Helpfile_PIA.txt";
	if ($options{h}){
							FileChecks::process_help($helpfile);														#FileChecks.pm	
	}				
							
	
##### Check megan filename and open file #####
	# Ensure that a filename has been given as an argument. -f flag
	# Get and clean input filename
	my $tophitfile =FileChecks::process_filename($options{f});									#FileChecks.pm
	chomp $tophitfile ; # $tophitfile is the header file.
	
##### Check blast filename and open file #####
	# Ensure that a filename has been given as an argument. -b flag
	# Get and clean input filename
 my $blastfile =FileChecks::process_filename($options{b});									#FileChecks.pm
	chomp $blastfile ;
	
	##### Check see if cap input #####
	# Ensure that a number has been given as an argument. -c flag
	my $capopt =($options{c});																	#FileChecks.pm
	my $cap = 100; # Set it to 100 as a default.
	if ($capopt) { $cap = $capopt;}
	
##### See if extended summary required #####
	#Create extended summary if -e flag called #		
	my $summary;
	if ($options{e}){
							$summary=1;																				#FileChecks.pm	
	}	
										
##### See if genome size adjustment required. If so, check genome size database file and open file #####
	#Do genome size adjustment process if -g flag called #		
	my $adjustment;
	my $genomefilename="Reference_files/All_Genomes_SS.txt";
	if ($options{g}){
							# Check that the filename has been included after -g flag
							# Get and clean input filename	
							$adjustment=1;		
	}	
	
##### Locate nodes.dmp and names.dmp files #####
	# Check that the path has been included after -l flag
	# Get and clean input path
	my $nodesfile="Reference_files/nodes.dmp";
	if ($options{n}){																			#FileChecks.pm
							$nodesfile=($options{n});
	}
	my $namesfile="Reference_files/names.dmp";
	if ($options{N}){																			#FileChecks.pm
							$namesfile=($options{N});
	}
	
##### Check genome size database file and open file #####
	# Check that the filename has been included after -G flag
	# Get and clean input filename	
	if ($options{G}){
							$adjustment=1;
							$genomefilename=FileChecks::process_filename($options{G});								#FileChecks.pm
							chomp $genomefilename;
							# Check file exists
							FileChecks::check_exists $genomefilename;												#FileChecks.pm
	}
	
	##### Check see if thread number input #####
	# Ensure that a number has been given as an argument. -t flag
	my $threadopt =($options{t});																	#FileChecks.pm
	my $thread = 2; # Set to 2 by default.
	if ($threadopt) { $thread = $threadopt;}	


# Start timer file
#-----------------
my $timer_filename = 'timer.txt';
open( my $timer_filehandle, '>', $timer_filename) or die "Cannot open $timer_filename for writing.\n$!\n"; # Note: this will overwrite old timer. I didn't see any point in appending old ones when you need to remove the other PIA outputs before running it again anyway.
use IO::Handle; # Enable autoflush.
$timer_filehandle -> autoflush(1); # Set autoflush to 1 for the log filehandle. This means that Perl won't buffer its output to the log, so the log should be updated in real time.
my $timestamp = localtime();
print $timer_filehandle "Input FASTA: $tophitfile\nInput BLAST: $blastfile\n\nPIA started at $timestamp\n\n";


# Generate a command file to run PIA_inner.pl on x threads
#---------------------------------------------------------
# How many headers will each thread have to process?
my $lines=`wc -l $tophitfile`; # wc returns the number of lines in the header file followed by the name of the header file.
$lines=~ s/^\s+//; # Remove any leading whitespace.
my @lines= split / /,$lines; # Then split on whitespace. So, there will be two elements: the number of lines and the filename.
$lines="$lines[0]"; # Redefine $lines to be just the number of lines. This is equivalent to the number of headers.

$lines=$lines/$thread; # Redefine $lines to be the number of headers to process per thread.
$lines=$lines+1; # Add 1. Does this ensure it's never 0. something?
chomp $lines; # Not sure what the chomp is for. If $lines was a number and a whitespace character, the whitespace would have been removed when we treated it as a numeric and added 1 to it. If it couldn't be treated as a numeric, there would have been an error.
$lines=~ s/\.\d+$//; # $lines is probably now a decimal. Remove any "." characters followed by digits. So, $lines gains 1 but loses whatever fraction it had. It's rounded up to the nearest integer.
#print "Number of headers each thread will process: $lines\n";

# Split the header file accordingly:
system("split -a 3 -l $lines $tophitfile $tophitfile"); # Split the header file such that each new file has $lines lines. The new files will be named [header file].aaa, [header file].aab, [header file].aac...
my @tophitfiles=split/\//,$tophitfile; # Split $tophitfile on "/" symbols in case it's not actually a filename, but a path.
$tophitfile=pop @tophitfiles; # Redefine $tophitfile as the final element, so the actual filename if $tophitfile is a path. If it isn't a path, $tophitfile stays the same.

my $splitfiles =`ls . | grep $tophitfile.`; # "ls ." lists the files and directories in the current directory. The grep narrows this list to files and directories starting with the header filename. It ends up being a list of the split header files.
my @splitfiles= split /\n/,$splitfiles; # Save that list in @splitfiles. To do: rename this to something more meaningful.

# Make the command file.
my $shellscript="shellscript.txt";
unless(open FILE, '>'."$shellscript") { die "\nUnable to create $shellscript\n"; } # If the $shellscript file can't be created, die with an error.

foreach my $file (@splitfiles){ # For each split header file, print to $shellscript the command to run the PIA_inner.pl file with the relevant options we checked earlier.
		print FILE "perl PIA_inner.pl -f $file -b $blastfile -n $nodesfile -N $namesfile &\n";
} # So, if we run these commands simultaneously, we'll be processing each of the split header files simultaneously, and that's threading.
# The fact that the command ends in "&" means that the command is told to run in the background, so new commands can be started on top. The commands should run alongside each other.

close FILE; # Close $shellscript. We're finished editing its contents.

chmod 0755, "$shellscript"; # But still need to change the permissions.

`./$shellscript `; # Finally, run it.


# Collate the PIA_inner.pl outputs from each thread
#---------------------------------------------
my @splitfiles2= `ls . | grep $tophitfile... | grep -v "out" `; # -v means "invert". This line searches for files or directories starting with $tophitfile but that don't contain "out".
my @splitfolders= `ls . | grep $tophitfile... | grep "out" `; # This searches for files or directories starting with $tophitfile that do contain "out". These will be directories that PIA_inner.pl has just made, such as test.headeraaa_out/.

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


mkdir "$tophitfile"."_out"; # Make a master output folder to collect all output in.
# Now make a variable for the filename of each of these PIA_inner.pl output files we've just looked at.
my $OI = $tophitfile."_out/".$tophitfile."_out.intersects.txt"; # So, $OI represents a collated, master out.intersects file inside the output folder. Note that this and $S get "out." because they are more of an output than the log.
my $S = $tophitfile."_out/".$tophitfile."_out.intersects.txt_Summary_Basic.txt";
my $L = $tophitfile."_out/".$tophitfile."_PIA_inner_logs.txt";

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
        next; # If the line contains a hash, which indicates the header line, skip it.
    }
    chomp $line;
	my @line = split("\t", $line);
	my $taxon = $line[0];
	my $hit_count = $line[1];
	if (exists $taxa_and_hits{$taxon}) {
		$taxa_and_hits{$taxon} = $taxa_and_hits{$taxon} + $hit_count;
	} else {
		$taxa_and_hits{$taxon} = $hit_count;
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
my $timer_final_destination = $tophitfile ."_out";
system ("mv $timer_filename $timer_final_destination");

exit;
